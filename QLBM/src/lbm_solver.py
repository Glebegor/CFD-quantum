"""Classical D2Q9 BGK reference solver for macroscopic CFD fields.

This solver provides the complete population field that the current
single-basis-state QPU streaming experiment does not yet encode.
"""

import json
import shutil
from pathlib import Path
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np

from animation import create_time_series_animation
from macrocomputations import CS2, CX, CY, reconstruct_fields


WEIGHTS = np.array(
    [4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36],
    dtype=float,
)


def equilibrium(density, velocity):
    """Return D2Q9 equilibrium populations with shape ``(9, ny, nx)``."""
    density = np.asarray(density, dtype=float)
    velocity = np.asarray(velocity, dtype=float)
    if velocity.shape != (*density.shape, 2):
        raise ValueError("velocity must have shape (*density.shape, 2)")

    cu = (
        CX[:, None, None] * velocity[None, ..., 0]
        + CY[:, None, None] * velocity[None, ..., 1]
    )
    u2 = np.sum(velocity * velocity, axis=-1)
    return (
        WEIGHTS[:, None, None]
        * density[None, ...]
        * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2[None, ...])
    )


class D2Q9Solver:
    """BGK solver with left velocity inlet, right outlet and periodic y."""

    def __init__(
        self,
        nx,
        ny,
        inlet_velocity_lattice=0.05,
        relaxation_time=0.6,
        initial_velocity_lattice=0.0,
    ):
        if nx < 3 or ny < 2:
            raise ValueError("the domain must be at least 3x2")
        if not 0.0 < inlet_velocity_lattice < 0.1:
            raise ValueError("inlet lattice velocity must be in (0, 0.1)")
        if relaxation_time <= 0.5:
            raise ValueError("relaxation_time must be greater than 0.5")

        self.nx = nx
        self.ny = ny
        self.inlet_velocity_lattice = inlet_velocity_lattice
        self.relaxation_time = relaxation_time
        self.omega = 1.0 / relaxation_time

        density = np.ones((ny, nx), dtype=float)
        velocity = np.zeros((ny, nx, 2), dtype=float)
        velocity[..., 0] = initial_velocity_lattice
        self.f = equilibrium(density, velocity)

    def step(self):
        """Apply one BGK collision and periodic streaming step."""
        density, velocity, _ = reconstruct_fields(self.f)
        f_equilibrium = equilibrium(density, velocity)
        post_collision = self.f - self.omega * (self.f - f_equilibrium)

        streamed = np.empty_like(post_collision)
        for direction, (dx, dy) in enumerate(zip(CX.astype(int), CY.astype(int))):
            streamed[direction] = np.roll(
                post_collision[direction],
                shift=(dy, dx),
                axis=(0, 1),
            )

        # Left velocity inlet: impose rho=1 and u=(u_in, 0) by equilibrium.
        inlet_density = np.ones((self.ny, 1), dtype=float)
        inlet_velocity = np.zeros((self.ny, 1, 2), dtype=float)
        inlet_velocity[..., 0] = self.inlet_velocity_lattice
        streamed[:, :, :1] = equilibrium(inlet_density, inlet_velocity)

        # Right zero-gradient outlet. The y direction remains periodic.
        streamed[:, :, -1] = streamed[:, :, -2]
        self.f = streamed
        return self.f

    def run_steps(self, count):
        if count <= 0:
            raise ValueError("step count must be positive")
        for _ in range(count):
            self.step()
        return self.f


def physical_fields(
    f,
    reference_density_kg_m3,
    reference_pressure_pa,
    velocity_scale_m_s,
):
    """Convert lattice populations into SI density, velocity and pressure."""
    rho_lattice, velocity_lattice, pressure_lattice = reconstruct_fields(f)
    density = reference_density_kg_m3 * rho_lattice
    velocity = velocity_scale_m_s * velocity_lattice
    gauge_pressure = (
        reference_density_kg_m3
        * velocity_scale_m_s**2
        * (pressure_lattice - CS2)
    )
    pressure = reference_pressure_pa + gauge_pressure
    return {
        "density": density,
        "velocity": velocity,
        "pressure": pressure,
        "gauge_pressure": gauge_pressure,
        "density_lattice": rho_lattice,
        "velocity_lattice": velocity_lattice,
        "pressure_lattice": pressure_lattice,
    }


def save_snapshot(
    output_dir,
    step,
    time_seconds,
    f,
    fields,
    metadata,
    history,
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = output_dir / f"step_{step:04d}"

    speed = np.linalg.norm(fields["velocity"], axis=-1)
    figure, axes = plt.subplots(3, 3, figsize=(18, 14))
    plots = (
        (fields["density"], "Density", "kg/m³", "viridis"),
        (fields["pressure"], "Absolute pressure", "Pa", "coolwarm"),
        (fields["gauge_pressure"], "Gauge pressure", "Pa", "coolwarm"),
        (speed, "Velocity magnitude", "m/s", "magma"),
        (fields["velocity"][..., 0], "Horizontal velocity $u_x$", "m/s", "coolwarm"),
        (fields["velocity"][..., 1], "Vertical velocity $u_y$", "m/s", "coolwarm"),
    )
    for axis, (values, title, label, cmap) in zip(axes.flat[:6], plots):
        image = axis.imshow(values, origin="lower", cmap=cmap, aspect="auto")
        axis.set_title(title)
        axis.set_xlabel("x cell")
        axis.set_ylabel("y cell")
        figure.colorbar(image, ax=axis, label=label)

    # Subsample the vector field to keep the plot readable on large grids.
    stride = max(1, min(f.shape[1], f.shape[2]) // 16)
    yy, xx = np.mgrid[0 : f.shape[1] : stride, 0 : f.shape[2] : stride]
    axes[2, 0].quiver(
        xx,
        yy,
        fields["velocity"][::stride, ::stride, 0],
        fields["velocity"][::stride, ::stride, 1],
        speed[::stride, ::stride],
        cmap="plasma",
        angles="xy",
        scale_units="xy",
        scale=1,
    )
    axes[2, 0].set_title("Velocity vector field")
    axes[2, 0].set_xlabel("x cell")
    axes[2, 0].set_ylabel("y cell")
    axes[2, 0].set_xlim(-0.5, f.shape[2] - 0.5)
    axes[2, 0].set_ylim(-0.5, f.shape[1] - 0.5)
    axes[2, 0].set_aspect("equal")
    axes[2, 0].grid(alpha=0.2)

    population_totals = np.sum(f, axis=(1, 2))
    axes[2, 1].bar(np.arange(9), population_totals, color="tab:blue")
    axes[2, 1].set_title("Total D2Q9 populations")
    axes[2, 1].set_xlabel("direction")
    axes[2, 1].set_ylabel(r"$\sum_{x,y} f_i$")
    axes[2, 1].set_xticks(np.arange(9))

    times = [item["time_seconds"] for item in history]
    mean_density = [item["mean_density_kg_m3"] for item in history]
    max_speed = [item["max_velocity_m_s"] for item in history]
    pressure_range = [item["pressure_range_pa"] for item in history]
    history_axis = axes[2, 2]
    density_line = history_axis.plot(
        times,
        mean_density,
        label="mean density [kg/m³]",
    )
    speed_line = history_axis.plot(
        times,
        max_speed,
        label="max speed [m/s]",
    )
    pressure_axis = history_axis.twinx()
    pressure_line = pressure_axis.plot(
        times,
        pressure_range,
        color="tab:green",
        label="pressure range [Pa]",
    )
    history_axis.set_title("Time-series diagnostics")
    history_axis.set_xlabel("time [s]")
    history_axis.set_ylabel("density / speed")
    pressure_axis.set_ylabel("pressure range [Pa]", color="tab:green")
    history_axis.set_xlim(0.0, max(time_seconds, 0.1))
    history_axis.grid(alpha=0.25)
    lines = density_line + speed_line + pressure_line
    history_axis.legend(
        lines,
        [line.get_label() for line in lines],
        fontsize=8,
    )

    figure.suptitle(f"D2Q9 BGK fields — t={time_seconds:.3f} s")
    figure.tight_layout()
    figure.savefig(stem.with_suffix(".png"), dpi=150, bbox_inches="tight")
    plt.close(figure)

    np.savez_compressed(
        stem.with_suffix(".npz"),
        f=f,
        time_seconds=time_seconds,
        **fields,
    )
    stem.with_suffix(".json").write_text(
        json.dumps(metadata, indent=2) + "\n",
        encoding="utf-8",
    )


def run_cfd_series(params):
    """Run the configured physical inlet-flow series and save all fields."""
    cfd_cells = int(getattr(params, "CFD_CELLS_SIZE", params.CELLS_SIZE))
    nx = cfd_cells
    ny = cfd_cells
    cell_size = float(getattr(params, "PHYSICAL_CELL_SIZE_M", 0.25))
    inlet_velocity = float(
        getattr(params, "PHYSICAL_INLET_VELOCITY_M_S", 1.0)
    )
    lattice_velocity = float(getattr(params, "LBM_INLET_VELOCITY", 0.05))
    output_dt = float(getattr(params, "TIME_STEP", 0.1))
    snapshot_count = int(getattr(params, "SNAPSHOT_COUNT", 40))
    relaxation_time = float(getattr(params, "LBM_RELAXATION_TIME", 0.6))
    reference_density = float(
        getattr(params, "REFERENCE_DENSITY_KG_M3", 1.225)
    )
    reference_pressure = float(
        getattr(params, "REFERENCE_PRESSURE_PA", 101325.0)
    )

    velocity_scale = inlet_velocity / lattice_velocity
    internal_dt = cell_size / velocity_scale
    substeps_float = output_dt / internal_dt
    substeps = round(substeps_float)
    if substeps < 1 or not np.isclose(substeps_float, substeps):
        raise ValueError(
            "TIME_STEP must be an integer multiple of the stable internal "
            f"LBM step ({internal_dt:g} s)"
        )

    output_dir = Path(
        getattr(params, "CFD_TIME_SERIES_DIR", "output/cfd_series")
    )
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)

    solver = D2Q9Solver(
        nx=nx,
        ny=ny,
        inlet_velocity_lattice=lattice_velocity,
        relaxation_time=relaxation_time,
    )
    started = perf_counter()
    snapshots = []
    for snapshot in range(1, snapshot_count + 1):
        solver.run_steps(substeps)
        time_seconds = snapshot * output_dt
        fields = physical_fields(
            solver.f,
            reference_density,
            reference_pressure,
            velocity_scale,
        )
        metadata = {
            "schema": "d2q9-bgk-physical-snapshot-v1",
            "snapshot": snapshot,
            "time_seconds": time_seconds,
            "grid": [ny, nx],
            "cell_size_m": cell_size,
            "inlet_velocity_m_s": inlet_velocity,
            "lattice_inlet_velocity": lattice_velocity,
            "internal_time_step_seconds": internal_dt,
            "lbm_substeps": substeps,
            "relaxation_time": relaxation_time,
            "reference_density_kg_m3": reference_density,
            "reference_pressure_pa": reference_pressure,
            "boundary_conditions": {
                "left": "equilibrium velocity inlet",
                "right": "zero-gradient outlet",
                "top_bottom": "periodic",
            },
            "density_min_max": [
                float(fields["density"].min()),
                float(fields["density"].max()),
            ],
            "pressure_min_max_pa": [
                float(fields["pressure"].min()),
                float(fields["pressure"].max()),
            ],
            "velocity_magnitude_min_max_m_s": [
                float(np.linalg.norm(fields["velocity"], axis=-1).min()),
                float(np.linalg.norm(fields["velocity"], axis=-1).max()),
            ],
            "mean_density_kg_m3": float(fields["density"].mean()),
            "max_velocity_m_s": float(
                np.linalg.norm(fields["velocity"], axis=-1).max()
            ),
            "pressure_range_pa": float(
                fields["pressure"].max() - fields["pressure"].min()
            ),
        }
        snapshots.append(metadata)
        save_snapshot(
            output_dir,
            snapshot,
            time_seconds,
            solver.f,
            fields,
            metadata,
            snapshots,
        )
        print(
            f"CFD snapshot {snapshot}/{snapshot_count}: "
            f"t={time_seconds:.3f} s, "
            f"max|u|={metadata['velocity_magnitude_min_max_m_s'][1]:.4f} m/s"
        )

    manifest = {
        "schema": "d2q9-bgk-physical-series-v1",
        "snapshot_count": snapshot_count,
        "output_time_step_seconds": output_dt,
        "duration_seconds": snapshot_count * output_dt,
        "internal_time_step_seconds": internal_dt,
        "substeps_per_snapshot": substeps,
        "wall_time_seconds": perf_counter() - started,
        "grid": [ny, nx],
        "snapshots": snapshots,
    }
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    gif, mp4 = create_time_series_animation(output_dir, output_dt)
    print(f"CFD manifest: {manifest_path}")
    print(f"CFD GIF: {gif}")
    if mp4:
        print(f"CFD MP4: {mp4}")
    return manifest


if __name__ == "__main__":
    import params

    run_cfd_series(params)
