"""Macroscopic D2Q9 field reconstruction shared by scripts and notebooks."""

import numpy as np


CX = np.array([0, 1, 0, -1, 0, 1, -1, -1, 1], dtype=float)
CY = np.array([0, 0, 1, 0, -1, 1, 1, -1, -1], dtype=float)
CS2 = 1.0 / 3.0


def reconstruct_fields(f):
    """Return density, velocity and isothermal pressure from D2Q9 populations."""
    populations = np.asarray(f, dtype=float)
    if populations.ndim != 3:
        raise ValueError(
            "f must be a 3D array with shape (9, ny, nx) or (ny, nx, 9)"
        )
    if populations.shape[0] == 9:
        fq = populations
    elif populations.shape[-1] == 9:
        fq = np.moveaxis(populations, -1, 0)
    else:
        raise ValueError("one axis of f must contain exactly 9 populations")
    if not np.all(np.isfinite(fq)):
        raise ValueError("f contains NaN or infinite values")

    density = np.sum(fq, axis=0)
    momentum_x = np.einsum("q,qyx->yx", CX, fq)
    momentum_y = np.einsum("q,qyx->yx", CY, fq)
    velocity = np.zeros((*density.shape, 2), dtype=float)
    np.divide(
        momentum_x,
        density,
        out=velocity[..., 0],
        where=density != 0.0,
    )
    np.divide(
        momentum_y,
        density,
        out=velocity[..., 1],
        where=density != 0.0,
    )
    pressure = CS2 * density
    return density, velocity, pressure
