"""Create portable animations from persisted QLBM snapshot plots."""

import json
import shutil
import subprocess
from pathlib import Path

from PIL import Image


def create_time_series_animation(output_dir, time_step=0.1):
    """Create GIF and, when ffmpeg is available, MP4 from step PNG files."""
    output_dir = Path(output_dir)
    frame_paths = sorted(output_dir.glob("step_*.png"))
    if not frame_paths:
        raise FileNotFoundError(f"no step_*.png frames found in {output_dir}")
    if time_step <= 0:
        raise ValueError("time_step must be positive")

    gif_path = output_dir / "qlbm_animation.gif"
    frames = []
    for frame_path in frame_paths:
        with Image.open(frame_path) as image:
            frame = image.convert("RGB")
            if frame.width > 1000:
                height = round(frame.height * 1000 / frame.width)
                frame = frame.resize(
                    (1000, height),
                    Image.Resampling.LANCZOS,
                )
            frames.append(
                frame.quantize(colors=128, method=Image.Quantize.MEDIANCUT)
            )
    frames[0].save(
        gif_path,
        save_all=True,
        append_images=frames[1:],
        duration=round(time_step * 1000),
        loop=0,
        optimize=False,
        disposal=2,
    )

    mp4_path = None
    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg:
        candidate = output_dir / "qlbm_animation.mp4"
        subprocess.run(
            [
                ffmpeg,
                "-y",
                "-framerate",
                f"{1.0 / time_step:g}",
                "-i",
                str(output_dir / "step_%04d.png"),
                "-vf",
                "pad=ceil(iw/2)*2:ceil(ih/2)*2",
                "-c:v",
                "libx264",
                "-pix_fmt",
                "yuv420p",
                "-movflags",
                "+faststart",
                str(candidate),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        mp4_path = candidate

    manifest_path = output_dir / "manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["animation"] = {
            "gif": gif_path.name,
            "mp4": None if mp4_path is None else mp4_path.name,
            "frames_per_second": 1.0 / time_step,
            "frame_count": len(frame_paths),
        }
        manifest_path.write_text(
            json.dumps(manifest, indent=2) + "\n",
            encoding="utf-8",
        )

    return gif_path, mp4_path


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", nargs="?", default="output/time_series")
    parser.add_argument("--time-step", type=float, default=0.1)
    arguments = parser.parse_args()
    gif, mp4 = create_time_series_animation(
        arguments.output_dir,
        arguments.time_step,
    )
    print(f"GIF saved to: {gif}")
    if mp4:
        print(f"MP4 saved to: {mp4}")
