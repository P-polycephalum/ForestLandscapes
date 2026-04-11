import os
import numpy as np
import zarr
import imageio.v2 as imageio

zarr_path = r"D:\BCI_50ha_timeseries\tiles\aligned_global\0_1\cube.zarr"
out_mp4 = r"D:\BCI_50ha_timeseries\tiles\aligned_global\0_1\timeseries.mp4"
fps = 6

z = zarr.open(zarr_path, mode="r")  # shape: (time, band, y, x)
n_frames = z.shape[0]

def to_uint8_rgb(frame_chw):
    # frame_chw: (band, y, x), use first 3 bands
    rgb = np.moveaxis(frame_chw[:3], 0, -1).astype(np.float32)  # (y, x, 3)

    # robust normalization to 8-bit
    lo = np.percentile(rgb, 2)
    hi = np.percentile(rgb, 98)
    if hi <= lo:
        hi = lo + 1.0
    rgb = np.clip((rgb - lo) / (hi - lo), 0, 1)
    return (rgb * 255).astype(np.uint8)

os.makedirs(os.path.dirname(out_mp4), exist_ok=True)

with imageio.get_writer(out_mp4, codec="libx264", quality=8) as writer:
    for i in range(n_frames):
        frame = z[i]  # (band, y, x)
        writer.append_data(to_uint8_rgb(frame))
        print(f"Frame {i+1}/{n_frames}")

print(f"\nSaved video: {out_mp4}")