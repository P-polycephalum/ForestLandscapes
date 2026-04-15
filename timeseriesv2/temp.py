import zarr
import dask.array as da
import napari
import os
import numpy as np

cube_root = r"D:\BCI_50ha_timeseries\tiles"
tile_type = "aligned_local"
tile_id = "0_0"
path = os.path.join(cube_root, tile_type, tile_id, "cube.zarr")
print(f"Opening {path}...")

os.path.exists(path)
z = zarr.open(path, mode="r")

times = list(z.attrs.get("time", []))
print(f"Shape: {z.shape}  —  {len(times)} time steps")

# (T, C, H, W) -> (T, H, W, C), keep RGB only
rgb = da.from_zarr(z)[:, :3, :, :]
rgb = da.moveaxis(rgb, 1, -1)

sample = rgb[0].compute()
p2, p98 = np.percentile(sample, [2, 98])

viewer = napari.Viewer()
viewer.add_image(
    rgb,
    rgb=True,
    name=f"{tile_type}/{tile_id}",
    contrast_limits=(float(p2), float(p98)),
)

# Show time label in title on slider move
def _on_step(event):
    t = viewer.dims.current_step[0]
    label = times[t] if t < len(times) else str(t)
    viewer.title = f"t={t}  {label}"

viewer.dims.events.current_step.connect(_on_step)
_on_step(None)

napari.run()
