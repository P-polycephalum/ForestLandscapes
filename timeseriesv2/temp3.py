import zarr
from tqdm import tqdm

path=r"D:\BCI_50ha_timeseries\tiles\aligned_global\0_0\cube.zarr"
path2=r"D:\BCI_50ha_timeseries\tiles\aligned_global\0_0\cube2.zarr"

z = zarr.open(path, mode='r')

cube = zarr.open(
        path2,
        mode="w",
        shape=(z.shape[0], z.shape[1], z.shape[2], z.shape[3]),
        chunks=(1, z.shape[1], z.chunks[2], z.chunks[3]),
        dtype=z.dtype,
        zarr_format=2,
    )

cube.attrs.update(z.attrs.asdict())

# Copy data one time step at a time
for t in tqdm(range(z.shape[0]), desc="Copying"):
    cube[t] = z[t]

print("Done.")
