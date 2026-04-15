import zarr, importlib.metadata, pathlib
print("zarr version:", zarr.__version__)
print("zarr location:", pathlib.Path(zarr.__file__).resolve())

# Check if it came from conda or pip
try:
    dist = importlib.metadata.distribution("zarr")
    installer = dist.read_text("INSTALLER") or "unknown"
    print("Installed via:", installer.strip())
except Exception as e:
    print("Could not determine installer:", e)
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.show()


###################