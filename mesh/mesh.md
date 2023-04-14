# Generate an ESMF mesh

About
=====

TODO: why this is needed?

Steps
=====

1) Uses the supergrid, ocean mask, and topography to create a netCDF file with the nominal grid of the model:

```
python gen_nc_grid.py
```

2) Convert grid to a SCRIP convention file

```
module load ncl
ncl gen_scrip.ncl
```

3) Create an ESMF mesh

```
qsub create__mesh
```

If these scripts run successfully, you should see file "*_mesh_*.nc" in your directory.

4) Convert the mesh file to CDF5

```
nccopy -k cdf5 tx2_3_mesh_YYMMDD.nc tx2_3_mesh_YYMMDD_cdf5.nc
```
