# Generate an ESMF mesh

Steps
=====

1) Uses the supergrid, ocean mask, and topography to create a netCDF file with the nominal grid of the model:

```
python gen_nc_grid.py
```
> **Note:** Make sure variable ```nc_topo``` in gen_nc_grid.py is pointing to the right topography file .


2) Convert grid to a SCRIP convention file

```
module load ncl
ncl gen_scrip.ncl
```
> **Note:** Make sure to set variables ncGrdFilePath, maskFilePath, and dstGridPath in ```gen_scrip.ncl```.

3) Create an ESMF mesh

```
qsub create__mesh
```
> **Note:** Make sure to set variables file_i and file_o ```create__mesh```.

If these scripts run successfully, you should see file "tx?_?_mesh_YYMMDD.nc" in your directory.

4) Convert the mesh file to CDF5

```
nccopy -k cdf5 tx2_3_mesh_YYMMDD.nc tx2_3_mesh_YYMMDD_cdf5.nc
```
