# Supergrid

About
=====
The supergrid is created using ORCA_gridgen which relies on the NEMO ocean modelling framework to generate tripolar grids for MOM6.
For a complete description of the code, please check [this user guide](ORCA_gridgen/tripole-user-guide.pdf). The
original code has been modified to create an supergrid (ocean_hgrid.nc), in addition to the
defaults files for NEMO (coordinates.nc and coordinates_north.nc).

TODO: describe what else was modified... parallel, u point at the Equator, etc

Usage
=====

File param.f90 has the modifications needed to generate the nominal 2/3 degree resolution grid (tx2_3);
The original (default) configuration used for ORCA1 is also included
for comparison (file param.f90.ori and trop.f90.ori). A simple ``diff`` of the files can be used to check
what has been changed to create tx2_3:

```
diff param.f90 param.f90.ori
```

The code is intended to compile under any compliant Fortran 90 compiler. It
must be compiled with a flag that promotes reals to 8 bytes. It must be linked
with the NetCDF F90 library. To compile ORCA_gridgen on Casper, load the following modules:

```
module load intel/19.1.1
module load netcdf/4.8.1
```

Then type:
```
cd ORCA_gridgen
make clean
make
```

This should create an executable called "tripole.exe". Next, type the following:

```
./tripole.exe
```

This should create three NetCDF files: ocean_hgrid.nc, coordinates.nc and coordinates_north.nc. For MOM6, we only use ocean_hgrid.nc.
