# Topography generation

In this section, we will walk through the steps required to generate an accurate representation of ocean topography. These steps can be summarized as follows: **1)** generate a land/ocean mask; **2)** visualy inspect and apply manual modifications to the mask via Jupyter Notebook; and **3)** generate a smoothed ocean topography using the land/sea mask and a high-resolution global topography dataset.

> code to generate mask:  create_model_topo.f90
> code to generate smoothed topography: interp_smooth.f90

a couple utility modules:
   -constants.f90
   -kinds.f90
   -mom6_grid.f90
   -ncdf_wrapper.f90

a make file : GNUmakefile
> **Note:** Expects gnu compiler as your default. On Casper, "module load gnu".

PBS script files to run the Fortran codes on Casper:
   -run_create_topo.tx1_4v2.pbs
   -run_intep_smooth_tx1_4v2.pbs

> **Note:** Namelist input for the program is set in these.

Notebook to hand edit mask: MaskEdit_tx2_3v2b.ipynb
   -utility script:  map_mask.py

The global bathymetry and topography at 15 arc sec derived from the Shuttle Radar Topography Mission (SRTM) can be found at:

    /glade/campaign/cgd/oce/datasets/obs/topography/SRTM/SRTM15_V2.4.nc

For more information on this dataset, see its associated publication:

Tozer, B, Sandwell, D. T., Smith, W. H. F., Olson, C., Beale, J. R., & Wessel, P. (2019). Global bathymetry and topography at 15 arc sec: SRTM15+. Earth and Space Science, 6, 1847. https://doi.org/10.1029/2019EA000658

## Compile Fortran codes
Modules that must be loaded:
- ncarenv/1.3
- gnu/10.1.0
- ncarcompilers/0.5.0
- netcdf/4.8.1
- mpt/2.19

Run the following:
```
gmake create_model_topo
gmake inter_smooth
```
## Land/sea mask and topo stats
This step will generate land/ocean mask information (mask - 0 if land, 1 if ocean - and fraction of cell area covered by ocean) and ocean topography statistics (mean, median, mean squared, minimum, and maximun) on MOM6 grid provided.

Run the following:

   ```qsub run_create_topo.tx2_3v2.pbs```

> **Note:** Open run_create_topo.tx2_3v2.pbs and make sure the namelist input is correct.

This should create file **topo.sub.{NSUB}.{GRID}.srtm.nc**, where NSUB and GRID are defined in run_create_topo.tx2_3v2.pbs.

## Check the land/sea mask
This step consists of visualy inspect the land/sea mask and possibly apply manual modifications. To do so,  open MaskEdit_tx2_3v2b.ipynb and execute all the cells. Modifications to the land/sea mask should be done in this Notebook. For example, if one wants to convert an ocean point around Antarctica to land, navigate to the "Antarctica" section and make the necessary changes there. Please make sure to add notes about the changes made.

## Generate smoothed topography
This step will generate a smoothed topography using the land/sea mask create in the previous step and a high-resolution global topography dataset.

Run the following:

   ```qsub run_intep_smooth_tx2_3v2.pbs```

> **Note:** Open run_intep_smooth_tx2_3v2.pbs and make sure the namelist input is correct.

This step will gerate file **topo.{GRID}.{TOPODATA}.{EDIT}.{SFNC}.nc**, where GRID,  TOPODATA, EDIT and SFNC are defined in run_intep_smooth_tx2_3v2.pbs. This is the final topography!

TODO: add instructions on how to modify the topography.

