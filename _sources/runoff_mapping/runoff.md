# Runoff mapping files

In this section, we will generage mapping files used to spread liquid and frozen water into the ocean. [CIME](https://github.com/ESMCI/cime) provides tools to help generate these maps. Please build the runoff_map executable by  following the instructions provided [here](https://github.com/ESMCI/cime/blob/master/tools/mapping/gen_mapping_files/runoff_to_ocn/INSTALL) before proceeding to the next steps.

## JRA-55 (used in C- and G-compsets)

Example of a namelist used to generate  mapping files from JRA-55 to tx2_3 (file map_JRA_to_tx2_3.nml):
```
&input_nml
   gridtype     = 'obs'
   file_roff    = '/glade/p/cesm/omwg_dev/JRA55/domain/domain.190212/domain.roff.JRA025m.190213.nc'
   file_ocn     = '/glade/work/gmarques/cesm/tx2_3/mesh/tx2_3_SCRIP_YYMMDD.nc'
   file_nn      = 'map_jra_to_tx2_3_nn.YYMMDD.nc'
   file_smooth  = 'map_jra_to_tx2_3_sm_e333r100_YYMMDD.nc'
   file_new     = 'map_jra_to_tx2_3_nnsm_e333r100_YYMMDD.nc'
   title        = 'runoff map: JRA-55 -> tx2_3, nearest neighbor and smoothed'
   eFold        = 1000000.0
   rMax         = 300000.0
   step1 = .true.
   step2 = .true.
   step3 = .true.
/
```
> **Note:** Please make sure to modify YYMMDD accordingly .

Where the variables can be divided into four categories:

1) Input grid files
  `gridtype` =  type of file_roff file, "rtm" or "obs" or "scrip"
              * rtm is a 720 x 360 grid ascii file
              * obs is a netcdf file with xc, yc, xv, yv, mask and area variable names
              * scrip is a scrip type grid file (must contain grid_area along with
                typical scrip grid variables)
`file_roff`  =  an ascii rdirc file OR an obs rtm file OR a scrip grid file
`file_ocn`   =  a scrip ocean grid file where the mask is 1 for all ocean grid
                cells (see note 3 below)
`file_ocn_coastal_mask`  =  a scrip ocean grid file where the mask is only 1
                            for coastal grid cells (see note 3 below)

NOTES:
a) gridtype, file_roff, and file_ocn MUST be specified in the namelist
b) if file_ocn_coastal_mask is not specified, file_ocn will be used
c) The file_ocn and file_ocn_coast_mask must be standard scrip grid files that
   include the cell area

2) Input parameters
  `eFold` = smoothing eFold distance in meters (default: 1000000)
  `rMax`  = maximum radius of effect in meters (default: 300000)

3) Settings
  title                 = ascii string to add to mapping files (default: 'unset')
  restrict_smooth_src_to_nn_dest = option to limit the source points for step2 to
                                   just the points that get mapped to in step1; if
                                   false, use all ocean points in
                                   file_ocn_coastal_mask instead (default: .true.)
  `step1`                 = computes nearest neighbor map (default: .true.)
  `step2`                 = computes smooth map (default: .true.)
  `step3`                 = multiply two maps together (default: .true.)

4) Output fields
  `file_nn`     = nearest neighbor mapping file (default: 'nn.nc')
  `file_smooth` = smoother mapping file (default: 'smooth.nc')
  `file_new`    = combined file (default: 'nnsm.nc') - this is the file that will be used.

Example of a scrip used to run **runoff_map** (file JRA_to_ocn.sh):
```
source /glade/work/gmarques/cesm.sandboxes/cesm2_3_beta08_carib/cime/tools/mapping/gen_mapping_files/runoff_to_ocn/src/.env_mach_specific.sh
export TMPDIR=/glade/scratch/gmarques
/glade/work/gmarques/cesm.sandboxes/cesm2_3_beta08_carib/cime/tools/mapping/gen_mapping_files/runoff_to_ocn/runoff_map < map_JRA_to_tx2_3.nml
```
Run the following:

   ```./JRA_to_ocn.sh```

## r05 (used in B-compsets)

Example of a namelist used to generate  mapping files from r05 to tx2_3 (file map_r05_to_tx2_3.nml):

    &input_nml
       gridtype     = 'rtm'
       file_roff    = '/glade/p/cesm/cseg/inputdata/lnd/clm2/rtmdata/rdirc.05.061026'
       file_ocn     = '/glade/work/gmarques/cesm/tx2_3/mesh/tx2_3_SCRIP_230415.nc'
       file_nn      = 'map_r05_to_tx2_3_nn.230415.nc'
       file_smooth  = 'map_r05_to_tx2_3_sm_e1000r1000_230415.nc'
       file_new     = 'map_r05_to_tx2_3_nnsm_e1000r1000_230415.nc'
       title        = 'runoff map: r05 -> tx2_3, nearest neighbor and smoothed'
       eFold        = 1000000.0
       rMax         = 1000000.0
       step1 = .true.
       step2 = .true.
       step3 = .true.
    /

Example of a scrip used to run **runoff_map** (file r05_to_ocn.sh):
```
source /glade/work/gmarques/cesm.sandboxes/cesm2_3_beta08_carib/cime/tools/mapping/gen_mapping_files/runoff_to_ocn/src/.env_mach_specific.sh
export TMPDIR=/glade/scratch/gmarques
/glade/work/gmarques/cesm.sandboxes/cesm2_3_beta08_carib/cime/tools/mapping/gen_mapping_files/runoff_to_ocn/runoff_map < map_r05_to_tx2_3.nml
```
Run the following:

   ```./r05_to_ocn.sh```

## Greenland 4km (used in B-compsets)

Example of a namelist used to generate  mapping files from gland4km to tx2_3 (file map_greenland_4km_to_tx2_3.nml):

    &input_nml
       gridtype     = 'scrip'
       file_roff    = '/glade/p/cesm/cseg/inputdata/glc/cism/griddata/SCRIPgrid_greenland_4km_epsg3413_c170414.nc'
       file_ocn     = '/glade/work/gmarques/cesm/tx2_3/mesh/tx2_3_SCRIP_230415.nc'
       file_nn      = 'map_gland4km_to_tx2_3_nn.230415.nc'
       file_smooth  = 'map_gland4km_to_tx2_3_sm_e1000r300_230415.nc'
       file_new     = 'map_gland4km_to_tx2_3_nnsm_e1000r300_230415.nc'
       title        = 'runoff map: gland4km -> tx2_3, nearest neighbor and smoothed'
       eFold        = 1000000.0
       rMax         = 300000.0
       step1 = .true.
       step2 = .true.
       step3 = .true.
    /

Example of a scrip used to run **runoff_map** (file gland4km_to_ocn.sh):
```
source /glade/work/gmarques/cesm.sandboxes/cesm2_3_beta08_carib/cime/tools/mapping/gen_mapping_files/runoff_to_ocn/src/.env_mach_specific.sh
export TMPDIR=/glade/scratch/gmarques
/glade/work/gmarques/cesm.sandboxes/cesm2_3_beta08_carib/cime/tools/mapping/gen_mapping_files/runoff_to_ocn/runoff_map < map_greenland_4km_to_tx2_3.nml
```
Run the following:

   ```./gland4km_to_ocn.sh```
