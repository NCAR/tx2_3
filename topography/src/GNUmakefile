#
# File:  GNUmakefile
#
#  The bourne shell is safest for gmake
SHELL = /bin/sh

TOP   = .
PARMS_DIR = $(TOP)
UTIL_DIR =  $(TOP)

SRCDIRS = $(TOP) $(PARMS_DIR) $(UTIL_DIR)

#  VPATH is the built-in symbol whose value is the path that gmake will
#  search for dependencies.
VPATH = $(SRCDIRS)

########################################################################
#                  Macro Definitions
########################################################################
#
INCLDIR =
#
FC=    gfortran
FFLAGS = -O2 -g -fopenmp -ffree-form -ffree-line-length-none -fconvert=big-endian $(INCLDIR)
FFLAGS_77 = -O2 -g -fopenmp fconvert=bigendian $(INCLDIR)
#
LDR =   $(FC)
LDFLAGS= $(FFLAGS)
#
LIBDIR =
#
LIBS = $(LIBDIR) -lnetcdff -lnetcdf
#
OBJECTS = kinds.o constants.o ncdf_wrapper.o mom6_grid.o sorting.o
#
########################################################################
#              Non-standard Suffix Rules
########################################################################
#
.SUFFIXES : .f90 .o
.f90.o :
	$(FC) -c $(FFLAGS) $*.f90
#
########################################################################
#            Main Program Build Rules
########################################################################
#
test_grid : $(OBJECTS) test_grid.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $@.o $(OBJECTS) $(LIBS)
create_model_topo : $(OBJECTS) create_model_topo.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $@.o $(OBJECTS) $(LIBS)
create_model_lmask : $(OBJECTS) create_model_lmask.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $@.o $(OBJECTS) $(LIBS)
create_model_ltype : $(OBJECTS) create_model_ltype.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $@.o $(OBJECTS) $(LIBS)
lake_fill : $(OBJECTS) lake_fill.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $@.o $(OBJECTS) $(LIBS)
interp_smooth : $(OBJECTS) interp_smooth.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $@.o $(OBJECTS) $(LIBS)
append_topo_interp_smooth : $(OBJECTS) append_topo_interp_smooth.o
	$(LDR) $(FFLAGS) $(LDFLAGS) -o $@ $@.o $(OBJECTS) $(LIBS)
#
########################################################################
#           Module Dependancies
########################################################################
# Local objects
test_grid.o : $(OBJECTS) test_grid.f90
create_model_topo.o : $(OBJECTS) create_model_topo.f90
create_model_lmask.o : $(OBJECTS) create_model_lmask.f90
create_model_ltype.o : $(OBJECTS) create_model_ltype.f90
lake_fill.o : $(OBJECTS) lake_fill.f90
append_topo_interp_smooth.o : $(OBJECTS) append_topo_interp_smooth.f90
#
# Utilities
kinds.o :
	$(FC) -c $(FFLAGS) $(UTIL_DIR)/kinds.f90
constants.o : kinds.o
	$(FC) -c $(FFLAGS) $(UTIL_DIR)/constants.f90
ncdf_wrapper.o :
	$(FC) -c $(FFLAGS) $(UTIL_DIR)/ncdf_wrapper.f90
mom6_grid.o :
	$(FC) -c $(FFLAGS) $(UTIL_DIR)/mom6_grid.f90
sorting.o :
	$(FC) -c $(FFLAGS) $(UTIL_DIR)/sorting.f90

###########################################################################
clean:
	rm -f *.o *.mod *.MOD
clobber: clean
	rm -f create_model_topo create_model_lmask lake_fill interp_smooth
