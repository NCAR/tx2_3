# ------------------Macro-Defs---------------------
F90=ifort -r8 -qopenmp
NC_CONFIG = nc-config
NETCDF_INC = -I $(shell $(NC_CONFIG) --prefix)/include
NETCDF_LIB = $(shell $(NC_CONFIG) --flibs)
LIBS = $(NETCDF_LIB)
# -------------------End-macro-Defs---------------------------

# Here is the link step 
tripole.exe:kinds.o nf90util.o param.o common.o functions.o nthcar.o hgrnth.o hgrctl.o trop.o shoots.o mshtrp.o mshnth.o mshsth.o gridio.o tripole.o 
	 $(F90) -o tripole.exe kinds.o nf90util.o param.o common.o functions.o nthcar.o hgrnth.o hgrctl.o trop.o shoots.o mshtrp.o mshnth.o mshsth.o gridio.o tripole.o   $(LIBS) 

# Here are the compile steps
 
kinds.o:./kinds.F90  
	 $(F90) -c -Dkey_64bit ./kinds.F90 $(NETCDF_INC) 

nf90util.o:./nf90util.f90 kinds.o 
	 $(F90) -c ./nf90util.f90 $(NETCDF_INC) 

param.o:./param.f90 kinds.o 
	 $(F90) -c ./param.f90 $(NETCDF_INC) 

common.o:./common.f90 kinds.o param.o 
	 $(F90) -c ./common.f90 $(NETCDF_INC) 

functions.o:./functions.f90 kinds.o param.o common.o 
	 $(F90) -c ./functions.f90 $(NETCDF_INC) 

nthcar.o:./nthcar.f90 kinds.o param.o common.o functions.o 
	 $(F90) -c ./nthcar.f90 $(NETCDF_INC) 

hgrnth.o:./hgrnth.f90 kinds.o param.o common.o functions.o nthcar.o 
	 $(F90) -c ./hgrnth.f90 $(NETCDF_INC) 

hgrctl.o:./hgrctl.f90 kinds.o param.o common.o functions.o 
	 $(F90) -c ./hgrctl.f90 $(NETCDF_INC) 

trop.o:./trop.f90 kinds.o param.o 
	 $(F90) -c ./trop.f90 $(NETCDF_INC) 

shoots.o:./shoots.f90 kinds.o trop.o 
	 $(F90) -c ./shoots.f90 $(NETCDF_INC) 

mshtrp.o:./mshtrp.f90 kinds.o nf90util.o shoots.o trop.o common.o param.o 
	 $(F90) -c ./mshtrp.f90 $(NETCDF_INC) 

mshnth.o:./mshnth.f90 kinds.o param.o common.o functions.o 
	 $(F90) -c ./mshnth.f90 $(NETCDF_INC) 

mshsth.o:./mshsth.f90 param.o common.o functions.o 
	 $(F90) -c ./mshsth.f90 $(NETCDF_INC) 

gridio.o:./gridio.f90 param.o common.o nf90util.o 
	 $(F90) -c ./gridio.f90 $(NETCDF_INC) 

tripole.o:tripole.f90 kinds.o nf90util.o param.o common.o functions.o hgrnth.o hgrctl.o mshtrp.o mshnth.o mshsth.o gridio.o 
	 $(F90) -c tripole.f90 $(NETCDF_INC) 
# This entry allows you to type " make clean " to get rid of
# all object and module files 
clean:
	rm -f -r f_{files,modd}* *.o *.mod *.M *.d V*.inc *.vo \
	V*.f *.dbg album F.err
  
