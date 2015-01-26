SHELL = /bin/bash
FC        := gfortran
ld        := $(FC)
moddir    := modules
lib       := ibridgelib.a
FCFLAGS   := -O3 -w -ffree-line-length-none -fopenmp -openmp
#FCFLAGS   := -O3 -ffree-line-length-none -fopenmp -openmp
FFTWFLG   := -lfftw3 -lfftw3_threads -lfftw3_omp -lm
NETCDFLIB := -L/usr/local/lib
NETCDFINC := -I/usr/local/include
NETCDFFLG := -lnetcdff -lnetcdf
ARFLAGS   := -crv
F77OBJS   := dverkb.f atmosdata.f iridata.f
F90OBJS   := fortgrapher.f90 mathsubs.f90 sysinfo.f90 sysgrids.f90 geoaxes.f90 doperators.f90 deset.f90 syssolver.f90 bridge.f90
%.o: %.f
	$(FC) $(FCFLAGS) -c $< -o $@
%.o: %.f90
	$(FC) $(FCFLAGS) $(FFTWFLG) $(NETCDFINC) $(NETCDFLIB) $(NETCDFFLG) -c $< -o $@
All: main
main: ibridge
ibridge : $(lib)(bridge.o)
	ar x $(lib) bridge.o
	$(FC) -p -o ibridge bridge.o $(lib) $(FFTWFLG) $(NETCDFLIB) $(NETCDFFLG) -liomp5 -lpthread
	rm bridge.o
$(lib)(bridge.o) : $(lib)(iridata.o) $(lib)(atmosdata.o) $(lib)(dverkb.o) \
                   $(lib)(rdsvnetcdf.o) $(lib)(mathsubs.o) $(lib)(fortgrapher.o) $(lib)(interpsubs.o) $(lib)(sysinfo.o) \
                   $(lib)(geoaxes.o) $(lib)(sysgrids.o) $(lib)(ionobg.o) \
		   $(lib)(doperators.o) $(lib)(deset.o) $(lib)(syssolver.o) $(lib)(sysgraphics.o) \
		   $(lib)(itest.o)

#$(moddir) :
#	mkdir -p $(moddir)
cleanall:
	rm -fr $(libobjs) *.o *.mod ibridge ibridgelib.a gmon.out

cleanpy:
	rm -fr *.dat *.py FES FLDENERGY* PPROFILE* GAMMA* OMEGA* *Field* GVProfile* *.pdf Data

cleanjob:
	rm -fr bridge.o* bridge.e* gmon.out 

cleanreset:
	rm -fr energySpectrum.out restartdata.out gmon.out
