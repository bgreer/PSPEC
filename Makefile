
COMPILER=ifort -O3 -ipo -openmp
COMPILER=gfortran -O3 -fopenmp

# FFTW library

FFTWLOC=/usr/
#FFTWLOC=/shared/fftw-3.2.1
FFTWLIB=$(FFTWLOC)/lib
FFTWINC=$(FFTWLOC)/include

# cfitsio library
FITSLOC=/usr/local
#FITSLOC=/shared/cfitsio
FITSLIB=$(FITSLOC)/lib
FITSINC=$(FITSLOC)/include

pspec : pspec.f90
	$(COMPILER) -o pspec pspec.f90 -I$(FFTWINC) -L$(FFTWLIB) -lfftw3 -L$(FITSLIB) -lcfitsio

clean :
	rm pspec

