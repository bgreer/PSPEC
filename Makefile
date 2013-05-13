
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

LINK=-I$(FFTWINC) -L$(FFTWLIB) -lfftw3 -L$(FITSLIB) -lcfitsio

pspec : pspec.f90 ParseInput.f90
	$(COMPILER) -c ParseInput.f90 $(LINK)
	$(COMPILER) pspec.f90 ParseInput.f90 -o pspec $(LINK)


clean :
	rm pspec *.mod

