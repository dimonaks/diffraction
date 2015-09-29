

CPP=g++
F90=gfortran
# CPP=gcc-mp-4.9
# F90=gfortran-mp-4.9
# CFLAGS = -I include -std=gnu99 -m64
CFLAGS = -undefined dynamic_lookup
FFLAGS = -ffixed-line-length-none -lgfortran
all: diffraction
	mv ./diffraction ..
	# cd ..; ./diffraction list

diffraction: diffraction.o spline.o

	${CPP} $(CFLAGS)  diffraction.o spline.o   -o diffraction

diffraction.o: diffraction.cpp diffraction.h
	
	${CPP} $(CFLAGS) -c  diffraction.cpp


spline.o: fortran/spline.f
	
	${F90} $(FFLAGS) -c  fortran/spline.f 

clean:
	rm *o


