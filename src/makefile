

CPP=g++
F90=gfortran
# CPP=gcc-mp-4.9
# F90=gfortran-mp-4.9
export LD_LIBRARY_PATH=/opt/local/lib/gcc49
# CFLAGS = -I include -std=gnu99 -m64
CFLAGS = 
all: diffraction
	mv ./diffraction ..
	# cd ..; ./diffraction list

diffraction: diffraction.o spline.o
	# export LD_LIBRARY_PATH=/opt/local/lib/gcc49

	${CPP} $(CFLAGS) -undefined dynamic_lookup diffraction.o spline.o   -o diffraction

diffraction.o: diffraction.cpp
	${CPP} $(CFLAGS)   -c diffraction.cpp

spline.o: fortran/spline.f
	${F90} -c -ffixed-line-length-none \
	fortran/spline.f -lgfortran

clean:
	rm *o

# cd src
# ${CPP} -c diffraction.cpp

# # cd fortran
# ${F90} -c -ffixed-line-length-none \
# fortran/spline.f

# # cd ..

 
# ${CPP} -o  diffraction diffraction.o spline.o

# cp diffraction ../


# VPATH   = src include
# CFLAGS  = -I include -std=gnu99 -m64
# CC      = gcc
# LDLIBS  =
# LDFLAGS =

# hello: hello.o
#     $(CC) $(CFLAGS) $^ -o $@

# hello.o: hello.c hello.h
#     $(CC) $(CFLAGS) -c $< -o $@ $(LDFLAGS) $(LDLIBS)
