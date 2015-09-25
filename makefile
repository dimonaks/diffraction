

CPP=g++
F90=gfortran


all: diffraction
	cp ./diffraction ..
	# cd ..; ./diffraction list

diffraction: diffraction.o spline.o
	${CPP} diffraction.o spline.o   -o diffraction

diffraction.o: diffraction.cpp
	${CPP} -c diffraction.cpp

spline.o: fortran/spline.f
	${F90} -c -ffixed-line-length-none \
	fortran/spline.f

clean:
	rm *o hello

# cd src
# ${CPP} -c diffraction.cpp

# # cd fortran
# ${F90} -c -ffixed-line-length-none \
# fortran/spline.f

# # cd ..

 
# ${CPP} -o  diffraction diffraction.o spline.o




# cp diffraction ../

