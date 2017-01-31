NAME = neighbors

.SUFFIXES: .f90

.f90.o:
	gfortran -O3 -ffast-math -finline-functions -finline-small-functions -fPIC -c $< -o $@
#	gfortran -O3 -pg -fbounds-check -fPIC -c $< -o $@

all: ${NAME}.o wrap clean

wrap:
	f90wrap -m ${NAME}_module ${NAME}.f90 -v -k kind_map
	f2py-f90wrap -c -m _${NAME}_module f90wrap*.f90 *.o

# fortran:
# #	gfortran -c globals.f90 ${NAME}.f90
# 	gfortran -c ${NAME}.f90

clean:
	-rm *.o *.mod #f90wrap_*

veryclean:
	-rm *.o *.mod *.so f90wrap_*
