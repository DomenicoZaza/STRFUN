
#MODEL
MODEL = -DSCAL 
#All possible models:    MODEL = -DSCAL -DFEEDBACK -DBUOYANCY

#COMPILER
COMP = mpif90 $(OPT)
#mpiifort 
#mpif90 $(OPT)
COMP_SERIAL = gfortran
#ifort
#gfortran

OPT = -g -fbacktrace  -fbounds-check -Waliasing -Wunderflow -Wsurprising -fbacktrace -fcheck=all -Wall -fcheck=all
#-O3 -funroll-loops -finline-functions -ftree-vectorize -fopt-info-vec-optimized -ffast-math -fwhole-file # -march=native 

#FLAG TO COMPILER (-cpp for pre processing of directives, INTEL -fpp)
FLAG =  -cpp  $(MODEL)

#LIBRARY PATHS
PATH2FFTW = $(FFTW_LIB)
#/cineca/prod/opt/libraries/fftw/3.3.9/intelmpi--oneapi-2021--binary/lib
#/opt/ohpc/pub/libs/gnu8/openmpi3/fftw/3.3.8/lib
PATH2OPENBLAS = $(OPENBLAS_LIB)
#/cineca/prod/opt/libraries/scalapack/2.1.0/intelmpi--oneapi-2021--binary/lib
#/opt/ohpc/pub/libs/gnu8/openblas/0.3.5/lib

#LINK LIBRARIES
LINK_LIB = -L$(PATH2FFTW) -lfftw3_mpi -lfftw3 -lm -L$(PATH2OPENBLAS) -lopenblas
#-L $(PATH2OPENBLAS) -lopenblas
#-lopenblas

#INCLUDE PATH
INC = -I$(FFTW_INCLUDE) -I$(OPENBLAS_INCLUDE)
#/cineca/prod/opt/libraries/fftw/3.3.9/intelmpi--oneapi-2021--binary/include -I /cineca/prod/opt/libraries/scalapack/2.1.0/intelmpi--oneapi-2021--binary/include 
#/opt/ohpc/pub/libs/gnu8/openmpi3/fftw/3.3.8/include



#OBJECT FILES USED BY POST
MODULES =  variables.o in_out.o transforms.o physics.o solvers.o  setup.o


#MAIN OBJECT
MAINOBJ = strfun.o

#EXECUTABLE
POST:		$(MODULES) $(MAINOBJ);	 $(COMP) $(FLAG) -o STRFUN.x $(MODULES) $(MAINOBJ) $(LINK_LIB) $(INC)


#COMPILE MODULES for POST
variables.o:		variables.f90;		$(COMP) $(FLAG) -c -o variables.o variables.f90
in_out.o:			in_out.f90;			$(COMP) $(FLAG) -c -o in_out.o in_out.f90
transforms.o:		transforms.f90;		$(COMP) $(FLAG) -c -o transforms.o transforms.f90 $(LINK_LIB) $(INC)
solvers.o:			solvers.f90;		$(COMP) $(FLAG) -c -o solvers.o solvers.f90
physics.o:			physics.f90;		$(COMP) $(FLAG) -c -o physics.o physics.f90
setup.o:			setup.f90;			$(COMP) $(FLAG) -c -o setup.o setup.f90

#POST-PROCESSING: 
strfun.o:		strfun.f90; 	$(COMP) $(FLAG) -c -o  strfun.o strfun.f90 $(LINK_LIB) $(INC)

# CLEAN
clean:
	rm -f *.o *.mod STRFUN.x

