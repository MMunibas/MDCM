# MDCM Makefile with serial, parallel (OpenMP) and GPU options.
#
# Serial: make serial
# OpenMP: make parallel
# GPU: make cuda
# All: make all
#
# Note that compilation with CUDA requires CUDA-Fortran to be installed and
# in your path

#FOR CUDA
F90_CUDA = nvfortran
#FFLAGS_PAR_CUDA = -fopenmp -D GPU -g #-Mbounds -Mchkptr -Mchkstk # nvfortran
#LFLAGS_PAR_CUDA = -fopenmp -D GPU -g -lcuda -Mcuda # nvfortran
FFLAGS_CUDA     = -g -cpp -D GPU #-Mbounds -Mchkptr -Mchkstk # nvfortran
#LFLAGS_CUDA     = -lcuda -D GPU -Mcuda
LFLAGS_CUDA     = -g -D GPU -lcuda -Mcuda
NVFLAGS = -D GPU -g -O0 #-Mbounds -Mchkptr -Mchkstk # nvfortran

# FOR CPU
F90 = gfortran
FFLAGS_PAR = -fopenmp -cpp -fcheck=bounds,mem -g #-pg
LFLAGS_PAR = -fopenmp -g #-pg
#FFLAGS     = -O3 -g 
FFLAGS     = -fcheck=bounds,mem -g -cpp #-pg
LFLAGS     = -g #-pg

.SUFFIXES: .f90 .f .o .cu .cuf

##########################
# Object Files for build #
##########################

OBJS = \
differential_evolution.o \
symmetry.o \
main.o \

CUDA_OBJS = \
cuda_differential_evolution.o \
cuda_symmetry.o \
rmse_gpu.o \
cuda_main.o \

OBJS_PAR = \
pdifferential_evolution.o \
psymmetry.o \
pmain.o \

# targets
parallel : pcubefit.x clean

serial : cubefit.x clean
cuda : cudacubefit.x clean
all : cubefit.x pcubefit.x cudacubefit.x clean

# rules
cubefit.x : $(OBJS) 
	${F90}  -o $@ $(LFLAGS) $(OBJS) 
	 
pcubefit.x : $(OBJS_PAR) 
	${F90}  -o $@ $(LFLAGS_PAR) $(OBJS_PAR) 

cudacubefit.x : $(CUDA_OBJS)
	${F90_CUDA}  -o $@ $(LFLAGS_CUDA) $(CUDA_OBJS)

#######################################
# Object dependencies and compilation #
#######################################
pdifferential_evolution.o : src/differential_evolution.F90
	$(F90) -c $(FFLAGS_PAR) $(INCLUDES) -o $@ src/differential_evolution.F90

psymmetry.o : src/symmetry.f90
	$(F90) -c $(FFLAGS_PAR) $(INCLUDES) -o $@ src/symmetry.f90

pmain.o : src/main.f90 \
pdifferential_evolution.o \
psymmetry.o 
	$(F90) -c $(FFLAGS_PAR) $(INCLUDES) -o $@ src/main.f90

differential_evolution.o : src/differential_evolution.F90
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ src/differential_evolution.F90
	
symmetry.o : src/symmetry.f90
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ src/symmetry.f90

main.o : src/main.f90 \
differential_evolution.o \
symmetry.o 
	$(F90) -c $(FFLAGS) $(INCLUDES) -o $@ src/main.f90

# CUDA binary
cuda_differential_evolution.o : src/differential_evolution.F90
	$(F90_CUDA) -c $(FFLAGS_CUDA) $(INCLUDES) -o $@ src/differential_evolution.F90
	
cuda_symmetry.o : src/symmetry.f90
	$(F90_CUDA) -c $(FFLAGS_CUDA) $(INCLUDES) -o $@ src/symmetry.f90

rmse_gpu.o : src/rmse_gpu.cuf
	$(F90_CUDA) -c $(NVFLAGS) $(INCLUDES) -o $@ src/rmse_gpu.cuf

cuda_main.o : src/main.f90 \
cuda_differential_evolution.o \
cuda_symmetry.o \
rmse_gpu.o
	$(F90_CUDA) -c $(FFLAGS_CUDA) $(INCLUDES) -o $@ src/main.f90

.PHONY: clean veryclean
clean:
	rm *.o *.mod

veryclean:
	rm *.x *.o *.mod



