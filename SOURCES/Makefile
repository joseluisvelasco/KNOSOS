##################################################################################
############################# KNOSOS Makefile ####################################
##################################################################################


###### Detect system #############################################################

HOSTNAME := $(shell uname -n|cut -f1 -d.)
NODENAME := $(shell uname -n|rev|cut -f2 -dl|rev)

ifeq ($(HOSTNAME), euler01)
CASE=mpi_euler
endif
ifeq ($(HOSTNAME), euler02)
CASE=mpi_euler
endif
ifeq ($(NODENAME), eu)
CASE=mpi_euler
endif
ifeq ($(HOSTNAME), dirac)
CASE=ifort_dirac
endif
ifeq ($(HOSTNAME), fusc3)
CASE=ifort_fusc3
endif
ifeq ($(HOSTNAME), macfus33)
CASE=gfortran_desktop
EXT=Darwin
endif
ifeq ($(HOSTNAME), MBdJLVG)
CASE=gfortran_desktop
EXT=Darwin
endif
ifeq ($(HOSTNAME), xUdJLVBG)
CASE=gfortran_desktop
EXT=Linux
endif
ifeq ($(HOSTNAME), task3d-l)
CASE=ifort_knl
endif
ifeq ($(HOSTNAME), draco01)
CASE=mpi_draco
endif
ifeq ($(HOSTNAME), draco02)
CASE=mpi_draco
endif


###### if running at draco #####################################

ifeq ($(CASE), mpi_draco)
COMPILER =mpiifort
#MPI_DIR  =
#PETSC_DIR=
#FFT_DIR  =
#MKL_DIR  =
#COMP_LIB =
FFLAGS   =-fpp -DDRACO -DMPIandPETSc -O3 -xHost -shared-intel -warn all -mcmodel large --with-debugging=0
IFLAGS   =-I${PETSC_DIR}/lib/petsc/conf/variables -I${PETSC_DIR}/lib/petsc/conf/rules -I${PETSC_DIR}/include -I${PETSC_DIR}/include/petsc -qopenmp -I${PETSC_DIR}/include
LFLAGS   =-L${PETSC_DIR}/lib -lpetsc -L/${MKL_HOME}/lib/intel64 -mkl -Wl,-rpath,${MKL_HOME}/lib/intel64,-rpath,${PETSC_DIR}/lib
endif

###### if running at TASK3D-L #####################################

ifeq ($(CASE), ifort_knl)
COMPILER =mpiifort
#MPI_DIR  =
PETSC_DIR=/usr/local/petsc-3.11.2
#FFT_DIR  =
#MKL_DIR  =
#COMP_LIB =
FFLAGS   =-fpp -DMPIandPETSc -DTASK3D -O3 -xHost -shared-intel -warn all -mcmodel large --with-debugging=0 
IFLAGS   =-I${PETSC_DIR}/include -I${PETSC_DIR}/include/petsc 
LFLAGS   =-mkl -lpetsc -L${PETSC_DIR}/lib 
endif


###### if running at euler cluster at CIEMAT #####################################

ifeq ($(CASE), mpi_euler)
COMPILER =mpif90 -traceback
MPI_DIR  =/opt/ofed_1.3.1/mpi/intel/openmpi-1.2.6/lib64/
PETSC_DIR=/home/localsoft/petsc-3.7.4/real
FFT_DIR  =/home/localsoft/intel-12/composerxe-2011.3.174/mkl/include/fft_w
MKL_DIR  =/home/roldan/gcc/520euler/lib64/
NETCDF_DIR=/soft/netcdf-4.0.1-intel/include
NETCDF_DIR2=/home/localsoft/netcdf-4.0.1-intel/lib
#COMP_LIB =
FFLAGS   =-cpp -DMPIandPETSc=.TRUE. -O3 -shared-intel -warn all -mcmodel large --with-debugging=0 
IFLAGS   =-I${PETSC_DIR}/include -I${PETSC_DIR}/include/petsc -I${FFT_DIR} -I${NETCDF_DIR} 
LFLAGS   =-mkl -L${MKL_DIR} -lpetsc -L${PETSC_DIR}/lib -L${MPI_DIR} -lnetcdf -L${NETCDF_DIR2} 
endif

###### if running at dirac computer at CIEMAT #####################################

ifeq ($(CASE), ifort_dirac)
COMPILER =ifort
#MPI_DIR  =
#PETSC_DIR=
#FFT_DIR  =
#MKL_DIR  =
#COMP_LIB =
FFLAGS   =-cpp -O3 -shared-intel -warn all
#IFLAGS   =
LFLAGS   =-mkl
endif

###### if running at fusc computer at CIEMAT #####################################

ifeq ($(CASE), ifort_fusc3)
COMPILER =ifort
#MPI_DIR  =
#PETSC_DIR=
#FFT_DIR  =
#MKL_DIR  =
COMP_LIB =LIBS/libnag_i4.a
FFLAGS   =-cpp -DNAG=.TRUE. -O3 -shared-intel -warn all
#IFLAGS   =
LFLAGS   =-mkl
endif

###### if running at Velasco's eMAC, Macbook or xUbuntu at CIEMAT ################

ifeq ($(CASE), gfortran_desktop)
COMPILER =gfortran
#MPI_DIR  =
#PETSC_DIR=
#FFT_DIR  =
#MKL_DIR  =
COMP_LIB =LIBS/libfftw3_${EXT}.a
FFLAGS   =-cpp --warn-all 
#IFLAGS   =
LFLAGS   =-llapack -lblas
endif

###### Targets, pre-requisites and commands ######################################

knosos.x: \
	global.o \
	amb_and_qn.o ana_ntv.o banana.o coefficients.o configuration.o dkes.o  \
	errors.o input.o knosos.o lagrange.o low_collisionality.o  mpi_and_end.o \
	plasmas.o plateau.o pfirsch_schlueter.o times.o trace_impurities.o transport.o

	${COMPILER} -o knosos.x global.o \
	amb_and_qn.o ana_ntv.o banana.o coefficients.o configuration.o dkes.o \
	errors.o input.o knosos.o lagrange.o low_collisionality.o mpi_and_end.o \
	plasmas.o plateau.o pfirsch_schlueter.o times.o trace_impurities.o transport.o \
	${COMP_LIB} ${LFLAGS} ${IFLAGS} ${FFLAGS} 

global.o: global.f90
	${COMPILER} ${FFLAGS} -c global.f90

amb_and_qn.o: amb_and_qn.f90
	${COMPILER} ${FFLAGS} -c amb_and_qn.f90

ana_ntv.o: ana_ntv.f90
	${COMPILER} ${FFLAGS} -c ana_ntv.f90

banana.o: banana.f90
	${COMPILER} ${FFLAGS} -c banana.f90

coefficients.o: coefficients.f90
	${COMPILER} ${FFLAGS} -c coefficients.f90

configuration.o: configuration.f90
	${COMPILER} ${FFLAGS} ${IFLAGS} -c configuration.f90

dkes.o: dkes.f90
	${COMPILER} ${FFLAGS} -c dkes.f90

errors.o: errors.f90
	${COMPILER} ${FFLAGS} -c errors.f90

input.o: input.f90
	${COMPILER} ${FFLAGS} -c input.f90

knosos.o: knosos.f90
	${COMPILER} ${FFLAGS} ${IFLAGS} -c knosos.f90

lagrange.o: lagrange.f90
	${COMPILER} ${FFLAGS} -c lagrange.f90

low_collisionality.o: low_collisionality.f90
	${COMPILER} ${FFLAGS} ${IFLAGS} -c low_collisionality.f90

mpi_and_end.o: mpi_and_end.f90
	${COMPILER} ${FFLAGS} -c mpi_and_end.f90

plasmas.o: plasmas.f90
	${COMPILER} ${FFLAGS} -c plasmas.f90

plateau.o: plateau.f90
	${COMPILER} ${FFLAGS} -c plateau.f90

pfirsch_schlueter.o: pfirsch_schlueter.f90
	${COMPILER} ${FFLAGS} ${IFLAGS} -c pfirsch_schlueter.f90

times.o: times.f90
	${COMPILER} ${FFLAGS} -c times.f90

trace_impurities.o: trace_impurities.f90
	${COMPILER} ${FFLAGS} -c trace_impurities.f90

transport.o: transport.f90
	${COMPILER} ${FFLAGS} -c transport.f90

clean:
	rm knosos.x *genmod* *.o global.mod
