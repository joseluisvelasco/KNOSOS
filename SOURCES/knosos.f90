!Main program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM KNOSOS

  USE GLOBAL
#ifdef MPIandPETSc
  USE petscsys
#endif
  IMPLICIT NONE
  !Parameters
  INTEGER, PARAMETER :: ntime=1
  REAL*8,  PARAMETER :: dt=1E-10
  !Allocatable
  INTEGER, ALLOCATABLE :: rank(:,:)
  REAL*8, ALLOCATABLE :: nb(:,:,:),dnbdpsi(:,:,:),Tb(:,:,:),dTbdpsi(:,:,:),Epsi(:,:)
  REAL*8, ALLOCATABLE :: Gb(:,:,:),Qb(:,:,:),Sb(:,:,:),Pb(:,:,:)
  !Other
  CHARACTER*100 serr
  INTEGER itime,is,ns,nbb,regb(nbx),jerr
  REAL*8 S(nsx),Zb(nbx),Ab(nbx),Zeff
#ifdef MPIandPETSc
  INTEGER ierr
! INCLUDE "mpif.h"
!#include <petsc/finclude/petscsys.h>
#endif

  !Initialize MPI
  CALL INITIALIZE_MPI()
  !Read input files (simulation parameters, models, flux-surfaces, species...)
  CALL READ_INPUT(ns,s,nbb,Zb,Ab,regb,Zeff)
  !Allocate some transport-related quantities
  ALLOCATE(nb(nbb,ns,nerr),dnbdpsi(nbb,ns,nerr),Tb(nbb,ns,nerr),dTbdpsi(nbb,ns,nerr),&
    & Epsi(ns,nerr),Gb(nbb,ns,nerr),Qb(nbb,ns,nerr),Sb(nbb,ns,nerr),Pb(nbb,ns,nerr),rank(ns,nerr))
  !Initialize PETSC and the random number generator according to the input, distribute MPI jobs
  CALL DISTRIBUTE_MPI(ns,rank)
#ifdef MPIandPETSc
  CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,myrank,myrank,PETSC_COMM_WORLD,ierr)
  CALL PETSCINITIALIZE(PETSC_NULL_CHARACTER,ierr)
#endif
  IF(nerr.GT.1) CALL INIT_RANDOM_SEED()
  CALL INIT_FILES()
  !Perform a loop in samples
  DO jerr=1,nerr
     !Follow time evolution  
     DO itime=1,ntime
        !Perform a loop in flux-surfaces
        DO is=1,ns
           !Êach process has been asigned some particular case
           IF(myrank.NE.rank(is,jerr)) CYCLE
           !Read the magnetic field
           IF(itime.EQ.1) THEN
              CALL READ_BFIELD(s(is))
              !Create (or read) a DKES-like database of monoenergetic transport coefficients
              CALL CALC_DATABASE(is,s(is))
              IF(ONLY_DB) CYCLE
              !Read the plasma profiles 
              CALL READ_PLASMAS(nbb,Zeff,s(is),Zb,Ab,nb(:,is,jerr),dnbdpsi(:,is,jerr),&
                   & Tb(:,is,jerr),dTbdpsi(:,is,jerr),Epsi(is,jerr))              
              !Read the particle and energy sources
              CALL READ_SOURCES(nbb,s(is),Sb(:,is,jerr),Pb(:,is,jerr))
           END IF
           !Calculate the flux for each surface
           !(incl. drift-kinetic eq., ambipolarity and quasineutrality)
           CALL SOLVE_DKE_QN_AMB(itime,nbb,Zb(1:nbb),Ab(1:nbb),regb(1:nbb),s(is),&
                & nb(:,is,jerr),dnbdpsi(:,is,jerr),Tb(:,is,jerr),dTbdpsi(:,is,jerr),&
                & Epsi(is,jerr),Gb(:,is,jerr),Qb(:,is,jerr))                
        END DO
        !Evolve the plasma profiles (without error propagation)
        CALL TRAnsPORT(nbb,ns,dt,S,Zb,Ab,nb,Gb(:,:,1),Sb,Tb,Qb(:,:,1),Pb,Epsi(:,1))
     END DO
  END DO
  IF(nerr.GT.1) CALL AVERAGE_SAMPLES(nbb,ns,s(1:ns),Epsi,Gb,Qb) 
  !End MPI and PETSC
  serr="Simulation complete"
  CALL END_ALL(serr,.TRUE.)

END PROGRAM KNOSOS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


