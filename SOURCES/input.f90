!Read input files

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE READ_INPUT(ns,s,nbb,Zb,Ab,regb,Zeff)

!-------------------------------------------------------------------------------------------------
!Read simulation input. This includes the number ns of surfaces s, the number of species nnb and
!their charge number Zb, mass number Ab and NC regime regb, and the effective charge Zeff.
!The rest are global variables.
!All of them are described in the user manual.  
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Namelist 'parameters' contains simulation parameters, including
  !decisions on how to solve some equations and what to plot
  NAMELIST /parameters/  GEN_FLAG,DEBUG,TIME,I0,PLOTG,&
       & USE_B1,USE_B0pB1,&
       & REMOVE_DIV,DELTA,&
       & MLAMBDA,MAL,TRUNCATE_B,PREC_B,PREC_EXTR,PREC_BINT,PREC_DQDV,PREC_INTV,&
       & NEFIELD,EFIELD,NCMUL,CMUL,NVMAG,VMAG,&
       & NER,ERMIN,ERMAX,ERACC,&
       & NERR,IPERR
  !Namelist 'model' contains physics parameters
  NAMELIST /model/ CALC_DB,ONLY_DB,INC_EXB,TANG_VM,CLASSICAL,ANISOTROPY,FRICTION,FACT_CON,&
       & SOLVE_AMB,TRIVIAL_AMB,FAST_AMB,SOLVE_QN,TRIVIAL_QN,ZERO_PHI1,ONLY_PHI1,D_AND_V,COMPARE_MODELS,&
       & FN,FI,FS,FP,FE,FR,FB,FNE,FTI,FTE,FER,&
       & NEOTRANSP,PENTA,TASK3D,SATAKE,ANA_NTV,JPP
  !Namelist 'surfaces' contains the list of flux-surfaces calculated
  NAMELIST /surfaces/ NS,S,DIRDB,DIRS
  !Namelist 'species' contains the list of species calculated
  NAMELIST /species/ NBB,ZB,AB,REGB,ZEFF
  !Namelist 'others' contains other variables
  NAMELIST /others/ QN,TRACE_IMP,PLATEAU_OR_PS,USE_B0,NTV,&
       & IMP1NU,CENTERED_ALPHA,CLOSEST_LAMBDA,CLOSEST_LAMBDA1
  !Output
  INTEGER NS,ib,NBB,REGB(nbx)
  REAL*8 S(nsx),ZB(nbx),AB(nbx),ZEFF
  !Other
  CHARACTER*100 serr,line
  INTEGER iline,iostat,nefield,ncmul,nvmag,iz,ia,is
  !DKES database
  INTEGER, PARAMETER :: ncmuld=18
  REAL*8 cmuld(ncmuld) /3E+2,1E+2,3E+1,1E+1,3E+0,1E+0,3E-1,1E-1,3E-2,1E-2,&
                      & 3E-3,1E-3,3E-4,1E-4,3E-5,1E-5,3E-6,1E-6/
  INTEGER, PARAMETER :: nefieldd=9
  REAL*8 efieldd(nefieldd) /0E-0,1E-5,3E-5,1E-4,3E-4,1E-3,3E-3,1E-2,3E-2/
  INTEGER, PARAMETER :: nvmagd=1
  REAL*8 vmagd(nvmagd) /0E+0/
  !3D database
  INTEGER, PARAMETER :: ncmulx=18
  REAL*8 cmulx(ncmulx) /3E+2,1E+2,3E+1,1E+1,3E+0,1E+0,3E-1,1E-1,3E-2,1E-2,&
                      & 3E-3,1E-3,3E-4,1E-4,3E-5,1E-5,3E-6,1E-6/
!  INTEGER, PARAMETER :: nefieldx=9
!  REAL*8 efieldx(nefieldx) /0E-0,1E-5,3E-5,1E-4,3E-4,1E-3,3E-3,1E-2,3E-2/
!  INTEGER, PARAMETER :: nvmagx=13
!  REAL*8 vmagx(nvmagx) /-3E-2,-1E-2,-3E-3,-1E-3,-3E-4,-1E-4,0E-0,&
!                      & +1E-4,+3E-4,+1E-3,+3E-3,+1E-2,+3E-2/
!  3D finer database 
!  INTEGER, PARAMETER :: ncmulx=52
!  REAL*8 cmulx(ncmulx) /5E+2,2E+2,1E+2,5E+1,2E+1,1E+1,5E+0,2E+0,1.0E+0,&
!                      & 5E-1,2E-1,1E-1,5E-2,2E-2,1E-2,5E-3,2E-3,1.5E-3,1.0E-3,&
!                      & 9E-4,8E-4,7E-4,6E-4,5E-4,4E-4,3E-4,2E-4,1.5E-4,1.2E-4,1E-4,&
!                      & 9E-5,8E-5,7E-5,6E-5,5E-5,4E-5,3E-5,2E-5,1.5E-5,1.2E-5,1E-5,&
!                      & 9E-6,8E-6,7E-6,6E-6,5E-6,4E-6,3E-6,2E-6,1.5E-6,1.2E-6,1E-6/
  INTEGER, PARAMETER :: nefieldx=50 
  REAL*8 efieldx(nefieldx) /0E-0,1E-5,1.2E-5,1.5E-5,2E-5,3E-5,4E-5,5E-5,6E-5,7E-5,8E-5,9E-5,&
                              &  1E-4,1.2E-4,1.5E-4,2E-4,3E-4,4E-4,5E-4,6E-4,7E-4,8E-4,9E-4,&
                              &  1E-3,1.2E-3,1.5E-3,2E-3,3E-3,4E-3,5E-3,6E-3,7E-3,8E-3,9E-3,&
                              &  1E-2,1.2E-2,1.5E-2,2E-2,3E-2,4E-2,5E-2,6E-2,7E-2,8E-2,9E-2,&
                              &  1E-1,1.2E-1,1.5E-1,2E-1,5E-5/
  INTEGER, PARAMETER :: nvmagx=66
  REAL*8 vmagx(nvmagx)/-9E-2,-8E-2,-7E-2,-8E-2,-5E-2,-4E-2,-3E-2,-2E-2,-1.5E-2,-1.2E-2,-1E-2,&
                     & -9E-3,-8E-3,-7E-3,-8E-3,-5E-3,-4E-3,-3E-3,-2E-3,-1.5E-3,-1.2E-3,-1E-3,&
                     & -9E-4,-8E-4,-7E-4,-8E-4,-5E-4,-4E-4,-3E-4,-2E-4,-1.5E-4,-1.2E-4,-1E-4,&
                     & +1E-4,+1.2E-4,+1.5E-4,+2E-4,+3E-4,+4E-4,+5E-4,+6E-4,+7E-4,+8E-4,+9E-4,&
                     & +1E-3,+1.2E-3,+1.5E-3,+2E-3,+3E-3,+4E-3,+5E-3,+6E-3,+7E-3,+8E-3,+9E-3,&
                     & +1E-2,+1.2E-2,+1.5E-2,+2E-2,+3E-2,+4E-2,+5E-2,+6E-2,+7E-2,+8E-2,+9E-2/
  !Database for neontrasp
  INTEGER, PARAMETER :: ncmuly=16
  REAL*8 cmuly(ncmuly) /1E+1,3E+0,1E+0,3E-1,1E-1,3E-2,1E-2,&
                      & 3E-3,1E-3,3E-4,1E-4,3E-5,1E-5,3E-6,1E-6,3E-7/
  INTEGER, PARAMETER :: nefieldy=27
  REAL*8 efieldy(nefieldy) /0.0,3E-7,1E-6,3E-6,1E-5,3E-5,1E-4,3E-4,1E-3,2E-3,5E-3,&
       & 1E-2,2E-2,3E-2,5E-2,1E-1,2E-1,3E-1,5E-1,7E-1,8E-1,1.0,1.2,1.5,2.0,3.0,5.0/
  INTEGER, PARAMETER :: nvmagy=1
  REAL*8 vmagy(nvmagy) /0.0/
!  INTEGER, PARAMETER :: nvmagy=66
!  REAL*8 vmagy(nvmagy)/-9E-2,-8E-2,-7E-2,-8E-2,-5E-2,-4E-2,-3E-2,-2E-2,-1.5E-2,-1.2E-2,-1E-2,&
!                     & -9E-3,-8E-3,-7E-3,-8E-3,-5E-3,-4E-3,-3E-3,-2E-3,-1.5E-3,-1.2E-3,-1E-3,&
!                     & -9E-4,-8E-4,-7E-4,-8E-4,-5E-4,-4E-4,-3E-4,-2E-4,-1.5E-4,-1.2E-4,-1E-4,&
!                     & +1E-4,+1.2E-4,+1.5E-4,+2E-4,+3E-4,+4E-4,+5E-4,+6E-4,+7E-4,+8E-4,+9E-4,&
!                     & +1E-3,+1.2E-3,+1.5E-3,+2E-3,+3E-3,+4E-3,+5E-3,+6E-3,+7E-3,+8E-3,+9E-3,&
!                     & +1E-2,+1.2E-2,+1.5E-2,+2E-2,+3E-2,+4E-2,+5E-2,+6E-2,+7E-2,+8E-2,+9E-2/
  REAL*8 dummy,cmul(MAX(ncmulx,ncmuly)),efield(MAX(nefieldx,nefieldy)),vmag(MAX(nvmagx,nvmagy))


  OPEN(unit=1000,file="STDOUT",form='formatted',action='write',iostat=iostat)

  !-------------------------------------------------------------------------------------------
  !Set values for namelist 'parameters'
  !-------------------------------------------------------------------------------------------

  !To be discontinued
  IMP1NU        =.FALSE.
  CENTERED_ALPHA=.FALSE.
  FIRST_ORDER_ALPHA=.FALSE.
  CLOSEST_LAMBDA=.TRUE.
  CLOSEST_LAMBDA1=.TRUE.

  !Default values
  !Debug
  GEN_FLAG=.FALSE. 
  DEBUG=   .FALSE.  
  TIME=    .TRUE.   
  I0=      -1        
  PLOTG=   .FALSE.  
  !How to describe magnetic field structure
  USE_B1=   .FALSE.  
  USE_B0pB1=.FALSE.  
  !Determine algorithms to be used
  REMOVE_DIV=    .FALSE.  
  DELTA=         .TRUE.   
  !Set numerical resolution
  MLAMBDA=128
  MAL    = 64
  TRUNCATE_B= -200     
  PREC_B=     1E-5     
  PREC_EXTR=  1E-7     
  PREC_BINT=  5E-2     
  PREC_DQDV=  5E-2     
  PREC_INTV=  1E-2     
  !Set details of the monoenergetic database
  NCMUL  =ncmuld       
  NEFIELD=nefieldd     
  NVMAG  =nvmagd
  CMUL(1:ncmuld)    =cmuld       
  EFIELD(1:nefieldd)=efieldd     
  VMAG(1:nvmagd)    =vmagd        
  !Set details of ambipolarity
  NER  = 21     
  ERMIN=+20.0  
  ERMAX=-20.0  
  ERACC= 0.1    
  !Estimate error bars
  NERR =1       
  IPERR=-1.0   

  !Read namelist 'parameters'
  OPEN(unit=1,file="input.parameters",form='formatted',action='read',iostat=iostat) 
  IF(iostat.NE.0) OPEN(unit=1,file="../input.parameters",form='formatted',action='read',iostat=iostat) 
  IF(iostat.EQ.0) THEN
     IF(myrank.EQ.0) WRITE(1000,*) 'File "input.parameters" found'
     READ (1,nml=parameters)
     CLOSE(1)
  END IF

  !-------------------------------------------------------------------------------------------
  !Set values for namelist 'model'
  !-------------------------------------------------------------------------------------------

  !Default values
  !Ambipolarity and quasineutrality
  FAST_AMB=   .FALSE.  
  SOLVE_AMB=  .TRUE.  
  TRIVIAL_AMB=.FALSE.
  SOLVE_QN=   .TRUE.  
  TRIVIAL_QN= .FALSE. 
  ZERO_PHI1 = .FALSE.
  ONLY_PHI1 = .FALSE.
  D_AND_V   =5
  COMPARE_MODELS=.FALSE.
  !Details of the drift-kinetic equation
  INC_EXB=   .FALSE. 
  TANG_VM=   .TRUE.  
  CLASSICAL= .TRUE.  
  ANISOTROPY=.TRUE. 
  FRICTION=  .TRUE. 
  FACT_CON=  -3      
  !Equations to be solved
  DIRDB  ="./"          
  DIRS   ="./"           
  CALC_DB=.FALSE.   
  ONLY_DB=.FALSE.   
  !Scan in parameters
  FN=1.0
  FI=1.0
  FS=1.0
  FP=1.0
  FE=1.0
  FR=1.0
  FB=1.0
  FNE=1.0
  FTE=1.0
  FTI=1.0
  FDLNE=1.0
  FDLTE=1.0
  FDLTI=1.0
  FER=1.0
  !Particular problems
  NEOTRANSP=.FALSE.
  TASK3D=   .FALSE. 
  PENTA=    .FALSE. 
  SATAKE=   .FALSE. 
  ANA_NTV=  .FALSE. 
  JPP=      .FALSE.
  
  !Read namelist 'model'
  OPEN(unit=1,file="input.model",form='formatted',action='read',iostat=iostat) 
  IF(iostat.NE.0) OPEN(unit=1,file="../input.model",form='formatted',action='read',iostat=iostat) 
  IF (iostat.EQ.0) THEN
     IF(myrank.EQ.0) WRITE(1000,*) 'File "input.model" found'
     READ (1,nml=model)
     CLOSE(1)
  END IF

  !-------------------------------------------------------------------------------------------
  !Set default values for namelists 'surfaces' and 'species'
  !-------------------------------------------------------------------------------------------

  ns=1
  s(1)=1
  !Default values: low collisionality hydrogen + adiabatic electrons, no impurities
  nbb    = 2
  Zb(1)  =-1.            
  Ab(1)  = 5.4858E-4
  regb(1)=-2           
  Zb(2)  =+1.            
  Ab(2)  = 1.0073
  regb(2)=+3
  IF(TASK3D) THEN
     !Could be overwritten by file 'input.surfaces'
     ns=40
     DO is=1,ns
        s(is)=(is-0.5)*(is-0.5)/REAL(ns*ns)
     END DO
     !Could be overwritten by file 'input.species'
     NBB=6
     REGB(1:2)=3
     ZB(1)=-1.
     AB(1)=5.48579909E-4
     OPEN(unit=1,file="input-prof.txt",action='read',iostat=iostat)
     IF(iostat.NE.0) OPEN(unit=1,file="../input-prof.txt",action='read',iostat=iostat)     
     IF(iostat.EQ.0) THEN
        IF(myrank.EQ.0) WRITE(1000,*) 'File "input-prof.txt" found'
        DO iline=1,7
           READ(1,*) line
        END DO
        READ(1,*) dummy,dummy,dummy,iz,ia
        CLOSE(1)
        Zb(2)=REAL(iz)
        Ab(2)=ia*1.00727647
     ELSE
        serr="Wrong data in input-prof.txt"
        CALL END_ALL(serr,.FALSE.)
     END IF
     Zb(3)  =+3.         
     Ab(3)  =6.941
     Zb(4)  =+19.         
     Ab(4)  =47.867
     Zb(5)  =+20.            
     Ab(5)  =50.9415
     Zb(6)  =+23.            
     Ab(6)  =55.847
     regb(3:nbb)=10
     ZEFF   =1.0000000001
     !Could be overwritten by 'input.model'
     SOLVE_QN  =.FALSE.
     TRIVIAL_QN=.TRUE.
     D_AND_V=5
  ELSE IF(NEOTRANSP) THEN
     ns=7
     DO is=1,ns
        s(is)=(is-0.5)*(is-0.5)/REAL(ns*ns)
     END DO
  END IF


  !-------------------------------------------------------------------------------------------
  !Set values for namelist 'surfaces'
  !-------------------------------------------------------------------------------------------

  !Read namelist 'surfaces'
  OPEN(unit=1,file="input.surfaces",form='formatted',action='read',iostat=iostat) 
  IF(iostat.NE.0) OPEN(unit=1,file="../input.surfaces",form='formatted',action='read',iostat=iostat) 
  IF (iostat.EQ.0) THEN
     IF(myrank.EQ.0) WRITE(1000,*) 'File "input.surfaces" found'
     READ (1,nml=surfaces)
     CLOSE(1)
  END IF
  
  !-------------------------------------------------------------------------------------------
  !Set values for namelist 'species'
  !-------------------------------------------------------------------------------------------

  !Read namelist 'species'
  OPEN(unit=1,file="input.species",form='formatted',action='read',iostat=iostat) 
  IF(iostat.NE.0) OPEN(unit=1,file="../input.species",form='formatted',action='read',iostat=iostat) 
  IF (iostat.EQ.0) THEN
     IF(myrank.EQ.0) WRITE(1000,*) 'File "input.species" found'
     READ (1,nml=species)
     CLOSE(1)
  END IF
  
  !-------------------------------------------------------------------------------------------

  !These variables are determined at some point of the run
  CALCULATED_INT=.FALSE.
  PHI1_READ=     .FALSE.
  DKES_READ=     .FALSE.
  CONVERGED=     .FALSE.
  TRACE_IMP=     .FALSE.
  PLATEAU_OR_PS= .FALSE.
  
  !Some of the choices above may be incompatible, have to be changed accordingly.
  DO ib=1,nbb !in case you don't remember the electron mass, just write a small number
     IF(Ab(ib).LT.1) Ab(ib)=5.48579909E-4 
     IF((ib.EQ.1.AND.Zb(ib).GT.0).OR.(ib.GT.1.AND.Zb(ib).LT.0)) THEN
        serr="Wrong data in input.species"
        CALL END_ALL(serr,.FALSE.)
     END IF
  END DO

  IF(FAST_AMB) THEN
     CALC_DB=.TRUE.
     NCMUL  =ncmulx  
     NEFIELD=nefieldx     
     NVMAG  =nvmagx
     CMUL(1:ncmulx)    =cmulx       
     EFIELD(1:nefieldx)=efieldx     
     VMAG(1:nvmagx)    =vmagx        
  END IF
  IF(NEOTRANSP.OR.PENTA) THEN
     CALC_DB=.TRUE.
     ONLY_DB=.TRUE.
     TANG_VM=.FALSE.
     INC_EXB=.TRUE.
     SOLVE_AMB=.FALSE.
     SOLVE_QN=.FALSE.       
     IF(NEOTRANSP) THEN
        ncmul  =ncmuly
        nefield=nefieldy
        nvmag  =nvmagy
        cmul(1:ncmuly)    =cmuly
        efield(1:nefieldy)=efieldy
        vmag(1:nvmagy)    =vmagy
     END IF
  END IF
  IF(.NOT.TANG_VM) THEN
     FS=0
     nvmag=1
     vmag=0.0
  END IF
  USE_B0=USE_B1.OR.USE_B0pB1
  DO ib=3,nbb
     IF(REGB(ib).GE.10) THEN
        TRACE_IMP    =.TRUE.
        IF(REGB(ib).LT.13) PLATEAU_OR_PS=.TRUE.
     END IF
  END DO
  IF(TRACE_IMP) THEN
     IF(REGB(2).LE.-1) REGB(2)=+3
     IF(.NOT.SOLVE_QN) TRIVIAL_QN=.TRUE.
  ELSE
     ANISOTROPY=.FALSE.
  END IF
  IF(TRIVIAL_QN)  SOLVE_QN =.FALSE.
  IF(TRIVIAL_AMB) SOLVE_AMB=.FALSE.
  QN=SOLVE_QN.OR.TRIVIAL_QN
  IF(SOLVE_QN) THEN
     DELTA=.TRUE.
     CLOSEST_LAMBDA=.FALSE.
  END IF
  IF(.NOT.SOLVE_QN) COMPARE_MODELS=.FALSE.
  IF(SATAKE) TRUNCATE_B=20
  NTV=SATAKE
  IF(PLOT_XYZ) THEN
     PREC_B=-1E-9
     TRUNCATE_B=-10
  END IF
  IF(MAL.GT.0) ALLOCATE(absnablar(MAL,MAL))

  IF(GEN_FLAG(1)) CLOSEST_LAMBDA=.NOT.CLOSEST_LAMBDA

  ncmult  =ncmul
  nefieldt=nefield
  nvmagt  =nvmag
  ALLOCATE( cmult(ncmult), efieldt(nefieldt), vmagt(nvmagt))
  ALLOCATE(lcmult(ncmult),lefieldt(nefieldt),lvmagt(nvmagt))
  ALLOCATE(lD11dkes1(ncmult,nefieldt),lD11dkes2(ncmult,nefieldt),lD11tab(ncmult,nefieldt,nvmagt))
  cmult  =cmul
  efieldt=efield
  vmagt  =vmag
  lcmult  =LOG(cmult)   
  lefieldt=LOG(efieldt)
  lvmagt  =LOG(ABS(lvmagt))
  lefieldt(1)=-1000 !dummy value to avoid log(0.0)


  
  !Write input parameters
  IF(myrank.EQ.0) THEN
     WRITE(1000,nml=model)
     WRITE(1000,nml=parameters)
     WRITE(1000,nml=surfaces)
     WRITE(1000,nml=species)
     WRITE(1000,nml=others)
     CLOSE(1000)
  END IF


END SUBROUTINE READ_INPUT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE INIT_FILES()

!-------------------------------------------------------------------------------------------------
!Initialize output files
!-------------------------------------------------------------------------------------------
  USE GLOBAL
  IMPLICIT NONE
  !Others
  CHARACTER*100 filename
  INTEGER iostat
  
  IF(DEBUG) THEN
     OPEN(unit=1300+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=1400+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=1500+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=2000+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=2100+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=2900+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=3000+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=3100+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=3200+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=3300+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=3400+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=4400+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=5000+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=5100+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=5200+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=5300+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=5400+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=5500+myrank,form='formatted',action='write',iostat=iostat)
     OPEN(unit=5600+myrank,form='formatted',action='write',iostat=iostat)
  END IF

  IF(numprocs.EQ.1) THEN
     filename="STDOUT"
     OPEN(unit=1000+myrank,file=filename,form='formatted',action='write',iostat=iostat,&
          access='append',status='old')
     filename="STDERR"
     OPEN(unit=1100+myrank,file=filename,form='formatted',action='write',iostat=iostat)
  ELSE
     WRITE(filename,'("STDOUT.",I2.2)') myrank
     OPEN(unit=1000+myrank,file=filename,form='formatted',action='write',iostat=iostat)
     WRITE(filename,'("STDERR.",I2.2)') myrank
     OPEN(unit=1100+myrank,file=filename,form='formatted',action='write',iostat=iostat)
  END IF
  
  IF(numprocs.EQ.1) filename="B.map"
  IF(numprocs.GT.1) WRITE(filename,'("B.map.",I2.2)') myrank
  OPEN(unit=100+myrank,file=filename,form='formatted',action='write',iostat=iostat)
  WRITE(100+myrank,&
       & '("s \zeta_{Boozer}  \theta_{Boozer}(right-handed)  B[T]  (v_B.\nabla\psi)[A.U.]")')

!  IF(numprocs.EQ.1) filename="results.knosos"
!  IF(numprocs.GT.1) WRITE(filename,'("results.knosos.",I2.2)') myrank
!  OPEN(unit=200+myrank,file=filename,form='formatted',action='write',iostat=iostat)
!  WRITE(200+myrank,'("cmul efield weov wtov L11m L11p L31m L31p L33m L33p scal11&
!       & scal13 scal33 max\_residual chip psip btheta bzeta vp vmag")')

  IF(SOLVE_AMB.OR.SOLVE_QN) THEN
     IF(numprocs.EQ.1) filename="flux.amb"
     IF(numprocs.GT.1) WRITE(filename,'("flux.amb.",I2.2)') myrank
     OPEN(unit=300+myrank,file=filename,form='formatted',action='write',iostat=iostat)
     WRITE(300+myrank,'("s E_r (\Gamma_b/n_b[m/s]  Q_b/n_b/T_b[m/s] L_1^b[m^2/s] L_2^b[m^2/s]&
          &  Z_b n_b[10^{19}m^{-3}] T_b[eV] dlnT_b/dr, b=1,NBB), size(e\varphi_1/T_i) ")')

     IF(TASK3D) THEN
        IF(numprocs.EQ.1) filename="knososTASK3D.flux"
        IF(numprocs.GT.1) WRITE(filename,'("knososTASK3D.flux.",I2.2)') myrank
        OPEN(unit=5300+myrank,file=filename,form='formatted',action='write',iostat=iostat)
        WRITE(5300+myrank,'("rho   E_r[kV/m]   Gamma_e[m^-2s^-1]   Gamma_i1[m^-2s^-1]   Gamma_i2[m^-2s^-1]&
             &  Gamma_i[m^-2s^-1]   Q_e[W/m^2]   Q_i1[W/m^2]   Q_i2[W/m^2]   Q_i[W/m^2]")')  
     END IF

     IF(COMPARE_MODELS) THEN
        IF(numprocs.EQ.1) filename="flux.amb.comp"
        IF(numprocs.GT.1) WRITE(filename,'("flux.amb.comp.",I2.2)') myrank
        OPEN(unit=4300+myrank,file=filename,form='formatted',action='write',iostat=iostat)
        WRITE(4300+myrank,'("s E_r[V/m] (\Gamma_b/n_b[m/s]  Q_b/n_b/T_b[m/s] L_1^b[m^2/s] L_2^b[m^2/s]&
          &  n_b[10^{19}m^{-3}] dlnn_b/dr T_b[eV] dlnT_b/dr Z_b, b=1,NBB), size(e\varphi_1/T_i) ")')
     END IF
     
     IF(SOLVE_QN) THEN
        IF(numprocs.EQ.1) filename="flux.modes"
        IF(numprocs.GT.1) WRITE(filename,'("flux.modes.",I2.2)') myrank
        OPEN(unit=700+myrank,file=filename,form='formatted',action='write',iostat=iostat)
        WRITE(700+myrank,'("s E_r[V/m] cosine(0)/sine(1) n  m &
	& (\Gamma_b/n_b[m/s]  Q_b/n_b/T_b[m/s] L_1^b[m^2/s] L_2^b[m^2/s], b=1,NBB)")')
     END IF
  END IF
  
  IF(QN) THEN
     IF(numprocs.EQ.1) filename="varphi1.modes"
     IF(numprocs.GT.1) WRITE(filename,'("varphi1.modes.",I2.2)') myrank
     OPEN(unit=400+myrank,file=filename,form='formatted',action='write',iostat=iostat)
     WRITE(400+myrank,&
     &'("s  cosine(0)/sine(1) n  m  \varphi_{nm} (angles are \zeta_{Boozer},\theta_{Boozer} right-handed)")')

     IF(numprocs.EQ.1) filename="varphi1.map"
     IF(numprocs.GT.1) WRITE(filename,'("varphi1.map.",I2.2)') myrank
     OPEN(unit=500+myrank,file=filename,form='formatted',action='write',iostat=iostat)
     WRITE(500+myrank,'("s  \zeta_{Boozer}  \theta_{Boozer}(right-handed) &
          & \varphi_1[V]  e\varphi_1/T_i  (v_E.\nabla\psi)[A.U.]  n1e/n0e  n1i/n0i ")')

     IF(COMPARE_MODELS) THEN

        IF(numprocs.EQ.1) filename="varphi1.modes.comp"
        IF(numprocs.GT.1) WRITE(filename,'("varphi1.modes.comp.",I2.2)') myrank
        OPEN(unit=4100+myrank,file=filename,form='formatted',action='write',iostat=iostat)
        WRITE(4100+myrank,&
        &'("s cosine(0)/sine(1) n  m  \varphi_{nm} (angles are \zeta_{Boozer},\theta_{Boozer} right-handed)")')
        
        IF(numprocs.EQ.1) filename="varphi1.map.comp"
        IF(numprocs.GT.1) WRITE(filename,'("varphi1.map.comp.",I2.2)') myrank
        OPEN(unit=4200+myrank,file=filename,form='formatted',action='write',iostat=iostat)
        WRITE(4200+myrank,'("s \zeta_{Boozer}  \theta_{Boozer}(right-handed) &
             & \varphi_1[V]  e\varphi_1/T_i  (v_E.\nabla\psi)[A.U.]  n1e/n0e  n1i/n0i ")')
        
     END IF

   END IF

   IF(numprocs.EQ.1) filename="flux.knosos"
   IF(numprocs.GT.1) WRITE(filename,'("flux.knosos.",I2.2)') myrank
   OPEN(unit=600+myrank,file=filename,form='formatted',action='write',iostat=iostat)
   WRITE(600+myrank,'("s E_r[V/m] (\Gamma_b/n_b[m/s]  Qb/n_b/T_b[m/s] L_1^b[m^2/s] L_2^b[m^2/s]&
        &  n_b[10^{19}m^{-3}] dlnn_b/dr T_b[eV] dlnT_b/dr Z_b, b=1,NBB), size(e\varphi_1/T_i) ")')
   
   IF(TASK3D) THEN
      IF(numprocs.EQ.1) filename="knososTASK3D.ambEr"
      IF(numprocs.GT.1) WRITE(filename,'("knososTASK3D.ambEr.",I2.2)') myrank
      OPEN(unit=5600+myrank,file=filename,form='formatted',action='write',iostat=iostat)
      WRITE(5600+myrank,'("rho   n_e[m^-3]   n_i1[m^-3   n_i2[m^-3]   T_e[eV]   T_i1[eV]   T_i2[eV]  &
           & E_r[kV/m]   Gamma[m^-2s^-1]   Q_e[W/m^2]   Q_i[W/m^2]   Chi_e[m^2/s]   Chi_i[m^2/s]")')
   END IF
   
   IF(nerr.GT.1) THEN
      IF(numprocs.EQ.1) THEN
         filename="flux.av"
         OPEN(unit=800+myrank,file=filename,form='formatted',action='write',iostat=iostat)
         WRITE(800+myrank,'("ohos E_r[V/m] err(E_r)[V/m] (\Gamma_b/n_b[m/s] err(\Ganna_b/n_b)[m/s] &
              & Qb/n_b/T_b[m/s], err(Qb/n_b/T_b)[m/s], b=1,NBB) ")')
      ELSE
         WRITE(filename,'("flux.av.",I2.2)') myrank
         OPEN(unit=800+myrank,file=filename,form='formatted',action='write',iostat=iostat)
         WRITE(800+myrank,'("s E_r[V/m] err(E_r)[V/m] (\Gamma_b/n_b[m/s] err(\Ganna_b/n_b)[m/s] &
              & Qb/n_b/T_b[m/s], err(Qb/n_b/T_b)[m/s], b=1,NBB) ")')
      END IF
   END IF
   
   IF(TRACE_IMP) THEN
      IF(numprocs.EQ.1) filename="imp.knosos"
      IF(numprocs.GT.1) WRITE(filename,'("imp.knosos.",I2.2)') myrank
      OPEN(unit=900+myrank,file=filename,form='formatted',action='write',iostat=iostat)
      WRITE(900+myrank,'("s  Z_z  A_z \Gamma_z/n_z[m/s]  V_z[m/s]  D_z[m^2/s] dlnn_z/dr[m^-1]  &
        & D_{E_r}[m^2/s] eE_r/T_z[m^-1]  D_T[m^2/s] dlnT_z/dr[m^-1]  D_n[m^2/s] dlnn_i/dr[m^-1]&
        & \Gamma_{anisotrp}/n_z[m/s]")')
   END IF
  
END SUBROUTINE INIT_FILES


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

