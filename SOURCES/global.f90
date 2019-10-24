!Global variables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE GLOBAL

  IMPLICIT NONE

  !Maximum number of points in (alpha,lambda), wells, points
  INTEGER, PARAMETER :: nlambdax  =1024
  INTEGER, PARAMETER :: nax       =128
  INTEGER, PARAMETER :: npointx   =400000
  INTEGER, PARAMETER :: nwx       =5000
  INTEGER, PARAMETER :: nsamp     =8

  !Maximum number of species and radial positions
  INTEGER, PARAMETER :: nsx       =141
  INTEGER, PARAMETER :: nbx       =11
  INTEGER, PARAMETER :: nbulk     =2

  !Bounce integrals
  INTEGER, PARAMETER :: nq0=6
  INTEGER, PARAMETER :: nqv=3

  !Algebraic constants
  REAL*8, PARAMETER :: PI          = 3.14159265358979
  REAL*8, PARAMETER :: TWOPI       = 6.28318530717959
  REAL*8, PARAMETER :: SQPI        = 1.77245385090552
  REAL*8, PARAMETER :: SQ2         = 1.4142135623731
  REAL*8, PARAMETER :: ZERO        = 0.00000000000000
  REAL*8, PARAMETER :: ONE         = 1.00000000000000
  REAL*8, PARAMETER :: MONE        =-1.00000000000000
  REAL*8, PARAMETER :: ALMOST_ZERO = 1.00000000000E-9
  REAL*8, PARAMETER :: SMALL       = 1.00000000000E-7
  COMPLEX*16, PARAMETER :: NUMI    = CMPLX(0.0,1.0)

  !Physics constants
  REAL*8, PARAMETER :: m_e =1.04396844E-08 !ratio between unitary mass and charge
  REAL*8, PARAMETER ::   e =1.60217662e-19 !electron charge

  !MPI constants
  INTEGER myrank,numprocs

  !Input parameters
  !Namellist /param/
  LOGICAL GEN_FLAG(10),DEBUG,TIME,PLOTG,USE_B1,USE_B0pB1,REMOVE_DIV,DELTA,FAST_AMB
  INTEGER I0,TRUNCATE_B,NER,NERR,MLAMBDA,MAL
  REAL*8 PREC_B,PREC_EXTR,PREC_BINT,PREC_DQDV,PREC_INTV,IPERR,ERMIN,ERMAX,ERACC
  !Namellist /model/
  CHARACTER*100 DIRDB
  CHARACTER*7 DIRS(nsx)
  LOGICAL CALC_DB,ONLY_DB,INC_EXB,TANG_VM,CLASSICAL,FRICTION,ANISOTROPY
  LOGICAL SOLVE_AMB,TRIVIAL_AMB,SOLVE_QN,TRIVIAL_QN,ZERO_PHI1,ONLY_PHI1,COMPARE_MODELS
  LOGICAL NEOTRANSP,TASK3D,PENTA,SATAKE,ANA_NTV,JPP
  INTEGER D_AND_V
  REAL*8 FN,FI,FS,FP,FE,FR,FB,FACT_CON,FNE,FTE,FTI,FDLNE,FDLTE,FDLTI,FER
  !Global
  LOGICAL DKES_READ,PHI1_READ,CALCULATED_INT,PRE_INTV,NTV,QN,USE_B0,CONVERGED,TRACE_IMP,PLATEAU_OR_PS
  LOGICAL IMP1NU,CLOSEST_LAMBDA,CLOSEST_LAMBDA1,CENTERED_ALPHA,FIRST_ORDER_ALPHA,STELL_ANTISYMMETRIC
  LOGICAL, PARAMETER :: PLOT_XYZ=.FALSE.

  !Flux-surface constants
  REAL*8 rad_a,rad_R,eps,eps32,avB2,etet
  REAL*8 atorflux,psip,iota,siota,aiota,iota2,diotadpsi,Bzeta,Btheta,torflux
  REAL*8 iBtpBz,aiBtpBz,sgnB,dBzdpsi,dBtdpsi,dB0dpsi

  !Magnetic field map
  INTEGER, PARAMETER :: mpolbd=128
  INTEGER, PARAMETER :: ntorbd=128
  INTEGER, PARAMETER :: nmd=2*mpolbd*ntorbd/nsamp
  INTEGER mpolb,ntorb,nzperiod,Nnm,Nnmp
  INTEGER, PARAMETER :: rprec   = SELECTED_REAL_KIND(12,100)
  REAL(rprec) borbi(-ntorbd:ntorbd,0:mpolbd)
  REAL(rprec) borbic(-ntorbd:ntorbd,0:mpolbd)  ,borbis(-ntorbd:ntorbd,0:mpolbd)
  REAL(rprec) borbic0(-ntorbd:ntorbd,0:mpolbd),borbis0(-ntorbd:ntorbd,0:mpolbd)
  REAL(rprec) dborbicdpsi(-ntorbd:ntorbd,0:mpolbd),dborbisdpsi(-ntorbd:ntorbd,0:mpolbd)
  REAL*8 rorbic(-ntorbd:ntorbd,0:mpolbd),rorbis(-ntorbd:ntorbd,0:mpolbd)
  REAL*8 zorbic(-ntorbd:ntorbd,0:mpolbd),zorbis(-ntorbd:ntorbd,0:mpolbd)
  REAL*8 porbic(-ntorbd:ntorbd,0:mpolbd),porbis(-ntorbd:ntorbd,0:mpolbd)
  REAL*8 bnmc(nmd),bnmc0(nmd),bnmc1(nmd),dbnmcdpsi(nmd)
  REAL*8 bnms(nmd),bnms0(nmd),bnms1(nmd),dbnmsdpsi(nmd)
  REAL*8 mp(nmd),np(nmd),ext_mp(nmd),ext_np(nmd)
  REAL*8, ALLOCATABLE :: absnablar(:,:)

  !Grid, resolution and experimetal points
  REAL*8 zmax,tmax,dzstep
  INTEGER array(nax,nax)

  !FFTs
  INTEGER*8, SAVE :: plan_fwd=0
  INTEGER*8, SAVE :: plan_bwd=0

  !Velocity integral
  INTEGER, PARAMETER :: nv=28
  INTEGER, PARAMETER :: iv0=6   !Representative velocity (v~=v_th) 
  REAL*8 v(nv),weight(nv),Sdke(nv),vdconst(nv),vmconst(nv),fdkes(nv),nu(nv),mu(nv),ftrace1nu(nv),mmuoT(nv)
  REAL*8 nuth(nbx),vth(nbx),nuzi(nbx)

  !Collisionality and normalized radial electric field
  REAL*8 cmul_PS,cmul_plateauPS,cmul_1NU
  REAL*8 D11onu,D11pla,D11nu

  !DKES/neotransp/PENTA-related variables
  INTEGER ncmult,nefieldt,nvmagt
  REAL*8, ALLOCATABLE ::  cmult(:), efieldt(:), vmagt(:)
  REAL*8, ALLOCATABLE :: lcmult(:),lefieldt(:),lvmagt(:)
  REAL*8, ALLOCATABLE :: lD11dkes1(:,:),lD11dkes2(:,:),lD11tab(:,:,:)


END MODULE GLOBAL
