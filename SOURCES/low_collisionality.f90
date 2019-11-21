
!Calculate neoclassical transport at low collisionalities
!Include 1/nu,sqrtnu and superbanana-plateau; do not need splitting B=B_0+delta*B_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_LOW_COLLISIONALITY(jv,Epsi,phi1c,Mbbnm,trMnm,D11,&
     & nalphab,zeta,theta,dn1dv,dn1nm)

!--------------------------------------------------------------------------------------------- 
!Calculate monoenergetic transport coefficient D11 and contribution quasineutrality dn1dv and 
!dn1nm at nalphabxnalphab grid in (zeta,theta)  for collisionality cmul=nu(jv)/v(jv) and 
!normalized radial electric field Epsi/v(jv), and in the presence of phi1c
!-----------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
!  Input
  INTEGER jv
  REAL*8 Epsi,phi1c(Nnmp),Mbbnm(Nnmp),trMnm(Nnmp)
  !Output
  INTEGER nalphab
  REAL*8 D11(Nnmp,Nnmp),zeta(nax),theta(nax),dn1dv(nax,nax),dn1nm(Nnmp,Nnmp)
  !Others
  CHARACTER*100 serr
  INTEGER ial,ilambda
  INTEGER, SAVE :: nal,nlambda
  REAL*8 D11r(100,100)

  IF(.NOT.QN) CONVERGED=.FALSE.

  IF(CONVERGED.OR.(MLAMBDA.GT.0.AND.MAL.GT.0)) THEN

     IF(MLAMBDA.GT.0.AND.MAL.GT.0) THEN
        nlambda=MLAMBDA
        nal=MAL
        CONVERGED=.TRUE.
     END IF
     !Write monoenergetic transport coefficients using DKES normalization
     CALL CALC_LOW_COLLISIONALITY_NANL(nal,nlambda,jv,Epsi,phi1c,Mbbnm,trMnm,&
          & D11,nalphab,zeta,theta,dn1dv,dn1nm)
     WRITE(200+myrank,'(3(1pe13.5)," NaN ",2(1pe13.5)," &
          & NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN")') &
          & nu(jv)/v(jv)/2.,Epsi/v(jv)*psip,vmconst(jv)/v(jv),fdkes(jv)*D11(1,1),fdkes(jv)*D11(1,1)
  ELSE

     D11=0
     CONVERGED=.FALSE.
     ial=1
     nal=32
     DO WHILE (nal.LE.nax)
        ilambda=1
        nlambda=32
        DO WHILE (nlambda.LE.nlambdax)
           CALL CALC_LOW_COLLISIONALITY_NANL(nal,nlambda,jv,Epsi,phi1c,Mbbnm,trMnm,&
                & D11r(ial,ilambda),nalphab,zeta,theta,dn1dv,dn1nm)
           IF(ilambda.GT.1.AND.ABS(D11r(ial,ilambda)/D11r(ial,ilambda-1)-1.0).LT.PREC_DQDV) THEN
              EXIT
           END IF
           ilambda=ilambda+1
           nlambda=nlambda*2  
        END DO
        D11r(ial,ilambda+1:100)=D11r(ial,ilambda)
        IF(ial.GT.1.AND.nlambda.LE.nlambdax.AND.&
             & ABS(D11r(ial,ilambda)/D11r(ial-1,100)-1.0).LT.PREC_DQDV) THEN
           WRITE(200+myrank,'(3(1pe13.5)," NaN ",2(1pe13.5),&
                & " NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN")') &
                & nu(jv)/v(jv)/2.,Epsi/v(jv)*psip,vmconst(jv)/v(jv),fdkes(jv)*D11r(ial,ilambda),fdkes(jv)*D11r(ial,ilambda)
           EXIT
        END IF
        ial=ial+1
        nal=nal*2
     END DO
     
     IF(nal.GT.nax.OR.nlambda.GT.nlambdax) THEN
        IF(nal.GE.nax) THEN
           WRITE(1100+myrank,*) 'Increase nax'
        ELSE IF(nlambda.GE.nlambdax) THEN
           WRITE(1100+myrank,*) 'Increase nlambdax'
        END IF
        serr="DKE not converged"
        IF(.NOT.ONLY_DB) CALL END_ALL(serr,.FALSE.)
     ELSE
        CONVERGED=.TRUE.
        D11=D11r(ial,ilambda)
        WRITE(1000+myrank,'(" DKE converged with nal=",I6," and nlambda=",I6)') nal,nlambda
        CALL CALC_LOW_COLLISIONALITY_NANL(nal,nlambda,jv,Epsi,phi1c,Mbbnm,trMnm,&
             D11,nalphab,zeta,theta,dn1dv,dn1nm)
     END IF

  END IF

END SUBROUTINE CALC_LOW_COLLISIONALITY


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_LOW_COLLISIONALITY_NANL(nal,nlambda,jv,Epsi,phi1c,Mbbnm,trMnm,&
     & D11,nalphab,zeta,theta,dn1dv,dn1nm)

!--------------------------------------------------------------------------------------------- 
!Calculate monoenergetic transport coefficient D11 and contribution quasineutrality dn1dv and 
!dn1nm at nalphabxnalphab grid in (zeta,theta), for collisionality cmul=nu(jv)/v(jv) and normalized 
!radial electric field Epsi/v(jv)m and in the presence of phi1c
!The DKE is solved in a nalxnlambda grid in (alpha,lambda)
!-----------------------------------------------------------------------------------------------

  USE GLOBAL
#ifdef IPPorNIFS
  USE petscsys
  USE petscksp
!  USE petscvec
!  USE petscmat
!  USE petscpc
#endif
  IMPLICIT NONE
!  Input
  INTEGER jv,nlambda,nal
  REAL*8 Epsi,phi1c(Nnmp),Mbbnm(Nnmp),trMnm(Nnmp),zeta(nax),theta(nax)
  !Output
  INTEGER nalphab
  REAL*8 D11(Nnmp,Nnmp),dn1dv(nax,nax),dn1nm(Nnmp,Nnmp)
  !Others
  INTEGER, SAVE :: nalpha,nalphab_save,npoint
  INTEGER, SAVE, ALLOCATABLE :: nbif(:),i_p_ap1(:),i_p_am1(:),i_p_ap2(:),i_p_am2(:)
  REAL*8, SAVE :: offset,theta_save(nax)
  REAL*8, SAVE, ALLOCATABLE :: lambda(:),dalpha_am1(:),dalpha_ap1(:),dalpha_am2(:),dalpha_ap2(:)
  REAL*8, SAVE, ALLOCATABLE :: dlambda_lm1(:),dlambda_lp1(:)
  !Wells and bounce points
  LOGICAL, ALLOCATABLE :: connected(:,:),bottom(:),ltemp(:),ltemp2(:,:)
  INTEGER nw,na
  REAL*8, ALLOCATABLE :: z1(:),t1(:),B1(:),hBpp1(:),vd1(:,:)
  REAL*8, ALLOCATABLE :: zb(:),tb(:),Bb(:),hBppb(:),vdb(:,:) 
  REAL*8, ALLOCATABLE :: z2(:),t2(:),B2(:),hBpp2(:),vd2(:,:) 
  REAL*8, ALLOCATABLE :: alphap_w(:),Bt(:),Btt(:),temp1(:),temp2(:,:)
  !Angular and lambda grid
  INTEGER, ALLOCATABLE :: i_w(:),i_a(:),itemp(:)
  INTEGER, SAVE, ALLOCATABLE :: i_l(:),i_p(:,:,:)
  REAL*8, SAVE, ALLOCATABLE :: thetap(:,:),zetap(:),zetax(:,:),B_al(:,:),vds_al(:,:,:)
  REAL*8, ALLOCATABLE :: one_o_lambda(:),alphap(:),lambdab_w(:),lambdac_w(:)
  REAL*8 dlambdap,dalphap
  !Lambda neighbours
  INTEGER, PARAMETER :: nbifx=5  !maximum number of bifurcations allowed
  INTEGER, SAVE, ALLOCATABLE :: i_p_lm1(:),i_p_lp1(:,:)
  !Alpha neighbours
  REAL*8, ALLOCATABLE, SAVE :: zlw(:),zrw(:)
  REAL*8, ALLOCATABLE       :: tlw(:),trw(:)
  INTEGER ipoint,ila
!!$  !Second adiabatic invariant
!!$  LOGICAL, ALLOCATABLE :: readpoint(:)
!!$  INTEGER ja,ia,il
!!$  REAL*8, ALLOCATABLE :: Jnorm(:,:)
  !Drift-kinetic equation
  INTEGER, SAVE, ALLOCATABLE :: nnz(:)
  REAL*8, SAVE, ALLOCATABLE :: BI1(:),BI2(:),BI3(:),BI4(:),BI5(:),BI6(:),BI7(:,:),gint(:,:)
  REAL*8 omega,cmul
#ifdef MPIandPETSc
  !Petsc
#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscmat.h>
#include <petsc/finclude/petscpc.h>
#include <petsc/finclude/petscksp.h>
!#include <petsc/finclude/petscviewer.h>
!  PetscViewer, SAVE :: viewer
  PetscErrorCode ierr
  PetscInt, SAVE :: innz(0:npointx-1)
  KSP, SAVE :: ksp
  Mat, SAVE :: matCOL,matVEAf,matVEAb,matVMAf,matVMAb
  INTEGER ipointm1,jpointm1,jpoint
#else
  REAL*8, SAVE, ALLOCATABLE :: COL(:,:),VEAf(:,:),VEAb(:,:),VMAf(:,:),VMAb(:,:)
  REAL*8, ALLOCATABLE :: mat(:,:)
#endif
  REAL*8, ALLOCATABLE :: rowCOL(:),rowVEAf(:),rowVEAb(:),rowVMAf(:),rowVMAb(:)
  !Time
  CHARACTER*30, PARAMETER :: routine="CALC_LOW_COLLISIONALITY"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart
!  CHARACTER*30, PARAMETER :: routine2="FILL_GRID"
!  INTEGER, SAVE :: ntotal2=0
!  REAL*8,  SAVE :: ttotal2=0
!  REAL*8,  SAVE :: t02=0
!  REAL*8 tstart2


  Mbbnm=Mbbnm !To be removed
  trMnm=trMnm

  CALL CPU_TIME(tstart)

  WRITE(1000+myrank,*) 'Calculating low collisionality regimes'

  !Set output to zero

  IF(CALCULATED_INT) GOTO 123

!  CALL CPU_TIME(tstart2)

  IF(PHI1_READ) bnmc0=bnmc0+2*borbic(0,0)*phi1c/vdconst(jv) 

 !Find and characterize wells
  ALLOCATE(connected(nwx,nwx),bottom(nwx),&
       & z1(nwx),t1(nwx),B1(nwx),hBpp1(nwx),vd1(nqv,nwx),&
       & zb(nwx),tb(nwx),Bb(nwx),hBppb(nwx),vdb(nqv,nwx),&
       & z2(nwx),t2(nwx),B2(nwx),hBpp2(nwx),vd2(nqv,nwx),& 
       & alphap_w(nwx),Bt(nwx),Btt(nwx),&
       & lambdab_w(nwx),lambdac_w(nwx))
  CALL CHARACTERIZE_WELLS(nal,na,nalpha,nw,z1,t1,B1,hBpp1,vd1, &
       & zb,tb,Bb,hBppb,vdb, &
       & z2,t2,B2,hBpp2,vd2, &
       & Bt,Btt,alphap_w,dalphap,bottom,connected,offset)
  !Resize arrays (nwx->nw)
  ALLOCATE(temp1(nw),temp2(nqv,nw),ltemp(nw))
  temp1=z1(1:nw);DEALLOCATE(z1);ALLOCATE(z1(nw));z1=temp1
  temp1=zb(1:nw);DEALLOCATE(zb);ALLOCATE(zb(nw));zb=temp1
  temp1=z2(1:nw);DEALLOCATE(z2);ALLOCATE(z2(nw));z2=temp1
  temp1=t1(1:nw);DEALLOCATE(t1);ALLOCATE(t1(nw));t1=temp1
  temp1=tb(1:nw);DEALLOCATE(tb);ALLOCATE(tb(nw));tb=temp1
  temp1=t2(1:nw);DEALLOCATE(t2);ALLOCATE(t2(nw));t2=temp1
  temp1=B1(1:nw);DEALLOCATE(B1);ALLOCATE(B1(nw));B1=temp1
  temp1=Bb(1:nw);DEALLOCATE(Bb);ALLOCATE(Bb(nw));Bb=temp1
  temp1=B2(1:nw);DEALLOCATE(B2);ALLOCATE(B2(nw));B2=temp1
  temp1=hBpp1(1:nw);DEALLOCATE(hBpp1);ALLOCATE(hBpp1(nw));hBpp1=temp1
  temp1=hBppb(1:nw);DEALLOCATE(hBppb);ALLOCATE(hBppb(nw));hBppb=temp1
  temp1=hBpp2(1:nw);DEALLOCATE(hBpp2);ALLOCATE(hBpp2(nw));hBpp2=temp1
  temp1=alphap_w(1:nw);DEALLOCATE(alphap_w);ALLOCATE(alphap_w(nw));alphap_w=temp1
  temp1=lambdab_w(1:nw);DEALLOCATE(lambdab_w);ALLOCATE(lambdab_w(nw));lambdab_w=temp1
  temp1=lambdac_w(1:nw);DEALLOCATE(lambdac_w);ALLOCATE(lambdac_w(nw));lambdac_w=temp1
  temp1= Bt(1:nw);DEALLOCATE(Bt) ;ALLOCATE(Bt(nw)) ;Bt=temp1
  temp1=Btt(1:nw);DEALLOCATE(Btt);ALLOCATE(Btt(nw));Btt=temp1
  ltemp=bottom(1:nw);DEALLOCATE(bottom);ALLOCATE(bottom(nw));bottom=ltemp
  temp2=vd1(1:nqv,1:nw);DEALLOCATE(vd1);ALLOCATE(vd1(nqv,nw));vd1=temp2
  temp2=vdb(1:nqv,1:nw);DEALLOCATE(vdb);ALLOCATE(vdb(nqv,nw));vdb=temp2
  temp2=vd2(1:nqv,1:nw);DEALLOCATE(vd2);ALLOCATE(vd2(nqv,nw));vd2=temp2
  DEALLOCATE(temp1,temp2,ltemp)  
  ALLOCATE(ltemp2(nw,nw))
  ltemp2=connected(1:nw,1:nw);DEALLOCATE(connected);ALLOCATE(connected(nw,nw));connected=ltemp2
  DEALLOCATE(ltemp2)
  !Create grid in alpha, then (zeta,theta)
  !Determine number of modes
  nalphab=1
  DO WHILE(nalphab.LT.nalpha)
     nalphab=nalphab*2
  END DO
  nalphab=nalphab/2
  nalphab_save=nalphab
  IF(ALLOCATED(zetap)) DEALLOCATE(zetap,thetap,zetax,B_al,vds_al)
  ALLOCATE(zetax(nalpha,nalphab),alphap(nalpha),&
         & zetap(nalphab),thetap(nalpha,nalphab),&
         & B_al(nalpha,nalphab),vds_al(Nnmp,nalpha,nalphab))
  CALL CREATE_ANGULAR_GRID(na,nalpha,nalphab,alphap,dalphap,offset,thetap,zetap,zetax,B_al,vds_al)
  zeta(1:nalphab) =zetap  !square grid
  theta(1:nalphab)=zetap*nzperiod 
  theta_save=theta

  IF(ALLOCATED(lambda)) DEALLOCATE(lambda,i_p)
  ALLOCATE(lambda(nlambda),one_o_lambda(nlambda),i_p(nlambda,nalpha,nalphab))

  !Set global grid in lambda
  CALL CREATE_LAMBDA_GRID(nlambda,nw,Bb,Btt,&
       & lambdab_w,lambdac_w,dlambdap,lambda,one_o_lambda)  !CHECK: Btt instead of Bt?

  !For each point in the (zeta,theta) grid, determine well and absolute point
  !For each absolute point, determine alpha, lambda and well number
  IF(ALLOCATED(i_l)) DEALLOCATE(i_l)
  ALLOCATE(i_l(npointx),i_a(npointx),i_w(npointx))
  CALL LABEL_GRIDPOINTS(na,nalpha,nalphab,nlambda,nw,bottom,connected,&
       & alphap_w,z1,zb,z2,Bb,Bt,lambda,&
       & thetap,zetap,zetax,B_al,npoint,i_a,i_l,i_w,i_p)
  ALLOCATE(itemp(npoint))
  itemp=i_l(1:npoint);DEALLOCATE(i_l);ALLOCATE(i_l(npoint));i_l=itemp
  itemp=i_a(1:npoint);DEALLOCATE(i_a);ALLOCATE(i_a(npoint));i_a=itemp
  itemp=i_w(1:npoint);DEALLOCATE(i_w);ALLOCATE(i_w(npoint));i_w=itemp
  DEALLOCATE(itemp)

  IF(ALLOCATED(BI1)) THEN
     DEALLOCATE(nbif,i_p_ap1,i_p_am1,i_p_ap2,i_p_am2,i_p_lp1,i_p_lm1,&
          & dalpha_am1,dalpha_ap1,dalpha_am2,dalpha_ap2,dlambda_lm1,dlambda_lp1,&
          & BI1,BI2,BI3,BI4,BI5,BI6,BI7,zlw,zrw,nnz)
#ifdef MPIandPETSc
     CALL MatDestroy(matCOL,ierr)
     CALL MatDestroy(matVEAf,ierr)
     CALL MatDestroy(matVEAb,ierr)
     IF(TANG_VM) THEN
        CALL MatDestroy(matVMAf,ierr)
        CALL MatDestroy(matVMAb,ierr)
     END IF
#endif
  END IF
  ALLOCATE(nbif(npoint),i_p_ap1(npoint),i_p_am1(npoint),i_p_ap2(npoint),i_p_am2(npoint))
  ALLOCATE(i_p_lm1(npoint),i_p_lp1(nbifx,npoint))
  ALLOCATE(dalpha_am1(npoint),dalpha_ap1(npoint),dalpha_am2(npoint),dalpha_ap2(npoint))
  ALLOCATE(dlambda_lm1(npoint),dlambda_lp1(npoint))
  ALLOCATE(BI1(npoint),BI2(npoint),BI3(npoint))
  ALLOCATE(BI4(npoint),BI5(npoint),BI6(npoint),BI7(npoint,Nnmp))
  ALLOCATE(zlw(npoint),zrw(npoint),tlw(npoint),trw(npoint)) 
  ALLOCATE(nnz(npoint))

  !Find neighbours in lambda
  CALL FIND_LAMBDA_NEIGHBOURS(npoint,nalpha,nalphab,nlambda,nw,nbifx,i_p,&
       & bottom,i_w,i_p_lm1,nbif,i_p_lp1)
  !Correct delta lambda at the bottom
  dlambda_lm1=dlambdap 
  dlambda_lp1=dlambdap 
  DO ipoint=1,npoint-1
     IF(i_p_lp1(1,ipoint).EQ.npoint) THEN 
        dlambda_lp1(ipoint)=lambdab_w(i_w(ipoint))-lambda(i_l(ipoint)) 
     ELSE IF(i_p_lm1(ipoint).EQ.npoint) THEN                        
        dlambda_lm1(ipoint)=lambda(i_l(ipoint))-MINVAL(lambdac_w) 
     END IF
  END DO

  !Calculate coefficients of the drift kinetic equation
  CALL COEFFICIENTS_DKE(npoint,i_w,i_l,nw,&
                 &  z1,t1,B1,hBpp1,vd1,&
                 &  zb,tb,Bb,hBppb,vdb,&
                 &  z2,t2,B2,hBpp2,vd2,&
                 &  nlambda,lambda,zlw,tlw,zrw,trw,& 
                 &  BI1,BI2,BI3,BI4,BI5,BI6,Nnmp,BI7)

!!$  ALLOCATE(Jnorm(nlambda,nalpha))
!!$  ALLOCATE(readpoint(npoint)) 
!!$  Jnorm=0
!!$  DO ila=1,nlambda
!!$     DO ia=1,na
!!$        readpoint=.FALSE.
!!$        DO ja=ia,nalpha,na
!!$           DO il=1,nalphab
!!$  
!!$              ipoint=i_p(ila,ja,il)
!!$              IF(readpoint(ipoint)) CYCLE
!!$              Jnorm(ila,ia)=Jnorm(ila,ia)+BI6(ipoint)
!!$              readpoint(ipoint)=.TRUE.
!!$           END DO
!!$        END DO
!!$        print *,'J',lambda(ila),ia,Jnorm(ila,ia)
!!$     END DO
!!$  END DO
!!$  stop

  !Find neighbours in alpha
  CALL FIND_ALPHA_NEIGHBOURS(npoint,i_w,i_a,i_l,i_p_lm1,zlw,tlw,zrw,trw,nw,&
  & z1,zb,z2,tb,Bb,Bt,alphap_w,bottom,connected,&
  & nlambda,one_o_lambda,na,nalpha,alphap,dalphap,offset, &
  & i_p_am1,i_p_ap1,dalpha_am1,dalpha_ap1,&
  & i_p_am2,i_p_ap2,dalpha_am2,dalpha_ap2)

  DEALLOCATE(connected,bottom,z1,t1,B1,hBpp1,vd1,zb,tb,Bb,hBppb,vdb,&
       &    z2,t2,B2,hBpp2,vd2,alphap_w,Bt,Btt,lambdab_w,lambdac_w) 

  !Find non-zero elements of the DKE matrix and initialize
  IF(ALLOCATED(gint)) DEALLOCATE(gint)
  ALLOCATE(rowCOL(npoint),rowVEAf(npoint),rowVEAb(npoint),rowVMAf(npoint),rowVMAb(npoint),gint(npoint,Nnmp))
#ifdef MPIandPETSc
  DO ipoint=1,npoint
     nnz(ipoint)=1
     IF(ipoint.NE.npoint.AND.i_l(ipoint).NE.0) THEN
        CALL FILL_DKE_ROW(ipoint,npoint,&
             & dalpha_am1(ipoint),dalpha_ap1(ipoint),dalpha_am2(ipoint),dalpha_ap2(ipoint),&
             & i_p_am1(ipoint),i_p_ap1(ipoint),i_p_am2(ipoint),i_p_ap2(ipoint),&
             & i_l,lambda(ila),dlambda_lm1,dlambda_lp1,i_p_lm1,nbif,nbifx,i_p_lp1,&
             & BI1(ipoint),BI2,BI4(ipoint),BI5(ipoint),&
             & rowCOL,rowVEAf,rowVEAb,rowVMAf,rowVMAb,nnz(ipoint),.TRUE.)
     END IF
     innz(ipoint-1)=nnz(ipoint)
  END DO
  innz(npoint-1)=1
  CALL INIT_LINEAR_PROBLEM(npoint,nnz,matCOL,matVEAf,matVEAb,matVMAf,matVMAb)
#else
  IF(ALLOCATED(COL)) DEALLOCATE(COL,VEAf,VEAb,VMAf,VMAb)
  ALLOCATE(COL(npoint,npoint),VEAf(npoint,npoint),VEAb(npoint,npoint),VMAf(npoint,npoint),VMAb(npoint,npoint))
#endif

  !For each point number, fill one row of the matrix with quantities that depend on the configuration only
  DO ipoint=1,npoint
     rowCOL=0
     rowVEAf=0
     rowVEAb=0
     rowVMAf=0
     rowVMAb=0
     ila =i_l(ipoint)
     IF(ipoint.EQ.npoint.OR.ila.EQ.0) THEN
        rowCOL(ipoint)=1./dlambdap/dlambdap
     ELSE
        CALL FILL_DKE_ROW(ipoint,npoint,&
             & dalpha_am1(ipoint),dalpha_ap1(ipoint),dalpha_am2(ipoint),dalpha_ap2(ipoint),&
             & i_p_am1(ipoint),i_p_ap1(ipoint),i_p_am2(ipoint),i_p_ap2(ipoint),&
             & i_l,lambda(ila),dlambda_lm1,dlambda_lp1,i_p_lm1,nbif,nbifx,i_p_lp1,&
             & BI1(ipoint),BI2,BI4(ipoint),BI5(ipoint),&
             & rowCOL,rowVEAf,rowVEAb,rowVMAf,rowVMAb,nnz(ipoint),.FALSE.)
     END IF

#ifdef MPIandPETSc 
     ipointm1=ipoint-1    
     DO jpoint=1,npoint
        IF((ABS(rowCOL(jpoint)).GT.ZERO).OR.&
             & (ABS(rowVEAf(jpoint)).GT.ZERO).OR. &
             & (ABS(rowVEAb(jpoint)).GT.ZERO).OR. &
             & (TANG_VM.AND.ABS(rowVMAf(jpoint)).GT.ZERO).OR. &
             & (TANG_VM.AND.ABS(rowVMAb(jpoint)).GT.ZERO)) THEN
           jpointm1=jpoint-1
           CALL MatSetValues(matCOL,1,ipointm1,1,jpointm1,rowCOL(jpoint),INSERT_VALUES,ierr)
        END IF
        IF(ABS(rowVEAf(jpoint)).GT.ZERO) &
             & CALL MatSetValues(matVEAf,1,ipointm1,1,jpointm1,rowVEAf(jpoint),INSERT_VALUES,ierr)
        IF(ABS(rowVEAb(jpoint)).GT.ZERO) &
             & CALL MatSetValues(matVEAb,1,ipointm1,1,jpointm1,rowVEAb(jpoint),INSERT_VALUES,ierr)
        IF(TANG_VM.AND.ABS(rowVMAf(jpoint)).GT.ZERO) &
             & CALL MatSetValues(matVMAf,1,ipointm1,1,jpointm1,rowVMAf(jpoint),INSERT_VALUES,ierr)
        IF(TANG_VM.AND.ABS(rowVMAb(jpoint)).GT.ZERO) &
             & CALL MatSetValues(matVMAb,1,ipointm1,1,jpointm1,rowVMAb(jpoint),INSERT_VALUES,ierr)
     END DO
#else
     COL(ipoint,:) =rowCOL
     VEAf(ipoint,:)=rowVEAf
     VEAb(ipoint,:)=rowVEAb
     IF(TANG_VM) THEN
        VMAf(ipoint,:)=rowVMAf
        VMAb(ipoint,:)=rowVMAb
     END IF
#endif

  END DO
#ifdef MPIandPETSc 
  
  CALL MatAssemblyBegin(matCOL,MAT_FINAL_ASSEMBLY,ierr)
  CALL MatAssemblyEnd(  matCOL,MAT_FINAL_ASSEMBLY,ierr)
  CALL MatAssemblyBegin(matVEAf,MAT_FINAL_ASSEMBLY,ierr)
  CALL MatAssemblyEnd(  matVEAf,MAT_FINAL_ASSEMBLY,ierr)
  CALL MatAssemblyBegin(matVEAb,MAT_FINAL_ASSEMBLY,ierr)
  CALL MatAssemblyEnd(  matVEAb,MAT_FINAL_ASSEMBLY,ierr)
  IF(TANG_VM) THEN
     CALL MatAssemblyBegin(matVMAf,MAT_FINAL_ASSEMBLY,ierr)
     CALL MatAssemblyEnd(  matVMAf,MAT_FINAL_ASSEMBLY,ierr)
     CALL MatAssemblyBegin(matVMAb,MAT_FINAL_ASSEMBLY,ierr)
     CALL MatAssemblyEnd(  matVMAb,MAT_FINAL_ASSEMBLY,ierr)
  END IF
  
#endif

!  CALL CALCULATE_TIME(routine2,ntotal2,t02,tstart2,ttotal2)

123 nalphab=nalphab_save 
  theta=theta_save
  IF(DEBUG) THEN
     DO ipoint=1,npoint
        WRITE(3300+myrank,'(I6,1(1pe13.5),2I6)') ipoint,BI3(ipoint)*vdconst(jv)/nu(jv),nal,nlambda
     END DO
  END IF

  !Use linear combinations of precalculated matrices to fill actual matrix for given values
  !the collisionality, radial electric field, etc
#ifdef MPIandPETSc 
  CALL FILL_MATRIX_PETSC(matCOL,jv,Epsi,matVEAf,matVEAb,matVMAf,matVMAb,ksp)
#else
  ALLOCATE(mat(npoint,npoint))
  CALL FILL_MATRIX(npoint,COL,jv,Epsi,VEAf,VEAb,VMAf,VMAb,mat)
#endif
!Invert linear system
#ifdef MPIandPETSc
  CALL INVERT_MATRIX_PETSC(nalphab,jv,npoint,BI3,BI7,phi1c,ksp,gint)
#else
  CALL INVERT_MATRIX(nalphab,jv,npoint,BI3,BI7,phi1c,mat,gint)
  DEALLOCATE(mat)
#endif  

  IF(PHI1_READ.AND.IMP1NU) THEN
     gint(:,2)=gint(:,1)*ftrace1nu(jv)*EXP(-mmuoT(jv)/lambda(i_l(:)))*mmuoT(jv)/lambda(i_l(:))
     gint(:,1)=gint(:,1)*ftrace1nu(jv)*EXP(-mmuoT(jv)/lambda(i_l(:)))
     WRITE(1000+myrank,*) 'Using mu and lambda'
  END IF

  !Calculate dn_1dv and D_{11} integrating it the velocity space and taking the flux-surface-average
  CALL INTEGRATE_G(nalpha,nalphab,nlambda,lambda,i_p,npoint,gint, &
       & zetap,thetap,theta(1:nalphab),B_al,vds_al,D11,dn1dv(1:nalphab,1:nalphab),dn1nm)

  IF(.NOT.PHI1_READ.OR..NOT.IMP1NU) THEN
     D11(1,:)=D11(1,:)*vdconst(jv)
     IF(PHI1_READ) D11     =D11     *weight(jv)
  END IF

!  ELSE
!     gint=gint*weight(jv)
!  END IF
  omega=ABS(Epsi)*psip/v(jv)
  cmul=nu(jv)/v(jv)/2.  

  !Connection with plateau regime
  IF(FACT_CON.GT.0.AND.cmul_1NU.GT.0.&
       & .AND.cmul.GT.cmul_1NU/FACT_CON.AND.omega.LT.1E-2) D11=D11+D11pla/fdkes(jv)
  
  !Write monoenergetic transport coefficients using DKES normalization
  IF(DEBUG.OR.ONLY_DB) THEN
     CALL FLUSH(10000+myrank)
     IF(cmul_1NU.GT.0) THEN
        WRITE(10000+myrank,'("4 ",6(1pe13.5),3I5,1pe13.5,I5)') nu(jv)/v(jv)/2.,Epsi/v(jv)*psip,vmconst(jv)/v(jv),&
             & fdkes(jv)*D11(1,1),weight(jv)/fdkes(jv),weight(jv)/fdkes(jv)*v(jv)*v(jv),nal,nlambda,nalphab,D11pla
     ELSE
        WRITE(10000+myrank,'("0 ",6(1pe13.5),3I5,1pe13.5,I5)') nu(jv)/v(jv)/2.,Epsi/v(jv)*psip,vmconst(jv)/v(jv),&
             & fdkes(jv)*D11(1,1),weight(jv)/fdkes(jv),weight(jv)/fdkes(jv)*v(jv)*v(jv),nal,nlambda,nalphab,D11pla
     END IF
     CALL FLUSH(10000+myrank)
  END IF
  
  IF(PHI1_READ) bnmc0=bnmc0-2*borbic(0,0)*phi1c/vdconst(jv)

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE CALC_LOW_COLLISIONALITY_NANL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CHARACTERIZE_WELLS(nal,na,nalpha,nw,&
                     & z1,t1,B1,hBpp1,vd1, &
                     & zb,tb,Bb,hBppb,vdb, &
                     & z2,t2,B2,hBpp2,vd2, &
                     & Bt,Btt,alphap_w,dalphap,bottom,connected,offset)

!-----------------------------------------------------------------------------------------------
!Find and characterize nw wells (including matched regions) characterized by the Boozer
!toroidal and poloidal angles z and t, and alpha, of its top points, 1 and 2, and bottom
!b, as well as the value of the half-second derivative of the magnetic field strength along
!the magnetic field line hBpp, the value of the magnetic drift vd, alphap_w, their top according
!to two definitions Bt and Btt, whether they are the the bottom.
!Determine a map of connections, connected, and offset TODO
!It exits when nalpha.GT.nal
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL  
  IMPLICIT NONE
  !Input
  INTEGER nal
  !Output
  LOGICAL matched(nwx),connected(nwx,nwx),bottom(nwx)
  INTEGER na,nalpha,nw
  REAL*8 z1(nwx),t1(nwx),B1(nwx),hBpp1(nwx),vd1(nqv,nwx)
  REAL*8 zb(nwx),tb(nwx),Bb(nwx),hBppb(nwx),vdb(nqv,nwx) 
  REAL*8 z2(nwx),t2(nwx),B2(nwx),hBpp2(nwx),vd2(nqv,nwx) 
  REAL*8 Bt(nwx),Btt(nwx),alphap_w(nwx),dalphap,offset
  !Others
  CHARACTER*100 serr
  INTEGER iw,jw,kw,lw,nwperiods,iturn,flag,fath(nwx),nfound,nfoundt,nw0,nw1,nw2,nwmax
  REAL*8 maxBt,maxBb
  !Time
  CHARACTER*30, PARAMETER :: routine="CHARACTERIZE_WELLS"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  na=3  !Follow na field lines along several periods and collapses them into one period;
        !need na>=3 when using i+/-1 and i+/-2 for the derivatives at i;
  nalpha=INT(nzperiod*na/aiota+1)
  IF(nalpha.GT.nal) THEN
     serr="increase nal"
!     CALL END_ALL(serr,.FALSE.)
  END IF 
  nalpha=0
  !Initialize variables
  connected=.FALSE.
  matched  =.FALSE.
  bottom   =.FALSE.
  iturn=0
  nw0=1
  nw2=0
  nwmax=0
  nfoundt=0
  fath=0
  alphap_w=-1000.0 
  !Maximum possible number of wells determined by size of arrays
  DO WHILE((nalpha.LT.nal.AND.(nwmax.LE.nwx.OR.(NTV.AND.nwmax.LE.nwx))).OR.iturn.NE.0)
     iturn=iturn+1
     !Find wells along the field lines
     IF(iturn.EQ.1) THEN        
        nw0=nw2+1
        !Characterize wells with extreme values
        CALL FIND_WELLS(na,nw0,z1,t1,B1,hBpp1,vd1, &
             &                 zb,tb,Bb,hBppb,vdb, &
             &                 z2,t2,B2,hBpp2,vd2, &
             &           nfound,alphap_w,offset)
        nw1=nw0
        nw2=nw0+nfound-1
        DO iw=nw1,nw2
           Bt(iw) =MIN(B1(iw),B2(iw))!smallest maximum
           Btt(iw)=MAX(B1(iw),B2(iw))!largest  maximum
           bottom(iw)=.TRUE.
        END DO
        nfoundt=nfoundt+nfound !to nfoundt local mimima, by matching
        nwmax=4*nfoundt        !we may obtain up to 2*nfoundt regions
                               !exit loop if nwmax.GT.NWX
     !Match wells  
     ELSE
        kw=nw2+1
        nw1=kw 
        !Try to match well iw...
        DO iw=nw0,nw2 
           IF(matched(iw)) CYCLE 
           !...with well jw
           DO jw=nw0,nw2
              IF(jw.EQ.iw.OR.matched(jw)) CYCLE 
              CALL MATCH_WELLS(kw, & 
                   & z1(iw),t1(iw),B1(iw),hBpp1(iw),vd1(:,iw), &
                   & zb(iw),tb(iw),Bb(iw),hBppb(iw),vdb(:,iw), &
                   & z2(iw),t2(iw),B2(iw),hBpp2(iw),vd2(:,iw), &
                   & z1(jw),t1(jw),B1(jw),hBpp1(jw),vd1(:,jw), &
                   & zb(jw),tb(jw),Bb(jw),hBppb(jw),vdb(:,jw), &
                   & z2(jw),t2(jw),B2(jw),hBpp2(jw),vd2(:,jw), &
                   & z1(kw),t1(kw),B1(kw),hBpp1(kw),vd1(:,kw), &
                   & zb(kw),tb(kw),Bb(kw),hBppb(kw),vdb(:,kw), &
                   & z2(kw),t2(kw),B2(kw),hBpp2(kw),vd2(:,kw), &
                   & 1,flag)
              !If they match, set limits, matrix of relations....
              IF(flag.EQ.1) THEN 
                 Bt(kw) =MIN(B1(kw),B2(kw))
                 Btt(kw)=MAX(B1(kw),B2(kw))
                 nwperiods=INT(ABS(z2(kw)-z1(kw))*nzperiod/TWOPI/1.5)+1
                 IF(DEBUG) WRITE(2900+myrank,'("    Region ",I6," (from ",I6," and" ,I6,   &
                      & "contained in",I6," periods")') kw,iw,jw,nwperiods
                 matched(iw)=.TRUE. 
                 matched(jw)=.TRUE. 
                 fath(iw)=kw
                 fath(jw)=kw
                 alphap_w(kw)=alphap_w(iw)
                 connected(kw,iw)=.TRUE. !needed for determing that well iw
                 connected(kw,jw)=.TRUE. !cover alphas that jw and kw covered
                 DO lw=1,nw            
                    IF(connected(iw,lw).OR.connected(jw,lw)) connected(kw,lw)=.TRUE.
                 END DO
                 kw=kw+1
              END IF
           END DO
        END DO
        nw2=kw-1 
     END IF
     nw=nw2
     
     !If there were matches, continue matching for the same alphas
     IF(.NOT.NTV.AND.nw2.GE.nw1) THEN

        IF(DEBUG) THEN
           maxBt=-1E3
           DO iw=1,nw
              IF(.NOT.matched(iw).AND.Bt(iw).GT.maxBt) maxBt=Bt(iw)
           END DO
           WRITE(2900+myrank,'(" Labelled regions ",I6," to ",I6," out of a maximum of ", &           
                & I6,":",1pe13.5," < 1/lambda <",1pe13.5)')                     &
                & nw1,nw2,nwmax,MINVAL(Bb),maxBt
        END IF
     !If there were no matches, turn to new alphas
     ELSE
        WRITE(1000+myrank,'(" Calculating a total of",I6," regions found following",I6," lines")') nw,na

        IF(DEBUG) THEN
           WRITE(2900+myrank,'(" Calculating a total of",I6," regions out of a maximum of",I6)') &
                & nw,nwmax
           maxBb=MAX(MAXVAL(B1),MAXVAL(B2))
           WRITE(2900+myrank,'(" Covering most of region   ",1pe13.5," < 1/lambda <",1pe13.5)') &
                & MINVAL(Bb),maxBt
           WRITE(2900+myrank,'(" Neglecting at least region",1pe13.5," < 1/lambda <",1pe13.5,   &
                & " (",1pe9.2," % of the lambda space)")')                            &
                & maxBb,maxBt,(1./maxBt-1./maxBb)/(1./MINVAL(Bb(nw0:nw))-1./maxBb)*100
        END IF
        IF(nw.GT.nwx) THEN
           serr="nw<nwx"
           CALL END_ALL(serr,.FALSE.)
        END IF
        !Double the number of alphas and start over
        iturn=0  
!        IF(SATAKE) THEN
!           dalphap=siota*TWOPI/na
!           nalpha=na
!        ELSE
        dalphap=iota*TWOPI/nzperiod/na  !dalphap has sign
        nalpha=INT(TWOPI/ABS(dalphap)+1)
!        END IF
        na=na*2  
     ENDIF

  END DO

!!$  DO lw=1,10
!!$     DO iw=1,nw
!!$        DO jw=1,nw
!!$           DO kw=1,nw              
!!$              IF(connected(kw,jw).AND.connected(jw,iw)) connected(kw,iw)=.TRUE.
!!$           END DO
!!$        END DO
!!$     END DO
!!$  END DO

  na=na/2

!  !If centered alpha derivatives, force to have odd number field lines
!  !so that periodicity couples odd and even nodes
!  IF(CENTERED_ALPHA.AND.MOD(nalpha,2).EQ.0) nalpha=nalpha-1 

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE CHARACTERIZE_WELLS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CREATE_ANGULAR_GRID(na,nalpha,nalphab,alphap,dalphap,offset,&
     & thetap,zetap,zetax,B_al,vds_al)

!-----------------------------------------------------------------------------------------------
!From na independent field lines, create:
!- a nalphaxnalpha (zetap,thetap) grid aligned with the magnetic field lines 
! and labelled with alphap (dalphap is the spacing);
!- a nax(nalphaxnalpha/na) grid (zetax,thetap) extended over several field periods.
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL  
  IMPLICIT NONE
  !Input
  INTEGER na,nalpha,nalphab
  REAL*8 alphap(nalpha),dalphap,offset
  !Output
  REAL*8 thetap(nalpha,nalphab),zetap(nalphab),zetax(nalpha,nalphab)
  REAL*8 B_al(nalpha,nalphab),vds_al(Nnmp,nalpha,nalphab)
  !Others
  INTEGER ia,il
  REAL*8 dzeta,dummy
!  INTEGER imin,np1,jp1,intnumn(nax),nmj
!  REAL*8 offsetb,doff,zetat(nalpha),thetat(nalphab),temp(nalphab,nalphab)
  !Time
  CHARACTER*30, PARAMETER :: routine="CREATE_ANGULAR_GRID"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  WRITE(1000+myrank,'(" Angular grid has size ",I4," in alpha and ",I4," in l. &
       & Maximum helicity in B is ",I4)') nalpha,nalphab,INT(MAX(MAXVAL(np),MAXVAL(mp)))
  
!  IF(SATAKE) THEN
!     offtheta=0
!     offzeta  =PI/nzperiod
!     dzeta=TWOPI/nalphab
!     DO ia=1,nalpha
!        DO il=1,nalphab
!           zetap(il) =offzeta            +(il-1)*dzeta
!           IF(iota.GT.0) THEN
!              thetap(ia,il)=offtheta+ia*dalphap+(il-1)*dzeta*iota
!           ELSE
!              thetap(ia,il)=TWOPI+offtheta+ia*dalphap+(il-1)*dzeta*iota
!           END IF
!           zetax(ia,il) =zetap(il)
!           CALL FILL_BNODE(zetap(il),thetap(ia,il),dummy,B_al(ia,il),vds_al(ia,il),dummy,.FALSE.)
!           IF(DEBUG) WRITE(3000+myrank,'(5(1pe13.5),2I5)') zetax(ia,il),&
!                &  zetap(il),thetap(ia,il),B_al(ia,il),ia,il
!        END DO
!     END DO

!  ELSE
  dzeta=TWOPI/nzperiod/nalphab
  DO il=1,nalphab
     zetap(il)=(il-1)*dzeta
     DO ia=1,nalpha !note that dalphap can be negative, if iota is, but...
        thetap(ia,il)=offset+ia*dalphap+(il-1)*iota*dzeta 
        zetax(ia,il) =zetap(il)+INT((ia-1)/na)*TWOPI/nzperiod              
        CALL FILL_BNODE(zetap(il),thetap(ia,il),dummy,B_al(ia,il),vds_al(:,ia,il),.FALSE.)
        IF(DEBUG) WRITE(3000+myrank,'(5(1pe13.5),2I5)') zetap(il),thetap(ia,il),&
             &  zetax(ia,il),thetap(ia,il),B_al(ia,il),ia,il
     END DO
  END DO
!  END IF

  DO ia=1,nalpha
     alphap(ia)=thetap(ia,1)-iota*zetap(1)  !...alpha has the correct sign
  END DO

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE CREATE_ANGULAR_GRID


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CREATE_LAMBDA_GRID(nlambda,nw,Bb,Bt,&
     & lambdab_w,lambdac_w,dlambdap,lambda,one_o_lambda)

!-----------------------------------------------------------------------------------------------
!Create a lambda and one_lambda grids with nlambda elements and dlambdap spacing
!from the values of the bottoms Bb and tops Bt of nw wells
!(the extremes in lambda, lambdac_w and lambdab_w for each well are calculated)
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL  
  IMPLICIT NONE
  !Input
  INTEGER nlambda,nw
  REAL*8 Bb(nw),Bt(nw)
  !Output
  REAL*8 lambdab_w(nw),lambdac_w(nw),dlambdap,lambda(nlambda),one_o_lambda(nlambda)
  !Others
  LOGICAL passing
  INTEGER ila,iw
  REAL*8 lambdac,lambdab
  !Time
  CHARACTER*30, PARAMETER :: routine="CREATE_LAMBDA_GRID"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  lambdab_w=1./Bb
  lambdab=MAXVAL(lambdab_w)
  lambdac_w=1./Bt
  lambdac=MINVAL(lambdac_w)


!  lambdac=1./2.795
  dlambdap=(lambdab-lambdac)/(nlambda+1)
  
  DO ila=1,nlambda
     lambda(ila)=lambdac+ila*dlambdap
  END DO
  one_o_lambda=1./lambda

  DO ila=1,nlambda
     passing=.TRUE.
     DO iw=1,nw
        IF(lambda(ila).GT.1/Bt(iw)) passing=.FALSE.
     END DO
     IF(passing) THEN
        WRITE(1000+myrank,*) 'Passing trajectories at ila=',ila
     END IF
  END DO

  WRITE(1000+myrank,'(" Velocity grid has size ",I4," in lambda")') nlambda

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
  
END SUBROUTINE CREATE_LAMBDA_GRID


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE LABEL_GRIDPOINTS(na,nalpha,nalphab,nlambda,nw,bottom,connected,&
     & alphap_w,z1,zb,z2,Bb,Bt,lambda,thetap,zetap,zetax,B_al,npoint,i_a,i_l,i_w,i_p)

!-----------------------------------------------------------------------------------------------
!For a nalpha x nlambda grid, and nw wells characterized by the matrix of connections 
!connected, alphap_w, z1, z2, Bt Bb, determine number of points npoint, alpha, lambda and well
!of each of them, i_a, i_l and i_w, and determine point number. i_p for each point of the extended
!grid for each value of lambda
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL  
  IMPLICIT NONE
  !Input
  LOGICAL bottom(nw),connected(nw,nw)
  INTEGER na,nalpha,nalphab,nlambda,nw
  REAL*8 alphap_w(nw),z1(nw),zb(nw),z2(nw),Bb(nw),Bt(nw),lambda(nlambda)
  REAL*8 zetap(nalphab),zetax(nalpha,nalphab)
  REAL*8 thetap(nalpha,nalphab),B_al(nalpha,nalphab)
  !Output
  INTEGER npoint,i_a(npointx),i_l(npointx),i_w(npointx),i_p(nlambda,nalpha,nalphab)
  !Others
  CHARACTER*100 serr
  INTEGER ja,ia,il,ila,iw,jw,assigned(nw),iwell(nalpha,nalphab,nlambda)
  REAL*8 one_o_lambda(nlambda),alpp,lambdab_w(nw),lambdac_w(nw),d1,d2,dist,distt
  REAL*8 MODANG
  !Time
  CHARACTER*30, PARAMETER :: routine="LABEL_GRIDPOINTS"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  !Ignore wells that are outside region of interest
  DO iw=1,nw
     IF(.NOT.bottom(iw)) CYCLE
     DO ia=1,na
        alpp=(thetap(ia,1)-iota*zetax(ia,1))
        IF(ABS(alphap_w(iw)-alpp).LT.PREC_EXTR.AND.zb(iw).LT.zetax(ia,1)) THEN
           Bb(iw)=1E10
           Bt(iw)=-1
           DO jw=iw+1,nw 
              IF(connected(jw,iw)) THEN
                 Bb(jw)=1E10
                 Bt(jw)=-1
              END IF
           END DO
        END IF
     END DO
     DO ia=nalpha-na+1,nalpha
        alpp=(thetap(ia,nalphab)-iota*zetax(ia,nalphab))
        IF(ABS(alphap_w(iw)-alpp).LT.PREC_EXTR.AND.zb(iw).GT.zetax(ia,nalphab)) THEN
           Bb(iw)=1E10
           Bt(iw)=-1
           DO jw=iw+1,nw 
              IF(connected(jw,iw)) THEN
                 Bb(jw)=1E10
                 Bt(jw)=-1
              END IF
           END DO
        END IF
     END DO
  END DO

  !Size of global grid
  one_o_lambda=1/lambda
  lambdab_w=1/Bb
  lambdac_w=1/Bt

  !For each point in the (zeta,theta) extended grid, determine point

  !Start by determine well label for each point
  iwell=0
  DO ia=1,nalpha
     DO il=1,nalphab
        alpp=(thetap(ia,il)-iota*zetax(ia,il))
        DO iw=1,nw
           IF(bottom(iw).AND.ABS(alphap_w(iw)-alpp)&
                & .LT.PREC_EXTR.AND.zetax(ia,il).GT.z1(iw).AND.zetax(ia,il).LE.z2(iw)) EXIT
        END DO
        DO jw=iw,nw
           IF(iw.EQ.jw.OR.connected(jw,iw)) THEN
              DO ila=1,nlambda
                 IF(MAX(bb(jw),B_al(ia,il)).LT.one_o_lambda(ila)&
                      & .AND.one_o_lambda(ila).LT.bt(jw)) iwell(ia,il,ila)=jw
              END DO

           END IF
        END DO
     END DO
  END DO

  npoint=0
  i_p=0
  i_l=0
  i_a=0
  i_w=0

  DO ila=1,nlambda  !CHECK leave out exactly "lambda_c"?
     assigned=0
     DO ia=1,nalpha
        DO il=1,nalphab
           iw=iwell(ia,il,ila) !CHECK mimimum lambda-size of region?
           IF(iw.EQ.0) CYCLE
           IF(assigned(iw).GT.0) THEN
              i_p(ila,ia,il)=assigned(iw)
           ELSE
              npoint=npoint+1
              i_p(ila,ia,il)=npoint
              i_a(npoint)   =ia
              i_l(npoint)   =ila
              i_w(npoint)   =iw
              assigned(iw)=npoint
           END IF
        END DO
     END DO
  END DO
  npoint=npoint+1    !ipoint=npoint is g=0
  IF(npoint.GT.npointx) THEN
     serr="npoint>npointx"
     CALL END_ALL(serr,.FALSE.)
  END IF
  
  !Contour conditions: for the points ignored in the first loop of the routine,
  !use periodicity to find equivalent points that have been labelled 
  DO ia=1,nalpha
     DO il=1,nalphab
        alpp=(thetap(ia,il)-iota*zetax(ia,il))
        DO iw=1,nw
           IF(bottom(iw).AND.ABS(alphap_w(iw)-alpp).LT.PREC_EXTR.AND.&
                & ((zetax(ia,il).GT.z1(iw)))) EXIT
        END DO
        IF(iw.EQ.nw+1) THEN           
           dist=1E5
           DO ja=1,nalpha
              IF(ja.LT.nalpha-na-1) CYCLE
              d1=MODANG( zetap(il)       -zetap(nalphab),TWOPI/nzperiod)
              d2=MODANG(thetap(ia,il)-thetap(ja,nalphab),TWOPI)
              distt=d1*d1+d2*d2
              IF(distt.LT.dist) THEN
                 DO ila=1,nlambda
                    IF(one_o_lambda(ila).GT.B_al(ia,il)) THEN
                       dist=distt
                       i_p(  ila,ia,il)=i_p(  ila,ja,nalphab)
                       iwell(ia,il,ila)=iwell(ja,nalphab,ila)
                    END IF
                 END DO
              END IF
           END DO
        END IF
        
        DO iw=1,nw
           IF(bottom(iw).AND.ABS(alphap_w(iw)-alpp).LT.PREC_EXTR.AND.&
                & ((zetax(ia,il).LT.z2(iw)))) EXIT
        END DO
        IF(iw.EQ.nw+1) THEN           
           dist=1E5
           DO ja=1,nalpha
              IF(ja.GT.na+1) CYCLE 
              d1=MODANG( zetap(il) -zetap(1),TWOPI/nzperiod)
              d2=MODANG(thetap(ia,il)-thetap(ja,1),TWOPI)
              distt=d1*d1+d2*d2
              IF(distt.LT.dist) THEN
                 DO ila=1,nlambda
                    IF(one_o_lambda(ila).GT.B_al(ia,il)) THEN
                       dist=distt
                       i_p(  ila,ia,il)=i_p(  ila,ja,1)
                       iwell(ia,il,ila)=iwell(ja,1,ila)
                    END IF
                 END DO
              END IF
           END DO
        END IF

     END DO
  END DO
  
  !i_p=npoint elsewhere, where g(npoint) is going to be 0
  DO ia=1,nalpha
     DO il=1,nalphab
        DO ila=1,nlambda
           IF(i_p(ila,ia,il).EQ.0) i_p(ila,ia,il)=npoint
        END DO
     END DO
  END DO

  WRITE(1000+myrank,'(" Global grid ",I6," points")') npoint

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE LABEL_GRIDPOINTS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 SUBROUTINE FIND_LAMBDA_NEIGHBOURS(npoint,nalpha,nalphab,nlambda,nw,nbifx,i_p,bottom,i_w,&
        & i_p_lm1,nbif,i_p_lp1)

!----------------------------------------------------------------------------------------------- 
!For each point i_p of npoint, find in nalpha x nlambda grid (with nbifx as maximum number of 
!inmediate neighbours at larger lambda), the neighbours in lambda i_p_lm1, nbif, i_p_lp1 are found
!using that there are nw wells characterized by i_w and bottom
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL  
  IMPLICIT NONE
  !Input
  LOGICAL bottom(nw)
  INTEGER npoint,nalpha,nalphab,nlambda,nw,nbifx,i_p(nlambda,nalpha,nalphab),i_w(npoint)
  !Output
  INTEGER i_p_lm1(npoint),i_p_lp1(nbifx,npoint),nbif(npoint)
  !Others
  INTEGER ia,il,ila,ipoint,jpoint
  !Time
  CHARACTER*30, PARAMETER :: routine="FIND_LAMBDA_NEIGHBOURS"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  !Default values are always npoint
  i_p_lm1=npoint
  !Decreasing lambda: one neighbour, since no bifurcations
  DO ia=1,nalpha
     DO il=1,nalphab
        DO ila=1,nlambda
           ipoint=i_p(ila,ia,il)
!           IF(ipoint.EQ.npoint) CYCLE           
!           IF(ila.GT.1) i_p_lm1(ipoint)=MIN(i_p(ila-1,ia,il),i_p_lm1(ipoint))
           IF(ipoint.NE.npoint.AND.ila.GT.1) i_p_lm1(ipoint)=MIN(i_p(ila-1,ia,il),i_p_lm1(ipoint))
        END DO 
     END DO
  END DO
  !Increasing lambda: several neighours are possible in a bifurcation
  nbif=0
  i_p_lp1=npoint
  DO ipoint=1,npoint-1
     DO jpoint=1,npoint-1
        IF (i_p_lm1(jpoint).EQ.ipoint) THEN
           nbif(ipoint)=nbif(ipoint)+1
           i_p_lp1(nbif(ipoint),ipoint)=jpoint
        END IF
     END DO

     IF(i_p_lp1(1,ipoint).EQ.npoint.AND..NOT.bottom(i_w(ipoint))) THEN
        bottom(i_w(ipoint))=.TRUE.
!        i_l(ipoint)=0
!        DO ia=1,nalpha
!           DO il=1,nalphab
!              DO ila=1,nlambda
!                 IF(i_p(ila,ia,il).EQ.ipoint) i_p(ila,ia,il)=npoint
!              END DO
!           END DO
!        END DO
!        DO jpoint=1,npoint
!           IF(i_p_lm1(jpoint).EQ.ipoint) i_p_lm1(jpoint)=npoint
!        END DO        
     END IF
  END DO

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE FIND_LAMBDA_NEIGHBOURS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE FIND_ALPHA_NEIGHBOURS(npoint,i_w,i_a,i_l,i_p_lm1,zlw,tlw,zrw,trw,&
     & nw,z1,zb,z2,tb,Bb,Bt,alphap_w,bottom,connected,              &
     & nlambda,one_o_lambda,na,nalpha,alphap,dalphap,offset,        &
     & i_p_am1,i_p_ap1,dalpha_am1,dalpha_ap1,&
     & i_p_am2,i_p_ap2,dalpha_am2,dalpha_ap2)
  
!----------------------------------------------------------------------------------------------- 
!For each point of npoint, located in nalpha x nlambda (alphap,one_o_lambda) grid (characterized 
!by dalphap, and offset), i_p_am1, i_p_am2, i_p_ap1, i_p_ap2 are found
!using that there are nw wells characterized by i_w, i_a, i_l, zlw, tlw, zrw, trw, z1, zb, z2, tb,
!Bb, Bt, alphap_w, bottom and connected
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL  
  IMPLICIT NONE
  !Input
  INTEGER npoint,i_w(npoint),i_a(npoint),i_l(npoint),i_p_lm1(npoint),nw,na,nalpha,nlambda
  LOGICAL bottom(nw),connected(nw,nw)
  REAL*8 z1(nw),zb(nw),z2(nw),tb(nw),Bb(nw),Bt(nw),one_o_lambda(nlambda)
  REAL*8 zlw(npoint),tlw(npoint),zrw(npoint),trw(npoint)
  REAL*8 alphap_w(nw),alphap(nalpha),dalphap,offset
  !Output
  INTEGER i_p_am1(npoint),i_p_ap1(npoint),i_p_am2(npoint),i_p_ap2(npoint)
  REAL*8 dalpha_am1(npoint),dalpha_ap1(npoint),dalpha_am2(npoint),dalpha_ap2(npoint)
  !Others
  INTEGER, PARAMETER :: nrep=30
  INTEGER ipoint,jpoint,iw,jw,kw,ila,jla,ial,irep,id,ka,i_wp(nlambda),ntot,nvec
  REAL*8 sign,da,alpp,zp,tp,zx,tx,zrep,trep
  REAL*8 MODANGHEL
  !Time
  CHARACTER*30, PARAMETER :: routine="FIND_ALPHA_NEIGHBOURS"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  sign=dalphap/ABS(dalphap)
  
  !For each point in the (zeta,theta) grid, determine two neighbours in alpha
  i_p_ap1=npoint
  i_p_am1=npoint
  dalpha_ap1=ABS(iota*TWOPI/nzperiod)/na
  dalpha_am1=dalpha_ap1
  DO ipoint=1,npoint-1

     !Determine representative (alpha,l) of a field line labelled by l
     iw =i_w(ipoint)
     ila=i_l(ipoint)
     ial=i_a(ipoint)

     IF(ila.EQ.0) CYCLE

     DO irep=0,nrep
        IF(i_p_am1(ipoint).NE.npoint.AND.i_p_ap1(ipoint).NE.npoint) CYCLE
        IF(irep.EQ.0) THEN
           trep=tb(iw)
           zrep=zb(iw)
        ELSE
           zrep=zlw(ipoint)+(zrw(ipoint)-zlw(ipoint))*(irep-1.)/(nrep-1.)
           trep=tlw(ipoint)+(trw(ipoint)-tlw(ipoint))*(irep-1.)/(nrep-1.)
        END IF
        zrep=MODANGHEL(zrep,offset,TWOPI/nzperiod)
        alpp=(trep-iota*zrep)
        IF((alpp-alphap(1).GT.PREC_EXTR.AND.alpp-alphap(nalpha).GT.PREC_EXTR).OR. &
             &  (alphap(1)-alpp.GT.PREC_EXTR.AND.alphap(nalpha)-alpp.GT.PREC_EXTR)) CYCLE
        
        !Loop in direction in alpha (i.e. alpha+dalpha and alpha-dalpha)
        DO id=-1,1,2
           IF((id.LT.0.AND.i_p_am1(ipoint).NE.npoint).OR.&
             &(id.GT.0.AND.i_p_ap1(ipoint).NE.npoint)) CYCLE
           da=ABS(iota*TWOPI/nzperiod)/na 
           tp=trep+id*da     /(1.+iota2)   ! the factor is cos(atan(iota))*d_alpha
           zp=zrep-id*da*iota/(1.+iota2)   ! the factor is sin(atan(iota))*d_alpha
           zp=MODANGHEL(zp,offset,TWOPI/nzperiod)
           alpp=tp-iota*zp
           !Ensure periodicity            
           IF((alpp-alphap(1).GT.PREC_EXTR.AND.alpp-alphap(nalpha).GT.PREC_EXTR).OR. &
                & (alphap(1)-alpp.GT.PREC_EXTR.AND.alphap(nalpha)-alpp.GT.PREC_EXTR)) THEN
              da=TWOPI-da*INT(TWOPI/da)
              tp=trep+id*da     /(1.+iota2)
              zp=zrep-id*da*iota/(1.+iota2)
              zp=MODANGHEL(zp,offset,TWOPI/nzperiod) 
              IF(    (iota.LT.0.AND.     alpp-alphap(1).GT.PREC_EXTR.AND.id.LT.0).OR.&
                   & (iota.LT.0.AND.alphap(nalpha)-alpp.GT.PREC_EXTR.AND.id.GT.0).OR.&
                   & (iota.GT.0.AND.     alpp-alphap(1).LT.PREC_EXTR.AND.id.GT.0).OR.&
                   & (iota.GT.0.AND.alphap(nalpha)-alpp.LT.PREC_EXTR.AND.id.LT.0)) THEN              
                 IF(sign.LT.0) THEN
                    IF(alpp-alphap(1)  .GT.PREC_EXTR)    tp=tp-TWOPI
                    IF(alphap(nalpha)-alpp.GT.PREC_EXTR) tp=tp+TWOPI
                 ELSE
                    IF(alphap(1)-alpp  .GT.PREC_EXTR)    tp=tp+TWOPI
                    IF(alpp-alphap(nalpha).GT.PREC_EXTR) tp=tp-TWOPI
                 END IF
              ELSE
                 IF(sign.LT.0) THEN
                    IF(alpp-alphap(1)  .GT.PREC_EXTR) tp=tp-TWOPI
                    IF(alphap(nalpha)-alpp.GT.PREC_EXTR) tp=tp+TWOPI
                 ELSE
                    IF(alphap(1)-alpp  .GT.PREC_EXTR) tp=tp+TWOPI
                    IF(alpp-alphap(nalpha).GT.PREC_EXTR) tp=tp-TWOPI
                 END IF
              END IF
              zp=MODANGHEL(zp,offset,TWOPI/nzperiod)
           END IF
           zp=MODANGHEL(zp,offset,TWOPI/nzperiod)
           ka=INT(((tp-iota*zp)-alphap(1)+sign*PREC_EXTR)/dalphap)
           zx=zp+INT(ka/na)*TWOPI/nzperiod
           tx=tp
           i_wp=-1

           alpp=(tx-iota*zx)
           DO kw=1,nw
              IF(.NOT.bottom(kw).OR.ABS(alphap_w(kw)-alpp).GT.PREC_EXTR &
                   & .OR.zx.LT.z1(kw).OR.zx.GE.z2(kw)) CYCLE
              EXIT
           END DO
           DO jw=kw,nw
              IF(kw.EQ.jw.OR.connected(jw,kw)) THEN
                 DO jla=1,nlambda
                    IF(bb(jw).LT.one_o_lambda(jla).AND.one_o_lambda(jla).LT.bt(jw)) i_wp(jla)=jw
                 END DO
              END IF
           END DO
           DO jpoint=1,npoint-1
              IF(i_l(jpoint).EQ.ila.AND.i_w(jpoint).EQ.i_wp(ila)) THEN
                 EXIT
              END IF
           END DO  
           IF(id.LT.0) THEN
              dalpha_am1(ipoint)=da
              i_p_am1(ipoint)=jpoint
           ELSE
              dalpha_ap1(ipoint)=da
              i_p_ap1(ipoint)=jpoint
           END IF
        END DO
     END DO
  END DO
  
  !Make sure each neighbours are consistent
  DO ipoint=1,npoint-1
     IF(i_p_ap1(ipoint).EQ.npoint) THEN
        DO jpoint=1,npoint-1
           IF(i_p_am1(jpoint).EQ.ipoint) THEN
              i_p_ap1(ipoint)=jpoint
              dalpha_ap1(ipoint)=dalpha_am1(jpoint)
           END IF
        END DO
     END IF
     IF(i_p_am1(ipoint).EQ.npoint) THEN
        DO jpoint=1,npoint-1
           IF(i_p_ap1(jpoint).EQ.ipoint) THEN
              i_p_am1(ipoint)=jpoint
              dalpha_am1(ipoint)=dalpha_ap1(jpoint)
           END IF
        END DO
     END IF
  END DO


  !When trajectories do not close in alpha, there may be problems, because no neighbour exist
  !at same value of lambda. Here, two different models are implemented
  
  IF(.NOT.CLOSEST_LAMBDA) THEN !Find largest lambda with closed trajectories in alpha
     jla=0
     DO ila=1,nlambda
        ntot=0
        nvec=0
        DO ipoint=1,npoint-1
           IF(i_l(ipoint).NE.ila) CYCLE
           ntot=ntot+1
           IF(i_p_ap1(ipoint).NE.npoint) nvec=nvec+1
        END DO
        IF(nvec.EQ.ntot.AND.ila.GT.jla) jla=ila
     END DO
  END IF

  DO ipoint=1,npoint-1
     IF(i_p_ap1(ipoint).EQ.npoint) THEN
        jpoint=ipoint
        DO jla=i_l(ipoint),1,-1
           jpoint=i_p_lm1(jpoint)
           IF(jpoint.EQ.npoint) EXIT
           IF((CLOSEST_LAMBDA.AND.i_p_ap1(jpoint).NE.npoint).OR.&
                (.NOT.CLOSEST_LAMBDA.AND.i_l(jpoint).EQ.ila)) THEN
              i_p_ap1(ipoint)   =i_p_ap1(jpoint)
              dalpha_ap1(ipoint)=dalpha_ap1(jpoint)
              EXIT
           END IF
        END DO
     END IF
     IF(i_p_am1(ipoint).EQ.npoint) THEN
        jpoint=ipoint
        DO jla=i_l(ipoint),1,-1
           jpoint=i_p_lm1(jpoint)
           IF(jpoint.EQ.npoint) EXIT
           IF((CLOSEST_LAMBDA.AND.i_p_am1(jpoint).NE.npoint).OR.&
                (.NOT.CLOSEST_LAMBDA.AND.i_l(jpoint).EQ.ila)) THEN
              i_p_am1(ipoint)   =i_p_am1(jpoint)
              dalpha_am1(ipoint)=dalpha_am1(jpoint)
              EXIT
           END IF
        END DO
     END IF
  END DO
  
  !Seconds neighbours
  i_p_ap2=npoint
  i_p_am2=npoint
  dalpha_ap2=dalpha_ap1
  dalpha_am2=dalpha_ap1
  DO ipoint=1,npoint-1
     IF(i_l(i_p_ap1(i_p_ap1(ipoint))).EQ.i_l(ipoint).OR..NOT.CLOSEST_LAMBDA1) THEN
        i_p_ap2(ipoint)=i_p_ap1(i_p_ap1(ipoint))
        dalpha_ap2(ipoint)=dalpha_ap1(i_p_ap1(ipoint))
     END IF
     IF(i_l(i_p_am1(i_p_am1(ipoint))).EQ.i_l(ipoint).OR..NOT.CLOSEST_LAMBDA1) THEN
        i_p_am2(ipoint)=i_p_am1(i_p_am1(ipoint))
        dalpha_am2(ipoint)=dalpha_am1(i_p_am1(ipoint))
     END IF
  END DO

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)


END SUBROUTINE FIND_ALPHA_NEIGHBOURS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE COEFFICIENTS_DKE(npoint,i_w,i_l,nw,z1,t1,B1,hBpp1,vd1,&
                         & zb,tb,Bb,hBppb,vdb,z2,t2,B2,hBpp2,vd2,&
                         & nlambda,lambda,zlw,tlw,zrw,trw, &
                         & BI1,BI2,BI3,BI4,BI5,Bi6,nmodes,BI7)
!----------------------------------------------------------------------------------------------- 
!For npoints characterized by i_w,i_l, with bounce points zlw, tlw, zrw and trw, located in nw 
!wells characterized by z,t,B,hBpp and vd at extremes 1, 2 and b, and in a lambda grid of size
!nlambda, calculate bounce integrals BI1, BI2, BI3, BI4, BI5, BI6 and, BI7, the latter with size
! nmodes
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL  
  IMPLICIT NONE
  !Input
  INTEGER npoint,i_w(npoint),i_l(npoint),nw,nlambda,nmodes
  REAL*8 z1(nw),t1(nw),B1(nw),hBpp1(nw),vd1(nqv,nw)
  REAL*8 zb(nw),tb(nw),Bb(nw),hBppb(nw),vdb(nqv,nw) 
  REAL*8 z2(nw),t2(nw),B2(nw),hBpp2(nw),vd2(nqv,nw) 
  REAL*8 lambda(nlambda),zlw(npoint),tlw(npoint),zrw(npoint),trw(npoint)
  !Output
  REAL*8 BI1(npoint),BI2(npoint),BI3(npoint),BI4(npoint),BI5(npoint),BI6(npoint)
  REAL*8 BI7(npoint,nmodes)
  !Others
  INTEGER ipoint,iw,ila,nq
  REAL*8, ALLOCATABLE :: Q(:)
  !Time
  CHARACTER*30, PARAMETER :: routine="COEFFICIENTS_DKE"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart
!  REAL*8 dJdpsi(nlambda)

  CALL CPU_TIME(tstart)

  nq=nq0
  IF(SOLVE_QN) nq=nq+nmodes
  ALLOCATE(Q(nq))

!  dJdpsi=0

  DO ipoint=1,npoint
     iw =i_w(ipoint)
     ila=i_l(ipoint)
     IF(ila.EQ.0.OR.ipoint.EQ.npoint) THEN
        BI1(ipoint)=0
        BI2(ipoint)=0
        BI3(ipoint)=0
        BI4(ipoint)=0
        BI5(ipoint)=0
        BI6(ipoint)=0
        IF(SOLVE_QN) BI7(ipoint,1:nmodes)=0
     ELSE
        CALL BOUNCES(iw,z1(iw),t1(iw),B1(iw),hBpp1(iw),vd1(:,iw), &
             &          zb(iw),tb(iw),Bb(iw),hBppb(iw),vdb(:,iw), &
             &          z2(iw),t2(iw),B2(iw),hBpp2(iw),vd2(:,iw), &
             & 1./lambda(ila),.FALSE.,nq,Q,&
             & zlw(ipoint),tlw(ipoint),zrw(ipoint),trw(ipoint))    
        BI1(ipoint)=Q(1)
        BI2(ipoint)=Q(2)
        BI3(ipoint)=Q(3)
        BI4(ipoint)=Q(4)
        BI5(ipoint)=Q(5)  
        BI6(ipoint)=Q(6)  
        IF(SOLVE_QN) BI7(ipoint,1:nmodes)=Q(7:nq)           
     END IF
     IF(DEBUG.AND.(iw.EQ.I0.OR.I0.EQ.0)) WRITE(3100+myrank,'(3I6,7(1pe13.5),I4,4(1pe13.5))') ipoint,ila,iw,&
          &  BI1(ipoint),BI2(ipoint),BI3(ipoint),&
          &  BI4(ipoint),BI5(ipoint),BI6(ipoint),BI7(ipoint,1),nlambda,zlw(ipoint),tlw(ipoint),zrw(ipoint),trw(ipoint)
!     dJdpsi(ila)=dJdpsi(ila)+BI4(ipoint)*ABS(BI2(ipoint))
!     WRITE(1000+myrank,'("dJdpsi_al ",1I6,7(1pe13.5))') ila,tlw(ipoint),BI4(ipoint),BI1(ipoint)
  END DO
!  stop
!  DO ila=1,nlambda
!     WRITE(1000,*) 'dJdpsi',ila,dJdpsi(ila)
!  END DO
!  stop
  
  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE COEFFICIENTS_DKE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef MPIandPETSc

SUBROUTINE INIT_LINEAR_PROBLEM(npoint,nnz,matCOL,matVEAf,matVEAb,matVMAf,matVMAb)
 
!----------------------------------------------------------------------------------------------- 
!Initialize linear problem of size npoint and nnz with PETSc: ksp and matrices matCOL, matVEAf
!matVEAb, matVMAf, matVMAb
!-----------------------------------------------------------------------------------------------    
                    
  USE GLOBAL  
  USE petscsys
  USE petscksp
  IMPLICIT NONE
!#include <petsc/finclude/petscksp.h>
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!#include <petsc/finclude/petscpc.h>
  !Input
  INTEGER npoint,nnz(npoint)
  !Output
!  KSP ksp
  Mat matCOL,matVEAf,matVEAb,matVMAf,matVMAb!,matA
  !Others
  PetscErrorCode ierr
  PetscInt iz,innz(0:npoint-1)
  INTEGER npalpha
!  INTEGER, PARAMETER :: MAXITS= 1000 !these choices could go in the parameters namelist
!  REAL, PARAMETER :: ATOL=1E-2
!  REAL, PARAMETER :: TOL =1E-1
!  PetscReal, PARAMETER :: factor=10.
!  PetscInt       Istart,Iend
!  PC  pc
  !Time
  CHARACTER*30, PARAMETER :: routine="INIT_LINEAR_PROBLEM"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  !Number of non-zero elements per row
  DO iz=1,npoint
     innz(iz-1)=nnz(iz)
  END DO
  !Collision operator
  CALL MatCreate(PETSC_COMM_WORLD,matCOL,ierr)
  CALL MatSetSizes(matCOL,PETSC_DECIDE,PETSC_DECIDE,npoint,npoint,ierr)
  CALL MatSetType(matCOL,MATAIJ,ierr)
  !  CALL MatSeqAIJSetPreallocation(matCOL,PETSC_NULL_INTEGER,innz,ierr)
  CALL MatSeqAIJSetPreallocation(matCOL,0,innz,ierr)
  CALL MatSetup(matCOL,ierr)

  IF(CENTERED_ALPHA) THEN
     npalpha=4
  ELSE
     npalpha=3
  END IF

  !Tangential ExB drift, forward derivative
  CALL MatCreate(PETSC_COMM_WORLD,matVEAf,ierr)
  CALL MatSetSizes(matVEAf,PETSC_DECIDE,PETSC_DECIDE,npoint,npoint,ierr)
  CALL MatSetType(matVEAf,MATAIJ,ierr)
  CALL MatSeqAIJSetPreallocation(matVEAf,npalpha,PETSC_NULL_INTEGER,ierr)
  CALL MatSetup(matVEAf,ierr)
  !Tangential ExB drift, backward derivative
  CALL MatCreate(PETSC_COMM_WORLD,matVEAb,ierr)
  CALL MatSetSizes(matVEAb,PETSC_DECIDE,PETSC_DECIDE,npoint,npoint,ierr)
  CALL MatSetType(matVEAb,MATAIJ,ierr)
  CALL MatSeqAIJSetPreallocation(matVEAb,npalpha,PETSC_NULL_INTEGER,ierr)
  CALL MatSetup(matVEAb,ierr)
  
  IF(TANG_VM) THEN
     !Tangential magnetic drift, forward derivative
     CALL MatCreate(PETSC_COMM_WORLD,matVMAf,ierr)
     CALL MatSetSizes(matVMAf,PETSC_DECIDE,PETSC_DECIDE,npoint,npoint,ierr)
     CALL MatSetType(matVMAf,MATAIJ,ierr)
     CALL MatSeqAIJSetPreallocation(matVMAf,npalpha,PETSC_NULL_INTEGER,ierr)
     CALL MatSetup(matVMAf,ierr)
     !Tangential magnetic drift, backward derivative
     CALL MatCreate(PETSC_COMM_WORLD,matVMAb,ierr)
     CALL MatSetSizes(matVMAb,PETSC_DECIDE,PETSC_DECIDE,npoint,npoint,ierr)
     CALL MatSetType(matVMAb,MATAIJ,ierr)
     CALL MatSeqAIJSetPreallocation(matVMAb,npalpha,PETSC_NULL_INTEGER,ierr)
     CALL MatSetup(matVMAb,ierr)
  END IF
  
!!$  CALL KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
!!$  CALL KSPSetTolerances(ksp,TOL,ATOL,PETSC_DEFAULT_REAL,MAXITS,ierr)
!!$!  CALL KSPSetType(ksp,KSPPREONLY,ierr)
!!$!  CALL KSPGetPC(ksp,pc,ierr)
!!$!  CALL PCSetType(pc,PCILU,ierr)
!!$!  CALL PCFactorSetLevels(pc,10,ierr)
!!$!  CALL PCFactorSetMatOrderingType(pc,"natural",ierr)
!!$
!!$  IF(DIRECT_SOLVER) THEN
!!$     CALL KSPSetType(ksp,KSPPREONLY,ierr)
!!$     CALL KSPGetPC(ksp,pc,ierr)
!!$     CALL PCSetType(pc,PCLU,ierr)
!!$     CALL PCFactorSetMatOrderingType(pc,"natural",ierr)
!!$!     CALL PCFactorSetMatOrderingType(pc,"nd",ierr)
!!$!     CALL PCFactorSetReuseOrdering(pc,PETSC_TRUE,ierr)
!!$!     CALL PCFactorSetReuseFill(pc,PETSC_TRUE,ierr)
!!$!#ifdef PETSC_HAVE_MUMPS
!!$!     call PCFactorSetMatSolverType(pc,MATSOLVERMUMPS,ierr)
!!$!     call PCFactorSetUpMatSolverType(pc,ierr)
!!$!#endif
!!$  ELSE
!!$     CALL KSPSetType(ksp,KSPGMRES,ierr)
!!$     CALL KSPGetPC(ksp,pc,ierr)     
!!$!     CALL PCSetType(pc,PCJACOBI,ierr)
!!$     CALL PCSetType(pc,PCILU,ierr)
!!$     CALL PCFactorSetLevels(pc,2,ierr)
!!$     CALL PCFactorSetMatOrderingType(pc,"natural",ierr)
!!$  END IF

!From more to less stable
!  CALL PCFactorSetMatOrderingType(pc,"rcm",ierr)
!  CALL PCFactorSetMatOrderingType(pc,"1wd",ierr)
!  CALL PCFactorSetMatOrderingType(pc,"nd",ierr) !default
!  CALL PCFactorSetMatOrderingType(pc,"qmd",ierr)
!  CALL PCFactorSetMatOrderingType(pc,"rowlength",ierr)
!  CALL PCFactorSetMatOrderingType(pc,"wbm",ierr)
!  CALL PCFactorSetMatOrderingType(pc,"spectral",ierr)
!  CALL PCFactorSetMatOrderingType(pc,"amd",ierr)
!  CALL PCFactorSetReuseOrdering(pc,.TRUE.,ierr)
!Other options
!  CALL PCFactorReorderForNonzeroDiagonal(pc,1E-3,ierr)
!  CALL PCFactorSetPivotInBlocks(pc,.TRUE.,ierr)
!  CALL PCFactorSetMatSolverType(pc,MATSOLVERMUMPS,ierr)
!  CALL PCFactorSetUpMatSolverType(pc,ierr)
!  CALL PCFactorSetFill(pc,factor,ierr)
  
  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE INIT_LINEAR_PROBLEM

#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE FILL_DKE_ROW(ipoint,npoint,&
     & dalpha_am1,dalpha_ap1,dalpha_am2,dalpha_ap2,i_p_am1,i_p_ap1,i_p_am2,i_p_ap2,&
     & i_l,lambda,dlambda_lm1,dlambda_lp1,i_p_lm1,nbif,nbifx,i_p_lp1,&
     & BI1,BI2,BI4,BI5,matCOL,matVEAf,matVEAb,matVMAf,matVMAb,nnz,flag_nnz)

!----------------------------------------------------------------------------------------------- 
!Fill row ipoint of npoint of several matrixes (matCOL, matVEAf, matVEAb, matVMAf, matVMAb), 
!with nnz non-zero elements, using spacing dlambda_lm1, dlambda_lp1, dalpha, factors
!factDKErhs, lambda, Epsi/vdconst, neighbours i_p_am1, i_p_ap1, i_p_lm1, nbif neighbours
!i_p_lp1(:) and coefficients BI[1-5]
!If flag_nnz, use this routine to calculate nnz only
!-----------------------------------------------------------------------------------------------  
                      
  USE GLOBAL  
  IMPLICIT NONE
  !Input
  LOGICAL flag_nnz
  INTEGER ipoint,npoint,nbifx
  INTEGER i_p_am1,i_p_ap1,i_p_am2,i_p_ap2,i_l(npoint),i_p_lm1(npoint),nbif(npoint),i_p_lp1(nbifx,npoint)
  REAL*8 lambda,dlambda_lm1(npoint),dlambda_lp1(npoint),dalpha_am1,dalpha_ap1,dalpha_am2,dalpha_ap2
  REAL*8 BI1,BI2(npoint),BI4,BI5
  !Output
  INTEGER nnz
  REAL*8 matCOL(npoint),matVEAf(npoint),matVEAb(npoint),matVMAf(npoint),matVMAb(npoint)
  !Others
  INTEGER ibif,jbif,j_p_lm1,j_p_lm2,j_p_lp1,j_p_lp2,jpoint
  REAL*8 denom,a,b,bm,bp,suma,sumb,d
  REAL*8 twodlambda,dlambda2,fourdlambda2
 !Time
!  CHARACTER*30, PARAMETER :: routine="FILL_DKE_ROW"
!  INTEGER, SAVE :: ntotal=0
!  REAL*8,  SAVE :: ttotal=0
!  REAL*8,  SAVE :: t0=0
!  REAL*8 tstart
!  
!  CALL CPU_TIME(tstart)

  j_p_lm1=i_p_lm1(ipoint)

  IF(j_p_lm1.EQ.npoint.AND.i_p_lp1(1,ipoint).EQ.npoint) THEN
     nnz=1
  ELSE IF(j_p_lm1.EQ.npoint) THEN
     nnz=2
  ELSE IF(nbif(ipoint).EQ.1.AND.nbif(j_p_lm1).EQ.1) THEN
     nnz=3
  ELSE
     nnz=2
     DO ibif=1,nbif(j_p_lm1)
        jpoint=i_p_lp1(ibif,j_p_lm1)
        DO jbif=1,nbif(jpoint)
           j_p_lp1=i_p_lp1(jbif,jpoint)
           j_p_lp2=i_p_lp1(1  ,j_p_lp1)
           IF(j_p_lp2.EQ.npoint) CYCLE
           nnz=nnz+1
           IF(jpoint.NE.ipoint) nnz=nnz+1
        END DO
     END DO
  END IF
  IF(nnz.NE.1) nnz=nnz+4
  IF(flag_nnz) RETURN

  !Collision operator
  twodlambda=dlambda_lm1(ipoint)+dlambda_lp1(ipoint)
  dlambda2  =dlambda_lm1(ipoint)*dlambda_lp1(ipoint)

  IF(.NOT.PHI1_READ.OR..NOT.GEN_FLAG(3)) THEN
     a=(BI2(ipoint)/lambda-BI1*lambda/2.)
     b=BI2(ipoint)
  ELSE
     a=0.5*(3*BI2(ipoint)/lambda-BI1*lambda*SQRT(lambda))
     b=BI2(ipoint)        
  END IF

  IF(j_p_lm1.EQ.npoint.AND.i_p_lp1(1,ipoint).EQ.npoint) THEN
     WRITE(1100+myrank,*) 'These points should not exist'     
     matCOL(ipoint)=1
  ELSE IF(j_p_lm1.EQ.npoint) THEN
     j_p_lp1=i_p_lp1(1,ipoint)
     denom=twodlambda*dlambda_lm1(ipoint)*dlambda_lp1(ipoint)/2
     matCOL(j_p_lp1)= b*dlambda_lm1(ipoint)/denom+ a/twodlambda
     matCOL(ipoint) =-b*twodlambda         /denom
  ELSE IF(nbif(ipoint).EQ.1.AND.nbif(j_p_lm1).EQ.1) THEN
     j_p_lp1=i_p_lp1(1,ipoint)
     suma=a/twodlambda
     sumb=b/dlambda2
     matCOL(j_p_lp1)=sumb+suma
     matCOL(ipoint) =-2*sumb
     matCOL(j_p_lm1)=sumb-suma
  ELSE
     j_p_lm2=i_p_lm1(j_p_lm1)
     bm=BI2(j_p_lm1)
     fourdlambda2=(dlambda_lp1(j_p_lm1)+dlambda_lm1(j_p_lm1))*twodlambda
     matCOL(j_p_lm2)=bm/fourdlambda2
     matCOL(ipoint)=-bm/fourdlambda2
     DO ibif=1,nbif(j_p_lm1)
        jpoint=i_p_lp1(ibif,j_p_lm1)
        DO jbif=1,nbif(jpoint)
           j_p_lp1=i_p_lp1(jbif,jpoint)
           j_p_lp2=i_p_lp1(1  ,j_p_lp1)
           IF(j_p_lp2.EQ.npoint) CYCLE
           bp=BI2(j_p_lp1)
           fourdlambda2=twodlambda*(dlambda_lp1(j_p_lp1)+dlambda_lm1(j_p_lp1))
           matCOL(j_p_lp2)=bp/fourdlambda2
           matCOL(jpoint)=matCOL(jpoint)-bp/fourdlambda2
        END DO
     END DO
  END IF

  !d_psi J d_alpha g
  IF(nnz.NE.1) THEN
     IF(INC_EXB) THEN
        d=BI5
     ELSE
        d=BI1
     END IF

!     nnz=nnz+4
     IF(CENTERED_ALPHA) THEN

        denom=7*(dalpha_ap1+dalpha_am1)-dalpha_ap2-dalpha_am2
        matVEAf(i_p_ap2)=  -d/denom
        matVEAf(i_p_ap1)=+8*d/denom
        matVEAf(i_p_am1)=-8*d/denom
        matVEAf(i_p_am2)=  +d/denom
        matVEAb(i_p_ap2)=  -d/denom
        matVEAb(i_p_ap1)=+8*d/denom
        matVEAb(i_p_am1)=-8*d/denom
        matVEAb(i_p_am2)=  +d/denom
        
     ELSE
        IF(i_p_ap2.NE.npoint.AND..NOT.FIRST_ORDER_ALPHA) THEN
           denom=dalpha_ap1*dalpha_ap2*(dalpha_ap1+dalpha_ap2)
           matVEAf(i_p_ap2)=-d*dalpha_ap1*dalpha_ap1/denom
           matVEAf(i_p_ap1)=+d*(dalpha_ap1+dalpha_ap2)*(dalpha_ap1+dalpha_ap2)/denom
           matVEAf(ipoint )=+d*(dalpha_ap1*dalpha_ap1-(dalpha_ap1+dalpha_ap2)*(dalpha_ap1+dalpha_ap2))&
                & /denom
        ELSE
           IF(i_l(i_p_ap1).EQ.i_l(ipoint)) THEN
              matVEAf(i_p_ap1)=+d/dalpha_ap1
           ELSE
              matVEAf(i_p_ap1)         =+(d/dalpha_ap1)*(1+i_l(ipoint)-i_l(i_p_ap1))
              matVEAf(i_p_lm1(i_p_ap1))=+(d/dalpha_ap1)*(i_l(i_p_ap1)-i_l(ipoint))
           END IF
           matVEAf(ipoint )=-d/dalpha_ap1
        END IF
        IF(i_p_am2.NE.npoint.AND..NOT.FIRST_ORDER_ALPHA) THEN
           denom=dalpha_am1*dalpha_am2*(dalpha_am1+dalpha_am2)
           matVEAb(i_p_am2)=+d*dalpha_am1*dalpha_am1/denom
           matVEAb(i_p_am1)=-d*(dalpha_am1+dalpha_am2)*(dalpha_am1+dalpha_am2)/denom
           matVEAb(ipoint )=-d*(dalpha_am1*dalpha_am1-(dalpha_am1+dalpha_am2)*(dalpha_am1+dalpha_am2))&
                & /denom
        ELSE
           IF(i_l(i_p_am1).EQ.i_l(ipoint)) THEN
              matVEAb(i_p_am1)=-d/dalpha_am1
           ELSE
              matVEAb(i_p_am1)         =-(d/dalpha_am1)*(1+i_l(ipoint)-i_l(i_p_am1))
              matVEAb(i_p_lm1(i_p_am1))=-(d/dalpha_am1)*(i_l(i_p_am1)-i_l(ipoint))
           END IF
           matVEAb(ipoint )=+d/dalpha_am1
        END IF
     END IF

     IF(TANG_VM) THEN

        d=-BI4

        IF(CENTERED_ALPHA) THEN

           denom=7*(dalpha_ap1+dalpha_am1)-dalpha_ap2-dalpha_am2
           matVMAf(i_p_ap2)=  -d/denom
           matVMAf(i_p_ap1)=+8*d/denom
           matVMAf(i_p_am1)=-8*d/denom
           matVMAf(i_p_am2)=  +d/denom
           matVMAb(i_p_ap2)=  -d/denom
           matVMAb(i_p_ap1)=+8*d/denom
           matVMAb(i_p_am1)=-8*d/denom
           matVMAb(i_p_am2)=  +d/denom

        ELSE

           IF(d.LT.0) THEN 
              IF(i_p_ap2.NE.npoint.AND..NOT.FIRST_ORDER_ALPHA) THEN
                 denom=dalpha_ap1*dalpha_ap2*(dalpha_ap1+dalpha_ap2) 
                 matVMAf(i_p_ap2)=-d*dalpha_ap1*dalpha_ap1/denom
                 matVMAf(i_p_ap1)=+d*(dalpha_ap1+dalpha_ap2)*(dalpha_ap1+dalpha_ap2)/denom
                 matVMAf(ipoint )=+d*(dalpha_ap1*dalpha_ap1-(dalpha_ap1+dalpha_ap2)*(dalpha_ap1+dalpha_ap2))/denom
              ELSE
                 matVMAf(i_p_ap1)=+d/dalpha_ap1
                 matVMAf(ipoint )=-d/dalpha_ap1
              END IF
              IF(i_p_am2.NE.npoint.AND..NOT.FIRST_ORDER_ALPHA) THEN
                 denom=dalpha_am1*dalpha_am2*(dalpha_am1+dalpha_am2) 
                 matVMAb(i_p_am2)=+d*dalpha_am1*dalpha_am1/denom
                 matVMAb(i_p_am1)=-d*(dalpha_am1+dalpha_am2)*(dalpha_am1+dalpha_am2)/denom
                 matVMAb(ipoint )=-d*(dalpha_am1*dalpha_am1-(dalpha_am1+dalpha_am2)*(dalpha_am1+dalpha_am2))/denom
              ELSE
                 matVMAb(i_p_am1)=-d/dalpha_am1
                 matVMAb(ipoint )=+d/dalpha_am1
              END IF
           ELSE IF(d.GT.0) THEN 
              IF(i_p_am2.NE.npoint.AND..NOT.FIRST_ORDER_ALPHA) THEN
                 denom=dalpha_am1*dalpha_am2*(dalpha_am1+dalpha_am2) 
                 matVMAf(i_p_am2)=+d*dalpha_am1*dalpha_am1/denom
                 matVMAf(i_p_am1)=-d*(dalpha_am1+dalpha_am2)*(dalpha_am1+dalpha_am2)/denom
                 matVMAf(ipoint )=-d*(dalpha_am1*dalpha_am1-(dalpha_am1+dalpha_am2)*(dalpha_am1+dalpha_am2))/denom
              ELSE
                 matVMAf(i_p_am1)=-d/dalpha_am1
                 matVMAf(ipoint )=+d/dalpha_am1
              END IF
              IF(i_p_ap2.NE.npoint.AND..NOT.FIRST_ORDER_ALPHA) THEN
                 denom=dalpha_ap1*dalpha_ap2*(dalpha_ap1+dalpha_ap2) 
                 matVMAb(i_p_ap2)=-d*dalpha_ap1*dalpha_ap1/denom
                 matVMAb(i_p_ap1)=+d*(dalpha_ap1+dalpha_ap2)*(dalpha_ap1+dalpha_ap2)/denom
                 matVMAb(ipoint )=+d*(dalpha_ap1*dalpha_ap1-(dalpha_ap1+dalpha_ap2)*(dalpha_ap1+dalpha_ap2))/denom
              ELSE
                 matVMAb(i_p_ap1)=+d/dalpha_ap1
                 matVMAb(ipoint )=-d/dalpha_ap1
              END IF

           END IF
           
        END IF
        
     END IF
  END IF
  
  IF(DEBUG) THEN
     WRITE(3150+myrank,'(10I6,5(1pe13.5),I4)') ipoint,i_p_am2,i_p_am1,i_p_ap1,i_p_ap2
     IF(ipoint.EQ.I0.OR.I0.EQ.0) THEN
        IF(ABS(matCOL(ipoint)).GT.ZERO) THEN
           DO jpoint=1,npoint
              WRITE(3200+myrank,'(2I6,5(1pe13.5),I4)') ipoint,jpoint,matCOL(jpoint),&
                   & matVEAb(jpoint),matVEAf(jpoint),matVMAb(jpoint),matVMAf(jpoint)
           END DO
        END IF
     ELSE
        WRITE(3200+myrank,'(2I6)') ipoint,nnz
     END IF
  END IF

END SUBROUTINE FILL_DKE_ROW


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef MPIandPETSc

SUBROUTINE FILL_MATRIX_PETSC(matCOL,jv,Epsi,matVEAf,matVEAb,matVMAf,matVMAb,ksp)

!-----------------------------------------------------------------------------------------------
!Use precalculated matrices matCOL, matVEAf, matVEAb, matVMAf, matVMAb and factors Epsi, nu(jv)
!and vdconst(jv) to create linear system ksp  
!-----------------------------------------------------------------------------------------------  

  USE GLOBAL
  USE petscsys
  USE petscksp
  IMPLICIT NONE
!#include <petsc/finclude/petscksp.h>
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!#include <petsc/finclude/petscpc.h>
  !Input
  INTEGER jv
  REAL*8 Epsi
  Mat matCOL,matVEAf,matVEAb,matVMAf,matVMAb
  !Output
  KSP ksp
  !Others
  PetscErrorCode ierr
  Mat matA!,matA2
  PC pc
  INTEGER, PARAMETER :: MAXITS= 1000 
  REAL*8, PARAMETER :: ATOL=1D-2
  REAL*8, PARAMETER :: TOL =1D-1
  !Time
  CHARACTER*30, PARAMETER :: routine="FILL_MATRIX"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)
!  CALL MatCreate(PETSC_COMM_WORLD,matA,ierr)
!  CALL MatSetSizes(matA,PETSC_DECIDE,PETSC_DECIDE,npoint,npoint,ierr)
!  CALL MatSetType(matA,MATAIJ,ierr)
!  CALL MatSeqAIJSetPreallocation(matA,PETSC_NULL_INTEGER,innz(0:npoint-1),ierr)
!  CALL MatSetup(matA,ierr)
  CALL MatDuplicate(matCOL,MAT_COPY_VALUES,matA,ierr)
!  CALL MatCOPY(matCOL,matA,ierr)
!  CALL MatAssemblyBegin(matA,MAT_FINAL_ASSEMBLY,ierr)
!  CALL MatAssemblyEnd(  matA,MAT_FINAL_ASSEMBLY,ierr)

  IF(sgnB*Epsi.LT.0) THEN
     CALL MatAXPY(matA,sgnB*Epsi/nu(jv),matVEAb,SUBSET_NONZERO_PATTERN,IERR)     
  ELSE IF(sgnB*Epsi.GT.0) THEN
     CALL MatAXPY(matA,sgnB*Epsi/nu(jv),matVEAf,SUBSET_NONZERO_PATTERN,IERR)
  END IF
  IF(TANG_VM) THEN
     IF(vmconst(jv).LT.0) THEN
        CALL MatAXPY(matA,vmconst(jv)/nu(jv),matVMAf,SUBSET_NONZERO_PATTERN,ierr)
     ELSE
        CALL MatAXPY(matA,vmconst(jv)/nu(jv),matVMAb,SUBSET_NONZERO_PATTERN,ierr)
     END IF
  END IF

  CALL KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
  CALL KSPSetTolerances(ksp,TOL,ATOL,PETSC_DEFAULT_REAL,MAXITS,ierr)
!  CALL KSPSetType(ksp,KSPPREONLY,ierr)
!  CALL KSPGetPC(ksp,pc,ierr)
!  CALL PCSetType(pc,PCILU,ierr)
!  CALL PCFactorSetLevels(pc,10,ierr)
!  CALL PCFactorSetMatOrderingType(pc,"natural",ierr)
  CALL KSPSetType(ksp,KSPPREONLY,ierr)
  CALL KSPGetPC(ksp,pc,ierr)
  CALL PCSetType(pc,PCLU,ierr)
!     CALL PCFactorSetMatOrderingType(pc,"natural",ierr)
  CALL PCFactorSetMatOrderingType(pc,"nd",ierr)
!     CALL PCFactorSetReuseOrdering(pc,PETSC_TRUE,ierr)
!     CALL PCFactorSetReuseFill(pc,PETSC_TRUE,ierr)
!#ifdef PETSC_HAVE_MUMPS
!     call PCFactorSetMatSolverType(pc,MATSOLVERMUMPS,ierr)
!     call PCFactorSetUpMatSolverType(pc,ierr)
!#endif

  CALL KSPSetOperators(ksp,matA,matA,ierr)

  CALL KSPSetUp(ksp,ierr)

  CALL MatDestroy(matA,ierr)

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE FILL_MATRIX_PETSC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#else 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE FILL_MATRIX(npoint,COL,jv,Epsi,VEAf,VEAb,VMAf,VMAb,mat)

!-----------------------------------------------------------------------------------------------
!Use precalculated matrices COL, VEAf, VEAb, VMAf, VMAb and factors Epsi, nu(jv) and vmconst(jv)
!to create matrix mat
!-----------------------------------------------------------------------------------------------  

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER jv,npoint
  REAL*8 Epsi,COL(npoint,npoint)
  REAL*8 VEAf(npoint,npoint),VEAb(npoint,npoint),VMAf(npoint,npoint),VMAb(npoint,npoint)
  !Output
  REAL*8 mat(npoint,npoint)
  !Time
  CHARACTER*30, PARAMETER :: routine="FILL_MATRIX"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)
  
  mat=COL
  IF(sgnB*Epsi.LT.0) THEN
     mat=mat+VEAb*sgnB*Epsi/nu(jv)
  ELSE IF(sgnB*Epsi.GT.0) THEN
     mat=mat+VEAf*sgnB*Epsi/nu(jv)
  END IF
  IF(TANG_VM) THEN
     IF(vmconst(jv).LT.0) THEN
        mat=mat+VMAf*vmconst(jv)/nu(jv)
     ELSE
        mat=mat+VMAb*vmconst(jv)/nu(jv)
     END IF
  END IF

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE FILL_MATRIX

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef MPIandPETSc

SUBROUTINE INVERT_MATRIX_PETSC(nalphab,jv,npoint,BI3,BI7,phi1c,ksp,gint)

!-----------------------------------------------------------------------------------------------
!Fill rhs with some linear combination (depending on phi1c,jv...) of arrays BI3 and BI7 of size
!npoint and invert linear system ksp to obtain g. Use nalphab to skip some helicities in phi1c
!-----------------------------------------------------------------------------------------------  

  USE GLOBAL
  use petscsys
  use petscksp
  IMPLICIT NONE
!#include <petsc/finclude/petscsys.h>
!#include <petsc/finclude/petscvec.h>
!#include <petsc/finclude/petscmat.h>
!#include <petsc/finclude/petscpc.h>
!#include <petsc/finclude/petscksp.h>
!#include <petsc/finclude/petscviewer.h>
  !Input
  INTEGER nalphab,jv,npoint
  REAL*8 BI3(npoint),BI7(npoint,Nnmp),phi1c(Nnmp)
  KSP ksp
  !Output
  REAL*8 gint(npoint,Nnmp)
  !Others
  CHARACTER*100 serr
  INTEGER ii,irhs!,jrhs
  REAL*8 c(npoint),g(npoint)
  PetscErrorCode ierr
  PetscInt indx(npoint)
  Vec vecb,vecx
  !Time
  CHARACTER*30, PARAMETER :: routine="INVERT_MATRIX"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)
  
  phi1c=phi1c !To be removed

  DO ii=1,npoint
     indx(ii)=ii-1
  END DO
  CALL VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,npoint,vecb,ierr)
  CALL VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,npoint,vecx,ierr)
  
  !Scan in possible rhs of the DKE·(corresponding to different contributions the the radial ExB from varphi1)
  DO irhs=1,Nnmp  !Skip helicities too large for the alpha precision
     IF(irhs.GT.1.AND.(ABS(ext_np(irhs)).GT.nalphab/nsamp.OR.ABS(ext_mp(irhs)).GT.nalphab/nsamp)) CYCLE
     CALL VecZeroEntries(vecb,ierr)
     CALL VecZeroEntries(vecx,ierr)
     !Fill rhs and invert
     IF(irhs.EQ.1) THEN
        c=BI3 !radial magnetic drift
!        IF(PHI1_READ) THEN
!           DO jrhs=2,Nnmp !total radial ExB drift added to the magnetic drift
!              c=c+phi1c(jrhs)*BI7(:,jrhs)*2*borbic(0,0)/vdconst(jv)
!           END DO
!        END IF
     ELSE !radial ExB drift
        c=BI7(:,irhs)/vdconst(jv)
     END IF
     CALL VecSetValues(vecb,npoint,indx,c,INSERT_VALUES,ierr)      
     CALL VecAssemblyBegin(vecb,ierr)
     CALL VecAssemblyEnd(vecb,ierr)
     !Solve
     CALL KSPSolve(ksp,vecb,vecx,ierr)
     g=0
     CALL VecGetValues(vecx,npoint,indx,g,ierr)
!     CALL KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD,ierr)
!     CALL PetscMemoryGetCurrentUsage(ierr)
!     CALL PetscViewerASCIIOpen(MPI_COMM_WORLD,'filename.xml',viewer,ierr)
!     CALL PetscViewerPushFormat(viewer,PETSC_VIEWER_DEFAULT,ierr)
!     CALL PetscLogView(viewer,ierr)
     !Distribution function at each point
     gint(:,irhs)=g*vdconst(jv)/nu(jv)
!     gint(:,irhs+1)=gint(:,irhs+1)+g*Sdke(jv)*weight(jv)
     IF(ierr.NE.0) THEN
        serr="Error when inverting the DKE"
        CALL END_ALL(serr,.FALSE.)
     END IF

     IF(.NOT.SOLVE_QN) EXIT
  END DO
  CALL KSPDestroy(ksp,ierr)

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)
  
END SUBROUTINE INVERT_MATRIX_PETSC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#else 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE INVERT_MATRIX(nalphab,jv,npoint,BI3,BI7,phi1c,mat,gint)

!-----------------------------------------------------------------------------------------------
!Fill rhs with some linear combination (depending on phi1c,jv...) of arrays BI3 and BI7 of size
!npoint and invert linear system with matrix mat to obtain g.
!Use nalphab to skip some helicities in phi1c
!-----------------------------------------------------------------------------------------------  

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER nalphab,jv,npoint
  REAL*8 BI3(npoint),BI7(npoint,Nnmp),phi1c(Nnmp),mat(npoint,npoint)
  !Output
  REAL*8 gint(npoint,Nnmp)
  !Others
  LOGICAL ierr
  CHARACTER*100 serr
  INTEGER irhs!,jrhs
  REAL*8 rhs(npoint),ipivot(npoint,npoint),g(npoint)
  !Time
  CHARACTER*30, PARAMETER :: routine="INVERT_MATRIX"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart

  CALL CPU_TIME(tstart)

  phi1c=phi1c !To be removed
    
  DO irhs=1,Nnmp  !Skip helicities too large for the alpha precision
     IF(irhs.GT.1.AND.(ABS(ext_np(irhs)).GT.nalphab/nsamp.OR.ABS(ext_mp(irhs)).GT.nalphab/nsamp)) CYCLE

     IF(irhs.EQ.1) THEN
        rhs=BI3 !radial magnetic drift
!        IF(PHI1_READ) THEN
!           DO jrhs=2,Nnmp !total radial ExB drift added to the magnetic drift
!              rhs=rhs+phi1c(jrhs)*BI7(:,jrhs)*2*borbic(0,0)/vdconst(jv)
!           END DO
!        END IF
     ELSE !radial ExB drift
        rhs=BI7(:,irhs)/vdconst(jv)
     END IF
     !Solve
     CALL DGETRF(npoint,npoint,mat,npoint,ipivot,ierr) 
     CALL DGETRS('No transpose',npoint,1,mat,npoint,ipivot,rhs,npoint,ierr)
     g=rhs
     !Distribution function at each point
     gint(:,irhs)=g*vdconst(jv)/nu(jv)     
!     gint(:,irhs+1)=gint(:,irhs+1)+g*Sdke(jv)*weight(jv)
     IF(ierr) THEN
        serr="Error when inverting the DKE"
        CALL END_ALL(serr,.FALSE.)
     END IF

     IF(.NOT.SOLVE_QN) EXIT
  END DO

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE INVERT_MATRIX

#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE INTEGRATE_G(nalpha,nalphab,nlambda,lambda,i_p,npoint,g,&
     & zetap,thetap,theta,B_al,vds_al,D11,dn1,dn1nm)

!-----------------------------------------------------------------------------------------------
!Calculate contribution D11, dn1 and dn1nm to the flux and to quasineutrality,
! by integrating in lambda grid of nlambda size the distribution function g known at npoints, i_p.
! The integral is calculated in (zetap,thetap) grid of size nalphabxnalpha using precalculated values
! of B_al and vds_al, and then interpolated to a (zetap,thetap) square grid, where flux-surface 
!average is done for the flux
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER nalpha,nalphab,nlambda,i_p(nlambda,nalpha,nalphab),npoint
  REAL*8 zetap(nalphab),thetap(nalpha,nalphab),theta(nalphab)
  REAL*8 lambda(nlambda),g(npoint,Nnmp),B_al(nalpha,nalphab),vds_al(Nnmp,nalpha,nalphab)
  !Output
  REAL*8 D11(Nnmp,Nnmp),dn1(nalphab,nalphab),dn1nm(Nnmp,Nnmp)
  !Others
  REAL*8, PARAMETER :: F5o12 =0.416666666666667
  REAL*8, PARAMETER :: F13o12=1.083333333333333
  INTEGER ia,il,ial,ila,jla,kla,ipoint,jpoint,kpoint,ind,nm
  REAL*8 fdlambda,fhdlambda2,lambdaB,d3vdlambdadK,dlambda,fint
  REAL*8 D11_alp(nalpha*nalphab),D11_ale(3*nalpha),D11_zt(nalphab,nalphab)
  REAL*8 dn1_alp(nalpha*nalphab),dn1_ale(3*nalpha)!,vds(Nnmp,nalpha*nalphab)
  REAL*8 one_o_modB,modB,lambda1,FSA,sqrt1mlB,dummy
  REAL*8 Jac(nalphab,nalphab),dn1_zt(nalphab,nalphab),dn1c_zt(nalphab,nalphab)
  COMPLEX*16 dn1c_nm(nalphab,nalphab)
  INTEGER, SAVE :: tnalpha
  REAL*8, SAVE, ALLOCATABLE :: D11p(:,:),dn1nmp(:,:)
  REAL*8, ALLOCATABLE :: thetape(:,:),vds_zt(:,:,:)
  !Time
  CHARACTER*30, PARAMETER :: routine="INTEGRATE_G"
  INTEGER, SAVE :: ntotal=0
  REAL*8,  SAVE :: ttotal=0
  REAL*8,  SAVE :: t0=0
  REAL*8 tstart
#ifdef MPIandPETSc
  !Others
!  INTEGER ierr
  INCLUDE "mpif.h"
#endif

  CALL CPU_TIME(tstart)
     
  PRE_INTV=(SOLVE_AMB.OR.(SOLVE_QN.AND..NOT..NOT.PHI1_READ.AND..NOT.ZERO_PHI1).OR.NERR.GT.1).AND.CONVERGED
  D11=0
  dn1=0
  dn1nm=0
  IF(DEBUG) THEN

     DO ipoint=1,npoint  
        WRITE(3400+myrank,'(I6,1(1pe13.5),2I6)') ipoint,g(ipoint,1),nalpha,nlambda
!        IF(nalphab.EQ.0) STOP
     END DO

  END IF

  IF(CALCULATED_INT.AND.PRE_INTV) GOTO 123
#ifndef NAG
  IF(plan_fwd.NE.0) THEN
     CALL DFFTW_DESTROY_PLAN(plan_fwd)
     plan_fwd=0!  IF(PRE_INTV.AND.QN)
  END IF
#endif

  IF(ALLOCATED(D11p)) THEN
     DEALLOCATE(D11p)
     IF(QN) DEALLOCATE(dn1nmp)
  END IF
  tnalpha=2*nalpha
  IF(aiota/nzperiod.GE.1) tnalpha=3*nalpha
  ALLOCATE(thetape(tnalpha,nalphab))
  IF(PRE_INTV) THEN
     ALLOCATE(D11p(Nnmp,npoint))
     IF(QN) ALLOCATE(dn1nmp(Nnmp,npoint))
  ELSE
     ALLOCATE(D11p(Nnmp,1))
     IF(QN) ALLOCATE(dn1nmp(Nnmp,1))
  END IF

  D11p=0
  IF(QN) dn1nmp=0
  DO il=1,nalphab
     IF(aiota/nzperiod.LT.1) THEN
        IF(iota.GT.0) THEN
           thetape(     1:nalpha,il)=thetap(1:nalpha,il)-TWOPI
        ELSE
           thetape(     1:nalpha,il)=thetap(1:nalpha,il)+TWOPI
        END IF
        thetape(nalpha+1:tnalpha,il)=thetap(1:nalpha,il)
     ELSE
        IF(iota.GT.0) THEN
           thetape(       1:  nalpha,il)=thetap(1:nalpha,il)-2*TWOPI
           thetape(nalpha+1:2*nalpha,il)=thetap(1:nalpha,il)-  TWOPI
        ELSE
           thetape(       1:  nalpha,il)=thetap(1:nalpha,il)+2*TWOPI
           thetape(nalpha+1:2*nalpha,il)=thetap(1:nalpha,il)+  TWOPI
        END IF
        thetape(2*nalpha+1: tnalpha,il)=thetap(1:nalpha,il)
     END IF
  END DO
!!$  doff=offset-(TWOPI/nalphab)*INT(offset*nalphab/TWOPI)
!!$  DO ia=1,nalphab
!!$!     IF(NTV) THEN 
!!$!        IF(siota.LT.0) THEN
!!$!           theta(ia)=TWOPI-(ia-1)*TWOPI/nalphab!+offsetb-doff
!!$!        ELSE
!!$!           theta(ia)=ia*TWOPI/nalphab!+offsetb-doff
!!$!        END IF
!!$!!        IF(siota.LT.0) THEN
!!$!!           theta(ia)=-TWOPI+(ia-1)*TWOPI/nalphab+offsetb-doff
!!$!!        ELSE
!!$!!           theta(ia)=ia*TWOPI/nalphab+offsetb-doff
!!$!!        END IF
!!$!     ELSE
!!$     IF(siota.LT.0) THEN
!!$        theta(ia)=-TWOPI+ia*TWOPI/nalphab+thetap(1,1)
!!$     ELSE
!!$        theta(ia)=(ia-1)*TWOPI/nalphab+thetap(1,1)
!!$     END IF
!!$!     END IF
!!$  END DO
!!$  imin=MINLOC(ABS(theta(:,1)),1)
!!$  temp(1:nalphab-imin+1,      :)=theta(imin:nalphab)
!!$  temp(nalphab-imin+2:nalphab)=theta(1:imin-1    )!+TWOPI
!!$  theta=temp
  !Precalculate quantities
  lambda1=lambda(1)
  dlambda=lambda(2)-lambda1
!  DO ia=1,nalpha
!     DO il=1,nalphab
!        ial=(il-1)*nalpha+ia
!        vds(:,ial)=vds_al(:,ia,il)
!     END DO
!  END DO
!     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
     ALLOCATE(vds_zt(Nnmp,nalphab,nalphab))
     DO ia=1,nalphab
        DO il=1,nalphab
           CALL FILL_BNODE(zetap(il),theta(ia),Jac(ia,il),dummy,vds_zt(:,ia,il),.FALSE.)
        END DO
     END DO

  DO ipoint=1,npoint

     D11_alp=0
     IF(QN) dn1_alp=0
     !Scan in the flux surface
     DO ia=1,nalpha 
        IF(.NOT.PRE_INTV.AND.ipoint.GT.1) EXIT
        !Scan in the flux surface
        DO il=1,nalphab 
           ial=(il-1)*nalpha+ia
           kla=-1
           modB=B_al(ia,il)
           one_o_modB=1./modB
           jla=INT((one_o_modB-lambda1)/dlambda)        
           IF(jla.LE.1) jla=2
           IF(jla.NE.0) fdlambda=(one_o_modB-lambda(jla))/dlambda
           kpoint=npoint
           !Integral in lambda
           DO ila=jla,2,-1 
              jpoint=i_p(ila,ia,il)
              fint=1
              IF(jpoint.EQ.npoint) THEN
                 jpoint=kpoint
                 IF(jpoint.EQ.npoint) CYCLE
                 fint=(ila-1.)/(kla-1.)
              ELSE
                 kpoint=jpoint
                 kla=ila
              END IF
              IF(PRE_INTV.AND.jpoint.NE.ipoint) CYCLE
              fhdlambda2=0.5*fdlambda*fdlambda
              IF(ila.EQ.jla) THEN
                 fint=fint*(F5o12+fdlambda+fhdlambda2)
              ELSE IF(ila.EQ.jla-1) THEN
                 fint=fint*(F13o12-fhdlambda2)
              ELSE IF (ila.EQ.2) THEN
                 fint=fint*F13o12
              END IF
              lambdaB=lambda(ila)*modB
              sqrt1mlB=SQRT(1.-lambdaB)
              d3vdlambdadK=fint*modB/sqrt1mlB*dlambda
              IF(PRE_INTV) THEN   !Calculate the linear contribution of g(ipoint)
                 D11_alp(ial)=D11_alp(ial)-d3vdlambdadK*vds_al(1,ia,il)*(1.0-0.5*lambdaB) /2.
                 IF(QN) dn1_alp(ial)=dn1_alp(ial)+d3vdlambdadK       
              ELSE                !Accumulate the total contribution g for each point jpoint
                 IF(.NOT.PHI1_READ.OR..NOT.IMP1NU) THEN
                    D11_alp(ial)=D11_alp(ial)-g(jpoint,1)*d3vdlambdadK*vds_al(1,ia,il)*(1.0-0.5*lambdaB)/2.
                 ELSE
                    D11_alp(ial)=D11_alp(ial)+g(jpoint,1)*&
                         & modB/(lambda(ila)*lambda(ila)*lambda(ila)*sqrt1mlb)*vds_al(1,ia,il)
                    
                 END IF
                 IF(QN) dn1_alp(ial)=dn1_alp(ial)+g(jpoint,1)*d3vdlambdadK        
              END IF
           END DO

           IF(DEBUG.AND.ipoint.EQ.1) WRITE(3451+myrank,'(4(1pe13.5),10I6)') &
                & zetap(il),thetap(ia,il),dn1_alp(ial),D11_alp(ial)

        END DO

     END DO

     !Copy values to extended grid and interpolate to square grid
     dn1_zt=0
     DO il=1,nalphab
        ial=(il-1)*nalpha
        IF(aiota/nzperiod.LT.1) THEN
           D11_ale(        1: nalpha)=D11_alp(ial+1:ial+nalpha)
           D11_ale( nalpha+1:tnalpha)=D11_alp(ial+1:ial+nalpha)
           IF(QN) THEN
              dn1_ale(        1: nalpha)=dn1_alp(ial+1:ial+nalpha)
              dn1_ale( nalpha+1:tnalpha)=dn1_alp(ial+1:ial+nalpha)
           END IF
        ELSE
           D11_ale(         1:  nalpha)=D11_alp(ial+1:ial+nalpha)
           D11_ale(  nalpha+1:2*nalpha)=D11_alp(ial+1:ial+nalpha)
           D11_ale(2*nalpha+1: tnalpha)=D11_alp(ial+1:ial+nalpha)
           IF(QN) THEN
              dn1_ale(         1:  nalpha)=dn1_alp(ial+1:ial+nalpha)
              dn1_ale( nalpha+1: 2*nalpha)=dn1_alp(ial+1:ial+nalpha)
              dn1_ale(2*nalpha+1: tnalpha)=dn1_alp(ial+1:ial+nalpha)
           END IF
        END IF

!        DO ia=1,tnalpha

!        END DO
        
        !Interpolation to square grid
        DO ia=1,nalphab
           ind=(ia-1)*nalphab+il
           D11_zt(ia,il)=0
           CALL LAGRANGE(thetape(1:tnalpha,il),D11_ale(1:tnalpha),tnalpha,&
                & theta(ia),D11_zt(ia,il),2)
           IF(QN.OR.SATAKE)THEN
              dn1_zt(ia,il)=0
              CALL LAGRANGE(thetape(1:tnalpha,il),dn1_ale(1:tnalpha),tnalpha,&
                 & theta(ia),dn1_zt(ia,il),2)
           END IF
           IF(DEBUG.AND.ipoint.EQ.1) WRITE(3450+myrank,'(10(1pe13.5),10I6)') &
                & zetap(il),theta(ia),dn1_zt(ia,il),D11_zt(ia,il)
        END DO
     END DO
     !Flux surface average
     D11p(1,ipoint)=FSA(nalphab,nalphab,D11_zt,Jac,1)

     IF(QN) THEN
        DO nm=2,Nnmp
           IF(ONLY_PHI1) EXIT 
           D11p(nm,ipoint)=FSA(nalphab,nalphab,dn1_zt*2*vds_zt(nm,:,:),Jac,1)
        END DO
        DO ia=1,nalphab
           DO il=1,nalphab
              dn1c_zt(il,ia)=dn1_zt(ia,il)              
           END DO
        END DO
        CALL FFTF_KN(nalphab,dn1c_zt,dn1c_nm)
        CALL FILL_ORBI(nalphab,nalphab,dn1c_nm,Nnmp,dn1nmp(:,ipoint))
     END IF
     IF(.NOT.PRE_INTV.AND.ipoint.EQ.1) EXIT
  END DO


123 IF(PRE_INTV) THEN

     DO nm=1,Nnmp
        IF(QN) CALL DGEMV('n',Nnmp,npoint,ONE,&
             & dn1nmp(:,:),Nnmp,g(:,nm),1,ZERO,dn1nm(:,nm),1)
        IF(.NOT.ONLY_PHI1) THEN
           CALL DGEMV('n',Nnmp,npoint,ONE,&
                & D11p(:,:),Nnmp,g(:,nm),1,ZERO,D11(:,nm),1)
        ELSE IF(nm.EQ.1) THEN
           CALL DGEMV('n',1,npoint,ONE,&
             & D11p(1,:),1,g(:,1),1,ZERO,D11(1,1),1)
        END IF
     END DO
!     CALL DGEMM('n','n',Nnmp,Nnmp,npoint-1,ONE,&
!          & D11p(:,1:npoint-1)  ,Nnmp,g(1:npoint-1,:),npoint-1,ZERO,D11,Nnmp)
!     IF(QN) CALL DGEMM('n','n',Nnmp,Nnmp,npoint-1,ONE,&
!          & dn1nmp(:,1:npoint-1),Nnmp,g(1:npoint-1,:),npoint-1,ZERO,dn1nm,Nnmp) 
  ELSE
     D11(1,1)=D11p(1,1)
     IF(QN) THEN
        dn1=dn1_zt
        dn1nm(:,1)=dn1nmp(:,1)
     END IF
 END IF

 CALCULATED_INT=.TRUE.

  CALL CALCULATE_TIME(routine,ntotal,t0,tstart,ttotal)

END SUBROUTINE INTEGRATE_G


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE FILL_ORBI(numz,numt,qnmc_in,num,qnm_out)

!-------------------------------------------------------------------------------------------------
!Read qnmc_in(numz,numt) from FFT and write num Fourier modes in 'ddkes2.data' format
!-------------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER numz,numt,num
  COMPLEX*16 qnmc_in(numz,numt)
  !Output
  REAL*8 qnm_out(num)
!  REAL qnm_out(num)
  !Others
  INTEGER n,m,nm,ind1,ind2,ind3,ind4
  REAL*8 qorbic(-ntorb:ntorb,0:mpolb),qorbis(-ntorb:ntorb,0:mpolb)
  REAL*8 qorbic_temp(-ntorb:ntorb,-mpolb:mpolb),qorbis_temp(-ntorb:ntorb,-mpolb:mpolb)
  
  qorbic_temp=0
  qorbis_temp=0
  qorbic    =0
  qorbis    =0
  
  DO m=-mpolb,mpolb
     DO n=-ntorb,ntorb
        IF(n.EQ.0.OR.m.EQ.0) THEN 
           IF(m.EQ.0.AND.n.EQ.0) THEN
              qorbic_temp(0,0)=REAL(qnmc_in(1,1))
              qorbis_temp(0,0)=AIMAG(qnmc_in(1,1))
           ELSE IF(n.EQ.0) THEN
              IF(m.GT.0) THEN
                 ind1=m+1
                 ind2=numt-m+1
                 IF(ind1.GT.numt.OR.ind2.LT.1) CYCLE
                 qorbic_temp(0,m)= REAL( qnmc_in(1,ind1)+qnmc_in(1,ind2))
                 qorbis_temp(0,m)=AIMAG(-qnmc_in(1,ind1)+qnmc_in(1,ind2))
              ELSE
                 ind1=1-m
                 ind2=numt+m+1
                 IF(ind1.GT.numt.OR.ind2.LT.1) CYCLE
                 qorbic_temp(0,-m)= REAL( qnmc_in(1,ind1)+qnmc_in(1,ind2))
                 qorbis_temp(0,-m)=AIMAG(-qnmc_in(1,ind1)+qnmc_in(1,ind2))
              END IF
           ELSE IF(m.EQ.0) THEN
              IF(n.GT.0) THEN
                 ind1=n+1
                 ind2=numz-n+1
                 IF(ind1.GT.numz.OR.ind2.LT.1) CYCLE
                 qorbic_temp(n,0)=qorbic_temp(n,0) +REAL( qnmc_in(ind1,1)-qnmc_in(ind2,1))
                 qorbis_temp(n,0)=qorbis_temp(n,0)+AIMAG(-qnmc_in(ind1,1)+qnmc_in(ind2,1))/2
              ELSE
                 ind1=1-n
                 ind2=numz+n+1
                 IF(ind1.GT.numz.OR.ind2.LT.1) CYCLE
                 qorbic_temp(-n,0)=qorbic_temp(-n,0)+ REAL( qnmc_in(ind1,1)+qnmc_in(ind2,1))
                 qorbis_temp(-n,0)=qorbis_temp(-n,0)+AIMAG(-qnmc_in(ind1,1)+qnmc_in(ind2,1))/2
              END IF
           END IF
        ELSE IF(n.GT.0) THEN
           IF(m.GT.0) THEN
              ind1=n+1
              ind2=m+1
              ind3=numz-n+1
              ind4=numt-m+1
              IF(ind1.GT.numz.OR.ind2.GT.numt.OR.ind3.LT.1.OR.ind4.LT.1) CYCLE
              qorbic_temp(n,m)= REAL(+qnmc_in(ind1,ind2)+qnmc_in(ind3,ind4))/2
              qorbis_temp(n,m)=AIMAG(-qnmc_in(ind1,ind2)+qnmc_in(ind3,ind4))
           ELSE
              ind1=numz-n+1
              ind2=1-m
              ind3=n+1
              ind4=numt+m+1
              IF(ind1.LT.1.OR.ind2.GT.numt.OR.ind3.GT.numz.OR.ind4.LT.1) CYCLE
              qorbic_temp(n,m)= REAL(+qnmc_in(ind1,ind2)+qnmc_in(ind3,ind4))/2
              qorbis_temp(n,m)=AIMAG(+qnmc_in(ind1,ind2)-qnmc_in(ind3,ind4))
           END IF
        ELSE IF(n.LT.0) THEN
           IF(m.GT.0) THEN
              ind1=1-n
              ind2=numt-m+1
              ind3=numz+n+1
              ind4=m+1
              IF(ind1.GT.numz.OR.ind2.LT.1.OR.ind3.LT.1.OR.ind4.GT.numt) CYCLE
              qorbic_temp(n,m)= REAL(+qnmc_in(ind1,ind2)+qnmc_in(ind3,ind4))/2
              qorbis_temp(n,m)=AIMAG(+qnmc_in(ind1,ind2)-qnmc_in(ind3,ind4))
           ELSE
              ind1=numz+n+1
              ind2=numt+m+1
              ind3=1-n
              ind4=1-m
              IF(ind1.LT.1.OR.ind2.LT.1.OR.ind3.GT.numz.OR.ind4.GT.numt) CYCLE
              qorbic_temp(n,m)= REAL(+qnmc_in(ind1,ind2)+qnmc_in(ind3,ind4))/2
              qorbis_temp(n,m)=AIMAG(+qnmc_in(ind1,ind2)-qnmc_in(ind3,ind4))
           END IF
        END IF
     END DO
  END DO

  DO m=-mpolb,mpolb
     DO n=-ntorb,ntorb
        IF(m.GE.0) THEN
           qorbic(n,m)=qorbic(n,m)+qorbic_temp(n,m)
           qorbis(n,m)=qorbis(n,m)+qorbis_temp(n,m)
        ELSE
           qorbic(-n,-m)=qorbic(-n,-m)+qorbic_temp(n,m)
!           qorbis(-n,-m)=qorbis(-n,-m)-qorbis_temp(n,m)
        END IF
     END DO
  END DO

  nm=0
  DO m=0,mpolb
     DO n=-ntorb,ntorb
        IF(m.EQ.0.AND.n.LT.0) CYCLE
        nm=nm+1
        qnm_out(nm)=qorbic(n,m)
     END DO
  END DO
  DO m=0,mpolb
     DO n=-ntorb,ntorb
        IF(m.EQ.0.AND.n.LT.0) CYCLE
        nm=nm+1
        qnm_out(nm)=qorbis(n,m)
     END DO
  END DO

END SUBROUTINE FILL_ORBI


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



