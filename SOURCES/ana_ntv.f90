!DEBUGGING Calculates NTV in a large aspect-ratio circular tokamak using analytical formulas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CALC_NTV(jt,s,ib,Zb,Ab,nb,dnbdpsi,Tb,dTbdpsi,Epsi,Gb)              

!----------------------------------------------------------------------------------------------- 
!Calculates elliptic integrals of first and second kind
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER jt,ib
  REAL*8 s,Zb,Ab,nb,dnbdpsi,Tb,dTbdpsi,Epsi
  !Output
  REAL*8 Gb
  !Others
  INTEGER, PARAMETER :: nk2=1000
  REAL*8, PARAMETER :: k0=0.827
!  REAL*8, PARAMETER :: k02=0.827
  INTEGER n,ik2
  REAL*8 eta1,eta2,k,ftheta,dk2!,eps
  REAL*8 nu_s,I1nu,Isqrtnu,Isbp,nustard,k2,detadv(nv),x1(nv),x2(nv)
!  REAL*8 alphan,betan,eK,eE,dummy,k0
  REAL*8 alphan,betan,eK,eE,dummy,vdummy(Nnmp),k02
  REAL*8 Gb_1nu,Gb_sbp,Gb_sqrtnu
  REAL*8 borbh(-ntorbd:ntorbd,0:mpolbd)
  INTEGER, PARAMETER :: npt=512
  INTEGER m,iz,it
  REAL*8 dzeta,dtheta,func,Gnm,dG,logn
  REAL*8 zeta(npt,npt),theta(npt,npt),B(npt,npt),Bh(npt,npt)!Bh(nax,nax)
  REAL*8 zeta_h(npt,npt),theta_h(npt,npt)
  REAL*8 mode(-ntorbd:ntorbd,-mpolbd:mpolbd)

  s=s
  
  dzeta =TWOPI/npt/nzperiod
  dtheta=TWOPI/npt
  mode=0
  DO iz=1,npt
     DO it=1,npt
        zeta(iz,it) =iz*dzeta
        theta(iz,it)=it*dtheta
        CALL CALCB(zeta(iz,it),theta(iz,it),0,.FALSE.,B(iz,it),dummy,dummy,&
             & dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,vdummy)
        func=avB2/B(iz,it)/B(iz,it)
!        func=(avB2/(B(iz,it)*B(iz,it))-1)
        DO n=-ntorb,ntorb
           DO m=0,mpolb
              mode(n,m)=mode(n,m)+&
                   & func*COS(m*theta(iz,it)+n*nzperiod*zeta(iz,it))
           END DO
        END DO
     END DO
  END DO
  mode=mode*dzeta*dtheta/PI/TWOPI
!  mode(0,0)=mode(0,0)/2.
  mode(0,0)=0!mode(0,0)/2.
  zeta_h =zeta
  theta_h=theta
  DO n=-ntorb,ntorb
     DO m=0,mpolb
        IF(n.EQ.0.AND.m.EQ.0) CYCLE
        Gnm=mode(n,m)/(iota*m+nzperiod*n)
        DO iz=1,npt
           DO it=1,npt
              dG=Gnm*SIN(m*theta(iz,it)+n*nzperiod*zeta(iz,it))
              zeta_h(iz,it) = zeta_h(iz,it)+dG
              theta_h(iz,it)=theta_h(iz,it)+dG*iota
           END DO
        END DO
     END DO
  END DO
  CALL CHANGE_COORDS(dummy,npt,npt,zeta_h,theta_h,B,0,npt,Bh)
  borbh=0
  DO iz=1,npt
     DO it=1,npt
        zeta(iz,it) =iz*dzeta
        theta(iz,it)=it*dtheta
        func=Bh(iz,it)
        DO n=-2*ntorb,2*ntorb
           DO m=0,2*mpolb
              borbh(n,m)=borbh(n,m)+&
                   & func*COS(m*theta(iz,it)+n*nzperiod*zeta(iz,it))
           END DO
        END DO
     END DO
  END DO
  borbh=borbh*dzeta*dtheta/PI/TWOPI
  borbh(0,0)=borbh(0,0)/2.

  
!  eps=SQRT(s)*rad_a/rad_R
  x1=v/vth(ib)
  x2=x1*x1
  nu_s=nuth(ib)/2

!  DO n=-ntorb,ntorb
!     DO m=0,mpolb
!        print *,n,m,borbh(n,m),borbic(n,m)
!     END DO
!  END DO
!  STOP
  borbic=borbh

  detadv=0.5*x2*x2*x1*nuth(ib)/nu
  detadv=detadv*1./x1
  CALL INT_VMAX(detadv,weight,eta1)
  detadv=detadv*(x2-2.5)
  CALL INT_VMAX(detadv,weight,eta2)
  dk2=1./nk2
  I1nu=0
  DO ik2=1,nk2
     k2=(ik2-0.5)*dk2
     k=SQRT(k2)
     ftheta=0
     DO n=-2*ntorb,2*ntorb
!        CALL CALC_ALPHABETA(k,n,alphan,betan,dummy,dummy)
        ftheta=ftheta+n*n*nzperiod*nzperiod*(alphan*alphan+betan*betan)
     END DO
     CALL ELLIPTIC(k,eK,eE)
     I1nu=I1nu+ftheta/(eE-(1-k2)*eK)
  END DO
  I1nu=I1nu*dk2
  Gb_1nu=-(1/4./SQ2/PI/SQPI)*(vth(ib)*vth(ib)*vth(ib)*vth(ib)/nu_s)*(Ab*m_e/Zb/iota)*(Ab*m_e/Zb/iota)*(eps)*SQRT(eps)*&
       & I1nu*(eta1*(dnbdpsi/nb+dTbdpsi/Tb-Zb*Epsi/Tb)+eta2*dTbdpsi/Tb)

!!!!!!!!!!!!!!!!

  detadv=0.5*x2*x2*x1*SQRT(nu/nuth(ib))
  detadv=detadv*1./x1
  CALL INT_VMAX(detadv,weight,eta1)
  detadv=detadv*(x2-2.5)
  CALL INT_VMAX(detadv,weight,eta2)
  nustard=ABS(4*nu_s*iota/(eps*Epsi))
  logn=LOG(16./SQRT(nustard))
  dk2=SQRT(nustard/logn)
  dk2=1e-9
  Isqrtnu=0
  DO n=-ntorb,ntorb
     CALL CALC_ALPHABETA(SQRT(1.-dk2),n,dummy,dummy,alphan,betan)
     Isqrtnu=Isqrtnu+SQRT(ABS(1.0*n*nzperiod))*(alphan*alphan+betan*betan)
!     IF(n.EQ.3) print *,'dn',Isqrtnu,alphan*alphan
  END DO
  CALL ELLIPTIC(SQRT(1.-dk2),eK,eE)
!  print *,'dk2',Isqrtnu,eK
  Isqrtnu=Isqrtnu/16./eK/eK
  Gb_sqrtnu=-(1./4./SQ2/PI/SQPI)*(vth(ib)*vth(ib)*vth(ib)*vth(ib))*(Ab*m_e/Zb/iota)*(Ab*m_e/Zb/iota)* &
       & ABS(iota/Epsi)*SQRT(nustard*eps*ABS(logn))*Isqrtnu* &
       & (eta1*(dnbdpsi/nb+dTbdpsi/Tb-Zb*Epsi/Tb)+eta2*dTbdpsi/Tb)

!  stop
!!!!!!!!!!!!!!!!1

!  k0=SQRT(k02)
  k02=k0*k0
  Isbp=0
  DO n=-2*ntorb,2*ntorb
!     CALL CALC_ALPHABETA(k0,n,dummy,dummy,alphan,betan)
     Isbp=Isbp+ABS(n*nzperiod)*(alphan*alphan+betan*betan)
  END DO
  CALL ELLIPTIC(k0,eK,dummy)
  Isbp=Isbp/16./eK/eK
  eta1=(3./4.)*SQPI*k02*(1-k02)*4.*eK
  Gb_sbp=-(1/SQ2/SQPI)*vth(ib)*vth(ib)*(Ab*m_e/ABS(Zb))*rad_R*SQRT(eps)*(psip/ABS(iota))* &
       & Isbp*eta1*(dnbdpsi/nb+dTbdpsi/Tb-Zb*Epsi/Tb)
 
  IF(ABS(Epsi).GT.1E-5) THEN
     Gb=Gb_sqrtnu
  ELSE
     IF(jt.EQ.1) THEN
        Gb=Gb_sbp
     ELSE IF(jt.EQ.2) THEN
        Gb=Gb_1nu
     END IF
  END IF
  
END SUBROUTINE CALC_NTV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_ALPHABETA(k,n,alpha_1nu,beta_1nu,alpha_sbp,beta_sbp)

!----------------------------------------------------------------------------------------------- 
!Calculates alpha_n and beta_n
!----------------------------------------------------------------------------------------------- 
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER n
  REAL*8 k
  !Output
  REAL*8 alpha_1nu,beta_1nu,alpha_sbp,beta_sbp
  !Others
  INTEGER, PARAMETER :: ntheta=10000
  INTEGER itheta
  REAL*8 k2,thetamin,thetamax,theta,dtheta,An,Bn,sine,sqr

  k2=k*k
  thetamin=-2*ASIN(k)
  thetamax=+2*ASIN(k)
  alpha_1nu=0
  alpha_sbp=0
  beta_1nu=0
  beta_sbp=0
  dtheta=(thetamax-thetamin)/ntheta
  DO itheta=1,ntheta
     theta=thetamin+(itheta-0.5)*dtheta
     sine=SIN(theta/2.)
     sqr=SQRT(k2-sine*sine)
     CALL CALC_AN_AND_BN(theta,n,An,Bn)
     alpha_1nu=alpha_1nu+An*sqr
     beta_1nu =beta_1nu +Bn*sqr
     alpha_sbp=alpha_sbp+An/sqr
     beta_sbp =beta_sbp +Bn/sqr
!     IF(n.EQ.3) print *,'ab',An*An,itheta,theta!borbh(3,6)*COS((6+nzperiod*3/iota)*theta),An,itheta
  END DO
  alpha_1nu=alpha_1nu*dtheta
  beta_1nu =beta_1nu *dtheta
  alpha_sbp=alpha_sbp*dtheta
  beta_sbp =beta_sbp *dtheta

END SUBROUTINE CALC_ALPHABETA


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_AN_AND_BN(theta,n,An,Bn)

!----------------------------------------------------------------------------------------------- 
!Calculates An(theta) and Bn(theta)
!----------------------------------------------------------------------------------------------- 
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER n
  REAL*8 theta
  !Output
  REAL*8 An,Bn
  !Others
  INTEGER m

  An=0
  Bn=0
  DO m=0,2*mpolb
     IF(n.EQ.0) THEN
        IF(m.EQ.0.OR.m.EQ.1) CYCLE
     END IF
     An=An-borbic(n,m)*COS((m+nzperiod*n/iota)*theta)
     Bn=Bn+borbic(n,m)*SIN((m+nzperiod*n/iota)*theta)
!     IF(theta.GT.1.947) print *,'ac',n,m,An*An
  END DO
  An=An/borbic(0,0)
  Bn=Bn/borbic(0,0)

END SUBROUTINE CALC_AN_AND_BN


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE ELLIPTIC(k,eK,eE)

!----------------------------------------------------------------------------------------------- 
!Calculates elliptic integrals of first and second kind
!----------------------------------------------------------------------------------------------- 

  USE GLOBAL
  IMPLICIT NONE
  !Input
  REAL*8 k,eK,eE
  !Others
  INTEGER, PARAMETER :: ntheta=10000
  INTEGER itheta
  REAL*8 k2,dtheta,theta,sine,sqroot

  eK=0.0
  eE=0.0
  k2=k*k

  dtheta=PI/2./ntheta
  DO itheta=1,ntheta
     theta=(itheta-0.5)*dtheta
     sine=SIN(theta)
     sqroot=SQRT(1.-k2*sine*sine)
     eK=eK+1./sqroot
     eE=eE+   sqroot
  END DO
  eK=eK*dtheta
  eE=eE*dtheta

  
END SUBROUTINE ELLIPTIC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE INT_VMAX(Q,weightb,intQ)

!----------------------------------------------------------------------------------------------- 
!Calculates the integral in v intQ of Q, weighed by weightb (convolution with Maxwellian)
!-----------------------------------------------------------------------------------------------
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  REAL*8 Q(nv),weightb(nv)
  !Output
  REAL*8 intQ
  !Others
  INTEGER iv

  intQ=0
  DO iv=1,nv
     intQ=intQ+weightb(iv)*Q(iv)
  END DO
  
END SUBROUTINE INT_VMAX


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
