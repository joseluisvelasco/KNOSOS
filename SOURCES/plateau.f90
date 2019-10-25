!Calculate neoclassical transport in the plateau regime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_PLATEAU(jv,Epsi,D11)

!----------------------------------------------------------------------------------------------- 
!Calculate monoenergetic transport coefficient D11 in the plateau regime for
!collisionality cmul=nu(jv)/v(jv) and normalized radial electric field efied=Epsi/v(jv)
!-----------------------------------------------------------------------------------------------
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER jv
  REAL*8 Epsi
  !Output
  REAL*8 D11
  !Others
  INTEGER m,n
  REAL*8 G11dkes,cmul
  COMPLEX*16 bmn_e(-ntorbd:ntorbd,-mpolbd:mpolbd)
  REAL*8, SAVE :: D11p=1

  WRITE(1000+myrank,*) 'Calculating PLATEAU'

  IF(D11p.LT.0) GO TO 999

  !Change Fourier representation
  bmn_e=0
  DO m=-mpolb,mpolb
     DO n=-ntorb,ntorb
        IF(m.LT.0) bmn_e(n,m)= borbic(-n,-m)/2.
        IF(m.GT.0) bmn_e(n,m)= borbic(+n,+m)/2.
        IF(m.EQ.0) bmn_e(n,0)=(borbic(+n,0)+borbic(-n,0))/2.
        IF(STELL_ANTISYMMETRIC) THEN
           IF(m.LT.0) bmn_e(n,m)=bmn_e(n,m)-NUMI*borbis(-n,-m)/2.
           IF(m.GT.0) bmn_e(n,m)=bmn_e(n,m)+NUMI*borbis(+n,+m)/2.
           IF(m.EQ.0) bmn_e(n,0)=bmn_e(n,0)+NUMI*(borbis(+n,0)-borbic(-n,0))/2.
        END IF
     END DO
  END DO
  D11=0 
  DO m=-mpolb,mpolb
     DO n=-ntorb,ntorb
        IF(n.EQ.0.AND.m.EQ.0) CYCLE
        D11=D11-(m*Bzeta-n*nzperiod*Btheta)*(m*Bzeta-n*nzperiod*Btheta)&
             & *REAL(bmn_e(n,m)*bmn_e(-n,-m))/ABS(iota*m+nzperiod*n) 
     END DO
  END DO 

  D11p=D11

999 D11=-0.125*pi*(vdconst(jv)*vdconst(jv)/v(jv))*(rad_R/borbic(0,0)/borbic(0,0)/aiBtpBz/aiBtpBz)*D11p 

  cmul=nu(jv)/v(jv)/2.
  !Connect with Pfirsch-Schlueter or 1/nu
  IF(FACT_CON.GT.0.AND.cmul_PS.GT.0) THEN
     IF(cmul.GT.cmul_PS /FACT_CON) D11=D11*(1+cmul/cmul_PS)
     IF(cmul.LT.cmul_1NU*FACT_CON) D11=D11*(1+cmul_1NU/cmul)
  END IF

  G11dkes=fdkes(jv)*D11
  WRITE(200+myrank,'(3(1pe13.5)," NaN ",2(1pe13.5),"  NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN")') &
       & nu(jv)/v(jv)/2,Epsi/v(jv)*psip,vmconst(jv)/v(jv),G11dkes,G11dkes
  IF(DEBUG) THEN
     IF(cmul_PS.GT.0) THEN !to be checked
        WRITE(10000+myrank,'("2 ",6(1pe13.5))') nu(jv)/v(jv)/2,Epsi/v(jv)*psip,vmconst(jv)/v(jv),&
             & G11dkes,weight(jv)/fdkes(jv),weight(jv)/fdkes(jv)*v(jv)*v(jv)
     ELSE
        WRITE(10000+myrank,'("0 ",6(1pe13.5))') nu(jv)/v(jv)/2,Epsi/v(jv)*psip,vmconst(jv)/v(jv),&
             & G11dkes,weight(jv)/fdkes(jv),weight(jv)/fdkes(jv)*v(jv)*v(jv)
     END IF
  END IF

END SUBROUTINE CALC_PLATEAU


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


