!Calculate neoclassical transport in the banana regime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE CALC_BANANA(jv,Epsi,D11)

!----------------------------------------------------------------------------------------------- 
!Calculate monoenergetic transport coefficient in the banana regime for collisionality
!cmul=nu(jv)/v(jv) and normalized radial electric field efied=Epsi/v(jv)
!-----------------------------------------------------------------------------------------------
  
  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER jv
  REAL*8 Epsi
  !Output
  REAL*8 D11
  !Others
  INTEGER, PARAMETER :: ntx    = 128
  INTEGER, PARAMETER :: nzx    = 128
  INTEGER iz,it
  REAL*8 G11dkes,cmul
  REAL*8, SAVE :: D11p=1
!  REAL*8 B_0
  REAL*8 dz, zeta(nzx)
  REAL*8 dt,theta(ntx)

  
  WRITE(1000+myrank,*) 'Calculating BANANA'

  IF(D11p.LT.0) GO TO 999

  dz=TWOPI/nzx/nzperiod
  dt=TWOPI/ntx
  DO iz=1,nzx
     zeta(iz)=(iz-1)*dz
  END DO
  DO it=1,ntx
     theta(it)=(it-1)*dt
  END DO
     
999 cmul=nu(jv)/v(jv)/2.
  !Connect with Pfirsch-Schlueter or 1/nu
!  IF(FACT_CON.GT.0.AND.cmul_PS.GT.0) THEN
!     IF(cmul.GT.cmul_PS /FACT_CON) D11=D11*(1+cmul/cmul_PS)
!     IF(cmul.LT.cmul_1NU*FACT_CON) D11=D11*(1+cmul_1NU/cmul)
!  END IF

  G11dkes=fdkes(jv)*D11
  WRITE(200+myrank,'(2(1pe13.5)," NaN NaN ",2(1pe13.5)," NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN")') &
       & nu(jv)/v(jv)/2.,Epsi/v(jv)*psip,fdkes(jv)*D11,fdkes(jv)*D11
  IF(DEBUG) THEN
     IF(cmul_PS.GT.0) THEN
        WRITE(10000+myrank,'("3 ",5(1pe13.5))') nu(jv)/v(jv),Epsi/v(jv)*psip,G11dkes,&
             & weight(jv)/fdkes(jv),weight(jv)/fdkes(jv)*v(jv)*v(jv)
     ELSE
        WRITE(10000+myrank,'("0 ",5(1pe13.5))') nu(jv)/v(jv),Epsi/v(jv)*psip,G11dkes,&
             & weight(jv)/fdkes(jv),weight(jv)/fdkes(jv)*v(jv)*v(jv)
     END IF
  END IF

END SUBROUTINE CALC_BANANA


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


