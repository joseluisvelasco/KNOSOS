
!Solve transport REVISE ALL

SUBROUTINE TRANSPORT(nbb,ns,dt,s,Zb,Ab,nb,dnbdpsi,Sb,Tb,dTbdpsi,Pb,Es)

!----------------------------------------------------------------------------------------------- 
!For densities and temperatures and electrostatic potential given by nb, dnbdpsi, Tb, dTbdpsi,
!and Epsi for nbb species of charge Zb and mass Ab at ns surfaces s, calculate evolution 
!after time step dt
!----------------------------------------------------------------------------------------------- 
  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER nbb,ns
  REAL*8 dt,s(ns),Zb(nbb),Ab(nbb),Gb(nbb,ns),Sb(nbb,ns),Qb(nbb,ns),Pb(nbb,ns)
  !Input/output
  REAL*8 nb(nbb,ns),dnbdpsi(nbb,ns),Tb(nbb,ns),dTbdpsi(nbb,ns),Es(ns)
  !Others
  REAL*8, PARAMETER  ::  prefact_Es=1213173.45142083 !e/(8pi^2m)
  INTEGER ib,is
  REAL*8 gradGb(nbb,ns),gradQb(nbb,ns),ohm(nbb,ns),coulomb(nbb,ns)
  REAL*8 dnbdt(nbb,ns),dTbdt(nbb,ns),dnbTbdt(nbb,ns),dEsdt(ns)
  REAL*8 fdummy,nu0(ns),fact_coulomb(ns),fact_Es(ns)

#ifdef MPIandPETSc

  INCLUDE "mpif.h"

  IF(numprocs.GT.1) THEN
     DO ib=1,nbb
        CALL REAL_ALLREDUCE(     nb(ib,:),ns)
        CALL REAL_ALLREDUCE(dnbdpsi(ib,:),ns)
        CALL REAL_ALLREDUCE(     Gb(ib,:),ns)
        CALL REAL_ALLREDUCE(     Sb(ib,:),ns)
        CALL REAL_ALLREDUCE(     Tb(ib,:),ns)
        CALL REAL_ALLREDUCE(dnbdpsi(ib,:),ns)
        CALL REAL_ALLREDUCE(dTbdpsi(ib,:),ns)
        CALL REAL_ALLREDUCE(     Qb(ib,:),ns)
        CALL REAL_ALLREDUCE(     Pb(ib,:),ns)
     END DO
     CALL REAL_ALLREDUCE( s,ns)
     CALL REAL_ALLREDUCE(Es,ns)
  END IF

#endif
  
  WRITE(1000+myrank,*) 'SUBROUTINE TRANSPORT not ready'
  RETURN

  DO ib=1,nbb
     CALL DERIVE(s(1:ns),Gb(ib,1:ns),ns,4,gradGb(ib,1:ns))
     CALL DERIVE(s(1:ns),Qb(ib,1:ns),ns,4,gradQb(ib,1:ns))
  END DO

  ohm(1,:)=-Es*Gb(1,:) !XXX
  ohm(2,:)=-ohm(1,:)
  DO is=1,ns
     CALL CALC_CTS(Zb(2),Ab(2),nb(2,is),Tb(2,is),fdummy,nu0(is))
  END DO
  fact_coulomb=4/sqpi*Zb(2)*Ab(1)/Ab(2)*nu0
  coulomb(1,:)=3*(Tb(1,:)-Tb(2,:))*fact_coulomb
  coulomb(2,:)=-coulomb(1,:)

  Sb=gradGb

  dnbdt  =-gradGb+Sb
  dnbTbdt=(-gradQb+Pb+ohm+coulomb)*2./3
  dTbdt  =(dnbTbdt-Tb*dnbdt)/nb

  dEsdt=0
  DO ib=1,nbb
     dEsdt=dEsdt+Zb(ib)*nb(ib,:)*Gb(ib,:)
  END DO
  fact_Es=prefact_Es/etet/SQRT(s(:))/nb(2,:)/Ab(2)
  dEsdt=dEsdt*fact_Es

  nb=nb+dt*dnbdt
  Tb=Tb+dt*dTbdt
  Es=Es+dt*dEsdt
  
  CALL DERIVE(s(1:ns),nb(ib,1:ns),ns,4,dnbdpsi(ib,1:ns))
  CALL DERIVE(s(1:ns),Tb(ib,1:ns),ns,4,dTbdpsi(ib,1:ns))

END SUBROUTINE TRANSPORT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE DERIVE(s,p,ns,order,dpdpsi)
 
  IMPLICIT NONE
  
  INTEGER order,ns
  REAL*8 s(ns),p(ns),dpdpsi(ns)

  INTEGER is
  REAL*8 dp(ns),dpsi(ns)
  REAL*8 c21,c41,c42,c61,c62,c63,c81,c82,c83,c84
  REAl f20,f21,f22,f40,f41,f42,f43,f44,f60,f61,f62,f63,f64,f65,f66

  c21=+1/2.
  c41=+2/3.
  c42=-1/12.
  c61=+3/4.
  c62=-3/20.
  c63=+1/60.
  c81=+4/5.
  c82=-1/5.
  c83=+4/105.
  c84=-1/280.

  f20=-3/2.
  f21=+2.
  f22=-1/2.
  f40=-25/12.
  f41=+4.
  f42=-3.
  f43=+4/3.
  f44=-1/4.
  f60=-49/20.
  f61=+6.
  f62=-15/2.
  f63=+20/3.
  f64=-15/4.
  f65=+6/5.
  f66=-1/6.

  dpsi=s(ns)-s(ns-1)
  dp(1) =p(2)-p(1)
  dp(ns)=p(ns)-p(ns-1)
  dpsi(1)=s(2)-s(1)
  dpsi(ns)=s(ns)-s(ns-1)
  DO is=1,ns
     IF(is.NE.ns) dpsi(is)=(s(is+1)-s(is-1))/2.
     IF(order.EQ.6) THEN
        IF((is.GE.4).AND.(is.LE.(ns-3))) THEN
           dp(is)=-c63*p(is-3)-c62*p(is-2)-c61*p(is-1)+c61*p(is+1)+c62*p(is+2)+c63*p(is+3)
        ELSE IF (is.LT.4) THEN
           dp(is)=+f60*p(is)+f61*p(is+1)+f62*p(is+2)+f63*p(is+3)+f64*p(is+4)+f65*p(is+5)+f66*p(is+6)
        ELSE
           dp(is)=-f60*p(is)-f61*p(is-1)-f62*p(is-2)-f63*p(is-3)-f64*p(is-4)-f65*p(is-5)-f66*p(is-6)
        END IF
     ELSE IF(order.EQ.4) THEN
        IF((is.GE.3).AND.(is.LE.(ns-2))) THEN
           dp(is)=-c42*p(is-2)-c41*p(is-1)+c41*p(is+1)+c42*p(is+2)
        ELSE IF (is.LT.3) THEN
           dp(is)=+f40*p(is)+f41*p(is+1)+f42*p(is+2)+f43*p(is+3)+f44*p(is+4)
        ELSE
           dp(is)=-f40*p(is)-f41*p(is-1)-f42*p(is-2)-f43*p(is-3)-f44*p(is-4)
        END IF
     ELSE IF(ABS(order).EQ.2) THEN
        IF((is.GE.2).AND.(is.LE.(ns-1))) THEN
           dp(is)=-c21*p(is-1)+c21*p(is+1)
        ELSE IF(is.EQ.1) THEN
           IF(order.LT.0) THEN
              dp(is)=p(is+1)-p(is)
           ELSE
              dp(is)=+f20*p(is)+f21*p(is+1)+f22*p(is+2)
           END IF
        ELSE
           dp(is)=-f20*p(is)-f21*p(is-1)-f22*p(is-2)
        END IF
     ELSE IF((order.EQ.0).AND.(is.LT.ns)) THEN
        dp(is)=p(is+1)-p(is)
     END IF
  END DO
  dpdpsi=dp/dpsi
  
END SUBROUTINE DERIVE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

