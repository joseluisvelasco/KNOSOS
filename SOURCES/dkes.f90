!Read and proces results from DKES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE READ_DKES_TABLE(is)

!----------------------------------------------------------------------------------------------- 
!Read monoenergetic transport coefficients from DKES for surface s(is)
!-----------------------------------------------------------------------------------------------   

  USE GLOBAL
  IMPLICIT NONE
  !Input
  INTEGER is
  !Others
  CHARACTER*100 line,file,dir_efield(nefieldt),dir_cmul(ncmult)
  INTEGER iefield,icmul,iostat,nefield,ncmul
  REAL*8 D11p,D11m,dummy

  !Values of EFIELD and CMUL are preset (see global.f90)
  lcmult  =LOG(cmult)   !JL
  lefieldt=LOG(efieldt) !JL
  lefieldt(1)=-1000 !dummy value to avoid log(0.0)

  !Folders for different values of EFIELD
  dir_efield(1) ="omega_0e-0/"
  dir_efield(2) ="omega_1e-5/"
  dir_efield(3) ="omega_3e-5/"
  dir_efield(4) ="omega_1e-4/"
  dir_efield(5) ="omega_3e-4/"
  dir_efield(6) ="omega_1e-3/"
  dir_efield(7) ="omega_3e-3/"   !these lines may produce warnings depending on the compiler
  IF(nefieldt.GT.7) dir_efield(8) ="omega_1e-2/"
  IF(nefieldt.GT.8) dir_efield(9) ="omega_3e-2/"
  IF(nefieldt.GT.9) dir_efield(10)="omega_1e-1/"  
  !Subfolders of omega_?e-? for different values of CMUL
  dir_cmul(1) ="cl_3e+2/"
  dir_cmul(2) ="cl_1e+2/"
  dir_cmul(3) ="cl_3e+1/"
  dir_cmul(4) ="cl_1e+1/"
  dir_cmul(5) ="cl_3e-0/"
  dir_cmul(6) ="cl_1e-0/"
  dir_cmul(7) ="cl_3e-1/"
  dir_cmul(8) ="cl_1e-1/"
  dir_cmul(9) ="cl_3e-2/"
  dir_cmul(10)="cl_1e-2/"
  dir_cmul(11)="cl_3e-3/"
  dir_cmul(12)="cl_1e-3/"
  dir_cmul(13)="cl_3e-4/"
  dir_cmul(14)="cl_1e-4/"
  dir_cmul(15)="cl_3e-5/"
  dir_cmul(16)="cl_1e-5/"
  IF(ncmult.GT.16) dir_cmul(17)="cl_3e-6/"
  IF(ncmult.GT.17) dir_cmul(18)="cl_1e-6/"

  file=TRIM(DIRDB)//TRIM(DIRS(is))
  WRITE(1000+myrank,*) 'Reading DKES output in folder ',file 

  !Read data
  nefield=0  
  ncmul=0
  D11pla=1e10
  DO iefield=1,nefieldt
     DO icmul=1,ncmult
        file=TRIM(DIRDB)//TRIM(DIRS(is))//TRIM(dir_efield(iefield))//TRIM(dir_cmul(icmul))//"results.data"
        OPEN(unit=1,file=TRIM(file),action='read',iostat=iostat) 
        IF (iostat.EQ.0) THEN 
           WRITE(1000+myrank,*) 'Reading subfolder ',TRIM(dir_efield(iefield))//TRIM(dir_cmul(icmul))
           IF(icmul.EQ.1) nefield=nefield+1
           IF(iefield.EQ.1) ncmul=ncmul+1
           READ(1,*) line
           READ(1,*) line
           READ(1,*) dummy,dummy,dummy,dummy,D11m,D11p,&
                & dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
           CLOSE(1)
           !Logarithms are stored
           IF(D11p/D11m.GT.1E4.OR.D11p.LT.0) D11p=EXP(lD11dkes1(icmul-1,iefield))*cmult(icmul-1)/cmult(icmul)
           IF(iefield.EQ.1.AND.0.5*(D11m+D11p).LT.D11pla) D11pla=0.5*(D11m+D11p)
           WRITE(10000+myrank,'("-100 ",3(1pe13.5))') &
                & cmult(icmul),efieldt(iefield),0.5*(D11m+D11p)
           lD11dkes1(icmul,iefield)=LOG(0.5*(D11p+D11m)) !both values are the same, which means
           lD11dkes2(icmul,iefield)=LOG(0.5*(D11p+D11m)) !ignoring error bars in D11
        END IF
     END DO
  END DO
  IF(nefield.EQ.0) THEN
     WRITE(1000+myrank,*) 'No DKES output found'
  ELSE
     nefieldt=nefield
     ncmult=ncmul
     DKES_READ=.TRUE.
  END IF

END SUBROUTINE READ_DKES_TABLE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE INTERP_DATABASE(it,jv,Epsi,D11,knososdb)

!----------------------------------------------------------------------------------------------- 
!Calculate monoenergetic transport coefficient D11 from database, interpolating
!at cmul=nu(jv)/v(jv) and efield=Epsi*psip/v(jv). The database comes from KNOSOS if
!knososdb.EQ..TRUE. and from DKES otherwise
!-----------------------------------------------------------------------------------------------

  USE GLOBAL
  IMPLICIT NONE
  !Input
  LOGICAL knososdb
  INTEGER it,jv
  REAL*8 Epsi
  !Output
  REAL*8 D11
  !Others
  REAL*8 efield,cmul,cmag,lD11,lD11dkes(ncmult,nefieldt)

  IF(knososDB) THEN
     WRITE(1000+myrank,*) 'Interpolating KNOSOS coefficients'
     cmag=vmconst(jv)*Epsi/ABS(Epsi)
  ELSE
     WRITE(1000+myrank,*) 'Interpolating DKES coefficients'
  END IF
  efield=ABS(Epsi*psip/v(jv))

  !Parameters of the monoenergetic calculation
  cmul=nu(jv)/v(jv)/2 !careful with the definition of collision frequency

  IF(it.EQ.-1) THEN !irrelevant, as lD11dkes1=lD11dkes2, but could be used for errorbars
     lD11dkes(1:ncmult,1:nefieldt)=lD11dkes1(1:ncmult,1:nefieldt)
  ELSE
     lD11dkes(1:ncmult,1:nefieldt)=lD11dkes2(1:ncmult,1:nefieldt)
  END IF

  IF(knososDB) THEN     !Interpolate from a database calculated with KNOSOS
     IF(efield.GT.1E-10) THEN
        CALL BILAGRANGE(lcmult(1:ncmult),lefieldt(2:nefieldt),lD11tab(1:ncmult,2:nefieldt,1),&
             & ncmult,nefieldt,LOG(cmul),LOG(efield),lD11,1)
!        CALL TRILAGRANGE(lcmult(1:ncmult),lefieldt(2:nefieldt),lcmagt(1:ncmagt),lD11tab(1:ncmult,2:nefieldt,1:ncmagt),&
!             & ncmult,nefieldt,ngmact,LOG(cmul),LOG(efield),LOG(ABS(cmag))*cmag/ABS(cmag),lD11,1)
     ELSE
        CALL LAGRANGE(lcmult(1:ncmult),lD11tab(1:ncmult,1,1),&
             & ncmult,LOG(cmul),lD11,1)
!        CALL BILAGRANGE(lcmult(1:ncmult),lcmagt(1:nmagt),lD11tab(1:ncmult,1,1:ncmagt),&
!             & ncmult,ncmagt,LOG(cmul),LOG(ABS(cmag))*cmag/ABS(cmag),lD11,1)
     END IF
  ELSE                  !Interpolate from DKES
     IF(efield.GT.0.1*efieldt(2)) THEN
        CALL BILAGRANGE(lcmult(1:ncmult),lefieldt(2:nefieldt),lD11dkes(1:ncmult,2:nefieldt),&
             & ncmult,nefieldt-1,LOG(cmul),LOG(efield),lD11,1)
     ELSE
        !If the radial electric field is close to zero, use only coefficients of efield=0
        CALL LAGRANGE(lcmult(1:ncmult),lD11dkes(1:ncmult,1),&
             & ncmult,LOG(cmul),lD11,1)     
     END IF
    
  END IF

  !After logarithmic interpolation, take exponential
  D11=EXP(lD11)

  WRITE(200+myrank,'(2(1pe13.5)," NaN NaN ",2(1pe13.5)," &
          & NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN")') &
          & cmul,efield,D11,D11
  WRITE(10000+myrank,'(I3,5(1pe13.5))') it,cmul,efield,D11,&
        & weight(jv)/fdkes(jv),weight(jv)/fdkes(jv)*v(jv)*v(jv)

  !Normalize
  D11=D11/fdkes(jv)
  
END SUBROUTINE INTERP_DATABASE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
