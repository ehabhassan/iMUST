MODULE SYSSOLVER
 USE MATHSUBS
 USE DESET
 USE DOPERATORS

!**************************************************************
!******************* BOUNDARY CONDITIONS **********************
!**************************************************************
! nfld         => NUMBER OF ORDINARY DIFFERENTIAL EQUATIONS
! fDIM1        => FIELDS FIRST DIMENSION
! fDIM2        => FIELDS SECOND DIMENSION
! fDIM3        => FIELDS THIRD DIMENSION
 INTEGER(4)              :: nfld
 INTEGER(4)              :: fDIM1,fDIM2,fDIM3

!**************************************************************
!************************* SPACE SOLVER ***********************
!**************************************************************
! RSPACESOLVER     => SOLVING THE DIFFERENTIAL EQUATIONS ON THE R-SPACE
! KSPACESOLVER     => SOLVING THE DIFFERENTIAL EQUATIONS ON THE K-SPACE
 LOGICAL                 :: RSPACESOLVER         =.FALSE.
 LOGICAL                 :: KSPACESOLVER         =.TRUE.

!***************************************************************
!**************** EIGENVALUES AND EIGENVECTORS *****************
!***************************************************************
 REAL(4),   ALLOCATABLE  :: eiKx1D(:),eiKy1D(:),eiKz1D(:)
 INTEGER(4),ALLOCATABLE  :: KxGMAX3D(:,:,:),KyGMAX3D(:,:,:),KzGMAX3D(:,:,:)
 COMPLEX(4),ALLOCATABLE  :: EIGENVAL7D(:,:,:,:,:,:,:)
 COMPLEX(4),ALLOCATABLE  :: EIGENVEC7D(:,:,:,:,:,:,:)



!***************************************************************
!                MODULE PROCEDURAL INTERFACES
!***************************************************************

 INTERFACE DESOLVER
   MODULE PROCEDURE r4r4DESOLVER1D,r4r8DESOLVER1D,r4c8DESOLVER1D,r8r8DESOLVER1D,r8c8DESOLVER1D, &
 &                  r4r4DESOLVER2D,r4r8DESOLVER2D,r4c8DESOLVER2D,r8r8DESOLVER2D,r8c8DESOLVER2D, &
 &                  r4r4DESOLVER3D,r4r8DESOLVER3D,r4c8DESOLVER3D,r8r8DESOLVER3D,r8c8DESOLVER3D
 END INTERFACE DESOLVER

 INTERFACE SYSEIGENVALS
   MODULE PROCEDURE r4SYSEIGENVALS,r8SYSEIGENVALS
 END INTERFACE SYSEIGENVALS

 INTERFACE SYSEIGENVECS
   MODULE PROCEDURE r4SYSEIGENVECS,r8SYSEIGENVECS
 END INTERFACE SYSEIGENVECS

 INTERFACE SYSEIGENVV
   MODULE PROCEDURE r4SYSEIGENVV,r8SYSEIGENVV
 END INTERFACE SYSEIGENVV

 INTERFACE SYSEIGENFLDS
   MODULE PROCEDURE r4r8c8SYSEIGENFLDS,r4c8r8SYSEIGENFLDS,r8r8c8SYSEIGENFLDS,r8c8r8SYSEIGENFLDS
 END INTERFACE SYSEIGENFLDS

 INTERFACE EPKEIGENSEVAL
   MODULE PROCEDURE R4EPKEIGENSEVAL,R8EPKEIGENSEVAL
 END INTERFACE EPKEIGENSEVAL


 CONTAINS


!****************************************************************
!********* DIFFERENTIAL EQUATIONS (DE) SPECTRAL SOLVER **********
!****************************************************************
 SUBROUTINE r4r4DESOLVER1D(ODESETPIN,LX1D,flds,t,tstp)
  implicit none
!********** DECLARATION OF INPUT VARIABLES AND ARRAYS ***********
  REAL(4),   INTENT(INOUT) :: flds(:,:)
  REAL(4),   INTENT(INOUT) :: tstp
  REAL(4),   INTENT(IN)    :: t
  REAL(4),   INTENT(IN)    :: LX1D(:)
  INTEGER(4),INTENT(IN)    :: ODESETPIN
!**** DECLARATION OF LOCAL VARIABLES WITHOUT INITIAL VALUES *****
  REAL(8)                  :: LocT,LocTSTP
  INTEGER(4)               :: MATDIMS(2)
  INTEGER(4)               :: neq,ind,nw
  INTEGER(4)               :: I,J,iD1,iD2
  INTEGER(4), SAVE         :: ACCESSID = 1
  INTEGER(4), SAVE         :: ODESETID = 1
!********** DECLARATION OF VARIABLES WITH INITIAL VALUES ********
  REAL(8)                  :: c24(24)=0.0d0
  REAL(8)                  :: tol=1.0d-3
  REAL(8)                  :: hmax,hmag
  REAL(8)                  :: sci=0.0d0
  REAL(8)                  :: sig=0.03d0
  REAL(8)                  :: tend

!*************** DECLARATION OF ALLOCATABLE ARRAYS **************
  REAL(8),   ALLOCATABLE   :: rwk(:,:,:,:)
  REAL(8),   ALLOCATABLE   :: Locflds(:,:)
  COMPLEX(8),ALLOCATABLE   :: kflds(:,:)

  MATDIMS=SHAPE(flds)

  ALLOCATE(Locflds(MATDIMS(1),MATDIMS(2)))
  LocT=REAL(t,8)
  LocTSTP=REAL(tstp,8)

  if(sci.eq.0d0)then
    hmax=2d0
  elseif(sci.ne.0d0)then
    hmax=2d0/abs(sci)
  endif
  hmag=hmax*tol**(1d0/6d0)

  c24(4)=hmag
  c24(5)=sci
  c24(6)=hmax
  ind=2

  tend=LocT+LocTSTP

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
  DO iD1=1,MATDIMS(1)
   DO iD2=1,MATDIMS(2)
    Locflds(iD1,iD2)=REAL(flds(iD1,iD2),8)
   END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

  IF(RSPACESOLVER) THEN
    fDIM1=MATDIMS(1)
    nfld =MATDIMS(2)

    neq=1*MATDIMS(1)*MATDIMS(2)
    nw=neq

    IF(ALLOCATED(rwk).EQV..FALSE.) ALLOCATE(rwk(1,MATDIMS(1),MATDIMS(2),9))
    CALL dverk(ODESETID,neq,r8fcn1D,LocT,Locflds,tend,tol,ind,c24,nw,rwk)
    hmag=c24(14)
  ELSE IF(KSPACESOLVER) THEN
    fDIM1=MATDIMS(1)/2+1
    nfld =MATDIMS(2)

    neq=2*(MATDIMS(1)/2+1)*MATDIMS(2)
    nw=neq

    IF(ALLOCATED(rwk).EQV..FALSE.)   ALLOCATE(rwk(2,MATDIMS(1)/2+1,MATDIMS(2),9))
    IF(ALLOCATED(kflds).EQV..FALSE.) ALLOCATE(kflds(MATDIMS(1)/2+1,MATDIMS(2)))
    DO I=1,MATDIMS(2)
      CALL FFTW(Locflds(:,I),kflds(:,I))
    END DO
    CALL dverk(ODESETID,neq,c8fcn1D,LocT,kflds,tend,tol,ind,c24,nw,rwk)
    hmag=c24(14)
    DO I=1,MATDIMS(2)
      CALL iFFTW(kflds(:,I),Locflds(:,I))
    END DO
  END IF

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
  DO iD1=1,MATDIMS(1)
   DO iD2=1,MATDIMS(2)
    flds(iD1,iD2)=REAL(Locflds(iD1,iD2),4)
   END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

 END SUBROUTINE r4r4DESOLVER1D

 SUBROUTINE r4r8DESOLVER1D(ODESETPIN,LX1D,flds,t,tstp)
  IMPLICIT NONE
!********** DECLARATION OF INCOMING VARIABLES AND ARRAYS*********
  REAL(4),   INTENT(IN)    :: LX1D(:)
  REAL(8),   INTENT(INOUT) :: flds(:,:)
  REAL(4),   INTENT(INOUT) :: tstp
  REAL(4),   INTENT(IN)    :: t
  INTEGER(4),INTENT(IN)    :: ODESETPIN

  CALL DESOLVER(ODESETPIN,1.0d0*LX1D,flds,t,tstp)

  RETURN
 END SUBROUTINE r4r8DESOLVER1D

 SUBROUTINE r8r8DESOLVER1D(ODESETPIN,LX1D,flds,t,tstp)
  IMPLICIT NONE
!********** DECLARATION OF INCOMING VARIABLES AND ARRAYS*********
  REAL(8),   INTENT(IN)    :: LX1D(:)
  REAL(8),   INTENT(INOUT) :: flds(:,:)
  REAL(4),   INTENT(INOUT) :: tstp
  REAL(4),   INTENT(IN)    :: t
  INTEGER(4),INTENT(IN)    :: ODESETPIN
!**** DECLARATION OF LOCAL VARIABLES WITHOUT INITIAL VALUES *****
  REAL(8)                  :: LocT,LocTSTP
  INTEGER(4)               :: MATDIMS(2)
  INTEGER(4)               :: neq,ind,nw
  INTEGER(4)               :: I,J
  INTEGER(4), SAVE         :: ACCESSID = 1
  INTEGER(4), SAVE         :: ODESETID = 1
!********** DECLARATION OF VARIABLES WITH INITIAL VALUES ********
  REAL(8)                  :: c24(24)=0.0d0
  REAL(8)                  :: tol=1.0d-3
  REAL(8)                  :: hmax,hmag
  REAL(8)                  :: sci=0.0d0
  REAL(8)                  :: sig=0.03d0
  REAL(8)                  :: tend

!*************** DECLARATION OF ALLOCATABLE ARRAYS **************
  REAL(8),   ALLOCATABLE   :: rwk(:,:,:,:)
  COMPLEX(8),ALLOCATABLE   :: kflds(:,:)

  MATDIMS=SHAPE(flds)

  if(sci.eq.0d0)then
    hmax=2d0
  elseif(sci.ne.0d0)then
    hmax=2d0/abs(sci)
  endif
  hmag=hmax*tol**(1d0/6d0)

  c24(4)=hmag
  c24(5)=sci
  c24(6)=hmax
  ind=2

  tend=LocT+LocTSTP

  IF(RSPACESOLVER) THEN
    fDIM1=MATDIMS(1)
    nfld =MATDIMS(2)

    neq=1*MATDIMS(1)*MATDIMS(2)
    nw=neq

    IF(ALLOCATED(rwk).EQV..FALSE.) ALLOCATE(rwk(1,MATDIMS(1),MATDIMS(2),9))
    CALL dverk(ODESETID,neq,r8fcn1D,LocT,flds,tend,tol,ind,c24,nw,rwk)
    hmag=c24(14)
  ELSE IF(KSPACESOLVER) THEN
    fDIM1=MATDIMS(1)/2+1
    nfld =MATDIMS(2)

    neq=2*(MATDIMS(1)/2+1)*MATDIMS(2)
    nw=neq

    IF(ALLOCATED(rwk).EQV..FALSE.)   ALLOCATE(rwk(2,MATDIMS(1)/2+1,MATDIMS(2),9))
    IF(ALLOCATED(kflds).EQV..FALSE.) ALLOCATE(kflds(MATDIMS(1)/2+1,MATDIMS(2)))
    DO I=1,MATDIMS(2)
      CALL FFTW(flds(:,I),kflds(:,I))
    END DO
    CALL dverk(ODESETID,neq,c8fcn1D,LocT,kflds,tend,tol,ind,c24,nw,rwk)
    hmag=c24(14)
    DO I=1,MATDIMS(2)
      CALL iFFTW(kflds(:,I),flds(:,I))
    END DO
  END IF

 END SUBROUTINE r8r8DESOLVER1D

 SUBROUTINE r4c8DESOLVER1D(ODESETPIN,LX1D,flds,t,tstp)
  IMPLICIT NONE
!********** DECLARATION OF INCOMING VARIABLES AND ARRAYS*********
  REAL(4),    INTENT(IN)    :: LX1D(:)
  COMPLEX(8), INTENT(INOUT) :: flds(:,:)
  REAL(4),    INTENT(INOUT) :: tstp
  REAL(4),    INTENT(IN)    :: t
  INTEGER(4), INTENT(IN)    :: ODESETPIN

  CALL DESOLVER(ODESETPIN,1.0d0*LX1D,flds,t,tstp)

  RETURN
 END SUBROUTINE r4c8DESOLVER1D

 SUBROUTINE r8c8DESOLVER1D(ODESETPIN,LX1D,flds,t,tstp)
  IMPLICIT NONE
!********** DECLARATION OF INCOMING VARIABLES AND ARRAYS*********
  REAL(8),    INTENT(IN)    :: LX1D(:)
  COMPLEX(8), INTENT(INOUT) :: flds(:,:)
  REAL(4),    INTENT(INOUT) :: tstp
  REAL(4),    INTENT(IN)    :: t
  INTEGER(4), INTENT(IN)    :: ODESETPIN
!**** DECLARATION OF LOCAL VARIABLES WITHOUT INITIAL VALUES *****
  REAL(8)                   :: LocT,LocTSTP
  INTEGER(4)                :: MATDIMS(2)
  INTEGER(4)                :: neq,ind,nw
  INTEGER(4)                :: I,J
  INTEGER(4), SAVE          :: ACCESSID = 1
  INTEGER(4), SAVE          :: ODESETID = 1
!********** DECLARATION OF VARIABLES WITH INITIAL VALUES ********
  REAL(8)                   :: c24(24)=0.0d0
  REAL(8)                   :: tol=1.0d-3
  REAL(8)                   :: hmax,hmag
  REAL(8)                   :: sci=0.0d0
  REAL(8)                   :: sig=0.03d0
  REAL(8)                   :: tend

!*************** DECLARATION OF ALLOCATABLE ARRAYS **************
  REAL(8),   ALLOCATABLE    :: rwk(:,:,:,:)

  MATDIMS=SHAPE(flds)
  fDIM1=MATDIMS(1)
  nfld =MATDIMS(2)

  if(sci.eq.0d0)then
    hmax=2d0
  elseif(sci.ne.0d0)then
    hmax=2d0/abs(sci)
  endif
  hmag=hmax*tol**(1d0/6d0)

  c24(4)=hmag
  c24(5)=sci
  c24(6)=hmax
  ind=2

  tend=LocT+LocTSTP

  neq=2*MATDIMS(1)*MATDIMS(2)
  nw=neq
  IF(ALLOCATED(rwk).EQV..FALSE.) ALLOCATE(rwk(2,MATDIMS(1),MATDIMS(2),9))
  CALL dverk(ODESETID,neq,c8fcn1D,LocT,flds,tend,tol,ind,c24,nw,rwk)
  hmag=c24(14)

  RETURN
 END SUBROUTINE r8c8DESOLVER1D

 SUBROUTINE r4r4DESOLVER2D(ODESETPIN,LX2D,LY2D,flds,t,tstp)
  implicit none
!********** DECLARATION OF INPUT VARIABLES AND ARRAYS ***********
  REAL(4),   INTENT(INOUT) :: flds(:,:,:)
  REAL(4),   INTENT(INOUT) :: tstp
  REAL(4),   INTENT(IN)    :: t
  REAL(4),   INTENT(IN)    :: LX2D(:,:),LY2D(:,:)
  INTEGER(4),INTENT(IN)    :: ODESETPIN
!**** DECLARATION OF LOCAL VARIABLES WITHOUT INITIAL VALUES *****
  REAL(8)                  :: LocT,LocTSTP
  INTEGER(4)               :: MATDIMS(3)
  INTEGER(4)               :: neq,ind,nw
  INTEGER(4)               :: I,J,iD1,iD2,iD3
  INTEGER(4), SAVE         :: ACCESSID = 1
  INTEGER(4), SAVE         :: ODESETID = 1
!********** DECLARATION OF VARIABLES WITH INITIAL VALUES ********
  REAL(8)                  :: c24(24)=0.0d0
  REAL(8)                  :: tol=1.0d-3
  REAL(8)                  :: hmax,hmag
  REAL(8)                  :: sci=0.0d0
  REAL(8)                  :: sig=0.03d0
  REAL(8)                  :: tend

!*************** DECLARATION OF ALLOCATABLE ARRAYS **************
  REAL(8),   ALLOCATABLE   :: rwk(:,:,:,:,:)
  COMPLEX(8),ALLOCATABLE   :: kflds(:,:,:)
  REAL(8),   ALLOCATABLE   :: Locflds(:,:,:)

  MATDIMS=SHAPE(flds)

  ALLOCATE(Locflds(MATDIMS(1),MATDIMS(2),MATDIMS(3)))
  LocT=REAL(t,8)
  LocTSTP=REAL(tstp,8)

  if(sci.eq.0d0)then
    hmax=2d0
  elseif(sci.ne.0d0)then
    hmax=2d0/abs(sci)
  endif
  hmag=hmax*tol**(1d0/6d0)

  c24(4)=hmag
  c24(5)=sci
  c24(6)=hmax
  ind=2

  tend=LocT+LocTSTP

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
  DO iD1=1,MATDIMS(1)
   DO iD2=1,MATDIMS(2)
    DO iD3=1,MATDIMS(3)
     Locflds(iD1,iD2,iD3)=REAL(flds(iD1,iD2,iD3),8)
    END DO
   END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

  IF(RSPACESOLVER) THEN
    fDIM1=MATDIMS(1)
    fDIM2=MATDIMS(2)
    nfld =MATDIMS(3)

    neq=1*MATDIMS(1)*MATDIMS(2)*MATDIMS(3)
    nw=neq

    IF(ALLOCATED(rwk).EQV..FALSE.) ALLOCATE(rwk(1,MATDIMS(1),MATDIMS(2),MATDIMS(3),9))
    CALL dverk(ODESETID,neq,r8fcn2D,LocT,Locflds,tend,tol,ind,c24,nw,rwk)
    hmag=c24(14)
  ELSE IF(KSPACESOLVER) THEN
    fDIM1=MATDIMS(1)/2+1
    fDIM2=MATDIMS(2)
    nfld =MATDIMS(3)

    neq=2*(MATDIMS(1)/2+1)*MATDIMS(2)*MATDIMS(3)
    nw=neq

    IF(ALLOCATED(rwk).EQV..FALSE.)   ALLOCATE(rwk(2,MATDIMS(1)/2+1,MATDIMS(2),MATDIMS(3),9))
    IF(ALLOCATED(kflds).EQV..FALSE.) ALLOCATE(kflds(MATDIMS(1)/2+1,MATDIMS(2),MATDIMS(3)))
    DO I=1,MATDIMS(3)
      CALL FFTW(Locflds(:,:,I),kflds(:,:,I))
    END DO
    CALL dverk(ODESETID,neq,c8fcn2D,LocT,kflds,tend,tol,ind,c24,nw,rwk)
    hmag=c24(14)
    DO I=1,MATDIMS(3)
      CALL iFFTW(kflds(:,:,I),Locflds(:,:,I))
    END DO
  END IF

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
  DO iD1=1,MATDIMS(1)
   DO iD2=1,MATDIMS(2)
    DO iD3=1,MATDIMS(3)
     flds(iD1,iD2,iD3)=REAL(Locflds(iD1,iD2,iD3),4)
    END DO
   END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

 END SUBROUTINE r4r4DESOLVER2D

 SUBROUTINE r4r8DESOLVER2D(ODESETPIN,LX2D,LY2D,flds,t,tstp)
  IMPLICIT NONE
!********** DECLARATION OF INCOMING VARIABLES AND ARRAYS*********
  REAL(4),   INTENT(IN)    :: LX2D(:,:),LY2D(:,:)
  REAL(8),   INTENT(INOUT) :: flds(:,:,:)
  REAL(4),   INTENT(INOUT) :: tstp
  REAL(4),   INTENT(IN)    :: t
  INTEGER(4),INTENT(IN)    :: ODESETPIN

  CALL DESOLVER(ODESETPIN,1.0d0*LX2D,1.0d0*LY2D,flds,t,tstp)

  RETURN
 END SUBROUTINE r4r8DESOLVER2D

 SUBROUTINE r8r8DESOLVER2D(ODESETPIN,LX2D,LY2D,flds,t,tstp)
  IMPLICIT NONE
!********** DECLARATION OF INCOMING VARIABLES AND ARRAYS*********
  REAL(8),   INTENT(IN)    :: LX2D(:,:),LY2D(:,:)
  REAL(8),   INTENT(INOUT) :: flds(:,:,:)
  REAL(4),   INTENT(INOUT) :: tstp
  REAL(4),   INTENT(IN)    :: t
  INTEGER(4),INTENT(IN)    :: ODESETPIN
!**** DECLARATION OF LOCAL VARIABLES WITHOUT INITIAL VALUES *****
  REAL(8)                  :: LocT,LocTSTP
  INTEGER(4)               :: MATDIMS(3)
  INTEGER(4)               :: neq,ind,nw
  INTEGER(4)               :: I,J
  INTEGER(4), SAVE         :: ACCESSID = 1
  INTEGER(4), SAVE         :: ODESETID = 1
!********** DECLARATION OF VARIABLES WITH INITIAL VALUES ********
  REAL(8)                  :: c24(24)=0.0d0
  REAL(8)                  :: tol=1.0d-3
  REAL(8)                  :: hmax,hmag
  REAL(8)                  :: sci=0.0d0
  REAL(8)                  :: sig=0.03d0
  REAL(8)                  :: tend

!*************** DECLARATION OF ALLOCATABLE ARRAYS **************
  REAL(8),   ALLOCATABLE   :: rwk(:,:,:,:,:)
  COMPLEX(8),ALLOCATABLE   :: kflds(:,:,:)

  MATDIMS=SHAPE(flds)

  if(sci.eq.0d0)then
    hmax=2d0
  elseif(sci.ne.0d0)then
    hmax=2d0/abs(sci)
  endif
  hmag=hmax*tol**(1d0/6d0)

  c24(4)=hmag
  c24(5)=sci
  c24(6)=hmax
  ind=2

  LocT=REAL(t,8)
  LocTSTP=REAL(tstp,8)
  tend=LocT+LocTSTP

  IF(RSPACESOLVER) THEN
    fDIM1=MATDIMS(1)
    fDIM2=MATDIMS(2)
    nfld =MATDIMS(3)

    neq=1*MATDIMS(1)*MATDIMS(2)*MATDIMS(3)
    nw=neq

    IF(ALLOCATED(rwk).EQV..FALSE.) ALLOCATE(rwk(1,MATDIMS(1),MATDIMS(2),MATDIMS(3),9))
    CALL dverk(ODESETPIN,neq,r8fcn2D,LocT,flds,tend,tol,ind,c24,nw,rwk)
    hmag=c24(14)
  ELSE IF(KSPACESOLVER) THEN
    fDIM1=MATDIMS(1)/2+1
    fDIM2=MATDIMS(2)
    nfld =MATDIMS(3)

    neq=2*(MATDIMS(1)/2+1)*MATDIMS(2)*MATDIMS(3)
    nw=neq

    IF(ALLOCATED(rwk).EQV..FALSE.)   ALLOCATE(rwk(2,MATDIMS(1)/2+1,MATDIMS(2),MATDIMS(3),9))
    IF(ALLOCATED(kflds).EQV..FALSE.) ALLOCATE(kflds(MATDIMS(1)/2+1,MATDIMS(2),MATDIMS(3)))
    DO I=1,MATDIMS(3)
      CALL FFTW(flds(:,:,I),kflds(:,:,I))
    END DO
    CALL dverk(ODESETPIN,neq,c8fcn2D,LocT,kflds,tend,tol,ind,c24,nw,rwk)
    hmag=c24(14)
    DO I=1,MATDIMS(3)
      CALL iFFTW(kflds(:,:,I),flds(:,:,I))
    END DO
  END IF

 END SUBROUTINE r8r8DESOLVER2D

 SUBROUTINE r4c8DESOLVER2D(ODESETPIN,LX2D,LY2D,flds,t,tstp)
  IMPLICIT NONE
!********** DECLARATION OF INCOMING VARIABLES AND ARRAYS*********
  REAL(4),    INTENT(IN)    :: LX2D(:,:),LY2D(:,:)
  COMPLEX(8), INTENT(INOUT) :: flds(:,:,:)
  REAL(4),    INTENT(INOUT) :: tstp
  REAL(4),    INTENT(IN)    :: t
  INTEGER(4), INTENT(IN)    :: ODESETPIN

  CALL DESOLVER(ODESETPIN,1.0d0*LX2D,1.0d0*LY2D,flds,t,tstp)

  RETURN
 END SUBROUTINE r4c8DESOLVER2D

 SUBROUTINE r8c8DESOLVER2D(ODESETPIN,LX2D,LY2D,flds,t,tstp)
  IMPLICIT NONE
!********** DECLARATION OF INCOMING VARIABLES AND ARRAYS*********
  REAL(8),    INTENT(IN)    :: LX2D(:,:),LY2D(:,:)
  COMPLEX(8), INTENT(INOUT) :: flds(:,:,:)
  REAL(4),    INTENT(INOUT) :: tstp
  REAL(4),    INTENT(IN)    :: t
  INTEGER(4), INTENT(IN)    :: ODESETPIN
!**** DECLARATION OF LOCAL VARIABLES WITHOUT INITIAL VALUES *****
  REAL(8)                   :: LocT,LocTSTP
  INTEGER(4)                :: MATDIMS(3)
  INTEGER(4)                :: neq,ind,nw
  INTEGER(4)                :: I,J
  INTEGER(4), SAVE          :: ACCESSID = 1
  INTEGER(4), SAVE          :: ODESETID = 1
!********** DECLARATION OF VARIABLES WITH INITIAL VALUES ********
  REAL(8)                   :: c24(24)=0.0d0
  REAL(8)                   :: tol=1.0d-3
  REAL(8)                   :: hmax,hmag
  REAL(8)                   :: sci=0.0d0
  REAL(8)                   :: sig=0.03d0
  REAL(8)                   :: tend

!*************** DECLARATION OF ALLOCATABLE ARRAYS **************
  REAL(8),   ALLOCATABLE    :: rwk(:,:,:,:,:)

  MATDIMS=SHAPE(flds)
  fDIM1=MATDIMS(1)
  fDIM2=MATDIMS(2)
  nfld =MATDIMS(3)

  if(sci.eq.0d0)then
    hmax=2d0
  elseif(sci.ne.0d0)then
    hmax=2d0/abs(sci)
  endif
  hmag=hmax*tol**(1d0/6d0)

  c24(4)=hmag
  c24(5)=sci
  c24(6)=hmax
  ind=2
 !ind=1

  LocT=REAL(t,8)
  LocTSTP=REAL(tstp,8)
  tend=LocT+LocTSTP

  neq=2*MATDIMS(1)*MATDIMS(2)*MATDIMS(3)
  nw=neq
  IF(ALLOCATED(rwk).EQV..FALSE.) ALLOCATE(rwk(2,MATDIMS(1),MATDIMS(2),MATDIMS(3),9))
  CALL dverk(ODESETID,neq,c8fcn2D,LocT,flds,tend,tol,ind,c24,nw,rwk)
  hmag=c24(14)

  RETURN
 END SUBROUTINE r8c8DESOLVER2D


 SUBROUTINE r4r4DESOLVER3D(ODESETPIN,LX3D,LY3D,LZ3D,flds,t,tstp)
  IMPLICIT NONE
!********** DECLARATION OF INPUT VARIABLES AND ARRAYS ***********
  REAL(4), INTENT(INOUT) :: flds(:,:,:,:)
  REAL(4), INTENT(INOUT) :: tstp
  REAL(4), INTENT(IN)    :: t
  REAL(4), INTENT(IN)    :: LX3D(:,:,:),LY3D(:,:,:),LZ3D(:,:,:)
  INTEGER(4),INTENT(IN)  :: ODESETPIN
!**** DECLARATION OF LOCAL VARIABLES WITHOUT INITIAL VALUES *****
  REAL(8)                :: LocT,LocTSTP
  INTEGER(4)             :: MATDIMS(4)
  INTEGER(4)             :: neq,ind,nw
  INTEGER(4)             :: I,J,iD1,iD2,iD3,iD4
  INTEGER(4), SAVE       :: ACCESSID = 1
  INTEGER(4), SAVE       :: ODESETID = 1
!********** DECLARATION OF VARIABLES WITH INITIAL VALUES ********
  REAL(8)                :: c24(24)=0.0d0
  REAL(8)                :: tol=1.0d-3
  REAL(8)                :: hmax,hmag
  REAL(8)                :: sci=0.0d0
  REAL(8)                :: sig=0.03d0
  REAL(8)                :: tend

!*************** DECLARATION OF ALLOCATABLE ARRAYS **************
  REAL(8),   ALLOCATABLE :: rwk(:,:,:,:,:,:)
  COMPLEX(8),ALLOCATABLE :: kflds(:,:,:,:)
  REAL(8),   ALLOCATABLE :: Locflds(:,:,:,:)

  MATDIMS=SHAPE(flds)

  ALLOCATE(Locflds(MATDIMS(1),MATDIMS(2),MATDIMS(3),MATDIMS(4)))
  LocT=REAL(t,8)
  LocTSTP=REAL(tstp,8)

  if(sci.eq.0d0)then
    hmax=2d0
  elseif(sci.ne.0d0)then
    hmax=2d0/abs(sci)
  endif
  hmag=hmax*tol**(1d0/6d0)

  c24(4)=hmag
  c24(5)=sci
  c24(6)=hmax
  ind=2

  tend=LocT+LocTSTP

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
  DO iD1=1,MATDIMS(1)
   DO iD2=1,MATDIMS(2)
    DO iD3=1,MATDIMS(3)
     DO iD4=1,MATDIMS(4)
      Locflds(iD1,iD2,iD3,iD4)=REAL(flds(iD1,iD2,iD3,iD4),8)
     END DO
    END DO
   END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

  IF(RSPACESOLVER) THEN
    fDIM1=MATDIMS(1)
    fDIM2=MATDIMS(2)
    fDIM3=MATDIMS(3)
    nfld =MATDIMS(4)

    neq=1*MATDIMS(1)*MATDIMS(2)*MATDIMS(3)*MATDIMS(4)
    nw=neq

    IF(ALLOCATED(rwk).EQV..FALSE.) ALLOCATE(rwk(1,MATDIMS(1),MATDIMS(2),MATDIMS(3),MATDIMS(4),9))
    CALL dverk(ODESETID,neq,r8fcn3D,LocT,Locflds,tend,tol,ind,c24,nw,rwk)
    hmag=c24(14)
  ELSE IF(KSPACESOLVER) THEN
    fDIM1=MATDIMS(1)/2+1
    fDIM2=MATDIMS(2)
    fDIM3=MATDIMS(3)
    nfld =MATDIMS(4)

    neq=2*(MATDIMS(1)/2+1)*MATDIMS(2)*MATDIMS(3)*MATDIMS(4)
    nw=neq

    IF(ALLOCATED(rwk).EQV..FALSE.)   ALLOCATE(rwk(2,MATDIMS(1)/2+1,MATDIMS(2),MATDIMS(3),MATDIMS(4),9))
    IF(ALLOCATED(kflds).EQV..FALSE.) ALLOCATE(kflds(MATDIMS(1)/2+1,MATDIMS(2),MATDIMS(3),MATDIMS(4)))
    DO I=1,MATDIMS(4)
      CALL FFTW(Locflds(:,:,1,I),kflds(:,:,1,I))
    END DO
    CALL dverk(ODESETID,neq,c8fcn3D,LocT,kflds,tend,tol,ind,c24,nw,rwk)
    hmag=c24(14)
    DO I=1,MATDIMS(4)
      CALL iFFTW(kflds(:,:,:,I),Locflds(:,:,:,I))
    END DO
  END IF

!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
  DO iD1=1,MATDIMS(1)
   DO iD2=1,MATDIMS(2)
    DO iD3=1,MATDIMS(3)
     DO iD4=1,MATDIMS(4)
      flds(iD1,iD2,iD3,iD4)=REAL(Locflds(iD1,iD2,iD3,iD4),4)
     END DO
    END DO
   END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

 END SUBROUTINE r4r4DESOLVER3D

 SUBROUTINE r4r8DESOLVER3D(ODESETPIN,LX3D,LY3D,LZ3D,flds,t,tstp)
  implicit none
!********** DECLARATION OF INCOMING VARIABLES AND ARRAYS*********
  REAL(8), INTENT(INOUT) :: flds(:,:,:,:)
  REAL(4), INTENT(INOUT) :: tstp
  REAL(4), INTENT(IN)    :: t
  REAL(4), INTENT(IN)    :: LX3D(:,:,:),LY3D(:,:,:),LZ3D(:,:,:)
  INTEGER(4),INTENT(IN)  :: ODESETPIN

  CALL DESOLVER(ODESETPIN,1.0d0*LX3D,1.0d0*LY3D,1.0d0*LZ3D,flds,t,tstp)

  RETURN
 END SUBROUTINE r4r8DESOLVER3D

 SUBROUTINE r8r8DESOLVER3D(ODESETPIN,LX3D,LY3D,LZ3D,flds,t,tstp)
  implicit none
!********** DECLARATION OF INCOMING VARIABLES AND ARRAYS*********
  REAL(8), INTENT(INOUT) :: flds(:,:,:,:)
  REAL(4), INTENT(INOUT) :: tstp
  REAL(4), INTENT(IN)    :: t
  REAL(8), INTENT(IN)    :: LX3D(:,:,:),LY3D(:,:,:),LZ3D(:,:,:)
  INTEGER(4),INTENT(IN)  :: ODESETPIN
!**** DECLARATION OF LOCAL VARIABLES WITHOUT INITIAL VALUES *****
  REAL(8)                :: LocT,LocTSTP
  INTEGER(4)             :: MATDIMS(4)
  INTEGER(4)             :: neq,ind,nw
  INTEGER(4)             :: I,J
  INTEGER(4), SAVE       :: ACCESSID = 1
  INTEGER(4), SAVE       :: ODESETID = 1
!********** DECLARATION OF VARIABLES WITH INITIAL VALUES ********
  REAL(8)                :: c24(24)=0.0d0
  REAL(8)                :: tol=1.0d-3
  REAL(8)                :: hmax,hmag
  REAL(8)                :: sci=0.0d0
  REAL(8)                :: sig=0.03d0
  REAL(8)                :: tend

!*************** DECLARATION OF ALLOCATABLE ARRAYS **************
  REAL(8),   ALLOCATABLE :: rwk(:,:,:,:,:,:)
  COMPLEX(8),ALLOCATABLE :: kflds(:,:,:,:)

  MATDIMS=SHAPE(flds)

  if(sci.eq.0d0)then
    hmax=2d0
  elseif(sci.ne.0d0)then
    hmax=2d0/abs(sci)
  endif
  hmag=hmax*tol**(1d0/6d0)

  c24(4)=hmag
  c24(5)=sci
  c24(6)=hmax
  ind=2

  tend=LocT+LocTSTP

  IF(RSPACESOLVER) THEN
    fDIM1=MATDIMS(1)
    fDIM2=MATDIMS(2)
    fDIM3=MATDIMS(3)
    nfld =MATDIMS(4)

    neq=1*MATDIMS(1)*MATDIMS(2)*MATDIMS(3)*MATDIMS(4)
    nw=neq

    IF(ALLOCATED(rwk)) DEALLOCATE(rwk)
    ALLOCATE(rwk(1,MATDIMS(1),MATDIMS(2),MATDIMS(3),MATDIMS(4),9))
    CALL dverk(ODESETID,neq,r8fcn3D,LocT,flds,tend,tol,ind,c24,nw,rwk)
    hmag=c24(14)
  ELSE IF(KSPACESOLVER) THEN
    fDIM1=MATDIMS(1)/2+1
    fDIM2=MATDIMS(2)
    fDIM3=MATDIMS(3)
    nfld =MATDIMS(4)

    neq=2*(MATDIMS(1)/2+1)*MATDIMS(2)*MATDIMS(3)*MATDIMS(4)
    nw=neq

    IF(ALLOCATED(rwk).EQV..FALSE.)   ALLOCATE(rwk(2,MATDIMS(1)/2+1,MATDIMS(2),MATDIMS(3),MATDIMS(4),9))
    IF(ALLOCATED(kflds).EQV..FALSE.) ALLOCATE(kflds(MATDIMS(1)/2+1,MATDIMS(2),MATDIMS(3),MATDIMS(4)))
    DO I=1,MATDIMS(4)
      CALL FFTW(flds(:,:,:,I),kflds(:,:,:,I))
    END DO
    CALL dverk(ODESETID,neq,c8fcn3D,LocT,kflds,tend,tol,ind,c24,nw,rwk)
    hmag=c24(14)
    DO I=1,MATDIMS(4)
      CALL iFFTW(kflds(:,:,:,I),flds(:,:,:,I))
    END DO
  END IF

 END SUBROUTINE r8r8DESOLVER3D

 SUBROUTINE r4c8DESOLVER3D(ODESETPIN,LX3D,LY3D,LZ3D,flds,t,tstp)
  IMPLICIT NONE
!********** DECLARATION OF INCOMING VARIABLES AND ARRAYS*********
  COMPLEX(8), INTENT(INOUT) :: flds(:,:,:,:)
  REAL(4),    INTENT(INOUT) :: tstp
  REAL(4),    INTENT(IN)    :: t
  REAL(4),    INTENT(IN)    :: LX3D(:,:,:),LY3D(:,:,:),LZ3D(:,:,:)
  INTEGER(4), INTENT(IN)    :: ODESETPIN

  CALL DESOLVER(ODESETPIN,1.0d0*LX3D,1.0d0*LY3D,1.0d0*LZ3D,flds,t,tstp)

  RETURN
 END SUBROUTINE r4c8DESOLVER3D

 SUBROUTINE r8c8DESOLVER3D(ODESETPIN,LX3D,LY3D,LZ3D,flds,t,tstp)
  implicit none
!********** DECLARATION OF INCOMING VARIABLES AND ARRAYS*********
  COMPLEX(8), INTENT(INOUT) :: flds(:,:,:,:)
  REAL(4),    INTENT(INOUT) :: tstp
  REAL(4),    INTENT(IN)    :: t
  REAL(8),    INTENT(IN)    :: LX3D(:,:,:),LY3D(:,:,:),LZ3D(:,:,:)
  INTEGER(4), INTENT(IN)    :: ODESETPIN
!**** DECLARATION OF LOCAL VARIABLES WITHOUT INITIAL VALUES *****
  REAL(8)                   :: LocT,LocTSTP
  INTEGER(4)                :: MATDIMS(4)
  INTEGER(4)                :: neq,ind,nw
  INTEGER(4)                :: I,J
  INTEGER(4), SAVE          :: ACCESSID = 1
  INTEGER(4), SAVE          :: ODESETID = 1
!********** DECLARATION OF VARIABLES WITH INITIAL VALUES ********
  REAL(8)                   :: c24(24)=0.0d0
  REAL(8)                   :: tol=1.0d-3
  REAL(8)                   :: hmax,hmag
  REAL(8)                   :: sci=0.0d0
  REAL(8)                   :: sig=0.03d0
  REAL(8)                   :: tend

!*************** DECLARATION OF ALLOCATABLE ARRAYS **************
  REAL(8),   ALLOCATABLE    :: rwk(:,:,:,:,:,:)

  MATDIMS=SHAPE(flds)
  fDIM1=MATDIMS(1)
  fDIM2=MATDIMS(2)
  fDIM3=MATDIMS(3)
  nfld =MATDIMS(4)

  if(sci.eq.0d0)then
    hmax=2d0
  elseif(sci.ne.0d0)then
    hmax=2d0/abs(sci)
  endif
  hmag=hmax*tol**(1d0/6d0)

  c24(4)=hmag
  c24(5)=sci
  c24(6)=hmax
  ind=2

  tend=LocT+LocTSTP

  neq=2*MATDIMS(1)*MATDIMS(2)*MATDIMS(3)*MATDIMS(4)
  nw=neq
  IF(ALLOCATED(rwk).EQV..FALSE.) ALLOCATE(rwk(2,MATDIMS(1),MATDIMS(2),MATDIMS(3),MATDIMS(4),9))
  CALL dverk(ODESETID,neq,c8fcn3D,LocT,flds,tend,tol,ind,c24,nw,rwk)
  hmag=c24(14)

  RETURN
 END SUBROUTINE r8c8DESOLVER3D


!****************************************************
!                      fcn
!****************************************************

 SUBROUTINE r8fcn1D(ODESETPIN,neq,t,fld,fldt)
  IMPLICIT NONE
  INTEGER(4),INTENT(IN)                       :: neq
  INTEGER(4),INTENT(IN)                       :: ODESETPIN
  INTEGER(4)                                  :: MATDIMS2D(2)
  INTEGER(4)                                  :: iDIM1,iDIM2

  REAL(8),INTENT(IN)                          :: t

  REAL(8),INTENT(INOUT),DIMENSION(fDIM1,nfld) :: fld
  REAL(8),INTENT(INOUT),DIMENSION(fDIM1,nfld) :: fldt

  REAL(8),ALLOCATABLE,DIMENSION(:,:,:)        ::  Lfld
  REAL(8),ALLOCATABLE,DIMENSION(:,:)          :: NLfld

  MATDIMS2D=SHAPE(fld)
  fldt=0.0d0
  IF(ALLOCATED( Lfld).EQV..FALSE.) ALLOCATE( Lfld(MATDIMS2D(1),MATDIMS2D(2),MATDIMS2D(2)))
  IF(ALLOCATED(NLfld).EQV..FALSE.) ALLOCATE(NLfld(MATDIMS2D(1),MATDIMS2D(2)))

  CALL ODESET(ODESETPIN,fld,Lfld,NLfld)

  DO iDIM1=1,MATDIMS2D(2)
   DO iDIM2=1,MATDIMS2D(2)
    fldt(:,iDIM1)=fldt(:,iDIM1)+ Lfld(:,iDIM1,iDIM2)
   END DO
    fldt(:,iDIM1)=fldt(:,iDIM1)+NLfld(:,iDIM1)
  END DO

  RETURN
 END SUBROUTINE r8fcn1D

 SUBROUTINE c8fcn1D(ODESETPIN,neq,t,fld,fldt)
  IMPLICIT NONE
  INTEGER(4),INTENT(IN)                          :: neq
  INTEGER(4),INTENT(IN)                          :: ODESETPIN
  INTEGER(4)                                     :: MATDIMS2D(2)
  INTEGER(4)                                     :: iDIM1,iDIM2

  REAL(8),INTENT(IN)                             :: t

  COMPLEX(8),INTENT(INOUT),DIMENSION(fDIM1,nfld) :: fld
  COMPLEX(8),INTENT(INOUT),DIMENSION(fDIM1,nfld) :: fldt

  COMPLEX(8),ALLOCATABLE,DIMENSION(:,:,:)        ::  Lfld
  COMPLEX(8),ALLOCATABLE,DIMENSION(:,:)          :: NLfld

  MATDIMS2D=SHAPE(fld)
  fldt=CMPLX(0.0d0,0.0d0)
  IF(ALLOCATED( Lfld).EQV..FALSE.) ALLOCATE( Lfld(MATDIMS2D(1),MATDIMS2D(2),MATDIMS2D(2)))
  IF(ALLOCATED(NLfld).EQV..FALSE.) ALLOCATE(NLfld(MATDIMS2D(1),MATDIMS2D(2)))

  CALL ODESET(ODESETPIN,fld,Lfld,NLfld)

  DO iDIM1=1,MATDIMS2D(2)
   DO iDIM2=1,MATDIMS2D(2)
    fldt(:,iDIM1)=fldt(:,iDIM1)+ Lfld(:,iDIM1,iDIM2)
   END DO
    fldt(:,iDIM1)=fldt(:,iDIM1)+NLfld(:,iDIM1)
  END DO

  RETURN
 END SUBROUTINE c8fcn1D

 SUBROUTINE r8fcn2D(ODESETPIN,neq,t,fld,fldt)
  IMPLICIT NONE
  INTEGER(4),INTENT(IN)                                :: neq
  INTEGER(4),INTENT(IN)                                :: ODESETPIN
  INTEGER(4)                                           :: MATDIMS3D(3)
  INTEGER(4)                                           :: iDIM1,iDIM2

  REAL(8),INTENT(IN)                                   :: t

  REAL(8),INTENT(INOUT),DIMENSION(fDIM1,fDIM2,nfld) :: fld
  REAL(8),INTENT(INOUT),DIMENSION(fDIM1,fDIM2,nfld) :: fldt

  REAL(8),ALLOCATABLE,DIMENSION(:,:,:,:)            ::  Lfld
  REAL(8),ALLOCATABLE,DIMENSION(:,:,:)              :: NLfld

  MATDIMS3D=SHAPE(fld)
  fldt=0.0d0
  IF(ALLOCATED( Lfld).EQV..FALSE.) ALLOCATE( Lfld(MATDIMS3D(1),MATDIMS3D(2),MATDIMS3D(3),MATDIMS3D(3)))
  IF(ALLOCATED(NLfld).EQV..FALSE.) ALLOCATE(NLfld(MATDIMS3D(1),MATDIMS3D(2),MATDIMS3D(3)))

  CALL ODESET(ODESETPIN,fld,Lfld,NLfld)

  DO iDIM1=1,MATDIMS3D(3)
   DO iDIM2=1,MATDIMS3D(3)
    fldt(:,:,iDIM1)=fldt(:,:,iDIM1)+ Lfld(:,:,iDIM1,iDIM2)
   END DO
    fldt(:,:,iDIM1)=fldt(:,:,iDIM1)+NLfld(:,:,iDIM1)
  END DO

  RETURN
 END SUBROUTINE r8fcn2D

 SUBROUTINE c8fcn2D(ODESETPIN,neq,t,fld,fldt)
  IMPLICIT NONE
  INTEGER(4),INTENT(IN)                                :: neq
  INTEGER(4),INTENT(IN)                                :: ODESETPIN
  INTEGER(4)                                           :: MATDIMS3D(3)
  INTEGER(4)                                           :: iDIM1,iDIM2

  REAL(8),INTENT(IN)                                   :: t

  COMPLEX(8),INTENT(INOUT),DIMENSION(fDIM1,fDIM2,nfld) :: fld
  COMPLEX(8),INTENT(INOUT),DIMENSION(fDIM1,fDIM2,nfld) :: fldt

  COMPLEX(8),ALLOCATABLE,DIMENSION(:,:,:,:)            ::  Lfld
  COMPLEX(8),ALLOCATABLE,DIMENSION(:,:,:)              :: NLfld

  MATDIMS3D=SHAPE(fld)
  fldt=CMPLX(0.0d0,0.0d0)
  IF(ALLOCATED( Lfld).EQV..FALSE.) ALLOCATE( Lfld(MATDIMS3D(1),MATDIMS3D(2),MATDIMS3D(3),MATDIMS3D(3)))
  IF(ALLOCATED(NLfld).EQV..FALSE.) ALLOCATE(NLfld(MATDIMS3D(1),MATDIMS3D(2),MATDIMS3D(3)))

  CALL ODESET(ODESETPIN,fld,Lfld,NLfld)

  DO iDIM1=1,MATDIMS3D(3)
   DO iDIM2=1,MATDIMS3D(3)
    fldt(:,:,iDIM1)=fldt(:,:,iDIM1)+ Lfld(:,:,iDIM1,iDIM2)
   END DO
    fldt(:,:,iDIM1)=fldt(:,:,iDIM1)+NLfld(:,:,iDIM1)
  END DO

  RETURN
 END SUBROUTINE c8fcn2D

 SUBROUTINE r8fcn3D(ODESETPIN,neq,t,fld,fldt)
  IMPLICIT NONE
  INTEGER(4),INTENT(IN)                                   :: neq
  INTEGER(4),INTENT(IN)                                   :: ODESETPIN
  INTEGER(4)                                              :: MATDIMS4D(4)
  INTEGER(4)                                              :: iDIM1,iDIM2

  REAL(8),INTENT(IN)                                      :: t

  REAL(8),INTENT(INOUT),DIMENSION(fDIM1,fDIM2,fDIM3,nfld) :: fld
  REAL(8),INTENT(INOUT),DIMENSION(fDIM1,fDIM2,fDIM3,nfld) :: fldt

  REAL(8),ALLOCATABLE,DIMENSION(:,:,:,:,:)                ::  Lfld
  REAL(8),ALLOCATABLE,DIMENSION(:,:,:,:)                  :: NLfld

  MATDIMS4D=SHAPE(fld)
  fldt=0.0d0
  IF(ALLOCATED( Lfld).EQV..FALSE.) ALLOCATE( Lfld(MATDIMS4D(1),MATDIMS4D(2),MATDIMS4D(3),MATDIMS4D(4),MATDIMS4D(4)))
  IF(ALLOCATED(NLfld).EQV..FALSE.) ALLOCATE(NLfld(MATDIMS4D(1),MATDIMS4D(2),MATDIMS4D(3),MATDIMS4D(4)))

  CALL ODESET(ODESETPIN,fld,Lfld,NLfld)

  DO iDIM1=1,MATDIMS4D(4)
   DO iDIM2=1,MATDIMS4D(4)
    fldt(:,:,:,iDIM1)=fldt(:,:,:,iDIM1)+ Lfld(:,:,:,iDIM1,iDIM2)
   END DO
    fldt(:,:,:,iDIM1)=fldt(:,:,:,iDIM1)+NLfld(:,:,:,iDIM1)
  END DO

  RETURN
 END SUBROUTINE r8fcn3D

 SUBROUTINE c8fcn3D(ODESETPIN,neq,t,fld,fldt)
  IMPLICIT NONE
  INTEGER(4),INTENT(IN)                                      :: neq
  INTEGER(4),INTENT(IN)                                      :: ODESETPIN
  INTEGER(4)                                                 :: MATDIMS4D(4)
  INTEGER(4)                                                 :: iDIM1,iDIM2

  REAL(8),INTENT(IN)                                         :: t

  COMPLEX(8),INTENT(INOUT),DIMENSION(fDIM1,fDIM2,fDIM3,nfld) :: fld
  COMPLEX(8),INTENT(INOUT),DIMENSION(fDIM1,fDIM2,fDIM3,nfld) :: fldt

  COMPLEX(8),ALLOCATABLE,DIMENSION(:,:,:,:,:)                ::  Lfld
  COMPLEX(8),ALLOCATABLE,DIMENSION(:,:,:,:)                  :: NLfld

  MATDIMS4D=SHAPE(fld)
  fldt=CMPLX(0.0d0,0.0d0)
  IF(ALLOCATED( Lfld).EQV..FALSE.) ALLOCATE( Lfld(MATDIMS4D(1),MATDIMS4D(2),MATDIMS4D(3),MATDIMS4D(4),MATDIMS4D(4)))
  IF(ALLOCATED(NLfld).EQV..FALSE.) ALLOCATE(NLfld(MATDIMS4D(1),MATDIMS4D(2),MATDIMS4D(3),MATDIMS4D(4)))

  CALL ODESET(ODESETPIN,fld,Lfld,NLfld)

  DO iDIM1=1,MATDIMS4D(4)
   DO iDIM2=1,MATDIMS4D(4)
    fldt(:,:,:,iDIM1)=fldt(:,:,:,iDIM1)+ Lfld(:,:,:,iDIM1,iDIM2)
   END DO
    fldt(:,:,:,iDIM1)=fldt(:,:,:,iDIM1)+NLfld(:,:,:,iDIM1)
  END DO

  RETURN
 END SUBROUTINE c8fcn3D



!********************************************************************************************
! TO FIND THE EIGENVALUES, EIGENVECTORS, GROWTH RATE AND FREQUESNCIES OF THE PDE SYSTEM
!********************************************************************************************
 SUBROUTINE r4SYSEIGENVALS(ODESETID,nDE,LOCALTYPE,LX3D,LY3D,LZ3D)
  IMPLICIT NONE
  REAL(4),     INTENT(IN)          :: LX3D(:,:,:),LY3D(:,:,:),LZ3D(:,:,:)
  INTEGER(4),  INTENT(IN)          :: nDE,ODESETID
  CHARACTER(*),INTENT(IN)          :: LOCALTYPE(:)

  CALL SYSEIGENVV(ODESETID,nDE,LOCALTYPE,LX3D,LY3D,LZ3D)
  IF(ALLOCATED(EIGENVEC7D)) DEALLOCATE(EIGENVEC7D)

 END SUBROUTINE r4SYSEIGENVALS

 SUBROUTINE r8SYSEIGENVALS(ODESETID,nDE,LOCALTYPE,LX3D,LY3D,LZ3D)
  IMPLICIT NONE
  REAL(8),     INTENT(IN)          :: LX3D(:,:,:),LY3D(:,:,:),LZ3D(:,:,:)
  INTEGER(4),  INTENT(IN)          :: nDE,ODESETID
  CHARACTER(*),INTENT(IN)          :: LOCALTYPE(:)

  CALL SYSEIGENVV(ODESETID,nDE,LOCALTYPE,LX3D,LY3D,LZ3D)
  IF(ALLOCATED(EIGENVEC7D)) DEALLOCATE(EIGENVEC7D)

 END SUBROUTINE r8SYSEIGENVALS

 SUBROUTINE r4SYSEIGENVECS(ODESETID,nDE,LOCALTYPE,LX3D,LY3D,LZ3D)
  IMPLICIT NONE
  REAL(4),     INTENT(IN)          :: LX3D(:,:,:),LY3D(:,:,:),LZ3D(:,:,:)
  INTEGER(4),  INTENT(IN)          :: nDE,ODESETID
  CHARACTER(*),INTENT(IN)          :: LOCALTYPE(:)

  CALL SYSEIGENVV(ODESETID,nDE,LOCALTYPE,LX3D,LY3D,LZ3D)
  IF(ALLOCATED(EIGENVAL7D)) DEALLOCATE(EIGENVAL7D)

 END SUBROUTINE r4SYSEIGENVECS

 SUBROUTINE r8SYSEIGENVECS(ODESETID,nDE,LOCALTYPE,LX3D,LY3D,LZ3D)
  IMPLICIT NONE
  REAL(8),     INTENT(IN)          :: LX3D(:,:,:),LY3D(:,:,:),LZ3D(:,:,:)
  INTEGER(4),  INTENT(IN)          :: nDE,ODESETID
  CHARACTER(*),INTENT(IN)          :: LOCALTYPE(:)

  CALL SYSEIGENVV(ODESETID,nDE,LOCALTYPE,LX3D,LY3D,LZ3D)
  IF(ALLOCATED(EIGENVAL7D)) DEALLOCATE(EIGENVAL7D)

 END SUBROUTINE r8SYSEIGENVECS

 SUBROUTINE r4SYSEIGENVV(ODESETID,nDE,LOCALTYPE,LX3D,LY3D,LZ3D)
  IMPLICIT NONE
  REAL(4),     INTENT(IN)          :: LX3D(:,:,:),LY3D(:,:,:),LZ3D(:,:,:)
  INTEGER(4),  INTENT(IN)          :: nDE,ODESETID
  CHARACTER(*),INTENT(IN)          :: LOCALTYPE(:)

  CALL SYSEIGENVV(ODESETID,nDE,LOCALTYPE,1.0d0*LX3D,1.0d0*LY3D,1.0d0*LZ3D)

  RETURN
 END SUBROUTINE r4SYSEIGENVV

 SUBROUTINE r8SYSEIGENVV(ODESETID,nDE,LOCALTYPE,LX3D,LY3D,LZ3D)
  IMPLICIT NONE
  REAL(8),     INTENT(IN)        :: LX3D(:,:,:),LY3D(:,:,:),LZ3D(:,:,:)
  INTEGER(4),  INTENT(IN)        :: nDE,ODESETID
  CHARACTER(*),INTENT(IN)        :: LOCALTYPE(:)

  CHARACTER(10)                  :: EIGENLOCALTYPE(3)
  REAL(8),   ALLOCATABLE         :: rfld1D(:,:),    rLfld1D(:,:,:),    rNLfld1D(:,:)
  REAL(8),   ALLOCATABLE         :: rfld2D(:,:,:),  rLfld2D(:,:,:,:),  rNLfld2D(:,:,:)
  REAL(8),   ALLOCATABLE         :: rfld3D(:,:,:,:),rLfld3D(:,:,:,:,:),rNLfld3D(:,:,:,:)
  COMPLEX(8),ALLOCATABLE         :: kfld1D(:,:),    kLfld1D(:,:,:),    kNLfld1D(:,:)
  COMPLEX(8),ALLOCATABLE         :: kfld2D(:,:,:),  kLfld2D(:,:,:,:),  kNLfld2D(:,:,:)
  COMPLEX(8),ALLOCATABLE         :: kfld3D(:,:,:,:),kLfld3D(:,:,:,:,:),kNLfld3D(:,:,:,:)
  COMPLEX(4),ALLOCATABLE         :: LEIGENVAL4D(:,:,:,:),LEIGENVEC5D(:,:,:,:,:)
  COMPLEX(4),ALLOCATABLE         :: EIGENMAT(:,:)
  INTEGER(4)                     :: LXDIM,LYDIM,LZDIM
  INTEGER(4)                     :: KxDIM,KyDIM,KzDIM
  INTEGER(4)                     :: iKx,iKy,iKz
  INTEGER(4)                     :: iXDIM,iYDIM,iZDIM
  INTEGER(4)                     :: iDE,iDE1,iDE2
  INTEGER(4)                     :: EIGENID,LMAXGINDEX(4)
  INTEGER(4)                     :: MATDIMS3D(3)

  PRINT*, "EIGENS CALCULATION >>> STARTED!"

  MATDIMS3D=SHAPE(LX3D)

  nfld=nDE
 !ALLOCATE THE BOUNDARY CONDITIONS ALLOCATABLE ARRAYS
  CALL RESETBC(nDE)
 !ALLOCATE THE DIFFERENTIATION SCHEMES ARRAYS
  CALL RESETDSCHEMES(nDE)
 !ALLOCATE THE INTEGRATION SCHEMES ARRAYS
  CALL RESETISCHEMES(nDE)
 !RESET AND SET THE FFTW VARIABLES
  IF((MATDIMS3D(2).EQ.1).AND.(MATDIMS3D(3).EQ.1)) THEN
   IF(KSPACESOLVER) THEN
    fDIM1=MATDIMS3D(1)/2+1
   ELSE IF(RSPACESOLVER) THEN
    fDIM1=MATDIMS3D(1)
   END IF
  ELSE IF((MATDIMS3D(1).EQ.1).AND.(MATDIMS3D(3).EQ.1)) THEN
   IF(KSPACESOLVER) THEN
    fDIM1=MATDIMS3D(2)/2+1
   ELSE IF(RSPACESOLVER) THEN
    fDIM1=MATDIMS3D(2)
   END IF
  ELSE IF((MATDIMS3D(1).EQ.1).AND.(MATDIMS3D(2).EQ.1)) THEN
   IF(KSPACESOLVER) THEN
    fDIM1=MATDIMS3D(3)/2+1
   ELSE IF(KSPACESOLVER) THEN
    fDIM1=MATDIMS3D(3)
   END IF
  ELSE IF(MATDIMS3D(1).EQ.1) THEN
   IF(KSPACESOLVER) THEN
    fDIM1=MATDIMS3D(2)/2+1
    fDIM2=MATDIMS3D(3)
   ELSE IF(RSPACESOLVER) THEN
    fDIM1=MATDIMS3D(2)
    fDIM2=MATDIMS3D(3)
   END IF
  ELSE IF(MATDIMS3D(2).EQ.1) THEN
   IF(KSPACESOLVER) THEN
    fDIM1=MATDIMS3D(1)/2+1
    fDIM2=MATDIMS3D(3)
   ELSE IF(RSPACESOLVER) THEN
    fDIM1=MATDIMS3D(1)
    fDIM2=MATDIMS3D(3)
   END IF
  ELSE IF(MATDIMS3D(3).EQ.1) THEN
   IF(KSPACESOLVER) THEN
    fDIM1=MATDIMS3D(1)/2+1
    fDIM2=MATDIMS3D(2)
   ELSE IF(RSPACESOLVER) THEN
    fDIM1=MATDIMS3D(1)
    fDIM2=MATDIMS3D(2)
   END IF
  ELSE
   IF(KSPACESOLVER) THEN
    fDIM1=MATDIMS3D(1)/2+1
    fDIM2=MATDIMS3D(2)
    fDIM3=MATDIMS3D(3)
   ELSE IF(RSPACESOLVER) THEN
    fDIM1=MATDIMS3D(1)
    fDIM2=MATDIMS3D(2)
    fDIM3=MATDIMS3D(3)
   END IF
  END IF

 !IF((MATDIMS3D(2).EQ.1).AND.(MATDIMS3D(2).EQ.1)) THEN
 ! IF(ALLOCATED(iKx1D)) DEALLOCATE(iKx1D)
 ! ALLOCATE(iKx1D(MATDIMS3D(1)/2+1))
 ! CALL KGRIDBUILDER(LX3D(:,1,1),iKx1D)
 !ELSE IF((MATDIMS3D(1).EQ.1).AND.(MATDIMS3D(3).EQ.1)) THEN
 ! IF(ALLOCATED(iKx1D)) DEALLOCATE(iKx1D)
 ! ALLOCATE(iKx1D(MATDIMS3D(2)/2+1))
 ! CALL KGRIDBUILDER(LY3D(1,:,1),iKx1D)
 !ELSE IF((MATDIMS3D(1).EQ.1).AND.(MATDIMS3D(2).EQ.1)) THEN
 ! IF(ALLOCATED(iKx1D)) DEALLOCATE(iKx1D)
 ! ALLOCATE(iKx1D(MATDIMS3D(3)/2+1))
 ! CALL KGRIDBUILDER(LZ3D(1,1,:),iKx1D)
 !ELSE IF(MATDIMS3D(3).EQ.1) THEN
 ! IF(ALLOCATED(iKx2D)) DEALLOCATE(iKx2D)
 ! IF(ALLOCATED(iKy2D)) DEALLOCATE(iKy2D)
 ! ALLOCATE(iKx2D(MATDIMS3D(1)/2+1,MATDIMS3D(2)))
 ! ALLOCATE(iKy2D(MATDIMS3D(1)/2+1,MATDIMS3D(2)))
 ! CALL KGRIDBUILDER(LX3D(:,:,1),iKx2D)
 ! CALL KGRIDBUILDER(LY3D(:,:,1),iKy2D)
 !ELSE IF(MATDIMS3D(2).EQ.1) THEN
 ! IF(ALLOCATED(iKx2D)) DEALLOCATE(iKx2D)
 ! IF(ALLOCATED(iKy2D)) DEALLOCATE(iKy2D)
 ! ALLOCATE(iKx2D(MATDIMS3D(1)/2+1,MATDIMS3D(3)))
 ! ALLOCATE(iKy2D(MATDIMS3D(1)/2+1,MATDIMS3D(3)))
 ! CALL KGRIDBUILDER(LX3D(:,1,:),iKx2D)
 ! CALL KGRIDBUILDER(LZ3D(:,1,:),iKy2D)
 !ELSE IF(MATDIMS3D(1).EQ.1) THEN
 ! IF(ALLOCATED(iKx2D)) DEALLOCATE(iKx2D)
 ! IF(ALLOCATED(iKy2D)) DEALLOCATE(iKy2D)
 ! ALLOCATE(iKx2D(MATDIMS3D(2)/2+1,MATDIMS3D(3)))
 ! ALLOCATE(iKy2D(MATDIMS3D(2)/2+1,MATDIMS3D(3)))
 ! CALL KGRIDBUILDER(LY3D(1,:,:),iKx2D)
 ! CALL KGRIDBUILDER(LZ3D(1,:,:),iKy2D)
 !ELSE
 ! IF(ALLOCATED(iKx3D)) DEALLOCATE(iKx3D)
 ! IF(ALLOCATED(iKy3D)) DEALLOCATE(iKy3D)
 ! IF(ALLOCATED(iKz3D)) DEALLOCATE(iKz3D)
 ! ALLOCATE(iKx3D(MATDIMS3D(1)/2+1,MATDIMS3D(2),MATDIMS3D(3)))
 ! ALLOCATE(iKy3D(MATDIMS3D(1)/2+1,MATDIMS3D(2),MATDIMS3D(3)))
 ! ALLOCATE(iKz3D(MATDIMS3D(1)/2+1,MATDIMS3D(2),MATDIMS3D(3)))
 ! CALL KGRIDBUILDER(LX3D,iKx3D)
 ! CALL KGRIDBUILDER(LY3D,iKy3D)
 ! CALL KGRIDBUILDER(LZ3D,iKz3D)
 !END IF

 !IF((MATDIMS3D(2).EQ.1).AND.(MATDIMS3D(2).EQ.1)) THEN
 ! IF(ALLOCATED(iKx1D)) DEALLOCATE(iKx1D)
 ! ALLOCATE(iKx1D(MATDIMS3D(1)/2+1))
 ! CALL KGRIDBUILDER(LX3D(:,1,1),iKx1D)
 !ELSE IF((MATDIMS3D(1).EQ.1).AND.(MATDIMS3D(3).EQ.1)) THEN
 ! IF(ALLOCATED(iKy1D)) DEALLOCATE(iKy1D)
 ! ALLOCATE(iKy1D(MATDIMS3D(2)/2+1))
 ! CALL KGRIDBUILDER(LY3D(1,:,1),iKy1D)
 !ELSE IF((MATDIMS3D(1).EQ.1).AND.(MATDIMS3D(2).EQ.1)) THEN
 ! IF(ALLOCATED(iKz1D)) DEALLOCATE(iKz1D)
 ! ALLOCATE(iKz1D(MATDIMS3D(3)/2+1))
 ! CALL KGRIDBUILDER(LZ3D(1,1,:),iKz1D)
 !ELSE IF(MATDIMS3D(3).EQ.1) THEN
 ! IF(ALLOCATED(iKx2D)) DEALLOCATE(iKx2D)
 ! IF(ALLOCATED(iKy2D)) DEALLOCATE(iKy2D)
 ! ALLOCATE(iKx2D(MATDIMS3D(1)/2+1,MATDIMS3D(2)))
 ! ALLOCATE(iKy2D(MATDIMS3D(1)/2+1,MATDIMS3D(2)))
 ! CALL KGRIDBUILDER(LX3D(:,:,1),iKx2D)
 ! CALL KGRIDBUILDER(LY3D(:,:,1),iKy2D)
 !ELSE IF(MATDIMS3D(2).EQ.1) THEN
 ! IF(ALLOCATED(iKx2D)) DEALLOCATE(iKx2D)
 ! IF(ALLOCATED(iKz2D)) DEALLOCATE(iKz2D)
 ! ALLOCATE(iKx2D(MATDIMS3D(1)/2+1,MATDIMS3D(3)))
 ! ALLOCATE(iKz2D(MATDIMS3D(1)/2+1,MATDIMS3D(3)))
 ! CALL KGRIDBUILDER(LX3D(:,1,:),iKx2D)
 ! CALL KGRIDBUILDER(LZ3D(:,1,:),iKz2D)
 !ELSE IF(MATDIMS3D(1).EQ.1) THEN
 ! IF(ALLOCATED(iKy2D)) DEALLOCATE(iKy2D)
 ! IF(ALLOCATED(iKz2D)) DEALLOCATE(iKz2D)
 ! ALLOCATE(iKy2D(MATDIMS3D(2)/2+1,MATDIMS3D(3)))
 ! ALLOCATE(iKz2D(MATDIMS3D(2)/2+1,MATDIMS3D(3)))
 ! CALL KGRIDBUILDER(LY3D(1,:,:),iKy2D)
 ! CALL KGRIDBUILDER(LZ3D(1,:,:),iKz2D)
 !ELSE
 ! IF(ALLOCATED(iKx3D)) DEALLOCATE(iKx3D)
 ! IF(ALLOCATED(iKy3D)) DEALLOCATE(iKy3D)
 ! IF(ALLOCATED(iKz3D)) DEALLOCATE(iKz3D)
 ! ALLOCATE(iKx3D(MATDIMS3D(1)/2+1,MATDIMS3D(2),MATDIMS3D(3)))
 ! ALLOCATE(iKy3D(MATDIMS3D(1)/2+1,MATDIMS3D(2),MATDIMS3D(3)))
 ! ALLOCATE(iKz3D(MATDIMS3D(1)/2+1,MATDIMS3D(2),MATDIMS3D(3)))
 ! CALL KGRIDBUILDER(LX3D,iKx3D)
 ! CALL KGRIDBUILDER(LY3D,iKy3D)
 ! CALL KGRIDBUILDER(LZ3D,iKz3D)
 !END IF


  IF(SIZE(LOCALTYPE).EQ.1) THEN
   EIGENLOCALTYPE(1)=TRIM(LOCALTYPE(1))
   EIGENLOCALTYPE(2)=TRIM(LOCALTYPE(1))
   EIGENLOCALTYPE(3)=TRIM(LOCALTYPE(1))
  ELSE IF(SIZE(LOCALTYPE).EQ.2) THEN
   EIGENLOCALTYPE(1)=TRIM(LOCALTYPE(1))
   EIGENLOCALTYPE(2)=TRIM(LOCALTYPE(2))
   EIGENLOCALTYPE(3)=TRIM(LOCALTYPE(2))
  ELSE IF(SIZE(LOCALTYPE).EQ.3) THEN
   EIGENLOCALTYPE(1)=TRIM(LOCALTYPE(1))
   EIGENLOCALTYPE(2)=TRIM(LOCALTYPE(2))
   EIGENLOCALTYPE(3)=TRIM(LOCALTYPE(3))
  ELSE IF(SIZE(LOCALTYPE).GT.3) THEN
   PRINT*, "WARNING: LOCALTYPE DIMENSIONS IS GREATER THAN THE SYSTEM DIMESNIONS"
   EIGENLOCALTYPE(1)=TRIM(LOCALTYPE(1))
   EIGENLOCALTYPE(2)=TRIM(LOCALTYPE(2))
   EIGENLOCALTYPE(3)=TRIM(LOCALTYPE(3))
  END IF

  IF(to_LowerCase(TRIM(EIGENLOCALTYPE(1))).EQ."") THEN
   LXDIM=1
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(1))).EQ."local") THEN
   LXDIM=1
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(1))).EQ."nonlocal") THEN
   LXDIM=MATDIMS3D(1)
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(1))).EQ."quasilocal") THEN
   LXDIM=MATDIMS3D(1)
  END IF

  IF(to_LowerCase(TRIM(EIGENLOCALTYPE(2))).EQ."") THEN
   LYDIM=1
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(2))).EQ."local") THEN
   LYDIM=1
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(2))).EQ."nonlocal") THEN
   LYDIM=MATDIMS3D(2)
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(2))).EQ."quasilocal") THEN
   LYDIM=MATDIMS3D(2)
  END IF

  IF(to_LowerCase(TRIM(EIGENLOCALTYPE(3))).EQ."") THEN
   LZDIM=1
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(3))).EQ."local") THEN
   LZDIM=1
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(3))).EQ."nonlocal") THEN
   LZDIM=MATDIMS3D(3)
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(3))).EQ."quasilocal") THEN
   LZDIM=MATDIMS3D(3)
  END IF

  IF(ALLOCATED(eiKx1D).EQV..FALSE.) THEN
   ALLOCATE(eiKx1D(MATDIMS3D(1)/2+1))
   CALL KGRIDBUILDER(LX3D(:,1,1),eiKx1D)
  END IF
  IF(ALLOCATED(eiKy1D).EQV..FALSE.) THEN
   ALLOCATE(eiKy1D(MATDIMS3D(2)/2+1))
   CALL KGRIDBUILDER(LY3D(1,:,1),eiKy1D)
  END IF
  IF(ALLOCATED(eiKz1D).EQV..FALSE.) THEN
   ALLOCATE(eiKz1D(MATDIMS3D(3)/2+1))
   CALL KGRIDBUILDER(LZ3D(1,1,:),eiKz1D)
  END IF
  KxDIM=SIZE(eiKx1D)
  KyDIM=SIZE(eiKy1D)
  KzDIM=SIZE(eiKz1D)
  IF((KxDIM.EQ.1).AND.(KyDIM.EQ.1)) THEN
   IF(ALLOCATED(iKz1D).EQV..TRUE.) THEN
    PRINT*,"WARNING: K-SPACE IS RE-ALLOCATED"
    DEALLOCATE(iKz1D)
   END IF
   ALLOCATE(iKz1D(KzDIM))
   iKz1D=(0.0d0,1.0d0)*eiKz1D
  ELSE IF((KyDIM.EQ.1).AND.(KzDIM.EQ.1)) THEN
   IF(ALLOCATED(iKx1D).EQV..TRUE.) THEN
    PRINT*,"WARNING: K-SPACE IS RE-ALLOCATED"
    DEALLOCATE(iKx1D)
   END IF
   ALLOCATE(iKx1D(KxDIM))
   iKx1D=(0.0d0,1.0d0)*eiKx1D
  ELSE IF((KxDIM.EQ.1).AND.(KzDIM.EQ.1)) THEN
   IF(ALLOCATED(iKy1D).EQV..TRUE.) THEN
    PRINT*,"WARNING: K-SPACE IS RE-ALLOCATED"
    DEALLOCATE(iKy1D)
   END IF
   ALLOCATE(iKy1D(KyDIM))
   iKy1D=(0.0d0,1.0d0)*eiKy1D
  ELSE IF(KxDIM.EQ.1) THEN
   IF((ALLOCATED(iKy2D).EQV..TRUE.).OR.(ALLOCATED(iKz2D).EQV..TRUE.)) THEN
    PRINT*,"WARNING: K-SPACE IS RE-ALLOCATED"
    DEALLOCATE(iKy2D)
    DEALLOCATE(iKz2D)
   END IF
   ALLOCATE(iKy2D(KyDIM,2*(KzDIM-1)))
   ALLOCATE(iKz2D(KyDIM,2*(KzDIM-1)))
   DO iKy=1,KyDIM
    iKz2D(iKy,1:KzDIM)=(0.0d0,1.0d0)*eiKz1D
    iKz2D(iKy,KzDIM+1:2*(KzDIM-1))=(0.0d0,-1.0d0)*ROT180(eiKz1D(1:KzDIM-1))
   END DO
   DO iKz=1,2*(KzDIM-1)
    iKy2D(:,iKz)=(0.0d0,1.0d0)*eiKy1D
   END DO
  ELSE IF(KyDIM.EQ.1) THEN
   IF((ALLOCATED(iKx2D).EQV..TRUE.).OR.(ALLOCATED(iKz2D).EQV..TRUE.)) THEN
    PRINT*,"WARNING: K-SPACE IS RE-ALLOCATED"
    DEALLOCATE(iKx2D)
    DEALLOCATE(iKz2D)
   END IF
   ALLOCATE(iKx2D(KxDIM,2*(KzDIM-1)))
   ALLOCATE(iKz2D(KxDIM,2*(KzDIM-1)))
   DO iKx=1,KxDIM
    iKz2D(iKx,1:KzDIM)=(0.0d0,1.0d0)*eiKz1D
    iKz2D(iKx,KzDIM+1:2*(KzDIM-1))=(0.0d0,-1.0d0)*ROT180(eiKz1D(1:KzDIM-1))
   END DO
   DO iKz=1,2*(KzDIM-1)
    iKx2D(:,iKz)=(0.0d0,1.0d0)*eiKx1D
   END DO
  ELSE IF(KzDIM.EQ.1) THEN
   IF((ALLOCATED(iKx2D).EQV..TRUE.).OR.(ALLOCATED(iKy2D).EQV..TRUE.)) THEN
    PRINT*,"WARNING: K-SPACE IS RE-ALLOCATED"
    DEALLOCATE(iKx2D); DEALLOCATE(iKy2D)
   END IF
   ALLOCATE(iKx2D(KxDIM,2*(KyDIM-1)))
   ALLOCATE(iKy2D(KxDIM,2*(KyDIM-1)))
   DO iKx=1,KxDIM
    iKy2D(iKx,1:KyDIM)=(0.0d0,1.0d0)*eiKy1D
    iKy2D(iKx,KyDIM+1:2*(KyDIM-1))=(0.0d0,-1.0d0)*ROT180(eiKy1D(1:KyDIM-1))
   END DO
   DO iKy=1,2*(KyDIM-1)
    iKx2D(:,iKy)=(0.0d0,1.0d0)*eiKx1D
   END DO
  ELSE
   IF((ALLOCATED(iKx3D).EQV..TRUE.).OR.(ALLOCATED(iKy3D).EQV..TRUE.).OR.(ALLOCATED(iKz3D).EQV..TRUE.)) THEN
    PRINT*,"WARNING: K-SPACE IS RE-ALLOCATED"
    DEALLOCATE(iKx3D); DEALLOCATE(iKy3D); DEALLOCATE(iKz3D)
   END IF
   ALLOCATE(iKx3D(KxDIM,2*(KyDIM-1),2*(KzDIM-1)))
   ALLOCATE(iKy3D(KxDIM,2*(KyDIM-1),2*(KzDIM-1)))
   ALLOCATE(iKz3D(KxDIM,2*(KyDIM-1),2*(KzDIM-1)))
   DO iKx=1,KxDIM
    DO iKz=1,2*(KzDIM-1)
     iKy3D(iKx,1:KyDIM,iKz)=(0.0d0,1.0d0)*eiKy1D
     iKy3D(iKx,KyDIM+1:2*(KyDIM-1),iKz)=(0.0d0,-1.0d0)*ROT180(eiKy1D(1:KyDIM-1))
    END DO
   END DO
   DO iKy=1,2*(KyDIM-1)
    DO iKz=1,2*(KzDIM-1)
     iKx3D(:,iKy,iKz)=(0.0d0,1.0d0)*eiKx1D
    END DO
   END DO
   DO iKx=1,KxDIM
    DO iKy=1,2*(KyDIM-1)
     iKz3D(iKx,iKy,1:KzDIM)=(0.0d0,1.0d0)*eiKz1D
     iKz3D(iKy,iKy,KzDIM+1:2*(KzDIM-1))=(0.0d0,-1.0d0)*ROT180(eiKz1D(1:KzDIM-1))
    END DO
   END DO
  END IF

 !IF(ALLOCATED(eiKx1D).EQV..FALSE.) THEN
 ! ALLOCATE(eiKx1D(KxDIM))
 ! DO iKx=1,KxDIM
 !  eiKx1D(iKx)=srtKx+(iKx-1)*((endKx-srtKx)/KxDIM)
 ! END DO
 !ELSE IF(ALLOCATED(eiKx1D)) THEN
 ! KxDIM=SIZE(eiKx1D)
 !END IF
 !IF(ALLOCATED(eiKy1D).EQV..FALSE.) THEN
 ! ALLOCATE(eiKy1D(KyDIM))
 ! DO iKy=1,KyDIM
 !  eiKy1D(iKy)=srtKy+(iKy-1)*((endKy-srtKy)/KyDIM)
 ! END DO
 !ELSE IF(ALLOCATED(eiKy1D)) THEN
 ! KyDIM=SIZE(eiKy1D)
 !END IF
 !IF(ALLOCATED(eiKz1D).EQV..FALSE.) THEN
 ! ALLOCATE(eiKz1D(KzDIM))
 ! DO iKz=1,KzDIM
 !  eiKz1D(iKz)=srtKz+(iKz-1)*((endKz-srtKz)/KzDIM)
 ! END DO
 !ELSE IF(ALLOCATED(eiKz1D)) THEN
 ! KzDIM=SIZE(eiKz1D)
 !END IF

  IF(ALLOCATED(EIGENMAT)) DEALLOCATE(EIGENMAT)
  ALLOCATE(EIGENMAT(nDE,nDE))

  IF(ALLOCATED(LEIGENVAL4D)) DEALLOCATE(LEIGENVAL4D)
  IF(ALLOCATED(LEIGENVEC5D)) DEALLOCATE(LEIGENVEC5D)
  ALLOCATE(LEIGENVAL4D(KxDIM,KyDIM,KzDIM,nDE))
  ALLOCATE(LEIGENVEC5D(KxDIM,KyDIM,KzDIM,nDE,nDE))

  IF(ALLOCATED(EIGENVAL7D)) DEALLOCATE(EIGENVAL7D)
  IF(ALLOCATED(EIGENVEC7D)) DEALLOCATE(EIGENVEC7D)
  ALLOCATE(EIGENVAL7D(LXDIM,LYDIM,LZDIM,KxDIM,KyDIM,KzDIM,nDE))
  ALLOCATE(EIGENVEC7D(LXDIM,LYDIM,LZDIM,KxDIM,KyDIM,KzDIM,nDE))

  IF(ALLOCATED(KxGMAX3D)) DEALLOCATE(KxGMAX3D)
  IF(ALLOCATED(KyGMAX3D)) DEALLOCATE(KyGMAX3D)
  IF(ALLOCATED(KzGMAX3D)) DEALLOCATE(KzGMAX3D)
  ALLOCATE(KxGMAX3D(LXDIM,LYDIM,LZDIM))
  ALLOCATE(KyGMAX3D(LXDIM,LYDIM,LZDIM))
  ALLOCATE(KzGMAX3D(LXDIM,LYDIM,LZDIM))

  IF(KSPACESOLVER.OR.RSPACESOLVER) THEN
   IF((KxDIM.EQ.1).AND.(KyDIM.EQ.1)) THEN
    IF(ALLOCATED(  kfld1D).EQV..FALSE.) ALLOCATE(  kfld1D(KzDIM,nDE))
    IF(ALLOCATED( kLfld1D).EQV..FALSE.) ALLOCATE( kLfld1D(KzDIM,nDE,nDE))
    IF(ALLOCATED(kNLfld1D).EQV..FALSE.) ALLOCATE(kNLfld1D(KzDIM,nDE))
      kfld1D=CMPLX(1.0d0,0.0d0)
     kLfld1D=CMPLX(0.0d0,0.0d0)
    kNLfld1D=CMPLX(0.0d0,0.0d0)
   ELSE IF((KyDIM.EQ.1).AND.(KzDIM.EQ.1)) THEN
    IF(ALLOCATED(  kfld1D).EQV..FALSE.) ALLOCATE(  kfld1D(KxDIM,nDE))
    IF(ALLOCATED( kLfld1D).EQV..FALSE.) ALLOCATE( kLfld1D(KxDIM,nDE,nDE))
    IF(ALLOCATED(kNLfld1D).EQV..FALSE.) ALLOCATE(kNLfld1D(KxDIM,nDE))
      kfld1D=CMPLX(1.0d0,0.0d0)
     kLfld1D=CMPLX(0.0d0,0.0d0)
    kNLfld1D=CMPLX(0.0d0,0.0d0)
   ELSE IF((KxDIM.EQ.1).AND.(KzDIM.EQ.1)) THEN
    IF(ALLOCATED(  kfld1D).EQV..FALSE.) ALLOCATE(  kfld1D(KyDIM,nDE))
    IF(ALLOCATED( kLfld1D).EQV..FALSE.) ALLOCATE( kLfld1D(KyDIM,nDE,nDE))
    IF(ALLOCATED(kNLfld1D).EQV..FALSE.) ALLOCATE(kNLfld1D(KyDIM,nDE))
      kfld1D=CMPLX(1.0d0,0.0d0)
     kLfld1D=CMPLX(0.0d0,0.0d0)
    kNLfld1D=CMPLX(0.0d0,0.0d0)
   ELSE IF(KxDIM.EQ.1) THEN
    IF(ALLOCATED(  kfld2D).EQV..FALSE.) ALLOCATE(  kfld2D(KyDIM,2*(KzDIM-1),nDE))
    IF(ALLOCATED( kLfld2D).EQV..FALSE.) ALLOCATE( kLfld2D(KyDIM,2*(KzDIM-1),nDE,nDE))
    IF(ALLOCATED(kNLfld2D).EQV..FALSE.) ALLOCATE(kNLfld2D(KyDIM,2*(KzDIM-1),nDE))
      kfld2D=CMPLX(1.0d0,0.0d0)
     kLfld2D=CMPLX(0.0d0,0.0d0)
    kNLfld2D=CMPLX(0.0d0,0.0d0)
   ELSE IF(KyDIM.EQ.1) THEN
    IF(ALLOCATED(  kfld2D).EQV..FALSE.) ALLOCATE(  kfld2D(KxDIM,2*(KzDIM-1),nDE))
    IF(ALLOCATED( kLfld2D).EQV..FALSE.) ALLOCATE( kLfld2D(KxDIM,2*(KzDIM-1),nDE,nDE))
    IF(ALLOCATED(kNLfld2D).EQV..FALSE.) ALLOCATE(kNLfld2D(KxDIM,2*(KzDIM-1),nDE))
      kfld2D=CMPLX(1.0d0,0.0d0)
     kLfld2D=CMPLX(0.0d0,0.0d0)
    kNLfld2D=CMPLX(0.0d0,0.0d0)
   ELSE IF(KzDIM.EQ.1) THEN
    IF(ALLOCATED(  kfld2D).EQV..FALSE.) ALLOCATE(  kfld2D(KxDIM,2*(KyDIM-1),nDE))
    IF(ALLOCATED( kLfld2D).EQV..FALSE.) ALLOCATE( kLfld2D(KxDIM,2*(KyDIM-1),nDE,nDE))
    IF(ALLOCATED(kNLfld2D).EQV..FALSE.) ALLOCATE(kNLfld2D(KxDIM,2*(KyDIM-1),nDE))
      kfld2D=CMPLX(1.0d0,0.0d0)
     kLfld2D=CMPLX(0.0d0,0.0d0)
    kNLfld2D=CMPLX(0.0d0,0.0d0)
   ELSE
    IF(ALLOCATED(  kfld3D).EQV..FALSE.) ALLOCATE(  kfld3D(KxDIM,2*(KyDIM-1),2*(KzDIM-1),nDE))
    IF(ALLOCATED( kLfld3D).EQV..FALSE.) ALLOCATE( kLfld3D(KxDIM,2*(KyDIM-1),2*(KzDIM-1),nDE,nDE))
    IF(ALLOCATED(kNLfld3D).EQV..FALSE.) ALLOCATE(kNLfld3D(KxDIM,2*(KyDIM-1),2*(KzDIM-1),nDE))
      kfld3D=CMPLX(1.0d0,0.0d0)
     kLfld3D=CMPLX(0.0d0,0.0d0)
    kNLfld3D=CMPLX(0.0d0,0.0d0)
   END IF
  END IF

  IF(RSPACESOLVER) THEN
   IF((KxDIM.EQ.1).AND.(KyDIM.EQ.1)) THEN
    IF(ALLOCATED(  rfld1D).EQV..FALSE.) ALLOCATE(  rfld1D(2*(KzDIM-1),nDE))
    IF(ALLOCATED( rLfld1D).EQV..FALSE.) ALLOCATE( rLfld1D(2*(KzDIM-1),nDE,nDE))
    IF(ALLOCATED(rNLfld1D).EQV..FALSE.) ALLOCATE(rNLfld1D(2*(KzDIM-1),nDE))
    DO iDE=1,nDE
     CALL iFFTW(kfld1D(:,iDE),rfld1D(:,iDE))
    END DO
   ELSE IF((KyDIM.EQ.1).AND.(KzDIM.EQ.1)) THEN
    IF(ALLOCATED(  rfld1D).EQV..FALSE.) ALLOCATE(  rfld1D(2*(KxDIM-1),nDE))
    IF(ALLOCATED( rLfld1D).EQV..FALSE.) ALLOCATE( rLfld1D(2*(KxDIM-1),nDE,nDE))
    IF(ALLOCATED(rNLfld1D).EQV..FALSE.) ALLOCATE(rNLfld1D(2*(KxDIM-1),nDE))
    DO iDE=1,nDE
     CALL iFFTW(kfld1D(:,iDE),rfld1D(:,iDE))
    END DO
   ELSE IF((KxDIM.EQ.1).AND.(KzDIM.EQ.1)) THEN
    IF(ALLOCATED(  rfld1D).EQV..FALSE.) ALLOCATE(  rfld1D(2*(KyDIM-1),nDE))
    IF(ALLOCATED( rLfld1D).EQV..FALSE.) ALLOCATE( rLfld1D(2*(KyDIM-1),nDE,nDE))
    IF(ALLOCATED(rNLfld1D).EQV..FALSE.) ALLOCATE(rNLfld1D(2*(KyDIM-1),nDE))
    DO iDE=1,nDE
     CALL iFFTW(kfld1D(:,iDE),rfld1D(:,iDE))
    END DO
   ELSE IF(KxDIM.EQ.1) THEN
    IF(ALLOCATED(  rfld2D).EQV..FALSE.) ALLOCATE(  rfld2D(2*(KyDIM-1),2*(KzDIM-1),nDE))
    IF(ALLOCATED( rLfld2D).EQV..FALSE.) ALLOCATE( rLfld2D(2*(KyDIM-1),2*(KzDIM-1),nDE,nDE))
    IF(ALLOCATED(rNLfld2D).EQV..FALSE.) ALLOCATE(rNLfld2D(2*(KyDIM-1),2*(KzDIM-1),nDE))
    DO iDE=1,nDE
     CALL iFFTW(kfld2D(:,:,iDE),rfld2D(:,:,iDE))
    END DO
   ELSE IF(KyDIM.EQ.1) THEN
    IF(ALLOCATED(  rfld2D).EQV..FALSE.) ALLOCATE(  rfld2D(2*(KxDIM-1),2*(KzDIM-1),nDE))
    IF(ALLOCATED( rLfld2D).EQV..FALSE.) ALLOCATE( rLfld2D(2*(KxDIM-1),2*(KzDIM-1),nDE,nDE))
    IF(ALLOCATED(rNLfld2D).EQV..FALSE.) ALLOCATE(rNLfld2D(2*(KxDIM-1),2*(KzDIM-1),nDE))
    DO iDE=1,nDE
     CALL iFFTW(kfld2D(:,:,iDE),rfld2D(:,:,iDE))
    END DO
   ELSE IF(KzDIM.EQ.1) THEN
    IF(ALLOCATED(  rfld2D).EQV..FALSE.) ALLOCATE(  rfld2D(2*(KxDIM-1),2*(KyDIM-1),nDE))
    IF(ALLOCATED( rLfld2D).EQV..FALSE.) ALLOCATE( rLfld2D(2*(KxDIM-1),2*(KyDIM-1),nDE,nDE))
    IF(ALLOCATED(rNLfld2D).EQV..FALSE.) ALLOCATE(rNLfld2D(2*(KxDIM-1),2*(KyDIM-1),nDE))
    DO iDE=1,nDE
     CALL iFFTW(kfld2D(:,:,iDE),rfld2D(:,:,iDE))
    END DO
   ELSE
    IF(ALLOCATED(  rfld3D).EQV..FALSE.) ALLOCATE(  rfld3D(2*(KxDIM-1),2*(KyDIM-1),2*(KzDIM-1),nDE))
    IF(ALLOCATED( rLfld3D).EQV..FALSE.) ALLOCATE( rLfld3D(2*(KxDIM-1),2*(KyDIM-1),2*(KzDIM-1),nDE,nDE))
    IF(ALLOCATED(rNLfld3D).EQV..FALSE.) ALLOCATE(rNLfld3D(2*(KxDIM-1),2*(KyDIM-1),2*(KzDIM-1),nDE))
    DO iDE=1,nDE
     CALL iFFTW(kfld3D(:,:,:,iDE),rfld3D(:,:,:,iDE))
    END DO
   END IF
  END IF

  DO iXDIM=1,LXDIM
   DO iYDIM=1,LYDIM
    DO iZDIM=1,LZDIM
     IF(((KxDIM.EQ.1).AND.(KyDIM.EQ.1)).OR.((KyDIM.EQ.1).AND.(KzDIM.EQ.1)).OR.((KxDIM.EQ.1).AND.(KzDIM.EQ.1))) THEN
      IF(KSPACESOLVER) THEN
       CALL ODESET(ODESETID,kfld1D,kLfld1D,kNLfld1D,iXDIM,iYDIM,iZDIM)
      ELSE IF(RSPACESOLVER) THEN
       CALL ODESET(ODESETID,rfld1D,rLfld1D,rNLfld1D,iXDIM,iYDIM,iZDIM)
       DO iDE1=1,nDE
        DO iDE2=1,nDE
         CALL FFTW(rLfld1D(:,iDE1,iDE2),kLfld1D(:,iDE1,iDE2))
        END DO
       END DO
      END IF
     ELSE IF((KxDIM.EQ.1).OR.(KyDIM.EQ.1).OR.(KzDIM.EQ.1)) THEN
      IF(KSPACESOLVER) THEN
       CALL ODESET(ODESETID,kfld2D,kLfld2D,kNLfld2D,iXDIM,iYDIM,iZDIM)
      ELSE IF(RSPACESOLVER) THEN
       CALL ODESET(ODESETID,rfld2D,rLfld2D,rNLfld2D,iXDIM,iYDIM,iZDIM)
       DO iDE1=1,nDE
        DO iDE2=1,nDE
         CALL FFTW(rLfld2D(:,:,iDE1,iDE2),kLfld2D(:,:,iDE1,iDE2))
        END DO
       END DO
      END IF
     ELSE
      IF(KSPACESOLVER) THEN
       CALL ODESET(ODESETID,kfld3D,kLfld3D,kNLfld3D,iXDIM,iYDIM,iZDIM)
      ELSE IF(RSPACESOLVER) THEN
       CALL ODESET(ODESETID,rfld3D,rLfld3D,rNLfld3D,iXDIM,iYDIM,iZDIM)
       DO iDE1=1,nDE
        DO iDE2=1,nDE
         CALL FFTW(rLfld3D(:,:,:,iDE1,iDE2),kLfld3D(:,:,:,iDE1,iDE2))
        END DO
       END DO
      END IF
     END IF
 !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iKx,iKy,iKz,EIGENMAT)
 !$OMP DO SCHEDULE(DYNAMIC,10)
     DO iKx=1,KxDIM
      DO iKy=1,KyDIM
       DO iKz=1,KzDIM
        IF((KxDIM.EQ.1).AND.(KyDIM.EQ.1)) THEN
         EIGENMAT=CMPLX(0.0,1.0)*kLfld1D(iKz,:,:)
        ELSE IF((KyDIM.EQ.1).AND.(KzDIM.EQ.1)) THEN
         EIGENMAT=CMPLX(0.0,1.0)*kLfld1D(iKx,:,:)
        ELSE IF((KxDIM.EQ.1).AND.(KzDIM.EQ.1)) THEN
         EIGENMAT=CMPLX(0.0,1.0)*kLfld1D(iKy,:,:)
        ELSE IF(KxDIM.EQ.1) THEN
         EIGENMAT=CMPLX(0.0,1.0)*kLfld2D(iKy,iKz,:,:)
        ELSE IF(KyDIM.EQ.1) THEN
         EIGENMAT=CMPLX(0.0,1.0)*kLfld2D(iKx,iKz,:,:)
        ELSE IF(KzDIM.EQ.1) THEN
         EIGENMAT=CMPLX(0.0,1.0)*kLfld2D(iKx,iKy,:,:)
        ELSE
         EIGENMAT=CMPLX(0.0,1.0)*kLfld3D(iKx,iKy,iKz,:,:)
        END IF
       !FINDING THE EIGENVALUES AND EIGENVECTORS USING EISPACK LIBRARY
       CALL EPKEIGENSEVAL(EIGENMAT,LEIGENVAL4D(iKx,iKy,iKz,:),LEIGENVEC5D(iKx,iKy,iKz,:,:))
       END DO
      END DO
     END DO
 !$OMP END DO
 !$OMP END PARALLEL

     CALL EIGENSMATADJUST(LEIGENVAL4D,LEIGENVEC5D)
     CALL EIGENSMATADJUST(LEIGENVAL4D,LEIGENVEC5D)

     LMAXGINDEX=MAXLOC(AIMAG(LEIGENVAL4D))
     EIGENVAL7D(iXDIM,iYDIM,iZDIM,:,:,:,1)=LEIGENVAL4D(:,:,:,LMAXGINDEX(4))
     EIGENVEC7D(iXDIM,iYDIM,iZDIM,:,:,:,:)=LEIGENVEC5D(:,:,:,:,LMAXGINDEX(4))
     KxGMAX3D(iXDIM,iYDIM,iZDIM)=LMAXGINDEX(1)
     KyGMAX3D(iXDIM,iYDIM,iZDIM)=LMAXGINDEX(2)
     KzGMAX3D(iXDIM,iYDIM,iZDIM)=LMAXGINDEX(3)
    END DO
   END DO
  END DO
 !GOLINEAR=iGOLINEAR
 !GONONLINEAR=iGONONLINEAR

  PRINT*, "EIGENS CALCULATION >>> FINISHED!"
  RETURN
 END SUBROUTINE r8SYSEIGENVV

 SUBROUTINE r4r8c8SYSEIGENFLDS(ODESETID,nDE,FLDTYPE,LOCALTYPE,LX3D,LY3D,LZ3D,rfld,ikfld)
  IMPLICIT NONE
  REAL(4),            INTENT(IN)  :: LX3D(:,:,:),LY3D(:,:,:),LZ3D(:,:,:)
  REAL(8),            INTENT(OUT) :: rfld(:,:,:,:)
  COMPLEX(8),OPTIONAL,INTENT(OUT) :: ikfld(:,:,:,:)
  COMPLEX(8),         ALLOCATABLE ::  kfld(:,:,:,:)
  INTEGER(4),         INTENT(IN)  :: nDE,ODESETID
  CHARACTER(*),       INTENT(IN)  :: FLDTYPE
  CHARACTER(*),       INTENT(IN)  :: LOCALTYPE(:)

  CALL SYSEIGENFLDS(ODESETID,nDE,FLDTYPE,LOCALTYPE,1.0d0*LX3D,1.0d0*LY3D,1.0d0*LZ3D,rfld,ikfld)

  RETURN
 END SUBROUTINE r4r8c8SYSEIGENFLDS

 SUBROUTINE r8r8c8SYSEIGENFLDS(ODESETID,nDE,FLDTYPE,LOCALTYPE,LX3D,LY3D,LZ3D,rfld,ikfld)
  IMPLICIT NONE
  REAL(8),            INTENT(IN)  :: LX3D(:,:,:),LY3D(:,:,:),LZ3D(:,:,:)
  REAL(8),            INTENT(OUT) :: rfld(:,:,:,:)
  COMPLEX(8),OPTIONAL,INTENT(OUT) :: ikfld(:,:,:,:)
  COMPLEX(8),         ALLOCATABLE ::  kfld(:,:,:,:)
  INTEGER(4),         INTENT(IN)  :: nDE,ODESETID
  CHARACTER(*),       INTENT(IN)  :: FLDTYPE
  CHARACTER(*),       INTENT(IN)  :: LOCALTYPE(:)

  CHARACTER(10)                   :: EIGENLOCALTYPE(3)
  INTEGER(4)                      :: LXDIM,LYDIM,LZDIM
  INTEGER(4)                      :: KxDIM,KyDIM,KzDIM
  INTEGER(4)                      :: nKx,nKy,nKz
  INTEGER(4)                      :: iKx,iKy,iKz
  INTEGER(4)                      :: iKxMAX,iKyMAX,iKzMAX
  INTEGER(4)                      :: iXDIM,iYDIM,iZDIM
  INTEGER(4)                      :: iD1,iD2,iCOUNTER,fldID
  INTEGER(4)                      :: iDE,iDE1,iDE2
  INTEGER(4)                      :: EIGENID
  INTEGER(4)                      :: MATDIMS3D(3)

  MATDIMS3D=SHAPE(LX3D)

  IF(SIZE(LOCALTYPE).EQ.1) THEN
   EIGENLOCALTYPE(1)=TRIM(LOCALTYPE(1))
   EIGENLOCALTYPE(2)=TRIM(LOCALTYPE(1))
   EIGENLOCALTYPE(3)=TRIM(LOCALTYPE(1))
  ELSE IF(SIZE(LOCALTYPE).EQ.2) THEN
   EIGENLOCALTYPE(1)=TRIM(LOCALTYPE(1))
   EIGENLOCALTYPE(2)=TRIM(LOCALTYPE(2))
   EIGENLOCALTYPE(3)=TRIM(LOCALTYPE(2))
  ELSE IF(SIZE(LOCALTYPE).EQ.3) THEN
   EIGENLOCALTYPE(1)=TRIM(LOCALTYPE(1))
   EIGENLOCALTYPE(2)=TRIM(LOCALTYPE(2))
   EIGENLOCALTYPE(3)=TRIM(LOCALTYPE(3))
  ELSE IF(SIZE(LOCALTYPE).GT.3) THEN
   PRINT*, "WARNING: LOCALTYPE DIMENSIONS IS GREATER THAN THE SYSTEM DIMESNIONS"
   EIGENLOCALTYPE(1)=TRIM(LOCALTYPE(1))
   EIGENLOCALTYPE(2)=TRIM(LOCALTYPE(2))
   EIGENLOCALTYPE(3)=TRIM(LOCALTYPE(3))
  END IF

  IF(to_LowerCase(TRIM(EIGENLOCALTYPE(1))).EQ."") THEN
   LXDIM=1
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(1))).EQ."local") THEN
   LXDIM=1
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(1))).EQ."nonlocal") THEN
   LXDIM=MATDIMS3D(1)
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(1))).EQ."quasilocal") THEN
   LXDIM=MATDIMS3D(1)
  END IF

  IF(to_LowerCase(TRIM(EIGENLOCALTYPE(2))).EQ."") THEN
   LYDIM=1
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(2))).EQ."local") THEN
   LYDIM=1
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(2))).EQ."nonlocal") THEN
   LYDIM=MATDIMS3D(2)
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(2))).EQ."quasilocal") THEN
   LYDIM=MATDIMS3D(2)
  END IF

  IF(to_LowerCase(TRIM(EIGENLOCALTYPE(3))).EQ."") THEN
   LZDIM=1
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(3))).EQ."local") THEN
   LZDIM=1
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(3))).EQ."nonlocal") THEN
   LZDIM=MATDIMS3D(3)
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(3))).EQ."quasilocal") THEN
   LZDIM=MATDIMS3D(3)
  END IF

  CALL SYSEIGENVV(ODESETID,nDE,LOCALTYPE,LX3D,LY3D,LZ3D)
  KxDIM=SIZE(eiKx1D)
  KyDIM=SIZE(eiKy1D)
  KzDIM=SIZE(eiKz1D)

  IF((KyDIM.EQ.1).AND.(KzDIM.EQ.1)) THEN
   IF(ALLOCATED(kfld).EQV..FALSE.) ALLOCATE(kfld(KxDIM,KyDIM,KzDIM,nDE))
  ELSE IF((KxDIM.EQ.1).AND.(KzDIM.EQ.1)) THEN
   IF(ALLOCATED(kfld).EQV..FALSE.) ALLOCATE(kfld(KxDIM,KyDIM,KzDIM,nDE))
  ELSE IF((KxDIM.EQ.1).AND.(KyDIM.EQ.1)) THEN
   IF(ALLOCATED(kfld).EQV..FALSE.) ALLOCATE(kfld(KxDIM,KyDIM,KzDIM,nDE))
  ELSE IF(KxDIM.EQ.1) THEN
   IF(ALLOCATED(kfld).EQV..FALSE.) ALLOCATE(kfld(KxDIM,KyDIM,2*(KzDIM-1),nDE))
  ELSE IF(KyDIM.EQ.1) THEN
   IF(ALLOCATED(kfld).EQV..FALSE.) ALLOCATE(kfld(KxDIM,KyDIM,2*(KzDIM-1),nDE))
  ELSE IF(KzDIM.EQ.1) THEN
   IF(ALLOCATED(kfld).EQV..FALSE.) ALLOCATE(kfld(KxDIM,2*(KyDIM-1),KzDIM,nDE))
  ELSE
   IF(ALLOCATED(kfld).EQV..FALSE.) ALLOCATE(kfld(KxDIM,2*(KyDIM-1),2*(KzDIM-1),nDE))
  END IF
  kfld=CMPLX(0.0d0,0.0d0)

  PRINT*,"FIELDS CALCULATIONS >>> STARTED!"
  DO iXDIM=1,LXDIM
   DO iYDIM=1,LYDIM
    DO iZDIM=1,LZDIM
     iKxMAX=KxGMAX3D(iXDIM,iYDIM,iZDIM)
     iKyMAX=KyGMAX3D(iXDIM,iYDIM,iZDIM)
     iKzMAX=KzGMAX3D(iXDIM,iYDIM,iZDIM)
     IF(to_LowerCase(TRIM(FLDTYPE)).EQ."max-mode") THEN
!******* INITIALIZE THE FIELDS WITH ONLY THE MAXIMUM GROWING MODE *********
      DO fldID=1,nDE
       kfld(iKxMAX,iKyMAX,iKzMAX,fldID)=1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKxMAX,iKyMAX,iKzMAX,fldID)
      END DO
     ELSE IF(to_LowerCase(TRIM(FLDTYPE)).EQ."positive-max-modes") THEN
!********* INITIALIZE THE FIELDS WITH ALL POSITIVE GROWING MODES **********
      IF((KyDIM.EQ.1).AND.(KzDIM.EQ.1)) THEN
       DO iKx=1,KxDIM
        IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,iKx,1,1,1)).GT.0.0d0) THEN
         DO fldID=1,nDE
          kfld(iKx,1,1,fldID)=kfld(iKx,1,1,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,1,1,fldID)
         END DO
        END IF
       END DO
      ELSE IF((KxDIM.EQ.1).AND.(KzDIM.EQ.1)) THEN
       DO iKy=1,KyDIM
        IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,1,iKy,1,1)).GT.0.0d0) THEN
         DO fldID=1,nDE
          kfld(1,iKy,1,fldID)=kfld(1,iKy,1,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,1,iKy,1,fldID)
         END DO
        END IF
       END DO
      ELSE IF((KxDIM.EQ.1).AND.(KyDIM.EQ.1)) THEN
       DO iKz=1,KzDIM
        IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,1,1,iKz,1)).GT.0.0d0) THEN
         DO fldID=1,nDE
          kfld(1,1,iKz,fldID)=kfld(1,1,iKz,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,1,1,iKz,fldID)
         END DO
        END IF
       END DO
      ELSE IF(KxDIM.EQ.1) THEN
       DO iKy=1,KyDIM
        IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,1,iKy,1,1)).GT.0.0d0) THEN
         DO fldID=1,nDE
          kfld(iKxMAX,iKy,1,fldID)=kfld(1,iKy,1,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,1,iKy,1,fldID)
         END DO
        END IF
        DO iKz=2,KzDIM
         IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,1,iKy,iKz,1)).GT.0.0d0) THEN
          DO fldID=1,nDE
           kfld(iKxMAX,iKy,iKz,fldID)=kfld(1,iKy,iKz,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,1,iKy,iKz,fldID)
           kfld(1,iKy,2*KzDIM-iKz,fldID)=kfld(1,iKy,2*KzDIM-iKz,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,1,iKy,iKz,fldID)
          END DO
         END IF
        END DO
       END DO
      ELSE IF(KyDIM.EQ.1) THEN
       DO iKx=1,KxDIM
        IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,iKx,1,1,1)).GT.0.0d0) THEN
         DO fldID=1,nDE
          kfld(iKx,1,1,fldID)=kfld(iKx,1,1,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,1,1,fldID)
         END DO
        END IF
        DO iKz=2,KzDIM
         IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,iKx,1,iKz,1)).GT.0.0d0) THEN
          DO fldID=1,nDE
           kfld(iKx,1,iKz,fldID)=kfld(iKx,iKyMAX,iKz,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,iKyMAX,iKz,fldID)
           kfld(iKx,1,2*KzDIM-iKz,fldID)=kfld(iKx,1,2*KzDIM-iKz,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,1,iKz,fldID)
          END DO
         END IF
        END DO
       END DO
      ELSE IF(KzDIM.EQ.1) THEN
       DO iKx=1,KxDIM
        IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,iKx,1,1,1)).GT.0.0d0) THEN
         DO fldID=1,nDE
          kfld(iKx,1,1,fldID)=kfld(iKx,1,1,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,1,1,fldID)
         END DO
        END IF
        DO iKy=2,KyDIM
         IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,iKx,iKy,1,1)).GT.0.0d0) THEN
          DO fldID=1,nDE
           kfld(iKx,iKy,1,fldID)=kfld(iKx,iKy,1,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,iKy,1,fldID)
           kfld(iKx,2*KyDIM-iKy,1,fldID)=kfld(iKx,2*KyDIM-iKy,1,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,iKy,1,fldID)
          END DO
         END IF
        END DO
       END DO
      ELSE
       DO iKx=1,KxDIM
        IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,iKx,1,1,1)).GT.0.0d0) THEN
         DO fldID=1,nDE
          kfld(iKx,1,1,fldID)=kfld(iKx,1,1,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,1,1,fldID)
         END DO
        END IF
        DO iKy=2,KyDIM
         IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,iKx,iKy,iKzMAX,1)).GT.0.0d0) THEN
          DO fldID=1,nDE
           kfld(iKx,iKy,iKzMAX,fldID)=kfld(iKx,iKy,iKzMAX,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,iKy,iKzMAX,fldID)
           kfld(iKx,2*KyDIM-iKy,iKzMAX,fldID)=kfld(iKx,2*KyDIM-iKy,iKzMAX,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,iKy,iKzMAX,fldID)
          END DO
         END IF
        END DO
        DO iKz=2,KzDIM
         IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,iKx,iKyMAX,iKz,1)).GT.0.0d0) THEN
          DO fldID=1,nDE
           kfld(iKx,iKyMAX,iKz,fldID)=kfld(iKx,iKyMAX,iKz,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,iKyMAX,iKz,fldID)
           kfld(iKx,iKyMAX,2*KzDIM-iKz,fldID)=kfld(iKx,iKyMAX,2*KzDIM-iKz,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,iKyMAX,iKz,fldID)
          END DO
         END IF
        END DO
       END DO
      END IF
     END IF
    END DO
   END DO
  END DO

  IF((KyDIM.EQ.1).AND.(KzDIM.EQ.1)) THEN
   DO iDE=1,nDE
    CALL iFFTW(kfld(:,1,1,iDE),rfld(:,1,1,iDE))
   END DO
  ELSE IF((KxDIM.EQ.1).AND.(KzDIM.EQ.1)) THEN
   DO iDE=1,nDE
    CALL iFFTW(kfld(1,:,1,iDE),rfld(1,:,1,iDE))
   END DO
  ELSE IF((KxDIM.EQ.1).AND.(KyDIM.EQ.1)) THEN
   DO iDE=1,nDE
    CALL iFFTW(kfld(1,1,:,iDE),rfld(1,1,:,iDE))
   END DO
  ELSE IF(KxDIM.EQ.1) THEN
   DO iDE=1,nDE
    CALL iFFTW(kfld(1,:,:,iDE),rfld(1,:,:,iDE))
   END DO
  ELSE IF(KyDIM.EQ.1) THEN
   DO iDE=1,nDE
    CALL iFFTW(kfld(:,1,:,iDE),rfld(:,1,:,iDE))
   END DO
  ELSE IF(KzDIM.EQ.1) THEN
   DO iDE=1,nDE
    CALL iFFTW(kfld(:,:,1,iDE),rfld(:,:,1,iDE))
   END DO
  ELSE
   DO iDE=1,nDE
    CALL iFFTW(kfld(:,:,:,iDE),rfld(:,:,:,iDE))
   END DO
  END IF

  IF(PRESENT(ikfld)) ikfld=kfld
  PRINT*,"FIELDS CALCULATIONS >>> FINISHED!"

  RETURN
 END SUBROUTINE r8r8c8SYSEIGENFLDS

 SUBROUTINE r4c8r8SYSEIGENFLDS(ODESETID,nDE,FLDTYPE,LOCALTYPE,LX3D,LY3D,LZ3D,kfld,rfld)
  IMPLICIT NONE
  REAL(4),         INTENT(IN)  :: LX3D(:,:,:),LY3D(:,:,:),LZ3D(:,:,:)
  REAL(8),OPTIONAL,INTENT(OUT) :: rfld(:,:,:,:)
  COMPLEX(8),      INTENT(OUT) :: kfld(:,:,:,:)
  INTEGER(4),      INTENT(IN)  :: nDE,ODESETID
  CHARACTER(*),    INTENT(IN)  :: FLDTYPE
  CHARACTER(*),    INTENT(IN)  :: LOCALTYPE(:)

  CALL SYSEIGENFLDS(ODESETID,nDE,FLDTYPE,LOCALTYPE,1.0d0*LX3D,1.0d0*LY3D,1.0d0*LZ3D,kfld,rfld)

  RETURN
 END SUBROUTINE r4c8r8SYSEIGENFLDS

 SUBROUTINE r8c8r8SYSEIGENFLDS(ODESETID,nDE,FLDTYPE,LOCALTYPE,LX3D,LY3D,LZ3D,kfld,rfld)
  IMPLICIT NONE
  REAL(8),         INTENT(IN)  :: LX3D(:,:,:),LY3D(:,:,:),LZ3D(:,:,:)
  REAL(8),OPTIONAL,INTENT(OUT) :: rfld(:,:,:,:)
  COMPLEX(8),      INTENT(OUT) :: kfld(:,:,:,:)
  INTEGER(4),      INTENT(IN)  :: nDE,ODESETID
  CHARACTER(*),    INTENT(IN)  :: FLDTYPE
  CHARACTER(*),    INTENT(IN)  :: LOCALTYPE(:)

  CHARACTER(10)                :: EIGENLOCALTYPE(3)
  INTEGER(4)                   :: LXDIM,LYDIM,LZDIM
  INTEGER(4)                   :: KxDIM,KyDIM,KzDIM
  INTEGER(4)                   :: nKx,nKy,nKz
  INTEGER(4)                   :: iKx,iKy,iKz
  INTEGER(4)                   :: iKxMAX,iKyMAX,iKzMAX
  INTEGER(4)                   :: iXDIM,iYDIM,iZDIM
  INTEGER(4)                   :: iD1,iD2,iCOUNTER,fldID
  INTEGER(4)                   :: iDE,iDE1,iDE2
  INTEGER(4)                   :: EIGENID
  INTEGER(4)                   :: MATDIMS3D(3)

  MATDIMS3D=SHAPE(LX3D)

  IF(SIZE(LOCALTYPE).EQ.1) THEN
   EIGENLOCALTYPE(1)=TRIM(LOCALTYPE(1))
   EIGENLOCALTYPE(2)=TRIM(LOCALTYPE(1))
   EIGENLOCALTYPE(3)=TRIM(LOCALTYPE(1))
  ELSE IF(SIZE(LOCALTYPE).EQ.2) THEN
   EIGENLOCALTYPE(1)=TRIM(LOCALTYPE(1))
   EIGENLOCALTYPE(2)=TRIM(LOCALTYPE(2))
   EIGENLOCALTYPE(3)=TRIM(LOCALTYPE(2))
  ELSE IF(SIZE(LOCALTYPE).EQ.3) THEN
   EIGENLOCALTYPE(1)=TRIM(LOCALTYPE(1))
   EIGENLOCALTYPE(2)=TRIM(LOCALTYPE(2))
   EIGENLOCALTYPE(3)=TRIM(LOCALTYPE(3))
  ELSE IF(SIZE(LOCALTYPE).GT.3) THEN
   PRINT*, "WARNING: LOCALTYPE DIMENSIONS IS GREATER THAN THE SYSTEM DIMESNIONS"
   EIGENLOCALTYPE(1)=TRIM(LOCALTYPE(1))
   EIGENLOCALTYPE(2)=TRIM(LOCALTYPE(2))
   EIGENLOCALTYPE(3)=TRIM(LOCALTYPE(3))
  END IF

  IF(to_LowerCase(TRIM(EIGENLOCALTYPE(1))).EQ."") THEN
   LXDIM=1
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(1))).EQ."local") THEN
   LXDIM=1
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(1))).EQ."nonlocal") THEN
   LXDIM=MATDIMS3D(1)
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(1))).EQ."quasilocal") THEN
   LXDIM=MATDIMS3D(1)
  END IF

  IF(to_LowerCase(TRIM(EIGENLOCALTYPE(2))).EQ."") THEN
   LYDIM=1
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(2))).EQ."local") THEN
   LYDIM=1
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(2))).EQ."nonlocal") THEN
   LYDIM=MATDIMS3D(2)
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(2))).EQ."quasilocal") THEN
   LYDIM=MATDIMS3D(2)
  END IF

  IF(to_LowerCase(TRIM(EIGENLOCALTYPE(3))).EQ."") THEN
   LZDIM=1
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(3))).EQ."local") THEN
   LZDIM=1
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(3))).EQ."nonlocal") THEN
   LZDIM=MATDIMS3D(3)
  ELSE IF(to_LowerCase(TRIM(EIGENLOCALTYPE(3))).EQ."quasilocal") THEN
   LZDIM=MATDIMS3D(3)
  END IF

  CALL SYSEIGENVV(ODESETID,nDE,LOCALTYPE,LX3D,LY3D,LZ3D)

  kfld=CMPLX(0.0d0,0.0d0)

  DO iXDIM=1,LXDIM
   DO iYDIM=1,LYDIM
    DO iZDIM=1,LZDIM
     iKxMAX=KxGMAX3D(iXDIM,iYDIM,iZDIM)
     iKyMAX=KyGMAX3D(iXDIM,iYDIM,iZDIM)
     iKzMAX=KzGMAX3D(iXDIM,iYDIM,iZDIM)
     IF(to_LowerCase(TRIM(FLDTYPE)).EQ."max-mode") THEN
!******* INITIALIZE THE FIELDS WITH ONLY THE MAXIMUM GROWING MODE *********
      DO fldID=1,nDE
       kfld(iKxMAX,iKyMAX,iKzMAX,fldID)=1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKxMAX,iKyMAX,iKzMAX,fldID)
      END DO
     ELSE IF(to_LowerCase(TRIM(FLDTYPE)).EQ."positive-max-modes") THEN
!********* INITIALIZE THE FIELDS WITH ALL POSITIVE GROWING MODES **********
      IF((KyDIM.EQ.1).AND.(KzDIM.EQ.1)) THEN
       DO iKx=1,KxDIM
        IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,iKx,1,1,1)).GT.0.0d0) THEN
         DO fldID=1,nDE
          kfld(iKx,1,1,fldID)=kfld(iKx,1,1,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,1,1,fldID)
         END DO
        END IF
       END DO
      ELSE IF((KxDIM.EQ.1).AND.(KzDIM.EQ.1)) THEN
       DO iKy=1,KyDIM
        IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,1,iKy,1,1)).GT.0.0d0) THEN
         DO fldID=1,nDE
          kfld(1,iKy,1,fldID)=kfld(1,iKy,1,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,1,iKy,1,fldID)
         END DO
        END IF
       END DO
      ELSE IF((KxDIM.EQ.1).AND.(KyDIM.EQ.1)) THEN
       DO iKz=1,KzDIM
        IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,1,1,iKz,1)).GT.0.0d0) THEN
         DO fldID=1,nDE
          kfld(1,1,iKz,fldID)=kfld(1,1,iKz,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,1,1,iKz,fldID)
         END DO
        END IF
       END DO
      ELSE IF(KxDIM.EQ.1) THEN
       DO iKy=1,KyDIM
        IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,1,iKy,1,1)).GT.0.0d0) THEN
         DO fldID=1,nDE
          kfld(iKxMAX,iKy,1,fldID)=kfld(1,iKy,1,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,1,iKy,1,fldID)
         END DO
        END IF
        DO iKz=2,KzDIM
         IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,1,iKy,iKz,1)).GT.0.0d0) THEN
          DO fldID=1,nDE
           kfld(iKxMAX,iKy,iKz,fldID)=kfld(1,iKy,iKz,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,1,iKy,iKz,fldID)
           kfld(1,iKy,2*KzDIM-iKz,fldID)=kfld(1,iKy,2*KzDIM-iKz,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,1,iKy,iKz,fldID)
          END DO
         END IF
        END DO
       END DO
      ELSE IF(KyDIM.EQ.1) THEN
       DO iKx=1,KxDIM
        IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,iKx,1,1,1)).GT.0.0d0) THEN
         DO fldID=1,nDE
          kfld(iKx,1,1,fldID)=kfld(iKx,1,1,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,1,1,fldID)
         END DO
        END IF
        DO iKz=2,KzDIM
         IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,iKx,1,iKz,1)).GT.0.0d0) THEN
          DO fldID=1,nDE
           kfld(iKx,1,iKz,fldID)=kfld(iKx,iKyMAX,iKz,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,iKyMAX,iKz,fldID)
           kfld(iKx,1,2*KzDIM-iKz,fldID)=kfld(iKx,1,2*KzDIM-iKz,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,1,iKz,fldID)
          END DO
         END IF
        END DO
       END DO
      ELSE IF(KzDIM.EQ.1) THEN
       DO iKx=1,KxDIM
        IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,iKx,1,1,1)).GT.0.0d0) THEN
         DO fldID=1,nDE
          kfld(iKx,1,1,fldID)=kfld(iKx,1,1,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,1,1,fldID)
         END DO
        END IF
        DO iKy=2,KyDIM
         IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,iKx,iKy,1,1)).GT.0.0d0) THEN
          DO fldID=1,nDE
           kfld(iKx,iKy,1,fldID)=kfld(iKx,iKy,1,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,iKy,1,fldID)
           kfld(iKx,2*KyDIM-iKy,1,fldID)=kfld(iKx,2*KyDIM-iKy,1,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,iKy,1,fldID)
          END DO
         END IF
        END DO
       END DO
      ELSE
       DO iKx=1,KxDIM
        IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,iKx,1,1,1)).GT.0.0d0) THEN
         DO fldID=1,nDE
          kfld(iKx,1,1,fldID)=kfld(iKx,1,1,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,1,1,fldID)
         END DO
        END IF
        DO iKy=2,KyDIM
         IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,iKx,iKy,iKzMAX,1)).GT.0.0d0) THEN
          DO fldID=1,nDE
           kfld(iKx,iKy,iKzMAX,fldID)=kfld(iKx,iKy,iKzMAX,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,iKy,iKzMAX,fldID)
           kfld(iKx,2*KyDIM-iKy,iKzMAX,fldID)=kfld(iKx,2*KyDIM-iKy,iKzMAX,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,iKy,iKzMAX,fldID)
          END DO
         END IF
        END DO
        DO iKz=2,KzDIM
         IF(AIMAG(EIGENVAL7D(iXDIM,iYDIM,iZDIM,iKx,iKyMAX,iKz,1)).GT.0.0d0) THEN
          DO fldID=1,nDE
           kfld(iKx,iKyMAX,iKz,fldID)=kfld(iKx,iKyMAX,iKz,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,iKyMAX,iKz,fldID)
           kfld(iKx,iKyMAX,2*KzDIM-iKz,fldID)=kfld(iKx,iKyMAX,2*KzDIM-iKz,fldID)+1.0d0*EIGENVEC7D(iXDIM,iYDIM,iZDIM,iKx,iKyMAX,iKz,fldID)
          END DO
         END IF
        END DO
       END DO
      END IF
     END IF
    END DO
   END DO
  END DO

  IF((KyDIM.EQ.1).AND.(KzDIM.EQ.1)) THEN
   IF(PRESENT(rfld)) THEN
    DO iDE=1,nDE
     CALL iFFTW(kfld(:,1,1,iDE),rfld(:,1,1,iDE))
    END DO
   END IF
  ELSE IF((KxDIM.EQ.1).AND.(KzDIM.EQ.1)) THEN
   IF(PRESENT(rfld)) THEN
    DO iDE=1,nDE
     CALL iFFTW(kfld(1,:,1,iDE),rfld(1,:,1,iDE))
    END DO
   END IF
  ELSE IF((KxDIM.EQ.1).AND.(KyDIM.EQ.1)) THEN
   IF(PRESENT(rfld)) THEN
    DO iDE=1,nDE
     CALL iFFTW(kfld(1,1,:,iDE),rfld(1,1,:,iDE))
    END DO
   END IF
  ELSE IF(KxDIM.EQ.1) THEN
   IF(PRESENT(rfld)) THEN
    DO iDE=1,nDE
     CALL iFFTW(kfld(1,:,:,iDE),rfld(1,:,:,iDE))
    END DO
   END IF
  ELSE IF(KyDIM.EQ.1) THEN
   IF(PRESENT(rfld)) THEN
    DO iDE=1,nDE
     CALL iFFTW(kfld(:,1,:,iDE),rfld(:,1,:,iDE))
    END DO
   END IF
  ELSE IF(KzDIM.EQ.1) THEN
   IF(PRESENT(rfld)) THEN
    DO iDE=1,nDE
     CALL iFFTW(kfld(:,:,1,iDE),rfld(:,:,1,iDE))
    END DO
   END IF
  ELSE
   IF(PRESENT(rfld)) THEN
    DO iDE=1,nDE
     CALL iFFTW(kfld(:,:,:,iDE),rfld(:,:,:,iDE))
    END DO
   END IF
  END IF

  RETURN
 END SUBROUTINE r8c8r8SYSEIGENFLDS

!************************************************************************************************
! A SUBROUTINE TO FIND THE EIGENVALUES AND EIGEVECTORS OF A SQUARE MATRIX WITHOUT SINGULARITIES
!************************************************************************************************
! USING EISPACK LIBRARY
!************************************************************************************************
 SUBROUTINE r4EPKEIGENSEVAL(TINMAT,TEIGENVAL,TEIGENVEC)
  IMPLICIT NONE
  COMPLEX(4),INTENT(IN)           :: TINMAT(:,:)
  COMPLEX(4),INTENT(OUT)          :: TEIGENVAL(:)
  COMPLEX(4),INTENT(OUT),OPTIONAL :: TEIGENVEC(:,:)
  INTEGER(4)                      :: WHICH,IERR,MM,NN
  INTEGER(4)                      :: MATDIMS(2)
  REAL(4),ALLOCATABLE             :: RINMAT(:,:),IINMAT(:,:)
  REAL(4),ALLOCATABLE             :: IEIGENVAL(:),REIGENVAL(:)
  REAL(4),ALLOCATABLE             :: IEIGENVEC(:,:),REIGENVEC(:,:)

  MATDIMS=SHAPE(TINMAT)
  MM=MATDIMS(1)
  NN=MATDIMS(2)
  IF(PRESENT(TEIGENVEC)) THEN
    WHICH=1
  ELSE
    WHICH=0
  END IF
  ALLOCATE(RINMAT(MM,NN),IINMAT(MM,NN))
  ALLOCATE(REIGENVAL(NN),IEIGENVAL(NN))
  IF(WHICH.EQ.1) ALLOCATE(REIGENVEC(MM,NN),IEIGENVEC(MM,NN))
  RINMAT=REAL(TINMAT)
  IINMAT=AIMAG(TINMAT)
  CALL CG(MM,NN,RINMAT,IINMAT,REIGENVAL,IEIGENVAL,WHICH,REIGENVEC,IEIGENVEC,IERR)
  IF(IERR.NE.0) THEN
    WRITE(*,*)  "ERROR IN CALCULATING THE EIGENVALUES AND EIGENVECTORS"
    RETURN
  END IF
  TEIGENVAL=CMPLX(REIGENVAL,IEIGENVAL)
  IF(WHICH.EQ.1) TEIGENVEC=CMPLX(REIGENVEC,IEIGENVEC)

  RETURN
 END SUBROUTINE r4EPKEIGENSEVAL

 SUBROUTINE r8EPKEIGENSEVAL(TINMAT,TEIGENVAL,TEIGENVEC)
  IMPLICIT NONE
  COMPLEX(8),INTENT(IN)           :: TINMAT(:,:)
  COMPLEX(8),INTENT(OUT)          :: TEIGENVAL(:)
  COMPLEX(8),INTENT(OUT),OPTIONAL :: TEIGENVEC(:,:)
  INTEGER(4)                      :: WHICH,IERR,MM,NN
  INTEGER(4)                      :: MATDIMS(2)
  REAL(8),ALLOCATABLE             :: RINMAT(:,:),IINMAT(:,:)
  REAL(8),ALLOCATABLE             :: IEIGENVAL(:),REIGENVAL(:)
  REAL(8),ALLOCATABLE             :: IEIGENVEC(:,:),REIGENVEC(:,:)

  MATDIMS=SHAPE(TINMAT)
  MM=MATDIMS(1)
  NN=MATDIMS(2)
  IF(PRESENT(TEIGENVEC)) THEN
    WHICH=1
  ELSE
    WHICH=0
  END IF
  IF(ALLOCATED(RINMAT)) DEALLOCATE(RINMAT)
  IF(ALLOCATED(REIGENVAL)) DEALLOCATE(REIGENVAL)
  IF(ALLOCATED(REIGENVEC)) DEALLOCATE(REIGENVEC)
  ALLOCATE(RINMAT(MM,NN),IINMAT(MM,NN))
  ALLOCATE(REIGENVAL(NN),IEIGENVAL(NN))
  IF(WHICH.EQ.1) ALLOCATE(REIGENVEC(MM,NN),IEIGENVEC(MM,NN))
  RINMAT=REAL(TINMAT)
  IINMAT=AIMAG(TINMAT)
! CALL CG(MM,NN,RINMAT,IINMAT,REIGENVAL,IEIGENVAL,WHICH,REIGENVEC,IEIGENVEC,IERR)
  IF(IERR.NE.0) THEN
    WRITE(*,*)  "ERROR IN CALCULATING THE EIGENVALUES AND EIGENVECTORS"
    RETURN
  END IF
  TEIGENVAL=REIGENVAL+CMPLX(0.0,1.0)*IEIGENVAL
  IF(WHICH.EQ.1) TEIGENVEC=REIGENVEC+CMPLX(0.0,1.0)*IEIGENVEC

  RETURN
 END SUBROUTINE r8EPKEIGENSEVAL


!******************************************************************************************
! A SUBROUTINE TO ADJUST THE EIGENVALUES, EIGENVECTORS OF THE SAME MODE AT THE SAME ARRAY
!******************************************************************************************
 SUBROUTINE EIGENSMATADJUST(EIGENVAL3D,EIGENVEC3D)
  IMPLICIT NONE
  REAL(4)                :: TEMPRNUM
  COMPLEX(4)             :: TEMPCNUM
  COMPLEX(4)             :: EIGENVAL3D(:,:,:,:),EIGENVEC3D(:,:,:,:,:)
  COMPLEX(4),ALLOCATABLE :: TEMPCMAT(:)
  INTEGER(4)             :: TEMPINUM
  INTEGER(4)             :: nKx,nKy,nKz,iKx,iKy,iKz
  INTEGER(4)             :: nfld,MATDIMS4D(4)

  MATDIMS4D=SHAPE(EIGENVAL3D)
  nKx =MATDIMS4D(1)
  nKy =MATDIMS4D(2)
  nKz =MATDIMS4D(3)
  nfld=MATDIMS4D(4)

  IF(ALLOCATED(TEMPCMAT)) DEALLOCATE(TEMPCMAT)
  ALLOCATE(TEMPCMAT(nfld))

  DO iKx=1,nKx
   DO iKz=1,nKz
    DO iKy=1,nKy
     TEMPRNUM=MAXVAL(AIMAG(EIGENVAL3D(iKx,iKy,iKz,:)))
     TEMPINUM=MAXLOC(AIMAG(EIGENVAL3D(iKx,iKy,iKz,:)),1)
     IF(TEMPINUM.NE.1) THEN
      TEMPCNUM=EIGENVAL3D(iKx,iKy,iKz,1)
      EIGENVAL3D(iKx,iKy,iKz,1)=EIGENVAL3D(iKx,iKy,iKz,TEMPINUM)
      EIGENVAL3D(iKx,iKy,iKz,TEMPINUM)=TEMPCNUM
      TEMPCMAT=EIGENVEC3D(iKx,iKy,iKz,:,1)
      EIGENVEC3D(iKx,iKy,iKz,:,1)=EIGENVEC3D(iKx,iKy,iKz,:,TEMPINUM)
      EIGENVEC3D(iKx,iKy,iKz,:,TEMPINUM)=TEMPCMAT
     END IF
    END DO
   END DO
  END DO

  RETURN
 END SUBROUTINE EIGENSMATADJUST



!********************************************************
!  CHANGE THE CASE OF A STRING TO LOWER AND UPPER CASE
!********************************************************

 FUNCTION to_UpperCase (str)
 !==============================
 !Changes a string to upper case
 !==============================
  IMPLICIT NONE
  Character(*), Intent(In) :: str
  Character(LEN(str))      :: to_UpperCase
  Character(LEN(str))      :: string

  Integer :: ic, i

  Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

 !Capitalize each letter if it is lowecase
  string = str
  do i = 1, LEN_TRIM(str)
   ic = INDEX(low, str(i:i))
   if (ic > 0) string(i:i) = cap(ic:ic)
  end do

  to_UpperCase=trim(string)

  RETURN
 End Function to_UpperCase

 FUNCTION to_LowerCase (str)
 !==============================
 !Changes a string to lower case
 !==============================
  IMPLICIT NONE
  Character(*), Intent(In) :: str
  Character(LEN(str))      :: to_LowerCase
  Character(LEN(str))      :: string

  Integer :: ic, i

  Character(26), Parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  Character(26), Parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

 !Capitalize each letter if it is lowecase
  string = str
  do i = 1, LEN_TRIM(str)
   ic = INDEX(cap, str(i:i))
   if (ic > 0) string(i:i) = low(ic:ic)
  end do

  to_LowerCase=TRIM(string)

  RETURN
 END FUNCTION to_LowerCase

END MODULE SYSSOLVER

