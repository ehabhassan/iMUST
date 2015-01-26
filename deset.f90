MODULE DESET
 USE DOPERATORS

!***************************************************************
!                MODULE PROCEDURAL INTERFACES
!***************************************************************
 INTERFACE ODESET
   MODULE PROCEDURE c8ODESET1D,c8ODESET2D,c8ODESET3D,r8ODESET1D,r8ODESET2D,r8ODESET3D
 END INTERFACE ODESET


 CONTAINS


!****************************************************************
!****************** SYSTEM OF FIRST-ORDER PDEs ******************
!****************************************************************

 SUBROUTINE r8ODESET1D(ODESETPIN,fld,Lfld,NLfld,iXDIM,iYDIM,iZDIM)
  IMPLICIT NONE
  INTEGER(4),INTENT(IN)                     :: ODESETPIN
  INTEGER(4),INTENT(IN),OPTIONAL            :: iXDIM,iYDIM,iZDIM

  REAL(8),INTENT(INOUT),DIMENSION(:,:)      ::   fld
  REAL(8),INTENT(INOUT),DIMENSION(:,:,:)    ::  Lfld
  REAL(8),INTENT(INOUT),DIMENSION(:,:)      :: NLfld

  INTEGER(4)                                :: iX   ,iY   ,iZ   
  INTEGER(4)                                :: ifld1,ifld2

  INTEGER(4)                                :: fDIMS1
  INTEGER(4)                                :: nflds
  INTEGER(4)                                :: fldMATDIMS(2)

  INTEGER(4), SAVE                          :: ACCESSID=1
  INTEGER(4), SAVE                          :: ODESETID=1

  REAL(8),ALLOCATABLE,DIMENSION(:)           :: tDy
  REAL(8),ALLOCATABLE,DIMENSION(:)           :: Dy

  fldMATDIMS=SHAPE(fld)
  fDIMS1=fldMATDIMS(1)
  nflds =fldMATDIMS(2)

  IF(PRESENT(iXDIM).EQV..TRUE.) THEN
   iX=iXDIM
  ELSE
   iX=1
  END IF

  IF(PRESENT(iYDIM).EQV..TRUE.) THEN
   iY=iYDIM
  ELSE
   iY=1
  END IF

  IF(PRESENT(iZDIM).EQV..TRUE.) THEN
   iZ=iZDIM
  ELSE
   iZ=1
  END IF

  IF(ALLOCATED( Dy).EQV..FALSE.) ALLOCATE( Dy(fDIMS1))
  IF(ALLOCATED(tDy).EQV..FALSE.) ALLOCATE(tDy(fDIMS1))

  IF(ACCESSID.GT.1.AND.ODESETID.NE.ODESETPIN) ACCESSID=1

  IF(ACCESSID.EQ.1) THEN
  !ALLOCATE THE BOUNDARY CONDITIONS ALLOCATABLE ARRAYS
   CALL RESETBC(nflds)
  !ALLOCATE THE DIFFERENTIATION SCHEMES ARRAYS
   CALL RESETDSCHEMES(nflds)
  !ALLOCATE THE INTEGRATION SCHEMES ARRAYS
   CALL RESETISCHEMES(nflds)
   IF(ODESETPIN.EQ.1) THEN
   !SET THE DIFFERENTIATION SCHEMES ARRAYS
    PSMDSCHEME=.TRUE.
   !SET THE INTEGERATION SCHEMES ARRAYS
    PSMISCHEME=.TRUE.
    ACCESSID=ACCESSID+1
    ODESETID=ODESETPIN
   END IF
  END IF

   Lfld=0.0d0
  NLfld=0.0d0

  IF(ODESETPIN.EQ.1) THEN

  END IF

 RETURN
 END SUBROUTINE r8ODESET1D

 SUBROUTINE r8ODESET2D(ODESETPIN,fld,Lfld,NLfld,iXDIM,iYDIM,iZDIM)
  USE SYSINFO
  USE IONOBG
  IMPLICIT NONE
  INTEGER(4),INTENT(IN)                       :: ODESETPIN
  INTEGER(4),INTENT(IN),OPTIONAL              :: iXDIM,iYDIM,iZDIM

  REAL(8),INTENT(INOUT),DIMENSION(:,:,:)      ::   fld
  REAL(8),INTENT(INOUT),DIMENSION(:,:,:,:)    ::  Lfld
  REAL(8),INTENT(INOUT),DIMENSION(:,:,:)      :: NLfld

  INTEGER(4)                                  :: iX,   iY,   iZ
  INTEGER(4)                                  :: ifld1,ifld2

  INTEGER(4)                                  :: fDIMS1,fDIMS2
  INTEGER(4)                                  :: nflds
  INTEGER(4)                                  :: fldMATDIMS(3)

  INTEGER(4), SAVE                            :: ACCESSID=1
  INTEGER(4), SAVE                            :: ODESETID=1

  REAL(8),ALLOCATABLE,DIMENSION(:,:)           :: tDy
  REAL(8),ALLOCATABLE,DIMENSION(:,:)           :: Dy


  fldMATDIMS=SHAPE(fld)
  fDIMS1=fldMATDIMS(1)
  fDIMS2=fldMATDIMS(2)
  nflds =fldMATDIMS(3)

  IF(PRESENT(iXDIM).EQV..TRUE.) THEN
   iX=iXDIM
  ELSE
   iX=1
  END IF

  IF(PRESENT(iYDIM).EQV..TRUE.) THEN
   iY=iYDIM
  ELSE
   iY=1
  END IF

  IF(PRESENT(iZDIM).EQV..TRUE.) THEN
   iZ=iZDIM
  ELSE
   iZ=1
  END IF

  IF(ALLOCATED( Dy).EQV..FALSE.) ALLOCATE( Dy(fDIMS1,fDIMS2))
  IF(ALLOCATED(tDy).EQV..FALSE.) ALLOCATE(tDy(fDIMS1,fDIMS2))

  IF(ACCESSID.GT.1.AND.ODESETID.NE.ODESETPIN) ACCESSID=1

  IF(ACCESSID.EQ.1) THEN
  !ALLOCATE THE BOUNDARY CONDITIONS ALLOCATABLE ARRAYS
   CALL RESETBC(nflds)
  !ALLOCATE THE DIFFERENTIATION SCHEMES ARRAYS
   CALL RESETDSCHEMES(nflds)
  !ALLOCATE THE INTEGRATION SCHEMES ARRAYS
   CALL RESETISCHEMES(nflds)
   IF(ODESETPIN.EQ.1) THEN
   !SET REQUIRED TYPE OF THE BOUNDARY CONDITIONS ARRAYS
    PERIODICBC=.TRUE.
   !SET THE DIFFERENTIATION SCHEMES ARRAYS
    PSMDSCHEME=.TRUE.
   !FDMDSCHEME=.TRUE.
   !SET THE INTEGERATION SCHEMES ARRAYS
    PSMISCHEME=.TRUE.
   !FDMISCHEME=.TRUE.
    ACCESSID=ACCESSID+1
    ODESETID=ODESETPIN
    PRINT*, "ODESET OF THE TYPE-I AND TYPE-II INSTABILITIES (R-SPACE)"
   ELSE IF(ODESETPIN.EQ.2) THEN
   !SET REQUIRED TYPE OF THE BOUNDARY CONDITIONS ARRAYS
    PERIODICBC=.TRUE.
   !DIRICHLETBC=.TRUE.
   !SET REQUIRED VALUE OF THE BOUNDARY CONDITIONS ARRAYS
   !DIRICHLETVAL=0
   !SET THE DIFFERENTIATION SCHEMES ARRAYS
    FDMDSCHEME=.TRUE.
   !PSMDSCHEME=.TRUE.
   !CSMDSCHEME=.TRUE.
   !SET THE INTEGERATION SCHEMES ARRAYS
    FDMISCHEME=.TRUE.
   !PSMISCHEME=.TRUE.
   !CSMISCHEME=.TRUE.
    ACCESSID=ACCESSID+1
    ODESETID=ODESETPIN
    PRINT*, "RSPACESOLVER: ODESET FOR THE WAVE EQUATION"
   ELSE IF(ODESETPIN.EQ.3) THEN
   !SET REQUIRED TYPE OF THE BOUNDARY CONDITIONS ARRAYS
   !PERIODICBC=.TRUE.
    DIRICHLETBC=.TRUE.
   !NEUMANNBC=.TRUE.
   !SET REQUIRED VALUE OF THE BOUNDARY CONDITIONS ARRAYS
    DIRICHLETVAL=0.0
   !NEUMANNVAL=0.0
   !SET THE DIFFERENTIATION SCHEMES ARRAYS
    CSMDSCHEME=.TRUE.
   !PSMDSCHEME=.TRUE.
   !FDMDSCHEME=.TRUE.
   !SET THE INTEGERATION SCHEMES ARRAYS
   !CSMISCHEME=.TRUE.
   !FDMISCHEME=.TRUE.
    ACCESSID=ACCESSID+1
    ODESETID=ODESETPIN
    PRINT*, "RSPACESOLVER: ODESET FOR THE HEAT EQUATION"
   END IF
  END IF

   Lfld=0.0d0
  NLfld=0.0d0

  IF(ODESETPIN.EQ.1) THEN
   !"RSPACESOLVER: ODESET OF THE TYPE-I AND TYPE-II INSTABILITIES"
!******* START OF NONLINEAR TERMS ********
!************** FIRST FIELD **************
    CALL GRADfGRADg(fld(:,:,1),fld(:,:,3),Y3D(1,:,:),Z3D(1,:,:),Dy)
    NLfld(:,:,1)=NLfld(:,:,1)+Dy

    CALL LPLCf(fld(:,:,3),Y3D(1,:,:),Z3D(1,:,:),tDy,1)
    NLfld(:,:,1)=NLfld(:,:,1)+fld(:,:,1)*tDy

!************** SECOND FIELD **************
    CALL GRADfGRADg(fld(:,:,1),fld(:,:,1),Y3D(1,:,:),Z3D(1,:,:),Dy)
    NLfld(:,:,2)=NLfld(:,:,2)+NUen(iX,iY,iZ)*Dy

    CALL GRADfGRADg(fld(:,:,1),fld(:,:,2),Y3D(1,:,:),Z3D(1,:,:),Dy)
    NLfld(:,:,2)=NLfld(:,:,2)-NUen(iX,iY,iZ)*Dy

    CALL GRADfGRADg(fld(:,:,1),fld(:,:,3),Y3D(1,:,:),Z3D(1,:,:),Dy)
    NLfld(:,:,2)=NLfld(:,:,2)-Dy/RHOe(iX,iY,iZ)**2

    CALL PBKfg(fld(:,:,2),fld(:,:,1),Y3D(1,:,:),Z3D(1,:,:),Dy)
    NLfld(:,:,2)=NLfld(:,:,2)-OMEGAce(iX,iY,iZ)*Dy

    CALL LPLCf(fld(:,:,2),Y3D(1,:,:),Z3D(1,:,:),tDy,1)
    CALL PBKfg(fld(:,:,2),tDy,Y3D(1,:,:),Z3D(1,:,:),Dy)
    NLfld(:,:,2)=NLfld(:,:,2)-(OMEGAce(iX,iY,iZ)*RHOe(iX,iY,iZ)**2)*Dy

    CALL LPLCf(fld(:,:,2),Y3D(1,:,:),Z3D(1,:,:),tDy,1)
    CALL PBKfg(fld(:,:,1),tDy,Y3D(1,:,:),Z3D(1,:,:),Dy)
    NLfld(:,:,2)=NLfld(:,:,2)+(OMEGAce(iX,iY,iZ)*RHOe(iX,iY,iZ)**2)*Dy

    CALL LPLCf(fld(:,:,1),Y3D(1,:,:),Z3D(1,:,:),tDy,1)
    CALL PBKfg(fld(:,:,2),tDy,Y3D(1,:,:),Z3D(1,:,:),Dy)
    NLfld(:,:,2)=NLfld(:,:,2)+(OMEGAce(iX,iY,iZ)*RHOe(iX,iY,iZ)**2)*Dy

!************** THIRD FIELD **************
    CALL GRADfGRADg(fld(:,:,3),fld(:,:,3),Y3D(1,:,:),Z3D(1,:,:),tDy)
    CALL LPLCf(tDy,Y3D(1,:,:),Z3D(1,:,:),Dy,1)
    NLfld(:,:,3)=NLfld(:,:,3)+0.5d0*Dy

   DO ifld1=2,nflds
    CALL intLPLCf(NLfld(:,:,ifld1),Y3D(1,:,:),Z3D(1,:,:),Dy)
    NLfld(:,:,ifld1)=Dy
   END DO

!********* START OF LINEAR TERMS *********
!************** FIRST FIELD **************
    CALL LPLCf(fld(:,:,3),Y3D(1,:,:),Z3D(1,:,:),Dy,1)
    Lfld(:,:,1,3)=Lfld(:,:,1,3)+Dy

!************** SECOND FIELD **************
    CALL LPLCf(fld(:,:,1),Y3D(1,:,:),Z3D(1,:,:),Dy,1)
    Lfld(:,:,2,1)=Lfld(:,:,2,1)+NUen(iX,iY,iZ)*Dy

    CALL Df(fld(:,:,1),Z3D(1,:,:),Dy,1)
    Lfld(:,:,2,1)=Lfld(:,:,2,1)+2.0*(NUen(iX,iY,iZ)/Lne(iX,iY,iZ))*Dy

    CALL Df(fld(:,:,1),Z3D(1,:,:),Dy,1)
    Lfld(:,:,2,1)=Lfld(:,:,2,1)+((Vexb(iX,iY,iZ)/RHOe(iX,iY,iZ)**2)*(NUen(iX,iY,iZ)/OMEGAce(iX,iY,iZ)))*Dy

    CALL Df(fld(:,:,1),Y3D(1,:,:),Dy,1)
    Lfld(:,:,2,1)=Lfld(:,:,2,1)-(Vexb(iX,iY,iZ)/RHOe(iX,iY,iZ)**2)*Dy

    CALL LPLCf(fld(:,:,1),Y3D(1,:,:),Z3D(1,:,:),tDy,1)
    CALL Df(tDy,Y3D(1,:,:),Dy,1)
    Lfld(:,:,2,1)=Lfld(:,:,2,1)+Vexb(iX,iY,iZ)*Dy

    CALL LPLCf(fld(:,:,2),Y3D(1,:,:),Z3D(1,:,:),Dy,1)
    Lfld(:,:,2,2)=Lfld(:,:,2,2)-NUen(iX,iY,iZ)*Dy

    CALL Df(fld(:,:,2),Z3D(1,:,:),Dy,1)
    Lfld(:,:,2,2)=Lfld(:,:,2,2)-(NUen(iX,iY,iZ)/Lne(iX,iY,iZ))*Dy

    CALL Df(fld(:,:,2),Y3D(1,:,:),Dy,1)
    Lfld(:,:,2,2)=Lfld(:,:,2,2)-(OMEGAce(iX,iY,iZ)/Lne(iX,iY,iZ))*Dy

    CALL LPLCf(fld(:,:,2),Y3D(1,:,:),Z3D(1,:,:),tDy,1)
    CALL Df(tDy,Y3D(1,:,:),Dy,1)
    Lfld(:,:,2,2)=Lfld(:,:,2,2)-Vexb(iX,iY,iZ)*Dy

    CALL LPLCf(fld(:,:,2),Y3D(1,:,:),Z3D(1,:,:),tDy,1)
    CALL Df(tDy,Y3D(1,:,:),Dy,1)
    Lfld(:,:,2,2)=Lfld(:,:,2,2)-(OMEGAce(iX,iY,iZ)*RHOe(iX,iY,iZ)**2/Lne(iX,iY,iZ))*Dy

    CALL LPLCf(fld(:,:,3),Y3D(1,:,:),Z3D(1,:,:),Dy,1)
    Lfld(:,:,2,3)=Lfld(:,:,2,3)-Dy/RHOe(iX,iY,iZ)**2

    CALL Df(fld(:,:,3),Z3D(1,:,:),Dy,1)
    Lfld(:,:,2,3)=Lfld(:,:,2,3)-Dy/(Lne(iX,iY,iZ)*RHOe(iX,iY,iZ)**2)

!************** THIRD FIELD **************
    CALL LPLCf(fld(:,:,1),Y3D(1,:,:),Z3D(1,:,:),Dy,1)
    Lfld(:,:,3,1)=Lfld(:,:,3,1)+Vthi(iX,iY,iZ)**2*Dy

    CALL LPLCf(fld(:,:,2),Y3D(1,:,:),Z3D(1,:,:),Dy,1)
    Lfld(:,:,3,2)=Lfld(:,:,3,2)+Vthi(iX,iY,iZ)**2*Dy

    CALL LPLCf(fld(:,:,3),Y3D(1,:,:),Z3D(1,:,:),Dy,1)
    Lfld(:,:,3,3)=Lfld(:,:,3,3)-NUin(iX,iY,iZ)*Dy

    CALL LPLCf(fld(:,:,3),Y3D(1,:,:),Z3D(1,:,:),tDy)
    CALL LPLCf(tDy,Y3D(1,:,:),Z3D(1,:,:),Dy)
    Lfld(:,:,3,3)=Lfld(:,:,3,3)+((4.0d0/3.0d0)*(Vthi(iX,iY,iZ)**2/NUin(iX,iY,iZ)))*Dy

   DO ifld1=2,nflds
    DO ifld2=1,nflds
     CALL intLPLCf(Lfld(:,:,ifld1,ifld2),Y3D(1,:,:),Z3D(1,:,:),Dy)
     Lfld(:,:,ifld1,ifld2)=Dy
    END DO
   END DO

  ELSE IF(ODESETPIN.EQ.2) THEN
   CALL FORCEBC(fld(:,:,1),1)
   CALL FORCEBC(fld(:,:,2),2)

   CALL LPLCf(fld(:,:,1),Y3D(1,:,:),Z3D(1,:,:),Dy)
   Lfld(:,:,1,1)=fld(:,:,2)
   Lfld(:,:,2,2)=Dy

  ELSE IF(ODESETPIN.EQ.3) THEN
   CALL FORCEBC(fld(:,:,1),1)

   CALL LPLCf(fld(:,:,1),Y3D(1,:,:),Z3D(1,:,:),Dy)
   Lfld(:,:,1,1)=Dy

  END IF

 RETURN
 END SUBROUTINE r8ODESET2D

 SUBROUTINE r8ODESET3D(ODESETPIN,fld,Lfld,NLfld,iXDIM,iYDIM,iZDIM)
  IMPLICIT NONE
  INTEGER(4),INTENT(IN)                         :: ODESETPIN
  INTEGER(4),INTENT(IN),OPTIONAL                :: iXDIM,iYDIM,iZDIM

  REAL(8),INTENT(INOUT),DIMENSION(:,:,:,:)      ::   fld
  REAL(8),INTENT(INOUT),DIMENSION(:,:,:,:,:)    ::  Lfld
  REAL(8),INTENT(INOUT),DIMENSION(:,:,:,:)      :: NLfld

  INTEGER(4)                                    :: iX   ,iY   ,iZ   
  INTEGER(4)                                    :: ifld1,ifld2

  INTEGER(4)                                    :: fDIMS1,fDIMS2,fDIMS3
  INTEGER(4)                                    :: nflds
  INTEGER(4)                                    :: fldMATDIMS(4)

  INTEGER(4), SAVE                              :: ACCESSID=1
  INTEGER(4), SAVE                              :: ODESETID=1

  REAL(8),ALLOCATABLE,DIMENSION(:,:,:)          :: tDy
  REAL(8),ALLOCATABLE,DIMENSION(:,:,:)          :: Dy


  fldMATDIMS=SHAPE(fld)
  fDIMS1=fldMATDIMS(1)
  fDIMS2=fldMATDIMS(2)
  fDIMS3=fldMATDIMS(3)
  nflds =fldMATDIMS(4)

  IF(PRESENT(iXDIM).EQV..TRUE.) THEN
   iX=iXDIM
  ELSE
   iX=1
  END IF

  IF(PRESENT(iYDIM).EQV..TRUE.) THEN
   iY=iYDIM
  ELSE
   iY=1
  END IF

  IF(PRESENT(iZDIM).EQV..TRUE.) THEN
   iZ=iZDIM
  ELSE
   iZ=1
  END IF

  IF(ALLOCATED( Dy).EQV..FALSE.) ALLOCATE( Dy(fDIMS1,fDIMS2,fDIMS3))
  IF(ALLOCATED(tDy).EQV..FALSE.) ALLOCATE(tDy(fDIMS1,fDIMS2,fDIMS3))

  IF(ACCESSID.GT.1.AND.ODESETID.NE.ODESETPIN) ACCESSID=1

  IF(ACCESSID.EQ.1) THEN
  !ALLOCATE THE BOUNDARY CONDITIONS ALLOCATABLE ARRAYS
   CALL RESETBC(nflds)
  !ALLOCATE THE DIFFERENTIATION SCHEMES ARRAYS
   CALL RESETDSCHEMES(nflds)
  !ALLOCATE THE INTEGRATION SCHEMES ARRAYS
   CALL RESETISCHEMES(nflds)
   IF(ODESETPIN.EQ.1) THEN
   !SET THE DIFFERENTIATION SCHEMES ARRAYS
    PSMDSCHEME=.TRUE.
   !SET THE INTEGERATION SCHEMES ARRAYS
    PSMISCHEME=.TRUE.
    ACCESSID=ACCESSID+1
    ODESETID=ODESETPIN
   END IF
  END IF

   Lfld=0.0d0
  NLfld=0.0d0

  IF(ODESETPIN.EQ.1) THEN

  END IF

  RETURN
 END SUBROUTINE r8ODESET3D


 SUBROUTINE c8ODESET1D(ODESETPIN,fld,Lfld,NLfld,iXDIM,iYDIM,iZDIM)
  IMPLICIT NONE
  INTEGER(4),INTENT(IN)                        :: ODESETPIN
  INTEGER(4),INTENT(IN),OPTIONAL               :: iXDIM,iYDIM,iZDIM

  INTEGER(4)                                   :: iX   ,iY   ,iZ   
  INTEGER(4)                                   :: ifld1,ifld2

  INTEGER(4)                                   :: fDIMS1
  INTEGER(4)                                   :: nflds
  INTEGER(4)                                   :: fldMATDIMS(2)

  INTEGER(4), SAVE                             :: ACCESSID=1
  INTEGER(4), SAVE                             :: ODESETID=1

  COMPLEX(8),INTENT(INOUT),DIMENSION(:,:)      ::   fld
  COMPLEX(8),INTENT(INOUT),DIMENSION(:,:,:)    ::  Lfld
  COMPLEX(8),INTENT(INOUT),DIMENSION(:,:)      :: NLfld

  COMPLEX(8),ALLOCATABLE,DIMENSION(:)           :: tDy
  COMPLEX(8),ALLOCATABLE,DIMENSION(:)           :: Dy


  fldMATDIMS=SHAPE(fld)
  fDIMS1=fldMATDIMS(1)
  nflds =fldMATDIMS(2)

  IF(PRESENT(iXDIM).EQV..TRUE.) THEN
   iX=iXDIM
  ELSE
   iX=1
  END IF

  IF(PRESENT(iYDIM).EQV..TRUE.) THEN
   iY=iYDIM
  ELSE
   iY=1
  END IF

  IF(PRESENT(iZDIM).EQV..TRUE.) THEN
   iZ=iZDIM
  ELSE
   iZ=1
  END IF

  IF(ALLOCATED( Dy).EQV..FALSE.) ALLOCATE( Dy(fDIMS1))
  IF(ALLOCATED(tDy).EQV..FALSE.) ALLOCATE(tDy(fDIMS1))

  IF(ACCESSID.GT.1.AND.ODESETID.NE.ODESETPIN) ACCESSID=1

  IF(ACCESSID.EQ.1) THEN
  !ALLOCATE THE BOUNDARY CONDITIONS ALLOCATABLE ARRAYS
   CALL RESETBC(nflds)
  !ALLOCATE THE DIFFERENTIATION SCHEMES ARRAYS
   CALL RESETDSCHEMES(nflds)
  !ALLOCATE THE INTEGRATION SCHEMES ARRAYS
   CALL RESETISCHEMES(nflds)

   IF(ODESETPIN.EQ.1) THEN
   !SET THE DIFFERENTIATION SCHEMES ARRAYS
    PSMDSCHEME=.TRUE.
   !SET THE INTEGERATION SCHEMES ARRAYS
    PSMISCHEME=.TRUE.
    ACCESSID=ACCESSID+1
    ODESETID=ODESETPIN
   END IF
  END IF

   Lfld=CMPLX(0.0d0,0.0d0)
  NLfld=CMPLX(0.0d0,0.0d0)

  IF(ODESETPIN.EQ.1) THEN

  END IF

 RETURN
 END SUBROUTINE c8ODESET1D

 SUBROUTINE c8ODESET2D(ODESETPIN,fld,Lfld,NLfld,iXDIM,iYDIM,iZDIM)
  USE IONOBG
  IMPLICIT NONE
  INTEGER(4),INTENT(IN)                          :: ODESETPIN
  INTEGER(4),INTENT(IN),OPTIONAL                 :: iXDIM,iYDIM,iZDIM

  COMPLEX(8),INTENT(INOUT),DIMENSION(:,:,:)      ::   fld
  COMPLEX(8),INTENT(INOUT),DIMENSION(:,:,:,:)    ::  Lfld
  COMPLEX(8),INTENT(INOUT),DIMENSION(:,:,:)      :: NLfld

  INTEGER(4)                                     :: iX,   iY,   iZ
  INTEGER(4)                                     :: ifld1,ifld2

  INTEGER(4)                                     :: fDIMS1,fDIMS2
  INTEGER(4)                                     :: nflds
  INTEGER(4)                                     :: fldMATDIMS(3)

  INTEGER(4), SAVE                               :: ACCESSID=1
  INTEGER(4), SAVE                               :: ODESETID=1

  COMPLEX(8),ALLOCATABLE,DIMENSION(:,:)           :: tDy
  COMPLEX(8),ALLOCATABLE,DIMENSION(:,:)           :: Dy


  fldMATDIMS=SHAPE(fld)
  fDIMS1=fldMATDIMS(1)
  fDIMS2=fldMATDIMS(2)
  nflds =fldMATDIMS(3)

  IF(PRESENT(iXDIM).EQV..TRUE.) THEN
   iX=iXDIM
  ELSE
   iX=1
  END IF

  IF(PRESENT(iYDIM).EQV..TRUE.) THEN
   iY=iYDIM
  ELSE
   iY=1
  END IF

  IF(PRESENT(iZDIM).EQV..TRUE.) THEN
   iZ=iZDIM
  ELSE
   iZ=1
  END IF

  IF(ALLOCATED( Dy).EQV..FALSE.) ALLOCATE( Dy(fDIMS1,fDIMS2))
  IF(ALLOCATED(tDy).EQV..FALSE.) ALLOCATE(tDy(fDIMS1,fDIMS2))

  IF(ACCESSID.GT.1.AND.ODESETID.NE.ODESETPIN) ACCESSID=1

  IF(ACCESSID.EQ.1) THEN
  !ALLOCATE THE BOUNDARY CONDITIONS ALLOCATABLE ARRAYS
   CALL RESETBC(nflds)
  !ALLOCATE THE DIFFERENTIATION SCHEMES ARRAYS
   CALL RESETDSCHEMES(nflds)
  !ALLOCATE THE INTEGRATION SCHEMES ARRAYS
   CALL RESETISCHEMES(nflds)

   IF(ODESETPIN.EQ.1) THEN
   !SET THE DIFFERENTIATION SCHEMES ARRAYS
    PSMDSCHEME=.TRUE.
   !SET THE INTEGERATION SCHEMES ARRAYS
    PSMISCHEME=.TRUE.
    ACCESSID=ACCESSID+1
    ODESETID=ODESETPIN
    PRINT*, "ODESET OF THE TYPE-I AND TYPE-II INSTABILITIES (K-SPACE)"
   END IF
  END IF

  IF(ALLOCATED(iKy2D).EQV..FALSE.) PRINT*,"iKy2D NOT ALLOCATED"
  IF(ALLOCATED(iKz2D).EQV..FALSE.) PRINT*,"iKz2D NOT ALLOCATED"
   Lfld=CMPLX(0.0d0,0.0d0)
  NLfld=CMPLX(0.0d0,0.0d0)

  IF(ODESETPIN.EQ.1) THEN
!************** FIRST FIELD **************
    CALL GRADfGRADg(fld(:,:,1),fld(:,:,3),iKy2D,iKz2D,Dy)
    NLfld(:,:,1)=NLfld(:,:,1)+Dy

    CALL LPLCf(fld(:,:,3),iKy2D,iKz2D,tDy,1)
    CALL fg(fld(:,:,1),tDy,Dy)
    NLfld(:,:,1)=NLfld(:,:,1)+Dy

!************** SECOND FIELD **************
    CALL GRADfGRADg(fld(:,:,1),fld(:,:,1),iKy2D,iKz2D,Dy)
    NLfld(:,:,2)=NLfld(:,:,2)+NUen(iX,iY,iZ)*Dy

    CALL GRADfGRADg(fld(:,:,1),fld(:,:,2),iKy2D,iKz2D,Dy)
    NLfld(:,:,2)=NLfld(:,:,2)-NUen(iX,iY,iZ)*Dy

    CALL GRADfGRADg(fld(:,:,1),fld(:,:,3),iKy2D,iKz2D,Dy)
    NLfld(:,:,2)=NLfld(:,:,2)-Dy/RHOe(iX,iY,iZ)**2

    CALL PBKfg(fld(:,:,2),fld(:,:,1),iKy2D,iKz2D,Dy)
    NLfld(:,:,2)=NLfld(:,:,2)-OMEGAce(iX,iY,iZ)*Dy

    CALL LPLCf(fld(:,:,2),iKy2D,iKz2D,tDy,1)
    CALL PBKfg(fld(:,:,2),tDy,iKy2D,iKz2D,Dy)
    NLfld(:,:,2)=NLfld(:,:,2)-(OMEGAce(iX,iY,iZ)*RHOe(iX,iY,iZ)**2)*Dy

    CALL LPLCf(fld(:,:,2),iKy2D,iKz2D,tDy,1)
    CALL PBKfg(fld(:,:,1),tDy,iKy2D,iKz2D,Dy)
    NLfld(:,:,2)=NLfld(:,:,2)+(OMEGAce(iX,iY,iZ)*RHOe(iX,iY,iZ)**2)*Dy

!   CALL LPLCf(fld(:,:,1),iKy2D,iKz2D,tDy,1)
!   CALL PBKfg(fld(:,:,2),tDy,iKy2D,iKz2D,Dy)
!   NLfld(:,:,2)=NLfld(:,:,2)+(OMEGAce(iX,iY,iZ)*RHOe(iX,iY,iZ)**2)*Dy

!************** THIRD FIELD **************
    CALL GRADfGRADg(fld(:,:,3),fld(:,:,3),iKy2D,iKz2D,tDy)
    CALL LPLCf(tDy,iKy2D,iKz2D,Dy,1)
    NLfld(:,:,3)=NLfld(:,:,3)+0.5d0*Dy

    DO ifld1=2,nflds
     CALL intLPLCf(NLfld(:,:,ifld1),iKy2D,iKz2D,Dy,1)
     NLfld(:,:,ifld1)=Dy
    END DO

!************** FIRST FIELD **************
    CALL LPLCf(fld(:,:,3),iKy2D,iKz2D,Dy,1)
    Lfld(:,:,1,3)=Lfld(:,:,1,3)+Dy

!************** SECOND FIELD **************
    CALL LPLCf(fld(:,:,1),iKy2D,iKz2D,Dy,1)
    Lfld(:,:,2,1)=Lfld(:,:,2,1)+NUen(iX,iY,iZ)*Dy

    CALL Df(fld(:,:,1),iKz2D,Dy,1)
    Lfld(:,:,2,1)=Lfld(:,:,2,1)+2.0*(NUen(iX,iY,iZ)/Lne(iX,iY,iZ))*Dy

    CALL Df(fld(:,:,1),iKz2D,Dy,1)
    Lfld(:,:,2,1)=Lfld(:,:,2,1)+((Vexb(iX,iY,iZ)/RHOe(iX,iY,iZ)**2)*(NUen(iX,iY,iZ)/OMEGAce(iX,iY,iZ)))*Dy

    CALL Df(fld(:,:,1),iKy2D,Dy,1)
    Lfld(:,:,2,1)=Lfld(:,:,2,1)-(Vexb(iX,iY,iZ)/RHOe(iX,iY,iZ)**2)*Dy

!   CALL LPLCf(fld(:,:,1),iKy2D,iKz2D,tDy,1)
!   CALL Df(tDy,iKy2D,Dy,1)
!   Lfld(:,:,2,1)=Lfld(:,:,2,1)+Vexb(iX,iY,iZ)*Dy

    CALL LPLCf(fld(:,:,2),iKy2D,iKz2D,Dy,1)
    Lfld(:,:,2,2)=Lfld(:,:,2,2)-NUen(iX,iY,iZ)*Dy

    CALL Df(fld(:,:,2),iKz2D,Dy,1)
    Lfld(:,:,2,2)=Lfld(:,:,2,2)-(NUen(iX,iY,iZ)/Lne(iX,iY,iZ))*Dy

    CALL Df(fld(:,:,2),iKy2D,Dy,1)
    Lfld(:,:,2,2)=Lfld(:,:,2,2)-(OMEGAce(iX,iY,iZ)/Lne(iX,iY,iZ))*Dy

    CALL LPLCf(fld(:,:,2),iKy2D,iKz2D,tDy,1)
    CALL Df(tDy,iKy2D,Dy,1)
    Lfld(:,:,2,2)=Lfld(:,:,2,2)-Vexb(iX,iY,iZ)*Dy

    CALL LPLCf(fld(:,:,2),iKy2D,iKz2D,tDy,1)
    CALL Df(tDy,iKy2D,Dy,1)
    Lfld(:,:,2,2)=Lfld(:,:,2,2)-(OMEGAce(iX,iY,iZ)*RHOe(iX,iY,iZ)**2/Lne(iX,iY,iZ))*Dy

    CALL LPLCf(fld(:,:,3),iKy2D,iKz2D,Dy,1)
    Lfld(:,:,2,3)=Lfld(:,:,2,3)-Dy/RHOe(iX,iY,iZ)**2

    CALL Df(fld(:,:,3),iKz2D,Dy,1)
    Lfld(:,:,2,3)=Lfld(:,:,2,3)-Dy/(Lne(iX,iY,iZ)*RHOe(iX,iY,iZ)**2)

!************** THIRD FIELD **************
    CALL LPLCf(fld(:,:,1),iKy2D,iKz2D,Dy,1)
    Lfld(:,:,3,1)=Lfld(:,:,3,1)+Vthi(iX,iY,iZ)**2*Dy

    CALL LPLCf(fld(:,:,2),iKy2D,iKz2D,Dy,1)
    Lfld(:,:,3,2)=Lfld(:,:,3,2)+Vthi(iX,iY,iZ)**2*Dy

    CALL LPLCf(fld(:,:,3),iKy2D,iKz2D,Dy,1)
    Lfld(:,:,3,3)=Lfld(:,:,3,3)-NUin(iX,iY,iZ)*Dy

    CALL LPLCf(fld(:,:,3),iKy2D,iKz2D,Dy,2)
    Lfld(:,:,3,3)=Lfld(:,:,3,3)+((4.0d0/3.0d0)*(Vthi(iX,iY,iZ)**2/NUin(iX,iY,iZ)))*Dy

    DO ifld1=2,nflds
     DO ifld2=1,nflds
      CALL intLPLCf(Lfld(:,:,ifld1,ifld2),iKy2D,iKz2D,Dy,1)
      Lfld(:,:,ifld1,ifld2)=Dy
     END DO
    END DO
  END IF

 RETURN
 END SUBROUTINE c8ODESET2D

 SUBROUTINE c8ODESET3D(ODESETPIN,fld,Lfld,NLfld,iXDIM,iYDIM,iZDIM)
  IMPLICIT NONE
  INTEGER(4),INTENT(IN)                            :: ODESETPIN
  INTEGER(4),INTENT(IN),OPTIONAL                   :: iXDIM,iYDIM,iZDIM

  COMPLEX(8),INTENT(INOUT),DIMENSION(:,:,:,:)      ::   fld
  COMPLEX(8),INTENT(INOUT),DIMENSION(:,:,:,:,:)    ::  Lfld
  COMPLEX(8),INTENT(INOUT),DIMENSION(:,:,:,:)      :: NLfld

  INTEGER(4)                                       :: iX   ,iY   ,iZ   
  INTEGER(4)                                       :: ifld1,ifld2

  INTEGER(4)                                       :: fDIMS1,fDIMS2,fDIMS3
  INTEGER(4)                                       :: nflds
  INTEGER(4)                                       :: fldMATDIMS(4)

  INTEGER(4), SAVE                                 :: ACCESSID=1
  INTEGER(4), SAVE                                 :: ODESETID=1

  COMPLEX(8),ALLOCATABLE,DIMENSION(:,:,:)           :: tDy
  COMPLEX(8),ALLOCATABLE,DIMENSION(:,:,:)           :: Dy


  fldMATDIMS=SHAPE(fld)
  fDIMS1=fldMATDIMS(1)
  fDIMS2=fldMATDIMS(2)
  fDIMS3=fldMATDIMS(3)
  nflds =fldMATDIMS(4)

  IF(PRESENT(iXDIM).EQV..TRUE.) THEN
   iX=iXDIM
  ELSE
   iX=1
  END IF

  IF(PRESENT(iYDIM).EQV..TRUE.) THEN
   iY=iYDIM
  ELSE
   iY=1
  END IF

  IF(PRESENT(iZDIM).EQV..TRUE.) THEN
   iZ=iZDIM
  ELSE
   iZ=1
  END IF

  IF(ALLOCATED( Dy).EQV..FALSE.) ALLOCATE( Dy(fDIMS1,fDIMS2,fDIMS3))
  IF(ALLOCATED(tDy).EQV..FALSE.) ALLOCATE(tDy(fDIMS1,fDIMS2,fDIMS3))

  IF(ACCESSID.GT.1.AND.ODESETID.NE.ODESETPIN) ACCESSID=1

  IF(ACCESSID.EQ.1) THEN
  !ALLOCATE THE BOUNDARY CONDITIONS ALLOCATABLE ARRAYS
   CALL RESETBC(nflds)
  !ALLOCATE THE DIFFERENTIATION SCHEMES ARRAYS
   CALL RESETDSCHEMES(nflds)
  !ALLOCATE THE INTEGRATION SCHEMES ARRAYS
   CALL RESETISCHEMES(nflds)
   IF(ODESETPIN.EQ.1) THEN
   !SET THE DIFFERENTIATION SCHEMES ARRAYS
    PSMDSCHEME=.TRUE.
   !SET THE INTEGERATION SCHEMES ARRAYS
    PSMISCHEME=.TRUE.
    ACCESSID=ACCESSID+1
    ODESETID=ODESETPIN
   END IF
  END IF

   Lfld=CMPLX(0.0d0,0.0d0)
  NLfld=CMPLX(0.0d0,0.0d0)

  IF(ODESETPIN.EQ.1) THEN

  END IF

 RETURN
 END SUBROUTINE c8ODESET3D


END MODULE DESET

