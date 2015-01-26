PROGRAM BRIDGE
 USE iTEST
 USE IONOBG
 USE SYSINFO
 USE GEOAXES
 USE SYSGRIDS
 USE MATHSUBS
 USE SYSSOLVER
 USE RDSVNETCDF
 USE DOPERATORS
 USE FORTGRAPHER
 USE SYSGRAPHICS
 IMPLICIT NONE
 INTEGER(4)            :: T_START,T_END
 INTEGER(4)            :: CLOCK_RATE,CLOCK_MAX
 INTEGER(4)            :: iTIMER,iCOUNTER
 INTEGER(4)            :: NETCDFiREC=1
 INTEGER(4)            :: fldID
 INTEGER(4)            :: MATDIMS4D(4)
 INTEGER(4)            :: MATDIMS1D(1)
 INTEGER(4)            :: nDE ! NUMBER OF DIFFERENTIAL EQUATIONS IN THE SYSTEM
 INTEGER(4)            :: ODESETID
 INTEGER(4),PARAMETER  :: ntstp=100000   
 COMPLEX(8),ALLOCATABLE:: kfld(:,:,:,:)
 REAL(8),   ALLOCATABLE:: rfld(:,:,:,:)
 REAL(8)               :: r8PI=4.0d0*ATAN(1.0d0)
 REAL(4)               :: t_crnt,t_nxt,t_strt=0.0
 REAL(4)               :: tstp=1.0e-5 
 REAL(4)               :: RECTIME
 REAL(4),   ALLOCATABLE:: T1D(:)

 INTEGER(4)            :: VARDIM(4)

 INTEGER(4)            :: GFC(10)
 CHARACTER(100)        :: GA(33)=""

 INTEGER(4)            :: NOWTIME(7)
 CHARACTER(20)         :: NETCDFFNAME="iBRIDGE.nc"
 CHARACTER(25)         :: FLDSNAME(3)=(/CHARACTER(25) :: "Density","ElectricPotential","IonVelocityPotential"/)
 CHARACTER(10)         :: FLDSUNIT(3)=(/CHARACTER(10) :: "meter^-3","volts","meter^2/s"/)
!CHARACTER(10)         :: DIMSNAME(2)=(/CHARACTER(10) :: "X","T"/)
!CHARACTER(10)         :: DIMSUNIT(2)=(/CHARACTER(10) :: "meter","second"/)
!CHARACTER(10)         :: DIMSNAME(3)=(/CHARACTER(10) :: "X","Y","T"/)
!CHARACTER(10)         :: DIMSUNIT(3)=(/CHARACTER(10) :: "meter","meter","second"/)
 CHARACTER(10)         :: DIMSNAME(4)=(/CHARACTER(10) :: "X","Y","Z","T"/)
 CHARACTER(10)         :: DIMSUNIT(4)=(/CHARACTER(10) :: "meter","meter","meter","second"/)

 REAL(4)               :: RGEOG,THETAGEOG,PHIGEOG
 REAL(4)               :: RGEOM,THETAGEOM,PHIGEOM
 REAL(4)               :: RGEOC,THETAGEOC,PHIGEOC
 REAL(4)               :: THETAN,PHIN
 REAL(4)               :: IGRFYEAR,Bo
 REAL(4)               :: H,Hx,Hy,Hz,Hxy,Hdip,Hvar
 REAL(4)               :: GEOGXYZ(3),GEOCXYZ(3)


!CALL iDIAGNOSE
!CALL NETCDF_VAR_NAMES('OPENMPtst.nc')

 NETCDFFNAME="testNETCDFmod.nc"

 ALLOCATE(X1D(10));X1D=(/ (iCOUNTER, iCOUNTER = 1,   10,   1) /)
 ALLOCATE(Y1D(10));Y1D=(/ (iCOUNTER, iCOUNTER = 10,  100,  10) /)
 ALLOCATE(Z1D(10));Z1D=(/ (iCOUNTER, iCOUNTER = 100, 1000, 100) /)
 CALL SAVEDIMS(TRIM(NETCDFFNAME),X1D,Y1D,Z1D,DIMSNAME,DIMSUNIT)
 CALL DEFNEWFLDS(TRIM(NETCDFFNAME),3,FLDSNAME,FLDSUNIT)

 ALLOCATE(rfld(10,10,10,6))
 rfld(:,:,:,1)=1.0d0
 rfld(:,:,:,2)=2.0d0
 rfld(:,:,:,3)=3.0d0

 NETCDFiREC=1
 RECTIME=1.0
 CALL SAVEFLDS(TRIM(NETCDFFNAME),NETCDFiREC,RECTIME,rfld(:,:,:,1),DIMSNAME(4),FLDSNAME(1))
 CALL SAVEFLDS(TRIM(NETCDFFNAME),NETCDFiREC,RECTIME,rfld(:,:,:,2),DIMSNAME(4),FLDSNAME(2))
 CALL SAVEFLDS(TRIM(NETCDFFNAME),NETCDFiREC,RECTIME,rfld(:,:,:,3),DIMSNAME(4),FLDSNAME(3))
  
 rfld(:,:,:,1)=4.0d0
 rfld(:,:,:,2)=5.0d0
 rfld(:,:,:,3)=6.0d0

 NETCDFiREC=2
 RECTIME=2.0
 CALL SAVEFLDS(TRIM(NETCDFFNAME),NETCDFiREC,RECTIME,rfld(:,:,:,1),DIMSNAME(4),FLDSNAME(1))
 CALL SAVEFLDS(TRIM(NETCDFFNAME),NETCDFiREC,RECTIME,rfld(:,:,:,2),DIMSNAME(4),FLDSNAME(2))
 CALL SAVEFLDS(TRIM(NETCDFFNAME),NETCDFiREC,RECTIME,rfld(:,:,:,3),DIMSNAME(4),FLDSNAME(3))
  
 CALL NETCDF_VAR_DIMS(TRIM(NETCDFFNAME),DIMSNAME(4),VARDIM)
 ALLOCATE(T1D(VARDIM(1)))

 CALL READDIMS(TRIM(NETCDFFNAME),X1D,Y1D,Z1D,T1D,DIMSNAME,DIMSUNIT)
 NETCDFiREC=1
 CALL READFLDS(TRIM(NETCDFFNAME),NETCDFiREC,rfld(:,:,:,4),FLDSNAME(1),FLDSUNIT(1))
 CALL READFLDS(TRIM(NETCDFFNAME),NETCDFiREC,rfld(:,:,:,5),FLDSNAME(2),FLDSUNIT(2))
 CALL READFLDS(TRIM(NETCDFFNAME),NETCDFiREC,rfld(:,:,:,6),FLDSNAME(3),FLDSUNIT(3))
 WRITE(*,'(10I4)'), NINT(rfld(:,1,1,4))
 WRITE(*,'(10I4)'), NINT(rfld(:,1,1,5))
 WRITE(*,'(10I4)'), NINT(rfld(:,1,1,6))

 NETCDFiREC=2
 CALL READFLDS(TRIM(NETCDFFNAME),NETCDFiREC,rfld(:,:,:,4),FLDSNAME(1),FLDSUNIT(1))
 CALL READFLDS(TRIM(NETCDFFNAME),NETCDFiREC,rfld(:,:,:,5),FLDSNAME(2),FLDSUNIT(2))
 CALL READFLDS(TRIM(NETCDFFNAME),NETCDFiREC,rfld(:,:,:,6),FLDSNAME(3),FLDSUNIT(3))
 WRITE(*,'(10I4)'), NINT(rfld(:,1,1,4))
 WRITE(*,'(10I4)'), NINT(rfld(:,1,1,5))
 WRITE(*,'(10I4)'), NINT(rfld(:,1,1,6))

 RETURN

 IGRFYEAR=1985.0
!CALL IGRF(IGRFYEAR,Bo,THETAN,PHIN,GEOGXYZ,GEOCXYZ)
!PRINT*, RAD2DEG(THETAN),RAD2DEG(PHIN)
!PRINT*, GEOGXYZ
!PRINT*, GEOCXYZ
!STOP

 RGEOG=6378.1
 THETAGEOG=DEG2RAD(40.0)
 PHIGEOG=DEG2RAD(0.0)

!PRINT*, RGEOG,RAD2DEG(THETAGEOG),RAD2DEG(PHIGEOG)
!CALL GEO2CD(THETAGEOG,PHIGEOG,THETAGEOC,PHIGEOC,IGRFYEAR)
!PRINT*, RGEOG,RAD2DEG(THETAGEOC),RAD2DEG(PHIGEOC)
!CALL CD2GEO(THETAGEOC,PHIGEOC,THETAGEOG,PHIGEOG,IGRFYEAR)
!PRINT*, RGEOG,RAD2DEG(THETAGEOG),RAD2DEG(PHIGEOG)
!STOP

!CALL EDBFIELD(RGEOG,THETAGEOG,PHIGEOG,IGRFYEAR,H,Hx,Hy,Hz,Hxy,Hdip,Hvar)
!PRINT*, "H = ",   H*1e9
!PRINT*, "Hx = ",  Hx*1e9
!PRINT*, "Hy = ",  Hy*1e9
!PRINT*, "Hz = ",  Hz*1e9
!PRINT*, "Hxy = ", Hxy*1e9
!PRINT*, "Hdip = ",RAD2DEG(Hdip)
!PRINT*, "Hvar = ",RAD2DEG(Hvar)
!RETURN

 PRINT*, RGEOG,RAD2DEG(THETAGEOG),RAD2DEG(PHIGEOG)
 CALL GEOG2GEOM(RGEOG,THETAGEOG,PHIGEOG,RGEOM,THETAGEOM,PHIGEOM,IGRFYEAR)
 PRINT*, RGEOM,RAD2DEG(THETAGEOM),RAD2DEG(PHIGEOM)
 CALL GEOM2GEOG(RGEOM,THETAGEOM,PHIGEOM,RGEOG,THETAGEOG,PHIGEOG,IGRFYEAR)
 PRINT*, RGEOG,RAD2DEG(THETAGEOG),RAD2DEG(PHIGEOG)
!CALL iDIAGNOSE
 RETURN

!******** Begining GFC INITIALIZATION ************
!GFC(PLTDPI)=50
!GFC(PLTFIGNUM)=1
!GFC(PLTSFIGNUM)=111
!GFC(PLTFIGHOLD)=GETVALOF("NO")
!GFC(PLTYYAXES)=GETVALOF("NO")
!GFC(DATSAVE)=GETVALOF("KEEP")
!GFC(DATTYPE)=GETVALOF("FORMATTED")
!GFC(MOVNFPS)=GETVALOF("NO")
!GFC(PNGSAVE)=GETVALOF("DELETE")
!GFC(PDFSAVE)=GETVALOF("SAVE")
!********** Ending GFC INITIALIZATION *************

!**************************!
!  INITIALIZE SYSTEM GRID  !
!**************************!
 XDIM=1
 YDIM=700
 ZDIM=200
 Xsrt=  0.0e0
 Xend=  0.0e0
 Ysrt=-16.0*(4.0*ATAN(1.0))
 Yend= 16.0*(4.0*ATAN(1.0))
 Zsrt=  0.0*(4.0*ATAN(1.0))
 Zend= 32.0*(4.0*ATAN(1.0))
 srtLAT=  0.000000
 endLAT=  0.000000
 srtLON= -0.000050
 endLON=  0.000050
 srtALT=105.000000
 endALT=105.100000
!**************************!
!  BUILDING PHYSICAL GRID  !
!**************************!
 CALL GRIDBUILDER

!**************************!
!  INITIALIZE SYSTEM TIME  !
!**************************!
!srtDT=20080320
!endDT=20080320
 srtDT=19870312
 endDT=19870312
 srtLT=120000
 endLT=130000
 stpLT=10000
 CALL DTLTMAT

!***************************************************!
!  INITIALIZE IONOSPHERE GRID WITH BACKGROUND DATA  !
!***************************************************!
 CALL RESETSYSBG
 CALL SYSBG((/"local"/),LAT1D,LON1D,ALT1D)

 IF(IonoLne) THEN
 !CALL RESETDSCHEMES(nDE)
 !FDMDSCHEME=.TRUE.
 !CALL Df(Ne(1,1,:),REAL(Z3D(1,1,:)),Lne(1,1,:),1)
 !Lne(1,1,:)=Ne(1,1,:)/Lne(1,1,:)
  Lne(1,1,1)=6000.0
 END IF

 IF(IonoVexb) THEN
 !CALL Df(EPOTEN(1,:,:),Z3D(1,:,:),Vexb(1,:,:),1)
 !Vexb(1,:,:)=-Vexb(1,:,:)/Babs(1,:,:)
  Vexb(1,1,1)=400.0
 END IF


 PRINT*, "Ne = ", Ne(1,1,1)
 PRINT*, "Te = ", Te(1,1,1)
 PRINT*, "Ti = ", Ti(1,1,1)
 PRINT*, "Bo = ", Babs(1,1,1)
 PRINT*, "Cs = ", Cs(1,1,1)
 PRINT*, "Vthe = ", Vthe(1,1,1)
 PRINT*, "Vthi = ", Vthi(1,1,1)
 PRINT*, "OMEGAce = ", OMEGAce(1,1,1)
 PRINT*, "OMEGAci = ", OMEGAci(1,1,1)
 PRINT*, "NUen = ", NUen(1,1,1)
 PRINT*, "NUin = ", NUin(1,1,1)
 PRINT*, "RHOe = ", RHOe(1,1,1)
 PRINT*, "RHOi = ", RHOi(1,1,1)
return

!******************************!
!  DECLARE THE NUMBER OF PDEs  !
!******************************!
 nDE=3

!******************************!
!  DECLARE THE SET-ID OF PDEs  !
!******************************!
 ODESETID=1

!**************************!
!  INITIALIZE FFTW SPACES  !
!**************************!
 CALL PSFFTWSET(YDIM,ZDIM)

!******************************************************!
!  SOLVING THE EIGENVALUE PROBLEM FOR THE SET OF PDEs  !
!******************************************************!
 CALL SYSEIGENVALS(ODESETID,nDE,(/"local"/),X3D,Y3D,Z3D)
!CALL SYSEIGENVV(ODESETID,nDE,(/"local"/),X3D,Y3D,Z3D)
!*****************************************************************!
!  PLOTING THE GROWTH-RATE AND PHASE VELOCITY OF THE SET OF PDEs  !
!*****************************************************************!
 CALL PLOTEIGENS("max-positive","yz","FIGNOHOLD",400.0,6000.0)
 RETURN


!**********************************************!
!  ALLOCATE THE FIELDS IN R-SPACE AND K-SPACE  !
!**********************************************!
 ALLOCATE(rfld(XDIM,YDIM,ZDIM,nDE))
 ALLOCATE(kfld(XDIM,YDIM/2+1,ZDIM,nDE))
 rfld=0.0d0
 kfld=CMPLX(0.0d0,0.0d0)

!*********************************************!
! INITIALIZING THE FIELDS FOR THE CURRENT RUN !
!*********************************************!
 IF(SYSRESTART) THEN
 !*******************!
 ! READ SAVED FIELDS !
 !*******************!
 ! FOR THE SYSTEMRESTART WE READ THE FIELDS AT THE LAST TIME STEP OF
 ! THE PREVIOUS RUN FROM THE netCDF FILE
  CALL NETCDF_VAR_DIMS(TRIM(NETCDFFNAME),TRIM(FLDSNAME(1)),MATDIMS4D)
  NETCDFiREC=MATDIMS4D(4)
! CALL READNETCDF(TRIM(NETCDFFNAME),NETCDFiREC,X3D(:,1,1),Y3D(1,:,1),Z3D(1,1,:),t_strt,rfld,DIMSNAME,DIMSUNIT,FLDSNAME,FLDSUNIT)

 !*******************************!
 ! CONVERT THE FIELDS TO K-SPACE !
 !*******************************!
  DO fldID=1,nDE
   CALL FFTW(rfld(1,:,:,fldID),kfld(1,:,:,fldID))
  END DO

 !*******************************************!
 !  ALLOCATE AND INITIALIZING SPECTRAL GRID  !
 !*******************************************!
  IF(ALLOCATED(iKy2D).EQV..FALSE.) ALLOCATE(iKy2D(YDIM/2+1,ZDIM))
  IF(ALLOCATED(iKz2D).EQV..FALSE.) ALLOCATE(iKz2D(YDIM/2+1,ZDIM))
  CALL KGRIDBUILDER(Y3D(1,:,:),iKy2D)
  CALL KGRIDBUILDER(Z3D(1,:,:),iKz2D)

 ELSE
 !***********************************************************************!
 ! FIND THE EIGENVECTORS OF THE THE SET OF PDEs TO INITIALIZE THE FIELDS !
 !***********************************************************************!
  CALL SYSEIGENFLDS(ODESETID,nDE,"positive-max-modes",(/"local"/),X3D,Y3D,Z3D,rfld,kfld)

 !*****************!
 ! SAVE THE FIELDS !
 !*****************!
! CALL SAVENETCDF("iBRIDGE.nc",NETCDFiREC,X3D(:,1,1),Y3D(1,:,1),Z3D(1,1,:),t_crnt,rfld,DIMSNAME,DIMSUNIT,FLDSNAME,FLDSUNIT)
 END IF

!*******************************!
! INITIALIZE THE TYPES OF PLOTS !
!*******************************!
 kfldPlot=.true.      ! SET TO .TRUE. WHEN WE NEED TO PLOT THE SPECTRAL CONTENTS OF THE FIELDS
 rfldPlot=.true.      ! SET TO .TRUE. WHEN WE NEED TO PLOT THE PHYSICAL CONTENTS OF THE FIELDS
 efldPlot=.true.      ! SET TO .TRUE. WHEN WE NEED TO PLOT THE ENERGY   CONTENTS OF THE FIELDS
 CALL FieldsPlot("yz",X3D,Y3D,Z3D,rfld,1,REAL(t_strt,4))            ! PLOT THE FIELDS DIRECTLY

 STOP

!***************************************!
! NUMERICAL SOLUTION OF THE SET OF PDEs !
!***************************************!
 DO iTIMER=1,ntstp
  t_crnt=t_strt+iTIMER*tstp
  CALL DESOLVER(ODESETID,Y3D(1,:,:),Z3D(1,:,:),kfld(1,:,:,:),t_crnt,tstp)
  IF(MOD(iTIMER,50).EQ.0) THEN
   CALL WHATISTIME(NOWTIME)
   WRITE(*,'(A9,I7,7I5)') "iTIMER = ", iTIMER,NOWTIME
  END IF
  IF(MOD(iTIMER,500).EQ.0) THEN
   DO fldID=1,nDE
    CALL iFFTW(kfld(1,:,:,fldID),rfld(1,:,:,fldID))
   END DO
   NETCDFiREC=NETCDFiREC+1
!  CALL SAVENETCDF("iBRIDGE.nc",NETCDFiREC,X3D(:,1,1),Y3D(1,:,1),Z3D(1,1,:),t_crnt,rfld,DIMSNAME,DIMSUNIT,FLDSNAME,FLDSUNIT)
   CALL FieldsPlot("yz",X3D,Y3D,Z3D,rfld,iTIMER+1,REAL(t_crnt,4))
  END IF
 END DO

!*******************************!
! INITIALIZE THE TYPES OF PLOTS !
!*******************************!
!kfldPlot=.true.             ! SET TO .TRUE. WHEN WE NEED TO PLOT THE SPECTRAL CONTENTS OF THE FIELDS
!rfldPlot=.true.             ! SET TO .TRUE. WHEN WE NEED TO PLOT THE PHYSICAL CONTENTS OF THE FIELDS
!efldPlot=.true.             ! SET TO .TRUE. WHEN WE NEED TO PLOT THE ENERGY   CONTENTS OF THE FIELDS
!CALL FieldsPlot("yz",X3D,Y3D,Z3D,"iBRIDGE.nc",FLDSNAME,DIMSNAME)   ! PLOT THE FIELDS FROM THE netCDF


 PRINT*,"RUN HAS BEEN DONE!"
END PROGRAM BRIDGE

