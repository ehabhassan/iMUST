MODULE iTEST

 CONTAINS

 SUBROUTINE iDIAGNOSE
  USE IONOBG
  USE SYSINFO
  USE GEOAXES
  USE SYSGRIDS
  USE MATHSUBS
  USE RDSVNETCDF
  USE INTERPSUBS
  USE DOPERATORS
  USE SYSGRAPHICS
  USE FORTGRAPHER
  IMPLICIT NONE
  INTEGER(4)            :: iD1,iD2,iD3
  INTEGER(4)            :: DIM1,DIM2,DIM3
  INTEGER(4)            :: iDE
  INTEGER(4)            :: iKx,iKy,iKz
  INTEGER(4)            :: iXDIM,iYDIM,iZDIM
  INTEGER(4)            :: iREC,iREC1,iREC2,iREC3
  INTEGER(4)            :: nREC,nREC1,nREC2,nREC3
  INTEGER(4)            :: inetCDF
  INTEGER(4)            :: LOOPSID,LOOPEID
  INTEGER(4)            :: nDIMS,nVARS
  INTEGER(4)            :: MATDIMS4D(4)
  INTEGER(4)            :: MATDIMS1D(1)
  INTEGER(4)            :: nDE  ! NUMBER OF DIFFERENTIAL EQUATIONS IN THE SYSTEM
  INTEGER(4)            :: nSAT ! NUMBER OF STATISTICS NEEDS TO BE DONE 

  COMPLEX(8),ALLOCATABLE:: kfld(:,:,:,:)
  REAL(8),   ALLOCATABLE:: rfld(:,:,:,:)
  REAL(8),   ALLOCATABLE:: rTIME(:)
  REAL(8),   ALLOCATABLE:: EyP(:),EyN(:)
  REAL(8),   ALLOCATABLE:: EzP(:),EzN(:)
  REAL(8),   ALLOCATABLE:: rNeSTAT(:,:),rNeSTAT2D(:,:,:)
  REAL(8),   ALLOCATABLE:: rEySTAT(:,:),rEySTAT2D(:,:,:)
  REAL(8),   ALLOCATABLE:: rEzSTAT(:,:),rEzSTAT2D(:,:,:)
  REAL(8),   ALLOCATABLE:: rVeySTAT(:,:),rVeySTAT2D(:,:,:)
  REAL(8),   ALLOCATABLE:: rVezSTAT(:,:),rVezSTAT2D(:,:,:)
  REAL(8),   ALLOCATABLE:: rViySTAT(:,:),rViySTAT2D(:,:,:)
  REAL(8),   ALLOCATABLE:: rVizSTAT(:,:),rVizSTAT2D(:,:,:)

  REAL(8),   ALLOCATABLE:: rNe(:,:,:,:)
  REAL(8),   ALLOCATABLE:: rEy(:,:,:,:),  rEz(:,:,:,:)
  REAL(8),   ALLOCATABLE:: rVey(:,:,:,:),rVez(:,:,:,:)
  REAL(8),   ALLOCATABLE:: rViy(:,:,:,:),rViz(:,:,:,:)

  COMPLEX(8),ALLOCATABLE:: kNe(:,:,:,:)
  COMPLEX(8),ALLOCATABLE:: kEy(:,:,:,:),  kEz(:,:,:,:)
  COMPLEX(8),ALLOCATABLE:: kVey(:,:,:,:),kVez(:,:,:,:)
  COMPLEX(8),ALLOCATABLE:: kViy(:,:,:,:),kViz(:,:,:,:)

  REAL(8),   ALLOCATABLE:: rFLXnvey(:,:,:,:),rFLXnvez(:,:,:,:)
  REAL(8),   ALLOCATABLE:: rFLXnviy(:,:,:,:),rFLXnviz(:,:,:,:)
  REAL(4)               :: RUNTIME

  LOGICAL               :: iNeSAT=.TRUE.
  LOGICAL               :: iEySAT=.FALSE.
  LOGICAL               :: iEzSAT=.FALSE.
  LOGICAL               :: iVeySAT=.FALSE.
  LOGICAL               :: iVezSAT=.FALSE.
  LOGICAL               :: iViySAT=.FALSE.
  LOGICAL               :: iVizSAT=.FALSE.

  INTEGER(4)            :: NOWTIME(7)
  CHARACTER(100)        :: TEXTPAD
  CHARACTER(100)        :: LocGA(8)=""
  CHARACTER(1)          :: LINECLR(7)=(/CHARACTER(1) :: "b","r","k","y","g","c","m"/)
  CHARACTER(1)          :: LINESHP(8)=(/CHARACTER(1) :: ":","--","-.","*","+","x","o","^"/)
! CHARACTER(35)         :: NETCDFFNAME(5)=(/CHARACTER(35) :: "keU120.nc","keU120noon19870312n.nc","keU120noon19870312.nc","keU120noon20080320.nc","keU123.nc"/)
! CHARACTER(25)         :: NETCDFFNAME(4)=(/CHARACTER(25) :: "keU120noon19870312T3.nc","keU120noon19870312n.nc","keU120noon19870312.nc","keU120noon20080320.nc"/)
  CHARACTER(40)         :: NETCDFFNAME(1)=""
! CHARACTER(30)         :: NETCDFFNAME(2)=(/CHARACTER(30) :: "keU120noon19870312Ln16.nc"/)
! CHARACTER(50)         :: NETCDFFNAME(3)=(/CHARACTER(50) :: "keU120noon19870312Ln1000.nc","keU120noon19870312Ln4000.nc","keU120noon19870312T3.nc"/)
  CHARACTER(25)         :: FLDSNAME(3)=(/CHARACTER(25) :: "Density","Electric Potential","Ion Velocity Potential"/)
  CHARACTER(10)         :: FLDSUNIT(3)=(/CHARACTER(10) :: "meter^-3","volts","meter^2/s"/)
  CHARACTER(10)         :: DIMSNAME(4)=(/CHARACTER(10) :: "X","Y","Z","T"/)
  CHARACTER(10)         :: DIMSUNIT(4)=(/CHARACTER(10) :: "meter","meter","meter","second"/)

!******** GRAPHICS VARIABLES DECLARATION *********
  INTEGER(4)            :: GFC(10)
  CHARACTER(100)        :: GA(100)=""
  CHARACTER(100)        :: GX(100)=""

!******** GFC INITIALIZATION ************
  GFC(PLTDPI)=200
  GFC(PLTFIGNUM)=1
  GFC(PLTSFIGNUM)=111
  GFC(PLTFIGHOLD)=GETVALOF("NO")
  GFC(PLTYYAXES)=GETVALOF("NO")
  GFC(DATSAVE)=GETVALOF("KEEP")
  GFC(DATTYPE)=GETVALOF("FORMATTED")
  GFC(MOVNFPS)=GETVALOF("NO")
  GFC(PNGSAVE)=GETVALOF("DELETE")
  GFC(PDFSAVE)=GETVALOF("SAVE")

 !CALL WHATISTIME(NOWTIME)
 !WRITE(*,'(A14,I7,7I5)') "RUN STARTS ON:", NOWTIME

! NETCDFFNAME(1)="keU120.nc"
  NETCDFFNAME(1)="OPENMPtst.nc"
 !NETCDFFNAME(1)="keU123.nc"
! NETCDFFNAME(1)="keU120noon20080320.nc"
! NETCDFFNAME(2)="keU120noon19870312n.nc"

! NETCDFFNAME(1)="keU120noon19870312Ln631.nc"
! NETCDFFNAME(2)="keU120noon19870312Ln16.nc"

  CALL NETCDF_INFO(TRIM(NETCDFFNAME(1)),nDIMS,nVARS)
  CALL NETCDF_VAR_DIMS(TRIM(NETCDFFNAME(1)),TRIM(FLDSNAME(1)),MATDIMS4D)
 !nREC1=MATDIMS4D(4)
 !CALL NETCDF_VAR_DIMS(TRIM(NETCDFFNAME(2)),TRIM(FLDSNAME(1)),MATDIMS4D)
 !nREC2=MATDIMS4D(4)
 !nREC=MIN(nREC1,nREC2)
  nREC=MATDIMS4D(4)
  nDE=nVARS-nDIMS
  XDIM=MATDIMS4D(1)
  YDIM=MATDIMS4D(2)
  ZDIM=MATDIMS4D(3)
  ALLOCATE(X1D(XDIM))
  ALLOCATE(Y1D(YDIM))
  ALLOCATE(Z1D(ZDIM))
  ALLOCATE(rfld(XDIM,YDIM,ZDIM,nDE))
  ALLOCATE(kfld(XDIM,YDIM/2+1,ZDIM,nDE))
  CALL READNETCDF(TRIM(NETCDFFNAME(1)),nREC,X1D,Y1D,Z1D,RUNTIME,rfld,DIMSNAME,DIMSUNIT,FLDSNAME,FLDSUNIT)
  ALLOCATE(X3D(XDIM,YDIM,ZDIM))
  ALLOCATE(Y3D(XDIM,YDIM,ZDIM))
  ALLOCATE(Z3D(XDIM,YDIM,ZDIM))
  CALL MESHGRID(X1D,Y1D,Z1D,X3D,Y3D,Z3D)

  Xsrt=X1D(1)
  Xend=X1D(XDIM)
  Ysrt=Y1D(1)
  Yend=Y1D(YDIM)
  Zsrt=Z1D(1)
  Zend=Z1D(ZDIM)
  srtLAT= -0.000000
  endLAT=  0.000000
  srtLON= -0.000050
  endLON=  0.000050
  srtALT=105.000000
  endALT=105.100000
  CALL GRIDBUILDER
  CALL DTLTMAT
  CALL SYSBG((/"local"/),LAT1D,LON1D,ALT1D)
  IF(IonoLne) THEN
   Lne(1,:,:)=6000.0
  END IF

  IF(IonoVexb) THEN
   Vexb(1,:,:)=400.0
  END IF

  IF(ALLOCATED(iKy2D).EQV..FALSE.) ALLOCATE(iKy2D(YDIM/2+1,ZDIM))
  IF(ALLOCATED(iKz2D).EQV..FALSE.) ALLOCATE(iKz2D(YDIM/2+1,ZDIM))
  CALL KGRIDBUILDER(Y3D(1,:,:),iKy2D)
  CALL KGRIDBUILDER(Z3D(1,:,:),iKz2D)

! nREC2=0
! DO iREC=1,nREC
!  CALL READNETCDF(TRIM(NETCDFFNAME(1)),iREC,X3D(:,1,1),Y3D(1,:,1),Z3D(1,1,:),RUNTIME,rfld,DIMSNAME,DIMSUNIT,FLDSNAME,FLDSUNIT)
!  IF(RUNTIME.GE.0.9) nREC2=nREC2+1
! END DO
! nREC1=nREC-nREC2+1

!!rfldPlot=.true.
!!kfldPlot=.true.
! GA(IPLOT)="contourf"
! GA(CBAR)="yes"
! GA(XLABEL)="East-West (m)"
! GA(SUPERTITLE)="Density Perturbation"
!!GA(FIGSIZE)="(40,32)"
! DO iREC=570,1001
!  CALL READNETCDF(TRIM(NETCDFFNAME(1)),iREC,X3D(:,1,1),Y3D(1,:,1),Z3D(1,1,:),RUNTIME,rfld,DIMSNAME,DIMSUNIT,FLDSNAME,FLDSUNIT)
! !CALL FORTPLOT(GFC,Y3D(1,:,:),Z3D(1,:,:),100.0*rfld(1,:,:,1),"",GA)
! !IF(NINT(1000.0*RUNTIME).EQ.600) THEN
!  IF(iREC.EQ.583) THEN
!  !WRITE(TEXTPAD,'(i3)'), NINT(1000.0*RUNTIME)
!   GA(TITLE)="Runtime = 600 ms"
!   GA(TITLEfs)="9"
!   GA(YLABEL)="Altitude (m)"
!   GFC(PLTFIGHOLD)=GETVALOF("YES")
!   GFC(PLTSFIGNUM)=221
!   CALL FORTPLOT(GFC,Y3D(1,:,:),Z3D(1,:,:),100.0*rfld(1,:,:,1),"",GA)
! !ELSEIF(NINT(1000.0*RUNTIME).EQ.1000) THEN
!  ELSEIF(iREC.EQ.1001) THEN
!   GA(TITLE)="Runtime = 1000 ms"
!   GA(TITLEfs)="9"
!   GA(YLABEL)=""
!   GA(CBARTITLE)="$\delta{n}$"
!   GFC(PLTFIGHOLD)=GETVALOF("NO")
!   GFC(PLTSFIGNUM)=222
!   CALL FORTPLOT(GFC,Y3D(1,:,:),Z3D(1,:,:),100.0*rfld(1,:,:,1),"",GA)
!  END IF
! END DO
!!CALL READNETCDF(TRIM(NETCDFFNAME(5)),420,X3D(:,1,1),Y3D(1,:,1),Z3D(1,1,:),RUNTIME,rfld,DIMSNAME,DIMSUNIT,FLDSNAME,FLDSUNIT)
!!CALL FieldsPlot("yz",X3D,Y3D,Z3D,rfld,420,RUNTIME)
! RETURN

  IF(ALLOCATED(rNe)    .EQV..FALSE.) ALLOCATE(rNe(XDIM,YDIM,ZDIM,nREC))
  IF(ALLOCATED(rEy)    .EQV..FALSE.) ALLOCATE(rEy(XDIM,YDIM,ZDIM,nREC))
  IF(ALLOCATED(rEz)    .EQV..FALSE.) ALLOCATE(rEz(XDIM,YDIM,ZDIM,nREC))
 !IF(ALLOCATED(rVey)   .EQV..FALSE.) ALLOCATE(rVey(XDIM,YDIM,ZDIM,nREC))
 !IF(ALLOCATED(rVez)   .EQV..FALSE.) ALLOCATE(rVez(XDIM,YDIM,ZDIM,nREC))
 !IF(ALLOCATED(rViy)   .EQV..FALSE.) ALLOCATE(rViy(XDIM,YDIM,ZDIM,nREC))
 !IF(ALLOCATED(rViz)   .EQV..FALSE.) ALLOCATE(rViz(XDIM,YDIM,ZDIM,nREC))

  IF(ALLOCATED(kNe)    .EQV..FALSE.) ALLOCATE(kNe(XDIM,YDIM/2+1,ZDIM,nREC))
  IF(ALLOCATED(kEy)    .EQV..FALSE.) ALLOCATE(kEy(XDIM,YDIM/2+1,ZDIM,nREC))
  IF(ALLOCATED(kEz)    .EQV..FALSE.) ALLOCATE(kEz(XDIM,YDIM/2+1,ZDIM,nREC))
 !IF(ALLOCATED(kVey)   .EQV..FALSE.) ALLOCATE(kVey(XDIM,YDIM/2+1,ZDIM,nREC))
 !IF(ALLOCATED(kVez)   .EQV..FALSE.) ALLOCATE(kVez(XDIM,YDIM/2+1,ZDIM,nREC))
 !IF(ALLOCATED(kViy)   .EQV..FALSE.) ALLOCATE(kViy(XDIM,YDIM/2+1,ZDIM,nREC))
 !IF(ALLOCATED(kViz)   .EQV..FALSE.) ALLOCATE(kViz(XDIM,YDIM/2+1,ZDIM,nREC))

  nSAT=5
  IF(ALLOCATED(rTIME)  .EQV..FALSE.) ALLOCATE(rTIME(nREC))
  IF(ALLOCATED(rNeSTAT) .EQV..FALSE.) ALLOCATE(rNeSTAT(nREC,nSAT))
  IF(ALLOCATED(rEySTAT) .EQV..FALSE.) ALLOCATE(rEySTAT(nREC,nSAT))
  IF(ALLOCATED(rEzSTAT) .EQV..FALSE.) ALLOCATE(rEzSTAT(nREC,nSAT))
 !IF(ALLOCATED(rVeySTAT).EQV..FALSE.) ALLOCATE(rVeySTAT(nREC,nSAT))
 !IF(ALLOCATED(rVezSTAT).EQV..FALSE.) ALLOCATE(rVezSTAT(nREC,nSAT))
 !IF(ALLOCATED(rViySTAT).EQV..FALSE.) ALLOCATE(rViySTAT(nREC,nSAT))
 !IF(ALLOCATED(rVizSTAT).EQV..FALSE.) ALLOCATE(rVizSTAT(nREC,nSAT))

 !IF(ALLOCATED(rNeSAT2D) .EQV..FALSE.) ALLOCATE(rNeSAT2D(nREC,YDIM,nSAT))
 !IF(ALLOCATED(rEySAT2D) .EQV..FALSE.) ALLOCATE(rEySAT2D(nREC,YDIM,nSAT))
 !IF(ALLOCATED(rEzSAT2D) .EQV..FALSE.) ALLOCATE(rEzSAT2D(nREC,YDIM,nSAT))
 !IF(ALLOCATED(rVeySAT2D).EQV..FALSE.) ALLOCATE(rVeySAT2D(nREC,YDIM,nSAT))
 !IF(ALLOCATED(rVezSAT2D).EQV..FALSE.) ALLOCATE(rVezSAT2D(nREC,YDIM,nSAT))
 !IF(ALLOCATED(rViySAT2D).EQV..FALSE.) ALLOCATE(rViySAT2D(nREC,YDIM,nSAT))
 !IF(ALLOCATED(rVizSAT2D).EQV..FALSE.) ALLOCATE(rVizSAT2D(nREC,YDIM,nSAT))

  LOOPSID=1
  LOOPEID=1
 !LOOPSID=1
 !LOOPEID=SIZE(NETCDFFNAME)
  DO inetCDF=LOOPSID,LOOPEID
   TEXTPAD="PROCESSING "//TRIM(NETCDFFNAME(inetCDF))//" ..."
   WRITE(*,TRIM(STRFMT(LEN(TRIM(TEXTPAD))))) TRIM(TEXTPAD) 
   DO iREC=1,nREC
  !DO iREC=1,nREC2
   !READ THE DATA IN THE netCDF FILE
    CALL READNETCDF(TRIM(NETCDFFNAME(inetCDF)),iREC,X3D(:,1,1),Y3D(1,:,1),Z3D(1,1,:),RUNTIME,rfld,DIMSNAME,DIMSUNIT,FLDSNAME,FLDSUNIT)
   !CALL READNETCDF(TRIM(NETCDFFNAME(inetCDF)),iREC+nREC1-1,X3D(:,1,1),Y3D(1,:,1),Z3D(1,1,:),RUNTIME,rfld,DIMSNAME,DIMSUNIT,FLDSNAME,FLDSUNIT)
   !FILLING THE TIME ARRAY
    rTIME(iREC)=RUNTIME
    DO iDE=1,nDE
    !CONVERTING THE FIELDS FROM THE R-SPACE TO K-SPACE
     CALL FFT(rfld(1,:,:,iDE),kfld(1,:,:,iDE))
    END DO
   !FILLING THE DENSITY MATRIX
    rNe(1,:,:,iREC)=rfld(1,:,:,1)
  ! kNe(1,:,:,iREC)=kfld(1,:,:,1)
   !CALCULATING THE ELECTRIC FIELD COMPONENTS
  ! kEy(1,:,:,iREC)=-iKy2D*((KB*Te(1,1,1)/Qe)*kfld(1,:,:,2))
  ! kEz(1,:,:,iREC)=-iKz2D*((KB*Te(1,1,1)/Qe)*kfld(1,:,:,2))
  ! CALL iFFT( kEy(1,:,:,iREC), rEy(1,:,:,iREC))
  ! CALL iFFT( kEz(1,:,:,iREC), rEz(1,:,:,iREC))
  ! kEy(1,:,:,1)=-iKy2D*((KB*Te(1,1,1)/Qe)*kfld(1,:,:,2))
  ! kEz(1,:,:,1)=-iKz2D*((KB*Te(1,1,1)/Qe)*kfld(1,:,:,2))
  ! CALL iFFT( kEy(1,:,:,1), rEy(1,:,:,1))
  ! CALL iFFT( kEz(1,:,:,1), rEz(1,:,:,1))

!   IF(ALLOCATED(EyP)) DEALLOCATE(EyP)
!   iD1=0;
!   DO iYDIM=1,YDIM
!    DO iZDIM=1,ZDIM
!     IF(rEy(1,iYDIM,iZDIM,1).GT.0) iD1=iD1+1
!    END DO
!   END DO
!   ALLOCATE(EyP(iD1))
!   iD1=1
!   DO iYDIM=1,YDIM
!    DO iZDIM=1,ZDIM
!     IF(rEy(1,iYDIM,iZDIM,1).GT.0) THEN
!      EyP(iD1)=rEy(1,iYDIM,iZDIM,1)
!      iD1=iD1+1
!     END IF
!    END DO
!   END DO
!   rEySTAT(iREC,3)=STDEV(EyP)
!   DEALLOCATE(EyP)

!   iD1=0
!   DO iYDIM=1,YDIM
!    DO iZDIM=1,ZDIM
!     IF(rEy(1,iYDIM,iZDIM,1).LT.0) iD1=iD1+1
!    END DO
!   END DO
!   ALLOCATE(EyN(iD1))
!   iD1=1
!   DO iYDIM=1,YDIM
!    DO iZDIM=1,ZDIM
!     IF(rEy(1,iYDIM,iZDIM,1).LT.0) THEN
!      EyN(iD1)=rEy(1,iYDIM,iZDIM,1)
!      iD1=iD1+1
!     END IF
!    END DO
!   END DO
!   rEySTAT(iREC,4)=STDEV(EyN)
!   DEALLOCATE(EyN)

!   iD1=0
!   DO iYDIM=1,YDIM
!    DO iZDIM=1,ZDIM
!     IF(rEz(1,iYDIM,iZDIM,1).GT.0) iD1=iD1+1
!    END DO
!   END DO
!   ALLOCATE(EzP(iD1))
!   iD1=1
!   DO iYDIM=1,YDIM
!    DO iZDIM=1,ZDIM
!     IF(rEz(1,iYDIM,iZDIM,1).GT.0) THEN
!      EzP(iD1)=rEz(1,iYDIM,iZDIM,1)
!      iD1=iD1+1
!     END IF
!    END DO
!   END DO
!   rEzSTAT(iREC,3)=STDEV(EzP)
!   DEALLOCATE(EzP)

!   iD1=0
!   DO iYDIM=1,YDIM
!    DO iZDIM=1,ZDIM
!     IF(rEz(1,iYDIM,iZDIM,1).LT.0) iD1=iD1+1
!    END DO
!   END DO
!   ALLOCATE(EzN(iD1))
!   iD1=1
!   DO iYDIM=1,YDIM
!    DO iZDIM=1,ZDIM
!     IF(rEz(1,iYDIM,iZDIM,1).LT.0) THEN
!      EzN(iD1)=rEz(1,iYDIM,iZDIM,1)
!      iD1=iD1+1
!     END IF
!    END DO
!   END DO
!   rEzSTAT(iREC,4)=STDEV(EzN)
!   DEALLOCATE(EzN)


   !CALCULATING THE ION VELOCITY COMPONENTS
  ! kViy(1,:,:,iREC)=-iKy2D*kfld(1,:,:,3)
  ! kViz(1,:,:,iREC)=-iKz2D*kfld(1,:,:,3)
  ! CALL iFFT( kViy(1,:,:,iREC), rViy(1,:,:,iREC))
  ! CALL iFFT( kViz(1,:,:,iREC), rViz(1,:,:,iREC))
   !CALCULATING THE FLUX OF ION DENSITY
   !rFLXnviy(1,:,:,iREC)=rNe(1,:,:,iREC)*rViy(1,:,:,iREC)
   !rFLXnviz(1,:,:,iREC)=rNe(1,:,:,iREC)*rViz(1,:,:,iREC)
   !CALCULATING THE ELECTRON VELOCITY COMPONENTS
  ! rVey(1,:,:,iREC)=-MUeH(1,1,1)*rEz(1,:,:,iREC)
  ! rVez(1,:,:,iREC)= MUeH(1,1,1)*rEy(1,:,:,iREC)
  ! CALL FFT( rVey(1,:,:,iREC), kVey(1,:,:,iREC))
  ! CALL FFT( rVez(1,:,:,iREC), kVez(1,:,:,iREC))
   !CALCULATING THE FLUX OF ELECTRON DENSITY
   !rFLXnvey(1,:,:,iREC)=ABS(kNe(1,:,:,iREC))*COS(ATAN2(AIMAG(kNe(1,:,:,iREC)),REAL(kNe(1,:,:,iREC))))*ABS(kVey(1,:,:,iREC))*COS(ATAN2(AIMAG(kVey(1,:,:,iREC)),REAL(kVey(1,:,:,iREC))))
   !rFLXnvez(1,:,:,iREC)=ABS(kNe(1,:,:,iREC))*COS(ATAN2(AIMAG(kNe(1,:,:,iREC)),REAL(kNe(1,:,:,iREC))))*ABS(kVez(1,:,:,iREC))*COS(ATAN2(AIMAG(kVez(1,:,:,iREC)),REAL(kVez(1,:,:,iREC))))
   !rFLXnvey(1,:,:,iREC)=rNe(1,:,:,iREC)*rVey(1,:,:,iREC)
   !rFLXnvez(1,:,:,iREC)=rNe(1,:,:,iREC)*rVez(1,:,:,iREC)

    IF(iNeSAT) THEN
    !CALCULATE THE MEAN VALUE OF THE PERTURBED DENSITY
   ! rNeSTAT(iREC,1)=MEAN(rNe(1,:,:,iREC))
    !CALCULATE THE STANDARD DEVIATION OF THE PERTURBED DENSITY
   ! rNeSTAT(iREC,2)=STDEV(rNe(1,:,:,iREC))
    !CALCULATE THE SKEWNESS OF THE PERTURBED DENSITY
   ! rNeSTAT(iREC,3)=SKEWNESS(rNe(1,:,:,iREC))
    !CALCULATE THE KURTOSIS OF THE PERTURBED DENSITY
   ! rNeSTAT(iREC,4)=KURTOSIS(rNe(1,:,:,iREC))
    !CALCULATE THE MAXIMUM OF THE PERTURBED DENSITY
     rNeSTAT(iREC,5)=MAXVAL(rNe(1,:,:,iREC))
    END IF
  
    IF(iEySAT) THEN
    !CALCULATE THE MEAN VALUE OF THE PERTURBED HORIZONTAL ELECTRIC FIELD
   ! rEySTAT(iREC,1)=MEAN(rEy(1,:,:,iREC))
    !CALCULATE THE STANDARD DEVIATION OF THE PERTURBED HORIZONTAL ELECTRIC FIELD
   ! rEySTAT(iREC,2)=STDEV(rEy(1,:,:,1))
   ! rEySTAT(iREC,2)=STDEV(rEy(1,:,:,iREC))
   ! rEySTAT(iREC,3)=STDEV(EyP)
   ! rEySTAT(iREC,4)=STDEV(EyN)
    !CALCULATE THE SKEWNESS OF THE PERTURBATION OF THE DENSITY
   ! rEySTAT(iREC,3)=SKEWNESS(rEy(1,:,:,iREC))
    !CALCULATE THE KURTOSIS OF THE PERTURBATION OF THE DENSITY
   ! rEySTAT(iREC,4)=KURTOSIS(rEy(1,:,:,iREC))
    !CALCULATE THE MAXIMUM OF THE PERTURBATION OF THE DENSITY
   ! rEySTAT(iREC,5)=MAXVAL(rEy(1,:,:,1))
    END IF

    IF(iEzSAT) THEN
    !CALCULATE THE MEAN VALUE OF THE PERTURBED HORIZONTAL ELECTRIC FIELD
   ! rEzSTAT(iREC,1)=MEAN(rEz(1,:,:,iREC))
    !CALCULATE THE STANDARD DEVIATION OF THE PERTURBED HORIZONTAL ELECTRIC FIELD
   ! rEzSTAT(iREC,2)=STDEV(rEz(1,:,:,1))
   ! rEzSTAT(iREC,2)=STDEV(rEz(1,:,:,iREC))
   ! rEzSTAT(iREC,3)=STDEV(EzP)
   ! rEzSTAT(iREC,4)=STDEV(EzN)
    !CALCULATE THE SKEWNESS OF THE PERTURBATION OF THE DENSITY
   ! rEzSTAT(iREC,3)=SKEWNESS(rEz(1,:,:,iREC))
    !CALCULATE THE KURTOSIS OF THE PERTURBATION OF THE DENSITY
   ! rEzSTAT(iREC,4)=KURTOSIS(rEz(1,:,:,iREC))
    !CALCULATE THE MAXIMUM OF THE PERTURBATION OF THE DENSITY
   ! rEzSTAT(iREC,5)=MAXVAL(rEz(1,:,:,1))
    END IF
   END DO

   IF(iNETCDF.EQ.LOOPEID) THEN
    GFC(PLTFIGHOLD)=GETVALOF("NO")
   ELSE
    GFC(PLTFIGHOLD)=GETVALOF("YES")
   END IF
   GA(IPLOT)="plot"
   GA(YTICKSFMT)="0.1e"
   GA(YTICKSfs)="9"
  !IF(inetCDF.EQ.1) GA(IPLOTLABEL)="$[\\varphi,\\nabla^2\\varphi]$"
  !IF(inetCDF.EQ.2) GA(IPLOTLABEL)="$[\\varphi,\\nabla^2\\varphi]+[\\varphi,\\nabla^2{n}]$"
  !IF(inetCDF.EQ.3) GA(IPLOTLABEL)="$[\\varphi,\\nabla^2\\varphi]+[n,\\nabla^2\\varphi]$"
  !IF(inetCDF.EQ.4) GA(IPLOTLABEL)="$[\\varphi,\\nabla^2\\varphi]+[\\varphi,\\nabla^2{n}]+[n,\\nabla^2\\varphi]$"
   IF(inetCDF.EQ.1) GA(IPLOTLABEL)="$\\delta{E_{y,Max}}$"
   IF(inetCDF.EQ.2) GA(IPLOTLABEL)="$\\delta{E_{y,Min}}$"
   GA(LINEPROP)=TRIM(LINECLR(inetCDF))
   GA(LEGEND)="YES"
   GA(LEGENDfs)="7"
   GA(LEGENDPOS)="best"

   GFC(PLTFIGHOLD)=GETVALOF("NO")
  !GFC(PLTFIGHOLD)=GETVALOF("YES")
   GFC(PLTSFIGNUM)=311
   IF(iNeSAT) THEN
   !GA(TITLE)="Perturbed Density"
   !GA(XLABEL)="Time (sec)"

  ! GFC(PLTFIGNUM)=1
  ! GA(TITLE)="Mean of Perturbed Density"
  ! GA(YLABEL)="$<\\delta{n}>$(%)"
  ! CALL FORTPLOT(GFC,rTIME,100.0*rNeSTAT(:,1),"",GA)
  !!CALL FORTPLOT(GFC,rTIME,100.0*rNeSTAT(:,1),"NeSTAT",GA)

  ! GFC(PLTFIGNUM)=2
  ! GA(TITLE)="Standard Deviation of Perturbed Density"
  ! GA(XLABEL)="Time (sec)"
  ! GA(YLABEL)="$\sigma_{\delta{n}}$(%)"
  ! GX(1)="ax.set_xticks([])"
  ! CALL FORTPLOT(GFC,rTIME,100.0*rNeSTAT(:,2),"",GA)
  ! GX(1)=""
   !CALL FORTPLOT(GFC,rTIME,100.0*rNeSTAT(:,2),"NeSTAT",GA)

  ! GFC(PLTFIGNUM)=3
  ! GA(TITLE)="Skewness of Perturbed Density"
  ! GA(YLABEL)="${\\delta{n}}$"
  !!GA(YLABEL)="${\\delta{n}}$ Skewness"
  !!CALL FORTPLOT(GFC,rTIME,1.0*rNeSTAT(:,3),"NeSTAT",GA)
  ! CALL FORTPLOT(GFC,rTIME,1.0*rNeSTAT(:,3),"",GA)

  ! GFC(PLTFIGNUM)=4
  ! GA(TITLE)="Kurtosis of Perturbed Density"
  ! GA(YLABEL)="${\\delta{n}}$"
  !!GA(YLABEL)="${\\delta{n}}$ Kurtosis"
  !!CALL FORTPLOT(GFC,rTIME,1.0*rNeSTAT(:,4),"NeSTAT",GA)
  ! CALL FORTPLOT(GFC,rTIME,1.0*rNeSTAT(:,4),"",GA)

    GFC(PLTFIGNUM)=5
    GA(TITLE)="Maximum Value of Perturbed Density"
   !GA(XLABEL)="Time (sec)"
    GA(YLABEL)="${\\delta{n}_{max}}$ (%)"
   !GA(YLABEL)="${\\delta{n}}$ Kurtosis"
   !CALL FORTPLOT(GFC,rTIME,1.0*rNeSTAT(:,4),"NeSTAT",GA)
    CALL FORTPLOT(GFC,rTIME,100.0*rNeSTAT(:,5),"",GA,GX)
   END IF

  !GFC(PLTSFIGNUM)=312
  !GA(LINESIZE)="2"
  !GA(GRID)="NO"
  !GFC(PLTSFIGNUM)=211
   IF(iEySAT) THEN
   !GA(TITLE)="Pertubed Horizontal E-Field"
  ! GA(TITLE)=""
   !GA(XLABEL)="Time (sec)"

  ! GFC(PLTFIGNUM)=1
  ! GA(YLABEL)="$<\\delta{E_y}>$(mV/m)"
  !!CALL FORTPLOT(GFC,rTIME,1000.0*rEySTAT(:,1),"EySTAT",GA)
  ! CALL FORTPLOT(GFC,rTIME,1000.0*rEySTAT(:,1),"",GA)

  ! GFC(PLTFIGNUM)=2
  ! GA(TITLE)="Standard Deviation of Perturbed Horizontal Electric Field"
  ! GA(XLABEL)=""
  ! GA(YLABEL)="$\sigma_{\delta{E_y}}$(mV/m)"
  !IF(iNETCDF.EQ.LOOPEID) THEN
  ! GA(IPLOTLABEL)="$\sigma_{\delta{E_y},Max}$"
  ! GA(LINEPROP)="r-"
  !ELSEIF(iNETCDF.EQ.LOOPSID) THEN
  ! GA(IPLOTLABEL)="$\sigma_{\delta{E_y},Min}$"
  ! GA(LINEPROP)="b-"
  !END IF
  ! CALL FORTPLOT(GFC,rTIME,1000.0*rEySTAT(:,2),"",GA)
  !IF(iNETCDF.EQ.LOOPEID) THEN
  ! GA(IPLOTLABEL)="$\sigma_{\delta{E_y^+},Max}$"
  ! GA(LINEPROP)="r:"
  !ELSEIF(iNETCDF.EQ.LOOPSID) THEN
  ! GA(IPLOTLABEL)="$\sigma_{\delta{E_y^+},Min}$"
  ! GA(LINEPROP)="b:"
  !END IF
  ! CALL FORTPLOT(GFC,rTIME,1000.0*rEySTAT(:,3),"",GA)
  !IF(iNETCDF.EQ.LOOPEID) THEN
  ! GA(IPLOTLABEL)="$\sigma_{\delta{E_y^-},Max}$"
  ! GA(LINEPROP)="r--"
  !ELSEIF(iNETCDF.EQ.LOOPSID) THEN
  ! GA(IPLOTLABEL)="$\sigma_{\delta{E_y^-},Min}$"
  ! GA(LINEPROP)="b--"
  !END IF
  ! CALL FORTPLOT(GFC,rTIME,1000.0*ABS(rEySTAT(:,4)),"",GA)
   !CALL FORTPLOT(GFC,rTIME,1000.0*rEySTAT(:,2),"EySTAT",GA)

  ! GFC(PLTFIGNUM)=3
  ! GA(YLABEL)="${\\delta{E_y}}$"
  !!GA(YLABEL)="${\\delta{E_y}}$ Skewness"
  !!CALL FORTPLOT(GFC,rTIME,1.0*rEySTAT(:,3),"EySTAT",GA)
  ! CALL FORTPLOT(GFC,rTIME,1.0*rEySTAT(:,3),"",GA)

  ! GFC(PLTFIGNUM)=4
  ! GA(YLABEL)="${\\delta{E_y}}$"
  !!GA(YLABEL)="${\\delta{E_y}}$ Kurtosis"
  !!CALL FORTPLOT(GFC,rTIME,1.0*rEySTAT(:,4),"EySTAT",GA)
  ! CALL FORTPLOT(GFC,rTIME,1.0*rEySTAT(:,4),"",GA)

  ! GFC(PLTFIGNUM)=5
  ! GA(TITLE)="Maximum of Perturbed Horizontal Electric Field"
  ! GA(YLABEL)="${\delta{E_y},max}$(mV/m)"
   !CALL FORTPLOT(GFC,rTIME,1000.0*rEySTAT(:,2),"EySTAT",GA)
  ! CALL FORTPLOT(GFC,rTIME,1000.0*rEySTAT(:,5),"",GA)
   END IF

!  GFC(PLTSFIGNUM)=313
  !GFC(PLTSFIGNUM)=212
  !IF(iNETCDF.EQ.LOOPEID) THEN
  ! GFC(PLTFIGHOLD)=GETVALOF("NO")
  !ELSE
  ! GFC(PLTFIGHOLD)=GETVALOF("YES")
  !END IF
   IF(iEzSAT) THEN
   !GA(TITLE)="Pertubed Vertical E-Field"
  ! GA(TITLE)=""
  ! GA(XLABEL)="Time (sec)"

  ! GFC(PLTFIGNUM)=1
  ! GA(YLABEL)="$<\\delta{E_z}>$(mV/m)"
  !!CALL FORTPLOT(GFC,rTIME,1000.0*rEzSTAT(:,1),"EzSTAT",GA)
  ! CALL FORTPLOT(GFC,rTIME,1000.0*rEzSTAT(:,1),"",GA)

  ! GFC(PLTFIGNUM)=2
  ! GA(TITLE)=""
  ! GA(TITLE)="Standard Deviation of Perturbed Vertical Electric Field"
  ! GA(XLABEL)="Time (sec)"
  ! GA(YLABEL)="$\sigma_{\delta{E_z}}$(mV/m)"
  !IF(iNETCDF.EQ.LOOPEID) THEN
  ! GA(IPLOTLABEL)="$\sigma_{\delta{E_z},Max}$"
  ! GA(LINEPROP)="r-"
  !ELSEIF(iNETCDF.EQ.LOOPSID) THEN
  ! GA(IPLOTLABEL)="$\sigma_{\delta{E_z},Min}$"
  ! GA(LINEPROP)="b-"
  !END IF
  ! CALL FORTPLOT(GFC,rTIME,1000.0*rEzSTAT(:,2),"",GA)
  !IF(iNETCDF.EQ.LOOPEID) THEN
  ! GA(IPLOTLABEL)="$\sigma_{\delta{E_z^+},Max}$"
  ! GA(LINEPROP)="r:"
  !ELSEIF(iNETCDF.EQ.LOOPSID) THEN
  ! GA(IPLOTLABEL)="$\sigma_{\delta{E_z^+},Min}$"
  ! GA(LINEPROP)="b:"
  !END IF
  ! CALL FORTPLOT(GFC,rTIME,1000.0*rEzSTAT(:,3),"",GA)
  !IF(iNETCDF.EQ.LOOPEID) THEN
  ! GA(IPLOTLABEL)="$\sigma_{\delta{E_z^-},Max}$"
  ! GFC(PLTFIGHOLD)=GETVALOF("NO")
  ! GA(LINEPROP)="r--"
  !ELSEIF(iNETCDF.EQ.LOOPSID) THEN
  ! GA(IPLOTLABEL)="$\sigma_{\delta{E_z^-},Min}$"
  ! GFC(PLTFIGHOLD)=GETVALOF("YES")
  ! GA(LINEPROP)="b--"
  !END IF
  ! CALL FORTPLOT(GFC,rTIME,1000.0*ABS(rEzSTAT(:,4)),"",GA)
   !CALL FORTPLOT(GFC,rTIME,1000.0*rEzSTAT(:,2),"EzSTAT",GA)

  ! GFC(PLTFIGNUM)=3
  ! GA(YLABEL)="${\\delta{E_z}}$"
  !!GA(YLABEL)="${\\delta{E_z}}$ Skewness"
  !!CALL FORTPLOT(GFC,rTIME,1.0*rEzSTAT(:,3),"EzSTAT",GA)
  ! CALL FORTPLOT(GFC,rTIME,1.0*rEzSTAT(:,3),"",GA)

  ! GFC(PLTFIGNUM)=4
  ! GA(YLABEL)="${\\delta{E_z}}$"
  !!GA(YLABEL)="${\\delta{E_z}}$ Kurtosis"
  !!CALL FORTPLOT(GFC,rTIME,1.0*rEzSTAT(:,4),"EzSTAT",GA)
  ! CALL FORTPLOT(GFC,rTIME,1.0*rEzSTAT(:,4),"",GA)

  ! GFC(PLTFIGNUM)=5
  ! GA(TITLE)="Maximum of Perturbed Vertical Electric Field"
  ! GA(YLABEL)="${\delta{E_z}}$(mV/m)"
   !CALL FORTPLOT(GFC,rTIME,1000.0*rEzSTAT(:,2),"EzSTAT",GA)
  ! CALL FORTPLOT(GFC,rTIME,1000.0*rEzSTAT(:,5),"",GA)
   END IF
  END DO

  RETURN

  GA=""
 !GA(IPLOT)="contourf"
  GA(XLABEL)="Zonal"
  GA(YLABEL)="Vertical"
 !GA(CBAR)="yes"
 !GA(CBARTITLE)="$N_e\\upsilon_{ez}$"
  GFC(PLTSFIGNUM)=111
  GFC(PLTFIGHOLD)=GETVALOF("NO")
  DO iREC=1,nREC
   rNeSTAT(iREC,5)=SUM(rFLXnvey(1,:,:,iREC))
  END DO
  CALL FORTPLOT(GFC,1.0d0*rTIME,rNeSTAT(:,5),"FLXnvey",GA)
  DO iREC=1,nREC
   rNeSTAT(iREC,5)=SUM(rFLXnvez(1,:,:,iREC))
  END DO
  CALL FORTPLOT(GFC,1.0d0*rTIME,rNeSTAT(:,5),"FLXnvez",GA)
 !DO iREC=701,1001,10
 ! GFC(PLTFIGNUM)=iREC
 ! CALL FORTPLOT(GFC,1.0d0*Y3D(1,:,1),SUM(rFLXnvez(1,:,:,iREC),2),"FLXnvey",GA)
 ! CALL FORTPLOT(GFC,1.0d0*Z3D(1,1,:),SUM(rFLXnvey(1,:,:,iREC),1),"FLXnvez",GA)
 ! CALL FORTPLOT(GFC,1.0d0*Y3D(1,:,:),1.0d0*Z3D(1,:,:),rFLXnvez(1,:,:,iREC),"FLXnvez",GA)
 ! CALL FORTPLOT(GFC,1.0d0*Y3D(1,:,:),1.0d0*Z3D(1,:,:),rFLXnvey(1,:,:,iREC),"FLXnvey",GA)
 ! CALL FORTPLOT(GFC,1.0d0*Y3D(1,:,:),1.0d0*Z3D(1,:,:),rFLXnviz(1,:,:,iREC),"FLXnviz",GA)
 ! CALL FORTPLOT(GFC,1.0d0*Y3D(1,:,:),1.0d0*Z3D(1,:,:),rFLXnviy(1,:,:,iREC),"FLXnviy",GA)
 ! CALL FORTPLOT(GFC,1.0d0*Y3D(1,:,:),1.0d0*Z3D(1,:,:),rNe(1,:,:,iREC),"rNe",GA)
 ! CALL FORTPLOT(GFC,1.0d0*Y3D(1,:,:),1.0d0*Z3D(1,:,:),rEy(1,:,:,iREC),"rEy",GA)
 ! CALL FORTPLOT(GFC,1.0d0*Y3D(1,:,:),1.0d0*Z3D(1,:,:),rEz(1,:,:,iREC),"rEz",GA)
 !END DO


 !CALL FieldsPlot("yz",X3D,Y3D,Z3D,rfld,1,REAL(t_strt,4))
 END SUBROUTINE iDIAGNOSE


END MODULE iTEST

