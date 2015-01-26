MODULE SYSINFO

!**********************************************************************
!**************** SHARED VARAIBLES AND MATERICES **********************
!**********************************************************************

!**************************************************************
!                  GENERAL CONSTANTS
!**************************************************************
 REAL(4),PARAMETER :: PI=4.0d0*ATAN(1.0d0)
 REAL(4),PARAMETER :: KB=1.3807E-23
 REAL(4),PARAMETER :: Vc=2.9979E+8
 REAL(4),PARAMETER :: MU0=4.0E-7*PI
 REAL(4),PARAMETER :: EPSILON0=8.8542E-12
 REAL(4),PARAMETER :: N_Avogaros=6.023E23 

!**************************************************************
!                    EARTH RADIUS
!**************************************************************
 REAL(4),PARAMETER :: Re=6378.1

!**************************************************************
!                    SYSTEM STATUS
!**************************************************************
 LOGICAL                :: SYSRESTART        =.FALSE.
!**************************************************************
!         SYSTEM GRID COORDINATES AND PARAMETERS
!**************************************************************
!********************** GRIDING SYSTEM ************************
 LOGICAL                 :: RECTGRID         =.TRUE.
 LOGICAL                 :: SPHRGRID         =.FALSE.
 LOGICAL                 :: GEOGRAPHICGRID   =.TRUE.
 LOGICAL                 :: GEOMAGNETICGRID  =.FALSE.
!**************** PRIMATRY GRIDING SYSTEM *********************
 LOGICAL                 :: pRECTGRID         =.FALSE.
 LOGICAL                 :: pSPHRGRID         =.FALSE.
 LOGICAL                 :: pGEOGRAPHICGRID   =.FALSE.
 LOGICAL                 :: pGEOMAGNETICGRID  =.FALSE.
!**************** SECONDARY GRIDING SYSTEM ********************
 LOGICAL                 :: sRECTGRID         =.FALSE.
 LOGICAL                 :: sSPHRGRID         =.FALSE.
 LOGICAL                 :: sGEOGRAPHICGRID   =.FALSE.
 LOGICAL                 :: sGEOMAGNETICGRID  =.FALSE.
!********************** GRID SPACINGS *************************
 LOGICAL                 :: EQUISPACED       =.TRUE.
 LOGICAL                 :: XEQUISPACED      =.FALSE.
 LOGICAL                 :: YEQUISPACED      =.FALSE.
 LOGICAL                 :: ZEQUISPACED      =.FALSE.
 LOGICAL                 :: CHEBYSHEVSPACED  =.FALSE.
 LOGICAL                 :: XCHEBYSHEVSPACED =.FALSE.
 LOGICAL                 :: YCHEBYSHEVSPACED =.FALSE.
 LOGICAL                 :: ZCHEBYSHEVSPACED =.FALSE.
 LOGICAL                 :: LEGENDERSPACED   =.FALSE.
 LOGICAL                 :: XLEGENDERSPACED  =.FALSE.
 LOGICAL                 :: YLEGENDERSPACED  =.FALSE.
 LOGICAL                 :: ZLEGENDERSPACED  =.FALSE.
!************ RECTANGLE NUMBER OF R-GRID POINTS ***************
 INTEGER(4)              :: XDIM=1 
 INTEGER(4)              :: YDIM=1
 INTEGER(4)              :: ZDIM=1
!*************** R-SPACE RECTANGLE DIMENSIONS RANGES ******************
 REAL(8)                 :: Xsrt=   0.0
 REAL(8)                 :: Xend=   0.0
 REAL(8)                 :: Ysrt= -16.0*PI
 REAL(8)                 :: Yend=  16.0*PI
 REAL(8)                 :: Zsrt=   0.0*PI
 REAL(8)                 :: Zend=  32.0*PI
!*************** GEO DIMENSIONS RANGES ******************
 REAL(4)                 :: srtLAT= - 0.0
 REAL(4)                 :: endLAT=   0.0
 REAL(4)                 :: srtLON= -0.000050
 REAL(4)                 :: endLON=  0.000050
 REAL(4)                 :: srtALT=  105.000  
 REAL(4)                 :: endALT=  105.100
!*************** RECTANGLE R-SPACE ARRAYS **********************
 REAL(8),    ALLOCATABLE :: X1D(:),X2D(:,:),X3D(:,:,:)
 REAL(8),    ALLOCATABLE :: Y1D(:),Y2D(:,:),Y3D(:,:,:)
 REAL(8),    ALLOCATABLE :: Z1D(:),Z2D(:,:),Z3D(:,:,:)
!*************** GEO R-SPACE ARRAYS ********************
 REAL(4),    ALLOCATABLE :: LAT1D(:)
 REAL(4),    ALLOCATABLE :: LON1D(:)
 REAL(4),    ALLOCATABLE :: ALT1D(:)

!**************************************************************
!************** LOCAL TIME AND DATE RANGE MATRICES ************
!**************************************************************
!INTEGER(4)              :: srtLT=120000,endLT=130000,stpLT=10000
!INTEGER(4)              :: srtDT=19870312,endDT=19870312
 INTEGER(4)              :: srtLT=10000,endLT=20000,stpLT=10000
 INTEGER(4)              :: srtDT=20080320,endDT=20080320
 INTEGER(4), ALLOCATABLE :: LTHURMAT(:),LTMINMAT(:),LTSECMAT(:)
 INTEGER(4), ALLOCATABLE :: DTYYYYMAT(:),DTMMMAT(:),DTDDMAT(:)
 INTEGER(4), ALLOCATABLE :: DTMMDDMAT(:)
 INTEGER(4), ALLOCATABLE :: LTMAT(:),DTMAT(:)

!**************************************************************
!                 MISCELLIOUS VARIABLES
!**************************************************************

 CONTAINS

 SUBROUTINE WHATISTIME(NOWTIME)
  IMPLICIT NONE
 !character(8)  :: date
 !character(10) :: time
 !character(5)  :: zone
  integer       :: values(8)
  INTEGER(4)    :: NOWTIME(*)
  !using keyword arguments
  call date_and_time(VALUES=values)
  NOWTIME(1:3)=values(1:3)
  NOWTIME(4:7)=values(5:8)
 !PRINT*, NOWTIME(1)
 !call date_and_time(date,time,zone,values)
 !call date_and_time(DATE=date,ZONE=zone)
 !call date_and_time(TIME=time)
 !call date_and_time(VALUES=values)
 !print '(a,2x,a,2x,a)', date, time, zone
 !print '(8i5)', values
  RETURN
 END SUBROUTINE WHATISTIME

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

END MODULE SYSINFO
