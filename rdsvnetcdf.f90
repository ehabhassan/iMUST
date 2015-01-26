MODULE RDSVNETCDF

!**********************************************************************
! BSD SOFTWARE LICENSE
!************************
! PERTAINS TO FORTGRAPHER
!
! COPYRIGHT(c) 2011-2012 EHAB ALI HASSAN:
!*****************************************
! THIS CODE IS DEVELOPED BY: EHAB ALI HASSAN; EHAB@UTEXAS.COM
! RESEARCH ASSISTANT AT UNIVERSITY OF TEXAS AT AUSTIN, TEXAS, US.
! RESEARCH ASSISTANT AT AIN SHAMS UNIVERSITY, CAIRO, EGYPT.
!
! IT IS AN OPEN SOURCE CODE AND YOU ARE ALLOWED TO MAKE ALL REQUIRED CHANGES
! YOU MAY NEED TO THIS CODE BUT YOU HAVE TO SEND ME ALL YOU MODIFICATION AND
! MAKE THESE MODIFICATIONS OPEN FOR ALL OTHERS TO USE FOR FREE.
!
! YOU ARE NOT ALLOWED TO CHANGE THE NAME OF THIS MODULE OR ANY OF ITS
! FUNCTIONS, SUBROUTINES, VARIABLES, ETC.
!
! YOU ARE ALLOWED TO REDISTRIBUTE THIS CODE TO OTHERS WITH THE SAME 
! RESITRICTIONS MENTIONED ABOVE.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
! - Redistributions of source code must retain the above copyright
!   notice, this list of conditions and the following disclaimer. 
!  
! - Redistributions in binary form must reproduce the above copyright
!   notice, this list of conditions and the following disclaimer listed
!   in this license in the documentation and/or other materials
!   provided with the distribution.
!  
! - Neither the name of the copyright holders nor the names of its
!   contributors may be used to endorse or promote products derived from
!   this software without specific prior written permission.
!  
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
!
!
! THIS MODULE CONTAINS SOME SUBROUTINES FOR READING A NETCDF FILE.
! REQUIREMENT:
! (1) FORTRAN 90 COMPILER
! (2) NETCDF LIBRARY FROM http://www.unidata.ucar.edu/downloads/netcdf/index.jsp
! (3) NETCDF FORTRAN F90 FROM http://www.unidata.ucar.edu/downloads/netcdf/index.jsp
!
! SUBROUTINE NETCDF_VAR_NAMES(NETCDF_FILENAME):
! THIS SUBROUTINE RETRIEVES ALL THE VARIABLE NAMES AND THEIR DIMENSIONS ON THE NETCDF FILE
! THIS SUBROUTINE TAKES THE NAME OF THE NETCDF FILE AS AN INPUT AND PRINT THE NAME OF 
! ALL VARIABLES ON THE SCREEN AS AN OUTPUT
!
! SUBROUTINE NETCDF_VAR_NDIMS(NETCDF_FILENAME,VAR_NAME,VAR_NDIMS):
! THIS SUBROUTINE RETRIEVES THE SIZE OF A VARIABLE.
! THIS SUBROUTINE TAKES THE NAME OF A NETCDF FILE AND THE NAME OF A VARIABLE
! (CASE-SENSITIVE) THAT CAN BE OBTAINED FROM SUBROUTINE NETCDF_VAR_NAMES.
! THE OUTPUT IS AN INTEGER(4) THAT HAS THE SIZE OF VARIABLE ARRAY.
!
! SUBROUTINE NETCDF_VAR_DIMS(NETCDF_FILENAME,VAR_NAME,VAR_DIMS):
! THIS SUBROUTINE RETRIEVES THE LENGTH OF EVERY DIMENSION OF A VARIABLE.
! THIS SUBROUTINE TAKES THE NAME OF THE NETCDF FILE AND THE NAME OF A VARIABLE
! (CASE-SENSITIVE) THAT CAN BE OBTAINED FROM SUBROUTINE NETCDF_VAR_NAMES.
! THE OUTPUT IS AN INTEGER(4) ARRAY OF DIMENSION NAME_VAR_NDIMS THAT CAN BE
! OBTAINED FROM SUBROUTINE NETCDF_VAR_NDIMS(NETCDF_FILENAME,VAR_NAME,VAR_NDIMS).
! YOU CAN ALLOCATE THE VAR_DIMS ARRAY USING THE VAR_NDIMS OUTPUT.
!
! SUBROUTINE NETCDF_VAR_1D_GET(NETCDF_FILENAME,VAR_NAME,VAR_MAT)
! THE SUBROUTINE RETRIEVES THE DATA MATRIX OF A VARIABLE.
! THIS SUBROUTINE TAKES THE NAME OF THE NETCDF FILE AND THE NAME OF A VARIABLE
! (CASE-SENSITIVE) THAT CAN BE OBTAINED FROM SUBROUTINE NETCDF_VAR_NAMES.
! THE OUTPUT IS A REAL(4) ARRAY OF DIMENSION NAME_VAR_NDIMS THAT CAN BE
! OBTAINED FROM SUBROUTINE NETCDF_VAR_NDIMS(NETCDF_FILENAME,VAR_NAME,VAR_NDIMS).
! THE LENGTH OF EACH DIMENSION IS OBTAINED FROM CALLING THE
! SUBROUTINE NETCDF_VAR_DIMS(NETCDF_FILENAME,VAR_NAME,VAR_DIMS).
! THEN YOU CAN ALLOCATE THE VAR_MAT USING THE VALUES IN THE VAR_DIMS ARRAY.
! subroutine handle_err(status,VarID)
! THIS SUBROUTINE RETRIEVES THE ERRORS IN CALLING ANY NETCDF SUBROUTINE.
! SO IT CAN HANDLE ALL THE ERRORS.
!
! SUBROUTINE SAVEDIMS3D(NETCDFFNAME,X1D,Y1D,Z1D,DIMNAMES,DIMUNITS)
! SUBROUTINE SAVEDIMS2D(NETCDFFNAME,X1D,Y1D,DIMNAMES,DIMUNITS)
! SUBROUTINE SAVEDIMS1D(NETCDFFNAME,X1D,DIMNAMES,DIMUNITS)
! THESE SUBROUTINES ARE DESIGNED TO SAVE THE DIMENSIONS IN CASE OF
! 1D, 2D, OR 3D.
! THE TIME IS ALWAYS THEIR AND YOU NEED TO INCLUDE ITS NAME AND UNIT
! IN DIMNAMES AND DIMUNITS, RESPECTIVELY.
! DIMNAMES AND DIMUNITS ARE CHARACTER ARRAYS OF DIMENSION EQUALS THE
! SYSTEM DIMENSION+1. FOR EXAMPLE: IF THE SYSTEM DIMENSION IS 1D, THE 
! DIMENSION OF BOTH DIMNAMES AND DIMUNITS IS 2D.
! NETCDFFNAME IS A CHARACTER VARIABLE FOR THE NETCDF FILENAME.
! THE SUBROUTINES TAKE REAL(4) OR REAL(8) FOR X1D, Y1D, AND Z1D.
!
!

  USE NETCDF

 INTERFACE NETCDF_VAR_GET
  MODULE PROCEDURE NETCDF_VAR_1D_GET,NETCDF_VAR_2D_GET,NETCDF_VAR_3D_GET,NETCDF_VAR_4D_GET
 END INTERFACE NETCDF_VAR_GET

 INTERFACE SAVEDIMS
  MODULE PROCEDURE r4SAVEDIMS3D,r8SAVEDIMS3D, &
                   r4SAVEDIMS2D,r8SAVEDIMS2D, &
                   r4SAVEDIMS1D,r8SAVEDIMS1D
 END INTERFACE SAVEDIMS

 INTERFACE SAVEFLDS
  MODULE PROCEDURE r4SAVEFLDS4D,r8SAVEFLDS4D, &
                   r4SAVEFLDS3D,r8SAVEFLDS3D, &
                   r4SAVEFLDS2D,r8SAVEFLDS2D, &
                   r4SAVEFLDS1D,r8SAVEFLDS1D
 END INTERFACE SAVEFLDS

 INTERFACE READDIMS
  MODULE PROCEDURE r4READDIMS3D,r8READDIMS3D, &
                   r4READDIMS2D,r8READDIMS2D, &
                   r4READDIMS1D,r8READDIMS1D
 END INTERFACE READDIMS

 INTERFACE READFLDS
  MODULE PROCEDURE r4READFLDS4D,r8READFLDS4D, &
                   r4READFLDS3D,r8READFLDS3D, &
                   r4READFLDS2D,r8READFLDS2D, &
                   r4READFLDS1D,r8READFLDS1D
 END INTERFACE READFLDS

 INTERFACE NUM2STR
  MODULE PROCEDURE I4NUM2STR
 END INTERFACE NUM2STR

 CONTAINS

SUBROUTINE NETCDF_INFO(NETCDF_FILENAME,nDIMS,nVARS,nGLOBALATTS)
 IMPLICIT NONE
!***************** INPUT VARIABLES DECLARATION ***************
 INTEGER(4),INTENT(OUT),OPTIONAL   :: nDims, nVars, nGlobalAtts
 CHARACTER(*),INTENT(IN) :: NETCDF_FILENAME
!**************** LOCAL VARIABLES DECLARATION *****************
 INTEGER(4)   :: status, ncid
!INTEGER(4)   :: unlimDimID

 status = NF90_OPEN(TRIM(NETCDF_FILENAME), nf90_NoWrite, ncid )
 !status = nf90_inquire(ncid, nDims, nVars, nGlobalAtts, unlimdimid)
 IF(PRESENT(nGLOBALATTS)) THEN
  status = nf90_inquire(ncid, nDims, nVars, nGlobalAtts)
 ELSE IF(PRESENT(nVARS)) THEN
  status = nf90_inquire(ncid, nDims, nVars)
 ELSE IF(PRESENT(nDIMS)) THEN
  status = nf90_inquire(ncid, nDims)
 END IF
 status = NF90_CLOSE(ncid )

 RETURN
END SUBROUTINE NETCDF_INFO


SUBROUTINE NETCDF_VAR_NAMES(NETCDF_FILENAME)
 IMPLICIT NONE
!***************** INPUT VARIABLES DECLARATION ***************
 CHARACTER(*),INTENT(IN) :: NETCDF_FILENAME
!**************** LOCAL VARIABLES DECLARATION *****************
 CHARACTER(100):: varNAME
 INTEGER(4)   :: NETCDF_ID
 INTEGER(4)   :: nDims, nVars, nGlobalAtts
 INTEGER(4)   :: status, unlimDimID
 INTEGER(4)   :: I

 status = NF90_OPEN(TRIM(NETCDF_FILENAME), nf90_NoWrite, ncid = NETCDF_ID)
 status = nf90_inquire(NETCDF_ID, nDims, nVars, nGlobalAtts, unlimdimid)
 do I=1,nVars
  status = nf90_inquire_variable(ncid = NETCDF_ID, varID = I, name = varNAME, ndims = nDims)
  print*, I, ":", TRIM(varName), "(", char(48+nDims) ,")"
 end do
 status = NF90_CLOSE(ncid = NETCDF_ID)

 RETURN
END SUBROUTINE NETCDF_VAR_NAMES


SUBROUTINE NETCDF_VAR_NDIMS(NETCDF_FILENAME,VAR_NAME,VAR_NDIMS)
 IMPLICIT NONE
 CHARACTER(*) :: NETCDF_FILENAME,VAR_NAME
 INTEGER(4)   :: NETCDF_ID,VAR_ID,VAR_NDIMS
 INTEGER(4)   :: status

 status = NF90_OPEN(TRIM(NETCDF_FILENAME), nf90_NoWrite, ncid = NETCDF_ID)
 status = nf90_inq_varid(ncid = NETCDF_ID, name = VAR_NAME, varid = VAR_ID)
 status = nf90_inquire_variable(ncid = NETCDF_ID, varID = VAR_ID , ndims = VAR_NDIMS)
 status = NF90_CLOSE(ncid = NETCDF_ID)

 RETURN
END SUBROUTINE NETCDF_VAR_NDIMS


SUBROUTINE NETCDF_VAR_DIMS(NETCDF_FILENAME,VAR_NAME,VAR_DIMS)
 IMPLICIT NONE
 CHARACTER(*) :: NETCDF_FILENAME,VAR_NAME
 INTEGER(4)   :: dimIDs(nf90_max_var_dims)
 INTEGER(4)   :: NETCDF_ID,VAR_ID,VAR_DIMS(:)
 INTEGER(4)   :: status, nDims
 INTEGER(4)   :: I

 status = NF90_OPEN(TRIM(NETCDF_FILENAME), nf90_NoWrite, ncid = NETCDF_ID)
 status = nf90_inq_varid(ncid = NETCDF_ID, name = VAR_NAME, varid = VAR_ID)
 status = nf90_inquire_variable(ncid = NETCDF_ID, varID = VAR_ID , dimids = dimIDs)
 do I=1,SIZE(VAR_DIMS)
  status = nf90_inquire_dimension(ncid = NETCDF_ID, dimid = dimIDs(I), len = VAR_DIMS(I))
 end do
 status = NF90_CLOSE(ncid = NETCDF_ID)

 RETURN
END SUBROUTINE NETCDF_VAR_DIMS

SUBROUTINE NETCDF_VAR_1D_GET(NETCDF_FILENAME,VAR_NAME,VAR_MAT)
 IMPLICIT NONE
 CHARACTER(*) :: NETCDF_FILENAME,VAR_NAME
 INTEGER(4)   :: NETCDF_ID,VAR_ID
 INTEGER(4)   :: status
 REAL(4)      :: VAR_MAT(:)

 status = NF90_OPEN(TRIM(NETCDF_FILENAME), nf90_NoWrite, ncid = NETCDF_ID)
 status = nf90_inq_varid(ncid = NETCDF_ID, name = VAR_NAME, varid = VAR_ID)
 status = nf90_get_var(NETCDF_ID, VAR_ID, VAR_MAT)
 status = NF90_CLOSE(ncid = NETCDF_ID)

 RETURN
END SUBROUTINE NETCDF_VAR_1D_GET


SUBROUTINE NETCDF_VAR_2D_GET(NETCDF_FILENAME,VAR_NAME,VAR_MAT)
 IMPLICIT NONE
 CHARACTER(*) :: NETCDF_FILENAME,VAR_NAME
 INTEGER(4)   :: NETCDF_ID,VAR_ID
 INTEGER(4)   :: status
 REAL(4)      :: VAR_MAT(:,:)

 status = NF90_OPEN(TRIM(NETCDF_FILENAME), nf90_NoWrite, ncid = NETCDF_ID)
 status = nf90_inq_varid(ncid = NETCDF_ID, name = VAR_NAME, varid = VAR_ID)
 status = nf90_get_var(NETCDF_ID, VAR_ID, VAR_MAT)
 status = NF90_CLOSE(ncid = NETCDF_ID)

 RETURN
END SUBROUTINE NETCDF_VAR_2D_GET


SUBROUTINE NETCDF_VAR_3D_GET(NETCDF_FILENAME,VAR_NAME,VAR_MAT)
 IMPLICIT NONE
 CHARACTER(*) :: NETCDF_FILENAME,VAR_NAME
 INTEGER(4)   :: NETCDF_ID,VAR_ID
 INTEGER(4)   :: status
 REAL(4)      :: VAR_MAT(:,:,:)

 status = NF90_OPEN(TRIM(NETCDF_FILENAME), nf90_NoWrite, ncid = NETCDF_ID)
 status = nf90_inq_varid(ncid = NETCDF_ID, name = VAR_NAME, varid = VAR_ID)
 status = nf90_get_var(NETCDF_ID, VAR_ID, VAR_MAT)
 status = NF90_CLOSE(ncid = NETCDF_ID)

 RETURN
END SUBROUTINE NETCDF_VAR_3D_GET

SUBROUTINE NETCDF_VAR_4D_GET(NETCDF_FILENAME,VAR_NAME,VAR_MAT)
 IMPLICIT NONE
 CHARACTER(*) :: NETCDF_FILENAME,VAR_NAME
 INTEGER(4)   :: NETCDF_ID,VAR_ID
 INTEGER(4)   :: status
 REAL(4)      :: VAR_MAT(:,:,:,:)

 status = NF90_OPEN(TRIM(NETCDF_FILENAME), nf90_NoWrite, ncid = NETCDF_ID)
 status = nf90_inq_varid(ncid = NETCDF_ID, name = VAR_NAME, varid = VAR_ID)
 status = nf90_get_var(NETCDF_ID, VAR_ID, VAR_MAT)
 status = NF90_CLOSE(ncid = NETCDF_ID)

 RETURN
END SUBROUTINE NETCDF_VAR_4D_GET



!********************************************!
! SAVE THE FIELDS AND DIMENSIONS TO A NETCDF !
!********************************************!
!********************************************!
! DEFINE and SAVE DIMENSIONS TO A NEW NETCDF !
!********************************************!
!******!
!  3D  !
!******!
 SUBROUTINE r4SAVEDIMS3D(NETCDFFNAME,X1D,Y1D,Z1D,DIMNAMES,DIMUNITS)
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(4),      INTENT(IN) :: X1D(:),Y1D(:),Z1D(:)
  CHARACTER(*), INTENT(IN) :: NETCDFFNAME
  CHARACTER(*)             :: DIMNAMES(:)
  CHARACTER(*)             :: DIMUNITS(:)

  CALL SAVEDIMS(NETCDFFNAME,1.0d0*X1D,1.0d0*Y1D,1.0d0*Z1D,DIMNAMES,DIMUNITS)

  RETURN
 END SUBROUTINE r4SAVEDIMS3D

 SUBROUTINE r8SAVEDIMS3D(NETCDFFNAME,X1D,Y1D,Z1D,DIMNAMES,DIMUNITS)
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(8),      INTENT(IN) :: X1D(:),Y1D(:),Z1D(:)
  CHARACTER(*), INTENT(IN) :: NETCDFFNAME
  CHARACTER(*)             :: DIMNAMES(:)
  CHARACTER(*)             :: DIMUNITS(:)

 !NETCDF INPUTS VARIABLES
  INTEGER(4)                  :: NCID
  INTEGER(4)                  :: NRECS    = NF90_UNLIMITED
  CHARACTER(5)                :: UNITS    = "units"

 !CHARACTER(50),ALLOCATABLE :: FLDNAMES(:),FLDUNITS(:)
  INTEGER(4)                  :: X_dimid,Y_dimid,Z_dimid,T_dimid
  INTEGER(4)                  :: X_varid,Y_varid,Z_varid,T_varid


  !DIMENSIONS INFORMATION HANDLING
  IF(DIMNAMES(1).EQ."") DIMNAMES(1)="X"
  IF(DIMNAMES(2).EQ."") DIMNAMES(2)="Y"
  IF(DIMNAMES(3).EQ."") DIMNAMES(3)="Z"
  IF(DIMNAMES(4).EQ."") DIMNAMES(4)="T"
  IF(DIMUNITS(1).EQ."") DIMUNITS(1)="meter"
  IF(DIMUNITS(2).EQ."") DIMUNITS(2)="meter"
  IF(DIMUNITS(3).EQ."") DIMUNITS(3)="meter"
  IF(DIMUNITS(4).EQ."") DIMUNITS(4)="second"

  !CREATE THE NETCDF FILE
  call check( nf90_create(NETCDFFNAME, nf90_clobber, ncid) )

  !Define the dimensions.
  call check( nf90_def_dim(ncid, TRIM(DIMNAMES(1)), SIZE(X1D),     X_dimid) )
  call check( nf90_def_dim(ncid, TRIM(DIMNAMES(2)), SIZE(Y1D),     Y_dimid) )
  call check( nf90_def_dim(ncid, TRIM(DIMNAMES(3)), SIZE(Z1D),     Z_dimid) )
  !The record dimension is defined to have unlimited length - it can grow as needed.
  call check( nf90_def_dim(ncid, TRIM(DIMNAMES(4)), NF90_UNLIMITED, T_dimid) )

  !Define the coordinate variables.
  call check( nf90_def_var(ncid, TRIM(DIMNAMES(1)), NF90_REAL, X_dimid, X_varid) )
  call check( nf90_def_var(ncid, TRIM(DIMNAMES(2)), NF90_REAL, Y_dimid, Y_varid) )
  call check( nf90_def_var(ncid, TRIM(DIMNAMES(3)), NF90_REAL, Z_dimid, Z_varid) )
  call check( nf90_def_var(ncid, TRIM(DIMNAMES(4)), NF90_REAL, T_dimid, T_varid) )

  !Assign units attributes to coordinate variables.
  call check( nf90_put_att(ncid, X_varid, UNITS, TRIM(DIMUNITS(1))) )
  call check( nf90_put_att(ncid, Y_varid, UNITS, TRIM(DIMUNITS(2))) )
  call check( nf90_put_att(ncid, Z_varid, UNITS, TRIM(DIMUNITS(3))) )
  call check( nf90_put_att(ncid, T_varid, UNITS, TRIM(DIMUNITS(4))) )

 !End define mode.
  call check( nf90_enddef(ncid) )

 !Write the coordinate variable data.
  call check( nf90_put_var(ncid, X_varid, X1D) )
  call check( nf90_put_var(ncid, Y_varid, Y1D) )
  call check( nf90_put_var(ncid, Z_varid, Z1D) )

  RETURN
 END SUBROUTINE r8SAVEDIMS3D

!******!
!  2D  !
!******!
 SUBROUTINE r4SAVEDIMS2D(NETCDFFNAME,X1D,Y1D,DIMNAMES,DIMUNITS)
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(4),      INTENT(IN)    :: X1D(:),Y1D(:)
  CHARACTER(*), INTENT(IN)    :: NETCDFFNAME
  CHARACTER(*)                :: DIMNAMES(:)
  CHARACTER(*)                :: DIMUNITS(:)

  CALL SAVEDIMS(NETCDFFNAME,1.0d0*X1D,1.0d0*Y1D,DIMNAMES,DIMUNITS)

  RETURN
 END SUBROUTINE r4SAVEDIMS2D

 SUBROUTINE r8SAVEDIMS2D(NETCDFFNAME,X1D,Y1D,DIMNAMES,DIMUNITS)
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(8),      INTENT(IN)    :: X1D(:),Y1D(:)
  CHARACTER(*), INTENT(IN)    :: NETCDFFNAME
  CHARACTER(*)                :: DIMNAMES(:)
  CHARACTER(*)                :: DIMUNITS(:)

 !NETCDF INPUTS VARIABLES
  INTEGER(4)                  :: NCID
  INTEGER(4)                  :: NRECS    = NF90_UNLIMITED
  CHARACTER(5)                :: UNITS    = "units"

 !CHARACTER(50),ALLOCATABLE :: FLDNAMES(:),FLDUNITS(:)
  INTEGER(4)                  :: X_dimid,Y_dimid,T_dimid
  INTEGER(4)                  :: X_varid,Y_varid,T_varid


  !DIMENSIONS INFORMATION HANDLING
  IF(DIMNAMES(1).EQ."") DIMNAMES(1)="X"
  IF(DIMNAMES(2).EQ."") DIMNAMES(2)="Y"
  IF(DIMNAMES(3).EQ."") DIMNAMES(3)="T"
  IF(DIMUNITS(1).EQ."") DIMUNITS(1)="meter"
  IF(DIMUNITS(2).EQ."") DIMUNITS(2)="meter"
  IF(DIMUNITS(3).EQ."") DIMUNITS(3)="second"

  !CREATE THE NETCDF FILE
  call check( nf90_create(NETCDFFNAME, nf90_clobber, ncid) )

  !Define the dimensions.
  call check( nf90_def_dim(ncid, TRIM(DIMNAMES(1)), SIZE(X1D),      X_dimid) )
  call check( nf90_def_dim(ncid, TRIM(DIMNAMES(2)), SIZE(Y1D),      Y_dimid) )
  !The record dimension is defined to have unlimited length - it can grow as needed.
  call check( nf90_def_dim(ncid, TRIM(DIMNAMES(3)), NF90_UNLIMITED, T_dimid) )

  !Define the coordinate variables.
  call check( nf90_def_var(ncid, TRIM(DIMNAMES(1)), NF90_REAL, X_dimid, X_varid) )
  call check( nf90_def_var(ncid, TRIM(DIMNAMES(2)), NF90_REAL, Y_dimid, Y_varid) )
  call check( nf90_def_var(ncid, TRIM(DIMNAMES(3)), NF90_REAL, T_dimid, T_varid) )

  !Assign units attributes to coordinate variables.
  call check( nf90_put_att(ncid, X_varid, UNITS, TRIM(DIMUNITS(1))) )
  call check( nf90_put_att(ncid, Y_varid, UNITS, TRIM(DIMUNITS(2))) )
  call check( nf90_put_att(ncid, T_varid, UNITS, TRIM(DIMUNITS(3))) )

 !End define mode.
  call check( nf90_enddef(ncid) )

 !Write the coordinate variable data.
  call check( nf90_put_var(ncid, X_varid, X1D) )
  call check( nf90_put_var(ncid, Y_varid, Y1D) )

  RETURN
 END SUBROUTINE r8SAVEDIMS2D

!******!
!  1D  !
!******!
 SUBROUTINE r4SAVEDIMS1D(NETCDFFNAME,X1D,DIMNAMES,DIMUNITS)
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(4),      INTENT(IN)    :: X1D(:)
  CHARACTER(*), INTENT(IN)    :: NETCDFFNAME
  CHARACTER(*)                :: DIMNAMES(:)
  CHARACTER(*)                :: DIMUNITS(:)

  CALL SAVEDIMS(NETCDFFNAME,1.0d0*X1D,DIMNAMES,DIMUNITS)

  RETURN
 END SUBROUTINE r4SAVEDIMS1D

 SUBROUTINE r8SAVEDIMS1D(NETCDFFNAME,X1D,DIMNAMES,DIMUNITS)
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(8),      INTENT(IN)    :: X1D(:)
  CHARACTER(*), INTENT(IN)    :: NETCDFFNAME
  CHARACTER(*)                :: DIMNAMES(:)
  CHARACTER(*)                :: DIMUNITS(:)

 !NETCDF INPUTS VARIABLES
  INTEGER(4)                  :: NCID
  INTEGER(4)                  :: NRECS    = NF90_UNLIMITED
  CHARACTER(5)                :: UNITS    = "units"

 !CHARACTER(50),ALLOCATABLE :: FLDNAMES(:),FLDUNITS(:)
  INTEGER(4)                  :: X_dimid,T_dimid
  INTEGER(4)                  :: X_varid,T_varid


  !DIMENSIONS INFORMATION HANDLING
  IF(DIMNAMES(1).EQ."") DIMNAMES(1)="X"
  IF(DIMNAMES(2).EQ."") DIMNAMES(2)="T"
  IF(DIMUNITS(1).EQ."") DIMUNITS(1)="meter"
  IF(DIMUNITS(2).EQ."") DIMUNITS(2)="second"

  !CREATE THE NETCDF FILE
  call check( nf90_create(NETCDFFNAME, nf90_clobber, ncid) )

  !Define the dimensions.
  call check( nf90_def_dim(ncid, TRIM(DIMNAMES(1)), SIZE(X1D),      X_dimid) )
  !The record dimension is defined to have unlimited length - it can grow as needed.
  call check( nf90_def_dim(ncid, TRIM(DIMNAMES(2)), NF90_UNLIMITED, T_dimid) )

  !Define the coordinate variables.
  call check( nf90_def_var(ncid, TRIM(DIMNAMES(1)), NF90_REAL, X_dimid, X_varid) )
  call check( nf90_def_var(ncid, TRIM(DIMNAMES(2)), NF90_REAL, T_dimid, T_varid) )

  !Assign units attributes to coordinate variables.
  call check( nf90_put_att(ncid, X_varid, UNITS, TRIM(DIMUNITS(1))) )
  call check( nf90_put_att(ncid, T_varid, UNITS, TRIM(DIMUNITS(2))) )

 !End define mode.
  call check( nf90_enddef(ncid) )

 !Write the coordinate variable data.
  call check( nf90_put_var(ncid, X_varid, X1D) )

  RETURN
 END SUBROUTINE r8SAVEDIMS1D


!*******************************!
! DEFINE NEW FIELDS TO A NETCDF !
!*******************************!
 SUBROUTINE DEFNEWFLDS(NETCDFFNAME,NDIMS,FLDSNAME,FLDSUNIT)
  IMPLICIT NONE
  INTEGER(4),   INTENT(IN)    :: NDIMS
  CHARACTER(*), INTENT(IN)    :: NETCDFFNAME
  CHARACTER(*), INTENT(INOUT) :: FLDSNAME(:)
  CHARACTER(*), INTENT(INOUT) :: FLDSUNIT(:)

 !NETCDF INPUTS VARIABLES
  INTEGER(4)                  :: NCID
  CHARACTER(5)                :: UNITS    = "units"

 !CHARACTER(50),ALLOCATABLE :: FLDNAMES(:),FLDUNITS(:)
  INTEGER(4),ALLOCATABLE      :: dimids(:)
  INTEGER(4)                  :: FLD_varid

  !LOOP INDICES
  INTEGER(4)                  :: iostate
  INTEGER(4)                  :: nflds,ifld

  nflds=SIZE(FLDSNAME)

  CALL check(nf90_open(NETCDFFNAME, nf90_write, ncid))
  !start re-Define mode.
  call check(nf90_redef(ncid) )
  !allocate the dimensions ids
  ALLOCATE(dimids(NDIMS+1))

  IF(NDIMS.EQ.3) THEN
   dimids = (/ 1, 2, 3, 4 /)
  ELSE IF(NDIMS.EQ.2) THEN
   dimids = (/ 1, 2, 3 /)
  ELSE IF(NDIMS.EQ.1) THEN
   dimids = (/ 1, 2 /)
  END IF
  DO ifld=1,nflds
   !check if the field is already defined before
   iostate = nf90_inq_varid(ncid, TRIM(FLDSNAME(ifld)), FLD_varid)
   IF(iostate.eq.nf90_noerr) THEN
    PRINT*, "A FIELD WITH THE SAME NAME ("//TRIM(FLDSNAME(ifld))// ") IS ALREADY EXIST!"
    CONTINUE
   END IF
   !Define the netCdF variables.
   call check( nf90_def_var(ncid, TRIM(FLDSNAME(ifld)), NF90_DOUBLE, dimids, FLD_varid) )
   !Assign units attributes to the netCDF variables.
   call check( nf90_put_att(ncid, FLD_varid, UNITS, TRIM(FLDSUNIT(ifld))) )
  END DO
  !End define mode.
  call check( nf90_enddef(ncid) )

  RETURN
 END SUBROUTINE DEFNEWFLDS

 SUBROUTINE DEFNEWFLD(NETCDFFNAME,NDIMS,FLDNAME,FLDUNIT)
  IMPLICIT NONE
  INTEGER(4),   INTENT(IN)    :: NDIMS
  CHARACTER(*), INTENT(IN)    :: NETCDFFNAME
  CHARACTER(*), INTENT(INOUT) :: FLDNAME
  CHARACTER(*), INTENT(INOUT) :: FLDUNIT

 !NETCDF INPUTS VARIABLES
  INTEGER(4)                  :: NCID
  CHARACTER(5)                :: UNITS    = "units"

 !CHARACTER(50),ALLOCATABLE :: FLDNAMES(:),FLDUNITS(:)
  INTEGER(4),ALLOCATABLE      :: dimids(:)
  INTEGER(4)                  :: FLD_varid

  !LOOP INDICES
  INTEGER(4)                  :: iostate


  CALL check(nf90_open(NETCDFFNAME, nf90_write, ncid))
  !start re-Define mode.
  call check(nf90_redef(ncid) )
  !allocate the dimensions ids
  ALLOCATE(dimids(NDIMS+1))
  
  IF(NDIMS.EQ.3) THEN
   dimids = (/ 1, 2, 3, 4 /)
  ELSE IF(NDIMS.EQ.2) THEN
   dimids = (/ 1, 2, 3 /)
  ELSE IF(NDIMS.EQ.1) THEN
   dimids = (/ 1, 2 /)
  END IF
  !check if the field is already defined before
  iostate = nf90_inq_varid(ncid, TRIM(FLDNAME), FLD_varid)
  IF(iostate.eq.nf90_noerr) THEN
   PRINT*, "A FIELD WITH THE SAME NAME ("//TRIM(FLDNAME)// ") IS ALREADY EXIST!"
  ELSE
   !Define the netCdF variables.
   call check( nf90_def_var(ncid, TRIM(FLDNAME), NF90_DOUBLE, dimids, FLD_varid) )
   !Assign units attributes to the netCDF variables.
   call check( nf90_put_att(ncid, FLD_varid, UNITS, TRIM(FLDUNIT)) )
  END IF
  !End define mode.
  call check( nf90_enddef(ncid) )

  RETURN
 END SUBROUTINE DEFNEWFLD


!*************************!
! SAVE 3D MULTIPLE FIELDS !
!*************************!
 SUBROUTINE r4SAVEFLDS4D(NETCDFFNAME,RECID,RECTIME,rfld,DIMNAME,FLDNAME)
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(4),      INTENT(IN)    :: rfld(:,:,:,:)
  REAL(4),      INTENT(IN)    :: RECTIME
  INTEGER(4),   INTENT(IN)    :: RECID
  CHARACTER(*), INTENT(IN)    :: NETCDFFNAME
  CHARACTER(*), INTENT(INOUT) :: DIMNAME
  CHARACTER(*), INTENT(INOUT) :: FLDNAME

  CALL SAVEFLDS(NETCDFFNAME,RECID,RECTIME,1.0d0*rfld,DIMNAME,FLDNAME)

  RETURN
 END SUBROUTINE r4SAVEFLDS4D

 SUBROUTINE r8SAVEFLDS4D(NETCDFFNAME,RECID,RECTIME,rfld,DIMNAME,FLDNAME)
  USE NETCDF
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(8),      INTENT(IN)    :: rfld(:,:,:,:)
  REAL(4),      INTENT(IN)    :: RECTIME
  INTEGER(4),   INTENT(IN)    :: RECID
  CHARACTER(*), INTENT(IN)    :: NETCDFFNAME
  CHARACTER(*), INTENT(INOUT) :: DIMNAME
  CHARACTER(*), INTENT(INOUT) :: FLDNAME

 !NETCDF INPUTS VARIABLES
  INTEGER(4)                  :: NCID
  INTEGER(4), PARAMETER       :: NDIMS    = 4

 !CHARACTER(50),ALLOCATABLE :: FLDNAMES(:),FLDUNITS(:)
  INTEGER(4)                  :: T_varid
  INTEGER(4)                  :: start(NDIMS), count(NDIMS)
  INTEGER(4),ALLOCATABLE      :: FLDS_varid(:)

 !LOOP INDICES
  INTEGER(4)                  :: MATDIMS4D(4),ifld


  MATDIMS4D=SHAPE(rfld)

  CALL check(nf90_open(TRIM(NETCDFFNAME), nf90_write, ncid))
  !These settings tell netcdf to write one timestep of data.
  count = (/ MATDIMS4D(1), MATDIMS4D(2), MATDIMS4D(3), 1 /)
  start = (/ 1, 1, 1, RECID /)

  !record the time.
  call check( nf90_inq_varid(ncid, TRIM(DIMNAME), T_varid) )
  call check( nf90_put_var(ncid, T_varid, RECTIME, start = (/RECID/)) )
  !record the data.
  DO ifld=1,MATDIMS4D(4)
   call check( nf90_inq_varid(ncid, TRIM(FLDNAME), FLDS_varid(ifld)) )
   call check( nf90_put_var(ncid, FLDS_varid(ifld), rfld(:,:,:,ifld), start = start, count = count) )
  END DO
  
  !Close the file.
  !This causes netCDF to flush all buffers and make sure your data are really written to disk.
  call check( nf90_close(ncid) )

  RETURN
 END SUBROUTINE r8SAVEFLDS4D

!**********************!
! SAVE 3D SINGLE FIELD !
!**********************!
 SUBROUTINE r4SAVEFLDS3D(NETCDFFNAME,RECID,RECTIME,rfld,DIMNAME,FLDNAME)
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(4),      INTENT(IN)    :: rfld(:,:,:)
  REAL(4),      INTENT(IN)    :: RECTIME
  INTEGER(4),   INTENT(IN)    :: RECID
  CHARACTER(*), INTENT(IN)    :: NETCDFFNAME
  CHARACTER(*), INTENT(INOUT) :: DIMNAME
  CHARACTER(*), INTENT(INOUT) :: FLDNAME

  CALL SAVEFLDS(NETCDFFNAME,RECID,RECTIME,1.0d0*rfld,DIMNAME,FLDNAME)

  RETURN
 END SUBROUTINE r4SAVEFLDS3D

 SUBROUTINE r8SAVEFLDS3D(NETCDFFNAME,RECID,RECTIME,rfld,DIMNAME,FLDNAME)
  USE NETCDF
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(8),      INTENT(IN)    :: rfld(:,:,:)
  REAL(4),      INTENT(IN)    :: RECTIME
  INTEGER(4),   INTENT(IN)    :: RECID
  CHARACTER(*), INTENT(IN)    :: NETCDFFNAME
  CHARACTER(*), INTENT(INOUT) :: DIMNAME
  CHARACTER(*), INTENT(INOUT) :: FLDNAME

 !NETCDF INPUTS VARIABLES
  INTEGER(4)                  :: NCID
  INTEGER(4), PARAMETER       :: NDIMS    = 4

 !CHARACTER(50),ALLOCATABLE :: FLDNAMES(:),FLDUNITS(:)
  INTEGER(4)                  :: start(NDIMS), count(NDIMS)
  INTEGER(4)                  :: T_varid,FLD_varid

 !LOOP INDICES
  INTEGER(4)                  :: MATDIMS3D(3)


  MATDIMS3D=SHAPE(rfld)
  !open the netCDF file
  CALL check(nf90_open(NETCDFFNAME, nf90_write, ncid))
  !These settings tell netcdf to write one timestep of data.
  count = (/ MATDIMS3D(1), MATDIMS3D(2), MATDIMS3D(3), 1 /)
  start = (/ 1, 1, 1, RECID /)
  !record the time.
  call check( nf90_inq_varid(ncid, TRIM(DIMNAME), T_varid) )
  call check( nf90_put_var(ncid, T_varid, RECTIME, start = (/RECID/)) )
  !record the data.
  call check( nf90_inq_varid(ncid, TRIM(FLDNAME), FLD_varid) )
  call check( nf90_put_var(ncid, FLD_varid, rfld, start = start, count = count) )
  !Close the file.
  call check( nf90_close(ncid) )

  RETURN
 END SUBROUTINE r8SAVEFLDS3D

!**********************!
! SAVE 2D SINGLE FIELD !
!**********************!
 SUBROUTINE r4SAVEFLDS2D(NETCDFFNAME,RECID,RECTIME,rfld,DIMNAME,FLDNAME)
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(4),      INTENT(IN)    :: rfld(:,:)
  REAL(4),      INTENT(IN)    :: RECTIME
  INTEGER(4),   INTENT(IN)    :: RECID
  CHARACTER(*), INTENT(IN)    :: NETCDFFNAME
  CHARACTER(*), INTENT(INOUT) :: DIMNAME
  CHARACTER(*), INTENT(INOUT) :: FLDNAME

  CALL SAVEFLDS(NETCDFFNAME,RECID,RECTIME,1.0d0*rfld,DIMNAME,FLDNAME)

  RETURN
 END SUBROUTINE r4SAVEFLDS2D

 SUBROUTINE r8SAVEFLDS2D(NETCDFFNAME,RECID,RECTIME,rfld,DIMNAME,FLDNAME)
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(8),      INTENT(IN)    :: rfld(:,:)
  REAL(4),      INTENT(IN)    :: RECTIME
  INTEGER(4),   INTENT(IN)    :: RECID
  CHARACTER(*), INTENT(IN)    :: NETCDFFNAME
  CHARACTER(*), INTENT(INOUT) :: DIMNAME
  CHARACTER(*), INTENT(INOUT) :: FLDNAME

 !NETCDF INPUTS VARIABLES
  INTEGER(4)                  :: NCID
  INTEGER(4), PARAMETER       :: NDIMS    = 3

 !CHARACTER(50),ALLOCATABLE :: FLDNAMES(:),FLDUNITS(:)
  INTEGER(4)                  :: start(NDIMS), count(NDIMS)
  INTEGER(4)                  :: T_varid,FLD_varid

 !LOOP INDICES
  INTEGER(4)                  :: MATDIMS2D(2)


  MATDIMS2D=SHAPE(rfld)
  !open netCDF file
  CALL check(nf90_open(NETCDFFNAME, nf90_write, ncid))
  !These settings tell netcdf to write one timestep of data.
  count = (/ MATDIMS2D(1), MATDIMS2D(2), 1 /)
  start = (/ 1, 1, RECID /)
  !record the time.
  call check( nf90_inq_varid(ncid, TRIM(DIMNAME), T_varid) )
  call check( nf90_put_var(ncid, T_varid, RECTIME, start = (/RECID/)) )
  !record the data.
  call check( nf90_inq_varid(ncid, TRIM(FLDNAME), FLD_varid) )
  call check( nf90_put_var(ncid, FLD_varid, rfld, start = start, count = count) )
  !Close the file.
  call check( nf90_close(ncid) )

  RETURN
 END SUBROUTINE r8SAVEFLDS2D

!**********************!
! SAVE 1D SINGLE FIELD !
!**********************!
 SUBROUTINE r4SAVEFLDS1D(NETCDFFNAME,RECID,RECTIME,rfld,DIMNAME,FLDNAME)
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(4),      INTENT(IN)    :: rfld(:)
  REAL(4),      INTENT(IN)    :: RECTIME
  INTEGER(4),   INTENT(IN)    :: RECID
  CHARACTER(*), INTENT(IN)    :: NETCDFFNAME
  CHARACTER(*), INTENT(INOUT) :: DIMNAME
  CHARACTER(*), INTENT(INOUT) :: FLDNAME

  CALL SAVEFLDS(NETCDFFNAME,RECID,RECTIME,1.0d0*rfld,DIMNAME,FLDNAME)

  RETURN
 END SUBROUTINE r4SAVEFLDS1D

 SUBROUTINE r8SAVEFLDS1D(NETCDFFNAME,RECID,RECTIME,rfld,DIMNAME,FLDNAME)
  USE NETCDF
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(8),      INTENT(IN)    :: rfld(:)
  REAL(4),      INTENT(IN)    :: RECTIME
  INTEGER(4),   INTENT(IN)    :: RECID
  CHARACTER(*), INTENT(IN)    :: NETCDFFNAME
  CHARACTER(*), INTENT(INOUT) :: DIMNAME
  CHARACTER(*), INTENT(INOUT) :: FLDNAME

 !NETCDF INPUTS VARIABLES
  INTEGER(4)                  :: NCID
  INTEGER(4), PARAMETER       :: NDIMS    = 2

 !CHARACTER(50),ALLOCATABLE :: FLDNAMES(:),FLDUNITS(:)
  INTEGER(4)                  :: start(NDIMS), count(NDIMS)
  INTEGER(4)                  :: T_varid,FLD_varid

 !LOOP INDICES
  INTEGER(4)                  :: MATDIMS1D(1)


  MATDIMS1D=SIZE(rfld)
  !open the netCDF file
  CALL check(nf90_open(TRIM(NETCDFFNAME), nf90_write, ncid))
  !These settings tell netcdf to write one timestep of data.
  count = (/ MATDIMS1D(1), 1 /)
  start = (/ 1, RECID /)

  !record the time.
  call check( nf90_inq_varid(ncid, TRIM(DIMNAME), T_varid) )
  call check( nf90_put_var(ncid, T_varid, RECTIME, start = (/RECID/)) )
  !record the data.
  call check( nf90_inq_varid(ncid, TRIM(FLDNAME), FLD_varid) )
  call check( nf90_put_var(ncid, FLD_varid, rfld, start = start, count = count) )
  !Close the file.
  call check( nf90_close(ncid) )

  RETURN
 END SUBROUTINE r8SAVEFLDS1D

!***************************************************!
! READ THE FIELDS AND DIMENSIONS FROM A NETCDF FILE !
!***************************************************!
!*****************!
! READ DIMENSIONS !
!*****************!
!******!
!  3D  !
!******!
 SUBROUTINE r4READDIMS3D(NETCDFFNAME,X1D,Y1D,Z1D,T1D,DIMNAMES,DIMUNITS)
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(4),      INTENT(OUT) :: X1D(:),Y1D(:),Z1D(:)
  REAL(4),      INTENT(OUT) :: T1D(:)
  CHARACTER(*), INTENT(OUT) :: DIMNAMES(:)
  CHARACTER(*), INTENT(OUT) :: DIMUNITS(:)
  CHARACTER(*), INTENT(IN)  :: NETCDFFNAME

 !LOCAL VARIABLES
  REAL(8),      ALLOCATABLE :: LX1D(:),LY1D(:),LZ1D(:)
  INTEGER(4)                :: iXDIM,iYDIM,iZDIM


  IF(ALLOCATED(LX1D).EQV..FALSE.) ALLOCATE(LX1D(SIZE(X1D)))
  IF(ALLOCATED(LY1D).EQV..FALSE.) ALLOCATE(LY1D(SIZE(Y1D)))
  IF(ALLOCATED(LZ1D).EQV..FALSE.) ALLOCATE(LZ1D(SIZE(Z1D)))
  CALL READDIMS(NETCDFFNAME,LX1D,LY1D,LZ1D,T1D,DIMNAMES,DIMUNITS)
  DO iXDIM=1,SIZE(X1D)
   X1D(iXDIM)=REAL(LX1D(iXDIM),4)
  END DO
  DO iYDIM=1,SIZE(Y1D)
   Y1D(iYDIM)=REAL(LY1D(iYDIM),4)
  END DO
  DO iZDIM=1,SIZE(Z1D)
   Z1D(iZDIM)=REAL(LZ1D(iZDIM),4)
  END DO

  RETURN
 END SUBROUTINE r4READDIMS3D

 SUBROUTINE r8READDIMS3D(NETCDFFNAME,X1D,Y1D,Z1D,T1D,DIMNAMES,DIMUNITS)
  USE NETCDF
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(8),      INTENT(OUT) :: X1D(:),Y1D(:),Z1D(:)
  REAL(4),      INTENT(OUT) :: T1D(:)
  CHARACTER(*), INTENT(OUT) :: DIMUNITS(:)
  CHARACTER(*), INTENT(IN)  :: DIMNAMES(:)
  CHARACTER(*), INTENT(IN)  :: NETCDFFNAME

 !NETCDF INPUTS VARIABLES
  INTEGER(4)                :: NCID
  INTEGER(4)                :: X_varid,Y_varid,Z_varid,T_varid
  CHARACTER(5)              :: UNITS   = "units"


  !Open the file. 
  call check(nf90_open(NETCDFFNAME, NF90_nowrite, ncid))

  !Get the varids of the latitude and longitude coordinate variables.
  call check( nf90_inq_varid(ncid, TRIM(DIMNAMES(1)), X_varid) )
  call check( nf90_inq_varid(ncid, TRIM(DIMNAMES(2)), Y_varid) )
  call check( nf90_inq_varid(ncid, TRIM(DIMNAMES(3)), Z_varid) )
  call check( nf90_inq_varid(ncid, TRIM(DIMNAMES(4)), T_varid) )

  !Read the latitude and longitude units.
  call check( nf90_get_att(ncid, X_varid, UNITS, DIMUNITS(1)) )
  call check( nf90_get_att(ncid, Y_varid, UNITS, DIMUNITS(2)) )
  call check( nf90_get_att(ncid, Z_varid, UNITS, DIMUNITS(3)) )
  call check( nf90_get_att(ncid, T_varid, UNITS, DIMUNITS(4)) )

  !Read the latitude and longitude data.
  call check( nf90_get_var(ncid, X_varid, X1D) )
  call check( nf90_get_var(ncid, Y_varid, Y1D) )
  call check( nf90_get_var(ncid, Z_varid, Z1D) )
  call check( nf90_get_var(ncid, T_varid, T1D) )
 
  !Close the file.
  !This causes netCDF to flush all buffers and make sure your data are really written to disk.
  call check( nf90_close(ncid) )

  RETURN
 END SUBROUTINE r8READDIMS3D

!******!
!  2D  !
!******!
 SUBROUTINE r4READDIMS2D(NETCDFFNAME,X1D,Y1D,T1D,DIMNAMES,DIMUNITS)
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(4),      INTENT(OUT) :: X1D(:),Y1D(:)
  REAL(4),      INTENT(OUT) :: T1D(:)
  CHARACTER(*), INTENT(OUT) :: DIMNAMES(:)
  CHARACTER(*), INTENT(OUT) :: DIMUNITS(:)
  CHARACTER(*), INTENT(IN)  :: NETCDFFNAME

 !LOCAL VARIABLES
  REAL(8),      ALLOCATABLE :: LX1D(:),LY1D(:)
  INTEGER(4)                :: iXDIM,iYDIM


  IF(ALLOCATED(LX1D).EQV..FALSE.) ALLOCATE(LX1D(SIZE(X1D)))
  IF(ALLOCATED(LY1D).EQV..FALSE.) ALLOCATE(LY1D(SIZE(Y1D)))
  CALL READDIMS(NETCDFFNAME,LX1D,LY1D,T1D,DIMNAMES,DIMUNITS)
  DO iXDIM=1,SIZE(X1D)
   X1D(iXDIM)=REAL(LX1D(iXDIM),4)
  END DO
  DO iYDIM=1,SIZE(Y1D)
   Y1D(iYDIM)=REAL(LY1D(iYDIM),4)
  END DO

  RETURN
 END SUBROUTINE r4READDIMS2D

 SUBROUTINE r8READDIMS2D(NETCDFFNAME,X1D,Y1D,T1D,DIMNAMES,DIMUNITS)
  USE NETCDF
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(8),      INTENT(OUT) :: X1D(:),Y1D(:)
  REAL(4),      INTENT(OUT) :: T1D(:)
  CHARACTER(*), INTENT(OUT) :: DIMUNITS(:)
  CHARACTER(*), INTENT(IN)  :: DIMNAMES(:)
  CHARACTER(*), INTENT(IN)  :: NETCDFFNAME

 !NETCDF INPUTS VARIABLES
  INTEGER(4)                :: NCID
  INTEGER(4)                :: X_varid,Y_varid,Z_varid,T_varid
  CHARACTER(5)              :: UNITS   = "units"


  !Open the file. 
  call check(nf90_open(NETCDFFNAME, NF90_nowrite, ncid))

  !Get the varids of the latitude and longitude coordinate variables.
  call check( nf90_inq_varid(ncid, TRIM(DIMNAMES(1)), X_varid) )
  call check( nf90_inq_varid(ncid, TRIM(DIMNAMES(2)), Y_varid) )
  call check( nf90_inq_varid(ncid, TRIM(DIMNAMES(3)), T_varid) )

  !Read the latitude and longitude units.
  call check( nf90_get_att(ncid, X_varid, UNITS, DIMUNITS(1)) )
  call check( nf90_get_att(ncid, Y_varid, UNITS, DIMUNITS(2)) )
  call check( nf90_get_att(ncid, T_varid, UNITS, DIMUNITS(3)) )

  !Read the latitude and longitude data.
  call check( nf90_get_var(ncid, X_varid, X1D) )
  call check( nf90_get_var(ncid, Y_varid, Y1D) )
  call check( nf90_get_var(ncid, T_varid, T1D) )
 
  !Close the file.
  !This causes netCDF to flush all buffers and make sure your data are really written to disk.
  call check( nf90_close(ncid) )

  RETURN
 END SUBROUTINE r8READDIMS2D

!******!
!  1D  !
!******!
 SUBROUTINE r4READDIMS1D(NETCDFFNAME,X1D,T1D,DIMNAMES,DIMUNITS)
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(4),      INTENT(OUT) :: X1D(:)
  REAL(4),      INTENT(OUT) :: T1D(:)
  CHARACTER(*), INTENT(OUT) :: DIMNAMES(:)
  CHARACTER(*), INTENT(OUT) :: DIMUNITS(:)
  CHARACTER(*), INTENT(IN)  :: NETCDFFNAME

 !LOCAL VARIABLES
  REAL(8),      ALLOCATABLE :: LX1D(:)
  INTEGER(4)                :: iXDIM


  IF(ALLOCATED(LX1D).EQV..FALSE.) ALLOCATE(LX1D(SIZE(X1D)))
  CALL READDIMS(NETCDFFNAME,LX1D,T1D,DIMNAMES,DIMUNITS)
  DO iXDIM=1,SIZE(X1D)
   X1D(iXDIM)=REAL(LX1D(iXDIM),4)
  END DO

  RETURN
 END SUBROUTINE r4READDIMS1D

 SUBROUTINE r8READDIMS1D(NETCDFFNAME,X1D,T1D,DIMNAMES,DIMUNITS)
  USE NETCDF
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(8),      INTENT(OUT) :: X1D(:)
  REAL(4),      INTENT(OUT) :: T1D(:)
  CHARACTER(*), INTENT(OUT) :: DIMUNITS(:)
  CHARACTER(*), INTENT(IN)  :: DIMNAMES(:)
  CHARACTER(*), INTENT(IN)  :: NETCDFFNAME

 !NETCDF INPUTS VARIABLES
  INTEGER(4)                :: NCID
  INTEGER(4)                :: X_varid,Y_varid,Z_varid,T_varid
  CHARACTER(5)              :: UNITS   = "units"


  !Open the file. 
  call check(nf90_open(NETCDFFNAME, NF90_nowrite, ncid))

  !Get the varids of the latitude and longitude coordinate variables.
  call check( nf90_inq_varid(ncid, TRIM(DIMNAMES(1)), X_varid) )
  call check( nf90_inq_varid(ncid, TRIM(DIMNAMES(2)), T_varid) )

  !Read the latitude and longitude units.
  call check( nf90_get_att(ncid, X_varid, UNITS, DIMUNITS(1)) )
  call check( nf90_get_att(ncid, T_varid, UNITS, DIMUNITS(2)) )

  !Read the latitude and longitude data.
  call check( nf90_get_var(ncid, X_varid, X1D) )
  call check( nf90_get_var(ncid, T_varid, T1D) )
 
  !Close the file.
  !This causes netCDF to flush all buffers and make sure your data are really written to disk.
  call check( nf90_close(ncid) )

  RETURN
 END SUBROUTINE r8READDIMS1D


!*************************!
! READ 3D MULTIPLE FIELDS !
!*************************!
 SUBROUTINE r4READFLDS4D(NETCDFFNAME,RECID,rfld,FLDNAMES,FLDUNITS)
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(4),      INTENT(OUT) :: rfld(:,:,:,:)
  INTEGER(4),   INTENT(IN)  :: RECID
  CHARACTER(*), INTENT(IN)  :: NETCDFFNAME
  CHARACTER(*), INTENT(IN)  :: FLDNAMES(:)
  CHARACTER(*), INTENT(OUT) :: FLDUNITS(:)

  REAL(8),      ALLOCATABLE   :: Lrfld(:,:,:,:)
  INTEGER(4)                  :: ifld,iXDIM,iYDIM,iZDIM
  INTEGER(4)                  :: MATDIMS4D(4)

  MATDIMS4D=SHAPE(rfld)
  if(ALLOCATED(Lrfld).EQV..FALSE.) ALLOCATE(Lrfld(MATDIMS4D(1),MATDIMS4D(2),MATDIMS4D(3),MATDIMS4D(4)))
  CALL READFLDS(NETCDFFNAME,RECID,Lrfld,FLDNAMES,FLDUNITS)
  DO iXDIM=1,MATDIMS4D(1)
   DO iYDIM=1,MATDIMS4D(2)
    DO iZDIM=1,MATDIMS4D(3)
     DO ifld=1,MATDIMS4D(4)
       rfld(iXDIM,iYDIM,iZDIM,ifld)=REAL(Lrfld(iXDIM,iYDIM,iZDIM,ifld),4)
     END DO
    END DO
   END DO
  END DO

  RETURN
 END SUBROUTINE r4READFLDS4D

 SUBROUTINE r8READFLDS4D(NETCDFFNAME,RECID,rfld,FLDNAMES,FLDUNITS)
  USE NETCDF
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(8),      INTENT(OUT) :: rfld(:,:,:,:)
  INTEGER(4),   INTENT(IN)  :: RECID
  CHARACTER(*), INTENT(IN)  :: NETCDFFNAME
  CHARACTER(*), INTENT(IN)  :: FLDNAMES(:)
  CHARACTER(*), INTENT(OUT) :: FLDUNITS(:)

 !INPUT VARIABLES FACTS
  INTEGER(4)                :: MATDIMS4D(4)

 !NETCDF INPUTS VARIABLES
  INTEGER(4)                :: NCID
  INTEGER(4), PARAMETER     :: NDIMS = 4
  INTEGER(4)                :: start(NDIMS),count(NDIMS)
  INTEGER(4),ALLOCATABLE    :: FLDS_varid(:)
  CHARACTER(5)              :: UNITS   = "units"

 !LOOP INDICES
  INTEGER(4)                :: ifld


  MATDIMS4D=SHAPE(rfld)

  !Open the file. 
  call check(nf90_open(NETCDFFNAME, NF90_nowrite, ncid))

  !These settings tell netcdf to write one timestep of data.
  !(The setting of start(4) inside the loop below tells netCDF which timestep to write.)
  count = (/ MATDIMS4D(1), MATDIMS4D(2), MATDIMS4D(3), 1 /)
  start = (/ 1, 1, 1, RECID /)

  !Write the pretend data.
  !The arrays only hold one timestep worth of data. We will just rewrite the same data for each timestep.
  DO ifld=1,MATDIMS4D(4)
   call check( nf90_inq_varid(ncid, TRIM(FLDNAMES(ifld)), FLDS_varid(ifld)) )
   call check( nf90_get_att(ncid, FLDS_varid(ifld), UNITS, FLDUNITS(ifld) ) )
   call check( nf90_get_var(ncid, FLDS_varid(ifld), rfld(:,:,:,ifld), start = start, count = count) )
  END DO
  
  !Close the file.
  call check( nf90_close(ncid) )

  RETURN
 END SUBROUTINE r8READFLDS4D


!**********************!
! READ 3D SINGLE FIELD !
!**********************!
 SUBROUTINE r4READFLDS3D(NETCDFFNAME,RECID,rfld,FLDNAME,FLDUNIT)
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(4),      INTENT(OUT) :: rfld(:,:,:)
  INTEGER(4),   INTENT(IN)  :: RECID
  CHARACTER(*), INTENT(IN)  :: NETCDFFNAME
  CHARACTER(*), INTENT(IN)  :: FLDNAME
  CHARACTER(*), INTENT(OUT) :: FLDUNIT

  REAL(8),      ALLOCATABLE   :: Lrfld(:,:,:)
  INTEGER(4)                  :: iXDIM,iYDIM,iZDIM
  INTEGER(4)                  :: MATDIMS3D(3)

  MATDIMS3D=SHAPE(rfld)
  if(ALLOCATED(Lrfld).EQV..FALSE.) ALLOCATE(Lrfld(MATDIMS3D(1),MATDIMS3D(2),MATDIMS3D(3)))
  CALL READFLDS(NETCDFFNAME,RECID,Lrfld,FLDNAME,FLDUNIT)
  DO iXDIM=1,MATDIMS3D(1)
   DO iYDIM=1,MATDIMS3D(2)
    DO iZDIM=1,MATDIMS3D(3)
     rfld(iXDIM,iYDIM,iZDIM)=REAL(Lrfld(iXDIM,iYDIM,iZDIM),4)
    END DO
   END DO
  END DO

  RETURN
 END SUBROUTINE r4READFLDS3D

 SUBROUTINE r8READFLDS3D(NETCDFFNAME,RECID,rfld,FLDNAME,FLDUNIT)
  USE NETCDF
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(8),      INTENT(OUT) :: rfld(:,:,:)
  INTEGER(4),   INTENT(IN)  :: RECID
  CHARACTER(*), INTENT(IN)  :: NETCDFFNAME
  CHARACTER(*), INTENT(IN)  :: FLDNAME
  CHARACTER(*), INTENT(OUT) :: FLDUNIT

 !INPUT VARIABLES FACTS
  INTEGER(4)                :: MATDIMS3D(3)

 !NETCDF INPUTS VARIABLES
  INTEGER(4)                :: NCID
  INTEGER(4), PARAMETER     :: NDIMS = 4
  INTEGER(4)                :: start(NDIMS),count(NDIMS)
  INTEGER(4)                :: FLD_varid
  CHARACTER(5)              :: UNITS   = "units"

 !LOOP INDICES
  INTEGER(4)                :: ifld


  MATDIMS3D=SHAPE(rfld)

  !Open the file. 
  call check(nf90_open(NETCDFFNAME, NF90_nowrite, ncid))

  !These settings tell netcdf to write one timestep of data.
  !(The setting of start(4) inside the loop below tells netCDF which timestep to write.)
  count = (/ MATDIMS3D(1), MATDIMS3D(2), MATDIMS3D(3), 1 /)
  start = (/ 1, 1, 1, RECID /)

  !Write the pretend data.
  !The arrays only hold one timestep worth of data. We will just rewrite the same data for each timestep.
  call check( nf90_inq_varid(ncid, TRIM(FLDNAME), FLD_varid) )
  call check( nf90_get_att(ncid, FLD_varid, UNITS, FLDUNIT ) )
  call check( nf90_get_var(ncid, FLD_varid, rfld, start = start, count = count) )
  
  !Close the file.
  call check( nf90_close(ncid) )

  RETURN
 END SUBROUTINE r8READFLDS3D


!**********************!
! READ 2D SINGLE FIELD !
!**********************!
 SUBROUTINE r4READFLDS2D(NETCDFFNAME,RECID,rfld,FLDNAME,FLDUNIT)
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(4),      INTENT(OUT) :: rfld(:,:)
  INTEGER(4),   INTENT(IN)  :: RECID
  CHARACTER(*), INTENT(IN)  :: NETCDFFNAME
  CHARACTER(*), INTENT(IN)  :: FLDNAME
  CHARACTER(*), INTENT(OUT) :: FLDUNIT

  REAL(8),      ALLOCATABLE   :: Lrfld(:,:)
  INTEGER(4)                  :: iXDIM,iYDIM
  INTEGER(4)                  :: MATDIMS2D(2)

  MATDIMS2D=SHAPE(rfld)
  if(ALLOCATED(Lrfld).EQV..FALSE.) ALLOCATE(Lrfld(MATDIMS2D(1),MATDIMS2D(2)))
  CALL READFLDS(NETCDFFNAME,RECID,Lrfld,FLDNAME,FLDUNIT)
  DO iXDIM=1,MATDIMS2D(1)
   DO iYDIM=1,MATDIMS2D(2)
     rfld(iXDIM,iYDIM)=REAL(Lrfld(iXDIM,iYDIM),4)
   END DO
  END DO

  RETURN
 END SUBROUTINE r4READFLDS2D

 SUBROUTINE r8READFLDS2D(NETCDFFNAME,RECID,rfld,FLDNAME,FLDUNIT)
  USE NETCDF
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(8),      INTENT(OUT) :: rfld(:,:)
  INTEGER(4),   INTENT(IN)  :: RECID
  CHARACTER(*), INTENT(IN)  :: NETCDFFNAME
  CHARACTER(*), INTENT(IN)  :: FLDNAME
  CHARACTER(*), INTENT(OUT) :: FLDUNIT

 !INPUT VARIABLES FACTS
  INTEGER(4)                :: MATDIMS2D(2)

 !NETCDF INPUTS VARIABLES
  INTEGER(4)                :: NCID
  INTEGER(4), PARAMETER     :: NDIMS = 3
  INTEGER(4)                :: start(NDIMS),count(NDIMS)
  INTEGER(4)                :: FLD_varid
  CHARACTER(5)              :: UNITS   = "units"

 !LOOP INDICES
  INTEGER(4)                :: ifld


  MATDIMS2D=SHAPE(rfld)

  !Open the file. 
  call check(nf90_open(NETCDFFNAME, NF90_nowrite, ncid))

  !These settings tell netcdf to write one timestep of data.
  !(The setting of start(4) inside the loop below tells netCDF which timestep to write.)
  count = (/ MATDIMS2D(1), MATDIMS2D(2), 1 /)
  start = (/ 1, 1, RECID /)

  !Write the pretend data.
  !The arrays only hold one timestep worth of data. We will just rewrite the same data for each timestep.
  call check( nf90_inq_varid(ncid, TRIM(FLDNAME), FLD_varid) )
  call check( nf90_get_att(ncid, FLD_varid, UNITS, FLDUNIT ) )
  call check( nf90_get_var(ncid, FLD_varid, rfld, start = start, count = count) )
  
  !Close the file.
  call check( nf90_close(ncid) )

  RETURN
 END SUBROUTINE r8READFLDS2D


!**********************!
! READ 1D SINGLE FIELD !
!**********************!
 SUBROUTINE r4READFLDS1D(NETCDFFNAME,RECID,rfld,FLDNAME,FLDUNIT)
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(4),      INTENT(OUT) :: rfld(:)
  INTEGER(4),   INTENT(IN)  :: RECID
  CHARACTER(*), INTENT(IN)  :: NETCDFFNAME
  CHARACTER(*), INTENT(IN)  :: FLDNAME
  CHARACTER(*), INTENT(OUT) :: FLDUNIT

  REAL(8),      ALLOCATABLE   :: Lrfld(:)
  INTEGER(4)                  :: iXDIM
  INTEGER(4)                  :: MATDIMS1D(1)

  MATDIMS1D=SIZE(rfld)
  if(ALLOCATED(Lrfld).EQV..FALSE.) ALLOCATE(Lrfld(MATDIMS1D(1)))
  CALL READFLDS(NETCDFFNAME,RECID,Lrfld,FLDNAME,FLDUNIT)
  DO iXDIM=1,MATDIMS1D(1)
   rfld(iXDIM)=REAL(Lrfld(iXDIM),4)
  END DO

  RETURN
 END SUBROUTINE r4READFLDS1D

 SUBROUTINE r8READFLDS1D(NETCDFFNAME,RECID,rfld,FLDNAME,FLDUNIT)
  USE NETCDF
  IMPLICIT NONE

 !INPUT VARIABLES
  REAL(8),      INTENT(OUT) :: rfld(:)
  INTEGER(4),   INTENT(IN)  :: RECID
  CHARACTER(*), INTENT(IN)  :: NETCDFFNAME
  CHARACTER(*), INTENT(IN)  :: FLDNAME
  CHARACTER(*), INTENT(OUT) :: FLDUNIT

 !INPUT VARIABLES FACTS
  INTEGER(4)                :: MATDIMS1D(1)

 !NETCDF INPUTS VARIABLES
  INTEGER(4)                :: NCID
  INTEGER(4), PARAMETER     :: NDIMS = 2
  INTEGER(4)                :: start(NDIMS),count(NDIMS)
  INTEGER(4)                :: FLD_varid
  CHARACTER(5)              :: UNITS   = "units"

 !LOOP INDICES
  INTEGER(4)                :: ifld


  MATDIMS1D=SIZE(rfld)

  !Open the file. 
  call check(nf90_open(NETCDFFNAME, NF90_nowrite, ncid))

  !These settings tell netcdf to write one timestep of data.
  !(The setting of start(4) inside the loop below tells netCDF which timestep to write.)
  count = (/ MATDIMS1D(1), 1 /)
  start = (/ 1, RECID /)

  !Write the pretend data.
  !The arrays only hold one timestep worth of data. We will just rewrite the same data for each timestep.
  call check( nf90_inq_varid(ncid, TRIM(FLDNAME), FLD_varid) )
  call check( nf90_get_att(ncid, FLD_varid, UNITS, FLDUNIT ) )
  call check( nf90_get_var(ncid, FLD_varid, rfld, start = start, count = count) )
  
  !Close the file.
  call check( nf90_close(ncid) )

  RETURN
 END SUBROUTINE r8READFLDS1D



!*************!
! ERROR CHECK !
!*************!
 SUBROUTINE CHECK(status)
    integer, intent ( in) :: status
    
    if(status.ne.nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
 END SUBROUTINE CHECK

 subroutine handle_err(status,VarID)
  integer, intent ( in) :: status
  integer, optional     :: VarID
  character(3)          :: cVarID

  if(present(VarID)) then
   if(VarID.LT.10) then
    cVarID=char(48+VarID)
   elseif (VarId.LT.100) then
    cVarID=char(48+(VarID/10))//char(48+MOD(VarID,10))
   elseif (VarId.GT.100) then
    cVarID=char(48+(VarID/100))//char(48+(MOD(VarID,100)/10))//char(48+MOD(MOD(VarID,100),10))
   end if
  end if
  if(status /= nf90_noerr) then
   if(present(VarID)) then
    print *, "(Error) Variable ID ", cVarID, " - ", trim(nf90_strerror(status))
   else
    print *, trim(nf90_strerror(status))
   end if
   !stop "Stopped"
    stop
  end if
 end subroutine handle_err

 FUNCTION I4NUM2STR(NUM)
  IMPLICIT NONE
  CHARACTER(10) :: I4NUM2STR
  INTEGER(4)    :: NUM

  IF((NUM/1000000000).GE.1) THEN
   WRITE(I4NUM2STR,'(I10)') NUM
  ELSEIF((NUM/100000000).GE.1) THEN
   WRITE(I4NUM2STR,'(I9)') NUM
  ELSEIF((NUM/10000000).GE.1) THEN
   WRITE(I4NUM2STR,'(I8)') NUM
  ELSEIF((NUM/1000000).GE.1) THEN
   WRITE(I4NUM2STR,'(I7)') NUM
  ELSEIF((NUM/100000).GE.1) THEN
   WRITE(I4NUM2STR,'(I6)') NUM
  ELSEIF((NUM/10000).GE.1) THEN
   WRITE(I4NUM2STR,'(I5)') NUM
  ELSEIF((NUM/1000).GE.1) THEN
   WRITE(I4NUM2STR,'(I4)') NUM
  ELSEIF((NUM/100).GE.1) THEN
   WRITE(I4NUM2STR,'(I3)') NUM
  ELSEIF((NUM/10).GE.1) THEN
   WRITE(I4NUM2STR,'(I2)') NUM
  ELSE
   WRITE(I4NUM2STR,'(I1)') NUM
  ENDIF

  RETURN
 END FUNCTION I4NUM2STR


END MODULE RDSVNETCDF
