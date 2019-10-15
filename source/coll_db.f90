! ================================================================================================ !
!  Collimator Database Module
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2018-03-22
!  Updated: 2019-09-10
! ================================================================================================ !
module coll_db

  use floatPrecision
  use parpro, only : mFileName
  use numerical_constants, only : zero

  implicit none

  public :: cdb_getFamilyID
  public :: cdb_getFamilyNSig

  integer,                  parameter     :: cdb_fNameLen  = 16     ! Length of collimator family name
  real(kind=fPrec),         parameter     :: cdb_defColGap = zero   ! Default collimator gap in sigma

  character(len=mFileName), public,  save :: cdb_fileName = " "     ! Database file
  logical,                  private, save :: cdb_dbOld    = .false. ! Old or new DB format
  logical,                  public,  save :: cdb_doNSig   = .false. ! Use the sigmas from fort.3 instead of DB
  integer,                  public,  save :: cdb_nColl    = 0       ! Number of collimators
  integer,                  public,  save :: cdb_nFam     = 0       ! Number of collimator families
  integer,                  public,  save :: cdb_setPos   = 0       ! The position in the DB file of the SETTINGS keyword

  ! Collimator Types (must be integer of power of 2)
  integer, parameter :: cdb_typPrimary   = 1
  integer, parameter :: cdb_typSecondary = 2
  integer, parameter :: cdb_typTertiary  = 4
  integer, parameter :: cdb_typOther     = 8
  integer, parameter :: cdb_typCrystal   = 16

  ! Main Database Arrays
  character(len=:), allocatable, public, save :: cdb_cName(:)       ! Collimator name
  character(len=:), allocatable, public, save :: cdb_cMaterial(:)   ! Collimator material
  integer,          allocatable, public, save :: cdb_cFamily(:)     ! Collimator family
  integer,          allocatable, public, save :: cdb_cType(:)       ! Collimator type
  real(kind=fPrec), allocatable, public, save :: cdb_cNSig(:)       ! Collimator sigma
  real(kind=fPrec), allocatable, public, save :: cdb_cNSigOrig(:)   ! Collimator sigma
  real(kind=fPrec), allocatable, public, save :: cdb_cLength(:)     ! Collimator length
  real(kind=fPrec), allocatable, public, save :: cdb_cOffset(:)     ! Collimator offset
  real(kind=fPrec), allocatable, public, save :: cdb_cRotation(:)   ! Collimator rotation
  real(kind=fPrec), allocatable, public, save :: cdb_cBx(:)         ! Collimator beta x
  real(kind=fPrec), allocatable, public, save :: cdb_cBy(:)         ! Collimator beta y
  logical,          allocatable, public, save :: cdb_cFound(:)      ! Found in lattice

  ! Additional Settings Arrays
  real(kind=fPrec), allocatable, public, save :: cdb_cTilt(:,:)     ! Collimator jaw tilt
  integer,          allocatable, public, save :: cdb_cMaterialID(:) ! Collimator material ID number
  integer,          allocatable, public, save :: cdb_cJawFit(:,:)   ! Collimator jaw fit index
  integer,          allocatable, public, save :: cdb_cSliced(:)     ! Collimator jaw fit sliced data index
  integer,          allocatable, public, save :: cdb_cSides(:)      ! 0 = two-sided, or 1,2 for single side 1 or 2

  ! Collimator Family Arrays
  character(len=:), allocatable, public, save :: cdb_famName(:)     ! Family name
  real(kind=fPrec), allocatable, public, save :: cdb_famNSig(:)     ! Family sigma
  real(kind=fPrec), allocatable, public, save :: cdb_famNSigOrig(:) ! Family sigma (original value from DB)
  integer,          allocatable, public, save :: cdb_famType(:)     ! Family type

  ! Element Map
  integer,          allocatable, public, save :: cdb_elemMap(:)     ! Map from single elements to DB

contains

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-03-19
!  Updated: 2019-08-30
!  Change the size of the collimator database arrays
! ================================================================================================ !
subroutine cdb_allocDB

  use parpro
  use mod_alloc
  use numerical_constants

  ! Main Database Arrays
  call alloc(cdb_cName,       mNameLen, cdb_nColl, " ",           "cdb_cName")
  call alloc(cdb_cMaterial,   4,        cdb_nColl, " ",           "cdb_cMaterial")
  call alloc(cdb_cFamily,               cdb_nColl, 0,             "cdb_cFamily")
  call alloc(cdb_cType,                 cdb_nColl, 0,             "cdb_cType")
  call alloc(cdb_cNSig,                 cdb_nColl, cdb_defColGap, "cdb_cNSig")
  call alloc(cdb_cNSigOrig,             cdb_nColl, cdb_defColGap, "cdb_cNSigOrig")
  call alloc(cdb_cLength,               cdb_nColl, zero,          "cdb_cLength")
  call alloc(cdb_cOffset,               cdb_nColl, zero,          "cdb_cOffset")
  call alloc(cdb_cRotation,             cdb_nColl, zero,          "cdb_cRotation")
  call alloc(cdb_cBx,                   cdb_nColl, zero,          "cdb_cBx")
  call alloc(cdb_cBy,                   cdb_nColl, zero,          "cdb_cBy")
  call alloc(cdb_cFound,                cdb_nColl, .false.,       "cdb_cFound")

  ! Additional Settings Arrays
  call alloc(cdb_cTilt,       2,        cdb_nColl, zero,          "cdb_cTilt")
  call alloc(cdb_cMaterialID,           cdb_nColl, 0,             "cdb_cMaterialID")
  call alloc(cdb_cJawFit,     2,        cdb_nColl, 0,             "cdb_cJawFit")
  call alloc(cdb_cSliced,               cdb_nColl, 0,             "cdb_cSliced")
  call alloc(cdb_cSides,                cdb_nColl, 0,             "cdb_cSides")

end subroutine cdb_allocDB

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-03-19
!  Updated: 2019-08-30
!  Change the size of the collimator family arrays
! ================================================================================================ !
subroutine cdb_allocFam

  use parpro
  use mod_alloc

  call alloc(cdb_famName, cdb_fNameLen, cdb_nFam, " ",           "cdb_famName")
  call alloc(cdb_famNSig,               cdb_nFam, cdb_defColGap, "cdb_famNSig")
  call alloc(cdb_famNSigOrig,           cdb_nFam, cdb_defColGap, "cdb_famNSigOrig")
  call alloc(cdb_famType,               cdb_nFam, 0,             "cdb_famType")

end subroutine cdb_allocFam

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-03-19
!  Updated: 2019-08-30
!  Change the size of other arrays depending on external size parameters
! ================================================================================================ !
subroutine cdb_expand_arrays(nele_new)
  use mod_alloc
  integer, intent(in) :: nele_new
  call alloc(cdb_elemMap,nele_new,0,"cdb_elemMap")
end subroutine cdb_expand_arrays

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-03-19
!  Updated: 2019-03-20
! ================================================================================================ !
subroutine cdb_readCollDB

  use parpro
  use crcoall
  use mod_units
  use mod_common
  use mod_settings
  use string_tools

  character(len=:), allocatable :: lnSplit(:)
  character(len=mInputLn) inLine
  character(len=cdb_fNameLen) cFam
  real(kind=fPrec)  cNSig
  integer dbUnit, ioStat, nSplit, collID, famID, elemEnd
  integer i, j, ix
  logical dbErr, spErr

  call f_requestUnit(cdb_fileName, dbUnit)
  call f_open(unit=dbUnit,file=cdb_fileName,formatted=.true.,mode="r",err=dbErr,status="old")
  if(dbErr) then
    write(lerr,"(a)") "COLLDB> ERROR Could not open the collimator database file '"//trim(cdb_fileName)//"'"
    call prror
  end if

! ============================================================================ !
!  Check DB Format: New or Old
! ============================================================================ !

10 continue
  read(dbUnit,"(a)",end=20) inLine
  if(inLine(1:1) == "#") goto 10
  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(nSplit > 1) then
    cdb_dbOld = .false. ! New style DB (multi-column)
  else
    cdb_dbOld = .true.  ! Old style DB (single-column)
  end if

20 continue
  call f_close(dbUnit)

  if(cdb_dbOld) then
    call cdb_readDB_oldFormat
  else
    call cdb_readDB_newFormat
  end if

  if(cdb_setPos > 0) then
    ! The DB has additional SETTINGS, parse them
    call cdb_readDBSettings
  end if

#ifdef ROOT
  call cdb_writeDB_ROOT
#endif
#ifdef HDF5
  call cdb_writeDB_HDF5
#endif

! ============================================================================ !
!  Post-Processing DB
! ============================================================================ !

  ! Set collimator types from family type, if we have any
  if(cdb_nFam > 0) then
    do i=1,cdb_nColl
      if(cdb_cFamily(i) > 0) then
        cdb_cType(i) = cdb_famType(cdb_cFamily(i))
      end if
    end do
  end if

  ! Map single elements to collimators
  do i=1,iu
    ix = ic(i)-nblo
    if(ix < 1) cycle

    collID = -1
    do j=1,cdb_nColl
      if(bez(ix) == cdb_cName(j)) then
        collID = j
        exit
      end if
    end do

    if(collID == -1) then
      if(bez(ix)(1:2) == "tc" .or. bez(ix)(1:2) == "td" .or. bez(ix)(1:3) == "col") then
        elemEnd = len_trim(bez(ix))
        if(bez(ix)(elemEnd-2:elemEnd) /= "_AP") then
          write(lout,"(a)") "COLLDB> WARNING Collimator not found in database: '"//trim(bez(ix))//"'"
        end if
      end if
    else
      cdb_elemMap(ix)    = collID
      cdb_cFound(collID) = .true.
    end if
  end do

  if(st_debug) then ! Dump a copy of the family and collimator databases
    call cdb_writeFam
    call cdb_writeDB
  end if

end subroutine cdb_readCollDB

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-03-19
!  Updated: 2019-08-30
!  Parsing the collimator section of the new database format. That is, the sigma settings and the
!  collimator descriptions. The parsing ends when it reaches the SETTINGS keyword.
! ================================================================================================ !
subroutine cdb_readDB_newFormat

  use parpro
  use crcoall
  use mod_alloc
  use mod_units
  use string_tools
  use coll_materials
  use numerical_constants

  character(len=:), allocatable :: lnSplit(:)
  character(len=mInputLn) inLine
  real(kind=fPrec) nSig
  integer i, dbUnit, ioStat, nSplit, iLine, famID, iColl, matID
  logical cErr, fErr, fExists

  fErr  = .false.
  cErr  = .false.
  iLine = 0
  iColl = 0

  write(lout,"(a)") "COLLDB> Reading collimator database, new format"

  call f_requestUnit(cdb_fileName, dbUnit)
  call f_open(unit=dbUnit,file=cdb_fileName,formatted=.true.,mode="r",status="old",err=fErr)
  if(fErr) then
    write(lerr,"(a)") "COLLDB> ERROR Cannot read from '"//trim(cdb_fileName)//"'"
    call prror
  end if

10 continue
  iLine = iLine + 1

  read(dbUnit,"(a)",end=20,iostat=ioStat) inLine
  if(ioStat /= 0) then
    write(lerr,"(a)") "COLLDB> ERROR Cannot read from '"//trim(cdb_fileName)//"'"
    call prror
  end if
  if(inLine(1:1) == "#") goto 10

  call chr_split(inLine, lnSplit, nSplit, cErr)
  if(cErr) then
    write(lerr,"(a,i0)") "COLLDB> ERROR Failed to parse database line ",iLine
    call prror
  end if
  if(nSplit == 0) goto 10 ! Skip empty lines

  if(lnSplit(1) == "SETTINGS") then
    write(lout,"(a,i0)") "COLLDB> SETTINGS flag encountered in collimator database on line ",iLine
    cdb_setPos = iLine
    goto 20
  end if

  if(lnSplit(1) == "NSIG_FAM") then ! Collimator Family
    if(nSplit /= 4) then
      write(lerr,"(a,i0,a)") "COLLDB> ERROR Collimator family description on line ",iLine," must be 4 values."
      write(lerr,"(a)")      "COLLDB>       NSIG_FAM famName sigmaSetting collType"
      call prror
    end if
    call chr_cast(lnSplit(3), nSig, cErr)
    call cdb_addFamily(trim(lnSplit(2)), nSig, famID, fExists)
    if(fExists .and. .not. cdb_doNSig) then
      ! If setting nsig in fort.3 is disabled, the DB values take precedence, so we overwrite them here
      cdb_famNSig(famID)     = nSig
      cdb_famNSigOrig(famID) = nSig
    end if
    select case(chr_toUpper(lnSplit(4)(1:3))) ! We only check the first three characters
    case("PRI")
      cdb_famType(famID) = cdb_typPrimary
    case("SEC")
      cdb_famType(famID) = cdb_typSecondary
    case("TER")
      cdb_famType(famID) = cdb_typTertiary
    case("OTH")
      cdb_famType(famID) = cdb_typOther
    case("CRY")
      cdb_famType(famID) = cdb_typCrystal
    case("UNK")
      cdb_famType(famID) = 0
    case default
      write(lerr,"(a,i0)") "COLLDB> ERROR Unknown collimator type '"//trim(lnSplit(4))//"' on line ",iLine
      call prror
    end select
    goto 10
  end if

  ! If not a family definition, it should be a collimator
  if(nSplit < 6) then
    write(lerr,"(a,i0,a)") "COLLDB> ERROR Collimator description on line ",iLine," has less than 6 values."
    call prror
  end if

  iColl = iColl + 1
  if(iColl > cdb_nColl) then
    cdb_nColl = cdb_nColl + 10
    call cdb_allocDB
  end if

  cdb_cName(iColl)     = lnSplit(1)
  cdb_cMaterial(iColl) = lnSplit(3)

  matID = collmat_getCollMatID(cdb_cMaterial(iColl))
  if(matID > 0) then
    cdb_cMaterialID(iColl) = matID
  else
    write(lerr,"(a)") "COLLDB> ERROR Material '"//trim(lnSplit(3))//"' not supported. Check your CollDB."
    call prror
  end if

  call chr_cast(lnSplit(4),cdb_cLength(iColl),  cErr)
  call chr_cast(lnSplit(5),cdb_cRotation(iColl),cErr)
  call chr_cast(lnSplit(6),cdb_cOffset(iColl),  cErr)

  cdb_cRotation(iColl) = cdb_cRotation(iColl)*rad

  if(nSplit > 6) call chr_cast(lnSplit(7),cdb_cBx(iColl),cErr)
  if(nSplit > 7) call chr_cast(lnSplit(8),cdb_cBy(iColl),cErr)

  if(chr_isNumeric(lnSplit(2))) then
    ! If column 3 is a number, we have no family assigned, so just use the value given.
    call chr_cast(lnSplit(2),cdb_cNSig(iColl),cErr)
    cdb_cNSigOrig(iColl) = cdb_cNSig(iColl)
  else
    ! Otherwise, we look up the name and require that it exists. We also save the family ID.
    famID = cdb_getFamilyID(lnSplit(2))
    if(famID == -1) then
      write(lerr,"(a,i0,a)") "COLLDB> ERROR Collimator opening '"//trim(adjustl(lnSplit(2)))//"' on line ",iLine,&
        " is not in family database"
      call prror
    else
      cdb_cFamily(iColl)   = famID
      cdb_cNSig(iColl)     = cdb_famNSig(famID)
      cdb_cNSigOrig(iColl) = cdb_famNSig(famID)
    end if
  end if

  goto 10

20 continue

  cdb_nColl = iColl
  call cdb_allocDB ! This should remove the unused lines in the DB

  call f_close(dbUnit)

end subroutine cdb_readDB_newFormat

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-03-19
!  Updated: 2019-09-02
!  Parses the old style database format with one value per line.
! ================================================================================================ !
subroutine cdb_readDB_oldFormat

  use parpro
  use crcoall
  use mod_units
  use string_tools
  use coll_materials
  use numerical_constants

  character(len=mInputLn) inLine
  character(len=cdb_fNameLen) famName
  character(len=mNameLen) collDummy
  logical cErr, fExists
  integer j, dbUnit, ioStat, iLine, famID, matID, collType

  cErr = .false.

  write(lout,"(a)") "COLLDB> Reading collimator database, old format"

  call f_requestUnit(cdb_fileName, dbUnit)
  call f_open(unit=dbUnit,file=cdb_fileName,formatted=.true.,mode="r",status="old")

  read(dbUnit,*,iostat=ioStat) inLine
  iLine = 1
  if(ioStat /= 0) goto 100

  read(dbUnit,*,iostat=ioStat) cdb_nColl
  iLine = 2
  if(ioStat /= 0) goto 100
  call cdb_allocDB

  do j=1,cdb_nColl

    ! Line 1: Hash, ignored
    read(dbUnit,*,iostat=ioStat) inLine
    iLine = iLine + 1
    if(ioStat /= 0) goto 100

    ! Line 2: Upper case name
    read(dbUnit,*,iostat=ioStat) collDummy
    iLine = iLine + 1
    if(ioStat /= 0) goto 100

    ! Line 3: Lower case name
    read(dbUnit,*,iostat=ioStat) cdb_cName(j)
    iLine = iLine + 1
    if(ioStat /= 0) goto 100

    ! Line 4: Collimator setting
    read(dbUnit,*,iostat=ioStat) inLine
    iLine = iLine + 1
    if(ioStat /= 0) goto 100
    call chr_cast(inLine, cdb_cNSigOrig(j), cErr)
    if(cErr) goto 100

    ! Line 5: Material
    read(dbUnit,*,iostat=ioStat) cdb_cMaterial(j)
    iLine = iLine + 1
    if(ioStat /= 0) goto 100

    ! Line 6: Length
    read(dbUnit,*,iostat=ioStat) inLine
    iLine = iLine + 1
    if(ioStat /= 0) goto 100
    call chr_cast(inLine, cdb_cLength(j), cErr)
    if(cErr) goto 100

    ! Line 7: Rotation
    read(dbUnit,*,iostat=ioStat) inLine
    iLine = iLine + 1
    if(ioStat /= 0) goto 100
    call chr_cast(inLine, cdb_cRotation(j), cErr)
    if(cErr) goto 100

    ! Line 8: Offset
    read(dbUnit,*,iostat=ioStat) inLine
    iLine = iLine + 1
    if(ioStat /= 0) goto 100
    call chr_cast(inLine, cdb_cOffset(j), cErr)
    if(cErr) goto 100

    ! Line 9: Beta X
    read(dbUnit,*,iostat=ioStat) inLine
    iLine = iLine + 1
    if(ioStat /= 0) goto 100
    call chr_cast(inLine, cdb_cBx(j), cErr)
    if(cErr) goto 100

    ! Line 10: Beta Y
    read(dbUnit,*,iostat=ioStat) inLine
    iLine = iLine + 1
    if(ioStat /= 0) goto 100
    call chr_cast(inLine, cdb_cBy(j), cErr)
    if(cErr) goto 100

    call cdb_generateFamName(cdb_cName(j), famName)
    famID = cdb_getFamilyID(famName)
    if(famID > 0 .and. cdb_doNSig) then
      cdb_cNSig(j) = cdb_famNSig(famID)
    else
      cdb_cNSig(j) = cdb_cNSigOrig(j)
      call cdb_addFamily(famName, cdb_cNSig(j), famID, fExists)
    end if
    cdb_cFamily(j) = famID

    call cdb_getCollType(cdb_cName(j), collType)
    cdb_famType(famID) = collType

    matID = collmat_getCollMatID(cdb_cMaterial(j))
    if(matID > 0) then
      cdb_cMaterialID(j) = matID
    else
      write(lerr,"(a)") "COLLDB> ERROR Material '"//trim(cdb_cMaterial(j))//"' not supported. Check your CollDB."
      call prror
    end if
  
  end do

  call f_freeUnit(dbUnit)

  return

100 continue
  write(lerr,"(2(a,i0))") "COLLDB> ERROR Cannot read DB file line ",iLine,", iostat = ",ioStat
  call prror

end subroutine cdb_readDB_oldFormat

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-03-19
!  Updated: 2019-09-02
!  Write a copy of the old style DB in new DB format.
! ================================================================================================ !
subroutine cdb_writeDB_newFromOld

  use crcoall
  use string_tools
  use mod_units
  use numerical_constants

  character(len=cdb_fNameLen) famName
  integer j, dbNew

  if(cdb_dbOld .eqv. .false.) then
    ! Already have a new format DB
    return
  end if

  write(lout,"(a)") "COLLDB> Converting old format DB to new format to file '"//trim(cdb_fileName)//".new'"

  call f_requestUnit(trim(cdb_fileName)//".new", dbNew)
  call f_open(unit=dbNew,file=trim(cdb_fileName)//".new",formatted=.true.,mode="w",status="replace")

  write(dbNew,"(a)") "# Automatically converted collimator DB from old format file '"//trim(cdb_fileName)//"'"
  write(dbNew,"(a)") "# Families"
  do j=1,cdb_nFam
    write(dbNew,"(a,1x,a16,1x,f13.6,1x,a)") "NSIG_FAM",cdb_famName(j),&
      cdb_famNSig(j),trim(cdb_getTypeName(cdb_famType(j)))
  end do
  write(dbNew,"(a)") "#"
  write(dbNew,"(a)") "# Collimators"

  write(dbNew,"(1a,a47,1x,a16,1x,a4,5(1x,a13))") "#",chr_rPad(" name",47),&
    "opening/fam","mat.","length[m]","angle[deg]","offset[m]","beta_x[m]","beta_y[m]"

  do j=1,cdb_nColl
    if(cdb_cFamily(j) > 0) then
      famName = cdb_famName(cdb_cFamily(j))
      write(dbNew,"(a48,1x,a16,1x,a4,5(1x,f13.6))") cdb_cName(j),&
      chr_lPad(trim(famName),16),cdb_cMaterial(j),cdb_cLength(j),cdb_cRotation(j)/rad,&
      cdb_cOffset(j),cdb_cBx(j),cdb_cBy(j)
    else
      write(dbNew,"(a48,1x,f16.3,1x,a4,5(1x,f13.6))") cdb_cName(j),&
      cdb_cNSig(j),cdb_cMaterial(j),cdb_cLength(j),cdb_cRotation(j)/rad,&
      cdb_cOffset(j),cdb_cBx(j),cdb_cBy(j)
    end if
  end do

  write(dbNew,"(a)") "#"
  write(dbNew,"(a)") "# Additional Collimator Settings"
  write(dbNew,"(a)") "SETTINGS"

  ! Onesided Collimators
  do j=1,cdb_nColl
    if(cdb_cSides(j) > 0) then
      write(dbNew,"(a,1x,a48,1x,i1)") "ONESIDED",cdb_cName(j),cdb_cSides(j)
    end if
  end do

  flush(dbNew)
  call f_freeUnit(dbNew)

end subroutine cdb_writeDB_newFromOld

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-08-30
!  Updated: 2019-08-30
!  Parse additional settings from the collimator database. This is treated separately since this
!  section is parsed in a standard name/value format like an input block in fort.3
! ================================================================================================ !
subroutine cdb_readDBSettings

  use parpro
  use crcoall
  use string_tools
  use mod_units
  use mod_alloc
  use coll_jawfit
  use numerical_constants

  character(len=:), allocatable :: lnSplit(:)
  character(len=mInputLn) inLine
  integer i, dbUnit, ioStat, nSplit, iLine, iColl, iFam, iTemp, iFit, fitID(2)
  logical cErr, fErr, isFam

  real(kind=fPrec) rParam(6)
  integer          iParam(3)
  logical          bParam(2)

  fErr  = .false.
  cErr  = .false.
  iLine = 0

  write(lout,"(a)") "COLLDB> Reading additional settings from collimator database"

  call f_requestUnit(cdb_fileName, dbUnit)
  call f_open(unit=dbUnit,file=cdb_fileName,formatted=.true.,mode="r",status="old",err=fErr)
  if(fErr) then
    write(lerr,"(a)") "COLLDB> ERROR Cannot read from '"//trim(cdb_fileName)//"'"
    call prror
  end if

10 continue
  iLine = iLine + 1

  read(dbUnit,"(a)",end=20,iostat=ioStat) inLine
  if(iLine <= cdb_setPos) goto 10 ! Skip already parsed lines

  if(ioStat /= 0) then
    write(lerr,"(a)") "COLLDB> ERROR Cannot read from '"//trim(cdb_fileName)//"'"
    call prror
  end if
  if(inLine(1:1) == "#") goto 10

  call chr_split(inLine, lnSplit, nSplit, cErr)
  if(cErr) then
    write(lerr,"(a)") "COLLDB> ERROR Failed to parse database line"
    goto 30
  end if
  if(nSplit == 0) goto 10 ! Skip empty lines

  iFam  = -1
  iColl = -1
  isFam = .false.

  ! Parse the keywords
  select case(lnSplit(1))

  case("ONESIDED") ! Treatment of one-sided collimators

    if(nSplit /= 3) then
      write(lerr,"(a,i0)") "COLLDB> ERROR ONESIDED expects 2 values, got ",nSplit-1
      write(lerr,"(a)")    "COLLDB>       ONESIDED collname|famname 1|2"
      goto 30
    end if
    call chr_cast(lnSplit(3), iTemp, cErr)
    if(iTemp /= 1 .and. iTemp /= 2) then
      write(lerr,"(a,i0)") "COLLDB> ERROR ONESIDED collimator value must be 1 or 2, got ",iTemp
      goto 30
    end if

    call cdb_getCollimatorOrFamilyID(lnSplit(2), iFam, iColl, isFam, cErr)
    if(cErr) goto 30
    if(isFam) then
      do i=1,cdb_nColl
        if(cdb_cFamily(i) == iFam) then
          if(iTemp == 2) then
            call cdb_rotateCollimator(i, pi)
            write(lout,"(a,i0,a)") "COLLDB> Collimator family '"//trim(lnSplit(2))//&
              "' is set as one-sided and rotated by pi (",iTemp,")"
          else
            write(lout,"(a,i0,a)") "COLLDB> Collimator family '"//trim(lnSplit(2))//"' is set as one-sided (",iTemp,")"
          end if
          cdb_cSides(i) = 1
        end if
      end do
    else
      if(iTemp == 2) then
        call cdb_rotateCollimator(iColl, pi)
        write(lout,"(a,i0,a)") "COLLDB> Collimator '"//trim(lnSplit(2))//"' is set as one-sided and rotated by pi (",iTemp,")"
      else
        write(lout,"(a,i0,a)") "COLLDB> Collimator '"//trim(lnSplit(2))//"' is set as one-sided (",iTemp,")"
      end if
      cdb_cSides(iColl) = 1
    end if

  case("JAW_PROFILE") ! Adding a Jaw Fit Profile
    if(nSplit < 3 .and. nSplit > 8) then
      write(lerr,"(a,i0)") "COLLDB> ERROR JAW_PROFILE expects 2 to 7 values, got ",nSplit-1
      write(lerr,"(a)")    "COLLDB>       JAW_PROFILE name fac0 [... fac5]"
      goto 30
    end if
    if(len_trim(lnSplit(2)) > jaw_fitNameLen) then
      write(lerr,"(2(a,i0))") "COLLDB> ERROR JAW_PROFILE name cannot be more than ",jaw_fitNameLen,&
        " characters, got ",len_trim(lnSplit(2))
      goto 30
    end if

    rParam(:) = zero
    if(nSplit > 2) call chr_cast(lnSplit(3), rParam(1), cErr)
    if(nSplit > 3) call chr_cast(lnSplit(4), rParam(2), cErr)
    if(nSplit > 4) call chr_cast(lnSplit(5), rParam(3), cErr)
    if(nSplit > 5) call chr_cast(lnSplit(6), rParam(4), cErr)
    if(nSplit > 6) call chr_cast(lnSplit(7), rParam(5), cErr)
    if(nSplit > 7) call chr_cast(lnSplit(8), rParam(6), cErr)

    call jaw_addJawFit(trim(lnSplit(2)), rParam(1:6), iFit, cErr)

    if(cErr) goto 30

  case("JAW_FIT") ! Apply Jaw Fit Profile
    if(nSplit /= 5 .and. nSplit /= 7 .and. nSplit /= 9) then
      write(lerr,"(a,i0)") "COLLDB> ERROR JAW_FIT expects 4, 6 ot 8 values, got ",nSplit-1
      write(lerr,"(a)")    "COLLDB>       JAW_FIT collname|famname nslices fit1 fit2 [scale1 scale2 [recentre1 recentre2]]"
      goto 30
    end if

    rParam(:) = one
    iParam(:) = -1
    bParam(:) = .false.
    if(nSplit > 2) call chr_cast(lnSplit(3), iParam(1), cErr)
    if(nSplit > 5) call chr_cast(lnSplit(6), rParam(1), cErr)
    if(nSplit > 6) call chr_cast(lnSplit(7), rParam(2), cErr)
    if(nSplit > 7) call chr_cast(lnSplit(8), bParam(1), cErr)
    if(nSplit > 8) call chr_cast(lnSplit(9), bParam(2), cErr)

    fitID(1) = jaw_getFitID(lnSplit(4))
    fitID(2) = jaw_getFitID(lnSplit(5))
    do i=1,2
      if(fitID(i) == -1) then
        write(lerr,"(a)") "COLLDB> ERROR Unknown fit profile '"//trim(lnSplit(3+1))//&
          "', did you forget to add it first?"
        goto 30
      end if
    end do

    call cdb_getCollimatorOrFamilyID(lnSplit(2), iFam, iColl, isFam, cErr)
    if(cErr) goto 30
    if(isFam) then
      do i=1,cdb_nColl
        if(cdb_cFamily(i) == iFam) then
          if(cdb_cSliced(i) /= 0) then
            write(lerr,"(a)") "COLLDB> ERROR Collimator '"//trim(cdb_cName(i))//"' has already been sliced"
            call prror
          end if
          cdb_cJawFit(:,i) = fitID
          call jaw_computeFit(cdb_cName(i), fitID, iParam(1), rParam(1:2), bParam(1:2), cdb_cLength(i), &
            cdb_cTilt(:,i), cdb_cOffset(i), iFit)
          cdb_cSliced(i) = iFit
        end if
      end do
    else
      if(cdb_cSliced(iColl) /= 0) then
        write(lerr,"(a)") "COLLDB> ERROR Collimator '"//trim(cdb_cName(iColl))//"' has already been sliced"
        call prror
      end if
      cdb_cJawFit(:,iColl) = fitID
      call jaw_computeFit(cdb_cName(iColl), fitID, iParam(1), rParam(1:2), bParam(1:2), cdb_cLength(iColl), &
        cdb_cTilt(:,iColl), cdb_cOffset(iColl), iFit)
      cdb_cSliced(iColl) = iFit
    end if

  case default
    write(lerr,"(a)") "COLLDB> ERROR Unknown keyword '"//trim(lnSplit(1))//"' in SETTINGS section"
    goto 30

  end select

  goto 10

20 continue

  call f_close(dbUnit)
  return

30 continue
  write(lerr,"(a,i0)") "COLLDB> ERROR Collimator DB '"//trim(cdb_fileName)//"' on line ",iLine
  call prror

end subroutine cdb_readDBSettings

#ifdef ROOT
subroutine cdb_writeDB_ROOT

  use parpro
  use iso_c_binding
  use root_output

  character(len=mNameLen+1) :: this_name     = C_NULL_CHAR
  character(len=5)          :: this_material = C_NULL_CHAR
  integer j

  if((root_flag .eqv. .false.) .or. root_CollimationDB /= 1) return

  do j=1,cdb_nColl
    this_name     = trim(adjustl(cdb_cName(j)))//C_NULL_CHAR
    this_material = trim(adjustl(cdb_cMaterial(j)))//C_NULL_CHAR
    call CollimatorDatabaseRootWrite(j, this_name, len_trim(this_name), this_material, len_trim(this_material), cdb_cNSig(j), &
      cdb_cLength(j), cdb_cRotation(j), cdb_cOffset(j))
  end do

end subroutine cdb_writeDB_ROOT
#endif

#ifdef HDF5
subroutine cdb_writeDB_HDF5

  use hdf5_output
  use string_tools

  type(h5_dataField), allocatable :: fldCollDB(:)
  character(len=:),   allocatable :: colNames(:)
  character(len=:),   allocatable :: colUnits(:)

  integer :: fmtCollDB, setCollDB, nSplit
  logical :: spErr

  if(h5_useForCOLL .eqv. .false.) return

  allocate(fldCollDB(8))

  fldCollDB(1) = h5_dataField(name="NAME",     type=h5_typeChar, size=mNameLen)
  fldCollDB(2) = h5_dataField(name="OPENING",  type=h5_typeReal)
  fldCollDB(3) = h5_dataField(name="MATERIAL", type=h5_typeChar, size=4)
  fldCollDB(4) = h5_dataField(name="LENGTH",   type=h5_typeReal)
  fldCollDB(5) = h5_dataField(name="ANGLE",    type=h5_typeReal)
  fldCollDB(6) = h5_dataField(name="OFFSET",   type=h5_typeReal)
  fldCollDB(7) = h5_dataField(name="BETAX",    type=h5_typeReal)
  fldCollDB(8) = h5_dataField(name="BETAY",    type=h5_typeReal)

  call h5_createFormat("collimation_db", fldCollDB, fmtCollDB)
  call h5_createDataSet("collimation_db", h5_collID, fmtCollDB, setCollDB, cdb_nColl)

  call chr_split("name opening material length angle offset beta_x beta_y",colNames,nSplit,spErr)
  call chr_split("text sigma text m rad m m m",colUnits,nSplit,spErr)

  call h5_writeDataSetAttr(setCollDB,"nColl",   cdb_nColl)
  call h5_writeDataSetAttr(setCollDB,"colNames",colNames)
  call h5_writeDataSetAttr(setCollDB,"colUnits",colUnits)

  call h5_prepareWrite(setCollDB, cdb_nColl)
  call h5_writeData(setCollDB, 1, cdb_nColl, cdb_cName(1:cdb_nColl))
  call h5_writeData(setCollDB, 2, cdb_nColl, cdb_cNSig(1:cdb_nColl))
  call h5_writeData(setCollDB, 3, cdb_nColl, cdb_cMaterial(1:cdb_nColl))
  call h5_writeData(setCollDB, 4, cdb_nColl, cdb_cLength(1:cdb_nColl))
  call h5_writeData(setCollDB, 5, cdb_nColl, cdb_cRotation(1:cdb_nColl))
  call h5_writeData(setCollDB, 6, cdb_nColl, cdb_cOffset(1:cdb_nColl))
  call h5_writeData(setCollDB, 7, cdb_nColl, cdb_cBx(1:cdb_nColl))
  call h5_writeData(setCollDB, 8, cdb_nColl, cdb_cBy(1:cdb_nColl))
  call h5_finaliseWrite(setCollDB)

  deallocate(fldCollDB)

end subroutine cdb_writeDB_HDF5
#endif

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-02
!  Updated: 2019-09-02
!  Rotate a collimater by rAngle radians
!  Note: does not check that zero <= rAngle <= twopi
! ================================================================================================ !
subroutine cdb_rotateCollimator(collID, rAngle)

  use numerical_constants

  integer,          intent(in) :: collID
  real(kind=fPrec), intent(in) :: rAngle

  cdb_cRotation(collID) = cdb_cRotation(collID) + rAngle
  if(cdb_cRotation(collID) >= twopi) then
    cdb_cRotation(collID) = modulo(cdb_cRotation(collID), twopi)
  end if
  if(cdb_cRotation(collID) < zero) then
    cdb_cRotation(collID) = cdb_cRotation(collID) + twopi
  end if

end subroutine cdb_rotateCollimator

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-03-22
!  Updated: 2019-03-22
!  Add a family to the family database
! ================================================================================================ !
subroutine cdb_addFamily(famName, nSig, famID, fExists)

  use crcoall

  character(len=*),  intent(in)  :: famName
  real(kind=fPrec),  intent(in)  :: nSig
  integer,           intent(out) :: famID
  logical,           intent(out) :: fExists

  famID = cdb_getFamilyID(famName)
  if(famID == -1) then
    cdb_nfam = cdb_nFam + 1
    call cdb_allocFam
    cdb_famName(cdb_nFam)     = famName
    cdb_famNSig(cdb_nFam)     = nSig
    cdb_famNSigOrig(cdb_nFam) = nSig
    famID   = cdb_nFam
    fExists = .false.
  else
    fExists = .true.
  end if

end subroutine cdb_addFamily

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-03-21
!  Updated: 2019-03-22
!  Find a family in the database and returns its ID
! ================================================================================================ !
integer function cdb_getFamilyID(famName) result(famID)

  character(len=*), intent(in) :: famName
  integer i

  famID = -1
  if(cdb_nFam > 0) then
    do i=1,cdb_nFam
      if(cdb_famName(i) == famName) then
        famID = i
        exit
      end if
    end do
  end if

end function cdb_getFamilyID

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-08-30
!  Updated: 2019-08-30
!  Find a collimator in the database and returns its ID
! ================================================================================================ !
integer function cdb_getCollimatorID(collName) result(collID)

  character(len=*), intent(in) :: collName
  integer i

  collID = -1
  if(cdb_nColl > 0) then
    do i=1,cdb_nColl
      if(cdb_cName(i) == collName) then
        collID = i
        exit
      end if
    end do
  end if

end function cdb_getCollimatorID

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-10
!  Updated: 2019-09-10
!  Look up a name in both family and collimator database
! ================================================================================================ !
subroutine cdb_getCollimatorOrFamilyID(itemName, iFam, iColl, isFam, fErr)

  use crcoall

  character(len=*), intent(in)    :: itemName
  integer,          intent(out)   :: iFam
  integer,          intent(out)   :: iColl
  logical,          intent(out)   :: isFam
  logical,          intent(inout) :: fErr

  iFam  = -1
  iColl = -1
  isFam = .false.

  iFam  = cdb_getFamilyID(itemName)
  iColl = cdb_getCollimatorID(itemName)
  if(iFam == -1 .and. iColl == -1) then
    write(lerr,"(a)") "COLLDB> ERROR Could not find '"//trim(itemName)//"' in neither collimator nor family database"
    fErr = .true.
    return
  end if
  if(iFam > 0 .and. iColl > 0) then
    write(lerr,"(a)") "COLLDB> ERROR Found '"//trim(itemName)//"' in both collimator and family database"
    fErr = .true.
    return
  end if

  isFam = iFam > 0

end subroutine cdb_getCollimatorOrFamilyID

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-03-21
!  Updated: 2019-03-22
!  Find a family in the database and returns its nsig
! ================================================================================================ !
real(kind=fPrec) function cdb_getFamilyNSig(famName) result(nSig)

  use crcoall

  character(len=*), intent(in) :: famName
  integer famID

  famID = cdb_getFamilyID(famName)
  if(famID > 0) then
    nSig = cdb_famNSig(famID)
  else
    write(lout,"(a)") "COLLDB> Warning No nsig value found for collimator family '"//trim(famName)//"'"
    nSig = 0.0_fPrec
  end if

end function cdb_getFamilyNSig

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-10-07
!  Updated: 2019-10-07
!  Get name from collimator type ID
! ================================================================================================ !
character(len=9) function cdb_getTypeName(typID)

  integer, intent(in) :: typID

  if(iand(typID,cdb_typPrimary) == cdb_typPrimary) then
    cdb_getTypeName = "PRIMARY"
  elseif(iand(typID,cdb_typSecondary) == cdb_typSecondary) then
    cdb_getTypeName = "SECONDARY"
  elseif(iand(typID,cdb_typTertiary) == cdb_typTertiary) then
    cdb_getTypeName = "TERTIARY"
  elseif(iand(typID,cdb_typOther) == cdb_typOther) then
    cdb_getTypeName = "OTHER"
  elseif(iand(typID,cdb_typCrystal) == cdb_typCrystal) then
    cdb_getTypeName = "CRYSTAL"
  else
    cdb_getTypeName = "UNKNOWN"
  end if

end function cdb_getTypeName

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-03-21
!  Updated: 2019-03-21
!  Write family database to file
! ================================================================================================ !
subroutine cdb_writeFam

  use mod_units

  integer j, famUnit

  call f_requestUnit("coll_families_dump.dat", famUnit)
  call f_open(unit=famUnit,file="coll_families_dump.dat",formatted=.true.,mode="w",status="replace")

  write(famUnit,"(a16,1x,a3,2(1x,a13))") "# famName       ","typ","nSig","nSigOrig"
  do j=1,cdb_nFam
    write(famUnit,"(a16,1x,a3,2(1x,f13.6))") cdb_famName(j),cdb_getTypeName(cdb_famType(j)), &
      cdb_famNSig(j),cdb_famNSigOrig(j)
  end do

  flush(famUnit)
  call f_close(famUnit)

end subroutine cdb_writeFam

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-04-09
!  Updated: 2019-04-09
!  Write collimator database to file
! ================================================================================================ !
subroutine cdb_writeDB

  use mod_units

  character(len=cdb_fNameLen) famName
  integer j, dbUnit

  call f_requestUnit("coll_db_dump.dat", dbUnit)
  call f_open(unit=dbUnit,file="coll_db_dump.dat",formatted=.true.,mode="w",status="replace")

  write(dbUnit,"(a20,1x,a16,1x,a3,2(1x,a13),1x,a4,5(1x,a13))") "# collName          ","famName         ",&
    "typ","nSig","nSigOrig","mat.","length","angle","offset","betax","betay"
  do j=1,cdb_nColl
    if(cdb_cFamily(j) > 0) then
      famName = cdb_famName(cdb_cFamily(j))
    else
      famName = " "
    end if
    write(dbUnit,"(a20,1x,a16,1x,a3,2(1x,f13.6),1x,a4,5(1x,f13.6))") cdb_cName(j),famName,          &
      cdb_getTypeName(cdb_cType(j)),cdb_cNSig(j),cdb_cNSigOrig(j),cdb_cMaterial(j),cdb_cLength(j),  &
      cdb_cRotation(j),cdb_cOffset(j),cdb_cBx(j),cdb_cBy(j)
  end do

  flush(dbUnit)
  call f_close(dbUnit)

end subroutine cdb_writeDB

! ================================================================================================ !
!  Compatibility Functions for old collimation code assuming LHC naming convention
! ================================================================================================ !

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-08-01
!  Updated: 2019-09-10
!  Set jaw fit from fort.3
! ================================================================================================ !
subroutine cdb_setMasterJawFit(nSlices, sMin, sMax, rc1, rc2, jawFit, fitScale)

  use parpro
  use crcoall
  use coll_jawfit
  use mod_common
  use mod_common_track
  use numerical_constants

  integer,          intent(in) :: nSlices
  real(kind=fPrec), intent(in) :: sMin, sMax
  real(kind=fPrec), intent(in) :: rc1, rc2
  real(kind=fPrec), intent(in) :: jawFit(2,6)
  real(kind=fPrec), intent(in) :: fitScale(2)

  integer i, ix, k, fitID(2), sliceID
  logical reCentre(2), fErr

  if(nSlices < 1) then
    return
  end if

  reCentre(:) = .false.
  if(rc1 /= zero) reCentre(1) = .true.
  if(rc2 /= zero) reCentre(2) = .true.

  fErr = .false.
  call jaw_addJawFit("FIT_1", jawFit(1,:), fitID(1), fErr)
  call jaw_addJawFit("FIT_2", jawFit(2,:), fitID(2), fErr)
  if(fErr) then
    write(lerr,"(a)") "COLLDB> ERROR While setting up jaw fit parameters"
    call prror
  end if

  do i=1,iu
    ix = ic(i)
    if(ix > nblo) then
      ix = ix-nblo
      k  = cdb_elemMap(ix)
      if(k > 0 .and. dcum(i) > sMin .and. dcum(i) < sMax) then
        if(cdb_cName(k)(1:4) == "tcsg" .or. cdb_cName(k)(1:3) == "tcp"  .or. &
           cdb_cName(k)(1:4) == "tcla" .or. cdb_cName(k)(1:3) == "tct"  .or. &
           cdb_cName(k)(1:4) == "tcli" .or. cdb_cName(k)(1:4) == "tcl." .or. &
           cdb_cName(k)(1:5) == "tcryo") then
          write(lout,"(a,f13.6)") "COLLDB> Will apply jaw fit to collimator '"//trim(bez(ix))//"' at position ",dcum(i)
          cdb_cJawFit(:,k) = fitID
          call jaw_computeFit(trim(bez(ix)), fitID, nSlices, fitScale, reCentre, cdb_cLength(k), cdb_cTilt(:,k), &
            cdb_cOffset(k), sliceID)
          cdb_cSliced(k) = sliceID
        end if
      end if
    end if
  end do

end subroutine cdb_setMasterJawFit

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-03-19
!  Updated: 2019-03-20
!  Extract Family Name from Old Format DB
! ================================================================================================ !
subroutine cdb_generateFamName(inElem, famName)

  use string_tools
  use numerical_constants

  character(len=mNameLen),     intent(in)  :: inElem
  character(len=cdb_fNameLen), intent(out) :: famName

  character(len=mNameLen) elemName

  famName  = " "
  elemName = chr_toLower(inElem)
  if(elemName(1:3) == "tcp") then
    if(elemName(7:9) == "3.b") then
      famName = "tcp3"
    else
      famName = "tcp7"
    end if
  else if(elemName(1:4) == "tcsg") then
    if(elemName(8:10) == "3.b" .or. elemName(9:11) == "3.b") then
      famName = "tcsg3"
    else
      famName = "tcsg7"
    end if
    if(elemName(5:6) == ".4" .and. elemName(8:9) == "6.") then
      famName = "tcstcdq"
    end if
  else if(elemName(1:4) == "tcsp") then
    if(elemName(9:11) == "6.b") then
      famName = "tcstcdq"
    end if
  else if(elemName(1:4) == "tcsm") then
    if(elemName(8:10) == "3.b" .or. elemName(9:11) == "3.b") then
      famName = "tcsm3"
    else
      famName = "tcsm7"
    end if
  else if(elemName(1:4) == "tcla") then
    if(elemName(9:11) == "7.b") then
      famName = "tcla7"
    else
      famName = "tcla3"
    endif
  else if(elemName(1:4) == "tcdq") then
    famName = "tcdq"
  else if(elemName(1:4) == "tcth" .or. elemName(1:5) == "tctxh" .or. elemName(1:5) == "tctph") then
! else if(elemName(1:4) == "tcth" .or. elemName(1:5) == "tctph") then
    if(elemName(8:8) == "1" .or. elemName(9:9) == "1" ) then
      famName = "tcth1"
    else if(elemName(8:8) == "2" .or. elemName(9:9) == "2") then
      famName = "tcth2"
    else if(elemName(8:8) == "5" .or. elemName(9:9) == "5") then
      famName = "tcth5"
    else if(elemName(8:8) == "8" .or. elemName(9:9) == "8") then
      famName = "tcth8"
    end if
  else if(elemName(1:4) == "tctv" .or. elemName(1:5) == "tctpv") then
    if(elemName(8:8) == "1" .or. elemName(9:9) == "1") then
      famName = "tctv1"
    else if(elemName(8:8) == "2" .or. elemName(9:9) == "2") then
      famName = "tctv2"
    else if(elemName(8:8) == "5" .or. elemName(9:9) == "5") then
      famName = "tctv5"
    else if(elemName(8:8) == "8" .or. elemName(9:9) == "8") then
      famName = "tctv8"
    end if
  else if(elemName(1:3) == "tdi") then
    famName = "tdi"
  else if(elemName(1:4) == "tclp" .or. elemName(1:4) == "tcl." .or. elemName(1:4) == "tclx") then
    famName = "tclp"
  else if(elemName(1:4) == "tcli") then
    famName = "tcli"
  else if(elemName(1:4) == "tcxr" .or. elemName(1:3) == "xrp") then
    famName = "tcxrp"
  else if(elemName(1:5) == "tcryo" .or. elemName(1:5) == "tcld.") then
    famName = "tcryo"
  else if(elemName(1:3) == "col") then
    if(elemName == "colm" .or. elemName(1:5) == "colh0") then
      famName = "tcth1"
    elseif(elemName(1:5) == "colv0") then
      famName = "tcth2"
    else if(elemName(1:5) == "colh1") then
      famName = "tcth5"
    else if(elemName(1:5) == "colv1") then
      famName = "tcth8"
    else if(elemName(1:5) == "colh2") then
      famName = "tctv1"
    end if
  else
    famName = "NONE"
  end if

end subroutine cdb_generateFamName

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-10-07
!  Updated: 2019-10-07
!  Checks the type of the collimator based on LHC naming convention (for support with old format)
! ================================================================================================ !
subroutine cdb_getCollType(inElem, collType)

  use string_tools

  character(len=mNameLen), intent(in)  :: inElem
  integer,                 intent(out) :: collType

  character(len=mNameLen) elemName

  collType = cdb_typOther
  elemName = chr_toLower(inElem)

  if(elemName(1:3) == "tcp") then
    collType = cdb_typPrimary
  else if(elemName(1:3) == "tcs") then
    collType = cdb_typSecondary
  else if(elemName(1:3) == "tcl" .or. elemName(1:3) == "tct" .or. &
          elemName(1:3) == "tcd" .or. elemName(1:3) == "tdi") then
    collType = cdb_typTertiary
  end if

end subroutine cdb_getCollType

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-09-02
!  Updated: 2019-09-02
!  Checks for one sided collimators from old input block and DB using the COLL block flag and
!  hardcoded onesided treatment for TCDQ and roman pots.
! ================================================================================================ !
subroutine cdb_setLHCOnesided(doOneSide)

  use crcoall

  logical, intent(in) :: doOneSide

  integer i

  if(cdb_dbOld .eqv. .false.) then
    ! Only do this if we're using the old DB format
    return
  end if

  do i=1,cdb_nColl
    cdb_cSides(i) = 0
    if(cdb_cName(i)(1:3) == "tcp" .and. doOneSide .or. cdb_cName(i)(1:4) == "tcdq" .or. cdb_cName(i)(1:5) == "tcxrp") then
      cdb_cSides(i) = 1
      write(lout,"(a)") "COLLDB> Collimator '"//trim(cdb_cName(i))//"' is treated as one-sided"
    end if
  end do

end subroutine cdb_setLHCOneSided

end module coll_db
