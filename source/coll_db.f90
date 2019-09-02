! ================================================================================================ !
!  Collimator Database Module
!  V.K. Berglyd Olsen, BE/ABP/HSS
!  Created: 2018-03-22
!  Updated: 2019-04-08
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
  integer,                  public,  save :: cdb_setPos   = 0       ! The position in the DB file where the settings start

  ! Main Database Arrays
  character(len=:), allocatable, public, save :: cdb_cName(:)       ! Collimator name
  character(len=:), allocatable, public, save :: cdb_cNameUC(:)     ! Collimator name upper case
  character(len=:), allocatable, public, save :: cdb_cMaterial(:)   ! Collimator material
  integer,          allocatable, public, save :: cdb_cFamily(:)     ! Collimator family
  real(kind=fPrec), allocatable, public, save :: cdb_cNSig(:)       ! Collimator sigma
  real(kind=fPrec), allocatable, public, save :: cdb_cNSigOrig(:)   ! Collimator sigma
  real(kind=fPrec), allocatable, public, save :: cdb_cLength(:)     ! Collimator length
  real(kind=fPrec), allocatable, public, save :: cdb_cOffset(:)     ! Collimator offset
  real(kind=fPrec), allocatable, public, save :: cdb_cRotation(:)   ! Collimator rotation
  real(kind=fPrec), allocatable, public, save :: cdb_cBx(:)         ! Collimator beta x
  real(kind=fPrec), allocatable, public, save :: cdb_cBy(:)         ! Collimator beta y
  logical,          allocatable, public, save :: cdb_cFound(:)      ! Found in lattice

  ! Additional Settings Arrays
  integer,          allocatable, public, save :: cdb_cSides(:)      ! 0 = two-sided, or 1,2 for single side 1 or 2

  ! Family Arrays
  character(len=:), allocatable, public, save :: cdb_famName(:)     ! Family name
  real(kind=fPrec), allocatable, public, save :: cdb_famNSig(:)     ! Family sigma
  real(kind=fPrec), allocatable, public, save :: cdb_famNSigOrig(:) ! Family sigma

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
  call alloc(cdb_cName,     mNameLen, cdb_nColl, " ",           "cdb_cName")
  call alloc(cdb_cNameUC,   mNameLen, cdb_nColl, " ",           "cdb_cNameUC")
  call alloc(cdb_cMaterial, 4,        cdb_nColl, " ",           "cdb_cMaterial")
  call alloc(cdb_cFamily,             cdb_nColl, 0,             "cdb_cFamily")
  call alloc(cdb_cNSig,               cdb_nColl, cdb_defColGap, "cdb_cNSig")
  call alloc(cdb_cNSigOrig,           cdb_nColl, cdb_defColGap, "cdb_cNSigOrig")
  call alloc(cdb_cLength,             cdb_nColl, zero,          "cdb_cLength")
  call alloc(cdb_cOffset,             cdb_nColl, zero,          "cdb_cOffset")
  call alloc(cdb_cRotation,           cdb_nColl, zero,          "cdb_cRotation")
  call alloc(cdb_cBx,                 cdb_nColl, zero,          "cdb_cBx")
  call alloc(cdb_cBy,                 cdb_nColl, zero,          "cdb_cBy")
  call alloc(cdb_cFound,              cdb_nColl, .false.,       "cdb_cFound")

  ! Additional Settings Arrays
  call alloc(cdb_cSides,              cdb_nColl, 0,             "cdb_cSides")

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
  if(st_debug) then
    call cdb_writeFam
    call cdb_writeDB
  end if

! ============================================================================ !
!  Post-Processing DB
! ============================================================================ !

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
  use string_tools
  use mod_units
  use mod_alloc
  use numerical_constants

  character(len=:), allocatable :: lnSplit(:)
  character(len=mInputLn) inLine
  real(kind=fPrec) nSig
  integer i, dbUnit, ioStat, nSplit, iLine, famID, iColl
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

  if(lnSplit(1) == "NSIG_FAM") then
    ! Collimator Family
    call chr_cast(lnSplit(3),nSig,cErr)
    call cdb_addFamily(trim(lnSplit(2)), nSig, famID, fExists)
    if(fExists .and. .not. cdb_doNSig) then
      ! If setting nsig in fort.3 is disabled, the DB values take precedence, so we overwrite them here
      cdb_famNSig(famID)     = nSig
      cdb_famNSigOrig(famID) = nSig
    end if
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
  cdb_cNameUC(iColl)   = chr_toUpper(lnSplit(1))
  cdb_cMaterial(iColl) = lnSplit(3)

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
!  Parses the old style database format with one calue per line.
! ================================================================================================ !
subroutine cdb_readDB_oldFormat

  use crcoall
  use parpro
  use string_tools
  use mod_units
  use numerical_constants

  character(len=mInputLn) inLine
  character(len=cdb_fNameLen) famName
  logical cErr
  integer j, dbUnit, ioStat, iLine, famID

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
    read(dbUnit,*,iostat=ioStat) cdb_cNameUC(j)
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
    end if
    cdb_cFamily(j) = famID

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
    write(dbNew,"(a,1x,a16,1x,f13.6)") "NSIG_FAM",cdb_famName(j),cdb_famNSig(j)
  end do
  write(dbNew,"(a)") "#"
  write(dbNew,"(a)") "# Collimators"

  write(dbNew,"(1a,a47,1x,a16,1x,a4,5(1x,a13))") "#",chr_rPad(" name",47),&
    "opening","mat.","length[m]","angle[deg]","offset[m]","beta_x[m]","beta_y[m]"

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
  use numerical_constants

  character(len=:), allocatable :: lnSplit(:)
  character(len=mInputLn) inLine
  integer i, dbUnit, ioStat, nSplit, iLine, iColl, iFam, iTemp
  logical cErr, fErr, isFam

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

  ! Look up the target collimator or family
  iFam  = -1
  iColl = -1
  isFam = .false.
  if(nSplit >= 2) then
    iFam  = cdb_getFamilyID(lnSplit(2))
    iColl = cdb_getCollimatorID(lnSplit(2))
    if(iFam == -1 .and. iColl == -1) then
      write(lerr,"(a)") "COLLDB> ERROR Could not find '"//trim(lnSplit(2))//"' in neither collimator nor family database"
      goto 30
    end if
    if(iFam > 0 .and. iColl > 0) then
      write(lerr,"(a)") "COLLDB> ERROR Found '"//trim(lnSplit(2))//"' in both collimator and family database"
      goto 30
    end if
    isFam = iFam > 0
  else
    write(lerr,"(a)") "COLLDB> ERROR Each keyword in the SETTINGS section of the collimator database requires "//&
      "a target collimator or family name"
    goto 30
  end if

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
    if(isFam) then
      do i=1,cdb_nColl
        if(cdb_cFamily(iFam) == iFam) then
          if(iTemp == 2) then
            call cdb_rotateCollimator(i, pi)
            write(lout,"(a,i0,a)") "COLLDB> Collimator family '"//trim(lnSplit(2))//"' is set as onesided and rotated (",iTemp,")"
          else
            write(lout,"(a,i0,a)") "COLLDB> Collimator family '"//trim(lnSplit(2))//"' is set as onesided (",iTemp,")"
          end if
          cdb_cSides(i) = 1
        end if
      end do
    else
      if(iTemp == 2) then
        call cdb_rotateCollimator(iColl, pi)
        write(lout,"(a,i0,a)") "COLLDB> Collimator '"//trim(lnSplit(2))//"' is set as onesided and rotated (",iTemp,")"
      else
        write(lout,"(a,i0,a)") "COLLDB> Collimator '"//trim(lnSplit(2))//"' is set as onesided (",iTemp,")"
      end if
      cdb_cSides(iColl) = 1
    end if


  ! ============= !

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
    this_name     = trim(adjustl(cdb_cNameUC(j)))//C_NULL_CHAR
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
  if(cdb_cRotation(collID) > twopi) then
    cdb_cRotation(collID) = cdb_cRotation(collID) - twopi
  else if(cdb_cRotation(collID) < zero) then
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
!  Created: 2019-03-21
!  Updated: 2019-03-21
!  Write family database to file
! ================================================================================================ !
subroutine cdb_writeFam

  use mod_units

  integer j, famUnit

  call f_requestUnit("coll_families_dump.dat", famUnit)
  call f_open(unit=famUnit,file="coll_families_dump.dat",formatted=.true.,mode="w",status="replace")

  write(famUnit,"(a16,2(1x,a13))") "# famName       ","nSig","nSigOrig"
  do j=1,cdb_nFam
    write(famUnit,"(a16,2(1x,f13.6))") cdb_famName(j),cdb_famNSig(j),cdb_famNSigOrig(j)
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

  write(dbUnit,"(a20,1x,a16,2(1x,a13),1x,a4,5(1x,a13))") "# collName          ","famName         ","nSig",&
    "nSigOrig","mat.","length","angle","offset","betax","betay"
  do j=1,cdb_nColl
    if(cdb_cFamily(j) > 0) then
      famName = cdb_famName(cdb_cFamily(j))
    else
      famName = " "
    end if
    write(dbUnit,"(a20,1x,a16,2(1x,f13.6),1x,a4,5(1x,f13.6))") cdb_cName(j),famName,cdb_cNSig(j),&
      cdb_cNSigOrig(j),cdb_cMaterial(j),cdb_cLength(j),cdb_cRotation(j),cdb_cOffset(j),cdb_cBx(j),cdb_cBy(j)
  end do

  flush(dbUnit)
  call f_close(dbUnit)

end subroutine cdb_writeDB

! ================================================================================================ !
!  Compatibility Functions for old collimation code assuming LHC naming convention
! ================================================================================================ !

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
    if(cdb_cNameUC(i)(1:3) == "TCP" .and. doOneSide .or. cdb_cNameUC(i)(1:4) == "TCDQ" .or. cdb_cNameUC(i)(1:5) == "TCXRP") then
      cdb_cSides(i) = 1
      write(lout,"(a)") "COLLDB> Collimator '"//trim(cdb_cName(i))//"' is treated as one sided"
    end if
  end do

end subroutine cdb_setLHCOneSided

end module coll_db
