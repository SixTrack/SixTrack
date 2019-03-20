! ================================================================================================ !
!  Collimator Database Module
! ================================================================================================ !
module collimation_db

  use parpro
  use floatPrecision

  implicit none

  character(len=mFileName), public,  save :: cdb_fileName = " "     ! Database file
  logical,                  private, save :: cdb_dbOld    = .false. ! Old or new DB format
  integer,                  public,  save :: cdb_nColl    = 0       ! Number of collimators
  integer,                  public,  save :: cdb_nfam     = 0       ! Number of collimator families

  ! Database arrays
  character(len=:), allocatable, public, save :: cdb_cName(:)     ! Collimator name
  character(len=:), allocatable, public, save :: cdb_cMaterial(:) ! Collimator material
  integer,          allocatable, public, save :: cdb_cFamily(:)   ! Collimator family
  real(kind=fPrec), allocatable, public, save :: cdb_cNSig(:)     ! Collimator sigma
  real(kind=fPrec), allocatable, public, save :: cdb_cLength(:)   ! Collimator length
  real(kind=fPrec), allocatable, public, save :: cdb_cOffset(:)   ! Collimator offset
  real(kind=fPrec), allocatable, public, save :: cdb_cRotation(:) ! Collimator rotation
  real(kind=fPrec), allocatable, public, save :: cdb_cBx(:)       ! Collimator beta x
  real(kind=fPrec), allocatable, public, save :: cdb_cBy(:)       ! Collimator beta y
  real(kind=fPrec), allocatable, public, save :: cdb_cTilt(:,:)   ! Collimator tilt
  logical,          allocatable, public, save :: cdb_cFound(:)    ! Found in lattice

  ! Family Arrays
  character(len=:), allocatable, public, save :: cdb_famName(:)  ! Family name
  real(kind=fPrec), allocatable, public, save :: cdb_famNSig(:)  ! Family sigma

contains

subroutine cdb_allocDB

  use mod_alloc
  use numerical_constants

  call alloc(cdb_cName,     mNameLen, cdb_nColl, " ",     "cdb_cName")
  call alloc(cdb_cMaterial, mNameLen, cdb_nColl, " ",     "cdb_cMaterial")
  call alloc(cdb_cFamily,   4,        cdb_nColl, " ",     "cdb_cFamily")
  call alloc(cdb_cNSig,               cdb_nColl, -1,      "cdb_cNSig")
  call alloc(cdb_cLength,             cdb_nColl, zero,    "cdb_cLength")
  call alloc(cdb_cOffset,             cdb_nColl, zero,    "cdb_cOffset")
  call alloc(cdb_cRotation,           cdb_nColl, zero,    "cdb_cRotation")
  call alloc(cdb_cBx,                 cdb_nColl, zero,    "cdb_cBx")
  call alloc(cdb_cBy,                 cdb_nColl, zero,    "cdb_cBy")
  call alloc(cdb_cTilt,               cdb_nColl, 2, zero, "cdb_cTilt")
  call alloc(cdb_cFound,              cdb_nColl, .false., "cdb_cFound")

end subroutine cdb_allocDB

subroutine cdb_allocFam

  use mod_alloc
  use numerical_constants

  call alloc(cdb_famName, 16, cdb_nFam, " ",  "cdb_famName")
  call alloc(cdb_famNSig,     cdb_nFam, zero, "cdb_famNSig")

end subroutine cdb_allocFam

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-03-19
!  Updated: 2019-03-20
! ================================================================================================ !
subroutine cdb_readCollDB(dbFile)

  use parpro
  use mod_units
  use mod_common
  use string_tools

  character(len=*), intent(in) :: dbFile

  character(len=:), allocatable :: lnSplit(:)
  character(len=mInputLn) inLine
  character(len=mNameLen) elemName
  character(len=16) cFam
  real(kind=fPrec)  cNSig
  integer dbUnit, ioStat, nLines, nSplit, collID, famID
  integer i, j, ix
  logical dbErr, spErr

  call f_requestUnit(dbFile, dbUnit)
  call f_open(unit=dbUnit,file=dbFile,formatted=.true.,mode="r",err=dbErr,status="old")
  if(dbErr) then
    write(lout,"(a)") "COLL> ERROR Could not open the collimator database file '"//trim(coll_db)//"'"
    call prror
  end if
  cdb_fileName = dbFile

! ============================================================================ !
!  Check DB Format: New or old, and count number of collimators
! ============================================================================ !

  nLines = 0

10 continue
  read(dbUnit,"(a)",end=20) inLine
  if(inLine(1:1) == "#") goto 10
  nLines = nLines + 1
  if(nLines == 1) then
    ! Check first line to count number of values
    call chr_split(inLine, lnSplit, nSplit, spErr)
    if(nSplit > 1) then
      cdb_dbOld = .false. ! New style DB (multi-column)
    else
      cdb_dbOld = .true.  ! Old style DB (single-column)
    end if
  end if
  goto 10

20 continue
  call f_close(dbUnit)

  if(cdb_dbOld) then
    call cdb_readDB_oldFormat
  else
    call cdb_readDB_newFormat
  end if

  ! Generate family names from DB

  goto 100

! ============================================================================ !
!  Parse new type DB
! ============================================================================ !
40 continue

! ============================================================================ !
!  Post-Processing DB
! ============================================================================ !
100 continue

  ! Map single elements to collimators
  do i=1,iu
    ix = ic(i)-nblo
    if(ix < 1) cycle

    elemName = chr_toLower(bez(ix))
    collID = -1
    do j=1,db_ncoll
      if(elemName == db_name2(j)) then
        collID = j
        exit
      end if
    end do

    if(collID == -1) then
      write(lout,"(a)") "COLL> WARNING Collimator not found in colldb: '"//trim(bez(ix))//"'"
    else
      db_elemMap(ix) = collID
    end if
  end do

end subroutine cdb_readCollDB

subroutine cdb_readDB_oldFormat

  use crcoall
  use parpro
  use string_tools
  use mod_units

  character(len=mInputLn) inLine
  logical cErr
  integer j, dbUnit, ioStat

  call f_requestUnit(cdb_fileName, dbUnit)
  call f_open(unit=dbUnit,file=cdb_fileName,formatted=.true.,mode="r",err=dbErr,status="old")

  read(coll_db_unit,*,iostat=ioStat) inLine
  if(ioStat /= 0) goto 100

  read(coll_db_unit,*,iostat=ioStat) cdb_nColl
  if(ioStat /= 0) goto 100
  call cdb_allocDB

  do j=1,cdb_nColl

    ! Line 1: Hash, ignored
    read(coll_db_unit,*,iostat=ioStat) inLine
    if(ioStat /= 0) goto 100

    ! Line 2: Upper case name, ignored
    read(coll_db_unit,*,iostat=ioStat) inLine
    if(ioStat /= 0) goto 100

    ! Line 3: Lower case name
    read(coll_db_unit,*,iostat=ioStat) cdb_cName(j)
    if(ioStat /= 0) goto 100

    ! Line 4: Sigma
    read(coll_db_unit,*,iostat=ioStat) inLine
    if(ioStat /= 0) goto 100
    call chr_cast(inLine, cdb_cNSig(j), cErr)
    if(cErr) goto 100

    ! Line 5: Material
    read(coll_db_unit,*,iostat=ioStat) cdb_cMaterial(j)
    if(ioStat /= 0) goto 100

    ! Line 6: Length
    read(coll_db_unit,*,iostat=ioStat) inLine
    if(ioStat /= 0) goto 100
    call chr_cast(inLine, cdb_cLength(j), cErr)
    if(cErr) goto 100

    ! Line 7: Rotation
    read(coll_db_unit,*,iostat=ioStat) inLine
    if(ioStat /= 0) goto 100
    call chr_cast(inLine, cdb_cRotation(j), cErr)
    if(cErr) goto 100

    ! Line 8: Offset
    read(coll_db_unit,*,iostat=ioStat) inLine
    if(ioStat /= 0) goto 100
    call chr_cast(inLine, cdb_cOffset(j), cErr)
    if(cErr) goto 100

    ! Line 9: Beta X
    read(coll_db_unit,*,iostat=ioStat) inLine
    if(ioStat /= 0) goto 100
    call chr_cast(inLine, cdb_cBx(j), cErr)
    if(cErr) goto 100

    ! Line 10: Beta Y
    read(coll_db_unit,*,iostat=ioStat) inLine
    if(ioStat /= 0) goto 100
    call chr_cast(inLine, cdb_cBy(j), cErr)
    if(cErr) goto 100

  end do

  call f_close(dbUnit)

  return

100 continue
  write(lout,"(a, i0)") "COLLDB> ERROR Cannot read DB file, iostat = ",ioStat
  call prror

end subroutine cdb_readDB_oldFormat

subroutine cdb_writeDB
! #ifdef ROOT
!   use iso_c_binding
!   use root_output
! #endif
! #ifdef HDF5
!   type(h5_dataField), allocatable :: fldCollDB(:)
!   character(len=:),   allocatable :: colNames(:)
!   character(len=:),   allocatable :: colUnits(:)
!   integer :: fmtCollDB, setCollDB, nSplit
!   logical :: spErr
! #endif

! #ifdef ROOT
! ! Temp variables to avoid fotran array -> C nightmares
!   character(len=mNameLen+1) :: this_name = C_NULL_CHAR
!   character(len=5) :: this_material = C_NULL_CHAR
! #endif
! #ifdef ROOT
!  do j=1,cdb_nColl
!     if(root_flag .and. root_CollimationDB.eq.1) then
!       this_name = trim(adjustl(db_name1(j))) // C_NULL_CHAR
!       this_material = trim(adjustl(db_material(j))) // C_NULL_CHAR
!       call CollimatorDatabaseRootWrite(j, this_name, len_trim(this_name), this_material, len_trim(this_material), db_nsig(j), &
!         db_length(j), db_rotation(j), db_offset(j))
!     end if
!   end do
! #endif
! #ifdef HDF5
!   if(h5_useForCOLL) then
!     allocate(fldCollDB(8))
!     fldCollDB(1) = h5_dataField(name="NAME",     type=h5_typeChar, size=mNameLen)
!     fldCollDB(2) = h5_dataField(name="OPENING",  type=h5_typeReal)
!     fldCollDB(3) = h5_dataField(name="MATERIAL", type=h5_typeChar, size=4)
!     fldCollDB(4) = h5_dataField(name="LENGTH",   type=h5_typeReal)
!     fldCollDB(5) = h5_dataField(name="ANGLE",    type=h5_typeReal)
!     fldCollDB(6) = h5_dataField(name="OFFSET",   type=h5_typeReal)
!     fldCollDB(7) = h5_dataField(name="BETAX",    type=h5_typeReal)
!     fldCollDB(8) = h5_dataField(name="BETAY",    type=h5_typeReal)
!     call h5_createFormat("collimation_db", fldCollDB, fmtCollDB)
!     call h5_createDataSet("collimation_db", h5_collID, fmtCollDB, setCollDB, db_ncoll)
!     call chr_split("name opening material length angle offset beta_x beta_y",colNames,nSplit,spErr)
!     call chr_split("text sigma text m rad m m m",colUnits,nSplit,spErr)
!     call h5_writeDataSetAttr(setCollDB,"nColl",   db_ncoll)
!     call h5_writeDataSetAttr(setCollDB,"colNames",colNames)
!     call h5_writeDataSetAttr(setCollDB,"colUnits",colUnits)
!     call h5_prepareWrite(setCollDB, db_ncoll)
!     call h5_writeData(setCollDB, 1, db_ncoll, db_name2(1:db_ncoll))
!     call h5_writeData(setCollDB, 2, db_ncoll, db_nsig(1:db_ncoll))
!     call h5_writeData(setCollDB, 3, db_ncoll, db_material(1:db_ncoll))
!     call h5_writeData(setCollDB, 4, db_ncoll, db_length(1:db_ncoll))
!     call h5_writeData(setCollDB, 5, db_ncoll, db_rotation(1:db_ncoll))
!     call h5_writeData(setCollDB, 6, db_ncoll, db_offset(1:db_ncoll))
!     call h5_writeData(setCollDB, 7, db_ncoll, db_bx(1:db_ncoll))
!     call h5_writeData(setCollDB, 8, db_ncoll, db_by(1:db_ncoll))
!     call h5_finaliseWrite(setCollDB)
!     deallocate(fldCollDB)
!   end if
! #endif
end subroutine cdb_writeDB

! ================================================================================================ !
!  Extract Family Name from Old Format DB
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Updated: 2019-03-20
! ================================================================================================ !
subroutine cdb_getFamily_oldDB

  do j=1,db_ncoll
    cFam     = " "
    cNSig    = zero
    elemName = chr_toLower(db_name2(j))
    if(elemName(1:3) == "tcp") then
      if(elemName(7:9) == "3.b") then
        cFam  = "tcp3"
        cNSig = nsig_tcp3
      else
        cFam  = "tcp7"
        cNSig = nsig_tcp7
      end if
    else if(elemName(1:4) == "tcsg") then
      if(elemName(8:10) == "3.b" .or. elemName(9:11) == "3.b") then
        cFam  = "tcsg3"
        cNSig = nsig_tcsg3
      else
        cFam  = "tcsg7"
        cNSig = nsig_tcsg7
      end if
      if(elemName(5:6) == ".4" .and. elemName(8:9) == "6.") then
        cFam  = "tcstcdq"
        cNSig = nsig_tcstcdq
      end if
    else if(elemName(1:4) == "tcsp") then
      if(elemName(9:11) == "6.b") then
        cFam  = "tcstcdq"
        cNSig = nsig_tcstcdq
      end if
    else if(elemName(1:4) == "tcsm") then
      if(elemName(8:10) == "3.b" .or. elemName(9:11) == "3.b") then
        cFam  = "tcsm3"
        cNSig = nsig_tcsm3
      else
        cFam  = "tcsm7"
        cNSig = nsig_tcsm7
      end if
    else if(elemName(1:4) == "tcla") then
      if(elemName(9:11) == "7.b") then
        cFam  = "tcla7"
        cNSig = nsig_tcla7
      else
        cFam  = "tcla3"
        cNSig = nsig_tcla3
      endif
    else if(elemName(1:4) == "tcdq") then
      cFam  = "tcdq"
      cNSig = nsig_tcdq
    else if(elemName(1:4) == "tcth" .or. elemName(1:5) == "tctph") then
      if(elemName(8:8) == "1" .or. elemName(9:9) == "1" ) then
        cFam  = "tcth1"
        cNSig = nsig_tcth1
      else if(elemName(8:8) == "2" .or. elemName(9:9) == "2") then
        cFam  = "tcth2"
        cNSig = nsig_tcth2
      else if(elemName(8:8) == "5" .or. elemName(9:9) == "5") then
        cFam  = "tcth5"
        cNSig = nsig_tcth5
      else if(elemName(8:8) == "8" .or. elemName(9:9) == "8") then
        cFam  = "tcth8"
        cNSig = nsig_tcth8
      end if
    else if(elemName(1:4) == "tctv" .or. elemName(1:5) == "tctpv") then
      if(elemName(8:8) == "1" .or. elemName(9:9) == "1") then
        cFam  = "tctv1"
        cNSig = nsig_tctv1
      else if(elemName(8:8) == "2" .or. elemName(9:9) == "2") then
        cFam  = "tctv2"
        cNSig = nsig_tctv2
      else if(elemName(8:8) == "5" .or. elemName(9:9) == "5") then
        cFam  = "tctv5"
        cNSig = nsig_tctv5
      else if(elemName(8:8) == "8" .or. elemName(9:9) == "8") then
        cFam  = "tctv8"
        cNSig = nsig_tctv8
      end if
    else if(elemName(1:3) == "tdi") then
      cFam  = "tdi"
      cNSig = nsig_tdi
    else if(elemName(1:4) == "tclp" .or. elemName(1:4) == "tcl." .or. elemName(1:4) == "tclx") then
      cFam  = "tclp"
      cNSig = nsig_tclp
    else if(elemName(1:4) == "tcli") then
      cFam  = "tcli"
      cNSig = nsig_tcli
    else if(elemName(1:4) == "tcxr") then
      cFam  = "tcxrp"
      cNSig = nsig_tcxrp
    else if(elemName(1:5) == "tcryo" .or. elemName(1:5) == "tcld.") then
      cFam  = "tcryo"
      cNSig = nsig_tcryo
    else
      write(lout,"(a)") "COLL> WARNING When setting opening for collimator '"//trim(elemName)//&
          "' from fort.3. Name not recognized. Setting nsig = 1000.0"
      cFam="default"
      cNSig=c1e3
    end if
    call collimate_getFamilyID(cFam,famID,.true.)
    coll_nSigFamily(famID) = cNSig
    db_family(j) = famID
  end do

end subroutine cdb_getFamily_oldDB

end module collimation_db
