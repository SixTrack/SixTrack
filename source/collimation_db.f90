! ================================================================================================ !
!  Collimator Database Module
! ================================================================================================ !
module collimation_db

  use parpro
  use floatPrecision

  implicit none

  integer, parameter :: cdb_fNameLen = 16

  character(len=mFileName), public,  save :: cdb_fileName = " "     ! Database file
  logical,                  private, save :: cdb_dbOld    = .false. ! Old or new DB format
  logical,                  public,  save :: cdb_doNSig   = .false. ! Use the sigmas from fort.3 isntead of DB
  integer,                  public,  save :: cdb_nColl    = 0       ! Number of collimators
  integer,                  public,  save :: cdb_nfam     = 0       ! Number of collimator families

  ! Database arrays
  character(len=:), allocatable, public, save :: cdb_cName(:)     ! Collimator name
  character(len=:), allocatable, public, save :: cdb_cNameUC(:)   ! Collimator name upper case
  character(len=:), allocatable, public, save :: cdb_cMaterial(:) ! Collimator material
  integer,          allocatable, public, save :: cdb_cFamily(:)   ! Collimator family
  real(kind=fPrec), allocatable, public, save :: cdb_cNSig(:)     ! Collimator sigma
  real(kind=fPrec), allocatable, public, save :: cdb_cNSigOrig(:) ! Collimator sigma
  real(kind=fPrec), allocatable, public, save :: cdb_cLength(:)   ! Collimator length
  real(kind=fPrec), allocatable, public, save :: cdb_cOffset(:)   ! Collimator offset
  real(kind=fPrec), allocatable, public, save :: cdb_cRotation(:) ! Collimator rotation
  real(kind=fPrec), allocatable, public, save :: cdb_cBx(:)       ! Collimator beta x
  real(kind=fPrec), allocatable, public, save :: cdb_cBy(:)       ! Collimator beta y
  logical,          allocatable, public, save :: cdb_cFound(:)    ! Found in lattice

  ! Family Arrays
  character(len=:), allocatable, public, save :: cdb_famName(:)  ! Family name
  real(kind=fPrec), allocatable, public, save :: cdb_famNSig(:)  ! Family sigma

  ! Element Map
  integer,          allocatable, public, save :: cdb_elemMap(:)  ! Map from single elements to DB

contains

subroutine cdb_allocDB

  use mod_alloc
  use numerical_constants

  call alloc(cdb_cName,     mNameLen, cdb_nColl, " ",     "cdb_cName")
  call alloc(cdb_cNameUC,   mNameLen, cdb_nColl, " ",     "cdb_cNameUC")
  call alloc(cdb_cMaterial, 4,        cdb_nColl, " ",     "cdb_cMaterial")
  call alloc(cdb_cFamily,             cdb_nColl, -1,      "cdb_cFamily")
  call alloc(cdb_cNSig,               cdb_nColl, zero,    "cdb_cNSig")
  call alloc(cdb_cNSigOrig,           cdb_nColl, zero,    "cdb_cNSigOrig")
  call alloc(cdb_cLength,             cdb_nColl, zero,    "cdb_cLength")
  call alloc(cdb_cOffset,             cdb_nColl, zero,    "cdb_cOffset")
  call alloc(cdb_cRotation,           cdb_nColl, zero,    "cdb_cRotation")
  call alloc(cdb_cBx,                 cdb_nColl, zero,    "cdb_cBx")
  call alloc(cdb_cBy,                 cdb_nColl, zero,    "cdb_cBy")
  call alloc(cdb_cFound,              cdb_nColl, .false., "cdb_cFound")

end subroutine cdb_allocDB

subroutine cdb_allocFam

  use mod_alloc
  use numerical_constants

  call alloc(cdb_famName, cdb_fNameLen, cdb_nFam, " ",  "cdb_famName")
  call alloc(cdb_famNSig,               cdb_nFam, c1e3, "cdb_famNSig")

end subroutine cdb_allocFam

subroutine cdb_expand_arrays(nele_new)

  use mod_alloc

  integer, intent(in) :: nele_new

  call alloc(cdb_elemMap,nele_new,-1,"cdb_elemMap")

end subroutine cdb_expand_arrays

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-03-19
!  Updated: 2019-03-20
! ================================================================================================ !
subroutine cdb_readCollDB(dbFile)

  use crcoall
  use mod_units
  use mod_common
  use string_tools

  character(len=*), intent(in) :: dbFile

  character(len=:), allocatable :: lnSplit(:)
  character(len=mInputLn) inLine
  character(len=cdb_fNameLen) cFam
  real(kind=fPrec)  cNSig
  integer dbUnit, ioStat, nLines, nSplit, collID, famID, elemEnd
  integer i, j, ix
  logical dbErr, spErr

  call f_requestUnit(dbFile, dbUnit)
  call f_open(unit=dbUnit,file=dbFile,formatted=.true.,mode="r",err=dbErr,status="old")
  if(dbErr) then
    write(lout,"(a)") "COLL> ERROR Could not open the collimator database file '"//trim(dbFile)//"'"
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

#ifdef ROOT
  call cdb_writeDB_ROOT
#endif
#ifdef HDF5
  call cdb_writeDB_HDF5
#endif
  call cdb_writeFam

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

subroutine cdb_readDB_newFormat
end subroutine cdb_readDB_newFormat

subroutine cdb_readDB_oldFormat

  use crcoall
  use parpro
  use string_tools
  use mod_units
  use numerical_constants

  character(len=mInputLn) inLine
  character(len=cdb_fNameLen) famName
  logical cErr
  integer j, dbUnit, dbNew, ioStat, iLine, famID

  cErr = .false.

  call f_requestUnit(cdb_fileName, dbUnit)
  call f_open(unit=dbUnit,file=cdb_fileName,formatted=.true.,mode="r",status="old")

  call f_requestUnit(trim(cdb_fileName)//".new", dbNew)
  call f_open(unit=dbNew,file=trim(cdb_fileName)//".new",formatted=.true.,mode="w",status="replace")
  write(dbNew,"(1a,a47,1x,a16,1x,a4,5(1x,a13))") "#",chr_rPad(" name",47),&
    "opening","mat.","length[m]","angle[deg]","offset[m]","beta_x[m]","beta_y[m]"

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

    ! Line 2: Upper case name, ignored
    read(dbUnit,*,iostat=ioStat) cdb_cNameUC(j)
    iLine = iLine + 1
    if(ioStat /= 0) goto 100

    ! Line 3: Lower case name
    read(dbUnit,*,iostat=ioStat) cdb_cName(j)
    iLine = iLine + 1
    if(ioStat /= 0) goto 100

    ! Line 4: Sigma
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
    call cdb_getFamilyID(famName, famID, .false.)
    if(famID > 0 .and. cdb_doNSig) then
      cdb_cNSig(j) = cdb_famNSig(famID)
    else
      cdb_cNSig(j) = cdb_cNSigOrig(j)
    end if
    cdb_cFamily(j) = famID
    if(famID > 0) then
      write(dbNew,"(a48,1x,a16,1x,a4,5(1x,f13.6))") cdb_cName(j),&
      chr_lPad(trim(famName),16),cdb_cMaterial(j),cdb_cLength(j),cdb_cRotation(j)/rad,&
      cdb_cOffset(j),cdb_cBx(j),cdb_cBy(j)
    else
      write(dbNew,"(a48,1x,f16.3,1x,a4,5(1x,f13.6))") cdb_cName(j),&
      cdb_cNSig(j),cdb_cMaterial(j),cdb_cLength(j),cdb_cRotation(j)/rad,&
      cdb_cOffset(j),cdb_cBx(j),cdb_cBy(j)
    end if

  end do

  call f_close(dbUnit)

  flush(dbNew)
  call f_close(dbNew)

  return

100 continue
  write(lout,"(2(a,i0))") "COLLDB> ERROR Cannot read DB file line ",iLine,", iostat = ",ioStat
  call prror

end subroutine cdb_readDB_oldFormat

#ifdef ROOT
subroutine cdb_writeDB_ROOT

  use parpro
  use iso_c_binding
  use root_output

  character(len=mNameLen+1) :: this_name     = C_NULL_CHAR
  character(len=5)          :: this_material = C_NULL_CHAR
  integer j

  if(root_flag .eqv. .false. .or. root_CollimationDB /= 1) return

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
!  Created: 2019-03-21
!  Updated: 2019-03-21
!  Find a family in the database and returns its ID.
!  Optionally, if it doesn't exist, it can be added
! ================================================================================================ !
subroutine cdb_getFamilyID(famName, famID, addIfNew)

  character(len=*),  intent(in)  :: famName
  integer,           intent(out) :: famID
  logical, optional, intent(in)  :: addIfNew

  integer i
  logical addNew

  if(present(addIfNew)) then
    addNew = addIfNew
  else
    addNew = .false.
  end if

  famID = -1
  if(cdb_nfam > 0) then
    do i=1,cdb_nfam
      if(cdb_famName(i) == famName) then
        famID = i
        exit
      end if
    end do
  end if

  if(famID == -1 .and. addNew) then
    cdb_nfam = cdb_nfam + 1
    call cdb_allocFam
    cdb_famName(cdb_nfam) = famName
    famID = cdb_nfam
  end if

end subroutine cdb_getFamilyID

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-03-21
!  Updated: 2019-03-21
!  Set the nsig for a specific family
! ================================================================================================ !
subroutine cdb_setFamilyNSig(famID, nSig)

  use crcoall

  integer,          intent(in) :: famID
  real(kind=fPrec), intent(in) :: nSig

  if(famID < 1 .or. famID > cdb_nFam) then
    write(lout,"(a,i0)") "COLLDB> ERROR Setting nsig for non-existing famID = ",famID
    write(lout,"(a,i0)") "COLLDB>       First ID is 1, last ID is ",cdb_nFam
    call prror
  end if

  if(nSig < 0.0_fPrec) then
    write(lout,"(a,i0)") "COLLDB> ERROR Setting nsig for famID = ",famID
    write(lout,"(a)")    "COLLDB>       Value must be larger than zero"
    call prror
  end if

  cdb_famNSig(famID) = nSig

end subroutine cdb_setFamilyNSig

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-03-21
!  Updated: 2019-03-21
!  Write family database to file
! ================================================================================================ !
subroutine cdb_writeFam

  use mod_units

  integer j, famUnit

  call f_requestUnit("coll_families.dat", famUnit)
  call f_open(unit=famUnit,file="coll_families.dat",formatted=.true.,mode="w",status="replace")

  write(famUnit,"(a16,1x,a13)") "# famName       ","nSig"
  do j=1,cdb_nFam
    write(famUnit,"(a16,1x,f13.6)") cdb_famName(j),cdb_famNSig(j)
  end do

  flush(famUnit)
  call f_close(famUnit)

end subroutine cdb_writeFam

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
! else if(elemName(1:4) == "tcth" .or. elemName(1:5) == "tctxh" .or. elemName(1:5) == "tctph") then
  else if(elemName(1:4) == "tcth" .or. elemName(1:5) == "tctph") then
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
  else if(elemName(1:4) == "tcxr") then
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

end module collimation_db
