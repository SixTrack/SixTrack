! ================================================================================================ !
!  Machine Geometry Module
! ~~~~~~~~~~~~~~~~~~~~~~~~~
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-03-28
!  Updated: 2019-03-28
! ================================================================================================ !
module mod_geometry

  use floatPrecision

  implicit none

  character(len=2), parameter,   private       :: geom_go     = "GO"
  character(len=3), parameter,   private       :: geom_cavity = "CAV"

  integer,                       public,  save :: geom_nSing  = 0
  integer,                       public,  save :: geom_nBloc  = 0
  integer,                       public,  save :: geom_nStru  = 0
  integer,                       public,  save :: geom_nBeam  = 0

  character(len=:), allocatable, public,  save :: geom_bez0(:)
  character(len=:), allocatable, public,  save :: geom_beze(:,:)
  character(len=:), allocatable, public,  save :: geom_ilm(:)

  public :: geom_checkSingElemUnique
  public :: geom_insertStruElem
  public :: geom_insertSingElem

contains

! ================================================================================================ !
!  Parse Single Element Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-05-xx
! ================================================================================================ !
subroutine geom_parseInputLineSING(inLine, iLine, iErr)

  use parpro
  use crcoall
  use mod_alloc
  use mod_common
  use mod_settings
  use string_tools
  use sixtrack_input

  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=:), allocatable   :: elemName
  integer nSplit
  logical spErr

  integer i

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lerr,"(a)") "GEOMETRY> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(nSplit <= 2) then
    write(lerr,"(a,i0)") "GEOMETRY> ERROR Single element line must have more than 2 values, got ",nSplit
    iErr = .true.
    return
  end if

  elemName = lnSplit(1)
  if(len(elemName) > mNameLen) then
    write(lerr,"(a,i0)") "GEOMETRY> ERROR Single element name too long. Max length is ",mNameLen
    iErr = .true.
    return
  end if

  ! Check that the name is unique
  do i=1,geom_nSing-1
    if(bez(i) == elemName) then
      write(lerr,"(a,i0)") "GEOMETRY> ERROR Single element '"//trim(elemName)//"' is not unique."
      iErr = .true.
      return
    end if
  end do

  ! Save Values
  if(nSplit > 1) call chr_cast(lnSplit(2),kz(geom_nSing),  iErr)
  if(nSplit > 2) call chr_cast(lnSplit(3),ed(geom_nSing),  iErr)
  if(nSplit > 3) call chr_cast(lnSplit(4),ek(geom_nSing),  iErr)
  if(nSplit > 4) call chr_cast(lnSplit(5),el(geom_nSing),  iErr)
  if(nSplit > 5) call chr_cast(lnSplit(6),bbbx(geom_nSing),iErr)
  if(nSplit > 6) call chr_cast(lnSplit(7),bbby(geom_nSing),iErr)
  if(nSplit > 7) call chr_cast(lnSplit(8),bbbs(geom_nSing),iErr)

  if(iErr) return

  if(kz(geom_nSing) == 25) then
    ed(geom_nSing) = ed(geom_nSing)/two
    ek(geom_nSing) = ek(geom_nSing)/two
  end if

  ! CHANGING SIGN OF CURVATURE OF VERTICAL THICK DIPOLE
  if((kz(geom_nSing) == 4 .or. kz(geom_nSing) == 5) .and. abs(el(geom_nSing)) > pieni) then
    ed(geom_nSing) = -one*ed(geom_nSing)
  end if

  ! Beam-Beam Elements
  if(kz(geom_nSing) == 20) then
    geom_nBeam = geom_nBeam+1
    if(geom_nBeam > nbb) then
      call expand_arrays(nele, npart, nblz, nblo, nbb+100)
      if(st_debug) then
        write(lout,"(a,i0)") "GEOMETRY> Increased beam-beam element storage to nbb = ",nbb
      end if
    end if
  end if

  !--------------------------------------------
  ! Handled by initialise_element subroutine:
  !--------------------------------------------
  ! CHANGING SIGN OF CURVATURE OF VERTICAL THICK DIPOLE
  ! THIN LENS (+/- 1-10)
  ! MULTIPOLES (11)
  ! CAVITY (+/- 12)
  ! CRABCAVITY (23/-23) / CC multipoles order 2/3/4 (+/- 23/26/27/28)
  ! ELECTRON LENSE (29)
  call initialise_element(geom_nSing,.true.)

  ! ACDIPOLE
  if(abs(kz(geom_nSing)) == 16) then
    if(abs(ed(geom_nSing)) <= pieni) then
      kz(geom_nSing) = 0
      ed(geom_nSing) = zero
      ek(geom_nSing) = zero
      el(geom_nSing) = zero
    else
      acdipph(geom_nSing) = el(geom_nSing)
      el(geom_nSing)      = zero
    end if
  end if

  ! General
  if(abs(el(geom_nSing)) > pieni .and. kz(geom_nSing) /= 0) then
    ithick = 1
  end if

  ! Expand Arrays
  if(geom_nSing > nele-2) then
    call expand_arrays(nele+100, npart, nblz, nblo, nbb)
    call alloc(geom_bez0, mNameLen, nele, " ", "geom_bez0")
  end if

  if(abs(kz(geom_nSing)) /= 12 .or. (abs(kz(geom_nSing)) == 12 .and. ncy2 == 0)) then
    kp(geom_nSing) = 0
  end if

  bez(geom_nSing)       = elemName
  geom_bez0(geom_nSing) = elemName

  !If no active RF cavities are seen so far in the single element list,
  ! add a CAV element to the end of the list.
  ! This is then overwritten when reading the next element, so that if
  ! and only if no active RF cavities are found, a CAV element can be
  ! used in the structure to enable 6D tracking using the parameters
  ! from the SYNC block.
  if(ncy2 == 0) then
    geom_nSing            = geom_nSing + 1
    il                    = geom_nSing
    bez(geom_nSing)       = geom_cavity
    geom_bez0(geom_nSing) = geom_cavity
    kp(geom_nSing)        = 6
  else
    il         = geom_nSing
    geom_nSing = geom_nSing + 1
  end if

end subroutine geom_parseInputLineSING

! ================================================================================================ !
!  Parse Block Definitions Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2018-05-xx
! ================================================================================================ !
subroutine geom_parseInputLineBLOC(inLine, iLine, iErr)

  use parpro
  use crcoall
  use mod_alloc
  use mod_common
  use string_tools

  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=:), allocatable   :: blocName
  integer nSplit
  logical spErr

  integer i, j, ka, ke, nInd
  logical eFound, isCont
  character(len=mNameLen) ilm0(40)

  integer :: k0 = 0

  call chr_split(inLine, lnSplit, nSplit, spErr, nIndent=nInd)
  if(spErr) then
    write(lerr,"(a)") "GEOMETRY> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if
  isCont = (nInd >= 5)

  if(nSplit < 2 .and. .not. isCont) then
    write(lerr,"(a,i0)") "GEOMETRY> ERROR Block definition line must be at least 2 values, got ",nSplit
    iErr = .true.
    return
  end if

  if(nSplit > 40) then
    write(lerr,"(a)") "GEOMETRY> ERROR Block definition cannot have more then 40 elements."
    iErr = .true.
    return
  end if

  ! If first line, read super period information
  if(iLine == 1) then
    call chr_cast(lnSplit(1),mper,iErr)
    if(mper > nper) then
      write(lerr,"(a,i0)") "GEOMETRY> ERROR Block definition number of super periods is too large. Max value is ",nper
      iErr = .true.
      return
    end if
    if(mper > nSplit+1) then
      write(lerr,"(a)") "GEOMETRY> ERROR Block definition number of super periods does not match the number of values."
      iErr = .true.
      return
    end if
    do i=1,mper
      call chr_cast(lnSplit(i+1),msym(i),iErr)
    end do
    if(iErr) return

    ! Init variables
    call alloc(geom_ilm,  mNameLen,       nelb, " ", "geom_ilm")
    call alloc(geom_beze, mNameLen, nblo, nelb, " ", "geom_beze")

    ! No need to parse anything more for this line
    return
  end if

  ! Parse normal line, iLine > 1
  ilm0(:) = " "

  if(isCont) then                             ! This line continues the previous BLOC
    blocName = " "                            ! No name returned, set an empty BLOC name
    ilm0(1:nSplit) = lnSplit(1:nSplit)        ! All elements are sub-elements. Save to buffer.
  else                                        ! This is a new BLOC
    blocName = lnSplit(1)
    ilm0(1:nSplit-1) = lnSplit(2:nSplit)      ! Shift the buffer
  end if

  if(blocName /= " ") then                    ! We have a new BLOC
    geom_nBloc = geom_nBloc + 1               ! Increment the BLOC number
    if(geom_nBloc > nblo-1) then              ! Expand arrays if needed
      call expand_arrays(nele, npart, nblz, nblo+50, nbb)
      call alloc(geom_beze, mNameLen, nblo, nelb, " ", "geom_beze")
    end if
    bezb(geom_nBloc) = blocName               ! Set the BLOC name in bezb
    k0               = 0                      ! Reset the single element counter
    mblo             = geom_nBloc             ! Update total number of BLOCs
  end if

  ka = k0 + 1
  ke = k0 + 40

  do i=ka, ke
    if(i > nelb) then
      write(lerr,"(a,2(i0,a))") "GEOMETRY> ERROR Block definitions can only have ",&
        nelb," elements. ",i," given."
      iErr = .true.
      return
    end if
    geom_ilm(i) = ilm0(i-k0)                  ! Append to sub-element buffer
    if(geom_ilm(i) == " ") exit               ! No more sub-elements to append
    mel(geom_nBloc)         = i               ! Update number of single elements in this block
    geom_beze(geom_nBloc,i) = geom_ilm(i)     ! Name of the current single element

    eFound = .false.
    do j=1,il                                 ! Search for the single element index
      if(geom_bez0(j) == geom_ilm(i)) then
        eFound = .true.
        exit
      end if
    end do
    if(eFound) then                           ! Handle element found
      mtyp(geom_nBloc,i) = j                  ! Save single element index
      if(kz(j) /= 8) then                     ! Count block length (kz=8 - edge focusing: skip)
        elbe(geom_nBloc) = elbe(geom_nBloc) + el(j)
      end if
    else                                      ! If not, the input is invalid
      write(lerr,"(a)") "GEOMETRY> ERROR Unknown single element '"//&
        chr_strip(geom_ilm(i))//"' in block definitions."
      iErr = .true.
      return
    end if
  end do

  k0 = i-1

end subroutine geom_parseInputLineBLOC

! ================================================================================================ !
!  Parse Structure Input Line
!  Rewritten from code from DATEN by VKBO
!  Last modified: 2019-03-27
! ================================================================================================ !
subroutine geom_parseInputLineSTRU(inLine, iLine, iErr)

  use parpro
  use crcoall
  use mod_common
  use string_tools
  use sixtrack_input

  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=:), allocatable   :: expLine
  integer nSplit
  logical spErr

  integer i, j

  expLine = chr_expandBrackets(inLine)
  call chr_split(expLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lerr,"(a)") "GEOMETRY> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(nSplit < 1 .or. nSplit > 40) then
    write(lerr,"(a)") "GEOMETRY> ERROR Structure input line cannot have less than 1 or more then 40 elements."
    iErr = .true.
    return
  end if

  do i=1,nSplit

    if(len_trim(lnSplit(1)) > mNameLen) then
      write(lerr,"(a)") "GEOMETRY> ERROR Structure element name cannot have more than ",mNameLen," characters."
      iErr = .true.
      return
    end if

    if(lnSplit(i) == " ") cycle
    if(lnSplit(i) == geom_go) then
      kanf = geom_nStru + 1
      cycle
    end if

    geom_nStru = geom_nStru + 1
    if(geom_nStru > nblz-3) then
      call expand_arrays(nele, npart, nblz+1000, nblo, nbb)
    end if

    do j=1,mblo ! is it a BLOC?
      if(bezb(j) == lnSplit(i)) then
        ic(geom_nStru)   = j
        bezs(geom_nStru) = bezb(j)
        exit
      end if
    end do

    do j=1,il ! is it a SINGLE ELEMENT?
      if(geom_bez0(j) == lnSplit(i)) then
        ic(geom_nStru)   = j+nblo
        bezs(geom_nStru) = geom_bez0(j)
        if(geom_bez0(j) == geom_cavity) then
          icy = icy+1
        end if
        exit
      end if
    end do
  end do

  mbloz = geom_nStru
  if(mbloz > nblz-3) then
    call expand_arrays(nele, npart, nblz+1000, nblo, nbb)
  end if

end subroutine geom_parseInputLineSTRU

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Parse Structure Input Line - Multi-Column Version
!  Created: 2019-03-27
!  Updated: 2019-03-28
! ================================================================================================ !
subroutine geom_parseInputLineSTRU_MULT(inLine, iLine, iErr)

  use parpro
  use crcoall
  use mod_common
  use string_tools
  use sixtrack_input

  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit, singID
  logical spErr, cErr

  integer j

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lerr,"(a)") "GEOMETRY> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(nSplit == 0) return

  if(lnSplit(1) == geom_go) then
    kanf = geom_nStru + 1
    return
  end if

  if(lnSplit(1) == "MULTICOL") then
    return
  end if

  if(nSplit < 3) then
    write(lerr,"(a)") "GEOMETRY> ERROR Multi-column structure input line cannot have less than 3 elements."
    iErr = .true.
    return
  end if

  if(len_trim(lnSplit(1)) > mNameLen .or. len_trim(lnSplit(2)) > mNameLen) then
    write(lerr,"(a)") "GEOMETRY> ERROR Structure element names in columns 1 and 2 cannot have more than ",mNameLen," characters."
    iErr = .true.
    return
  end if

  geom_nStru = geom_nStru + 1
  if(geom_nStru > nblz-3) then
    call expand_arrays(nele, npart, nblz+1000, nblo, nbb)
  end if

  bezs(geom_nStru) = trim(lnSplit(1))
  call chr_cast(lnSplit(3), elpos(geom_nStru), cErr)
  if(elpos(geom_nStru) < elpos(geom_nStru-1)) then ! Note: elpos(0) does exist, and should be zero
    write(lerr,"(a)") "GEOMETRY> ERROR Structure element '"//trim(bezs(geom_nStru))//&
      "' cannot be positioned ahead of previous element"
    iErr = .true.
    return
  end if

  singID = -1
  do j=1,mblo ! is it a BLOC?
    if(bezb(j) == lnSplit(2)) then
      singID = j
      exit
    end if
  end do

  if(singID == -1) then
    do j=1,il ! is it a SINGLE ELEMENT?
      if(geom_bez0(j) == lnSplit(2)) then
        singID = j+nblo
        if(geom_bez0(j) == geom_cavity) then
          icy = icy+1
        end if
        exit
      end if
    end do
  end if

  if(singID == -1) then
    write(lerr,"(a)") "GEOMETRY> ERROR Unknown element '"//trim(lnSplit(2))//"' in STRUCTURE INPUT."
    iErr = .true.
    return
  else
    ic(geom_nStru) = singID
  end if

  mbloz = geom_nStru
  if(mbloz > nblz-3) then
    call expand_arrays(nele, npart, nblz+1000, nblo, nbb)
  end if

end subroutine geom_parseInputLineSTRU_MULT

! ================================================================================================ !
!  A.Mereghetti, V.K. Berglyd Olsen, BE-ABP-HSS
!  Insert a New Empty Element (empty) in SINGLE ELEMENTS
!  Updated: 2019-03-28
! ================================================================================================ !
integer function geom_insertSingElem()

  use parpro
  use mod_common, only : il, ithick

  il = il + 1
  if(il > nele-2) then
    call expand_arrays(nele+50, npart, nblz, nblo, nbb)
    if(ithick == 1) then
      call expand_thickarrays(nele, npart, nblz, nblo )
    end if
  end if
  geom_insertSingElem = il

end function geom_insertSingElem

! ================================================================================================ !
!  A.Mereghetti, V.K. Berglyd Olsen, BE-ABP-HSS
!  Insert a New Empty Element (empty) in STRUCTURE ELEMENTS
!  Updated: 2019-03-28
!  iEl : Index in lattice to insert new element. 0 = append, negative count from last index
! ================================================================================================ !
integer function geom_insertStruElem(iEl)

  use parpro
#ifdef DEBUG
  use crcoall
#endif
  use mod_common

  integer, intent(in) :: iEl

  integer i,iIns
  character(len=mNameLen) tmpC

  if(iu > nblz-3) then
    call expand_arrays(nele, npart, nblz+100, nblo, nbb)
  end if

  iu = iu + 1
  if(iEl <= 0) then
    iIns = iu+iEl
  else
    iIns = iEl
  end if

#ifdef DEBUG
  write(lout,*) "geom_insertStruElem: iEl,iu,iIns:",iEl,iu,iIns
#endif

  ic(iIns:iu)      = cshift(ic(iIns:iu),      -1)
  elpos(iIns:iu)   = cshift(elpos(iIns:iu),   -1)
! bezs(iIns:iu)    = cshift(bezs(iIns:iu),    -1)
  icext(iIns:iu)   = cshift(icext(iIns:iu),   -1)
  icextal(iIns:iu) = cshift(icextal(iIns:iu), -1)
  dcum(iIns:iu)    = cshift(dcum(iIns:iu),    -1)

  ! Do string arrays with a loop due to a gfrotran bug in at least 8.3
  tmpC = bezs(iu)
  do i=iu-1,iIns,-1
    bezs(i+1) = bezs(i)
  end do
  bezs(iIns) = tmpC

  ! Update s coordinate of added element
  elpos(iIns) = elpos(iIns-1)
  dcum(iIns)  = dcum(iIns-1)

  geom_insertStruElem = iu

end function geom_insertStruElem

! ================================================================================================ !
!  A.Mereghetti, V.K. Berglyd Olsen, BE-ABP-HSS
!  Check that a given entry in the sequence is unique
!  Updated: 2019-03-28
!  iEl  : Index in lattice structure to be checked
!  ixEl : Index in array of SINGLE ELEMENTs of the element to be checked
! ================================================================================================ !
integer function geom_checkSingElemUnique(iEl, ixEl)

  use parpro
  use mod_common, only : iu, ic

  integer, intent(in) :: iEl
  integer, intent(in) :: ixEl

  integer i, ix

  geom_checkSingElemUnique = -1

  do i=1,iu
    ix = ic(i)-nblo
    if(ix > 0) then ! SINGLE ELEMENT
      if(i /= iEl .and. ix == ixEl) then
        geom_checkSingElemUnique = i
        return
      end if
    end if
  end do

end function geom_checkSingElemUnique

! ================================================================================================ !
!  A.Mereghetti, V.K. Berglyd Olsen, BE-ABP-HSS
!  Check that a given entry in the sequence is unique
!  Updated: 2019-03-28
!  sLoc     : S-coordinate where element should be
!  isLast   : True if return last lens at sLoc
!  iEl      : Index in lattice structure of found element
!  ixEl     : Index in array of SINGLE ELEMENTs of found element
!  wasFound : True if element was found
! ================================================================================================ !
subroutine geom_findElemAtLoc(sLoc, isLast, iEl, ixEl, wasFound)

  use parpro
  use crcoall
  use floatPrecision
  use mod_common, only : iu, tlen, ic, dcum
  use numerical_constants, only : zero

  real(kind=fPrec), intent(out) :: sLoc
  logical,          intent(in)  :: isLast
  integer,          intent(out) :: iEl
  integer,          intent(out) :: ixEl
  logical,          intent(out) :: wasFound

  integer i, iDelta, iCheck, iMax, iStep
  logical lSlide

  iEl  = -1
  ixEl = -1
  sLoc = zero

  if(sLoc > tlen .or. sLoc < zero) then
    write(lerr,"(a,2(f11.4,a))") "GEOMETRY> ERROR Find Element: "//&
      "Requested s-location: ",sLoc," is not inside accelerator range [0:",tlen,"]"
    call prror
  endif

  ! Fast search
  iCheck = iu
  iDelta = iu
  do while(iDelta > 1 .or. iCheck > 0 .or. iCheck < iu)
    if(dcum(iCheck) == sLoc) exit
    iDelta = nint(real(iDelta/2))
    if(dcum(iCheck) < sLoc) then
      iCheck = iCheck + iDelta
    else
      iCheck = iCheck - iDelta
    end if
  end do

  ! Finalise search
  if(dcum(iCheck) < sLoc .or. dcum(iCheck) == sLoc .and. isLast) then
    iMax   = iu
    iStep  = 1
    lslide = isLast
  else
    iMax   = 1
    iStep  = -1
    lSlide = .not.isLast
  endif
  do i=iCheck,iMax,iStep
    if(dcum(i) < sLoc) cycle
    if(dcum(i) > sLoc) then
      if(iEl == -1) iEl = i
    else
      iEl = i
      wasFound = .true.
      if(lSlide) cycle
    end if
    exit
  end do

  if(wasFound) then
    ixEl = ic(iEl)-nblo
    if(ixEl < 0) ixEl=ic(iEl) ! drift
  else
    write(lout,"(a,2(f11.4,a))") "GEOMETRY> Find Element: s-location: ",sLoc," was not found in acclerator range [0:",tlen,"]"
  end if

end subroutine geom_findElemAtLoc

! ================================================================================================ !
!  A. Mereghetti, D. Sinuela Pastor, for the FLUKA Team
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Original: 2014-06-13
!  Updated:  2019-03-29
!  Calculate dcum, as done in linopt and when parsing BLOCs (daten):
!    lengths of thick lens elements are taken on the curvilinear
!    reference system; thus, no difference between the length
!    of SBENDs and the one of RBENDs, as they are both the ARC one;
!  For future needs:  ds=two/ed(ix)*asin(el(ix)*ed(ix)/two)
! ================================================================================================ !
subroutine geom_calcDcum

  use parpro
  use crcoall
  use mod_meta
  use mod_units
  use mod_common
  use string_tools
  use floatPrecision
  use numerical_constants

  character(len=24) tmpS, tmpE
  character(len=99) fmtH, fmtC
  real(kind=fPrec) tmpDcum, delS, sGo, sEnd
  integer i, j, k, ix, outUnit

  write(lout,"(a)") "GEOMETRY> Calculating machine length"

  tmpDcum = zero

  ! Loop all over the entries in the accelerator structure
  do i=1,iu
    ix = ic(i)
    if(ix > nblo) then ! SINGLE ELEMENT
      ix = ix-nblo
      if(el(ix) > zero) then
        tmpDcum = tmpDcum  + el(ix)
      end if
      if(strumcol) then
        elpos(i) = elpos(i) + el(ix)/2 ! Change from centre of element to end of element
      else
        elpos(i) = tmpDcum ! Just copy dcum(i)
      end if
    else ! BLOC
      ! Iterate over elements
      do j=1,mel(ix)
        k = mtyp(ix,j)
        if(el(k) > zero) then
          tmpDcum = tmpDcum + el(k)
        end if
        if(strumcol) then
          elpos(i) = elpos(i) + el(k)/2 ! Change from centre of element to end of element
        else
          elpos(i) = tmpDcum ! Just copy dcum(i)
        end if
      end do
    end if
    dcum(i) = tmpDcum
  end do

  ! Correct the elpos array after a GO reshuffle
  if(kanf > 1 .and. strumcol) then
    sGo  = elpos(1)
    sEnd = elpos(iu-kanf+1)-sGo
    elpos(1:iu-kanf+1)  = elpos(1:iu-kanf+1)  - sGo
    elpos(iu-kanf+2:iu) = elpos(iu-kanf+2:iu) + sEnd
  end if

  ! Assign the last value to the closing MARKER:
  dcum(iu+1)  = tmpDcum
  elpos(iu+1) = elpos(iu)

  write(lout,"(a,f17.10,a)") "GEOMETRY> Machine length was ",dcum(iu+1)," [m]"

  call meta_write("MultiColStructBlock", strumcol)
  call meta_write("MachineLengthInput",  tlen)
  call meta_write("MachineLengthAccum",  dcum(iu+1))
  call meta_write("MachineLengthStruct", elpos(iu+1))

  if(idp /= 0.and. ition /= 0) then ! 6D tracking
    if(abs(dcum(iu+1) - tlen) > eps_dcum) then
      write(lout,"(a)")          ""
      write(lout,"(a)")          "    WARNING Problem with SYNC block detected"
      write(lout,"(a,f17.10)")   "            TLEN in SYNC block = ",tlen
      write(lout,"(a,f17.10)")   "            Length from DCUM   = ",dcum(iu+1)
      write(lout,"(a,f17.10)")   "            Difference         = ",dcum(iu+1)-tlen
      write(lout,"(a,e27.16,a)") "            Relative error     = ",2*(dcum(iu+1)-tlen)/(dcum(iu+1)+tlen)," [m]"
      write(lout,"(a,f17.10,a)") "            Tolerance eps_dcum = ",eps_dcum," [m]"
      write(lout,"(a)")          "    Please fix the TLEN parameter in your SYNC block"
      write(lout,"(a)")          "    so that it matches the calculated machine length from DCUM."
      write(lout,"(a)")          "    If incorrect, the RF frequency may be (slightly) wrong."
      write(lout,"(a)")          ""
      write(lout,"(a)")          str_divLine
      ! It's a warning not an error, and the consequences seem relatively small.
      ! Ideally, tlen should be calculated automatically based on the sequence.
    end if
  else
    tlen = dcum(iu+1)
  endif

  ! Enabled by the PRINT_DCUM flag in the SETTINGS block
  if(print_dcum) then

    call f_requestUnit("machine_length.dat",outUnit)
    call f_open(unit=outUnit,file="machine_length.dat",formatted=.true.,mode="w",status="replace")

    tmpS = "START"
    tmpE = "END"
    fmtH = "(a1,1x,a29,1x,a31,2(1x,a17),1x,a10)"
    fmtC = "(i6,1x,a24,1x,i6,1x,a24,2(1x,f17.9),1x,f10.3)"

    if(strumcol) then
      write(outUnit,"(a)") "# Structure Element MULTICOL flag is ON. s_madx column is read from input."
    else
      write(outUnit,"(a)") "# Structure Element MULTICOL flag is OFF. s_madx column is a copy of s_dcum."
    end if

    write(outUnit,fmtH) "#",chr_rPad("idST StructureElement",29),chr_rPad("idSE SingleElement",29),&
      "s_dcum[m]","s_madx[m]","delta[nm]"
    delS = (elpos(0) - dcum(0))*1e9
    write(outUnit,fmtC) 0,tmpS,-1,tmpS,dcum(0),elpos(0),delS
    do i=1,iu
      ix   = ic(i)
      delS = (elpos(i) - dcum(i))*1e9
      if(ix > nblo) then ! SINGLE ELEMENT
        ix = ix-nblo
        write(outUnit,fmtC) i,bezs(i)(1:24),ix,bez(ix)(1:24),dcum(i),elpos(i),delS
      else ! BLOC
        write(outUnit,fmtC) i,bezs(i)(1:24),ix,bezb(ix)(1:24),dcum(i),elpos(i),delS
      end if
    end do
    delS = (elpos(iu+1) - dcum(iu+1))*1e9
    write(outUnit,fmtC) iu+1,tmpE,-1,tmpE,dcum(iu+1),elpos(iu+1),delS

    flush(outUnit)
    call f_freeUnit(outUnit)
  end if

end subroutine geom_calcDcum

! ================================================================================================ !
!  A. Mereghetti, V.K. Berglyd Olsen
!  Original:  2016-03-14 (orglat)
!  Rewritten: 2019-04-03
!  Updated:   2019-04-03
!  Reshuffle the lattice
! ================================================================================================ !
subroutine geom_reshuffleLattice

  use parpro
  use crcoall
  use floatPrecision
  ! use numerical_constants
  use mod_common
  use mod_commons
  use mod_common_track

  implicit none

  integer i,j,ii,jj
  character(len=mNameLen), allocatable :: tmpC(:)

  if(mper > 1) then
    do i=2,mper
      do j=1,mbloz
        ii = (i-1)*mbloz+j
        jj = j
        if(msym(i) < 0) then
          jj = mbloz-j+1
        end if
        ic(ii) = msym(i)*ic(jj)
        if(ic(ii) < -nblo) then
          ic(ii) = -ic(ii)
        end if
      end do
    end do
  end if
  iu = mper*mbloz

  ! "GO" is the first structure element, we're done.
  if(kanf == 1) return

  ! Otherwise, reshuffle the structure
  write(lout,"(a)") "GEOMETRY> GO keyword not the first lattice element: Reshuffling lattice structure."

  if(iorg >= 0) then
    mzu(1:iu)   = cshift(mzu(1:iu),     kanf-1)
  end if
  ic(1:iu)      = cshift(ic(1:iu),      kanf-1)
  icext(1:iu)   = cshift(icext(1:iu),   kanf-1)
  icextal(1:iu) = cshift(icextal(1:iu), kanf-1)
! bezs(1:iu)    = cshift(bezs(1:iu),    kanf-1)
  elpos(1:iu)   = cshift(elpos(1:iu),   kanf-1)

  ! Do string arrays manually due to a gfortran bug in at least 8.3
  allocate(tmpC(kanf))
  tmpC(1:kanf-1)     = bezs(1:kanf-1)
  bezs(1:iu-kanf+1)  = bezs(kanf:iu)
  bezs(iu-kanf+2:iu) = tmpC(1:kanf-1)

end subroutine geom_reshuffleLattice

end module mod_geometry

! ================================================================================================ !
!  K. Sjobak, A. Santamaria, V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2016-12-23
!  Updated: 2019-08-29
!
!  Initialise a lattice element with index elIdx, such as done when reading fort.2 and in DYNK.
!
!  Never delete an element from the lattice, even if it is not making a kick.
!  If the element is not recognized, do nothing (for now).
!  If trying to initialize an element (not lfirst) which is disabled, print an error and exit.
!
!  This subroutine belongs with mod_geometry, but is outside the module due to circular dependency
!  with the DYNK module. If further issues with dependencies are encountered, it can also be moved
!  to the mainarrays.f90 file.
! ================================================================================================ !
subroutine initialise_element(ix,lfirst)

  use crcoall
  use floatPrecision
  use mathlib_bouncer
  use numerical_constants

  use parpro
  use parbeam
  use tracking
  use mod_common
  use mod_common_main
  use mod_common_track
  use mod_utils

  use cheby, only : cheby_kz
  use dynk,  only : dynk_elemData

  implicit none

  integer, intent(in) :: ix
  logical, intent(in) :: lfirst

  integer i,m,k,im,nmz,izu,ibb,ii,j
  real(kind=fPrec) r0,r0a,bkitemp,sfac1,sfac2,sfac2s,sfac3,sfac4,sfac5,crkveb_d,cikveb_d,rho2b_d,   &
    tkb_d,r2b_d,rb_d,rkb_d,xrb_d,zrb_d,cbxb_d,cbzb_d,crxb_d,crzb_d,xbb_d,zbb_d,napx0
  real(kind=fPrec) crkveb(npart),cikveb(npart),rho2b(npart),tkb(npart),r2b(npart),rb(npart),        &
    rkb(npart),xrb(npart),zrb(npart),xbb(npart),zbb(npart),crxb(npart),crzb(npart),cbxb(npart),     &
    cbzb(npart)

  ! Nonlinear Elements
  if(abs(kz(ix)) >= 1 .and. abs(kz(ix)) <= 10) then
    if(.not.lfirst) then
      do i=1,iu
        if(ic(i)-nblo == ix) then
          if(ktrack(i) == 31) goto 100 !ERROR
          sm(ix)  = ed(ix)          ! Also done in envar() which is called from clorb()
          smiv(i) = sm(ix)+smizf(i) ! Also done in program maincr
          smi(i)  = smiv(i)         ! Also done in program maincr
          select case(abs(kz(ix)))
          case(1)
            call setStrack(1,i)
          case(2)
            call setStrack(2,i)
          case(3)
            call setStrack(3,i)
          case(4)
            call setStrack(4,i)
          case(5)
            call setStrack(5,i)
          case(6)
            call setStrack(6,i)
          case(7)
            call setStrack(7,i)
          case(8)
            call setStrack(8,i)
          case(9)
            call setStrack(9,i)
          case(10)
            call setStrack(10,i)
          end select
        end if
      end do
    end if

  ! Multipoles
  elseif(kz(ix) == 11) then
    if(lfirst) then
      if(abs(el(ix)+one) <= pieni) then
        dki(ix,1) = ed(ix)
        dki(ix,3) = ek(ix)
        ed(ix) = one
        ek(ix) = one
        el(ix) = zero
      else if(abs(el(ix)+two) <= pieni) then
        dki(ix,2) = ed(ix)
        dki(ix,3) = ek(ix)
        ed(ix) = one
        ek(ix) = one
        el(ix) = zero
      end if
    else
      do i=1,iu
        if(ic(i)-nblo == ix) then
          nmz = nmu(ix)
          im  = irm(ix)
          do k=1,nmz
            aaiv(k,i) = scalemu(im)*(ak0(im,k)+amultip(k,i)*aka(im,k))
            bbiv(k,i) = scalemu(im)*(bk0(im,k)+bmultip(k,i)*bka(im,k))
          end do
        end if
      end do
    end if

  ! Cavities (ktrack = 2 for thin)
  elseif(abs(kz(ix)) == 12) then
    dynk_elemData(ix,3) = el(ix)
    phasc(ix) = el(ix)*rad
    el(ix) = zero
    if(lfirst) then
      if(abs(ed(ix)) > pieni .and. abs(ek(ix)) > pieni) then
        ncy2   = ncy2 + 1
        kp(ix) = 6
      end if
    else
      hsyc(ix) = ((twopi)*ek(ix))/tlen                             ! SYNC block
      hsyc(ix) = (c1m3*hsyc(ix)) * real(sign(1,kz(ix)),kind=fPrec) ! trauthin/trauthck
    end if

  ! Wire
  else if(kz(ix) == 15) then
    ed(ix) = zero
    ek(ix) = zero
    el(ix) = zero

  ! BEAM-BEAM
  elseif(kz(ix) == 20) then

    if(nbeam == 0 .and. .not. lfirst) then
      write(lerr,"(a)") "BEAMBEAM> ERROR Beam-beam element encountered, but no BEAM block in '"//trim(fort3)//"'"
      call prror
    end if

    if(lfirst) then
      ptnfac(ix)  = el(ix)
      el(ix)      = zero
      parbe(ix,5) = ed(ix)
      ed(ix)      = zero
      parbe(ix,6) = ek(ix)
      ek(ix)      = zero
    end if
    ! This is to inialize all the beam-beam element before the tracking (or to update it for DYNK).
    if(.not.lfirst) then
      do i=1,iu
        if(ic(i)-nblo == ix) then
          ibb = imbb(i)
          if(parbe(ix,2) > zero) then
            if(beam_expflag == 1) then
              bbcu(ibb,1)  = parbe(ix,7)
              bbcu(ibb,4)  = parbe(ix,8)
              bbcu(ibb,6)  = parbe(ix,9)
              bbcu(ibb,2)  = parbe(ix,10)
              bbcu(ibb,9)  = parbe(ix,11)
              bbcu(ibb,10) = parbe(ix,12)
              bbcu(ibb,3)  = parbe(ix,13)
              bbcu(ibb,5)  = parbe(ix,14)
              bbcu(ibb,7)  = parbe(ix,15)
              bbcu(ibb,8)  = parbe(ix,16)
              do ii=1,10
                bbcu(ibb,ii) = bbcu(ibb,ii)*c1m6
              end do
            end if
            ktrack(i)   = 44
            parbe(ix,4) = (((-one*crad)*ptnfac(ix))*half)*c1m6
            if(ibeco == 1) then
              track6d(1,1) = parbe(ix,5)*c1m3
              track6d(2,1) = zero
              track6d(3,1) = parbe(ix,6)*c1m3
              track6d(4,1) = zero
              track6d(5,1) = zero
              track6d(6,1) = zero
              napx0 = napx
              napx  = 1
              call beamint(napx,track6d,parbe,sigz,bbcu,imbb(i),ix,ibtyp,ibbc, mtc)
              beamoff(1,imbb(i)) = track6d(1,1)*c1e3
              beamoff(2,imbb(i)) = track6d(3,1)*c1e3
              beamoff(3,imbb(i)) = track6d(5,1)*c1e3
              beamoff(4,imbb(i)) = track6d(2,1)*c1e3
              beamoff(5,imbb(i)) = track6d(4,1)*c1e3
              beamoff(6,imbb(i)) = track6d(6,1)
              napx = napx0
            end if

          else if(parbe(ix,2) == zero) then
            if(beam_expflag == 1) then
              bbcu(ibb,1) = parbe(ix,1)
              bbcu(ibb,2) = parbe(ix,3)
              bbcu(ibb,3) = parbe(ix,13)
            end if
            if(ibbc == 1) then
              sfac1  = bbcu(ibb,1)+bbcu(ibb,2)
              sfac2  = bbcu(ibb,1)-bbcu(ibb,2)
              sfac2s = one
              if(sfac2 < zero) sfac2s = -one
              sfac3 = sqrt(sfac2**2+(four*bbcu(ibb,3))*bbcu(ibb,3))
              if(sfac3 > sfac1) then
                write(lerr,"(a)") "BEAMBEAM> ERROR 6D beam-beam with tilt not possible."
                call prror
              end if
              sfac4 = (sfac2s*sfac2)/sfac3
              sfac5 = (((-one*sfac2s)*two)*bbcu(ibb,3))/sfac3
              sigman(1,ibb) = sqrt(((sfac1+sfac2*sfac4)+(two*bbcu(ibb,3))*sfac5)*half)
              sigman(2,ibb) = sqrt(((sfac1-sfac2*sfac4)-(two*bbcu(ibb,3))*sfac5)*half)
              bbcu(ibb,11)  = sqrt(half*(one+sfac4))
              bbcu(ibb,12)  = (-one*sfac2s)*sqrt(half*(one-sfac4))
              if(bbcu(ibb,3) < zero) bbcu(ibb,12) = -one*bbcu(ibb,12)
            else
              bbcu(ibb,11)  = one
              sigman(1,ibb) = sqrt(bbcu(ibb,1))
              sigman(2,ibb) = sqrt(bbcu(ibb,2))
            end if

            ! Round beam
            nbeaux(imbb(i)) = 0
            if(sigman(1,imbb(i)) == sigman(2,imbb(i))) then
              if(nbeaux(imbb(i)) == 2 .or. nbeaux(imbb(i)) == 3) then
                write(lerr,"(a)") "BEAMBEAM> ERROR At each interaction point the beam must be either "//&
                  "round or elliptical for all particles"
                call prror
              else
                nbeaux(imbb(i)) = 1
                sigman2(1,imbb(i)) = sigman(1,imbb(i))**2
              end if
            end if

            ! Elliptic beam x>z
            if(sigman(1,imbb(i)) > sigman(2,imbb(i))) then
              if(nbeaux(imbb(i)) == 1 .or. nbeaux(imbb(i)) == 3) then
                write(lerr,"(a)") "BEAMBEAM> ERROR At each interaction point the beam must be either "//&
                  "round or elliptical for all particles"
                call prror
              else
                nbeaux(imbb(i)) = 2
                ktrack(i)       = 42
                sigman2(1,imbb(i)) = sigman(1,imbb(i))**2
                sigman2(2,imbb(i)) = sigman(2,imbb(i))**2
                sigmanq(1,imbb(i)) = sigman(1,imbb(i))/sigman(2,imbb(i))
                sigmanq(2,imbb(i)) = sigman(2,imbb(i))/sigman(1,imbb(i))
              end if
            end if

            ! Elliptic beam z>x
            if(sigman(1,imbb(i)) < sigman(2,imbb(i))) then
              if(nbeaux(imbb(i)) == 1 .or. nbeaux(imbb(i)) == 2) then
                write(lerr,"(a)") "BEAMBEAM> ERROR At each interaction point the beam must be either "//&
                  "round or elliptical for all particles"
                call prror
              else
                nbeaux(imbb(i)) = 3
                ktrack(i)       = 43
                sigman2(1,imbb(i)) = sigman(1,imbb(i))**2
                sigman2(2,imbb(i)) = sigman(2,imbb(i))**2
                sigmanq(1,imbb(i)) = sigman(1,imbb(i))/sigman(2,imbb(i))
                sigmanq(2,imbb(i)) = sigman(2,imbb(i))/sigman(1,imbb(i))
              end if
            end if

            strack(i) = crad*ptnfac(ix)
            if(ibbc.eq.0) then
              crkveb_d = parbe(ix,5)
              cikveb_d = parbe(ix,6)
            else
              crkveb_d = parbe(ix,5)*bbcu(imbb(i),11)+parbe(ix,6)*bbcu(imbb(i),12)
              cikveb_d = parbe(ix,6)*bbcu(imbb(i),11)-parbe(ix,5)*bbcu(imbb(i),12)
            end if

            if(nbeaux(imbb(i)) == 1) then
              ktrack(i) = 41
              if(ibeco == 1) then
                rho2b_d = crkveb_d**2+cikveb_d**2
                tkb_d   = rho2b_d/(two*sigman2(1,imbb(i)))
                beamoff(4,imbb(i)) = ((strack(i)*crkveb_d)/rho2b_d)*(one-exp_mb(-one*tkb_d))
                beamoff(5,imbb(i)) = ((strack(i)*cikveb_d)/rho2b_d)*(one-exp_mb(-one*tkb_d))
              end if
            end if

            if(ktrack(i) == 42) then
              if(ibeco == 1) then
                r2b_d = two*(sigman2(1,imbb(i))-sigman2(2,imbb(i)))
                rb_d  = sqrt(r2b_d)
                rkb_d = (strack(i)*pisqrt)/rb_d
                xrb_d = abs(crkveb_d)/rb_d
                zrb_d = abs(cikveb_d)/rb_d
                if(ibtyp == 0) then
                  call errf(xrb_d,zrb_d,crxb_d,crzb_d)
                  tkb_d = (crkveb_d**2/sigman2(1,imbb(i))+cikveb_d**2/sigman2(2,imbb(i)))*half
                  xbb_d = sigmanq(2,imbb(i))*xrb_d
                  zbb_d = sigmanq(1,imbb(i))*zrb_d
                  call errf(xbb_d,zbb_d,cbxb_d,cbzb_d)
                else if(ibtyp == 1) then
                  tkb_d = (crkveb_d**2/sigman2(1,imbb(i))+cikveb_d**2/sigman2(2,imbb(i)))*half
                  xbb_d = sigmanq(2,imbb(i))*xrb_d
                  zbb_d = sigmanq(1,imbb(i))*zrb_d
                else
                  tkb_d = zero ! -Wmaybe-uninitialized
                end if
              else
                rkb_d = zero ! -Wmaybe-uninitialized
                tkb_d = zero ! -Wmaybe-uninitialized
              end if
              beamoff(4,imbb(i))=(rkb_d*(crzb_d-exp_mb(-one*tkb_d)*cbzb_d))*sign(one,crkveb_d)
              beamoff(5,imbb(i))=(rkb_d*(crxb_d-exp_mb(-one*tkb_d)*cbxb_d))*sign(one,cikveb_d)
            end if

            if(ktrack(i) == 43) then
              if(ibeco == 1) then
                r2b_d = two*(sigman2(2,imbb(i))-sigman2(1,imbb(i)))
                rb_d  = sqrt(r2b_d)
                rkb_d = (strack(i)*pisqrt)/rb_d
                xrb_d = abs(crkveb_d)/rb_d
                zrb_d = abs(cikveb_d)/rb_d
                if(ibtyp == 0) then
                  call errf(zrb_d,xrb_d,crzb_d,crxb_d)
                  tkb_d = (crkveb_d**2/sigman2(1,imbb(i))+cikveb_d**2/sigman2(2,imbb(i)))*half
                  xbb_d = sigmanq(2,imbb(i))*xrb_d
                  zbb_d = sigmanq(1,imbb(i))*zrb_d
                  call errf(zbb_d,xbb_d,cbzb_d,cbxb_d)
                else if(ibtyp == 1) then
                  tkb_d = (crkveb_d**2/sigman2(1,imbb(i))+cikveb_d**2/sigman2(2,imbb(i)))*half
                  xbb_d = sigmanq(2,imbb(i))*xrb_d
                  zbb_d = sigmanq(1,imbb(i))*zrb_d
                else
                  tkb_d = zero ! -Wmaybe-uninitialized
                end if
              else
                rkb_d = zero ! -Wmaybe-uninitialized
                tkb_d = zero ! -Wmaybe-uninitialized
              end if
              beamoff(4,imbb(i)) = (rkb_d*(crzb_d-exp_mb(-one*tkb_d)*cbzb_d))*sign(one,crkveb_d)
              beamoff(5,imbb(i)) = (rkb_d*(crxb_d-exp_mb(-one*tkb_d)*cbxb_d))*sign(one,cikveb_d)
            end if
          end if
        end if
      end do
    end if

  ! Crab Cavities
  ! Note: If setting something else than el(),
  ! DON'T call initialise_element on a crab, it will reset the phase to 0.
  elseif(abs(kz(ix)) == 23) then
    crabph(ix) = el(ix)
    el(ix)     = zero

  ! CC Mult kick order 2
  elseif(abs(kz(ix)) == 26) then
    crabph2(ix) = el(ix)
    el(ix)      = zero

  ! CC Mult kick order 3
  elseif(abs(kz(ix)) == 27) then
    crabph3(ix) = el(ix)
    el(ix)      = zero

  ! CC Mult kick order 4
  else if(abs(kz(ix)) == 28) then
    crabph4(ix) = el(ix)
    el(ix)      = zero

  ! e-lens
  else if(kz(ix) == 29) then
    ed(ix) = zero
    ek(ix) = zero
    el(ix) = zero

  ! Chebyshev lens
  else if(kz(ix) == cheby_kz) then
    ed(ix) = zero
    ek(ix) = zero
    el(ix) = zero
  end if

  return

  ! Error handlers
100 continue
  write(lerr,"(a,i0)") "INITELEM> ERROR Tried to set the strength of an element which is disabled. bez = ", bez(ix)
  call prror

end subroutine initialise_element
