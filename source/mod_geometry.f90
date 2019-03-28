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
  use mod_common
  use string_tools
  use sixtrack_input

  implicit none

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
    write(lout,"(a)") "GEOMETRY> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(nSplit <= 2) then
    write(lout,"(a,i0)") "GEOMETRY> ERROR Single element line must have more than 2 values, got ",nSplit
    iErr = .true.
    return
  end if

  elemName = lnSplit(1)
  if(len(elemName) > mNameLen) then
    write(lout,"(a,i0)") "GEOMETRY> ERROR Single element name too long. Max length is ",mNameLen
    iErr = .true.
    return
  end if

  ! Check that the name is unique
  do i=1,geom_nSing-1
    if(bez(i) == elemName) then
      write(lout,"(a,i0)") "GEOMETRY> ERROR Single element '"//trim(elemName)//"' is not unique."
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
  endif

  ! CHANGING SIGN OF CURVATURE OF VERTICAL THICK DIPOLE
  if((kz(geom_nSing) == 4 .or. kz(geom_nSing) == 5) .and. abs(el(geom_nSing)) > pieni) then
    ed(geom_nSing) = -one*ed(geom_nSing)
  end if

  ! CAVITIES
  if(abs(kz(geom_nSing)) == 12) then
    if(abs(ed(geom_nSing)) > pieni .and. abs(ek(geom_nSing)) > pieni) then
      sixin_ncy2          = sixin_ncy2 + 1
      itionc(geom_nSing) = kz(geom_nSing)/abs(kz(geom_nSing))
      kp(geom_nSing)     = 6
    end if
  end if

  !--------------------------------------------
  ! Handled by initialize_element subroutine:
  !--------------------------------------------
  ! CHANGING SIGN OF CURVATURE OF VERTICAL THICK DIPOLE
  ! THIN LENS (+/- 1-10)
  ! MULTIPOLES (11)
  ! CAVITY (+/- 12)
  ! CRABCAVITY (23/-23) / CC multipoles order 2/3/4 (+/- 23/26/27/28)
  ! ELECTRON LENSE (29)
  call initialize_element(geom_nSing,.true.)

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
    call expand_arrays(nele+100, npart, nblz, nblo)
    call alloc(geom_bez0, mNameLen, nele, " ", "geom_bez0")
  end if

  if(abs(kz(geom_nSing)) /= 12 .or. (abs(kz(geom_nSing)) == 12 .and. sixin_ncy2 == 0)) then
    kp(geom_nSing) = 0
  end if

  bez(geom_nSing)        = elemName
  geom_bez0(geom_nSing) = elemName

  !If no active RF cavities are seen so far in the single element list,
  ! add a CAV element to the end of the list.
  ! This is then overwritten when reading the next element, so that if
  ! and only if no active RF cavities are found, a CAV element can be
  ! used in the structure to enable 6D tracking using the parameters
  ! from the SYNC block.
  if(sixin_ncy2 == 0) then
    geom_nSing             = geom_nSing + 1
    il                     = geom_nSing
    bez(geom_nSing)        = geom_cavity
    geom_bez0(geom_nSing)  = geom_cavity
    kp(geom_nSing)         = 6
  else
    il          = geom_nSing
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

  implicit none

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
    write(lout,"(a)") "GEOMETRY> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if
  isCont = (nInd >= 5)

  if(nSplit < 2 .and. .not. isCont) then
    write(lout,"(a,i0)") "GEOMETRY> ERROR Block definition line must be at least 2 values, got ",nSplit
    iErr = .true.
    return
  end if

  if(nSplit > 40) then
    write(lout,"(a)") "GEOMETRY> ERROR Block definition cannot have more then 40 elements."
    iErr = .true.
    return
  end if

  ! If first line, read super period information
  if(iLine == 1) then
    call chr_cast(lnSplit(1),mper,iErr)
    if(mper > nper) then
      write(lout,"(a,i0)") "GEOMETRY> ERROR Block definition number of super periods is too large. Max value is ",nper
      iErr = .true.
      return
    end if
    if(mper > nSplit+1) then
      write(lout,"(a)") "GEOMETRY> ERROR Block definition number of super periods does not match the number of values."
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
  do i=1,40
    ilm0(i) = " "
  end do

  if(isCont) then                             ! This line continues the previous BLOC
    blocName = " "                            ! No name returned, set an empty BLOC name
    do i=1,nSplit                             ! All elements are sub-elements. Save to buffer.
      ilm0(i) = lnSplit(i)
    end do
  else                                        ! This is a new BLOC
    blocName = lnSplit(1)
    do i=1,nSplit-1                           ! Save the rest to buffer
      ilm0(i) = lnSplit(i+1)
    end do
  end if

  if(blocName /= " ") then                    ! We have a new BLOC
    geom_nBloc = geom_nBloc + 1               ! Increment the BLOC number
    if(geom_nBloc > nblo-1) then              ! Expand arrays if needed
      call expand_arrays(nele, npart, nblz, nblo+50)
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
      write(lout,"(a,2(i0,a))") "GEOMETRY> ERROR Block definitions can only have ",&
        nelb," elements. ",i," given."
      iErr = .true.
      return
    end if
    geom_ilm(i) = ilm0(i-k0)                  ! Append to sub-element buffer
    if(geom_ilm(i) == " ") exit               ! No more sub-elements to append
    mel(geom_nBloc)          = i              ! Update number of single elements in this block
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
      write(lout,"(a)") "GEOMETRY> ERROR Unknown single element '"//&
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

  implicit none

  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=:), allocatable   :: expLine
  integer nSplit
  logical spErr

  integer i, j
  character(len=mNameLen) ilm0(40)

  ilm0(:) = " "

  expLine = chr_expandBrackets(inLine)
  call chr_split(expLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "GEOMETRY> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(nSplit < 1 .or. nSplit > 40) then
    write(lout,"(a)") "GEOMETRY> ERROR Structure input line cannot have less than 1 or more then 40 elements."
    iErr = .true.
    return
  end if

  do i=1,nSplit
    ilm0(i) = lnSplit(i)
  end do

  do i=1,40

    if(ilm0(i) == " ") cycle
    if(ilm0(i) == geom_go) then
      kanf = geom_nStru + 1
      cycle
    end if

    geom_nStru = geom_nStru + 1
    if(geom_nStru > nblz-3) then
      call expand_arrays(nele,npart,nblz+1000,nblo)
    end if

    do j=1,mblo ! is it a BLOC?
      if(bezb(j) == ilm0(i)) then
        ic(geom_nStru)     = j
        icname(geom_nStru) = bezb(j)
        cycle
      end if
    end do

    do j=1,il ! is it a SINGLE ELEMENT?
      if(geom_bez0(j) == ilm0(i)) then
        ic(geom_nStru)     = j+nblo
        icname(geom_nStru) = geom_bez0(j)
        if(geom_bez0(j) == geom_cavity) then
          sixin_icy = sixin_icy+1
        end if
        cycle
      end if
    end do
  end do

  mbloz = geom_nStru
  if(mbloz > nblz-3) then
    call expand_arrays(nele,npart,nblz+1000,nblo)
  end if

end subroutine geom_parseInputLineSTRU

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Parse Structure Input Line - Multi-Column Version
!  Created: 2019-03-27
!  Updated: 2019-03-28
! ================================================================================================ !
subroutine geom_parseInputLineSTRU_MULT(inLine, iLine, iErr)

  use crcoall
  use mod_common
  use string_tools
  use sixtrack_input

  implicit none

  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit, singID
  logical spErr, cErr

  integer j

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "GEOMETRY> ERROR Failed to parse input line."
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
    write(lout,"(a)") "GEOMETRY> ERROR Multi-column structure input line cannot have less than 3 elements."
    iErr = .true.
    return
  end if

  geom_nStru = geom_nStru + 1
  if(geom_nStru > nblz-3) then
    call expand_arrays(nele,npart,nblz+1000,nblo)
  end if

  icname(geom_nStru) = lnSplit(1)
  call chr_cast(lnSplit(3), icpos(geom_nStru), cErr)

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
          sixin_icy = sixin_icy+1
        end if
        exit
      end if
    end do
  end if

  if(singID == -1) then
    write(lout,"(a)") "GEOMETRY> ERROR Unknown element '"//trim(lnSplit(2))//"' in STRUCTURE INPUT."
    iErr = .true.
    return
  else
    ic(geom_nStru) = singID
  end if

  mbloz = geom_nStru
  if(mbloz > nblz-3) then
    call expand_arrays(nele,npart,nblz+1000,nblo)
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

  implicit none

  il = il + 1
  if(il > nele-2) then
    call expand_arrays(nele+50, npart, nblz, nblo )
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
  use mod_common

  implicit none

  integer, intent(in) :: iEl

  integer iIns

  if(iu > nblz-3) then
    call expand_arrays(nele, npart, nblz+100, nblo)
  end if

  iu = iu + 1
  if(iEl == 0) then
    iIns = iu
  elseif(iEl < 0) then
    iIns = iu+iEl
  else
    iIns = iEl
  end if

  ic(iIns:iu)      = cshift(ic(iIns:iu),      1)
  icpos(iIns:iu)   = cshift(icpos(iIns:iu),   1)
  icname(iIns:iu)  = cshift(icname(iIns:iu),  1)
  icext(iIns:iu)   = cshift(icext(iIns:iu),   1)
  icextal(iIns:iu) = cshift(icextal(iIns:iu), 1)
  dcum(iIns:iu)    = cshift(dcum(iIns:iu),    1)

  ! Update s coordinate of added element
  icpos(iIns) = icpos(iIns-1)
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

  implicit none

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

  implicit none

! interface variables
  real(kind=fPrec), intent(out) :: sLoc
  logical,          intent(in)  :: isLast
  integer,          intent(out) :: iEl
  integer,          intent(out) :: ixEl
  logical,          intent(out) :: wasFound

  integer i, iDelta, iCheck, iMax, iStep
  logical lSlide

  iEl  = -1
  ixEl = -1

  if(sLoc > tlen .or. sLoc < zero) then
    write(lout,"(a,2(f11.4,a))") "GEOMETRY> ERROR Find Element: "//&
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

end module mod_geometry
