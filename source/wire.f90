module wire

  use floatPrecision
  use parpro, only : nele, nblz
  use numerical_constants, only : zero

  implicit none

  ! A. Patapenka (NIU), M. Fitterer (FNAL)
  ! Common block for wire definition
  ! variables to save wire parameters for tracking etc.

  ! wire current [A]
  real(kind=fPrec), allocatable, save :: wire_current(:) !(nele)
  ! integrated length of the wire [m]
  real(kind=fPrec), allocatable, save :: wire_lint(:) !(nele)
  ! physical length of the wire [m]
  real(kind=fPrec), allocatable, save :: wire_lphys(:) !(nele)

  ! integer to include or not closed orbit in the separation between beam and wire
  ! 0  : Un-initialized if wire element not found
  ! +1 : dispx is the distance between x0=y0=0 and the wire
  ! -1 : dispx is the distance between the closed orbit and the wire
  !    x=y=0    <->   xco     <->    xwire
  !               closed orbit    wire position
  ! wire_flagco = +1: dispx = xwire -> rx = x + xsep
  ! wire_flagco = -1: dispx = xwire - xco -> rx = x - xco + xsep
  ! -> rx = x + xwire
  integer, allocatable, save          :: wire_flagco(:) !(nele)
  ! hor./vert. displacement of the wire [mm]
  real(kind=fPrec), allocatable, save :: wire_dispx(:),wire_dispy(:) !(nele)
  ! hor./vert. tilt of the wire [degrees] -90 < tilty < 90,
  !  uses the same definition as the DISP block
  real(kind=fPrec), allocatable, save :: wire_tiltx(:), wire_tilty(:) !(nele)


  ! wire parameters for closed orbit calculation (FOX part)
  ! for FOX length of variable names must be smaller 8
  integer, parameter :: wire_max = 350 ! max. number of wires (same as BB interactions)
  real(kind=fPrec), save :: wire_clo(6,wire_max) = zero ! closed orbit at wire

  ! wire number for each structure element (default = 0 if no wire)
  integer, allocatable, save :: wire_num(:) ! (nblz)

contains

subroutine wire_expand_arrays(nele_new, nblz_new)

  use mod_alloc
  use numerical_constants, only : zero

  implicit none

  integer, intent(in) :: nele_new
  integer, intent(in) :: nblz_new

  call alloc(wire_current,nele_new,zero,'wire_current')
  call alloc(wire_lint,nele_new,zero,'wire_lint')
  call alloc(wire_lphys,nele_new,zero,'wire_lphys')

  call alloc(wire_flagco,nele_new,0,'wire_flagco')
  call alloc(wire_num,nblz_new,0,"wire_num")

  call alloc(wire_dispx,nele_new,zero,'wire_dispx')
  call alloc(wire_dispy,nele_new,zero,'wire_dispy')
  call alloc(wire_tiltx,nele_new,zero,'wire_tiltx')
  call alloc(wire_tilty,nele_new,zero,'wire_tilty')

end subroutine wire_expand_arrays

! ================================================================================================ !
!  Parse WIRE Input Line kz=+/-15,ktrack=45
!  A. Patapenka, NIU
!  M. Fitterer, FNAL
!  V. K. Berglyd Olsen, CERN, BE-ABP-HSS
!  Last modified: 2018-06-25
! ================================================================================================ !
subroutine wire_parseInputLine(inLine, iLine, iErr)

  use mod_common
  use mod_settings
  use sixtrack_input
  use string_tools

  implicit none

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit, j, iElem
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "WIRE> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(nSplit /= 9) then
    write(lout,"(a,i0)") "WIRE> ERROR Expected 9 input va;ues, got ",nSplit
    iErr = .true.
    return
  end if

  ! Loop over single elements and set parameters of wire
  iElem = -1
  do j=1,nele
    if(bez(j) == trim(lnSplit(1))) then
      iElem = j
      exit
    end if
  end do
  if(iElem == -1) then
    write(lout,"(a)") "WIRE> ERROR Element '"//trim(lnSplit(1))//"' not found in single element list."
    iErr = .true.
    return
  end if

  if(kz(iElem) /= 15) then
    write(lout,"(2(a,i0),a)") "WIRE> ERROR Element type kz(",iElem,") = ",kz(iElem)," != +15"
    iErr = .true.
    return
  end if
  if(el(iElem) /= 0 .or. ek(iElem) /= 0 .or. ed(iElem) /= 0) then
    ! Check the element type (kz(iElem)_wire=+/-15)
    write(lout,"(a)")       "WIRE> ERROR Length el(iElem) (wire is treated as thin element), "//&
      "and first and second field have to be zero:"
    write(lout,"(2(a,i0))") "WIRE>       el(",iElem,") = ",el(iElem)," != 0"
    write(lout,"(2(a,i0))") "WIRE>       ed(",iElem,") = ",ed(iElem)," != 0"
    write(lout,"(2(a,i0))") "WIRE>       ek(",iElem,") = ",ek(iElem)," != 0"
    iErr = .true.
    return
  end if
  if(wire_flagco(iElem) /= 0) then
    write(lout,"(a)") "WIRE> ERROR The element '"//trim(bez(iElem))//"' was defined twice."
    iErr = .true.
    return
  end if

  ! Parse the element
  call chr_cast(lnSplit(2),wire_flagco(iElem), iErr)
  call chr_cast(lnSplit(3),wire_current(iElem),iErr)
  call chr_cast(lnSplit(4),wire_lint(iElem),   iErr)
  call chr_cast(lnSplit(5),wire_lphys(iElem),  iErr)
  call chr_cast(lnSplit(6),wire_dispx(iElem),  iErr)
  call chr_cast(lnSplit(7),wire_dispy(iElem),  iErr)
  call chr_cast(lnSplit(8),wire_tiltx(iElem),  iErr)
  call chr_cast(lnSplit(9),wire_tilty(iElem),  iErr)

  ! Make checks for the wire parameters
  if(wire_flagco(iElem) /= 1 .and. wire_flagco(iElem) /= -1) then
    write(lout,"(a)")    "WIRE> ERROR Flag for defining the wire separation must be -1 (disp* = distance closed orbit and beam)"
    write(lout,"(a,i0)") "WIRE>       or 1 (disp* = distance from x=y=0 <-> beam), but wire_flagco = ",wire_flagco(iElem)
    iErr = .true.
    return
  end if
  if((wire_lint(iElem) < 0) .or. (wire_lphys(iElem) < 0)) then
    write(lout,"(a)")          "WIRE> ERROR Integrated and physical length must larger than 0."
    write(lout,"(2(a,e15.8))") "WIRE>       wire_lint = ",wire_lint(iElem),", wire_lphys = ",wire_lphys(iElem)
    iErr = .true.
    return
  end if
  if((abs(wire_tiltx(iElem)) >= 90) .or. (abs(wire_tilty(iElem)) >= 90)) then
    write(lout,"(a)")          "WIRE> ERROR Tilt angle must be within [-90,90] degrees."
    write(lout,"(2(a,e15.8))") "WIRE>       wire_tiltx = ",wire_tiltx(iElem),", wire_tilty = ",wire_tilty(iElem)
    iErr = .true.
    return
  end if

  if(st_debug) then
    call sixin_echoVal("name",    bez(iElem),         "WIRE",iLine)
    call sixin_echoVal("flagco",  wire_flagco(iElem), "WIRE",iLine)
    call sixin_echoVal("current", wire_current(iElem),"WIRE",iLine)
    call sixin_echoVal("int_len", wire_lint(iElem),   "WIRE",iLine)
    call sixin_echoVal("phys_len",wire_lphys(iElem),  "WIRE",iLine)
    call sixin_echoVal("disp_x",  wire_dispx(iElem),  "WIRE",iLine)
    call sixin_echoVal("disp_y",  wire_dispy(iElem),  "WIRE",iLine)
    call sixin_echoVal("tilt_x",  wire_tiltx(iElem),  "WIRE",iLine)
    call sixin_echoVal("tilt_y",  wire_tilty(iElem),  "WIRE",iLine)
  end if
  if(iErr) return

  if(abs(wire_flagco(iElem)*(wire_current(iElem)*(wire_lint(iElem)*&
        (wire_lphys(iElem)*(wire_dispx(iElem)+wire_dispy(iElem)))))) <= pieni) then
    kz(iElem) = 0 ! treat element as marker
  end if

end subroutine wire_parseInputLine

subroutine wire_parseInputDone(iErr)

  use crcoall
  use mod_common, only : kz,bez

  implicit none

  logical, intent(inout) :: iErr

  integer j

  ! Loop over single elements to check that they have been defined in the fort.3 block
  do j=1,nele
    if(kz(j) == 15) then
      if(wire_flagco(j) == 0) then
        write(lout,"(a)") "WIRE> ERROR Wire element '"//trim(bez(j))//"'not defined in fort.3."
        write(lout,"(a)") "WIRE>       You must define every wire in the WIRE block."
        iErr = .true.
        return
      end if
    end if
  end do

end subroutine wire_parseInputDone

end module wire
