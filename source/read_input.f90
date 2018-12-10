! ================================================================================================ !
!  Read Input
! ~~~~~~~~~~~~
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-18
!  Module for parsing special input files like fort.33, fort.13, fort.10, etc.
! ================================================================================================ !

module read_input

  use crcoall

  implicit none

contains

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-10-31
!  Rewritten parsing of initial coordinates from fort.13. Moved from maincr.
! ================================================================================================ !
subroutine readFort13

  use parpro
  use mod_hions
  use mod_common
  use mod_commonmn
  use string_tools
  use mod_units

  implicit none

  integer            j, ioStat, nPair
  logical            cErr, rErr, fErr
  character(len=100) inLine

  nPair = 0
  rErr  = .false.
  cErr  = .false.
  fErr  = .false.

  call units_openUnit(unit=13,file="fort.13",formatted=.true.,mode="r",err=fErr)
  if(fErr) then
    write(lout,"(a)") "FORT13> ERROR Could not open fort.13."
    call prror(-1)
  end if

  do j=1,napx,2
    nPair = nPair + 1

    ! Particle 1

    read(13,"(a)", iostat=ioStat) inLine
    rErr = ioStat /= 0
    call chr_cast(trim(inLine), xv1(j), cErr)

    read(13,"(a)", iostat=ioStat) inLine
    rErr = ioStat /= 0
    call chr_cast(trim(inLine), yv1(j), cErr)

    read(13,"(a)", iostat=ioStat) inLine
    rErr = ioStat /= 0
    call chr_cast(trim(inLine), xv2(j), cErr)

    read(13,"(a)", iostat=ioStat) inLine
    rErr = ioStat /= 0
    call chr_cast(trim(inLine), yv2(j), cErr)

    read(13,"(a)", iostat=ioStat) inLine
    rErr = ioStat /= 0
    call chr_cast(trim(inLine), sigmv(j), cErr)

    read(13,"(a)", iostat=ioStat) inLine
    rErr = ioStat /= 0
    call chr_cast(trim(inLine), dpsv(j), cErr)

    ! Particle 2

    read(13,"(a)", iostat=ioStat) inLine
    rErr = ioStat /= 0
    call chr_cast(trim(inLine), xv1(j+1), cErr)

    read(13,"(a)", iostat=ioStat) inLine
    rErr = ioStat /= 0
    call chr_cast(trim(inLine), yv1(j+1), cErr)

    read(13,"(a)", iostat=ioStat) inLine
    rErr = ioStat /= 0
    call chr_cast(trim(inLine), xv2(j+1), cErr)

    read(13,"(a)", iostat=ioStat) inLine
    rErr = ioStat /= 0
    call chr_cast(trim(inLine), yv2(j+1), cErr)

    read(13,"(a)", iostat=ioStat) inLine
    rErr = ioStat /= 0
    call chr_cast(trim(inLine), sigmv(j+1), cErr)

    read(13,"(a)", iostat=ioStat) inLine
    rErr = ioStat /= 0
    call chr_cast(trim(inLine), dpsv(j+1), cErr)

    ! Reference Energy + Energy of Particle 1 & 2

    read(13,"(a)", iostat=ioStat) inLine
    rErr = ioStat /= 0
    call chr_cast(trim(inLine), e0, cErr)

    read(13,"(a)", iostat=ioStat) inLine
    rErr = ioStat /= 0
    call chr_cast(trim(inLine), ejv(j), cErr)

    read(13,"(a)", iostat=ioStat) inLine
    rErr = ioStat /= 0
    call chr_cast(trim(inLine), ejv(j+1), cErr)

    ! Report Errors
    if(rErr) then
      write(lout,"(a,i0,a)") "FORT13> ERROR While reading particle pair ",nPair," from file."
      call prror(-1)
    end if
    if(cErr) then
      write(lout,"(a,i0,a)") "FORT13> ERROR While converting particle pair ",nPair," to float."
      call prror(-1)
    end if

    mtc(j)        = one
    mtc(j+1)      = one
    nucm(j)       = nucm0
    nucm(j+1)     = nucm0

    ejfv(j)       = sqrt(ejv(j)**2   - nucm(j)**2)
    ejfv(j+1)     = sqrt(ejv(j+1)**2 - nucm(j+1)**2)
    oidpsv(j)     = one/(one + dpsv(j))
    oidpsv(j+1)   = one/(one + dpsv(j+1))
    moidpsv(j)    = mtc(j)   / (one+dpsv(j))
    moidpsv(j+1)  = mtc(j+1) / (one+dpsv(j+1))
    omoidpsv(j)   = c1e3*((one-mtc(j))   * oidpsv(j))
    omoidpsv(j+1) = c1e3*((one-mtc(j+1)) * oidpsv(j+1))

  end do

  e0f = sqrt(e0**2 - nucm0**2)

end subroutine readFort13

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-18
!  Reading fort.33. Moved from maincr/mainda.
! ================================================================================================ !
subroutine readFort33

  use string_tools
  use mod_common,     only : clo6, clop6
  use sixtrack_input, only : sixin_echoVal
  use mod_settings
  use mod_units
  use parpro

  implicit none

  character(len=:), allocatable :: lnSplit(:)
  character(len=mInputLn)       :: inLine
  integer nSplit, ioStat
  logical spErr, fErr

  call units_openUnit(unit=33,file="fort.33",formatted=.true.,mode="r",err=fErr)

  read(33,"(a)",iostat=ioStat) inLine
  if(ioStat /= 0) then
    write(lout,"(a,i0)") "READ33> ERROR Failed to read line from 'fort.33'. iostat = ",ioStat
    call prror(-1)
  end if

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "READ33> ERROR Failed to parse line from 'fort.33'"
    call prror(-1)
  end if

  if(nSplit > 0) call chr_cast(lnSplit(1),clo6(1), spErr)
  if(nSplit > 1) call chr_cast(lnSplit(2),clop6(1),spErr)
  if(nSplit > 2) call chr_cast(lnSplit(3),clo6(2), spErr)
  if(nSplit > 3) call chr_cast(lnSplit(4),clop6(2),spErr)
  if(nSplit > 4) call chr_cast(lnSplit(5),clo6(3), spErr)
  if(nSplit > 5) call chr_cast(lnSplit(6),clop6(3),spErr)

  if(st_debug) then
    call sixin_echoVal("clo6(1)", clo6(1), "READ33",1)
    call sixin_echoVal("clop6(1)",clop6(1),"READ33",1)
    call sixin_echoVal("clo6(2)", clo6(2), "READ33",1)
    call sixin_echoVal("clop6(2)",clop6(2),"READ33",1)
    call sixin_echoVal("clo6(3)", clo6(3), "READ33",1)
    call sixin_echoVal("clop6(3)",clop6(3),"READ33",1)
  end if

end subroutine readFort33

end module read_input
