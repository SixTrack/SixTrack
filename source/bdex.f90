! ================================================================================================ !
! Beam Distribution EXchange
! At one or more elements, exchange the current beam distribution
! for one given by an external program.
! K. Sjobak BE/ABP-HSS, 2016
! Based on FLUKA coupling version by
! A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team, 2014.
! ================================================================================================ !
module bdex

  use floatPrecision
  use crcoall
  use numerical_constants
  use parpro

  use mod_alloc
  use string_tools

  implicit none

  ! Is BDEX in use?
  logical, save :: bdex_enable = .false.
  ! Debug mode?
  logical, save :: bdex_debug  = .false.

  ! BDEX in use for this element?
  ! 0: No.
  ! 1: Do a particle exchange at this element.
  ! Other values are reserved for future use.
  integer, allocatable, save :: bdex_elementAction(:) !(nele)
  ! Which BDEX channel does this element use?
  ! This points in the bdex_channel arrays.
  integer, allocatable, save :: bdex_elementChannel(:) !(nele)

  integer, parameter :: bdex_maxchannels = 16
  integer, save      :: bdex_nchannels   = 0

  ! Basic data for the bdex_channels, one row/channel
  ! Column 1: Type of channel. Values:
  !           0: Channel not in use
  !           1: PIPE channel
  !           2: TCPIP channel (not implemented)
  ! Column 2: Meaning varies, based on the value of col. 1:
  !           If col 1 is PIPE, then it is the output format.
  ! Column 3: Meaning varies, based on the value of col. 1:
  !           If col 1 is PIPE, it points to the first (of two) files in bdex_stringStorage.
  ! Column 4: Meaning varies, based on the value of col. 1:
  !           If col 1 is PIPE, then it is the unit number to use (first of two consecutive).
  integer, save :: bdex_channels(bdex_maxchannels,4) = 0
  ! The names of the BDEX channel
  character(len=mStrLen), save :: bdex_channelNames(bdex_maxchannels) = " "

  !Number of places in the bdex_xxStorage arrays
  integer, parameter :: bdex_maxStore       = 20
  integer, save      :: bdex_nstringStorage = 0
  character(len=mStrLen), save :: bdex_stringStorage(bdex_maxStore) = " "

contains

subroutine bdex_allocate_arrays
  use crcoall
  implicit none
  call alloc(bdex_elementAction,nele,0,'bdex_elementAction')
  call alloc(bdex_elementChannel,nele,0,'bdex_elementChannel')
end subroutine bdex_allocate_arrays

subroutine bdex_expand_arrays(nele_new)
  use crcoall
  implicit none
  integer, intent(in) :: nele_new
  call alloc(bdex_elementAction,nele_new,0,'bdex_elementAction')
  call alloc(bdex_elementChannel,nele_new,0,'bdex_elementChannel')
end subroutine bdex_expand_arrays

subroutine bdex_closeFiles
  integer i
  if(bdex_enable) then
      do i=0,bdex_nchannels
        if (bdex_channels(i,1).eq.1) then
            close(bdex_channels(i,4))   !inPipe
            write(bdex_channels(i,4)+1,"(a)") "CLOSEUNITS"
            close(bdex_channels(i,4)+1) !outPipe
        endif
      enddo
  endif
end subroutine bdex_closeFiles

! ================================================================================================ !
!  Parse Input Line
!  K. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-06-25
!  Based on FLUKA coupling version by
!  A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team, 2014.
!  Rewritten from code from DATEN by VKBO
! ================================================================================================ !
subroutine bdex_parseInputLine(inLine, iLine, iErr)

  use string_tools

  implicit none

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lerr,"(a)") "BDEX> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if
  if(nSplit == 0) return

  select case(lnSplit(1)(1:4))

  case("DEBU")
    bdex_debug = .true.
    write (lout,"(a)") "BDEX> DEBUG enabled"

  case("ELEM")
    call bdex_parseElem(inLine,iErr)
    if(iErr) return

  case("CHAN")
    call bdex_parseChan(inLine,iErr)
    if(iErr) return

  case default
    write (lerr,"(a)") "BDEX> ERROR Expected keywords DEBU, NEXT, ELEM, or CHAN"
    iErr = .true.
    return

  end select

end subroutine bdex_parseInputLine

subroutine bdex_parseElem(inLine,iErr)

  use mod_common
  use string_tools

  implicit none

  character(len=*), intent(in)    :: inLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit, ii, jj
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lerr,"(a)") "BDEX> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(bdex_debug)then
    write(lout,"(a,i0,a)") "BDEX> DEBUG Got an ELEM block, len=",len(inLine),": '"//trim(inLine)//"'"
    do ii=1,nSplit
      write(lout,"(a,i3,a)") "BDEX> DEBUG Field(",ii,") = '"//trim(lnSplit(ii))//"'"
    end do
  end if

  ! Parse ELEM
  if(nSplit /= 4) then
    write(lerr,"(a)") "BDEX> ERROR ELEM expects the following arguments:"
    write(lerr,"(a)") "BDEX>       ELEM chanName elemName action"
    iErr = .true.
    return
  end if

  jj = -1
  do ii=1,il ! Match the single element
    if(bez(ii) == lnSplit(3)) then
      jj = ii
      exit
    end if
  end do
  if(jj == -1) then
    write(lerr,"(a)") "BDEX> ERROR The element '"//trim(lnSplit(3))//"' was not found in the single element list."
    iErr = .true.
    return
  end if

  if(kz(jj) /= 0 .or. el(jj) > pieni) then
    write(lerr,"(2(a,i0))") "BDEX> ERROR The element '"//trim(bez(jj))//"' is not a marker. kz=",kz(jj),", el=",el(jj)
    iErr = .true.
    return
  end if

  ! Action
  call chr_cast(lnSplit(4),bdex_elementAction(jj),iErr)
  if(bdex_elementAction(jj) /= 1) then
    write(lerr,"(a,i0)") "BDEX> ERROR Only action 1 (exchange) is currently supported, got",bdex_elementAction(jj)
    iErr = .true.
    return
  end if

  bdex_elementChannel(jj) = -1
  do ii=1,bdex_nchannels ! Match channel name
    if(bdex_channelNames(ii) == lnSplit(2)) then
      bdex_elementChannel(jj) = ii
      exit
    end if
  end do
  if(bdex_elementChannel(jj) == -1) then
    write(lerr,"(a)") "BDEX> ERROR The channel '"//trim(lnSplit(2))//"' was not found."
    iErr = .true.
    return
  end if

end subroutine bdex_parseElem

subroutine bdex_parseChan(inLine,iErr)

  use mod_common
  use string_tools

  implicit none

  character(len=*), intent(in)    :: inLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit, ii, jj
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lerr,"(a)") "BDEX> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(bdex_debug) then
    write (lout,"(a,i0,a)") "BDEX> DEBUG Got a CHAN block, len=",len(inLine), ": '"//trim(inLine)//"'"
    do ii=1,nSplit
      write (lout,"(a,i3,a)") "BDEX> DEBUG Field(",ii,") ='"//trim(lnSplit(ii))//"'"
    end do
  end if

  if(nSplit < 3) then
    write(lerr,"(a,i0)") "BDEX> ERROR CHAN expects at least 3 arguments, got ",nSplit
    iErr = .true.
    return
  end if

  ! Parse CHAN
  select case(trim(lnSplit(3)))
  case("PIPE")
    call bdex_initialisePipe(inLine,iErr)
    if(iErr) return
  case("TCPIP")
    call bdex_initialiseTCPIP(inLine,iErr)
    if(iErr) return
  case default
    write(lerr,"(a)") "BDEX> ERROR Unknown keyword in CHAN: '"//trim(lnSplit(3))//"'. Expected PIPE or TCPIP"
    iErr = .true.
    return
  end select

end subroutine bdex_parseChan

subroutine bdex_parseInputDone

  use mod_common
  use string_tools

  implicit none

  integer ii

  if(bdex_debug) then
    write(lout,"(a)") "BDEX> DEBUG Finished parsing BDEX block"
  end if
  bdex_enable = .true.

  if(bdex_debug) then
    write(lout,"(a)")    "BDEX> DEBUG Done parsing block, data dump:"
    write(lout,"(a,l1)") "BDEX> DEBUG * bdex_enable = ",bdex_enable
    write(lout,"(a,l1)") "BDEX> DEBUG * bdex_debug  = ",bdex_debug
    do ii=1,il
      if(bdex_elementAction(ii) /= 0) then
        write(lout,"(3(a,i0))") "BDEX> DEBUG Single element number",ii,"named '"//trim(bez(ii))//"' "//&
          "bdex_elementAction(#)=",bdex_elementAction(ii)," bdex_elementChannel(#)=",bdex_elementChannel(ii)
      end if
    end do
    write(lout,"(2(a,i0))") "BDEX> DEBUG bdex_nchannels=",bdex_nchannels," >= ",bdex_maxchannels
    do ii=1,bdex_nchannels
      write(lout,"(a,i0,a,4(i0,1x))") "BDEX> DEBUG Channel #",ii," bdex_channelNames(#)='"//trim(bdex_channelNames(ii))//"'"//&
        "bdex_channels(#,:)=",bdex_channels(ii,1),bdex_channels(ii,2),bdex_channels(ii,3),bdex_channels(ii,4)
    end do
    write(lout,"(2(a,i0))") "BDEX> DEBUG bdex_nstringStorage=",bdex_nstringStorage," >= ",bdex_maxStore
    do ii=1,bdex_nstringStorage
      write(lout,"(a,i0,a)") "BDEXDEBUG> #",ii,"= '"//trim(bdex_stringStorage(ii))//"'"
    end do
    write(lout,"(a)") "BDEX> DEBUG Dump completed."
  end if

end subroutine bdex_parseInputDone

! The following subroutines where extracted from deck bdexancil:
! Deck with the initialization etc. routines for BDEX
subroutine bdex_initialisePIPE(inLine,iErr)

  use mod_common
  use parpro
  use string_tools

  implicit none

  character(len=*), intent(in)    :: inLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit
  logical spErr

  logical lopen
  integer stat

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lerr,"(a)") "BDEX> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  ! PIPE: Use a pair of pipes to communicate the particle distributions
  ! Arguments: InFileName OutFileName format fileUnit
  if(nSplit /= 7) then
    write(lerr,"(a)") "BDEX> ERROR CHAN PIPE expects the following arguments:"
    write(lerr,"(a)") "BDEX>       CHAN chanName PIPE InFileName OutFileName format fileUnit"
    iErr = .true.
    return
  end if

  bdex_nchannels = bdex_nchannels+1
  if(bdex_nchannels > bdex_maxchannels) then
    write(lerr,"(a)") "BDEX> ERROR Max channels exceeded!"
    iErr = .true.
    return
  end if

  if(bdex_nStringStorage+2 > bdex_maxStore) then
    write(lerr,"(a)") "BDEX> ERROR maxStore exceeded for strings!"
    iErr = .true.
    return
  end if

  ! Store config data

  bdex_channelNames(bdex_nchannels) = trim(lnSplit(2)) ! channelName
  bdex_channels(bdex_nchannels,1) = 1                  ! TYPE is PIPE
  bdex_nstringStorage = bdex_nstringStorage+1
  bdex_channels(bdex_nchannels,3) = bdex_nstringStorage

  bdex_stringStorage(bdex_nstringStorage) = trim(lnSplit(4)) ! inPipe
  bdex_nstringStorage = bdex_nstringStorage+1
  bdex_stringStorage(bdex_nstringStorage) = trim(lnSplit(5)) ! outPipe

  call chr_cast(lnSplit(6),bdex_channels(bdex_nchannels,2),iErr) ! Output Format
  call chr_cast(lnSplit(7),bdex_channels(bdex_nchannels,4),iErr) ! fileUnit

  ! Open the inPipe
  inquire(unit=bdex_channels(bdex_nchannels,4),opened=lopen)
  if(lopen) then
    write(lerr,"(a,i0,a)")"BDEX> ERROR CHAN:PIPE unit=",bdex_channels(bdex_nchannels,4),&
      " for file '"//bdex_stringStorage(bdex_channels(bdex_nchannels,3))//"' was already taken"
    iErr = .true.
    return
  end if

  write(lout,"(a)") "BDEX> Opening input pipe '"//trim(bdex_stringStorage(bdex_channels(bdex_nchannels,3)))//"'"
  open(unit=bdex_channels(bdex_nchannels,4), file=bdex_stringStorage(bdex_channels(bdex_nchannels,3) ),&
    action='read',iostat=stat,status="OLD")
  if(stat /= 0) then
    write(lerr,"(a,i0)") "BDEX> ERROR opening file '",bdex_stringStorage( bdex_channels(bdex_nchannels,3)),"', stat=",stat
    iErr = .true.
    return
  end if

  ! Open the outPipe
  inquire(unit=bdex_channels(bdex_nchannels,4)+1,opened=lopen)
  if(lopen) then
    write(lerr,"(a,i0,a)")"BDEX> ERROR CHAN:PIPE unit=",bdex_channels(bdex_nchannels,4)+1,&
      " for file '"//bdex_stringStorage(bdex_channels(bdex_nchannels,3)+1 )//"' was already taken"
    iErr = .true.
    return
  end if

  write(lout,"(a)") "BDEX> Opening output pipe '"//trim(bdex_stringStorage(bdex_channels(bdex_nchannels,3)+1))//"'"
  open(unit=bdex_channels(bdex_nchannels,4)+1,file=bdex_stringStorage( bdex_channels(bdex_nchannels,3)+1),&
    action='write',iostat=stat,status="OLD")
  if(stat /= 0) then
    write(lerr,"(a,i0)") "BDEX> ERROR opening file '",bdex_stringStorage( bdex_channels(bdex_nchannels,3)+1 ),"' stat=",stat
    iErr = .true.
    return
  end if
  write(bdex_channels(bdex_nchannels,4)+1,"(a)") "BDEX-PIPE !******************!"

end subroutine bdex_initialisePipe

subroutine bdex_initialiseTCPIP(inLine,iErr)

  use parpro
  use string_tools

  implicit none

  character(len=*), intent(in)    :: inLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lerr,"(a)") "BDEX> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  ! TCPIP: Communicate over a TCP/IP port, like the old FLUKA coupling version did.
  ! Currently not implemented.
  write(lerr,"(a)") "BDEX> ERROR CHAN:TCPIP currently not supported in BDEX."
  iErr = .true.
  return

end subroutine bdex_initialiseTCPIP

! The following subroutines where extracted from deck bdexancil:
! Deck with the routines used by BDEX during tracking
! i  : current structure element
! ix : current single element
! n  : turn number
subroutine bdex_track(i,ix,n)

  use parpro
  use string_tools
  use mod_common
  use mod_common_track
  use mod_common_main

  implicit none

  integer, intent(in) :: i, ix, n

  character(len=25) tmpStr(11)
  logical rErr

  !Temp variables
  integer j, k, ii

  if(bdex_elementAction(ix) == 1) then ! Particle exchange
    if (bdex_debug) then
      write(lout,"(a,i0)") "BDEX> DEBUG Doing particle exchange in bez = ",bez(ix)
    endif

    if(bdex_channels(bdex_elementChannel(ix),1) == 1) then ! PIPE channel
      ! TODO: Fix the format!
      write(bdex_channels(bdex_elementChannel(ix),4)+1,"(a,i10,1x,a,a,1x,a,i10,1x,a,i5)") &
        "BDEX TURN=",n,"BEZ=",bez(ix),"I=",i,"NAPX=",napx

      ! Write out particles
      do j=1,napx
        call chr_fromReal(xv1(j),  tmpStr(1), 19,2,rErr)
        call chr_fromReal(yv1(j),  tmpStr(2), 19,2,rErr)
        call chr_fromReal(xv2(j),  tmpStr(3), 19,2,rErr)
        call chr_fromReal(yv2(j),  tmpStr(4), 19,2,rErr)
        call chr_fromReal(sigmv(j), tmpStr(5), 19,2,rErr)
        call chr_fromReal(ejv(j),   tmpStr(6), 19,2,rErr)
        call chr_fromReal(ejfv(j),  tmpStr(7), 19,2,rErr)
        call chr_fromReal(rvv(j),   tmpStr(8), 19,2,rErr)
        call chr_fromReal(dpsv(j),  tmpStr(9), 19,2,rErr)
        call chr_fromReal(oidpsv(j),tmpStr(10),19,2,rErr)
        call chr_fromReal(dpsv1(j), tmpStr(11),19,2,rErr)
        write(bdex_channels(bdex_elementChannel(ix),4)+1,"(11(a25,1x),i24)") tmpStr(1),tmpStr(2),tmpStr(3),tmpStr(4),tmpStr(5),&
          tmpStr(6),tmpStr(7),tmpStr(8),tmpStr(9),tmpStr(10),tmpStr(11),partID(j)
      end do
      write(bdex_channels(bdex_elementChannel(ix),4)+1,'(a)') "BDEX WAITING..."

      !Read back particles
      read(bdex_channels(bdex_elementChannel(ix),4),*) j

      if ( j == -1 ) then
        !Don't change the distribution at all
        if (bdex_debug) then
          write(lout,"(a)") "BDEX> DEBUG No change in distribution."
        endif
      else
        ! TODO: Handle j < npart
        if (j > npart) then
          call expand_arrays(nele, j, nblz, nblo)
          ! TODO: Initialise per-particle optics variables
        endif
        napx=j
        if (bdex_debug) then
          write(lout,"(a,i0,a)") "BDEX> DEBUG Reading ",napx," particles back."
        endif
        do j=1,napx
          read(bdex_channels(bdex_elementChannel(ix),4),*) &
                xv1(j),yv1(j),xv2(j),yv2(j),sigmv(j), &
                ejv(j),ejfv(j),rvv(j),dpsv(j),oidpsv(j), &
                dpsv1(j),partID(j)
          ! TODO: Handle secondary particle tracking arrays properly
          ! TODO: CRLIBM
        enddo
      endif

      write(bdex_channels(bdex_elementChannel(ix),4)+1,'(a)') "BDEX TRACKING..."
    endif

  else
    write(lerr,"(a,i0,a)") "BDEX> ERROR elementAction = ",bdex_elementAction(i)," not understood."
    call prror
  endif

end subroutine bdex_track

end module bdex
