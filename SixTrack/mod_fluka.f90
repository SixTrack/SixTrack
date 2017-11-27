module mod_fluka

  use floatPrecision

  use, intrinsic :: ISO_FORTRAN_ENV, only : int8, int16, int32, int64

  ! A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
  ! last modified: 18-01-2016
  ! fortran 90 module for coupling SixTrack to FLUKA
  ! NOTA BENE: 
  !    napx  (SixTrack) -> npart     (mod_fluka)
  !    npart (SixTrack) -> max_npart (mod_fluka)

  implicit none
  private

  public :: fluka_mod_init
  public :: fluka_mod_end

  public :: fluka_connect
  public :: fluka_end

  public :: fluka_send_receive
  public :: fluka_send
  public :: fluka_receive
  public :: fluka_lostpart
  public :: fluka_set_synch_part
  public :: fluka_init_max_uid
  public :: fluka_is_running


  ! FlukaIO Connection parameters
  character(len = 255), public  :: fluka_host
  integer, public :: fluka_port
  character(len = 255), parameter :: fluka_net_nfo_file = 'network.nfo'
  
  ! FlukaIO interface
  external ntinit, ntconnect, ntend
  external ntsendp,     &
           ntsendeob,   &
           ntsendeoc,   &
           ntsendipt,   &
           ntrecv,      &
           ntwait,      &
           ntsendnpart, &
           ntsendbrhono

  integer :: ntconnect,   &
             ntsendp,     &
             ntsendeob,   &
             ntsendeoc,   &
             ntsendipt,   &
             ntrecv,      &
             ntwait,      &
             ntsendnpart, &
             ntsendbrhono,&
             ntend

  ! FlukaIO Message types
  integer(kind=int8), parameter :: FLUKA_PART = 1, &
                                   FLUKA_EOB  = 2, &
                                   FLUKA_EOC  = 3, &
                                   FLUKA_CONF = 4, &
                                   FLUKA_IPT  = 5, &
                                   FLUKA_HSK  = 6, &
                                   FLUKA_NPART= 7, &
                                   FLUKA_BRHO = 8
  ! connection ID
  integer :: fluka_cid


  ! FLUK input block
  logical, public :: fluka_enable    = .false.                     ! enable coupling
  logical, public :: fluka_connected = .false.                     ! fluka is connected
  logical, public :: fluka_debug     = .false.                     ! write debug messages
  integer, public :: fluka_log_unit  = 888                         ! logical unit for log messages
  ! fluka insertions
  logical, public :: fluka_inside = .false.                        ! Are we in a fluka insertion?
  integer, public, allocatable :: fluka_type(:)                    ! type of insertion (one per SINGLE ELEMENT)
  integer, public, allocatable :: fluka_geo_index(:)               ! index of insertion (one per SINGLE ELEMENT)
  real(kind=fPrec), public, allocatable :: fluka_synch_length(:)   ! length of insertion [m] (one per SINGLE ELEMENT)
  ! recognised insertion types
  integer, parameter, public :: FLUKA_NONE    = 0, & ! no insertion
                                FLUKA_ELEMENT = 1, & ! insertion covers only the present SINGLE ELEMENT
                                FLUKA_ENTRY   = 2, & ! SINGLE ELEMENT marking the start of the insertion
                                FLUKA_EXIT    = 3    ! SINGLE ELEMENT marking the end   of the insertion
  ! ancillary tracking values
  integer, public :: fluka_max_npart                          ! Maximum number of particles (array size)
  integer, public :: fluka_max_uid                            ! Highest particle ID
  integer(kind=int32), public, allocatable :: fluka_uid(:)    ! particle ID
  integer(kind=int32), public, allocatable :: fluka_gen(:)    ! ID of parent particle
  real(kind=fPrec), public, allocatable    :: fluka_weight(:) ! statistical weight (>0.0)

  ! Useful values
  integer :: fluka_nsent     ! Temporary count of sent particles
  integer :: fluka_nrecv     ! Temporary count of received particles
  real(kind=fPrec), public :: fluka_clight ! [m/s]

  ! Reference particle
  real(kind=fPrec), public :: fluka_e0     ! [GeV]
  real(kind=fPrec), public :: fluka_pc0    ! [GeV/c]
  real(kind=fPrec), public :: fluka_mass0  ! [GeV/c2]
  real(kind=fPrec), public :: fluka_brho0  ! [Tm]
  integer,          public :: fluka_chrg0  ! []

  ! hack for Li-7
  !real(kind=fPrec), parameter :: MLI = 6.5338351903884577D0
  !real(kind=fPrec), parameter :: RLI = 3D0/1D0

  contains

  !----------------------------------------------------------------------------
  ! set the module up
  subroutine fluka_mod_init(maxp, nele, clight)
    implicit none

    ! interface variables
    integer :: maxp, nele
    real(kind=fPrec) :: clight

    ! temporary variables
    integer :: j

    fluka_max_npart = maxp
    fluka_clight    = clight

    allocate(fluka_uid(fluka_max_npart))
    allocate(fluka_gen(fluka_max_npart))
    allocate(fluka_weight(fluka_max_npart))
    allocate(fluka_type(nele))
    allocate(fluka_geo_index(nele))
    allocate(fluka_synch_length(nele))

    do j = 1, fluka_max_npart
      fluka_uid(j) = j
      fluka_gen(j) = j
    end do

    fluka_weight       = 1.0D+00
    fluka_type         = FLUKA_NONE
    fluka_geo_index    = 0
    fluka_synch_length = 0.0D+00

    open(fluka_log_unit)

  end subroutine fluka_mod_init

  !----------------------------------------------------------------------------
  ! un-set the module
  subroutine fluka_mod_end()
    implicit none
    deallocate(fluka_uid)
    deallocate(fluka_gen)
    deallocate(fluka_weight)
    deallocate(fluka_type)
    deallocate(fluka_geo_index)
    deallocate(fluka_synch_length)
    close(fluka_log_unit)
  end subroutine fluka_mod_end

  !----------------------------------------------------------------------------
  ! acquire info for network communication
  subroutine fluka_read_config(net_nfo_file, host, port)
    implicit none

    ! interface variables
    character(len=255) :: net_nfo_file
    character(len=255) :: host
    integer :: port

    open(90, file=net_nfo_file, status='old')
    read(90, *) host
    read(90, *) port
    close(90)

  end subroutine fluka_read_config

  !----------------------------------------------------------------------------
  ! start communication with fluka
  integer function fluka_connect()
    implicit none

    call fluka_read_config(fluka_net_nfo_file, fluka_host, fluka_port)

    write(fluka_log_unit,*) '# Connecting to host: ', fluka_host, ', in port: ', fluka_port
    write(fluka_log_unit,*) '# Maximum number of particles: ', fluka_max_npart
    call ntinit()
    fluka_cid = ntconnect(fluka_host, fluka_port)
    fluka_connect = fluka_cid

  end function fluka_connect

  !----------------------------------------------------------------------------
  ! close communication with fluka
  subroutine fluka_end()
    implicit none

    ! Finish connection
    integer :: n

    ! Fluka I/O parameters
    integer(kind=int32)         :: flid, flgen
    real(kind=fPrec)  :: flwgt, flx, fly, flz, flxp, flyp, flpc, flm, flt
    integer(kind=int16)         :: flaa, flzz
    integer(kind=int8)          :: mtype

    write(fluka_log_unit,*) "# FlukaIO: sending End of Computation signal"

    ! Send end of computation
    n = ntsendeoc(fluka_cid)
    if(n.lt.0) then
      write(fluka_log_unit,*) "# FlukaIO error: Error sending End of Computation"
      flush(fluka_log_unit)
      return
    end if

    ! Wait end of comp
    n = ntwait(fluka_cid, mtype, &
          flid, flgen, flwgt, flx, fly, flz, flxp, flyp, flaa, flzz, &
          flm, flpc, flt)
    if(n.eq.-1) then
      write(fluka_log_unit,*) "# FlukaIO error: Server timed out while waiting End of Computation"
      flush(fluka_log_unit)
      return
    end if
    if(mtype.ne.FLUKA_EOC) then
      write(fluka_log_unit,*) "# FlukaIO warning: Received unexpected message at shutdown"
    end if

    ! At this point both ends agreed to disconnect

    ! Close connection
    n = ntend(fluka_cid)

    return
  end subroutine fluka_end

  !----------------------------------------------------------------------------
  ! send and receive particles from Fluka
  integer function fluka_send_receive(turn, ipt, el, npart, x, xp, y, yp, s, etot)
    implicit none

    ! Parameters
    integer           :: turn, ipt, npart
    real(kind=fPrec)  :: el
    real(kind=fPrec)  :: x(fluka_max_npart), xp(fluka_max_npart)
    real(kind=fPrec)  :: y(fluka_max_npart), yp(fluka_max_npart)
    real(kind=fPrec)  :: s(fluka_max_npart), etot(fluka_max_npart)

    fluka_send_receive = fluka_send(turn, ipt, el, npart, x, xp, y, yp, s, etot)
    if(fluka_send_receive.eq.-1) return

    fluka_send_receive = fluka_receive(turn, ipt, el, npart, x, xp, y, yp, s, etot)
  end function fluka_send_receive

  !----------------------------------------------------------------------------
  ! just send particles to Fluka
  integer function fluka_send(turn, ipt, el, npart, x, xp, y, yp, s, etot)
    implicit none

    ! Interface variables
    integer           :: turn, ipt, npart
    real(kind=fPrec)  :: el
    real(kind=fPrec)  :: x(fluka_max_npart), xp(fluka_max_npart)
    real(kind=fPrec)  :: y(fluka_max_npart), yp(fluka_max_npart)
    real(kind=fPrec)  :: s(fluka_max_npart), etot(fluka_max_npart)

    ! Fluka I/O parameters
    integer(kind=int32)         :: flid, flgen
    real(kind=fPrec)  :: flwgt, flx, fly, flz, flxp, flyp, flzp, flet, flm, flt
    integer(kind=int16)         :: flaa, flzz
    integer(kind=int8)          :: mtype

    ! Auxiliary variables
    integer :: j
    integer :: n

    flush(fluka_log_unit)

    fluka_send = 0

    n = ntsendipt(fluka_cid, turn, ipt)
    if(n.eq.-1) then
      write(fluka_log_unit,*) "# FlukaIO error: Error sending Insertion Point"
      fluka_cid = -1
      fluka_send = -1
      return
    end if

    fluka_nsent = 0
    fluka_nrecv = 0
    mtype = 0

!   atomic number:
    flzz = 1
!   mass number:
    flaa = 1
!   particle mass [GeV/c2]:
    flm  = fluka_mass0

    do j=1, npart

      flid  = fluka_uid(j)
      flgen = fluka_gen(j)
      flwgt = fluka_weight(j)

      flx   = x(j) * 1.0D-01  ! from [mm] to [cm]
      fly   = y(j) * 1.0D-01  ! from [mm] to [cm]
      flz   = 0.0D+00

      flxp  = xp(j) * 1.0D-03 ! from [1.0E-03] to [1.0]
      flyp  = yp(j) * 1.0D-03 ! from [1.0E-03] to [1.0]
      ! director cosines:
      ! full transformation:
      flzp  = sqrt( 1.0D+00 / ( flxp**2 + flyp**2 + 1.0D+00 ) )
!      ! taylor expansion, for consistency with drifts in SixTrack:
!      flzp  = 1d0 / ( 1d0 + ( flxp**2+flyp**2 )/2d0 )
      flxp  = flxp * flzp
      flyp  = flyp * flzp

      ! total energy:
      flet  = etot(j) * 1.0D-03 ! from [MeV] to [GeV]
      ! longitudinal phase:
      flt   = -s(j) * 1.0D-03 / ( (fluka_pc0/fluka_e0)*fluka_clight ) ! from [mm] to [s]

      if(fluka_debug) then
        write(fluka_log_unit, '(">",2I8,7(1X,1PE25.18))') flid, flgen, &
             flx, fly, flxp, flyp, flm, flet, flt
      end if

! Hack for lithium-7
!      flm     = MLI
!      flpc    = flpc * RLI
!      flaa    = 7
!      flzz    = 3

      ! Send particle
      n = ntsendp(fluka_cid, &
            flid, flgen, flwgt, &
            flx, fly, flz, &
            flxp, flyp, flzp, &
            flaa, flzz, flm, flet, flt)

      if(n.eq.-1) then
        write(fluka_log_unit,*) "# FlukaIO error: Error sending Particle"
        fluka_cid = -1
        fluka_send = -1
        return
      end if

      fluka_nsent = fluka_nsent + 1

    end do

    ! Send end of batch
    n = ntsendeob(fluka_cid)

    if(n.lt.0) then
      write(fluka_log_unit,*) "# FlukaIO error: Error sending End of Batch"
      fluka_cid = -1
      fluka_send = -1
      return
    end if

  end function fluka_send

  !----------------------------------------------------------------------------
  ! just receive particles from Fluka
  integer function fluka_receive(turn, ipt, el, npart, x, xp, y, yp, s, etot)
    implicit none

    ! Interface variables
    integer           :: turn, ipt, npart
    real(kind=fPrec)  :: el
    real(kind=fPrec)  :: x(fluka_max_npart), xp(fluka_max_npart)
    real(kind=fPrec)  :: y(fluka_max_npart), yp(fluka_max_npart)
    real(kind=fPrec)  :: s(fluka_max_npart), etot(fluka_max_npart)

    ! Fluka I/O parameters
    integer(kind=int32)         :: flid, flgen
    real(kind=fPrec)  :: flwgt, flx, fly, flz, flxp, flyp, flzp, flet, flm, flt
    integer(kind=int16)         :: flaa, flzz
    integer(kind=int8)          :: mtype

    ! Auxiliary variables
    integer :: n, j

    fluka_receive = 0

    fluka_nrecv = 0
    mtype = 0

    ! assign default values
    do j = 1, fluka_max_npart
      fluka_uid(j) = j
      fluka_gen(j) = j

      fluka_weight(j) = 1.0D+00

      x   (j) = 0.0D+00
      y   (j) = 0.0D+00
      xp  (j) = 0.0D+00
      yp  (j) = 0.0D+00
      etot(j) = 0.0D+00
      s   (j) = 0.0D+00
    end do

    ! Wait until end of turn (Synchronize)
    do while(mtype.ne.FLUKA_EOB)
      n = ntwait(fluka_cid, mtype, &
              flid, flgen, flwgt, &
              flx, fly, flz, &
              flxp, flyp, flzp, &
              flaa, flzz, flm, flet, flt)

      if(n.eq.-1) then
        write(fluka_log_unit,*) "# FlukaIO error: Server timed out while waiting for message"
        fluka_cid = -1
        fluka_receive = -1
        return
      end if

      if(mtype.eq.FLUKA_PART) then

         fluka_nrecv = fluka_nrecv + 1

         if(fluka_nrecv.gt.fluka_max_npart) then

            write(fluka_log_unit, *) &
                 '# FlukaIO error: reached maximum number of particles, ', &
                 'no space left to store other incoming particles'
            fluka_cid = -1
            fluka_receive = -1
            return

         else

            if(fluka_debug) then
               write(fluka_log_unit, '("<",2I8,7(1X,1PE25.18))') flid, flgen, &
                    flx, fly, flxp, flyp, flm, flet, flt
            end if

            fluka_uid(fluka_nrecv)    = flid
            fluka_gen(fluka_nrecv)    = flgen
            if (fluka_uid(fluka_nrecv).gt.fluka_max_uid) then
               fluka_max_uid = fluka_uid(fluka_nrecv)
! AM ->                ! generate a new uid
! AM ->                fluka_max_uid = fluka_max_uid + 1
! AM ->                fluka_uid(fluka_nrecv) = fluka_max_uid
            end if
            fluka_weight(fluka_nrecv) = flwgt
            x(fluka_nrecv)            = flx * 1.0D+01   ! from [cm]  to [mm]
            y(fluka_nrecv)            = fly * 1.0D+01   ! from [cm]  to [mm]
            xp(fluka_nrecv)           = flxp / flzp * 1.0D+03 ! from director cosine to x' [1.0E-03] 
            yp(fluka_nrecv)           = flyp / flzp * 1.0D+03 ! from director cosine to x' [1.0E-03]
            etot(fluka_nrecv)         = flet * 1.0D+03  ! from [GeV] to [MeV]
            s(fluka_nrecv)            = ( el - (fluka_pc0/fluka_e0)*(flt*fluka_clight) ) * 1.0D+03 ! from [s] to [mm]

         end if
      end if

      !Finished waiting end of turn
    end do

    npart = fluka_nrecv

    write(fluka_log_unit,*) "# FlukaIO: turn = ", turn, &
      " ipt = ", ipt, &
      " sent = ", fluka_nsent, &
      " received = ", fluka_nrecv, &
      " max_uid = ", fluka_max_uid
    flush(fluka_log_unit)

  end function fluka_receive

  !----------------------------------------------------------------------------
  ! compact ancillary tracking arrays
  subroutine fluka_lostpart(npart, i)
    integer, intent(in) :: npart, i

    if(fluka_debug) then
      write(fluka_log_unit, *) '# fluka_lostpart called with npart (lnapx for SixTrack) = ', npart, ', i = ', i
      flush(fluka_log_unit)
    end if

    fluka_uid(i:npart-1) = fluka_uid(i+1:npart)
    fluka_gen(i:npart-1) = fluka_gen(i+1:npart)
    fluka_weight(i:npart-1) = fluka_weight(i+1:npart)

  end subroutine fluka_lostpart

  !----------------------------------------------------------------------------
  ! set reference particle properties (mainly for longitudinal dynamics)
  integer function fluka_set_synch_part( e0, pc0, mass0, chrg0 )
    implicit none

    ! interface variables
    real(kind=fPrec) :: e0, pc0, mass0
    integer          :: chrg0

    ! Auxiliary variables
    integer :: n

    fluka_set_synch_part = 0

    fluka_e0    = e0    *1.0D-03 ! from  [MeV]    to [GeV]   
    fluka_pc0   = pc0   *1.0D-03 ! from  [MeV/c]  to [GeV/c] 
    fluka_mass0 = mass0 *1.0D-03 ! from  [MeV/c2] to [GeV/c2]
    fluka_chrg0 = chrg0
    write(fluka_log_unit,*) ' updated synch part:'
    write(fluka_log_unit,*) ' - total en    [GeV]:',fluka_e0
    write(fluka_log_unit,*) ' - momentum  [GeV/c]:',fluka_pc0
    write(fluka_log_unit,*) ' - mass     [GeV/c2]:',fluka_mass0
    write(fluka_log_unit,*) ' - charge        [e]:',fluka_chrg0
    flush(fluka_log_unit)

    ! update magnetic rigidity, unless division by clight and 10^-9
    fluka_brho0 = fluka_pc0 / dble(fluka_chrg0)

    ! inform Fluka about the new magnetic rigidity
    n = ntsendbrhono(fluka_cid, fluka_brho0)
    if (n .lt. 0) then
      fluka_set_synch_part = -1
      return
    end if
    write(fluka_log_unit,*) ' synchronised magnetic rigidity with Fluka'
    write(fluka_log_unit,*) '    transmitted value [Tm/0.3]:', fluka_brho0
    write(fluka_log_unit,*) '    in proper units       [Tm]:', fluka_brho0 / ( fluka_clight*1.0D-09 )
    flush(fluka_log_unit)

  end function fluka_set_synch_part

  !----------------------------------------------------------------------------
  ! set max ID
  integer function fluka_init_max_uid( npart )
    implicit none

    ! interface variables
    integer :: npart

    ! Auxiliary variables
    integer :: n

    fluka_init_max_uid = 0

    fluka_max_uid = npart

    n = ntsendnpart(fluka_cid, npart)
    if (n .lt. 0) then
      fluka_init_max_uid = -1
      return
    end if

  end function fluka_init_max_uid

  !----------------------------------------------------------------------------
  ! check if fluka is running, ie if it created the 
  integer function fluka_is_running()
    implicit none

    ! temporary variables
    logical :: lexist

    fluka_is_running = 0
    inquire( file=fluka_net_nfo_file, exist=lexist)

    if (.not.lexist) then
       write(fluka_log_unit,*) '# Error: file containing network infos ', fluka_net_nfo_file
       write(fluka_log_unit,*) '#        does not exist!!'
       fluka_is_running = -1
    endif

  end function fluka_is_running

end module mod_fluka
