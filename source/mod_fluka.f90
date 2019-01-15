module mod_fluka

  use floatPrecision
  use numerical_constants
  use mod_alloc
  use mod_units, only : f_requestUnit

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
  public :: fluka_mod_expand_arrays
  public :: fluka_mod_end

  public :: fluka_connect
  public :: fluka_end

  public :: fluka_send_receive
  public :: fluka_send
  public :: fluka_receive
  public :: fluka_shuffleLostParticles
  public :: fluka_set_synch_part
  public :: fluka_init_max_uid
  public :: fluka_is_running
  public :: fluka_init_brhono

  public :: fluka_close

  public :: fluka_parsingDone
  public :: fluka_parseInputLine

#ifdef ROOT
  public :: root_FLUKA_DumpInsertions
#endif

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

  integer(kind=int32) :: ntconnect,   &
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
  integer(kind=int32) :: fluka_cid

  ! FLUK input block
  logical, public :: fluka_enable    = .false.                     ! enable coupling
  logical, public :: fluka_connected = .false.                     ! fluka is connected
  logical, public :: fluka_debug     = .false.                     ! write debug messages
  integer, public :: fluka_log_unit                    ! logical unit for log messages (was 888)
  ! hisix: write isotope info
  integer, public :: isotope_log_unit                  ! logical unit for isotope-id output (was 822)

  ! fluka insertions
  logical, public :: fluka_inside = .false.                        ! Are we in a fluka insertion?
  integer(kind=int32), public, allocatable :: fluka_type(:)        ! type of insertion (one per SINGLE ELEMENT)
  integer(kind=int32), public, allocatable :: fluka_geo_index(:)   ! index of insertion (one per SINGLE ELEMENT)
  real(kind=fPrec), public, allocatable :: fluka_synch_length(:)   ! length of insertion [m] (one per SINGLE ELEMENT)
  ! recognised insertion types
  integer(kind=int32), parameter, public :: FLUKA_NONE    = 0, & ! no insertion
                                            FLUKA_ELEMENT = 1, & ! insertion covers only the present SINGLE ELEMENT
                                            FLUKA_ENTRY   = 2, & ! SINGLE ELEMENT marking the start of the insertion
                                            FLUKA_EXIT    = 3    ! SINGLE ELEMENT marking the end   of the insertion
  ! ancillary tracking values
  integer(kind=int32), public :: fluka_max_npart                          ! Maximum number of particles (array size)
  integer(kind=int32), public :: fluka_max_uid                            ! Highest particle ID
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
  integer(kind=int16),          public :: fluka_chrg0  ! []
  integer(kind=int16),          public :: fluka_a0     ! nucelon number (hisix)
  integer(kind=int16),          public :: fluka_z0     ! charge multiplicity (hisix)

  save

  contains

  !----------------------------------------------------------------------------
  ! set the module up
  subroutine fluka_mod_init(npart, nele, clight)

    implicit none

    ! interface variables
    integer :: npart, nele
    real(kind=fPrec) :: clight

    ! temporary variables
    integer :: j

    fluka_max_npart = npart
    fluka_clight    = clight

    call alloc(fluka_uid,          npart, 0, 'fluka_uid')
    call alloc(fluka_gen,          npart, 0, 'fluka_gen')
    call alloc(fluka_weight,       npart, one, 'fluka_weight')
    call alloc(fluka_type,         nele, FLUKA_NONE, 'fluka_type')
    call alloc(fluka_geo_index,    nele, 0, 'fluka_geo_index')
    call alloc(fluka_synch_length, nele, zero, 'fluka_synch_length')

    do j = 1, npart
      fluka_uid(j) = j
      fluka_gen(j) = j
    end do

!    fluka_weight       = one
!    fluka_type         = FLUKA_NONE
!    fluka_geo_index    = 0
!    fluka_synch_length = zero

    call f_requestUnit('fluka.log', fluka_log_unit)
    call f_requestUnit('fluka_isotope.log', isotope_log_unit)
    open(unit=fluka_log_unit, file='fluka.log')
    open(unit=isotope_log_unit, file='fluka_isotope.log')

  end subroutine fluka_mod_init

  subroutine fluka_mod_expand_arrays(npart_new, nele_new)

    use parpro, only : npart

    implicit none

    integer :: npart_new, nele_new, j

    call alloc(fluka_uid,          npart_new, 0, 'fluka_uid')
    call alloc(fluka_gen,          npart_new, 0, 'fluka_gen')
    call alloc(fluka_weight,       npart_new, one, 'fluka_weight')
    call alloc(fluka_type,         nele_new, FLUKA_NONE, 'fluka_type')
    call alloc(fluka_geo_index,    nele_new, 0, 'fluka_geo_index')
    call alloc(fluka_synch_length, nele_new, zero, 'fluka_synch_length')

    do j = npart+1, npart_new
      fluka_uid(j) = j
      fluka_gen(j) = j
    end do

    fluka_max_npart = npart_new

  end subroutine fluka_mod_expand_arrays

  !----------------------------------------------------------------------------
  ! un-set the module
  subroutine fluka_mod_end()
    implicit none
    call dealloc(fluka_uid,'fluka_uid')
    call dealloc(fluka_gen,'fluka_gen')
    call dealloc(fluka_weight,'fluka_weight')
    call dealloc(fluka_type,'fluka_type')
    call dealloc(fluka_geo_index,'fluka_geo_index')
    call dealloc(fluka_synch_length,'fluka_synch_length')

    close(fluka_log_unit)
    close(isotope_log_unit)
  end subroutine fluka_mod_end

  !----------------------------------------------------------------------------
  ! acquire info for network communication
  subroutine fluka_read_config(net_nfo_file, host, port)
    implicit none

    ! interface variables
    character(len=255) :: net_nfo_file
    character(len=255) :: host
    integer :: port
    integer :: net_nfo_unit
    integer :: ios

    call f_requestUnit(net_nfo_file, net_nfo_unit)
    open(net_nfo_unit, file=net_nfo_file, status='old')
    read(unit=net_nfo_unit, fmt=*, iostat=ios) host
    if(ios .ne. 0) then
      write(lout,*)
      write(lout,*) 'FLUKA> Could not read the host name from network.nfo'
      write(lout,*)
      call prror(-1)
    end if

    read(unit=net_nfo_unit, fmt=*, iostat=ios) port
    if(ios .ne. 0) then
      write(lout,*)
      write(lout,*) 'FLUKA> Could not read the port number from network.nfo'
      write(lout,*) 'FLUKA> Is the FLUKA server running and has it had time to write the port number?'
      write(lout,*)
      call prror(-1)
    end if

    close(net_nfo_unit)

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
    integer(kind=int32) :: n

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
  integer function fluka_send_receive(turn, ipt, el, npart, xv1, xv2, yv1, yv2, s, etot, aa, zz, mass)
    implicit none

    ! Parameters
    integer(kind=int32) :: turn, ipt
    integer           ::  npart
    real(kind=fPrec)  :: el

    real(kind=fPrec), allocatable :: xv1(:)
    real(kind=fPrec), allocatable :: yv1(:)
    real(kind=fPrec), allocatable :: xv2(:)
    real(kind=fPrec), allocatable :: yv2(:)
    real(kind=fPrec), allocatable :: s(:)
    real(kind=fPrec), allocatable :: etot(:)

    real(kind=fPrec), allocatable :: mass(:)
    integer(kind=int16), allocatable :: aa(:)
    integer(kind=int16), allocatable :: zz(:)

    fluka_send_receive = fluka_send(turn, ipt, el, npart, xv1, xv2, yv1, yv2, s, etot, aa, zz, mass)
    if(fluka_send_receive.eq.-1) return

    fluka_send_receive = fluka_receive(turn, ipt, el, npart, xv1, xv2, yv1, yv2, s, etot, aa, zz, mass)
  end function fluka_send_receive

  !----------------------------------------------------------------------------
  ! just send particles to Fluka
  integer function fluka_send(turn, ipt, el, npart, xv1, xv2, yv1, yv2, s, etot, aa, zz, mass)
    implicit none

    ! Interface variables
    integer(kind=int32) :: turn, ipt
    integer           :: npart
    real(kind=fPrec)  :: el

    real(kind=fPrec), allocatable :: xv1(:)
    real(kind=fPrec), allocatable :: yv1(:)
    real(kind=fPrec), allocatable :: xv2(:)
    real(kind=fPrec), allocatable :: yv2(:)
    real(kind=fPrec), allocatable :: s(:)
    real(kind=fPrec), allocatable :: etot(:)

    real(kind=fPrec), allocatable :: mass(:)
    integer(kind=int16), allocatable :: aa(:)
    integer(kind=int16), allocatable :: zz(:)

    ! Fluka I/O parameters
    integer(kind=int32) :: flid, flgen
    real(kind=fPrec)    :: flwgt, flx, fly, flz, flxp, flyp, flzp, flet, flm, flt
    integer(kind=int16) :: flaa, flzz
    integer(kind=int8)  :: mtype

    ! Auxiliary variables
    integer :: j
    integer(kind=int32) :: n

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
!    flzz = 1
!   mass number:
!    flaa = 1
!   particle mass [GeV/c2]:
!    flm  = fluka_mass0

    do j=1, npart

      flid  = fluka_uid(j)
      flgen = fluka_gen(j)
      flwgt = fluka_weight(j)

      flx   = xv1(j) * c1m1  ! from [mm] to [cm]
      fly   = xv2(j) * c1m1  ! from [mm] to [cm]
      flz   = zero

      flxp  = yv1(j) * c1m3 ! from [1.0E-03] to [1.0]
      flyp  = yv2(j) * c1m3 ! from [1.0E-03] to [1.0]
      ! director cosines:
      ! full transformation:
      flzp  = sqrt( one / ( flxp**2 + flyp**2 + one ) )
!      ! taylor expansion, for consistency with drifts in SixTrack:
!      flzp  = 1d0 / ( 1d0 + ( flxp**2+flyp**2 )/2d0 )
      flxp  = flxp * flzp
      flyp  = flyp * flzp

      ! total energy:
      flet  = etot(j) * c1m3 ! from [MeV] to [GeV]
      ! longitudinal phase:
      flt   = -s(j) * c1m3 / ( (fluka_pc0/fluka_e0)*fluka_clight ) ! from [mm] to [s]


      ! Ion properties (PH for hiSix)
      flm   = mass(j) * c1m3      ! unit is [GeV]
      flaa  = aa(j)
      flzz  = zz(j)

      if(fluka_debug) then
        write(fluka_log_unit, '(">",2I8,7(1X,1PE25.18),2I8)') flid, flgen, &
             flx, fly, flxp, flyp, flm, flet, flt, flaa, flzz             !PH: added flaa,flzz
        flush(fluka_log_unit)
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
  ! The call from fluka.s90 is:
  ! fluka_receive( nturn, fluka_geo_index(ix), eltot, napx, xv1(:), yv1(:), xv2(:), yv2(:), sigmv, ejv, naa(:), nzz(:), nucm(:))
  ! When the above arrays are made allocatable, the below variables will need updating - see mod_common_main and mod_hions
  integer function fluka_receive(turn, ipt, el, napx, xv1, xv2, yv1, yv2, s, etot, aa, zz, mass)

    use parpro

    implicit none

    ! Interface variables
    integer(kind=int32) :: turn, ipt
    integer           :: napx
    real(kind=fPrec)  :: el

    real(kind=fPrec), allocatable :: xv1(:)
    real(kind=fPrec), allocatable :: yv1(:)
    real(kind=fPrec), allocatable :: xv2(:)
    real(kind=fPrec), allocatable :: yv2(:)
    real(kind=fPrec), allocatable :: s(:)
    real(kind=fPrec), allocatable :: etot(:)

    real(kind=fPrec), allocatable :: mass(:)
    integer(kind=int16), allocatable :: aa(:)
    integer(kind=int16), allocatable :: zz(:)

    ! Fluka I/O parameters
    integer(kind=int32) :: flid, flgen
    real(kind=fPrec)    :: flwgt, flx, fly, flz, flxp, flyp, flzp, flet, flm, flt
    integer(kind=int16) :: flaa, flzz
    integer(kind=int8)  :: mtype

    ! Auxiliary variables
    integer(kind=int32) :: n, j

    fluka_receive = 0

    fluka_nrecv = 0
    mtype = 0

    ! assign default values
    do j = 1, npart
      fluka_uid(j) = j
      fluka_gen(j) = j

      fluka_weight(j) = one

      xv1 (j) = zero
      xv2 (j) = zero
      yv1 (j) = zero
      yv2 (j) = zero
      etot(j) = zero
      s   (j) = zero
! hisix: we should also parse m0,A0,Z0
      aa  (j) = 1
      zz  (j) = 1
      mass(j) = zero
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

         if(fluka_nrecv .gt. npart) then

            !If we hit the particle limit, we will need to  do a global array expand on npart, lets increase by 50 for now
            call expand_arrays(nele, npart+50, nblz, nblo)

!            write(fluka_log_unit, *) &
!                 '# FlukaIO error: reached maximum number of particles, ', &
!                 'no space left to store other incoming particles'
!            fluka_cid = -1
!            fluka_receive = -1
!            return
         end if

            if(fluka_debug) then
               write(fluka_log_unit, '("<",2I8,7(1X,1PE25.18),2I8)') flid, flgen, &
                    flx, fly, flxp, flyp, flm, flet, flt, flaa, flzz ! PH for hiSix
               flush(fluka_log_unit)
            end if

            fluka_uid(fluka_nrecv)    = flid
            fluka_gen(fluka_nrecv)    = flgen
            if (fluka_uid(fluka_nrecv).gt.fluka_max_uid) then
               fluka_max_uid = fluka_uid(fluka_nrecv)
! AM ->                ! generate a new uid
! AM ->                fluka_max_uid = fluka_max_uid + 1
! AM ->                fluka_uid(fluka_nrecv) = fluka_max_uid
!
! PH for hisix: write the particle species and their initial conditions to fort.822
!
               write(isotope_log_unit,*) fluka_uid(fluka_nrecv),flgen, ipt, flaa, flzz, flet * c1e3

            end if

            fluka_weight(fluka_nrecv) = flwgt
            xv1(fluka_nrecv)         = flx * c1e1   ! from [cm]  to [mm]
            xv2(fluka_nrecv)         = fly * c1e1   ! from [cm]  to [mm]
            yv1(fluka_nrecv)         = flxp / flzp * c1e3 ! from director cosine to x' [1.0E-03]
            yv2(fluka_nrecv)         = flyp / flzp * c1e3 ! from director cosine to x' [1.0E-03]
            etot(fluka_nrecv)         = flet * c1e3  ! from [GeV] to [MeV]
            s(fluka_nrecv)            = ( el - (fluka_pc0/fluka_e0)*(flt*fluka_clight) ) * c1e3 ! from [s] to [mm]
            aa(fluka_nrecv)           = flaa          !PH for hiSix
            zz(fluka_nrecv)           = flzz          !PH for hiSix
            mass(fluka_nrecv)         = flm  * c1e3  ! from [GeV] to [MeV]         !PH for hiSix
      end if

      !Finished waiting end of turn
    end do

    napx = fluka_nrecv

    write(fluka_log_unit,*) "# FlukaIO: turn = ", turn, &
      " ipt = ", ipt, &
      " sent = ", fluka_nsent, &
      " received = ", fluka_nrecv, &
      " max_uid = ", fluka_max_uid
    flush(fluka_log_unit)

  end function fluka_receive

  !----------------------------------------------------------------------------
  ! compact ancillary tracking arrays
  subroutine fluka_shuffleLostParticles(tnapx, j)

    integer, intent(in) :: tnapx
    integer, intent(in) :: j

    if(fluka_debug) then
      write(fluka_log_unit, *) '# fluka_shuffleLostParticles called with napx (lnapx for SixTrack) = ', tnapx, ', j = ', j
      flush(fluka_log_unit)
    end if

    fluka_uid(j:tnapx)    = cshift(fluka_uid(j:tnapx),    1)
    fluka_gen(j:tnapx)    = cshift(fluka_gen(j:tnapx),    1)
    fluka_weight(j:tnapx) = cshift(fluka_weight(j:tnapx), 1)

  end subroutine fluka_shuffleLostParticles

  !----------------------------------------------------------------------------
  ! set reference particle properties (mainly for longitudinal dynamics)
  integer function fluka_set_synch_part( e0, pc0, mass0, a0, z0 )
    implicit none

    ! interface variables
    real(kind=fPrec) :: e0, pc0, mass0
    integer(kind=int16)          :: a0, z0

    ! Auxiliary variables
    integer(kind=int32) :: n

    fluka_set_synch_part = 0

    fluka_e0    = e0    *c1m3 ! from  [MeV]    to [GeV]
    fluka_pc0   = pc0   *c1m3 ! from  [MeV/c]  to [GeV/c]
    fluka_mass0 = mass0 *c1m3 ! from  [MeV/c2] to [GeV/c2]
!    fluka_chrg0 = chrg0
    fluka_a0 = a0
    fluka_z0 = z0

    write(fluka_log_unit,*) ' updated synch part:'
    write(fluka_log_unit,*) ' - total en    [GeV]:',fluka_e0
    write(fluka_log_unit,*) ' - momentum  [GeV/c]:',fluka_pc0
    write(fluka_log_unit,*) ' - mass     [GeV/c2]:',fluka_mass0
!    write(fluka_log_unit,*) ' - charge        [e]:',fluka_chrg0
    write(fluka_log_unit,*) ' - mass number    []:',fluka_a0
    write(fluka_log_unit,*) ' - charge        [e]:',fluka_z0
    flush(fluka_log_unit)

    ! update magnetic rigidity, unless division by clight and 10^-9
    fluka_brho0 = fluka_pc0 / real(fluka_z0,real64)

    ! inform Fluka about the new magnetic rigidity
    n = ntsendbrhono(fluka_cid, fluka_brho0)
    if (n .lt. 0) then
      fluka_set_synch_part = -1
      return
    end if
    write(fluka_log_unit,*) ' synchronised magnetic rigidity with Fluka'
    write(fluka_log_unit,*) '    transmitted value [Tm/0.3]:', fluka_brho0
    write(fluka_log_unit,*) '    in proper units       [Tm]:', fluka_brho0 / ( fluka_clight*c1m9 )
    flush(fluka_log_unit)

  end function fluka_set_synch_part

  !----------------------------------------------------------------------------
  ! set max ID
  integer function fluka_init_max_uid( npart )
    implicit none

    ! interface variables
    integer(kind=int32) :: npart

    ! Auxiliary variables
    integer(kind=int32) :: n

    fluka_init_max_uid = 0

    fluka_max_uid = npart

    n = ntsendnpart(fluka_cid, npart)
    if (n .lt. 0) then
      fluka_init_max_uid = -1
      return
    end if

  end function fluka_init_max_uid


  !----------------------------------------------------------------------------
  ! set Brho nominal
  integer function fluka_init_brhono( brhono )
    implicit none

    ! interface variables
    real(kind=real64) :: brhono

    ! Auxiliary variables
    integer :: n

    n = ntsendbrhono(fluka_cid, brhono)
    if (n .lt. 0) then
      fluka_init_brhono = -1
      return
    end if

  end function fluka_init_brhono


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

!     A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!     last modified: 17-07-2013
!     clean closure of communication with fluka and un-set mod_fluka
!     inserted in main code by the 'fluka' compilation flag
subroutine fluka_close

  use crcoall

  implicit none

  integer fluka_con

  fluka_con = fluka_is_running()
  if(fluka_con.eq.0) then
    if( .not. fluka_connected ) then
!         temporarily connect to fluka, to properly terminate the run
      fluka_con = fluka_connect()
      if(fluka_con.eq.-1) then
!           no hope to properly close the run
        write(lout,*) '[Fluka] unable to connect to fluka while'
        write(lout,*) '        closing the simulation: please,'
        write(lout,*) '        manually kill all its instances'
        write(fluka_log_unit,*) '# unable to connect to fluka while'
        write(fluka_log_unit,*) '#  closing the simulation: please,'
        write(fluka_log_unit,*) '#  manually kill all its instances'
        goto 1982
      endif
      write(lout,*) '[Fluka] Successfully connected to Fluka server'
      write(lout,*) '[Fluka]     (only temporarily)'
      write(fluka_log_unit,*) '# Successfully connected to Fluka server'
      write(fluka_log_unit,*) '#     (only temporarily)'
    end if
    call fluka_end
  end if
1982 call fluka_mod_end
  flush(lout)
!      flush(fluka_log_unit)
end subroutine fluka_close

! ================================================================================================ !
!  Parse Fluka Coupling Input Line
!  A. Mereghetti, D. Sinuela Pastor, for the FLUKA Team
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-06-25
! ================================================================================================ !
subroutine fluka_parseInputLine(inLine, iLine, iErr)

  use string_tools
  use mod_common, only : il,bez

  implicit none

  character(len=*), intent(in)    :: inLine
  integer,          intent(in)    :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mNameLen) entrElem, exitElem
  real(kind=fPrec) tmplen
  integer nSplit, i, entrIdx, exitIdx, ii
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "FLUKA> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  select case(lnSplit(1)(1:4))

  case("DEBU")
    write(lout,"(a)") "FLUKA> DEBUG enabled"
    fluka_debug = .true.

  case("LOGU")
    write(lout,"(a,i0,a)") "FLUKA> NOTE The LOGU flag is deprecated. A log unit is assigned automatically."

  case default

    if(nSplit /= 4) then
      write(lout,"(a,i0)") "FLUKA> ERROR Exected 4 values in input line,got ",nSplit
      iErr = .true.
      return
    end if

    entrElem = " "
    exitElem = " "

    if(nSplit > 0) entrElem = trim(lnSplit(1))
    if(nSplit > 1) exitElem = trim(lnSplit(2))
    if(nSplit > 2) call chr_cast(lnSplit(3),ii,    iErr)
    if(nSplit > 3) call chr_cast(lnSplit(4),tmplen,iErr)

    ! 1. find name of entrance element in the list of SINGLE ELEMENTs:
    entrIdx = -1
    do i=1,il
      if(bez(i) == entrElem) then
        entrIdx = i
        exit
      end if
    end do
    if(entrIdx == -1) then
      write(lout,"(a)") "FLUKA> ERROR Un-identified SINGLE ELEMENT '"//trim(entrElem)//"'"
      iErr = .true.
      return
    end if

    ! 2. find exit element in the list of SINGLE ELEMENTs:
    exitIdx = -1
    do i=1,il
      if(bez(i) == exitElem) then
        exitIdx = i
        exit
      end if
    end do
    if(exitIdx == -1) then
      write(lout,"(a)") "FLUKA> ERROR Un-identified SINGLE ELEMENT '"//trim(exitElem)//"'"
      iErr = .true.
      return
    end if

    ! 3. check that the current markers have not been already flagged
    if(fluka_type(entrIdx) /= FLUKA_NONE ) then
      write(lout,"(a)")       "FLUKA> ERROR Single element '"//trim(bez(entrIdx))//"' was alredy labelled as fluka marker."
      write(lout,"(2(a,i0))") "FLUKA> ERROR fluka_type(entrance) = ",fluka_type(entrIdx)," at position = ",entrIdx
      iErr = .true.
      return
    end if
    if(fluka_type(exitIdx) /= FLUKA_NONE ) then
      write(lout,"(a)")       "FLUKA> ERROR Single element '"//trim(bez(exitIdx))//"' was alredy labelled as fluka marker."
      write(lout,"(2(a,i0))") "FLUKA> ERROR fluka_type(exit) = ",fluka_type(exitIdx)," at position = ",exitIdx
      iErr = .true.
      return
    end if

    ! 4. disentangle between just a simple element or an interval of elements
    !    in the accelerator structure, labelled as Fluka insertion:
    if(entrIdx == exitIdx) then
      fluka_type(entrIdx)         = FLUKA_ELEMENT
      fluka_geo_index(entrIdx)    = ii
      fluka_synch_length(entrIdx) = tmplen
      write(fluka_log_unit,"(a,i0)") "# Found         Fluka element as SING EL num ",entrIdx
    else
      fluka_type(entrIdx)         = FLUKA_ENTRY
      fluka_geo_index(entrIdx)    = ii
      fluka_type(exitIdx)         = FLUKA_EXIT
      fluka_geo_index(exitIdx)    = ii
      fluka_synch_length(exitIdx) = tmplen
      write(fluka_log_unit,"(a,i0)") "# Found entrance Fluka element as SING EL num ",entrIdx
      write(fluka_log_unit,"(a,i0)") "# Found exit     Fluka element as SING EL num ",exitIdx
    end if

    ! Wait to find at least one FLUKA insertion before actually enabling the coupling
    if(.not.fluka_enable) fluka_enable = .true.

  end select

end subroutine fluka_parseInputLine

subroutine fluka_parsingDone

  use mod_common, only : il,bez

  integer ii

  if(fluka_enable) then
    ! Dump all elements found:
    write(lout,"(a)") "FLUKA>  Name                 | Type | Insertion Point | Synch Length [m] "
    write(lout,"(a)") "FLUKA> ----------------------+------+-----------------+------------------"
    do ii=1,il
      if(fluka_type(ii) /= FLUKA_NONE) then
        write(lout,"(a,a20,a,i4,a,i15,a,e15.8)") "FLUKA>  ",bez(ii)(1:20)," | ",fluka_type(ii)," | ",&
          fluka_geo_index(ii)," | ",fluka_synch_length(ii)
      end if
    end do
    write(lout,"(a)")    "FLUKA> ----------------------+------+-----------------+------------------"
    write(lout,"(a)")    "FLUKA> Keys to FLUKA types:"
    write(lout,"(a,i0)") "FLUKA> Simple element: ",FLUKA_ELEMENT
    write(lout,"(a,i0)") "FLUKA> Entrance point: ",FLUKA_ENTRY
    write(lout,"(a,i0)") "FLUKA> Exit point:     ",FLUKA_EXIT
  else
    write(lout,"(a)") "FLUKA> WARNING No elements flagged for coupling!"
    write(lout,"(a)") "FLUKA>         Disabling coupling flags/labelling."
    fluka_enable = .false.
    fluka_debug  = .false.
    do ii=1,il
      fluka_type(ii) = FLUKA_NONE
    end do
  end if

end subroutine fluka_parsingDone

#ifdef ROOT
subroutine root_FLUKA_DumpInsertions

  use root_output
  use mod_common, only : bez
  use parpro, only : mNameLen, nele

  implicit none

! loop index
  integer :: k

! fluka id
  integer(kind=int32) :: ii

! name to go with the fluka id
  character(len=mNameLen+1) :: this_name

! loop over each element entry
  do k=0, nele
!   extract the fluka geo index value, which usually will be 0 for non-insertions
    ii = fluka_geo_index(k)
    if(ii .eq. 0) then
      continue
    else

      if(fluka_type(k) .eq. FLUKA_ENTRY) then
!       this entry exists, so add it to root
        this_name = trim(adjustl(bez(k))) // C_NULL_CHAR
        call root_FLUKA_Names(ii, this_name, len_trim(this_name))
      end if

    end if
  end do

end subroutine root_FLUKA_DumpInsertions
#endif

end module mod_fluka

