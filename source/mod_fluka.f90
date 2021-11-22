module mod_fluka

  use floatPrecision
  use numerical_constants
  use mod_alloc

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

  public :: fluka_close

  public :: fluka_parsingDone
  public :: fluka_parseInputLine

#ifdef ROOT
  public :: root_FLUKA_DumpInsertions
#endif

  public :: kernel_fluka_element
  public :: kernel_fluka_entrance
  public :: kernel_fluka_exit

  public :: check_coupling_integrity
  public :: check_coupling_start_point

  ! HION variables that are only used for FLUKA
  ! ien0,ien1: ion energy entering/leaving the collimator
  real(kind=fPrec),    public :: ien0, ien1
  integer(kind=int16), public :: nnuc0,nnuc1

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
  ! current fluka insertion:
  integer :: fluka_i=-1              ! lattice entry (FLUKA_ELEMENT/FLUKA_EXIT)
  integer :: fluka_ix=-1             ! single element entry (FLUKA_ELEMENT/FLUKA_EXIT)
  integer :: fluka_nturn=-1          ! turn
  integer :: fluka_last_sent_mess=-1 ! last sent message
  integer :: fluka_last_rcvd_mess=-1 ! last received message
  ! recognised insertion types
  integer(kind=int32), parameter, public :: FLUKA_NONE    = 0, & ! no insertion
                                            FLUKA_ELEMENT = 1, & ! insertion covers only the present SINGLE ELEMENT
                                            FLUKA_ENTRY   = 2, & ! SINGLE ELEMENT marking the start of the insertion
                                            FLUKA_EXIT    = 3    ! SINGLE ELEMENT marking the end   of the insertion
  ! ancillary tracking values
  integer(kind=int32), public :: fluka_max_npart                          ! Maximum number of particles (array size)
  integer,          public, allocatable    :: pids(:)         ! Particle ID moved from hisixtrack, to be harmonised

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
  integer(kind=int16),          public :: fluka_z0     ! nuclear charge

  save

contains

  !----------------------------------------------------------------------------
  ! set the module up
  subroutine fluka_mod_init(npart, nele, clight)

    use mod_units
    use mod_common, only : fort208, unit208
    implicit none

    ! interface variables
    integer :: npart, nele
    real(kind=fPrec) :: clight

    ! temporary variables
    integer :: j

    fluka_max_npart = npart
    fluka_clight    = clight

    call alloc(pids,               npart, 0, "pids")
    call alloc(fluka_type,         nele, FLUKA_NONE, 'fluka_type')
    call alloc(fluka_geo_index,    nele, 0, 'fluka_geo_index')
    call alloc(fluka_synch_length, nele, zero, 'fluka_synch_length')

    if(unit208 == -1) then
      call f_requestUnit(fort208,unit208)
      call f_open(unit=unit208,file=fort208,formatted=.true.,mode="w")
#ifdef CR
      fort208Pos = 0
#endif
    end if

    call f_requestUnit("fluka.log",        fluka_log_unit)
    call f_requestUnit("fluka_isotope.log",isotope_log_unit)
    call f_open(unit=fluka_log_unit,  file="fluka.log",        formatted=.true.,mode="w")
    call f_open(unit=isotope_log_unit,file="fluka_isotope.log",formatted=.true.,mode="w")

  end subroutine fluka_mod_init

  subroutine fluka_mod_expand_arrays(npart_new, nele_new)

    use parpro, only : npart

    implicit none

    integer :: npart_new, nele_new, j

    call alloc(pids,               npart_new, 0, "pids")
    call alloc(fluka_type,         nele_new, FLUKA_NONE, 'fluka_type')
    call alloc(fluka_geo_index,    nele_new, 0, 'fluka_geo_index')
    call alloc(fluka_synch_length, nele_new, zero, 'fluka_synch_length')

    fluka_max_npart = npart_new

  end subroutine fluka_mod_expand_arrays

  !----------------------------------------------------------------------------
  ! un-set the module
  subroutine fluka_mod_end()
    implicit none
    call dealloc(pids,"pids")
    call dealloc(fluka_type,'fluka_type')
    call dealloc(fluka_geo_index,'fluka_geo_index')
    call dealloc(fluka_synch_length,'fluka_synch_length')

    close(fluka_log_unit)
    close(isotope_log_unit)
  end subroutine fluka_mod_end

  !----------------------------------------------------------------------------
  ! acquire info for network communication
  subroutine fluka_read_config(net_nfo_file, host, port)

    use mod_units

    implicit none

    ! interface variables
    character(len=255) :: net_nfo_file
    character(len=255) :: host
    integer :: port
    integer :: net_nfo_unit
    integer :: ios

    call f_requestUnit(net_nfo_file, net_nfo_unit)
    call f_open(net_nfo_unit, file=net_nfo_file, formatted=.true., mode="rw", status='old')
    read(unit=net_nfo_unit, fmt=*, iostat=ios) host
    if(ios .ne. 0) then
      write(lerr,'(A)') 'FLUKA> ERROR Could not read the host name from network.nfo'
      call prror
    end if

    read(unit=net_nfo_unit, fmt=*, iostat=ios) port
    if(ios .ne. 0) then
      write(lerr,'(A)') 'FLUKA> ERROR Could not read the port number from network.nfo'
      write(lerr,'(A)') 'FLUKA>       Is the FLUKA server running and has it had time to write the port number?'
      call prror
    end if

    call f_close(net_nfo_unit)

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
    integer(kind=int16)         :: flaa, flzz, flq
    integer(kind=int8)          :: mtype

    integer(kind=int32)         :: flpdgid
    real(kind=fPrec)            :: flsx, flsy, flsz

    write(lout,'(A)') 'FLUKA> call to fluka_end'
    write(fluka_log_unit,'(A)') "# FlukaIO: sending End of Computation signal"

    ! Send end of computation
    n = ntsendeoc(fluka_cid)
    if(n.lt.0) then
      write(fluka_log_unit,'(A,i0,A)') "# FlukaIO error: ", n, " - Error sending End of Computation"
      flush(fluka_log_unit)
      return
    end if

    ! Wait end of comp
    n = ntwait(fluka_cid, mtype, &
          flid, flgen, flwgt, flx, fly, flz, flxp, flyp, &
          flm, flpc, flt, flpdgid, flq, flsx, flsy, flsz)
    if(n.lt.0) then
      write(fluka_log_unit,'(A,i0,A)') "# FlukaIO error: ", n, " - Server timed out while waiting End of Computation"
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
  integer function fluka_send_receive(turn, ipt, el, npart, xv1, xv2, yv1, yv2, s, etot, aa, zz, mass, qq, pdg_id, &
                   partID, parentID, partWeight, spinx, spiny, spinz)
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
    integer(kind=int16), allocatable :: qq(:)
    integer(kind=int32), allocatable :: pdg_id(:)
    integer(kind=int32), allocatable :: partID(:)
    integer(kind=int32), allocatable :: parentID(:)
    real(kind=fPrec), allocatable :: partWeight(:)
    real(kind=fPrec), allocatable :: spinx(:)
    real(kind=fPrec), allocatable :: spiny(:)
    real(kind=fPrec), allocatable :: spinz(:)

    fluka_send_receive = fluka_send(turn, ipt, el, npart, xv1, xv2, yv1, yv2, s, etot, aa, zz, mass, qq, pdg_id, &
                         partID, parentID, partWeight, spinx, spiny, spinz)
    if(fluka_send_receive.lt.0) return

    fluka_send_receive = fluka_receive(turn, ipt, el, npart, xv1, xv2, yv1, yv2, s, etot, aa, zz, mass, qq, pdg_id, &
                         partID, parentID, partWeight, spinx, spiny, spinz)
  end function fluka_send_receive

  !----------------------------------------------------------------------------
  ! just send particles to Fluka
  integer function fluka_send(turn, ipt, el, npart, xv1, xv2, yv1, yv2, s, etot, aa, zz, mass, qq, pdg_id, &
                   partID, parentID, partWeight, spinx, spiny, spinz)
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
    integer(kind=int16), allocatable :: qq(:)
    integer(kind=int32), allocatable :: pdg_id(:)
    integer(kind=int32), allocatable :: partID(:)
    integer(kind=int32), allocatable :: parentID(:)
    real(kind=fPrec), allocatable :: partWeight(:)
    real(kind=fPrec), allocatable :: spinx(:)
    real(kind=fPrec), allocatable :: spiny(:)
    real(kind=fPrec), allocatable :: spinz(:)

    ! Fluka I/O parameters
    integer(kind=int32) :: flid, flgen
    real(kind=fPrec)    :: flwgt, flx, fly, flz, flxp, flyp, flzp, flet, flm, flt
    integer(kind=int16) :: flaa, flzz, flq
    integer(kind=int8)  :: mtype

    integer(kind=int32) :: flpdgid
    real(kind=fPrec)    :: flsx, flsy, flsz

    ! Auxiliary variables
    integer :: j
    integer(kind=int32) :: n

    flush(fluka_log_unit)

    fluka_send = 0
    fluka_last_rcvd_mess = -1

    n = ntsendipt(fluka_cid, turn, ipt)
    if(n.lt.0) then
      write(fluka_log_unit,'(A,i0,A)') "# FlukaIO error: ", n, " - Error sending Insertion Point"
      fluka_cid = -1
      fluka_send = n
      return
    end if
    fluka_last_sent_mess=FLUKA_IPT

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

      flid  = partID(j)
      flgen = parentID(j)
      flwgt = partWeight(j)

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

      flpdgid = pdg_id(j)
      flq  = qq(j)
      flsx = spinx(j)
      flsy = spiny(j)
      flsz = spinz(j)

      if(fluka_debug) then
        write(fluka_log_unit, '(">",2I8,12(1X,1PE25.18),4I11)') flid, flgen, &
             flx, fly, flz, flxp, flyp, flzp, flm, flet, flt, flsx, flsy, flsz, flaa, flzz, flq, flpdgid
        flush(fluka_log_unit)
      end if

      ! Send particle
      n = ntsendp(fluka_cid, &
            flid, flgen, flwgt, &
            flx, fly, flz, &
            flxp, flyp, flzp, &
            flm, flet, flt, &
            flpdgid, flq, flsx, flsy, flsz)

      if(n.lt.0) then
        write(fluka_log_unit,'(A,i0,A)') "# FlukaIO error: ", n, " - Error sending Particle"
        fluka_cid = -1
        fluka_send = -1
        return
      end if

      fluka_nsent = fluka_nsent + 1
      fluka_last_sent_mess=FLUKA_PART

    end do

    ! Send end of batch
    n = ntsendeob(fluka_cid)

    if(n.lt.0) then
      write(fluka_log_unit,'(A,i0,A)') "# FlukaIO error: ", n, " - Error sending End of Batch"
      fluka_cid = -1
      fluka_send = -1
      return
    end if
    fluka_last_sent_mess=FLUKA_EOB

  end function fluka_send

  !----------------------------------------------------------------------------
  ! just receive particles from Fluka
  ! The call from fluka.s90 is:
  ! fluka_receive( nturn, fluka_geo_index(ix), eltot, napx, xv1(:), yv1(:), xv2(:), yv2(:), sigmv, ejv, naa(:), nzz(:), nucm(:))
  ! When the above arrays are made allocatable, the below variables will need updating - see mod_commonmn and mod_hions
  integer function fluka_receive(turn, ipt, el, napx, xv1, xv2, yv1, yv2, s, etot, aa, zz, mass, qq, pdg_id, &
                                 partID, parentID, partWeight, spinx, spiny, spinz)

    use parpro
    use mod_pdgid
    use mod_common_main, only : MaximumPartID

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
    integer(kind=int16), allocatable :: qq(:)
    integer(kind=int32), allocatable :: pdg_id(:)
    integer(kind=int32), allocatable :: partID(:)
    integer(kind=int32), allocatable :: parentID(:)
    real(kind=fPrec), allocatable :: partWeight(:)
    real(kind=fPrec), allocatable :: spinx(:)
    real(kind=fPrec), allocatable :: spiny(:)
    real(kind=fPrec), allocatable :: spinz(:)

    ! Fluka I/O parameters
    integer(kind=int32) :: flid, flgen
    real(kind=fPrec)    :: flwgt, flx, fly, flz, flxp, flyp, flzp, flet, flm, flt
    integer(kind=int16) :: flaa, flzz, flq
    integer(kind=int8)  :: mtype

    integer(kind=int32)         :: flpdgid
    real(kind=fPrec)            :: flsx, flsy, flsz

    ! Auxiliary variables
    integer(kind=int32) :: n, j

    fluka_receive = 0
    fluka_last_sent_mess = -1

    fluka_nrecv = 0
    mtype = 0

    ! assign default values
    do j = 1, npart
      partID(j) = j
      parentID(j) = j

      partWeight(j) = one

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
      qq  (j) = 1
      pdg_id(j) = 0
      spinx = zero
      spiny = zero
      spinz = zero
    end do

    ! Wait until end of turn (Synchronize)
    do while(mtype.ne.FLUKA_EOB)
      n = ntwait(fluka_cid, mtype, &
              flid, flgen, flwgt, &
              flx, fly, flz, &
              flxp, flyp, flzp, &
              flm, flet, flt, &
              flpdgid, flq, flsx, flsy, flsz)

      if(n.lt.0) then
        write(fluka_log_unit,'(A,i0,A)') "# FlukaIO error: ", n ," - Server timed out while waiting for message"
        fluka_cid = -1
        fluka_receive = n
        return
      end if

      if(mtype.eq.FLUKA_PART) then

         fluka_nrecv = fluka_nrecv + 1
         fluka_last_rcvd_mess = FLUKA_PART

         if(fluka_nrecv .gt. npart) then

            !If we hit the particle limit, we will need to  do a global array expand on npart, lets increase by 50 for now
            call expand_arrays(nele, npart+50, nblz, nblo, nbb)

!            write(fluka_log_unit, *) &
!                 '# FlukaIO error: reached maximum number of particles, ', &
!                 'no space left to store other incoming particles'
!            fluka_cid = -1
!            fluka_receive = -1
!            return
         end if

            call CalculateAZ(flpdgid, flaa, flzz)

            if(fluka_debug) then
               write(fluka_log_unit, '("<",2I8,12(1X,1PE25.18),4I11)') flid, flgen, &
                    flx, fly, flz, flxp, flyp, flzp, flm, flet, flt, flsx, flsy, flsz, flaa, flzz, flq, flpdgid
               flush(fluka_log_unit)
            end if

            partID(fluka_nrecv)    = flid
            parentID(fluka_nrecv)    = flgen
            if (partID(fluka_nrecv).gt.MaximumPartID) then
               MaximumPartID = partID(fluka_nrecv)

! PH for hisix: write the particle species and their initial conditions to fort.822
               write(isotope_log_unit,*) partID(fluka_nrecv),flgen, ipt, flaa, flzz, flet * c1e3

            end if

            partWeight(fluka_nrecv)  = flwgt
            xv1(fluka_nrecv)         = flx * c1e1   ! from [cm]  to [mm]
            xv2(fluka_nrecv)         = fly * c1e1   ! from [cm]  to [mm]
            yv1(fluka_nrecv)         = flxp / flzp * c1e3 ! from director cosine to x' [1.0E-03]
            yv2(fluka_nrecv)         = flyp / flzp * c1e3 ! from director cosine to x' [1.0E-03]
            etot(fluka_nrecv)         = flet * c1e3  ! from [GeV] to [MeV]
            s(fluka_nrecv)            = ( el - (fluka_pc0/fluka_e0)*(flt*fluka_clight) ) * c1e3 ! from [s] to [mm]
            aa(fluka_nrecv)           = flaa          !PH for hiSix
            zz(fluka_nrecv)           = flzz          !PH for hiSix
            mass(fluka_nrecv)         = flm  * c1e3  ! from [GeV] to [MeV]         !PH for hiSix
            qq(fluka_nrecv)           = flq
            pdg_id(fluka_nrecv)       = flpdgid
            spinx(fluka_nrecv)        = flsx
            spiny(fluka_nrecv)        = flsy
            spinz(fluka_nrecv)        = flsz

!            The conversion is now done inside the coupling server
!            call GetPDGid_fromFLUKA(-2, pdg_id(fluka_nrecv), flaa, flzz)
      end if

      !Finished waiting end of turn
    end do

    napx = fluka_nrecv
    fluka_last_rcvd_mess = FLUKA_EOB

    write(fluka_log_unit,*) "# FlukaIO: turn = ", turn, &
      " ipt = ", ipt, &
      " sent = ", fluka_nsent, &
      " received = ", fluka_nrecv, &
      " max_uid = ", MaximumPartID
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

  end subroutine fluka_shuffleLostParticles

  !----------------------------------------------------------------------------
  ! set reference particle properties (mainly for longitudinal dynamics)
  integer function fluka_set_synch_part( e0, pc0, mass0, a0, z0, q0 )
    implicit none

    ! interface variables
    real(kind=fPrec) :: e0, pc0, mass0
    integer(kind=int16) :: a0, z0, q0

    ! Auxiliary variables
    integer(kind=int32) :: n

    fluka_set_synch_part = 0

    fluka_e0    = e0    *c1m3 ! from  [MeV]    to [GeV]
    fluka_pc0   = pc0   *c1m3 ! from  [MeV/c]  to [GeV/c]
    fluka_mass0 = mass0 *c1m3 ! from  [MeV/c2] to [GeV/c2]
    fluka_a0 = a0
    fluka_z0 = z0
    fluka_chrg0 = q0

    write(fluka_log_unit,*) ' updated synch part:'
    write(fluka_log_unit,*) ' - total en    [GeV]:',fluka_e0
    write(fluka_log_unit,*) ' - momentum  [GeV/c]:',fluka_pc0
    write(fluka_log_unit,*) ' - mass     [GeV/c2]:',fluka_mass0
    write(fluka_log_unit,*) ' - A mass number  []:',fluka_a0
    write(fluka_log_unit,*) ' - Z number       []:',fluka_z0
    write(fluka_log_unit,*) ' - charge state  [e]:',fluka_chrg0
    flush(fluka_log_unit)

    ! update magnetic rigidity, unless division by clight and 10^-9
    fluka_brho0 = fluka_pc0 / real(fluka_chrg0,real64)

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

    use mod_common_main, only : MaximumPartID

    implicit none

    ! interface variables
    integer(kind=int32) :: npart

    ! Auxiliary variables
    integer(kind=int32) :: n

    fluka_init_max_uid = 0

    MaximumPartID = npart

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

!     A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!     last modified: 17-07-2013
!     clean closure of communication with fluka and un-set mod_fluka
!     inserted in main code by the 'fluka' compilation flag
subroutine fluka_close

     use crcoall

     implicit none

     integer fluka_con

     write(lout,'(A)') 'FLUKA> call to fluka_close'
     if (fluka_enable) then
       if (fluka_last_sent_mess==FLUKA_EOB .and. fluka_last_rcvd_mess.eq.-1) then
         ! FLUKA is still crunching stuff
         write(lout,'(A)') 'FLUKA> fluka_close called while particles are still on the Fluka side'
         write(lout,'(A)') 'FLUKA>    dummy wait to have a clean closure'
         write(fluka_log_unit,'(A)') '# fluka_close called while particles are still on the Fluka side'
         write(fluka_log_unit,'(A)') '#   dummy wait to have a clean closure'
         call kernel_fluka_exit
       end if
       fluka_con = fluka_is_running()
       if(fluka_con.eq.0) then
         if( .not. fluka_connected ) then
!              temporarily connect to fluka, to properly terminate the run
           fluka_con = fluka_connect()
           if(fluka_con.lt.0) then
!                no hope to properly close the run
             write(lerr,'(A)') 'FLUKA> ERROR Unable to connect to fluka while closing the simulation:'
             write(lerr,'(A)') 'FLUKA>       please, manually kill all its instances'
             write(fluka_log_unit,'(A)') '# unable to connect to fluka while closing the simulation:'
             write(fluka_log_unit,'(A)') '#   please, manually kill all its instances'
             goto 1982
           endif
           write(lout,'(A)') 'FLUKA> Successfully connected to Fluka server (only temporarily)'
           write(fluka_log_unit,'(A)') '# Successfully connected to Fluka server (only temporarily)'
         end if
         call fluka_end
       end if
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
    write(lerr,"(a)") "FLUKA> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if
  if(nSplit == 0) return

  select case(lnSplit(1)(1:4))

  case("DEBU")
    write(lout,"(a)") "FLUKA> DEBUG enabled"
    fluka_debug = .true.

  case("LOGU")
    write(lout,"(a,i0,a)") "FLUKA> NOTE The LOGU flag is deprecated. A log unit is assigned automatically."

  case default

    if(nSplit /= 4) then
      write(lerr,"(a,i0)") "FLUKA> ERROR Exected 4 values in input line,got ",nSplit
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
      write(lerr,"(a)") "FLUKA> ERROR Un-identified SINGLE ELEMENT '"//trim(entrElem)//"'"
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
      write(lerr,"(a)") "FLUKA> ERROR Un-identified SINGLE ELEMENT '"//trim(exitElem)//"'"
      iErr = .true.
      return
    end if

    ! 3. check that the current markers have not been already flagged
    if(fluka_type(entrIdx) /= FLUKA_NONE ) then
      write(lerr,"(a)")       "FLUKA> ERROR Single element '"//trim(bez(entrIdx))//"' was alredy labelled as fluka marker."
      write(lerr,"(2(a,i0))") "FLUKA> ERROR fluka_type(entrance) = ",fluka_type(entrIdx)," at position = ",entrIdx
      iErr = .true.
      return
    end if
    if(fluka_type(exitIdx) /= FLUKA_NONE ) then
      write(lerr,"(a)")       "FLUKA> ERROR Single element '"//trim(bez(exitIdx))//"' was alredy labelled as fluka marker."
      write(lerr,"(2(a,i0))") "FLUKA> ERROR fluka_type(exit) = ",fluka_type(exitIdx)," at position = ",exitIdx
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
  do k=1, nele
!   extract the fluka geo index value, which usually will be 0 for non-insertions
    ii = fluka_geo_index(k)
    if(ii .eq. 0) then
      continue
    else

      if(fluka_type(k) /= FLUKA_NONE) then
!       this entry exists, so add it to root
        this_name = trim(adjustl(bez(k))) // C_NULL_CHAR
        call root_FLUKA_Names(ii, this_name, len_trim(this_name), fluka_type(k))
      end if

    end if
  end do

end subroutine root_FLUKA_DumpInsertions
#endif

subroutine kernel_fluka_element( nturn, i, ix )
!
!-----------------------------------------------------------------------
!     A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!     last modified: 07-02-2014
!     'transport' subroutine, for a Fluka insertion corresponding to a
!       single SINGLE ELEMENT, of non-zero length
!     inserted in main code by the 'fluka' compilation flag
!-----------------------------------------------------------------------
!
      use floatPrecision
      use numerical_constants, only : zero, one, c1e3, c1m3
      use crcoall
      use parpro
      use mod_common
      use mod_common_track
      use mod_common_main
#ifdef CR
      use coll_common, only : fort208Pos
#endif

#ifdef ROOT
      use root_output
#endif

      implicit none

!     interface variables:
      integer nturn, i, ix

!     temporary variables
      integer ret, j, k
      integer pid_q               ! ph: hisix
      save

      fluka_i = i
      fluka_ix = ix
      fluka_nturn = nturn

      if (fluka_debug) then
!        where are we?
         write(fluka_log_unit,*)'# In fluka element of type ', fluka_type(fluka_ix)
         write(fluka_log_unit,*)'#   i=', fluka_i
         write(fluka_log_unit,*)'#   ix=', fluka_ix
         write(fluka_log_unit,*)'#   fluka_geo_index=',fluka_geo_index(fluka_ix)
         write(fluka_log_unit,*)'#   eltot=',fluka_synch_length( fluka_ix )
      end if

!     PH hisix compute the number of nucleons sent to FLUKA
!     PH hisix compute the total ion energy sent to FLUKA
      nnuc0 = 0
      ien0  = zero
      do j=1,napx
        nnuc0   = nnuc0 + naa(j)
        ien0    = ien0 + ejv(j)
        ! array of particle ids sent to FLUKA
        pids(j) = partID(j)
      end do


      ret = fluka_send_receive( fluka_nturn, fluka_geo_index(fluka_ix), fluka_synch_length( fluka_ix ), &
           napx, xv1, xv2, yv1, yv2, sigmv, ejv, naa, nzz, nucm, nqq, pdgid, &
           partID, parentID, partWeight, spin_x, spin_y, spin_z )

      if (ret.lt.0) then
         write(lerr,'(A,i0,A)')'FLUKA> ERROR ', ret, ' in Fluka communication returned by fluka_send_receive...'
         write(fluka_log_unit,'(A,i0,A)')'# Error ', ret, ' in Fluka communication returned by fluka_send_receive...'
         call prror
      end if

      nnuc1 = 0                 ! hisix: number of nucleons leaving the collimato
      ien1  = zero              ! hisix: total energy leaving the collimator
!     particles to be tracked
      do j=1,napx
!        Update values related to losses
         partID(j) = j
         pstop (j) = .false.
!        Update variables depending on total energy
!         ejfv  (j) = sqrt((ejv(j)-pma)*(ejv(j)+pma))
!         rvv   (j) = (ejv(j)*e0f)/(e0*ejfv(j))
!         dpsv  (j) = (ejfv(j)-e0f)/e0f
!         oidpsv(j) = one/(one+dpsv(j))
         ejfv  (j) = sqrt((ejv(j)-nucm(j))*(ejv(j)+nucm(j)))   ! hisix: ion mass
         rvv   (j) = (ejv(j)*e0f)/(e0*ejfv(j))                 ! hisix: remains unchanged
         dpsv  (j) = (ejfv(j)*(nucm0/nucm(j))-e0f)/e0f         ! hisix: new delta
         oidpsv(j) = one/(one+dpsv(j))
         dpsv1 (j) = (dpsv(j)*c1e3)*oidpsv(j)
         if(nqq(j) .eq. 0) then
           mtc (j) = zero
         else
           mtc (j) = (nqq(j)*nucm0)/(qq0*nucm(j))              ! hisix: mass to charge
         endif
         moidpsv (j) = mtc(j)*oidpsv(j)                        ! hisix
         omoidpsv(j) = c1e3*((one-mtc(j))*oidpsv(j))           ! hisix
         nnuc1       = nnuc1 + naa(j)                          ! outcoming nucleons
         ien1        = ien1  + ejv(j)                          ! outcoming energy
      end do

!     hisix: compute the nucleon and energy difference
!              reduce by factor 1e-3 to get the energy in GeV
      if((ien0-ien1).gt.one) then
        write(unit208,"(2(i6,1x),e24.16)") fluka_geo_index(fluka_ix), nnuc0-nnuc1, c1m3*(ien0-ien1)
#ifdef CR
        fort208Pos = fort208Pos + 1
#endif
#ifdef ROOT
        if(root_flag .and. root_FLUKA .eq. 1) then
          call root_EnergyDeposition(fluka_geo_index(fluka_ix), nnuc0-nnuc1, c1m3*(ien0-ien1))
        end if
#endif
      end if

!     hisix: check which particle ids have not been sent back
!            write their ids to fort.209
      do j=1,npart                                             ! loop over all pids possible
        pid_q = zero

        do k=1,napx                                            ! loop over pids received
          if(pids(j).eq.partID(k)) then
            pid_q = one
          end if
        end do
      end do

!     empty places
      do j=napx+1,npart
!        Update values related to losses
         partID(j) = j
         pstop (j) = .true.
!        Update values related to momentum
         ejv   (j) = zero
         rvv   (j) = one
         ejfv  (j) = zero
         dpsv  (j) = zero
         oidpsv(j) = one
         dpsv1 (j) = zero
         mtc   (j) = one            ! hiSix
         naa   (j) = aa0            ! hiSix
         nzz   (j) = zz0            ! hiSix
         nqq   (j) = qq0
         pdgid (j) = pdgid0
         spin_x(j) = zero
         spin_y(j) = zero
         spin_z(j) = zero
         nucm  (j) = nucm0          ! hiSix
         moidpsv (j) = one          ! hiSix
         omoidpsv(j) = zero         ! hiSix
      end do

!     au revoir:
      fluka_i = -1
      fluka_ix = -1
      fluka_nturn = -1
      flush(unit208)
      return
end subroutine kernel_fluka_element

subroutine kernel_fluka_entrance( nturn, i, ix )
!
!-----------------------------------------------------------------------
!     A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!     last modified: 07-02-2014
!     'transport' subroutine, for the marker starting a Fluka insertion
!     inserted in main code by the 'fluka' compilation flag
!-----------------------------------------------------------------------
!
      use floatPrecision
      use numerical_constants, only : zero
      use crcoall
      use parpro
      use mod_common
      use mod_common_track
      use mod_common_main

      implicit none


!     interface variables:
      integer nturn, i, ix, ixt

!     temporary variables
      integer ret, j
      integer k                   ! ph: hisix
      integer pid_q               ! ph: hisix

      save

!     keep track of exit element
!     NB: check_coupling_integrity and check_coupling_start_point have
!         already verify the conditions for which this search is always successful
      do j=i+1,iu
        if(ktrack(j).ne.1.and.ic(j).gt.nblo) then
          ixt=ic(j)-nblo
          if(fluka_geo_index(ix).eq.fluka_geo_index(ixt))then
             fluka_i = j
             fluka_ix = ixt
             fluka_nturn = nturn
             exit
          end if
        end if
      end do

      if (fluka_debug) then
!        where are we?
         write(fluka_log_unit,*)'# In fluka element of type ', fluka_type(ix)
         write(fluka_log_unit,*)'#   i=', i
         write(fluka_log_unit,*)'#   ix=', ix
         write(fluka_log_unit,*)'#   fluka_geo_index=',fluka_geo_index(ix)
         write(fluka_log_unit,*)'#   eltot=',zero
      end if

      ! P. HERMES for hisix
      ! send also A,Z,m to FLUKA

!     PH hisix compute the number of nucleons sent to FLUKA
!     PH hisix compute ion energy sent to FLUKA
!     PH initialize array of particle ids
      nnuc0 = 0
      ien0  = zero
      do j=1,npart
        pids(j) = 0
      end do

      do j=1,napx
        nnuc0   = nnuc0 + naa(j)
        ien0    = ien0  + ejv(j)
        pids(j) = partID(j)   ! array of particle ids sent to FLUKA
!    write(*,*),'PH:',pids(j)
      end do

      ret = fluka_send( fluka_nturn, fluka_geo_index(fluka_ix), zero, &
           napx, xv1, xv2, yv1, yv2, sigmv, ejv, naa, nzz, nucm, nqq, pdgid, &
           partID, parentID, partWeight, spin_x, spin_y, spin_z )

      if (ret.lt.0) then
         write(lerr,'(A,i0,A)')'FLUKA> ERROR ', ret,' in Fluka communication returned by fluka_send...'
         write(fluka_log_unit,'(A,i0,A)')'# Error ', ret, ' in Fluka communication returned by fluka_send...'
         call prror
      end if

!     au revoir:
      return
end subroutine kernel_fluka_entrance

subroutine kernel_fluka_exit
!
!-----------------------------------------------------------------------
!     A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!     last modified: 07-02-2014
!     'transport' subroutine, for the marker closing a Fluka insertion
!     inserted in main code by the 'fluka' compilation flag
!-----------------------------------------------------------------------
!
      use floatPrecision
      use numerical_constants, only : zero, one, c1e3, c1m3
      use crcoall
      use parpro
      use mod_common
      use mod_common_track
      use mod_common_main
#ifdef CR
      use coll_common, only : fort208Pos
#endif
#ifdef ROOT
      use root_output
#endif

      implicit none

!     interface variables:
      integer nturn, i, ix

!     temporary variables
      integer ret, j, k
      integer pid_q               ! ph: hisix

      save

      if (fluka_debug) then
!        where are we?
         write(fluka_log_unit,*)'# In fluka element of type ', fluka_type(fluka_ix)
         write(fluka_log_unit,*)'#   i=', fluka_i
         write(fluka_log_unit,*)'#   ix=', fluka_ix
         write(fluka_log_unit,*)'#   fluka_geo_index=',fluka_geo_index(fluka_ix)
         write(fluka_log_unit,*)'#   eltot=',fluka_synch_length( fluka_ix )
      end if

      ret = fluka_receive( fluka_nturn, fluka_geo_index(fluka_ix), fluka_synch_length( fluka_ix ), &
           napx, xv1, xv2, yv1, yv2, sigmv, ejv, naa, nzz, nucm, nqq, pdgid, partID, parentID,&
           partWeight, spin_x, spin_y, spin_z )

      if (ret.lt.0) then
         write(lerr,'(A,i0,A)')'FLUKA> ERROR ', ret, ' in Fluka communication returned by fluka_receive...'
         write(fluka_log_unit,'(A,i0,A)')'# Error ',ret, ' in Fluka communication returned by fluka_receive...'
         call prror
      end if

      nnuc1 = 0                 ! hisix: number of nucleons leaving the collimator
      ien1  = zero              ! hisix: total energy leaving the collimator
!     particles to be tracked
      do j=1,napx
!        Update values related to losses
         partID(j) = j
         pstop (j) = .false.
!        Update variables depending on total energy
!         ejfv  (j) = sqrt((ejv(j)-pma)*(ejv(j)+pma))
!         rvv   (j) = (ejv(j)*e0f)/(e0*ejfv(j))
!         dpsv  (j) = (ejfv(j)-e0f)/e0f
!         oidpsv(j) = one/(one+dpsv(j))
!         dpsv1 (j) = (dpsv(j)*c1e3)*oidpsv(j)
         ejfv  (j) = sqrt((ejv(j)-nucm(j))*(ejv(j)+nucm(j)))   ! hisix: ion mass
         rvv   (j) = (ejv(j)*e0f)/(e0*ejfv(j))                 ! hisix: remains unchanged
         dpsv  (j) = (ejfv(j)*(nucm0/nucm(j))-e0f)/e0f         ! hisix: new delta
         oidpsv(j) = one/(one+dpsv(j))
         dpsv1 (j) = (dpsv(j)*c1e3)*oidpsv(j)
         mtc     (j) = (nqq(j)*nucm0)/(qq0*nucm(j))            ! hisix: mass to charge
         moidpsv (j) = mtc(j)*oidpsv(j)                        ! hisix
         omoidpsv(j) = c1e3*((one-mtc(j))*oidpsv(j))           ! hisix
         nnuc1       = nnuc1 + naa(j)                          ! outcoming nucleons
         ien1        = ien1  + ejv(j)                          ! outcoming energy
      end do

!       hisix: compute the nucleon and energy difference
!              reduce by factor 1e-3 to get the energy in GeV
        if((ien0-ien1).gt.one) then
          write(unit208,"(2(i6,1x),e24.16)") fluka_geo_index(fluka_ix), nnuc0-nnuc1, c1m3*(ien0-ien1)
#ifdef CR
          fort208Pos = fort208Pos + 1
#endif
#ifdef ROOT
          if(root_flag .and. root_FLUKA .eq. 1) then
            call root_EnergyDeposition(fluka_geo_index(fluka_ix), nnuc0-nnuc1, c1m3*(ien0-ien1))
          end if
#endif
        end if
!
!     hisix: check which particle ids have not been sent back
!            write their ids to fort.209
      do j=1,npart                                       ! loop over all pids possible
        pid_q = zero
        do k=1,napx                                    ! loop over pids received
          if(pids(j).eq.partID(k)) then
            pid_q = one
          end if
        end do
      end do

!     empty places
      do j=napx+1,npart
!        Update values related to losses
         partID(j) = j
         pstop (j) = .true.
!        Update values related to momentum
         ejv   (j) = zero
         rvv   (j) = one
         ejfv  (j) = zero
         dpsv  (j) = zero
         oidpsv(j) = one
         dpsv1 (j) = zero
         mtc   (j) = one            ! hiSix
         naa   (j) = aa0            ! hiSix
         nzz   (j) = zz0            ! hiSix
         nqq   (j) = qq0
         pdgid (j) = pdgid0
         spin_x(j) = zero
         spin_y(j) = zero
         spin_z(j) = zero
         nucm  (j) = nucm0          ! hiSix
         moidpsv (j) = one          ! hiSix
         omoidpsv(j) = zero         ! hiSix
      end do

!     au revoir:
      fluka_i = -1
      fluka_ix = -1
      fluka_nturn = -1
      flush(unit208)
      return
end subroutine kernel_fluka_exit

subroutine check_coupling_integrity
!
!-----------------------------------------------------------------------
!     A.Mereghetti, for the FLUKA Team
!     last modified: 23-05-2019
!     check that an entrance MARKER is followed by an exit one in the
!        accelerator sequence, even though not immediately
!     NB: together with check_coupling_start_point, this subroutine is
!           fundamental to verify that every entrance point has an exit
!           one, downstream of it, and the FLUKA insertion is not accross
!           the lattice extremities
!     inserted in main code by the 'fluka' compilation flag
!-----------------------------------------------------------------------

      use floatPrecision
      use crcoall
      use parpro
      use mod_common
      use mod_common_track
      implicit none

!     temporary variables
      integer i1 , i2
      integer ix1, ix2
      integer istart, istop
      logical lerror, lfound, lcurturn

      lerror = .false.

      write(lout,'(a)') 'FLUKA> '
      write(lout,10040)
      write(lout,'(a)') 'FLUKA> CALL TO CHECK_COUPLING_INTEGRITY '
      write(lout,'(a)') 'FLUKA> NB: only entrance/exit markers are listed;'
      write(lout,'(a)') 'FLUKA>     a single entry is by definition righteous'
      write(lout,10040)
      write(lout,'(a)') 'FLUKA> '
      write(lout,'(a)') 'FLUKA>         keys to FLUKA types:'
      write(lout,'(a,i0,a)') 'FLUKA> ',FLUKA_ELEMENT,' --> simple element'
      write(lout,'(a,i0,a)') 'FLUKA> ',FLUKA_ENTRY,' --> entrance point'
      write(lout,'(a,i0,a)') 'FLUKA> ',FLUKA_EXIT,' --> exit point'
      write(lout,'(a)') 'FLUKA> '
      write(lout,10040)
      write(lout,10020) 'entry type', 'name', 'ID SING EL', 'ID struct', 'ID geom'

      i1=1
      do while ( i1.le.iu )
        if(ktrack(i1).ne.1.and.ic(i1).gt.nblo) then
!         SINGLE ELEMENT
          ix1=ic(i1)-nblo
          if ( fluka_type(ix1).eq.FLUKA_ENTRY ) then
            write(lout,10040)
            write(lout,10030) fluka_type(ix1), bez(ix1), ix1, i1, fluka_geo_index(ix1)
            istart = i1+1
            istop  = iu
            lcurturn = .true.
            lfound = .false.
 1982       continue
            do i2=istart,istop
              if(ktrack(i2).ne.1.and.ic(i2).gt.nblo) then
!               SINGLE ELEMENT
                ix2=ic(i2)-nblo
                if ( fluka_type(ix2).eq.FLUKA_EXIT ) then
                  if(fluka_geo_index(ix1).eq.fluka_geo_index(ix2))then
                    write(lout,10030) fluka_type(ix2), bez(ix2), ix2, i2, fluka_geo_index(ix2)
                    i1 = i2
                    lfound = .true.
                    if ( lcurturn ) then
                      exit
                    else
                      goto 1983
                    endif
                  else
                    write(lerr,"(a)") "FLUKA> ERROR Un-matched geo index"
                    write(lerr,10030) fluka_type(ix2), bez(ix2), ix2, i2, fluka_geo_index(ix2)
                    lerror = .true.
                  endif
                elseif ( fluka_type(ix2).ne.FLUKA_NONE ) then
                  write(lerr,"(a)") "FLUKA> ERROR Non-exit point when entrance is on"
                  write(lerr,10030) fluka_type(ix2), bez(ix2), ix2, i2, fluka_geo_index(ix2)
                  lerror = .true.
                endif
              endif
            enddo
            if ( .not. lfound ) then
              if ( lcurturn ) then
!               the exit marker was not found: maybe the fluka insertion
!                  is across the end/beginning of the accelerator structure;
!               -> restart the research, upstream of the entrance marker:
                istart = 1
                istop  = i1
                lcurturn = .false.
                goto 1982
              else
!               failing research:
!               NB: in principle, this should never happen, but let's be picky
                write(lerr,"(a)") "FLUKA> ERROR Entrance point does not have the exit"
                lerror = .true.
              endif
            endif
          endif
        endif

!       go to next accelerator entry
        i1 = i1+1
      enddo
      write(lout,10040)

 1983 continue
      if ( lerror ) then
        write(lout,'(a)') ' at least one inconsistency in flagging elements'
        write(lout,'(a)') '    for coupling: please check carefully...'
        call prror
      endif

!     au revoir:
      return
10020 format('FLUKA> ',1X,A10,1X,A4,12X,3(1X,A10))
10030 format('FLUKA> ',1X,I10,1X,A16,3(1X,I10))
10040 format('FLUKA> ',62('-'))
end subroutine check_coupling_integrity

subroutine check_coupling_start_point()
!-----------------------------------------------------------------------
!     A.Mereghetti, CERN BE-ABP-HSS
!     last modified: 20-03-2019
!     check that the lattice structure (after re-shiffle due to GO
!        statement) does not start inside a FLUKA insertion region
!     NB: together with check_coupling_integrity, this subroutine is
!           fundamental to verify that every entrance point has an exit
!           one, downstream of it, and the FLUKA insertion is not accross
!           the lattice extremities
!-----------------------------------------------------------------------

  use parpro, only : nblo, nele
  use mod_common, only : iu, ic, bez
  use crcoall, only : lout, lerr
  use mod_common_track, only : ktrack
!  use mod_fluka, only : FLUKA_ELEMENT, FLUKA_ENTRY, FLUKA_EXIT, &
!           fluka_type, fluka_geo_index

  implicit none

! temporary variables
  integer ii, ix, iInside, jj

  iInside=-1
  do ii=1,iu
    if(ktrack(ii).ne.1.and.ic(ii).gt.nblo) then
      ! SINGLE ELEMENT
      ix=ic(ii)-nblo
      if ( fluka_type(ix).eq.FLUKA_EXIT ) then
        write(lerr,"(a,i0)") "FLUKA> ERROR Lattice structure starts inside FLUKA insertion region # ",fluka_geo_index(ix)
        do jj=1,nele
          if ( fluka_geo_index(ix).eq.fluka_geo_index(jj).and.fluka_type(jj).eq.FLUKA_ENTRY ) then
            write(lerr,"(a,i0)") "FLUKA>       entrance marker: "//trim(bez(jj))//" - exit marker: "//trim(bez(ix))
            exit
          end if
        end do
        write(lerr,"(a,i0)") "FLUKA>       The actual lattice starting point should be outside a FLUKA insergion region"
        write(lerr,"(a,i0)") "FLUKA>       Please update your lattice structure or set the GO in a sensible position"
        iInside=fluka_geo_index(ix)
        call prror
        exit
      elseif ( fluka_type(ix).eq.FLUKA_ENTRY .or. fluka_type(ix).eq.FLUKA_ELEMENT ) then
        write(lout,"(a)") ""
        write(lout,"(a,i0)") "FLUKA> Lattice structure starts upstream of FLUKA insertion region #",fluka_geo_index(ix)
        write(lout,"(a)") ""
        iInside=fluka_geo_index(ix)
        exit
      end if
    end if
  end do
  if ( iInside==-1 ) then
    write(lout,"(a)") ""
    write(lout,"(a,i0)") "FLUKA> No FLUKA insertion region found!"
    write(lout,"(a)") ""
  end if

! au revoir:
  return

end subroutine check_coupling_start_point

end module mod_fluka

