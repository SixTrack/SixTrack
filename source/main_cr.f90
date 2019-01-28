!-----------------------------------------------------------------------
!
!  SIXTRACK
!
!  SIXDIMENSIONAL PARTICLE-TRACKING
!
!-----------------------------------------------------------------------
!
!  F. SCHMIDT, M. VANTTINEN
!
!  COLLIMATION VERSION, NOVEMBER 2004
!
!  G. ROBERT-DEMOLAIZE
!
!  COLLIMATION UPGRADE, JUNE 2005
!
!  G. ROBERT-DEMOLAIZE, S. REDAELLI
!
!  UPGRADED FOR COUPLING TO FLUKA, JULY 2013
!
!  A. MEREGHETTI, D. SINUELA PASTOR
!
!  FURTHER UPGRADE FOR COUPLING TO FLUKA, MAY-JUNE 2014
!
!  A. MEREGHETTI, P. GARCIA ORTEGA
!
!-----------------------------------------------------------------------
!  SIXTRACR CHECKPOINT/RESTART and CRLIBM (ENS Lyon)
!
!  E. MCINTOSH FEBRUARY 2005
!-----------------------------------------------------------------------
!  USED DISKS:
!
!  GEOMETRY AND STRENGTH OF THE ACCELERATOR : UNIT  2
!  TRACKING PARAMETER                       : UNIT  3
!  NORMAL PRINTOUT                          : UNIT  6
!  TRACKING DATA                            : UNIT  8
!  DATA FOR SUMMARY OF THE POSTPROCESSING   : UNIT 10
!  AUXILIARY FILE FOR THE INPUT             : UNIT 11
!  ASCII FILE WITH THE HORIZONTAL FFT DATA  : UNIT 14
!  ASCII FILE WITH THE VERTICAL FFT DATA    : UNIT 15
!  METAFILE FOR PLOTTING WITH GKS           : UNIT 20
!
!  FOR CR VERSION:
!  CHECKPOINT/RESTART FILES                 : UNIT 95,96
!  OPTIONAL DUMP.DEBUG FILE                 : UNIT 99
!  PROGRESS FILE                            : UNIT 91
!  INTERMEDIATE OUTPUT FILE (LOUT)          : UNIT 92
!  CHECKPOINT/RESTART LOGFILE               : UNIT 93
!  TEMPORARY SCRATCH FILE for C/R           : UNIT 94
!-----------------------------------------------------------------------

program maincr

  use floatPrecision
  use file_units
  use string_tools
  use mathlib_bouncer
  use physical_constants
  use numerical_constants

  use scatter, only : scatter_active, scatter_initialise
  use dynk,    only : dynk_izuIndex
  use fma,     only : fma_postpr, fma_flag
  use dump,    only : dump_initialise, dumpclo,dumptas,dumptasinv
  use zipf,    only : zipf_numfiles, zipf_dozip

  use, intrinsic :: iso_fortran_env, only : output_unit
  use mod_units
  use aperture
  use mod_ranecu
  use mod_alloc,      only : alloc_init
  use mod_fluc,       only : fluc_randomReport, fluc_errAlign, fluc_errZFZ
  use postprocessing, only : postpr, writebin_header, writebin
  use read_input,     only : readFort33

#ifdef FLUKA
  use mod_fluka
#endif
#ifdef FFIELD
  ! Modification by B.DALENA and T.PUGNAT
  use mod_ffield,     only :ffield_mod_init,ffield_mod_end
#endif
#ifdef HDF5
  use hdf5_output
#endif
#ifdef ROOT
  use root_output
#endif
#ifdef CR
  use checkpoint_restart
#endif

  use crcoall
  use parpro
  use mod_settings
  use mod_common
  use mod_commonmn
  use mod_commons
  use mod_commont

  use mod_hions
  use mod_dist
  use matrix_inv
  use aperture
  use wire

  implicit none

interface

  subroutine envarsv(dpsv,oidpsv,rvv,ekv)

    use floatPrecision
    use parpro
    use mod_commond

    implicit none

    real(kind=fPrec) :: dpsv(npart)
    real(kind=fPrec) :: oidpsv(npart)
    real(kind=fPrec) :: rvv(npart)
    real(kind=fPrec), allocatable, intent(inout) :: ekv(:,:)
  end subroutine envarsv

end interface

  ! "Old" variables
  integer i,itiono,i2,i3,ia,ia2,iar,iation,ib,ib0,ib1,ib2,ib3,id,ie,ig,ii,im,iposc,ix,izu,j,j2,jj,  &
    k,kpz,kzz,l,ll,m,ncorruo,ncrr,nd,nd2,ndafi2,nerror,nlino,nlinoo,nmz,nthinerr
  real(kind=fPrec) alf0s1,alf0s2,alf0s3,alf0x2,alf0x3,alf0z2,alf0z3,amp00,bet0s1,bet0s2,bet0s3,     &
    bet0x2,bet0x3,bet0z2,bet0z3,chi,coc,dam1,dchi,ddp1,dp0,dp00,dp10,dpsic,dps0,dsign,gam0s1,gam0s2,&
    gam0s3,gam0x1,gam0x2,gam0x3,gam0z1,gam0z2,gam0z3,phag,r0,r0a,rat0,sic,tasia56,tasiar16,tasiar26,&
    tasiar36,tasiar46,tasiar56,tasiar61,tasiar62,tasiar63,tasiar64,tasiar65,taus,x11,x13
  integer idummy(6)
  character(len=4) cpto
#ifdef CR
  logical isOpen
#endif
  character(len=8) cdate,ctime,progrm !Note: Keep in sync with writebin_header and more
                                      !DANGER: If the len changes, CRCHECK will break.
#ifdef BOINC
  character(len=256) filename
#endif
#ifdef CRLIBM
  integer nchars
  parameter (nchars=160)
  character(len=nchars) ch
  character(len=nchars+nchars) ch1
  real(kind=fPrec) round_near
#endif
#ifdef FLUKA
  integer fluka_con
#endif

  ! New Variables
  character(len=:), allocatable :: featList
#ifndef STF
  character(len=7)  tmpFile
#endif
  character(len=23) timeStamp
  character(len=8)  tsDate
  character(len=10) tsTime

  logical fErr ! For file units

#include "version.f90"

  ! ---------------------------------------------------------------------------------------------- !
  errout_status = 0 ! Set to nonzero before calling abend in case of error.
  lout = 92
#ifndef CR
  lout = output_unit
#endif

  call funit_initUnits ! This one has to be first
  call units_initUnits
  call alloc_init      ! Initialise tmod_alloc
  call allocate_arrays ! Initial allocation of memory

  ! Set napx,napxo,trtime for error handling
  napx   = 0
  napxo  = 0
  trtime = 0.0
  napxto = 0

  !----------------------------------------------------------------------------------------------- !
  ! Features
  featList = ""
#ifdef TILT
  featList = featList//" TILT"
#endif
#ifdef FAST
  featList = featList//" FAST"
#endif
#ifdef STF
  featList = featList//" STF"
#endif
#ifdef COLLIMAT
  featList = featList//" COLLIMAT"
#endif
#ifdef CRLIBM
  featList = featList//" CRLIBM"
  call disable_xp()
#endif
#ifdef FIO
  featList = featList//" FIO"
#endif
#ifdef CR
  featList = featList//" CR"
  stxt = ""
#endif
#ifdef ROOT
  featList = featList//" ROOT"
#endif
#ifdef HDF5
  featList = featList//" HDF5"
  call h5_initHDF5()
#endif
#ifdef BOINC
  featList = featList//" BOINC"
  call boinc_init()
! call boinc_init_graphics()
#endif
#ifdef LIBARCHIVE
  featList = featList//" LIBARCHIVE"
#endif

#ifdef CR
  ! Main start for Checkpoint/Restart
  sythckcr = .false.
  numlcr   = 1
  rerun    = .false.
  start    = .true.
  restart  = .false.
  checkp   = .false.
  fort95   = .false.
  fort96   = .false.
  sixrecs  = 0
  binrec   = 0
  bnlrec   = 0
  bllrec   = 0
  crtime3  = 0.0
  ! do i=1,(npart+1)/2
  !   binrecs(i) = 0
  ! end do

#ifdef BOINC
611 continue
#endif
  ! Very first get rid of any previous partial output
  inquire(unit=lout, opened=isOpen)
  if(isOpen) close(lout)
  call units_openUnit(unit=lout,fileName="fort.92",formatted=.true.,mode="w",err=fErr,status="replace")

  ! Now position the checkpoint/restart logfile=93
  call units_openUnit(unit=93,fileName="fort.93",formatted=.true.,mode="w",err=fErr)
606 continue
  read(93,"(a1024)",end=607) arecord
  goto 606
607 continue
  backspace(93,iostat=ierro)
#ifdef BOINC
  ! and if BOINC issue an informatory message
  if(start) then
    write(93,"(a)") "SIXTRACR starts for the very first time"
  else
    write(93,"(a)") "SIXTRACR retry after unzip of Sixin.zip"
  end if
#endif
  ! Now we see if we have a fort.6 which implies that we can perhaps just restart using all exisiting files
  ! including the last checkpoints. If not, we just do a start (with an unzip for BOINC)
  ! call units_openUnit(unit=6,fileName="fort.6",formatted=.true.,mode="w",err=fErr,status="old")
  ! if(fErr) goto 602
  ! stxt = "SIXTRACR reruns on: "
  call units_openUnit(unit=output_unit,fileName="fort.6",formatted=.true.,mode="w",err=fErr,status="old")
  if(fErr) then
#ifdef BOINC
    ! No fort.6 so we do an unzip of Sixin.zip
    ! BUT ONLY IF WE HAVE NOT DONE IT ALREADY
    ! and CLOSE 92 and 93
    if(start) then
      start=.false.
      close(92)
      close(93)
      ! Now, if BOINC, after no fort.6, call UNZIP Sixin.zip
      ! Name hard-wired in our boinc_unzip_.
      ! Either it is only the fort.* input data or it is a restart.
      call boincrf("Sixin.zip",filename)
      ! This function expects a normal, trimmed fortran string; it will do the zero-padding internally.
      call f_read_archive(trim(filename),".")
      goto 611
    end if
    call units_openUnit(unit=output_unit,fileName="fort.6",formatted=.true.,mode="w",err=fErr)
#else
    call units_openUnit(unit=output_unit,fileName="fort.6",formatted=.true.,mode="w",err=fErr,status="new")
#endif
    ! Set up start message depending on fort.6 or not
    stxt = "SIXTRACR starts on: "
  else
    ! Set up start message depending on fort.6 or not
    stxt = "SIXTRACR reruns on: "
    rerun=.true.
  end if
  call units_openUnit(unit=95,fileName="fort.95",formatted=.false.,mode="w",err=fErr,status="old")
  if(fErr) then
    call units_openUnit(unit=95,fileName="fort.95",formatted=.false.,mode="w",err=fErr,status="new")
  else
    fort95 = .true.
  end if
  call units_openUnit(unit=96,fileName="fort.96",formatted=.false.,mode="w",err=fErr,status="old")
  if(fErr) then
    call units_openUnit(unit=96,fileName="fort.96",formatted=.false.,mode="w",err=fErr,status="new")
  else
    fort96 = .true.
  end if
  call units_openUnit(unit=91,fileName="fort.91",formatted=.true.,mode="w",err=fErr)
#else
  lout = output_unit
#endif

  ! Open Regular File Units
  call units_openUnit(unit=2, fileName="fort.2", formatted=.true., mode="r", err=fErr) ! Should be opened in DATEN
  call units_openUnit(unit=3, fileName="fort.3", formatted=.true., mode="r", err=fErr) ! Should be opened in DATEN
! call units_openUnit(unit=4, fileName="fort.4", formatted=.true., mode="w", err=fErr) ! Handled by mod_fluc
  call units_openUnit(unit=7, fileName="fort.7", formatted=.true., mode="w", err=fErr,recl=303)
! call units_openUnit(unit=8, fileName="fort.8", formatted=.true., mode="r", err=fErr) ! Handled by mod_fluc
  call units_openUnit(unit=9, fileName="fort.9", formatted=.true., mode="w", err=fErr)
  call units_openUnit(unit=11,fileName="fort.11",formatted=.true., mode="w", err=fErr)
  call units_openUnit(unit=12,fileName="fort.12",formatted=.true., mode="w", err=fErr)
  call units_openUnit(unit=13,fileName="fort.13",formatted=.true., mode="r", err=fErr) ! Should only be opened when reading
  call units_openUnit(unit=14,fileName="fort.14",formatted=.true., mode="w", err=fErr)
  call units_openUnit(unit=15,fileName="fort.15",formatted=.true., mode="w", err=fErr)
! call units_openUnit(unit=16,fileName="fort.16",formatted=.true., mode="r", err=fErr) ! Handled by mod_fluc
! call units_openUnit(unit=17,fileName="fort.17",formatted=.true., mode="w", err=fErr) ! Not in use? Should mirror fort.16
  call units_openUnit(unit=18,fileName="fort.18",formatted=.true., mode="w", err=fErr)
! call units_openUnit(unit=19,fileName="fort.19",formatted=.true., mode="rw",err=fErr) ! Not in use?
  call units_openUnit(unit=20,fileName="fort.20",formatted=.true., mode="w", err=fErr)
  call units_openUnit(unit=21,fileName="fort.21",formatted=.true., mode="w", err=fErr)
! call units_openUnit(unit=22,fileName="fort.22",formatted=.true. ,mode="w", err=fErr) ! Not in use?
! call units_openUnit(unit=23,fileName="fort.23",formatted=.true., mode="w", err=fErr) ! Not in use?
! call units_openUnit(unit=24,fileName="fort.24",formatted=.true., mode="w", err=fErr) ! Not in use?
! call units_openUnit(unit=25,fileName="fort.25",formatted=.true., mode="w", err=fErr) ! Not in use?
! call units_openUnit(unit=26,fileName="fort.26",formatted=.true., mode="w", err=fErr) ! Not in use?
  call units_openUnit(unit=27,fileName="fort.27",formatted=.true., mode="w", err=fErr)
  call units_openUnit(unit=28,fileName="fort.28",formatted=.true., mode="w", err=fErr)
  call units_openUnit(unit=29,fileName="fort.29",formatted=.true., mode="w", err=fErr)
! call units_openUnit(unit=30,fileName="fort.30",formatted=.true., mode="r", err=fErr) ! Handled by mod_fluc
  call units_openUnit(unit=31,fileName="fort.31",formatted=.true., mode="w", err=fErr)
  call units_openUnit(unit=32,fileName="fort.32",formatted=.false.,mode="w", err=fErr)
  call units_openUnit(unit=34,fileName="fort.34",formatted=.true., mode="w", err=fErr)
! call units_openUnit(unit=35,fileName="fort.35",formatted=.true., mode="w", err=fErr) ! Not in use?

#ifdef STF
  ! Open Single Track File
  call units_openUnit(unit=90,fileName="singletrackfile.dat",formatted=.false.,mode="w",err=fErr)
#else
  ! Open binary files 59 to 90 for particle pair 1 to 32
  do i=59,90
    write(tmpFile,"(a5,i2)") "fort.",i
    call units_openUnit(unit=i,fileName=tmpFile,formatted=.false.,mode="w",err=fErr)
  end do
#endif

  call units_openUnit(unit=98,fileName="fort.98",formatted=.true.,mode="w",err=fErr)

  ! Eric for the DA coefficients in BINARY
  call units_openUnit(unit=110,fileName="fort.110",formatted=.false.,mode="w",err=fErr)
  call units_openUnit(unit=111,fileName="fort.111",formatted=.false.,mode="w",err=fErr)

#ifdef DEBUG
  call units_openUnit(unit=99 ,fileName="dump",  formatted=.false.,mode="w",err=fErr)
  call units_openUnit(unit=100,fileName="arrays",formatted=.false.,mode="w",err=fErr)
#endif

  ! Heavy Ion Output
  call units_openUnit(unit=208,fileName="fort.208",formatted=.true.,mode="w",err=fErr) ! coll losses (energy)
  call units_openUnit(unit=209,fileName="fort.209",formatted=.true.,mode="w",err=fErr) ! coll losses in function of particle i
  call units_openUnit(unit=210,fileName="fort.210",formatted=.true.,mode="w",err=fErr) ! mtc after each collimator interaction

  ! ---------------------------------------------------------------------------------------------- !
  ! Write Header

  ! TimeStamp
  call date_and_time(tsDate,tsTime)
  timeStamp = tsDate(1:4)//"-"//tsDate(5:6)//"-"//tsDate(7:8)//" "//&
              tsTime(1:2)//":"//tsTime(3:4)//":"//tsTime(5:10)

  write(lout,"(a)") ""
  write(lout,"(a)") "    SixTrack :: Version "//trim(version)//" :: Released "//trim(moddate)
  write(lout,"(a)") "  "//repeat("=",128)
  write(lout,"(a)") "    Git SHA Hash: "//trim(git_revision)
  write(lout,"(a)") "    Built With:  "//featList
  write(lout,"(a)") "    Start Time:   "//timeStamp
  write(lout,"(a)") ""
  write(lout,"(a)") str_divLine

#ifdef CR
  ! Log start messages
  write(93,"(a)") ""
  write(93,"(a)") "SIXTRACR MAINCR"
  write(93,"(a)") stxt//timeStamp
  endfile(93,iostat=ierro)
  backspace(93,iostat=ierro)
#endif

!     A normal start, time0 is beginning
      pretime=0.0
      trtime=0.0
      posttime=0.0
      tottime=0.0
      time0=0.0
      time1=0.0
      time2=0.0
      time3=0.0
      tlim=1e7
      call timest
      call timex(time0)
      do 20 i=1,mmul
        cr(i)=zero
        ci(i)=zero
   20 continue
      do 30 i=1,2
        eps(i)=zero
        epsa(i)=zero
        ekk(i)=zero
        qw(i)=zero
        qwc(i)=zero
   30 continue
      qwc(3)=zero
      call comnul
      commen=' '
      progrm='SIXTRACK'

#ifdef ROOT
      call SixTrackRootFortranInit
#endif

#ifdef FLUKA
  call fluka_mod_init(npart_initial, nele_initial, clight)
#endif

#ifdef FFIELD
  ! Modification by B.DALENA and T.PUGNAT
  call ffield_mod_init(npart_initial, nele_initial)
#endif

  call daten

#ifdef HDF5
  if(h5_isActive) then
    call h5_openFile()
    call h5_writeSimInfo()
  end if
#endif
  call aperture_init

  if (ithick.eq.1) call allocate_thickarrays

#ifdef DEBUG
!     call dumpbin('adaten',999,9999)
!     call abend('after  daten                                      ')
#endif
#if defined(DEBUG) && defined(CR)
!     write(93,*) 'ERIC IL= ',il
!     endfile (93,iostat=ierro)
!     backspace (93,iostat=ierro)
#endif
#ifdef CR
      checkp=.true.
      call crcheck
#endif
      if(ithick.eq.1) write(lout,"(a)") "MAINCR> Structure input file has -thick- linear elements"
      if(ithick.eq.0) write(lout,"(a)") "MAINCR> Structure input file has -thin- linear elements"
      if(ibidu.eq.2) then
        write(lout,10025)
        goto 550
      endif

#ifndef FLUKA
  ! SETTING UP THE PLOTTING
  if(ipos.eq.1.and.(idis.ne.0.or.icow.ne.0.or.istw.ne.0.or.iffw.ne.0)) then
    call hlimit(nplo)
    call hplint(kwtype)
    call igmeta(-20,-111)
    cpto='NPTO'
    if(icr.eq.1) cpto='PTO '
    call hplopt(cpto,1)
    call hplopt('DATE',1)
    call hplset('DATE',1.)
    call hplset('CSIZ',.15)
  endif

  ! Postprocessing is on, but there are no particles
  if(ipos.eq.1.and.napx.eq.0) then
    ! Now we open fort.10 unless already opened for BOINC
    call units_openUnit(unit=10,fileName="fort.10",formatted=.true.,mode="w",err=fErr,recl=8195)

#ifndef STF
    do i=1,ndafi !ndafi = number of files to postprocess (set by fort.3)
#ifndef CR
      call postpr(91-i)
#else
      write(93,"(a,i0)") "MAINCR> Calling POSTPR nnuml = ",nnuml
      endfile(93,iostat=ierro)
      backspace(93,iostat=ierro)
      call postpr(91-i,nnuml)
#endif
    end do
#else
    ! ndafi normally set in fort.3 to be "number of files to postprocess"
    ! Inside the postpr subroutine ndafi is modified as:
    ! ndafi=itopa(total particles) if once particle per header i.e ntwin=1,
    ! ndafi=itopa/2 if 2 particle per header i.e ntwin=2
    do i=1,(2*ndafi),2
#ifndef CR
      call postpr(i)
#else
      write(93,"(a,i0)") "MAINCR> Calling POSTPR nnuml = ",nnuml
      endfile(93,iostat=ierro)
      backspace(93,iostat=ierro)
      call postpr(i,nnuml)
#endif
    end do
#endif
! END ifndef STF

    call sumpos
    goto 520 ! Jump to after particle&optics initialization, and also after tracking.
  end if !if(ipos.eq.1.and.napx.eq.0)
#endif
! END ifndef FLUKA

  do i=1,20
    fake(1,i)=zero
    fake(2,i)=zero
  end do

  itra=2
  amp00=amp(1)
  if(napx.ne.1) damp=((amp00-amp0)/real(napx-1,fPrec))/two                 !hr05
  napx=2*napx
  aperture_napxStart=napx
  iation=abs(ition)
  ib0=0
  dp00=dp1
  if(napx.le.0.or.imc.le.0) goto 490
  do 260 m=1,mmac
#ifdef DEBUG
!       call warr('mmac and m',0d0,nmac,m,0,0)
!       write(*,*) 'do 260 mmac/m',mmac,m
#endif
    !--MULTIPOLE WITH THEIR RANDOM VALUES ADDED
    ! mmac is currently not allowed to be larger than 1
    ! zfz array is n ow handled by mod_fluc, and using the below code
    ! will break tests
    ! if(m.ge.2) then
    !   call recuin(m*izu0,irecuin)
    !   call ranecu(zfz,nzfz,mcut)
    !   rsum=zero
    !   do i=1,nzfz
    !     rsum=rsum+zfz(i)
    !   end do
    !   rmean=rsum/real(nzfz,fPrec)
    !   rsqsum=zero
    !   do i=1,nzfz
    !     rsqsum=rsqsum+(zfz(i)-rmean)*(zfz(i)-rmean)
    !   end do
    !   rdev=sqrt(rsqsum/real(nzfz,fPrec))
    !   write(lout,10320) m*izu0,nzfz,rmean,rdev
    !   write(lout,10070)
    ! endif

    ! A.Mereghetti (CERN, BE-ABP-HSS), 06-03-2018
    ! possible to re-shuffle lattice structure
    if(m.eq.1) call orglat

    ! A.Mereghetti, P. G. Ortega and D.Sinuela Pastor, for the FLUKA Team
    ! last modified: 01-07-2014
    ! call routine for calculating dcum, necessary for the online
    !    aperture check and in case of dumping particle population
    !    or statistics or beam matrix
    call cadcum
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

    ! A.Mereghetti (CERN, BE-ABP-HSS), 16-12-2016
    ! initialise aperture of first and last elements of sequence
    if (limifound) then
      write(lout,"(a)") "MAINCR> Check that beginning/end of lattice structure is assigned aperture markers."
      call contour_aperture_markers( iu, 1, .false. )
    end if

#ifdef FLUKA
    if (fluka_enable) call check_coupling_integrity
#endif

    ! dump aperture model
    if (ldmpaper) call dump_aperture_model
    ! dump x-sections at specific locations
    if (mxsec.gt.0) call dump_aperture_xsecs
    ! map errors, now that the sequence is no longer going to change
    if(m.eq.1) then
      call ord
      if(allocated(zfz)) call fluc_randomReport
    end if

    call clorb(ded)

#ifdef ROOT
    if(root_flag) then
      call SixTrackRootInit()
      call ConfigurationOutputRootSet_npart(napx)
      call ConfigurationOutputRootSet_nturns(nnuml)
      call ConfigurationRootWrite()

      ! Dump the accelerator lattice
      if(root_flag .and. root_Accelerator == 1) then
        ! loop all over the entries in the accelerator structure
        do i=1,iu
          ix=ic(i)
          if(ix.gt.nblo) then
            ix=ix-nblo
            call AcceleratorRootWrite(trim(adjustl(bez(ix)))//C_NULL_CHAR,&
              len_trim(trim(adjustl(bez(ix)))//C_NULL_CHAR), kz(ix), ed(ix), ek(ix), el(ix))
          else
            do j=1,mel(ix)
              k=mtyp(ix,j)
              call AcceleratorRootWrite(trim(adjustl(bez(k)))//C_NULL_CHAR, &
                len_trim(trim(adjustl(bez(k)))//C_NULL_CHAR), kz(k), ed(k), ek(k), el(k))
            end do
          end if
        end do
      end if

#ifdef FLUKA
     !Must be called after input parsing and root configuration/init is finished.
     if(root_flag .and. root_FLUKA.eq.1) then
       call root_FLUKA_DumpInsertions
     end if
#endif

   end if
#endif

#ifdef DEBUG
!     call dumpbin('aclorb',1,1)
!     call abend('after  clorb                                      ')
#endif
    do l=1,2
      clo0(l)=clo(l)
      clop0(l)=clop(l)
    end do
    call clorb(zero)
#ifdef DEBUG
!     call dumpbin('aclorb',1,1)
!     call abend('after  clorb                                      ')
#endif
    do l=1,2
      ll=2*l
      di0(l)=(clo0(l)-clo(l))/ded
      dip0(l)=(clop0(l)-clop(l))/ded
    end do
    call corrorb

    if(irmod2.eq.1) call rmod(dp1)
    if(iqmod.ne.0) call qmod0
    if(ichrom.eq.1.or.ichrom.eq.3) call chroma
    if(iskew.ne.0) call decoup
    if(ilin.eq.1.or.ilin.eq.3) then
      call linopt(dp1)
    end if
#ifdef DEBUG
!     call dumpbin('bbb',96,996)
!     call abend('bbb                                               ')
#endif
    ! beam-beam element
    nlino = nlin
    nlin  = 0
    if(nbeam.ge.1) then
      do i=1,nele
        if(kz(i).eq.20) then
          nlin=nlin+1
          if(nlin.gt.nele) call prror(81)
          bezl(nlin)=bez(i)
        end if
      end do
    end if
    if(isub == 1) call subre(dp1)
    if(ise  == 1) call search(dp1)
#ifdef DEBUG
!     call dumpbin('asearch',95,995)
!     call abend('asearch                                           ')
#endif
    !! Initialize kicks
    izu=0
    do i=1,iu
#ifdef DEBUG
!       call warr('i/iu',0d0,i,iu,0,0)
!       write(*,*) 'do 150 i/iu',i,iu
#endif
      ix=ic(i)
      if(ix.le.nblo) cycle
      ix=ix-nblo
      kpz=kp(ix)
      kzz=kz(ix)
      if(kpz.eq.6.or.kzz.eq.0.or.kzz.eq.20.or.kzz.eq.22) cycle
      if(kzz.eq.15) cycle
      if(iorg.lt.0) mzu(i)=izu
      izu=mzu(i)+1
      smizf(i)=zfz(izu)*ek(ix)
      smiv(m,i)=sm(ix)+smizf(i) ! Also in initalize_element!
      smi(i)=smiv(m,i)          ! Also in initalize_element!
#ifdef DEBUG
!         call warr('smizf(i)',smizf(i),i,0,0,0)
!         call warr('smiv(m,i)',smiv(m,i),m,i,0,0)
!         call warr('smi(i)',smi(i),i,0,0,0)
#endif
      izu=izu+1
      xsiv(m,i)=xpl(ix)+zfz(izu)*xrms(ix)
      xsi(i)=xsiv(m,i)
      izu=izu+1
      zsiv(m,i)=zpl(ix)+zfz(izu)*zrms(ix)
      zsi(i)=zsiv(m,i)
      if(mout2.eq.1) then
        if(kzz.eq.11) zfz(izu-2)=zero
        if(abs(ek(ix)).le.pieni) zfz(izu-2)=zero
        if(abs(xrms(ix)).le.pieni) zfz(izu-1)=zero
        if(abs(zrms(ix)).le.pieni) zfz(izu)=zero
        if(icextal(i) > 0) then
          write(31,"(a48,1p,d19.11,2d14.6,d17.9)") bez(ix),zfz(izu-2),zfz(izu-1),zfz(izu),fluc_errAlign(3,icextal(i))
        else if(icextal(i) < 0) then
          write(31,"(a48,1p,d19.11,2d14.6,d17.9)") bez(ix),zfz(izu-2),zfz(izu-1),zfz(izu),fluc_errZFZ(4,-icextal(i))
        else
          write(31,"(a48,1p,d19.11,2d14.6,d17.9)") bez(ix),zfz(izu-2),zfz(izu-1),zfz(izu),zero
        end if
      endif

!-- MULTIPOLE BLOCK
      if(kzz.eq.11) then
        dynk_izuIndex(ix)=izu
!-- Initialize multipoles, combining settings from fort.2 with
!-- coefficients from MULT and random values from FLUC.
!-- Used in program maincr and from initialize_element.
        r0=ek(ix)
        if(abs(r0).le.pieni) cycle
        nmz=nmu(ix)
        if(nmz.eq.0) then
          izu=izu+2*mmul
          cycle
        end if
        im=irm(ix)
        r0a=one
        do k=1,nmz
          izu=izu+1
          aaiv(k,m,i)=(ed(ix)*(ak0(im,k)+zfz(izu)*aka(im,k)))/r0a !hr05
          aai(i,k)=aaiv(k,m,i)
          izu=izu+1
          bbiv(k,m,i)=(ed(ix)*(bk0(im,k)+zfz(izu)*bka(im,k)))/r0a !hr05
          bbi(i,k)=bbiv(k,m,i)
          r0a=r0a*r0
        end do
        izu=izu+2*mmul-2*nmz
      end if
    end do
#ifdef DEBUG
!     call dumpbin('ado 150',150,150)
!     call abend('ado 150                                           ')
#endif
        dp1=zero
        if(ichrom.gt.1) then
          itiono=ition
          ition=0
          call chromda
          ition=itiono
          do ncrr=1,iu
            ix=ic(ncrr)
            if(ix.gt.nblo) ix=ix-nblo
            if(ix.eq.is(1).or.iratioe(ix).eq.is(1)) then
              smiv(m,ncrr)=smi(ncrr)
            else if(ix.eq.is(2).or.iratioe(ix).eq.is(2)) then
              smiv(m,ncrr)=smi(ncrr)
            endif
          enddo
        endif
        dp1=dp00
        dp0=dp00
        if(imc.gt.1) then
          ddp1=(two*dp0)/(real(imc,fPrec)-one)                                 !hr05
        endif
        do 250 ib=1,imc
          if(imc.gt.1) then
            dp1=dp0-(real(ib,fPrec)-one)*ddp1                                  !hr05
          endif
          dp10=dp1
!-----------------------------------------------------------------------
          if(idp /= 1 .or. iation /= 1) iclo6=0
          if(iclo6 == 1 .or. iclo6 == 2) then
            if(ib == 1) then
              if(iclo6r == 0) then
                clo6(1)  = clo(1)
                clop6(1) = clop(1)
                clo6(2)  = clo(2)
                clop6(2) = clop(2)
                clo6(3)  = zero
                clop6(3) = zero
              else
                write(lout,"(a)") "MAINCR> Reading closed orbit guess from fort.33"
                call readFort33
              end if
              call clorb(zero)
              call betalf(zero,qw)
              call phasad(zero,qwc)
              sigm(1) = clo6(3)
              dps(1)  = clop6(3)
              call qmodda(3,qwc)
              if(ilin >= 2) then
                nlinoo = nlin
                nlin   = nlino
                ilinc  = 1
                call mydaini(2,2,6,3,6,1)
                nlin   = nlinoo
              end if
              dp1 = dp10+clop6(3)
            end if
            if(iqmod6 == 1) then
              do ncrr=1,iu
                ix=ic(ncrr)
                if(ix.gt.nblo) ix=ix-nblo
                if(ix.eq.iq(1).or.iratioe(ix).eq.iq(1)) then
                  smiv(m,ncrr)=smi(ncrr)
                else if(ix.eq.iq(2).or.iratioe(ix).eq.iq(2)) then
                  smiv(m,ncrr)=smi(ncrr)
                endif
              enddo
            endif

            do 190 ib1=1,napx
              ib3=ib1+(m+ib-2)*napx
!--beam-beam element
              clo6v(1,ib3)=clo6(1)
              clo6v(2,ib3)=clo6(2)
              clo6v(3,ib3)=clo6(3)
              clop6v(1,ib3)=clop6(1)
              clop6v(2,ib3)=clop6(2)
              clop6v(3,ib3)=clop6(3)
              di0xs(ib3)=di0(1)
              di0zs(ib3)=di0(2)
              dip0xs(ib3)=dip0(1)
              dip0zs(ib3)=dip0(2)
              qwcs(ib3,1)=qwc(1)
              qwcs(ib3,2)=qwc(2)
              qwcs(ib3,3)=qwc(3)

              do i2=1,6
                do j2=1,6
                  tas(ib3,i2,j2)=tasm(i2,j2)
                end do
              end do

  190       continue
          else
            if(idp.eq.1.and.iation.eq.1) then
              ncorruo=ncorru
              ncorru=1
              call clorb(zero)
#ifdef DEBUG
!     call dumpbin('aclorb',1,1)
!     call abend('after  clorb                                      ')
#endif
              call betalf(zero,qw)
              call phasad(zero,qwc)
#ifdef DEBUG
!     call dumpbin('abetphas',1,1)
!     call abend('after  abetphas                                   ')
#endif
!--beam-beam element
              if(nbeam.ge.1) then
              nd=3
              nd2=6
#include "include/beamcou.f90"
              endif
              ncorru=ncorruo
              iqmodc=3
              call mydaini(2,2,6,3,6,1)
#ifdef DEBUG
!     call dumpbin('bmydaini',999,9999)
!     call abend('before mydaini                                    ')
#endif
              do i=1,2
                qwc(i)=real(int(qwc(i)),fPrec)+wxys(i)
              enddo
              if(ilin.ge.2) then
#ifdef DEBUG
!     call dumpbin('bmydaini',999,9999)
!     call abend('before mydaini                                    ')
#endif
                nlinoo=nlin
                nlin=nlino
                ilinc=1
                call mydaini(2,2,6,3,6,1)
#ifdef DEBUG
!     call dumpbin('amydaini',999,9999)
!     call abend('after  mydaini                                    ')
#endif
                nlin=nlinoo
              endif
            else
              dps(1)=dp1
              ncorruo=ncorru
              ncorru=1
              call clorb(dp1)
              call betalf(dp1,qw)
              call phasad(dp1,qwc)
#ifdef DEBUG
!     call dumpbin('abetphas',1,1)
!     call abend('after  abetphas                                   ')
#endif
              dp1=zero
!--beam-beam element
#ifdef DEBUG
!     call dumpbin('bbeam',1,1)
!     call abend('after bbeam                                       ')
!     write(*,*) 'call qmodda at beam-beam'
#endif
              dp1=dps(1)
              ncorru=ncorruo
              if(nvar2.le.5) then
                itiono=ition
                ition=0
              endif
              call qmodda(2,qwc)
#ifdef DEBUG
!     call dumpbin('aqmodda',3,2)
!     call abend('after  qmodda 3 2                                 ')
#endif
              if(nvar2.le.5) ition=itiono
              if(nvar2.le.4.and.ithick.eq.1) call envar(dp1)

              if(ilin.ge.2) then
                nlinoo=nlin
                nlin=nlino
                iqmodc=2
                call mydaini(1,2,5,2,5,1)
                ilinc=1
                call mydaini(2,2,5,2,5,1)
                nlin=nlinoo
              endif

              do ncrr=1,iu
                ix=ic(ncrr)
                if(ix.gt.nblo) ix=ix-nblo
                if(ix.eq.iq(1).or.iratioe(ix).eq.iq(1)) then
                  smiv(m,ncrr)=smi(ncrr)
                else if(ix.eq.iq(2).or.iratioe(ix).eq.iq(2)) then
                  smiv(m,ncrr)=smi(ncrr)
                endif
              enddo
            endif

            do 170 i=1,napx
              iar=(m+ib-2)*napx+i
              clo6v(1,iar)=clo(1)
              clop6v(1,iar)=clop(1)
              clo6v(2,iar)=clo(2)
              clop6v(2,iar)=clop(2)
              di0xs(iar)=di0(1)
              di0zs(iar)=di0(2)
              dip0xs(iar)=dip0(1)
              dip0zs(iar)=dip0(2)
              qwcs(iar,1)=qwc(1)
              qwcs(iar,2)=qwc(2)
              qwcs(iar,3)=zero

              do i2=1,4
                do j2=1,4
                  tas(iar,i2,j2)=tasm(i2,j2)
                end do
              end do

  170       continue
          endif
          iar=(m+ib-2)*napx+1

! save tas matrix and closed orbit for later dumping of the beam
! distribution at the first element (i=-1)
! dumptas(*,*) [mm,mrad,mm,mrad,1] canonical variables
! tas(iar,*,*) [mm,mrad,mm,mrad,1] canonical variables
! clo6v,clop6v [mm,mrad,mm,mrad,1] canonical variables (x' or px?)
! for the initialization of the particles. Only in 5D thick the ta
! matrix is different for each particle.
! -> implement a check for this!
! In 4d,6d thin+thick and 5d thin we have:
!   tas(ia,*,*) = tas(1,*,*) for all particles ia
          if (iar .eq. 1) then
             do i3=1,3
                dumpclo(-1,i3*2-1) = clo6v(i3,1)
                dumpclo(-1,i3*2)   = clop6v(i3,1)
             enddo
             dumptas(-1,:,:) = tas(1,:,:)
!     invert the tas matrix
             call invert_tas(dumptasinv(-1,:,:),dumptas(-1,:,:))
!     dumptas and dumptasinv are now in units [mm,mrad,mm,mrad,1]
          endif
!     tas(iar,*,*) [mm,mrad,mm,mrad,1]

! convert to [mm,mrad,mm,mrad,1.e-3] for optics calculation
          tasiar16=tas(iar,1,6)*c1m3
          tasiar26=tas(iar,2,6)*c1m3
          tasiar36=tas(iar,3,6)*c1m3
          tasiar46=tas(iar,4,6)*c1m3
          tasiar56=tas(iar,5,6)*c1m3
          tasiar61=tas(iar,6,1)*c1e3
          tasiar62=tas(iar,6,2)*c1e3
          tasiar63=tas(iar,6,3)*c1e3
          tasiar64=tas(iar,6,4)*c1e3
          tasiar65=tas(iar,6,5)*c1e3
          bet0(1)=tas(iar,1,1)**2+tas(iar,1,2)**2                        !hr05
          bet0x2 =tas(iar,1,3)**2+tas(iar,1,4)**2                        !hr05
          bet0x3 =tas(iar,1,5)**2+tasiar16**2                            !hr05
          gam0x1 =tas(iar,2,1)**2+tas(iar,2,2)**2                        !hr05
          gam0x2 =tas(iar,2,3)**2+tas(iar,2,4)**2                        !hr05
          gam0x3 =tas(iar,2,5)**2+tasiar26**2                            !hr05
      alf0(1)=-one*(tas(iar,1,1)*tas(iar,2,1)+tas(iar,1,2)*tas(iar,2,2)) !hr05
      alf0x2 =-one*(tas(iar,1,3)*tas(iar,2,3)+tas(iar,1,4)*tas(iar,2,4)) !hr05
      alf0x3 =-one*(tas(iar,1,5)*tas(iar,2,5)+tasiar16*tasiar26)         !hr05
          bet0(2)=tas(iar,3,3)**2+tas(iar,3,4)**2                        !hr05
          bet0z2 =tas(iar,3,1)**2+tas(iar,3,2)**2                        !hr05
          bet0z3 =tas(iar,3,5)**2+tasiar36**2                            !hr05
          gam0z1 =tas(iar,4,3)**2+tas(iar,4,4)**2                        !hr05
          gam0z2 =tas(iar,4,1)**2+tas(iar,4,2)**2                        !hr05
          gam0z3 =tas(iar,4,5)**2+tasiar46**2                            !hr05
      alf0(2)=-one*(tas(iar,3,3)*tas(iar,4,3)+tas(iar,3,4)*tas(iar,4,4)) !hr05
      alf0z2 =-one*(tas(iar,3,1)*tas(iar,4,1)+tas(iar,3,2)*tas(iar,4,2)) !hr05
      alf0z3 =-one*(tas(iar,3,5)*tas(iar,4,5)+tasiar36*tasiar46)         !hr05
          bet0s1 =tas(iar,5,5)**2+tasiar56**2                            !hr05
          bet0s2 =tas(iar,5,1)**2+tas(iar,5,2)**2                        !hr05
          bet0s3 =tas(iar,5,3)**2+tas(iar,5,4)**2                        !hr05
          gam0s1 =tasiar65**2+tas(iar,6,6)**2                            !hr05
          gam0s2 =tasiar61**2+tasiar62**2                                !hr05
          gam0s3 =tasiar63**2+tasiar64**2                                !hr05
          alf0s1 =-one*(tas(iar,5,5)*tasiar65+tasiar56*tas(iar,6,6))     !hr05
          alf0s2 =-one*(tas(iar,5,1)*tasiar61+tas(iar,5,2)*tasiar62)     !hr05
          alf0s3 =-one*(tas(iar,5,3)*tasiar63+tas(iar,5,4)*tasiar64)     !hr05
#ifdef DEBUG
!     call dumpbin('abib1',1,1)
!     call abend('after bib1                                        ')
#endif
          do 220 ib1=1,napx
            iar=ib1+(m+ib-2)*napx

            do ib2=1,6
              do ib3=1,6
                tau(ib2,ib3)=tas(iar,ib3,ib2)
              end do
            end do

            if(abs(tau(1,1)).le.pieni.and.abs(tau(2,2)).le.pieni) then
              tau(1,1)=one
              tau(2,2)=one
            endif
            if(abs(tau(3,3)).le.pieni.and.abs(tau(4,4)).le.pieni) then
              tau(3,3)=one
              tau(4,4)=one
            endif
            if(abs(tau(5,5)).le.pieni.and.abs(tau(6,6)).le.pieni) then
              tau(5,5)=one
              tau(6,6)=one
              call dinv(6,tau,6,idummy,nerror)
              its6d=0
              if(ntwin.ne.2) then
                taus=(((((((((((((((((((                                &!hr05
     &abs(tau(5,1))+abs(tau(5,2)))+abs(tau(5,3)))+abs                   &!hr05
     &(tau(5,4)))+abs(tau(5,5)))+abs(tau(5,6)))+abs(tau(6,1)))          &!hr05
     &+abs(tau(6,2)))+abs(tau(6,3)))+abs(tau(6,4)))+abs                 &!hr05
     &(tau(6,5)))+abs(tau(6,6)))+abs(tau(1,5)))+abs(tau(2,5)))          &!hr05
     &+abs(tau(3,5)))+abs(tau(4,5)))+abs(tau(1,6)))+abs                 &!hr05
     &(tau(2,6)))+abs(tau(3,6)))+abs(tau(4,6)))-two                      !hr05
                if(abs(taus).ge.pieni) its6d=1
              endif
              do ib2=1,6
                do ib3=1,6
                  tasau(iar,ib2,ib3)=tau(ib2,ib3)
                end do
              end do
            endif
  220     continue
          if(ierro.ne.0) then
            write(lout,10230) dp1
            goto 520
          endif
          write(lout,10070)
          phag=(phas*c180e0)/pi                                           !hr05
          if((idp.eq.0).or.(abs(phas).le.pieni.and.ition.eq.0))         &
     &write(lout,10170)                                                 &
     &qwc(1),clo(1),clop(1),                                            &
     &bet0(1),alf0(1),gam0x1,bet0x2,alf0x2,gam0x2,                      &
     &qwc(2),clo(2),clop(2),                                            &
     &bet0(2),alf0(2),gam0z1,bet0z2,alf0z2,gam0z2
          if(idp.eq.1.and.iation.eq.1.and.abs(phas).gt.pieni) then
            if(iclo6.eq.0) then
              write(lout,10150) phag,                                   &
     &qwc(1),clo(1),clop(1),                                            &
     &bet0(1),alf0(1),gam0x1,bet0x2,alf0x2,gam0x2,                      &
     &qwc(2),clo(2),clop(2),                                            &
     &bet0(2),alf0(2),gam0z1,bet0z2,alf0z2,gam0z2
            else
              write(lout,10160) phag,                                   &
     &qwc(1),clo6(1),clop6(1),                                          &
     &bet0(1),alf0(1),gam0x1,bet0x2,alf0x2,gam0x2,                      &
     &bet0x3,alf0x3,gam0x3,                                             &
     &qwc(2),clo6(2),clop6(2),                                          &
     &bet0(2),alf0(2),gam0z1,bet0z2,alf0z2,gam0z2,                      &
     &bet0z3,alf0z3,gam0z3,                                             &
     &qwc(3),clo6(3),clop6(3),                                          &
     &bet0s1,alf0s1,gam0s1,bet0s2,alf0s2,gam0s2,                        &
     &bet0s3,alf0s3,gam0s3
            endif
          endif
          if(idp.eq.1.and.ition.eq.0.and.abs(phas).gt.pieni)            &
     &write(lout,10190) phag,                                           &
     &qwc(1),clo(1),clop(1),                                            &
     &bet0(1),alf0(1),gam0x1,bet0x2,alf0x2,gam0x2,                      &
     &qwc(2),clo(2),clop(2),                                            &
     &bet0(2),alf0(2),gam0z1,bet0z2,alf0z2,gam0z2
          if(idp.eq.1.and.abs(phas).le.pieni.and.iation.eq.1) then
            if(iclo6.eq.0) then
              write(lout,10210)                                         &
     &qwc(1),clo(1),clop(1),                                            &
     &bet0(1),alf0(1),gam0x1,bet0x2,alf0x2,gam0x2,                      &
     &qwc(2),clo(2),clop(2),                                            &
     &bet0(2),alf0(2),gam0z1,bet0z2,alf0z2,gam0z2
            else
              write(lout,10220)                                         &
     &qwc(1),clo6(1),clop6(1),                                          &
     &bet0(1),alf0(1),gam0x1,bet0x2,alf0x2,gam0x2,                      &
     &bet0x3,alf0x3,gam0x3,                                             &
     &qwc(2),clo6(2),clop6(2),                                          &
     &bet0(2),alf0(2),gam0z1,bet0z2,alf0z2,gam0z2,                      &
     &bet0z3,alf0z3,gam0z3,                                             &
     &qwc(3),clo6(3),clop6(3),                                          &
     &bet0s1,alf0s1,gam0s1,bet0s2,alf0s2,gam0s2,                        &
     &bet0s3,alf0s3,gam0s3
            endif
          endif
          write(lout,10080) dp1
          e0f=sqrt(e0**2-nucm0**2)                                         !hr05
          if(iclo6.eq.0) then
            write(lout,10110) clo(1),clop(1),clo(2),clop(2),idz(1),     &
     &idz(2),                                                           &
     &iver, idfor,iclo6,ition
          else
            write(lout,10120) clo6(1),clop6(1),clo6(2),clop6(2),clo6(3),&
     &clop6(3), idz(1),idz(2),iver,idfor,iclo6,ition
          endif

          do ib1=1,napx
            !Loop over all particles (not pairs) with the same
            ! ib (momentum variation, 1..imc ) and
            !  m (seed,               1..mmac).
            !It appears that only the odd (1,3,5,..) indices are actually used?
            ib2=ib0+ib1        ! ib0 is fixed to 0 => ib2 equals ib1
            clov(1,ib2)=clo(1)
            clov(2,ib2)=clo(2)
            clopv(1,ib2)=clop(1)
            clopv(2,ib2)=clop(2)
            bet0v(ib2,1)=bet0(1)
            bet0v(ib2,2)=bet0(2)
            alf0v(ib2,1)=alf0(1)
            alf0v(ib2,2)=alf0(2)
            ampv(ib2)=amp(1)-damp*real(ib1-1,fPrec) !hr05

            if(ib1.eq.napx-1 .and. ib1.ne.1) then
              !Make sure that last amplitude EXACTLY corresponds to the end amplitude amp0
              ! This is helpfull when doing DA studies and checking the "overlap"
              ampv(ib2)=amp0
            end if

            dp0v(ib2)=dp10
            dpsv(ib2)=dp10
            oidpsv(ib2)=one/(one+dp1)
! Heavy ion variable
            moidpsv(ib2)=mtc(ib2)/(one+dp1)
            nms(ib2)=m

            if(ithick.eq.1) then
              do i=1,nele
                ekv(ib2,i)=ek(i)
              end do
            end if

          end do

          ib0=ib0+napx
  250   continue
  260 continue
#ifdef DEBUG
!     call dumpbin('ado 260',260,260)
!     call abend('ado 260                                           ')
#endif

      napx=(napx*imc)*mmac                                               !hr05

#ifdef FLUKA

!     A.Mereghetti, P. Garcia Ortega, D.Sinuela Pastor, V. Vlachoudis
!             for the FLUKA Team
!     last modified: 11-06-2014
!     start connection to FLUKA and initialise max ID
!     inserted in main code by the 'fluka' compilation flag
      if(fluka_enable) then
        fluka_con = fluka_is_running()
        if(fluka_con.eq.-1) then
          write(lout,*) '[Fluka] Error: Fluka is expected to run but it is'
           write(lout,*) '               NOT actually the case'
          write(fluka_log_unit,*) '# Fluka is expected to run but it is'
          write(fluka_log_unit,*) '               NOT actually the case'
          call prror(-1)
        endif
        write(lout,*) '[Fluka] Initializing FlukaIO interface...'
        write(fluka_log_unit,*) '# Initializing FlukaIO interface...'
        fluka_con = fluka_connect()
        if(fluka_con.eq.-1) then
          write(lout,*) '[Fluka] Error connecting to Fluka server'
          write(fluka_log_unit,*) '# Error connecting to Fluka server'
          call prror(-1)
        endif
        write(lout,*) '[Fluka] Successfully connected to Fluka server'
        write(fluka_log_unit,*) '# Successfully connected to Fluka server'
        fluka_connected = .true.
      endif

#endif

#ifdef CR
      write(93,*) 'MAINCR setting napxo=',napx
      endfile (93,iostat=ierro)
      backspace (93,iostat=ierro)
#endif
      napxo=napx
      if(ibidu.eq.1) then
        ! Note: Keep in sync with read(32) below
        write(32) &
        ierro,erbez,pi2,pisqrt,rad,il,mper,mblo,mbloz,msym,kanf,iu,ic,    &
        ed,el,ek,sm,kz,kp,xpl,xrms,zpl,zrms,mel,mtyp,mstr,a,bl1,bl2,rvf,  &
        idfor,napx,napxo,numlr,nde,nwr,ird,imc,irew,ntwin,iclo6,iclo6r,   &
        iver,ibidu,qs,e0,pma,ej,ejf,phas0,phas,hsy,crad,                  &
        hsyc,phasc,dppoff,sigmoff,tlen,                                   &
        iicav,itionc,ition,idp,ncy,ixcav,dpscor,                          &
        sigcor,icode,idam,its6d,bk0,ak0,bka,aka,benki,benkc,r00,irm,nmu,  &
        zfz,iorg,mzu,bezr,izu0,mmac,mcut,tiltc,tilts,                     &
        mout2,icext,icextal,aper,di0,dip0,ta,dma,dmap,dkq,dqq,de0,ded,dsi,&
        dech,dsm0,itco,itcro,itqv,qw0,iq,iqmod,kpa,iqmod6,bez,            &
        elbe,bezb,ilin,nt,iprint,ntco,eui,euii,nlin,bezl,betam,pam,betac, &
        pac,bclorb,nhmoni,nhcorr,nvmoni,nvcorr,ncororb,sigma0,iclo,       &
        ncorru,ncorrep,icomb0,icomb,ratio,ratioe,iratioe,                 &
        icoe,ise,mesa,mp,m21,m22,m23,                                     &
        ise1,ise2,ise3,isea,qxt,qzt,tam1,tam2,isub,nta,nte,ipt,totl,rtc,  &
        rts,ire,ipr,irmod2,dtr,nre,nur,nch,nqc,npp,nrr,nu,dphix,dphiz,qx0,&
        qz0,dres,dfft,cma1,cma2,nstart,nstop,iskip,iconv,imad,ipos,iav,   &
        iwg,ivox,ivoz,ires,ifh,toptit,kwtype,itf,icr,idis,icow,istw,iffw, &
        nprint,ndafi,qwsk,betx,betz,pttemp,temptr,kxxa,                   &
        alfx,alfz,iskew,nskew,hmal,sixtit,commen,ithick,clo6,clop6,dki,   &
        sigman,sigman2,sigmanq,clobeam,beamoff,parbe,track6d,ptnfac,      &
        sigz,sige,partnum,parbe14,emitx,emity,emitz,gammar,nbeam,ibbc,    &
        ibeco,ibtyp,lhc,cotr,rrtr,imtr,bbcu,ibb6d,imbb,wire_num,          &
        as,al,sigm,dps,idz,dp1,itra,                                      &
        x,y,bet0,alf0,clo,clop,cro,is,ichrom,nnumxv,xsi,zsi,smi,aai,      &
        bbi,ampt,tlim,tasm,preda,idial,nord,nvar,                         &
        nvar2,nsix,ncor,ipar,nordf,                                       &
        nvarf,nord1,ndimf,idptr,inorm,imod1,imod2,                        &
        ekv,fokqv,aaiv,bbiv,smiv,zsiv,xsiv,xsv,zsv,qw,qwc,clo0,           &
        clop0,eps,epsa,ekk,cr,ci,xv,yv,dam,ekkv,sigmv,dpsv,dp0v,sigmv6,   &
        dpsv6,ejv,ejfv,xlv,zlv,pstop,rvv,                                 &
        ejf0v,numxv,nms,nlostp,dpd,                                       &
        dpsq,fok,rho,fok1,si,co,g,gl,sm1,sm2,sm3,sm12,as3,as4,as6,sm23,   &
        rhoc,siq,aek,afok,hp,hm,hc,hs,wf,wfa,wfhi,rhoi,hi,fi,hi1,xvl,yvl, &
        ejvl,dpsvl,oidpsv,sigmvl,iv,aperv,ixv,clov,clopv,alf0v,bet0v,ampv,&
        clo6v,clop6v,hv,bl1v,tas,qwcs,di0xs,di0zs,dip0xs,dip0zs,xau,cloau,&
        di0au,tau,tasau,wx,x1,x2,fake,e0f,numx,cotr,rrtr,imtr

      endif
  550 continue


      if (idp.eq.0.or.ition.eq.0) then
         !4D tracking
         if (iclo6 .ne. 0) then
            write(lout,*) "ERROR: Doing 4D tracking but iclo6=",iclo6
            write(lout,*) "Expected iclo6.eq.0. for 4D tracking."
            call prror(-1)
         endif
      else
         !6D tracking
         if (iclo6 .eq. 0) then
            write(lout,*) "ERROR: Doing 6D tracking but iclo6=",iclo6
            write(lout,*) "Expected iclo6.ne.0. for 6D tracking."
            call prror(-1)
         endif
      endif


!!!   GENERATE THE INITIAL DISTRIBUTION
      if(ibidu.eq.2) then
        ! Note: Keep in sync with write(32) above
        read(32) &
        ierro,erbez,pi2,pisqrt,rad,il,mper,mblo,mbloz,msym,kanf,iu,ic,    &
        ed,el,ek,sm,kz,kp,xpl,xrms,zpl,zrms,mel,mtyp,mstr,a,bl1,bl2,rvf,  &
        idfor,napx,napxo,numlr,nde,nwr,ird,imc,irew,ntwin,iclo6,iclo6r,   &
        iver,ibidu,qs,e0,pma,ej,ejf,phas0,phas,hsy,crad,                  &
        hsyc,phasc,dppoff,sigmoff,tlen,                                   &
        iicav,itionc,ition,idp,ncy,ixcav,dpscor,                          &
        sigcor,icode,idam,its6d,bk0,ak0,bka,aka,benki,benkc,r00,irm,nmu,  &
        zfz,iorg,mzu,bezr,izu0,mmac,mcut,tiltc,tilts,                     &
        mout2,icext,icextal,aper,di0,dip0,ta,dma,dmap,dkq,dqq,de0,ded,dsi,&
        dech,dsm0,itco,itcro,itqv,qw0,iq,iqmod,kpa,iqmod6,bez,            &
        elbe,bezb,ilin,nt,iprint,ntco,eui,euii,nlin,bezl,betam,pam,betac, &
        pac,bclorb,nhmoni,nhcorr,nvmoni,nvcorr,ncororb,sigma0,iclo,       &
        ncorru,ncorrep,icomb0,icomb,ratio,ratioe,iratioe,                 &
        icoe,ise,mesa,mp,m21,m22,m23,                                     &
        ise1,ise2,ise3,isea,qxt,qzt,tam1,tam2,isub,nta,nte,ipt,totl,rtc,  &
        rts,ire,ipr,irmod2,dtr,nre,nur,nch,nqc,npp,nrr,nu,dphix,dphiz,qx0,&
        qz0,dres,dfft,cma1,cma2,nstart,nstop,iskip,iconv,imad,ipos,iav,   &
        iwg,ivox,ivoz,ires,ifh,toptit,kwtype,itf,icr,idis,icow,istw,iffw, &
        nprint,ndafi,qwsk,betx,betz,pttemp,temptr,kxxa,                   &
        alfx,alfz,iskew,nskew,hmal,sixtit,commen,ithick,clo6,clop6,dki,   &
        sigman,sigman2,sigmanq,clobeam,beamoff,parbe,track6d,ptnfac,      &
        sigz,sige,partnum,parbe14,emitx,emity,emitz,gammar,nbeam,ibbc,    &
        ibeco,ibtyp,lhc,cotr,rrtr,imtr,bbcu,ibb6d,imbb,wire_num,          &
        as,al,sigm,dps,idz,dp1,itra,                                      &
        x,y,bet0,alf0,clo,clop,cro,is,ichrom,nnumxv,xsi,zsi,smi,aai,      &
        bbi,ampt,tlim,tasm,preda,idial,nord,nvar,                         &
        nvar2,nsix,ncor,ipar,nordf,                                       &
        nvarf,nord1,ndimf,idptr,inorm,imod1,imod2,                        &
        ekv,fokqv,aaiv,bbiv,smiv,zsiv,xsiv,xsv,zsv,qw,qwc,clo0,           &
        clop0,eps,epsa,ekk,cr,ci,xv,yv,dam,ekkv,sigmv,dpsv,dp0v,sigmv6,   &
        dpsv6,ejv,ejfv,xlv,zlv,pstop,rvv,                                 &
        ejf0v,numxv,nms,nlostp,dpd,                                       &
        dpsq,fok,rho,fok1,si,co,g,gl,sm1,sm2,sm3,sm12,as3,as4,as6,sm23,   &
        rhoc,siq,aek,afok,hp,hm,hc,hs,wf,wfa,wfhi,rhoi,hi,fi,hi1,xvl,yvl, &
        ejvl,dpsvl,oidpsv,sigmvl,iv,aperv,ixv,clov,clopv,alf0v,bet0v,ampv,&
        clo6v,clop6v,hv,bl1v,tas,qwcs,di0xs,di0zs,dip0xs,dip0zs,xau,cloau,&
        di0au,tau,tasau,wx,x1,x2,fake,e0f,numx,cotr,rrtr,imtr

        damp=((amp(1)-amp0)/real(napx/2-1,fPrec))/two                          !hr05
      endif

      do i=1,npart
        pstop(i)=.false.
        nnumxv(i)=numl
        numxv(i)=numl
      end do

      rat0=rat


!----- Initial distribution creation

!     A.Mereghetti, for the FLUKA Team
!     last modified: 14-06-2014
!     acquisition of initial distribution moved out of loop
!     always in main code

      if ( idfor.eq.3 ) then
!       A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!       last modified: 17-07-2013
!       initialize particle distribution, read from file
!       always in main code

        if(.not. dist_enable) then
          write(lout,"(a)") "MAINCR> ERROR idfor set to 3 but DIST block not present."
          call prror(-1)
        endif


        e0f=sqrt(e0**2-nucm0**2)       ! hisix

        call dist_readdis( napx, npart, e0, e0f, clight, xv(1,:), xv(2,:), yv(1,:), yv(2,:), sigmv(:), ejfv(:) &
& ,naa(:), nzz(:), nucm(:) )      ! hisix

!       finalise beam distribution creation
        do j=1, napx
!         values related to losses
          nlostp(j) = j
          pstop (j) = .false.

!         values related to momentum
!         old proton only terms:
!          ejv   (j) = sqrt(ejfv(j)**2+pma**2)
!          dpsv  (j) = (ejfv(j)-e0f)/e0f
!          oidpsv(j) = one/(one+dpsv(j))

          ejv   (j)   = sqrt(ejfv(j)**2+nucm(j)**2)              ! hiSix
          dpsv  (j)   = (ejfv(j)*(nucm0/nucm(j))-e0f)/e0f         ! hiSix
          oidpsv(j)   = one/(one+dpsv(j))
          mtc     (j) = (nzz(j)*nucm0)/(zz0*nucm(j))
          moidpsv (j) = mtc(j)*oidpsv(j)
          omoidpsv(j) = c1e3*((one-mtc(j))*oidpsv(j))

!         check existence of on-momentum particles in the distribution
          if ( abs(dpsv(j)).lt.c1m15 .or.  abs( (ejv(j)-e0)/e0 ) .lt.c1m15 ) then

!           warning with old infos:
            write(lout,*)''
            write(lout,'(5X,A22)') 'on-momentum particle!!'
            write(lout,'(5X,10X,4(1X,A25))') "momentum [MeV/c]","total energy [MeV]","Dp/p","1/(1+Dp/p)"
            write(lout,'(5X,"ORIGINAL: ",4(1X,1PE25.18))') ejfv(j), ejv(j), dpsv(j), oidpsv(j)

!            ejfv(j)   = e0f
!            ejv(j)    = e0
            ejfv(j)   = e0f*(nucm(j)/nucm0)          ! P. HERMES for hiSix
            ejv(j)    = sqrt(ejfv(j)**2+nucm(j)**2)  ! P. HERMES for hiSix
            dpsv(j)   = zero
            oidpsv(j) = one

!           warning with new infos:
            write(lout,'(5X,"CORRECTED:",4(1X,1PE25.18))') ejfv(j), ejv(j), dpsv(j), oidpsv(j)
            write(lout,*)''
          endif
        end do

! hisix
        write(lout,*) 'Heavy-Ion SixTrack'
        write(lout,*) '------------------'
        write(lout,*) 'Reference ion species: [A,Z,M]', aa0, zz0, nucm0
        write(lout,*) 'Reference energy [Z TeV]: ', c1m6*e0/zz0

! hisix - debugging
!        write(lout,*) 'Properties of tracked ion bunch [A,Z,E(MeV)], etc'
!        do j=1,napx
!          write(lout,*) naa(j),nzz(j),e0f*(nucm(j)/nucm0), ejfv(j), mtc(j), dpsv(j), ejv(j)
!        end do

!       A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!       last modified: 07-02-2014
!       in principle there is no need to fill in the unused places:
!       - nlostp(j) = j        with j=1,npart    filled in trauthin/trauthck
!       - pstop (j) = .false.  with j=1,npart    filled in maincr
!       - ejv   (j) = zero     with j=1,npart    filled in maincr
!       - dpsv  (j) = zero     with j=1,npart    filled in maincr
!       - oidpsv(j) = one      with j=1,npart    filled in maincr
!       nevertheless, let's do it, to be fully sure:
        do j=napx+1,npart
!         values related to losses
          nlostp(j) = j
          pstop (j) = .true.
!         values related to momentum
          ejv   (j) = zero
          dpsv  (j) = zero
          oidpsv(j) = one

          mtc   (j) = one         ! P. HERMES for hiSix
          naa   (j) = aa0
          nzz   (j) = zz0
          nucm  (j) = nucm0
          moidpsv (j) = one
          omoidpsv(j) = zero      ! P. HERMES for hiSix
        enddo

!       add closed orbit
        if(iclo6.eq.2) then
          do j=1, napx
            xv(1,j)=xv(1,j)+clo6v(1,j)
            yv(1,j)=yv(1,j)+clop6v(1,j)
            xv(2,j)=xv(2,j)+clo6v(2,j)
            yv(2,j)=yv(2,j)+clop6v(2,j)
            sigmv(j)=sigmv(j)+clo6v(3,j)
            dpsv(j)=dpsv(j)+clop6v(3,j)
            oidpsv(j)=one/(one+dpsv(j))
            moidpsv(j)=mtc(j)/(one+dpsv(j))
            omoidpsv(j) = c1e3*((one-mtc(j))*oidpsv(j))
          end do
        end if

!       echo
        if ( dist_echo ) then
           open(unit=dist_echo_unit)
           rewind(dist_echo_unit)
           write(dist_echo_unit,'(" # ",A40,1PE25.18)') " total energy of synch part [MeV]: ", e0
           write(dist_echo_unit,'(" # ",A40,1PE25.18)') " momentum of synch part [MeV/c]: ", e0f
           write(dist_echo_unit,*) '#'
           write(dist_echo_unit,*) '# for every particle (j)'
           write(dist_echo_unit,*) '# xv(1), yv(1), xv(2), yv(2), sigmv, ejfv'
           do j = 1, napx
             write(dist_echo_unit,'(6(1X,1PE25.18))') xv(1, j), yv(1, j), xv(2,j), yv(2,j), sigmv(j), ejfv(j)
           end do
           close(dist_echo_unit)
        endif

      endif

      do 340 ia=1,napx,2
        if(idfor.ne.2.and.idfor.ne.3) then
!---------------------------------------  SUBROUTINE 'ANFB' IN-LINE
          if(st_quiet==0) write(lout,10050)
          tasia56=tas(ia,5,6)*c1m3
          bet0x2=tas(ia,1,3)**2+tas(ia,1,4)**2                           !hr05
          bet0z2=tas(ia,3,1)**2+tas(ia,3,2)**2                           !hr05
          bet0s1=tas(ia,5,5)**2+tasia56**2                               !hr05
          dsign=one
          rat=rat0
          if(tas(ia,3,3).lt.(-one*pieni)) rat=-one*rat                   !hr05
          if(rat.lt.(-one*pieni)) dsign=-one*one                         !hr05
          x11=ampv(ia)/(sqrt(bet0v(ia,1))+sqrt(abs(rat)*bet0x2))
          x13=(x11*dsign)*sqrt(abs(rat))                                 !hr05
          amp(2)=(dsign*real(1-iver,fPrec))*(abs(x11)*sqrt(bet0z2)+abs(x13)*sqrt(bet0v(ia,2)))                 !hr05
          x1(5)=zero
          x1(6)=dpsv(ia)*sqrt(bet0s1)
          chi=chi0*rad
          dchi=chid*rad
          do 320 i2=1,2
            i3=ia+i2-1
            sic=sin_mb(chi)
            coc=cos_mb(chi)
            x1(1)=x11*coc
            x1(2)=x11*sic
            x1(3)=x13*coc
            x1(4)=x13*sic
            do 300 ii=1,6
              x2(ii)=zero
              do 290 jj=1,6
                x2(ii)=x2(ii)+tas(ia,ii,jj)*x1(jj)
  290         continue
  300       continue
            if(iclo6.eq.1.or.iclo6.eq.2) then
              x2(2)=x2(2)/((one+x2(6))+clop6v(3,ia))                     !hr05
              x2(4)=x2(4)/((one+x2(6))+clop6v(3,ia))                     !hr05
            endif
            if(abs(bet0s1).le.pieni) x2(6)=dpsv(ia)
            if(iver.eq.1) then
              x2(3)=zero
              x2(4)=zero
            endif
            do 310 l=1,2
              ll=(l-1)*2
              xv(l,i3)=x2(1+ll)+exz(i2,1+ll)
              yv(l,i3)=x2(2+ll)+exz(i2,2+ll)
  310       continue
            sigmv(i3)=x2(5)+exz(i2,5)
            dpsv(i3)=x2(6)
            dpsic=dpsv(i3)+clop6v(3,ia)
            if(idp.eq.1.and.abs(ition).eq.1.and.iclo6.eq.0) then
              xv(1,i3)=xv(1,i3)+di0xs(ia)*dpsic
              xv(2,i3)=xv(2,i3)+di0zs(ia)*dpsic
              yv(1,i3)=yv(1,i3)+dip0xs(ia)*dpsic
              yv(2,i3)=yv(2,i3)+dip0zs(ia)*dpsic
            endif
            chi=chi+dchi
  320     continue
          if(st_quiet==0) write(lout,10260) ia,nms(ia)*izu0,dpsv(ia)
          if(st_quiet == 0) then
            write(lout,10060) xv(1,ia),yv(1,ia),xv(2,ia),yv(2,ia),sigmv(ia),dpsv(ia), &
                              xv(1,ia+1),yv(1,ia+1),xv(2,ia+1),yv(2,ia+1),sigmv(ia+1),dpsv(ia+1)
          end if
!---------------------------------------  END OF 'ANFB'
          if(iclo6.eq.2) then
            xv(1,ia)=xv(1,ia)+clo6v(1,ia)
            yv(1,ia)=yv(1,ia)+clop6v(1,ia)
            xv(2,ia)=xv(2,ia)+clo6v(2,ia)
            yv(2,ia)=yv(2,ia)+clop6v(2,ia)
            sigmv(ia)=sigmv(ia)+clo6v(3,ia)
            dpsv(ia)=dpsv(ia)+clop6v(3,ia)
            xv(1,ia+1)=xv(1,ia+1)+clo6v(1,ia)
            yv(1,ia+1)=yv(1,ia+1)+clop6v(1,ia)
            xv(2,ia+1)=xv(2,ia+1)+clo6v(2,ia)
            yv(2,ia+1)=yv(2,ia+1)+clop6v(2,ia)
            sigmv(ia+1)=sigmv(ia+1)+clo6v(3,ia)
            dpsv(ia+1)=dpsv(ia+1)+clop6v(3,ia)
            oidpsv(ia)=one/(one+dpsv(ia))
            oidpsv(ia+1)=one/(one+dpsv(ia+1))
          else
            xv(1,ia)=xv(1,ia)+(clov(1,ia)*real(idz(1),fPrec))*          &
     &real(1-idfor,fPrec)    !hr05
            yv(1,ia)=yv(1,ia)+(clopv(1,ia)*real(idz(1),fPrec))*         &
     &real(1-idfor,fPrec)   !hr05
            xv(2,ia)=xv(2,ia)+(clov(2,ia)*real(idz(2),fPrec))*          &
     &real(1-idfor,fPrec)    !hr05
            yv(2,ia)=yv(2,ia)+(clopv(2,ia)*real(idz(2),fPrec))*         &
     &real(1-idfor,fPrec)   !hr05
            xv(1,ia+1)=xv(1,ia+1)+(clov(1,ia)*real(idz(1),fPrec))*      &
     &real(1-idfor,fPrec)  !hr05
            yv(1,ia+1)=yv(1,ia+1)+(clopv(1,ia)*real(idz(1),fPrec))*     &
     &real(1-idfor,fPrec) !hr05
            xv(2,ia+1)=xv(2,ia+1)+(clov(2,ia)*real(idz(2),fPrec))*      &
     &real(1-idfor,fPrec)  !hr05
            yv(2,ia+1)=yv(2,ia+1)+(clopv(2,ia)*real(idz(2),fPrec))*     &
     &real(1-idfor,fPrec) !hr05
          endif
          ejfv(ia)=e0f*(one+dpsv(ia))
          ejfv(ia+1)=e0f*(one+dpsv(ia+1))
          ejv(ia)=sqrt(ejfv(ia)**2+nucm0**2)                               !hr05
          ejv(ia+1)=sqrt(ejfv(ia+1)**2+nucm0**2)                           !hr05
          epsa(1)=(ampv(ia)**2/bet0v(ia,1))                              !hr05
          epsa(2)=(amp(2)**2/bet0v(ia,2))                                !hr05

          moidpsv(ia)=mtc(ia)/(one+dpsv(ia))
          moidpsv(ia+1)=mtc(ia+1)/(one+dpsv(ia+1))
          omoidpsv(ia)=c1e3*((one-mtc(ia))*oidpsv(ia))
          omoidpsv(ia+1)=c1e3*((one-mtc(ia+1))*oidpsv(ia+1))
          nucm(ia)=nucm0
          nucm(ia+1)=nucm0

          if(st_quiet==0) write(lout,10020) ampv(ia),amp(2),epsa
        else if(idfor.eq.2) then
#ifndef CRLIBM
          read(13,*,iostat=ierro) xv(1,ia),yv(1,ia),xv(2,ia),yv(2,ia),  &
     &sigmv(ia),dpsv(ia),xv(1,ia+1),yv(1,ia+1),xv(2,ia+1),yv            &
     &(2,ia+1), sigmv(ia+1),dpsv(ia+1),e0,ejv(ia),ejv(ia+1)
#endif
#ifdef CRLIBM
          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ xv(1,ia)]"
             call prror(-1)
          endif
          xv(1,ia) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV xv(1,ia)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ yv(1,ia)]"
             call prror(-1)
          endif
          yv(1,ia) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV yv(1,ia)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ xv(2,ia)]"
             call prror(-1)
          endif
          xv(2,ia) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV xv(2,ia)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ yv(2,ia)]"
             call prror(-1)
          endif
          yv(2,ia) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV yv(2,ia)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ sigmv(ia)]"
             call prror(-1)
          endif
          sigmv(ia) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV sigmv(ia)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ dpsv(ia)]"
             call prror(-1)
          endif
          dpsv(ia) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV dpsv(ia)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ xv(1,ia+1)]"
             call prror(-1)
          endif
          xv(1,ia+1) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV xv(1,ia+1)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ yv(1,ia+1)]"
             call prror(-1)
          endif
          yv(1,ia+1) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV yv(1,ia+1)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ xv(2,ia+1)]"
             call prror(-1)
          endif
          xv(2,ia+1) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV xv(2,ia+1)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ yv(2,ia+1)]"
             call prror(-1)
          endif
          yv(2,ia+1) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV yv(2,ia+1)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)                                              &
     &            "Error when reading fort.13 [READ sigmv(ia+1)]"
             call prror(-1)
          endif
          sigmv(ia+1) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)                                              &
     &            "Error when reading fort.13 [CONV sigmv(ia+1)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)                                              &
     &            "Error when reading fort.13 [READ dpsv(ia+1)]"
             call prror(-1)
          endif
          dpsv(ia+1) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)                                              &
     &            "Error when reading fort.13 [CONV dpsv(ia+1)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ e0]"
             call prror(-1)
          endif
          e0 = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV e0]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ ejv(ia)]"
             call prror(-1)
          endif
          ejv(ia) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV ejv(ia)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ ejv(ia+1)]"
             call prror(-1)
          endif
          ejv(ia+1) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV ejv(ia+1)]"
             call prror(-1)
          endif

#endif
          if(ierro.ne.0) call prror(56)
          mtc(ia)=one
          mtc(ia+1)=one
          nucm(ia)=nucm0
          nucm(ia+1)=nucm0
          e0f=sqrt(e0**2-nucm0**2)                                         !hr05
          ejfv(ia)=sqrt(ejv(ia)**2-nucm(ia)**2)                               !hr05
          ejfv(ia+1)=sqrt(ejv(ia+1)**2-nucm(ia+1)**2)                           !hr05
          oidpsv(ia)=one/(one+dpsv(ia))
          oidpsv(ia+1)=one/(one+dpsv(ia+1))
          moidpsv(ia)=mtc(ia)/(one+dpsv(ia))
          moidpsv(ia+1)=mtc(ia+1)/(one+dpsv(ia+1))
          omoidpsv(ia)=c1e3*((one-mtc(ia))*oidpsv(ia))
          omoidpsv(ia+1)=c1e3*((one-mtc(ia+1))*oidpsv(ia+1))
        endif
        if (idfor /= 3 .and. st_quiet == 0) then
          write(lout,10090) xv(1,ia),yv(1,ia),xv(2,ia),yv(2,ia),sigmv(ia),dpsv(ia),xv(1,ia+1),&
                            yv(1,ia+1),xv(2,ia+1),yv(2,ia+1),sigmv(ia+1),dpsv(ia+1),e0,ejv(ia),ejv(ia+1)
        end if
        idam=3
        icode=0
        if(abs(xv(1,ia)).le.pieni.and.abs(yv(1,ia)).le.pieni) then
          idam=idam-1
        else
          icode=icode+1
        endif
        if(abs(xv(2,ia)).le.pieni.and.abs(yv(2,ia)).le.pieni) then
          idam=idam-1
        else
          icode=icode+2
        endif
        if(idp.eq.0.or.abs(ition).eq.0) then
          idam=idam-1
        else
          icode=icode+4
        endif
        if(idam.le.0) idam=1
        if(icode.le.0) icode=1
        ia2=(ia+1)/2
        if(ntwin.ne.2) then
          if(mod(ia+1,2).eq.0) then
            xau(1,1)= xv(1,ia)
            xau(1,2)= yv(1,ia)
            xau(1,3)= xv(2,ia)
            xau(1,4)= yv(2,ia)
            xau(1,5)=sigmv(ia)
            xau(1,6)= dpsv(ia)
            xau(2,1)= xv(1,ia+1)
            xau(2,2)= yv(1,ia+1)
            xau(2,3)= xv(2,ia+1)
            xau(2,4)= yv(2,ia+1)
            xau(2,5)=sigmv(ia+1)
            xau(2,6)= dpsv(ia+1)
            cloau(1)= clo6v(1,ia)
            cloau(2)=clop6v(1,ia)
            cloau(3)= clo6v(2,ia)
            cloau(4)=clop6v(2,ia)
            cloau(5)= clo6v(3,ia)
            cloau(6)=clop6v(3,ia)
            di0au(1)= di0xs(ia)
            di0au(2)=dip0xs(ia)
            di0au(3)= di0zs(ia)
            di0au(4)=dip0zs(ia)

            do ib2=1,6
              do ib3=1,6
                tau(ib2,ib3)=tasau(ia,ib2,ib3)
              end do
            end do

            call distance(xau,cloau,di0au,tau,dam1)
            dam(ia)=dam1
            dam(ia+1)=dam1
          endif !endif(mod(ia+1,2).eq.0)


!     Write header of track output file(s) used by postprocessing
!     for case ntwin.ne.2
#ifdef CR
          if (.not.restart) then
#endif
#ifndef STF
            call writebin_header(ia,ia,91-ia2,ierro,                    &
     &        cdate,ctime,progrm)
#ifdef CR
            flush(91-ia2)
            binrecs(ia2)=1
          endif
#endif
#endif
#ifdef STF
            call writebin_header(ia,ia,90,ierro,                        &
     &        cdate,ctime,progrm)
#ifdef CR
            flush(90)
            binrecs(ia2)=1
          endif
#endif
#endif
        else !ELSE for "if(ntwin.ne.2)"

!     Write header of track output file(s) used by postprocessing
!     for case ntwin.eq.2

#ifdef CR
          if (.not.restart) then
#endif
#ifndef STF
            call writebin_header(ia,ia+1,91-ia2,ierro,                  &
     &        cdate,ctime,progrm)
#ifdef CR
            flush(91-ia2)
            binrecs(ia2)=1
          endif
#endif
#endif
#ifdef STF
            call writebin_header(ia,ia+1,90,ierro,                      &
     &        cdate,ctime,progrm)
#ifdef CR
            flush(90)
            binrecs(ia2)=1
          endif
#endif
#endif
        endif !ENDIF (ntwin.ne.2)
        if(ierro.ne.0) then
          write(lout,*)
          write(lout,*) '*** ERROR ***,PROBLEMS WRITING TO FILE # : ',91&
     &-ia2
          write(lout,*) 'ERROR CODE : ',ierro
          write(lout,*)
          goto 520
        endif
  340 continue
#ifdef CR
      if (lhc.ne.9) binrec=1    ! binrec:
                                ! The maximum number of reccords writen for all tracking data files
                                ! Thus crbinrecs(:) .le. binrec
#endif
      if(e0.gt.pieni) then
        do j=1,napx
          rvv(j)=(ejv(j)*e0f)/(e0*ejfv(j))
        end do
      else
        call prror(79)
      endif

!-----/ End of initial distribution

  if(ithick.eq.1) then
!------ Compute matrices for linear tracking
    call envarsv(dpsv,moidpsv,rvv,ekv)
    if(idp.eq.0 .or. ition.eq.0) then
! ------- Only in case of thck4d
      call blocksv
    end if
  end if

#ifdef FLUKA
!     P.Garcia Ortega, A.Mereghetti and V.Vlachoudis, for the FLUKA Team
!     last modified: 26-08-2014
!     send napx to fluka
!     inserted in main code by the 'fluka' compilation flag
      if(fluka_enable) then
        write(lout,*) '[Fluka] Sending napx: ', napx
        write(fluka_log_unit,*) '# Sending napx: ', napx
        fluka_con = fluka_init_max_uid( napx )

        if (fluka_con .lt. 0) then
           write(lout,*) '[Fluka] Error: failed to send napx to fluka ',&
     &  napx
           write(fluka_log_unit, *) '# failed to send napx to fluka ',  &
     &  napx
           call prror(-1)
        end if

        write(lout,*) '[Fluka] Sending napx successful;'
        write(fluka_log_unit,*) '# Sending napx successful;'
        flush(lout)
        flush(fluka_log_unit)
      endif

!     A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!     last modified: 18-01-2016
!     initialise energy/momentum/rest mass of reference particle in mod_fluka
!         and synch magnetic rigidity with Fluka (for the time being, consider
!         only protons);
!     inserted in main code by the 'fluka' compilation flag
      if(fluka_enable) then
        write(lout,*) '[Fluka] Updating ref particle'
        write(fluka_log_unit,*) '# Updating ref particle'
        flush(lout)
        flush(fluka_log_unit)

        fluka_con = fluka_set_synch_part( e0, e0f, nucm0, aa0, zz0 )

        if (fluka_con .lt. 0) then
          write(lout, *) '[Fluka] Error: failed to update ref particle'
          write(fluka_log_unit, *) '# failed to update ref particle'
          call prror(-1)
        end if

        write(lout,*) '[Fluka] Updating ref successful;'
        write(fluka_log_unit,*) '# Updating ref particle successful;'
        flush(lout)
        flush(fluka_log_unit)
      endif

#endif

!     A.Mereghetti, P.Garcia Ortega and D.Sinuela Pastor, for the FLUKA Team
!     K. Sjobak, for BE/ABP-HSS
!     M. Fitterer, for FNAL
!     last modified: 21/02-2016
!     open units for dumping particle population or statistics
!     always in main code

      ! Initialise DUMP
      call dump_initialise

      ! ! ! Initialize SCATTER ! ! !
      if (scatter_active) then
         call scatter_initialise
      endif

#ifdef ROOT
! flush the root file
!  call SixTrackRootWrite()
#endif

!                                !
!     ****** TRACKING ******     !
!                                !
      write(lout,10200)
#ifdef DEBUG
!     call dumpbin('btrack',1,1)
!     call abend('btrack                                            ')
#endif
#ifdef DEBUG
                   !call system('../crmain  >> crlog')
#endif
      time1=0.
      call timex(time1)
! time1 is now pre-processing CPU
! note that this will be reset evry restart as we redo pre-processing
      pretime=time1-time0
!---------------------------------------  LOOP OVER TURNS TO BE TRACKED
      if(ithick.eq.0) call trauthin(nthinerr)
      if(ithick.eq.1) call trauthck(nthinerr)
#ifdef DEBUG
!     call dumpbin('atrack',1,1)
!     call abend('atrack                                            ')
#endif
      time2=0.
      call timex(time2)
! trtime is now the tracking time, BUT we must add other time for C/R
      trtime=time2-time1
#ifdef CR
! because now crpoint will write tracking time
! using time3 as a temp
! and crcheck/crstart will reset crtime3
      trtime=trtime+crtime3
#endif
      if(nthinerr.eq.3000) goto 520
      if(nthinerr.eq.3001) goto 460
!---------------------------------------  END OF LOOP OVER TURNS
  460 continue
#ifndef FLUKA
      napxto=0
#endif
! and set numx=nnuml (for writebin) NOT for LOST particles
! because all lost set nnuml=numl
      numx=nnuml
      id=0

#ifndef FLUKA
#ifdef CR
      if (.not.restart) then
! If restart is true , we haven't done any tracking
! and must be running from very last checkpoint
        write(93,*) 'Very last call to WRITEBIN?'
        write(93,*) 'numlmax,nnuml,numl',numlmax,nnuml,numl
        endfile (93,iostat=ierro)
        backspace (93,iostat=ierro)
        if (nnuml.eq.numl) then
! We REALLY have finished (or all particles lost)
! When all lost, nthinerr=3001, we set nnuml=numl
! and make sure we do the last WRITEBIN
          write(93,*) 'Very last call to WRITEBIN'
          endfile (93,iostat=ierro)
          backspace (93,iostat=ierro)
          call writebin(nthinerr)
          if(nthinerr.eq.3000) goto 520
        else
! I assume we are stopping because we have done nnuml turns
! which should be numlmax and do a writebin only if time
          write(93,*) 'Very last call to WRITEBIN?'
          write(93,*) 'numlmax,nnuml,nwri',numlmax,nnuml,nwri
          endfile (93,iostat=ierro)
          backspace (93,iostat=ierro)
          if(mod(nnuml,nwri).eq.0) then
            write(93,*) 'Very last call to WRITEBIN'
            endfile (93,iostat=ierro)
            backspace (93,iostat=ierro)
            call writebin(nthinerr)
            if(nthinerr.eq.3000) goto 520
          endif
        endif
! and do the very last checkpoint
        call callcrp()
      endif
#endif
#ifndef CR
      call writebin(nthinerr)
      if(nthinerr.eq.3000) goto 520
#endif
      ! If CR we have to worry about turns printed in fort.6
      ! If lost should be OK, otherwise we need to use nnuml instead
      ! of the numl in numxv/nnumxv???? Eric.
      ! where we reset [n]numxv to nnuml UNLESS particle lost
      ! Now we shall try using that fix at start of tracking
      write(lout,"(a)") str_divLine
      write(lout,"(a)") ""
      write(lout,"(a)") "    OOOOOOOOOOOOOOOOOOOOOOOOO"
      write(lout,"(a)") "    OO                     OO"
      write(lout,"(a)") "    OO  TRACKING COMPLETE  OO"
      write(lout,"(a)") "    OO                     OO"
      write(lout,"(a)") "    OOOOOOOOOOOOOOOOOOOOOOOOO"
      write(lout,"(a)") ""
      write(lout,"(a)") str_divLine
      write(lout,"(a)") ""
      write(lout,"(a)") "    PARTICLE SUMMARY:"
      write(lout,"(a)") ""

      do ia=1,napxo,2
        ie=ia+1
        ia2=(ie)/2
        napxto = napxto+numxv(ia)+numxv(ie)

        if(pstop(ia).and.pstop(ie)) then !-- BOTH PARTICLES LOST
          write(lout,10000) ia,nms(ia)*izu0,dp0v(ia),numxv(ia),abs(xvl(1,ia)),aperv(ia,1),abs(xvl(2,ia)),aperv(ia,2)
          write(lout,10000) ie,nms(ia)*izu0,dp0v(ia),numxv(ie),abs(xvl(1,ie)),aperv(ie,1),abs(xvl(2,ie)),aperv(ie,2)
          if(st_quiet == 0) write(lout,10280) xvl(1,ia),yvl(1,ia),xvl(2,ia),yvl(2,ia),sigmvl(ia),dpsvl(ia), &
            xvl(1,ie),yvl(1,ie),xvl(2,ie),yvl(2,ie),sigmvl(ie),dpsvl(ie),e0,ejvl(ia),ejvl(ie)
          write(12,10280,iostat=ierro) xvl(1,ia),yvl(1,ia),xvl(2,ia),yvl(2,ia),sigmvl(ia),dpsvl(ia), &
            xvl(1,ie),yvl(1,ie),xvl(2,ie),yvl(2,ie),sigmvl(ie),dpsvl(ie),e0,ejvl(ia),ejvl(ie)
          if(ierro /= 0) write(lout,"(2(a,i0))") "MAINCR> WARNING fort.12 has corrupted output probably due to lost particle ",&
            ia," or ",ie
        end if

        if(.not.pstop(ia).and.pstop(ie)) then !-- SECOND PARTICLE LOST
          id=id+1
          if(st_quiet == 0) then
            write(lout,10240) ia,nms(ia)*izu0,dp0v(ia),numxv(ia)
          else if(st_quiet == 1) then
            write(lout,10241) ia,nms(ia)*izu0,dp0v(ia),numxv(ia)
          end if
          write(lout,10000) ie,nms(ia)*izu0,dp0v(ia),numxv(ie),abs(xvl(1,ie)),aperv(ie,1),abs(xvl(2,ie)),aperv(ie,2)
          if(st_quiet==0) write(lout,10280) xv(1,id),yv(1,id),xv(2,id),yv(2,id),sigmv(id),dpsv(id), &
            xvl(1,ie),yvl(1,ie),xvl(2,ie),yvl(2,ie),sigmvl(ie),dpsvl(ie),e0,ejv(id),ejvl(ie)
          write(12,10280,iostat=ierro) xv(1,id),yv(1,id),xv(2,id),yv(2,id),sigmv(id),dpsv(id), &
            xvl(1,ie),yvl(1,ie),xvl(2,ie),yvl(2,ie),sigmvl(ie),dpsvl(ie),e0,ejv(id),ejvl(ie)
          if(ierro.ne.0) write(lout,"(a,i0)") "MAINCR> WARNING fort.12 has corrupted output probably due to lost particle ",ie
        end if

        if(pstop(ia).and..not.pstop(ie)) then !-- FIRST PARTICLE LOST
          id=id+1
          write(lout,10000) ia,nms(ia)*izu0,dp0v(ia),numxv(ia),abs(xvl(1,ia)),aperv(ia,1),abs(xvl(2,ia)),aperv(ia,2)
          if(st_quiet == 0) then
            write(lout,10240) ie,nms(ia)*izu0,dp0v(ia),numxv(ie)
          else if(st_quiet == 1) then
            write(lout,10241) ie,nms(ia)*izu0,dp0v(ia),numxv(ie)
          end if
          if(st_quiet==0) write(lout,10280) xvl(1,ia),yvl(1,ia),xvl(2,ia),yvl(2,ia),sigmvl(ia),dpsvl(ia), &
            xv(1,id),yv(1,id),xv(2,id),yv(2,id),sigmv(id),dpsv(id),e0,ejvl(ia),ejv(id)
          write(12,10280,iostat=ierro) xvl(1,ia),yvl(1,ia),xvl(2,ia),yvl(2,ia),sigmvl(ia),dpsvl(ia), &
            xv(1,id),yv(1,id),xv(2,id),yv(2,id),sigmv(id),dpsv(id),e0,ejvl(ia),ejv(id)
          if(ierro.ne.0) write(lout,"(a,i0)") "MAINCR> WARNING fort.12 has corrupted output probably due to lost particle ",ia
        end if

        if(.not.pstop(ia).and..not.pstop(ie)) then !-- BOTH PARTICLES STABLE
          id=id+1
          ig=id+1
          if(st_quiet == 0) then
            write(lout,10270) ia,ie,nms(ia)*izu0,dp0v(ia),numxv(ia)
          else if(st_quiet == 1) then
            write(lout,10271) ia,ie,nms(ia)*izu0,dp0v(ia),numxv(ia)
          end if
          if(st_quiet==0) write(lout,10280) xv(1,id),yv(1,id),xv(2,id),yv(2,id),sigmv(id),dpsv(id), &
            xv(1,ig),yv(1,ig),xv(2,ig),yv(2,ig),sigmv(ig),dpsv(ig),e0,ejv(id),ejv(ig)
          write(12,10280,iostat=ierro) xv(1,id),yv(1,id),xv(2,id),yv(2,id),sigmv(id),dpsv(id), &
            xv(1,ig),yv(1,ig),xv(2,ig),yv(2,ig),sigmv(ig),dpsv(ig),e0,ejv(id),ejv(ig)
          if(ierro.ne.0) write(lout,"(a)") "MAINCR> WARNING fort.12 has corrupted output although particles stable"
          id=ig
        end if
      end do
#endif
#ifdef FLUKA
      ! A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
      ! last modified: 17-07-2013
      ! print stable particles only
      ! inserted in main code by the 'fluka' compilation flag
      write(lout,"(a)") ""
      write(lout,"(a)") str_divLine
      if ( napxo .gt. 0 ) then
        write(lout,"(a)") ""
        write(lout,10350) napxo
        write(lout,"(a)") ""
        write(lout,10360) 'ID', 'GEN', 'WEIGHT', 'X [m]', 'XP []', 'Y [m]', 'YP[]', 'PC [GeV]', 'DE [eV]', 'DT [s]'
        write(lout,"(a)") ""
        do ia=1,napxo
          if(.not.pstop(ia)) then
            write(lout,10370) fluka_uid(ia),fluka_gen(ia),fluka_weight(ia), &
              xv(1,ia)*c1m3, yv(1,ia)*c1m3, xv(2,ia)*c1m3, yv(2,ia)*c1m3, &
              ejfv(ia)*c1m3,(ejv(ia)-e0)*c1e6,-c1m3*(sigmv(ia)/clight)*(e0/e0f)
          end if
        end do
      end if
#endif

! POSTPROCESSING (POSTPR)

! and we need to open fort.10 unless already opened for BOINC
#ifdef NAGFOR
#ifdef BOINC
  call boincrf('fort.10',filename)
#ifdef FIO
  open(10,file=filename,form='formatted',status='unknown',round='nearest',recl=8195)
#else
  open(10,file=filename,form='formatted',status='unknown',recl=8195)
#endif
#else
#ifdef FIO
  open(10,file='fort.10',form='formatted',status='unknown',round='nearest',recl=8195)
#else
  open(10,file='fort.10',form='formatted',status='unknown',recl=8195)
#endif
#endif
#else
#ifdef BOINC
  call boincrf('fort.10',filename)
#ifdef FIO
  open(10,file=filename,form='formatted',status='unknown',round='nearest')
#else
  open(10,file=filename,form='formatted',status='unknown')
#endif
#else
#ifdef FIO
  open(10,file='fort.10',form='formatted',status='unknown',round='nearest')
#else
  open(10,file='fort.10',form='formatted',status='unknown')
#endif
#endif
#endif

#ifndef FLUKA
#ifndef STF
        iposc=0
        if(ipos.eq.1) then !Variable IPOS=1 -> postprocessing block present in fort.3
          do 480 ia=1,napxo,2
            ia2=(ia+1)/2
            iposc=iposc+1
#ifndef CR
            call postpr(91-ia2) !Postprocess file "fort.(91-ia2)"
#endif
#ifdef CR
            write(93,*) 'Calling POSTPR nnuml=',nnuml
            endfile (93,iostat=ierro)
            backspace (93,iostat=ierro)
            call postpr(91-ia2,nnuml)
#endif
  480     continue
          if(iposc.ge.1) call sumpos
        endif !END if(ipos.eq.1)
        goto 520 !Done postprocessing

  490   if(ipos.eq.1) then !GOTO here if(napx.le.0.or.imc.le.0) (skipping tracking)
          ndafi2=ndafi
          do 500 ia=1,ndafi2
            if(ia.gt.ndafi) goto 510
#ifndef CR
            call postpr(91-ia)
#endif
#ifdef CR
            write(93,*) 'Calling POSTPR nnuml=',nnuml
            endfile (93,iostat=ierro)
            backspace (93,iostat=ierro)
            call postpr(91-ia,nnuml)
#endif
  500     continue
  510     if(ndafi.ge.1) call sumpos
        endif
#endif
#ifdef STF
        iposc=0
        if(ipos.eq.1) then !Variable IPOS=1 -> postprocessing block present in fort.3
           do 480 ia=1,napxo,2
              iposc=iposc+1
#ifndef CR
              call postpr(ia) !Postprocess particle ia (and ia+1 if ntwin=2)
#endif
#ifdef CR
              write(93,*) 'Calling POSTPR nnuml=',nnuml
              endfile (93,iostat=ierro)
              backspace (93,iostat=ierro)
              call postpr(ia,nnuml)
#endif
  480      continue
          if(iposc.ge.1) call sumpos
        endif
        goto 520 !Done postprocessing

  490   if(ipos.eq.1) then !GOTO here if(napx.le.0.or.imc.le.0) (skipping tracking)
          ndafi2=ndafi
          do 500 ia=1,(2*ndafi2),2
            if(ia.gt.ndafi) goto 510
#ifndef CR
            call postpr(ia)
#endif
#ifdef CR
            write(93,*) 'Calling POSTPR nnuml=',nnuml
            endfile (93,iostat=ierro)
            backspace (93,iostat=ierro)
            call postpr(ia,nnuml)
#endif
  500     continue
  510     if(ndafi.ge.1) call sumpos
        endif
#endif

 520  continue !Finished postprocessing (POST in fort.3)

!     start fma
      if(fma_flag) then
        write(lout,*)'Calling FMA_POSTPR'
        call fma_postpr
      endif
!--HPLOTTING END
      if(ipos.eq.1.and.                                                 &
     &(idis.ne.0.or.icow.ne.0.or.istw.ne.0.or.iffw.ne.0)) then
        call igmeta(999,0)
        call hplend
      endif
#endif

#ifdef FLUKA
!     A.Mereghetti, for the FLUKA Team
!     last modified: 28-05-2014
!     collect a couple of goto statements, sending code flow
!       to different plotting points, which are not actually
!       inserted
!     inserted in main code by the 'fluka' compilation flag
 490  continue
 520  continue
  call fluka_close
#endif
#ifdef FFIELD
  ! Modification by B.DALENA and T.PUGNAT
  call ffield_mod_end()
#endif
      time3=0.
      call timex(time3)
! Note that crpoint no longer destroys time2
      posttime=time3-time2
#ifdef DEBUG
      write(lout,*) 'BUG:',time3,time2,pretime,trtime,posttime
#ifdef CR
      write(93,*)   'BUG:',time3,time2,pretime,trtime,posttime
#endif
#endif
#ifdef CR
! and TRY a FIX for napxto
!     if (nnuml.ne.numl) then
!       napxto=0
!       write(lout,*) 'numl=',numl,' nnuml=',nnuml
! We may have stopped because of numlmax
!       do ia=1,napxo
!         if (numxv(ia).eq.numl) then
! assumed stable
!     write(lout,*) 'ia=',ia,nnuml
!           napxto=napxto+nnuml
!         else
! assumed lost
!     write(lout,*) 'ia=',ia,' numxv=',numxv
!           napxto=napxto+numxv(ia)
!         endif
!       enddo
!     endif
#endif

  ! Get grand total including post-processing
  tottime = (pretime+trtime)+posttime
  write(lout,"(a)")         ""
  write(lout,"(a)")         str_divLine
  write(lout,"(a)")         ""
  write(lout,"(a)")         "    Computing Time Summary"
  write(lout,"(a)")         "  =========================="
  write(lout,"(a,f12.3,a)") "    Preparating Calculations: ",pretime, " second(s)"
  write(lout,"(a,f12.3,a)") "    Particle Tracking:        ",trtime,  " second(s)"
  write(lout,"(a,f12.3,a)") "    Post Processing:          ",posttime," second(s)"
  write(lout,"(a,f12.3,a)") "    Total Time Used:          ",tottime, " second(s)"
  write(lout,"(a,i8)")      "    Particle Turns:           ",napxto
  write(lout,"(a)")         ""
  write(lout,"(a)")         str_divLine

  if (zipf_numfiles.gt.0) then
    call zipf_dozip
  endif
#ifdef HDF5
  if(h5_isReady) then
    call h5_writeAttr(h5_rootID,"PreTime",  pretime)
    call h5_writeAttr(h5_rootID,"TrackTime",trtime)
    call h5_writeAttr(h5_rootID,"PostTime", posttime)
    call h5_writeAttr(h5_rootID,"TotalTime",tottime)
  end if
#endif
#ifdef ROOT
  if(root_flag) then
    call RunTimeRootWrite(pretime, trtime, posttime)
    call SixTrackRootExit()
  end if
#endif
  call alloc_exit
  call closeUnits ! Must be last as it also closes fort.6
! ----------------------------------------------------------------------
!   We're done in maincr, no error :)
! ----------------------------------------------------------------------
#ifdef CR
  call abend('                                                  ')
#endif
#ifndef CR
  stop
#endif
10000 format(/4x,"Tracking ended abnormally for particle: ",i0,         &
             /4x,"Random seed:        ",i8,                             &
             /4x,"Momentum deviation: ",f14.5,                          &
             /4x,"Lost in revolution: ",i8,                             &
             /4x,"Horiz:  Amplitude = ",ES23.16,"  Aperture = ",f15.3   &
             /4x,"Vert:   Amplitude = ",ES23.16,"  Aperture = ",f15.3/)
10020 format(/t10,'UNCOUPLED AMPLITUDES AND EMITTANCES:', /t10,         &
     &'AMPLITUDE-X = ',f15.3,10x,'AMPLITUDE-Y = ',f15.3, '  MM'/t10,    &
     &'EMITTANCE-X = ',f15.3,10x,'EMITTANCE-Y = ',f15.3, '  PI*MRAD*MM')
10025 format(/t10,'Run started from binary dump file # 32')
10050 format(//131('-')//t10,27('O')/t10,2('O'),23x,2('O')/t10,         &
     &'OO  INITIAL COORDINATES  OO'/ t10,2('O'),23x,2('O')/t10,27('O')  &
     &//131('-')//)
10060 format(/5x,'---- TWIN-TRAJECTORIES NO CL.ORBIT ADDED'/ 5x,'/X1  /'&
     &,f47.33/5x,'/XP1 /',f47.33/ 5x,'/Y1  /',f47.33/5x,'/YP1 /',f47.33/&
     &5x,'/SIG1/',f47.33/5x,'/DP1 /',f47.33/ 5x,'/X2  /',f47.33/5x,     &
     &'/XP2 /',f47.33/ 5x,'/Y2  /',f47.33/5x,'/YP2 /',f47.33/ 5x,       &
     &'/SIG2/',f47.33/5x,'/DP2 /',f47.33/)
10070 format(/131('-'))
10080 format(/t10,'REL. MOMENTUM DEVIATION=',f19.16/ t8,                &
     &'========================================')
10090 format(/5x,'---- INITIAL COORD. OF TWIN-TRAJECTORIES'/ 15(10x,f47.&
     &33/))
10110 format(/5x,'---- CLOSED ORBIT AND DECOUPLING (1=COU,0=DECOU)'/ 5x,&
     &'/CLX /',f47.33/5x,'/CLXP/',f47.33/ 5x,'/CLY /',f47.33/5x,'/CLYP/'&
     &,f47.33/ 5x,'/DCX / ',i4/5x,'/DCY / ',i4/ 5x,'/IVER /',i4/ 5x,    &
     &'/IDFOR/',i4/ 5x,'/ICLO6/',i4/ 5x,'/ITION/',i4/5x/)
10120 format(/5x,'---- CLOSED ORBIT AND DECOUPLING (1=COU,0=DECOU)'/ 5x,&
     &'/CLX /',f47.33/5x,'/CLXP/',f47.33/ 5x,'/CLY /',f47.33/5x,'/CLYP/'&
     &,f47.33/ 5x,'/CLS /',f47.33/5x,'/CLSP/',f47.33/ 5x,'/DCX / ',i4/5 &
     &x,'/DCY / ',i4/ 5x,'/IVER /',i4/ 5x,'/IDFOR/',i4/ 5x,'/ICLO6/',i4/&
     &5x,'/ITION/',i4/5x/)
10150 format(/t10,'TRACKING WITH SYNCHROTRON OSCILLATIONS'/ 15x,        &
     &'ACCELERATION WITH PHASE = ',f8.4/ t15,                           &
     &'       TUNE             CLO                CLOP           ',     &
     &'     BET0             ALF0           GAMMA      '//              &
     &t10,'  X  ',6(1x,ES17.10)/                                        &
     &t69,3(1x,ES17.10)/t69,3(1x,ES17.10)/                              &
     &t10,'  Y  ',6(1x,ES17.10))
10160 format(/t10,'TRACKING WITH SYNCHROTRON OSCILLATIONS'/ 15x,        &
     &'ACCELERATION WITH PHASE = ',f8.4/ t15,                           &
     &'       TUNE             CLO                CLOP           ',     &
     &'     BET0             ALF0           GAMMA      '//              &
     &t10,'  X  ',6(1x,ES17.10)/                                        &
     &t69,3(1x,ES17.10)/t69,3(1x,ES17.10)/                              &
     &t10,'  Y  ',6(1x,ES17.10)/                                        &
     &t69,3(1x,ES17.10)/t69,3(1x,ES17.10)/                              &
     &t10,'  S  ',3(1x,ES17.10),3(1x,ES17.10)/                          &
     &t69,3(1x,ES17.10)/t69,3(1x,ES17.10))
10170 format(/t10,'TRACKING FOR CONSTANT MOMENTUM DEVIATION'// 15x,     &
     &'------ NO ACCELERATION ------'// t15,                            &
     &'       TUNE         CLO            CLOP           ',             &
     &'   BET0           ALF0           GAMMA      '//                  &
     &t10,'  X  ',f14.10,2(1x,g15.8),1x,f15.9,1x,f15.10,f15.9/          &
     &t62,f15.9,1x,f15.10,f15.9/                                        &
     &t10,'  Y  ',f14.10,2(1x,g15.8),1x,f15.9,1x,f15.10,f15.9/          &
     &t62,f15.9,1x,f15.10,f15.9/)
10180 format(t5//t5,'BACK-TRACKING'/ t5, '============='//)
10190 format(t10,'TRACKING FOR CONSTANT MOMENTUM DEVIATION'// 15x,      &
     &'ACCELERATION WITH PHASE = ',f8.4/ t15,                           &
     &'       TUNE         CLO            CLOP           ',             &
     &'   BET0           ALF0           GAMMA      '//                  &
     &t10,'  X  ',f14.10,2(1x,g15.8),1x,f15.9,1x,f15.10,f15.9/          &
     &t62,f15.9,1x,f15.10,f15.9/                                        &
     &t10,'  Y  ',f14.10,2(1x,g15.8),1x,f15.9,1x,f15.10,f15.9/          &
     &t62,f15.9,1x,f15.10,f15.9/)
10200 format(//131('-')//t10,16('O')/t10,2('O'),12x,2('O')/t10,         &
     &'OO  TRACKING  OO', /t10,2('O'),12x,2('O')/t10,16('O')//131('-')//&
     &)
10210 format(/t10,'TRACKING WITH SYNCHROTRON OSCILLATIONS'/ 15x,        &
     &'------ NO ACCELERATION ------'// t15,                            &
     &'       TUNE             CLO                CLOP           ',     &
     &'     BET0             ALF0           GAMMA      '//              &
     &t10,'  X  ',6(1x,ES17.10)/                                        &
     &t69,3(1x,ES17.10)/t69,3(1x,ES17.10)/                              &
     &t10,'  Y  ',6(1x,ES17.10))
10220 format(/t10,'TRACKING WITH SYNCHROTRON OSCILLATIONS'/ 15x,        &
     &'------ NO ACCELERATION ------'// t15,                            &
     &'       TUNE             CLO                CLOP           ',     &
     &'     BET0             ALF0           GAMMA      '//              &
     &t10,'  X  ',6(1x,ES17.10)/                                        &
     &t69,3(1x,ES17.10)/t69,3(1x,ES17.10)/                              &
     &t10,'  Y  ',6(1x,ES17.10)/                                        &
     &t69,3(1x,ES17.10)/t69,3(1x,ES17.10)/                              &
     &t10,'  S  ',3(1x,ES17.10),3(1x,ES17.10)/                          &
     &t69,3(1x,ES17.10)/t69,3(1x,ES17.10))
10230 format(t10,'NO OPTICAL SOLUTION FOR',2x,f19.16,2x,'RELATIVE MOMENTUM DEVIATION')
10240 format(/4x,"Particle ",i7,"             Stable. Random seed: ",i0," Momentum deviation: ",g12.5," Revolution: ",i0/)
10241 format( 4x,"Particle ",i7,"             Stable. Random seed: ",i0," Momentum deviation: ",g12.5," Revolution: ",i0)
10260 format(/4x,"Particle ",i7,"                     Random seed: ",i0," Momentum deviation: ",g12.5/)
10270 format(/4x,"Particle ",i7," and ",i7, " Stable. Random seed: ",i0," Momentum deviation: ",g12.5," Revolution: ",i0/)
10271 format( 4x,"Particle ",i7," and ",i7, " Stable. Random seed: ",i0," Momentum deviation: ",g12.5," Revolution: ",i0)
10280 format(4x,f47.33)
10320 format(//131('-')//t10,'DATA BLOCK FLUCTUATIONS OF MULTIPOLES'//  &
     &t10,'RANDOM STARTING NUMBER=  ',i20/ t10,                         &
     &'RANDOM NUMBERS GENERATED:',i20/ t10,'MEAN VALUE=',f15.7,         &
     &'  -   DEVIATION=',f15.7)
10330 format(/10x,'ERROR IN OPENING FILES')
#ifdef FLUKA
!     A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!     last modified: 17-07-2013
!     print stable particles only (format directives)
!     inserted in main code by the 'fluka' compilation flag
10350 format(4X,I8,1X,'SURVIVING PARTICLES:')
10360 format(2(1X,A8),8(1X,A16))
10370 format(2(1X,I8),8(1X,1PE16.9))
10380 format(10x,f47.33)
#endif
end program maincr
