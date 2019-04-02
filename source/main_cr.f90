! ============================================================================ !
!
!  SIXTRACK MAIN PROGRAM
! =======================
!  Authors:       See README.md
!  Change Log:    See CHANGELOG.md
!  Licence:       LGPL v2.1, See COPYING and LICENCE
!  Documentation: See doc folder
!  Website:       http://sixtrack.web.cern.ch
!
! ============================================================================ !
program maincr

  use floatPrecision
  use mod_units
  use string_tools
  use mathlib_bouncer
  use physical_constants
  use numerical_constants

  use dynk,    only : dynk_izuIndex
  use fma,     only : fma_postpr, fma_flag
  use dump,    only : dump_initialise, dumpclo,dumptas,dumptasinv
  use zipf,    only : zipf_numfiles, zipf_dozip
  use scatter, only : scatter_init

  use, intrinsic :: iso_fortran_env, only : output_unit
  use mod_meta
  use mod_time
  use aperture
  use mod_ranecu
  use mod_particles
  use mod_geometry,   only : geom_calcDcum
  use mod_alloc,      only : alloc_init
  use mod_fluc,       only : fluc_randomReport, fluc_errAlign, fluc_errZFZ
  use postprocessing, only : postpr, writebin_header, writebin
  use read_write,     only : writeFort12, readFort13, readFort33

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
  use mod_common_main
  use mod_commons
  use mod_common_track

  use mod_hions
  use mod_dist
  use matrix_inv
  use aperture
  use wire
  use mod_version
#ifdef HASHLIB
  use mod_hash
#endif

  implicit none

interface

  subroutine envarsv(dpsv,oidpsv,rvv,ekv)

    use floatPrecision
    use parpro
    use mod_common_da

    implicit none

    real(kind=fPrec) :: dpsv(npart)
    real(kind=fPrec) :: oidpsv(npart)
    real(kind=fPrec) :: rvv(npart)
    real(kind=fPrec), allocatable, intent(inout) :: ekv(:,:)
  end subroutine envarsv

end interface

  ! "Old" variables
  integer i,itiono,i2,i3,ia,ia2,iar,iation,ib,ib0,ib1,ib2,ib3,id,ie,ig,ii,im,iposc,ix,izu,j,j2,jj,  &
    k,kpz,kzz,l,ll,ncorruo,ncrr,nd,nd2,ndafi2,nerror,nlino,nlinoo,nmz,nthinerr
  real(kind=fPrec) alf0s1,alf0s2,alf0s3,alf0x2,alf0x3,alf0z2,alf0z3,amp00,bet0s1,bet0s2,bet0s3,     &
    bet0x2,bet0x3,bet0z2,bet0z3,chi,coc,dam1,dchi,ddp1,dp0,dp00,dp10,dpsic,dps0,dsign,gam0s1,gam0s2,&
    gam0s3,gam0x1,gam0x2,gam0x3,gam0z1,gam0z2,gam0z3,phag,r0,r0a,rat0,sic,tasia56,tasiar16,tasiar26,&
    tasiar36,tasiar46,tasiar56,tasiar61,tasiar62,tasiar63,tasiar64,tasiar65,taus,x11,x13,damp
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
#ifdef FLUKA
  integer fluka_con
#endif

  ! New Variables
  character(len=:), allocatable :: featList, compName
#ifndef STF
  character(len=7)  tmpFile
#endif
  character(len=23) timeStamp
  character(len=8)  tsDate
  character(len=10) tsTime

  logical fErr ! For file units

  ! ---------------------------------------------------------------------------------------------- !
  errout = 0 ! Set to nonzero before calling abend in case of error.
#ifdef CR
  lout = 92
#else
  lout = output_unit
#endif

#ifdef BOINC
  call boinc_init
! call boinc_init_graphics
#endif
  call f_initUnits
  call meta_initialise ! The meta data file.
  call time_initialise ! The time data file. Need to be as early as possible as it sets cpu time 0.
  call alloc_init      ! Initialise mod_alloc
  call allocate_arrays ! Initial allocation of memory
#ifdef HASHLIB
  call hash_initialise
#endif

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
  call h5_initHDF5
#endif
#ifdef BOINC
  featList = featList//" BOINC"
#endif
#ifdef LIBARCHIVE
  featList = featList//" LIBARCHIVE"
#endif
#ifdef PYTHIA
  featList = featList//" PYTHIA"
#endif

  compName = "default"
#ifdef GFORTRAN
  compName = "gfortran"
#endif
#ifdef IFORT
  compName = "ifort"
#endif
#ifdef NAGFOR
  compName = "nagfor"
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
  ! Goes here after unzip for BOINC
#endif
  ! Very first get rid of any previous partial output
  call f_close(lout)
  call f_open(unit=lout,file="fort.92",formatted=.true.,mode="rw",err=fErr,status="replace")

  ! Now position the checkpoint/restart logfile=93
  call f_open(unit=93,file="fort.93",formatted=.true.,mode="rw",err=fErr)
606 continue
  read(93,"(a1024)",end=607) arecord
  goto 606
607 continue
  backspace(93,iostat=ierro)
#ifdef BOINC
  ! and if BOINC issue an informatory message
  if(start) then
    write(93,"(a)") "SIXTRACR> starts for the very first time"
  else
    write(93,"(a)") "SIXTRACR> retry after unzip of Sixin.zip"
  end if
#endif
  ! Now we see if we have a fort.6 which implies that we can perhaps just restart using all exisiting files
  ! including the last checkpoints. If not, we just do a start (with an unzip for BOINC)
  ! call f_open(unit=6,file="fort.6",formatted=.true.,mode="w",err=fErr,status="old")
  ! if(fErr) goto 602
  ! stxt = "SIXTRACR reruns on: "
  call f_open(unit=output_unit,file="fort.6",formatted=.true.,mode="rw",err=fErr,status="old")
  if(fErr) then
#ifdef BOINC
    ! No fort.6 so we do an unzip of Sixin.zip
    ! BUT ONLY IF WE HAVE NOT DONE IT ALREADY
    ! and CLOSE 92 and 93
    if(start) then
      start=.false.
      call f_close(92)
      call f_close(93)
      ! Now, if BOINC, after no fort.6, call UNZIP Sixin.zip
      ! Name hard-wired in our boinc_unzip_.
      ! Either it is only the fort.* input data or it is a restart.
      call boincrf("Sixin.zip",filename)
      ! This function expects a normal, trimmed fortran string; it will do the zero-padding internally.
      call f_read_archive(trim(filename),".")
      goto 611
    end if
    call f_open(unit=output_unit,file="fort.6",formatted=.true.,mode="rw",err=fErr)
#else
    call f_open(unit=output_unit,file="fort.6",formatted=.true.,mode="rw",err=fErr,status="new")
#endif
    ! Set up start message depending on fort.6 or not
    stxt = "SIXTRACR> starts on: "
  else
    ! Set up start message depending on fort.6 or not
    stxt = "SIXTRACR> reruns on: "
    rerun=.true.
  end if
  call f_open(unit=95,file="fort.95",formatted=.false.,mode="rw",err=fErr,status="old")
  if(fErr) then
    call f_open(unit=95,file="fort.95",formatted=.false.,mode="rw",err=fErr,status="new")
  else
    fort95 = .true.
  end if
  call f_open(unit=96,file="fort.96",formatted=.false.,mode="rw",err=fErr,status="old")
  if(fErr) then
    call f_open(unit=96,file="fort.96",formatted=.false.,mode="rw",err=fErr,status="new")
  else
    fort96 = .true.
  end if
  call f_open(unit=91,file="fort.91",formatted=.true.,mode="rw",err=fErr)
#else
  lout = output_unit
#endif

  ! Open Regular File Units
  call f_open(unit=18,file="fort.18",formatted=.true., mode="rw",err=fErr) ! DA file
  call f_open(unit=19,file="fort.19",formatted=.true., mode="r", err=fErr) ! DA file
  call f_open(unit=20,file="fort.20",formatted=.true., mode="w", err=fErr) ! DA file
  call f_open(unit=21,file="fort.21",formatted=.true., mode="w", err=fErr) ! DA file
  call f_open(unit=31,file="fort.31",formatted=.true., mode="w", err=fErr)

#ifdef STF
  ! Open Single Track File
  call f_open(unit=90,file="singletrackfile.dat",formatted=.false.,mode="rw",err=fErr)
#else
  ! Open binary files 59 to 90 for particle pair 1 to 32
  do i=59,90
    write(tmpFile,"(a5,i0)") "fort.",i
    call f_open(unit=i,file=tmpFile,formatted=.false.,mode="rw",err=fErr)
  end do
#endif

  call f_open(unit=111,file="fort.111",formatted=.false.,mode="rw",err=fErr) ! DA file, binary

#ifdef DEBUG
  ! call f_open(unit=99 ,file="dump",  formatted=.false.,mode="rw",err=fErr)
  ! call f_open(unit=100,file="arrays",formatted=.false.,mode="rw",err=fErr)
#endif

  call time_timeStamp(time_afterFileUnits)

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
  write(lout,"(a)") "    Compiler:     "//trim(compName)
  write(lout,"(a)") "    Built With:   "//trim(adjustl(featList))
  write(lout,"(a)") "    Start Time:   "//timeStamp
  write(lout,"(a)") ""
  write(lout,"(a)") str_divLine

  call meta_write("SixTrackVersion", trim(version))
  call meta_write("ReleaseDate",     trim(moddate))
  call meta_write("GitHash",         trim(git_revision))
  call meta_write("Compiler",        trim(compName))
  call meta_write("Features",        trim(adjustl(featList)))
  call meta_write("StartTime",       timeStamp)

#ifdef CR
  ! Log start messages
  write(93,"(a)") ""
  write(93,"(a)") "SIXTRACR> MAINCR Starting"
  write(93,"(a)") stxt//timeStamp
  flush(93)
#endif

  call time_timerStart
  call time_timerCheck(time0)
  progrm = "SIXTRACK"

#ifdef ROOT
  call SixTrackRootFortranInit
#endif

#ifdef FLUKA
  call fluka_mod_init(npart_initial, nele_initial, clight)
#endif

#ifdef FFIELD
  ! Modification by B.DALENA and T.PUGNAT
  call ffield_mod_init
#endif

  call daten
  call time_timeStamp(time_afterDaten)

#ifdef HDF5
  if(h5_isActive) then
    call h5_openFile
    call h5_writeSimInfo
  end if
#endif

  if (ithick == 1) call allocate_thickarrays

#ifdef CR
  checkp=.true.
  call crcheck
  call time_timeStamp(time_afterCRCheck)
#endif
  if(ithick == 1) write(lout,"(a)") "MAINCR> Structure input file has -thick- linear elements"
  if(ithick == 0) write(lout,"(a)") "MAINCR> Structure input file has -thin- linear elements"

  call scatter_init
  call aperture_init

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

  if(nprint == 1 .and. ipos == 1) then
    ! Open fort.14 and fort.15 for postprocessing
    call f_open(unit=14,file="fort.14",formatted=.true.,mode="w")
    call f_open(unit=15,file="fort.15",formatted=.true.,mode="w")
  end if

  ! Postprocessing is on, but there are no particles
  if(ipos.eq.1.and.napx.eq.0) then
    ! Now we open fort.10 unless already opened for BOINC
    call f_open(unit=10, file="fort.10", formatted=.true., mode="rw",err=fErr,recl=8195)
    call f_open(unit=110,file="fort.110",formatted=.false.,mode="w", err=fErr)

#ifndef STF
    do i=1,ndafi !ndafi = number of files to postprocess (set by fort.3)
#ifndef CR
      call postpr(91-i)
#else
      write(93,"(a,i0)") "MAINCR> Calling POSTPR nnuml = ",nnuml
      flush(93)
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
      flush(93)
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

  itra  = 2
  amp00 = amp(1)
  damp  = zero
  if(napx /= 1) damp=((amp00-amp0)/real(napx-1,fPrec))/two
  napx  = 2*napx
  call expand_arrays(nele, napx*imc, nblz, nblo)

  ! Log some meta data
  meta_nPartInit = napx*imc
  call meta_write("NumParticles",         napx*imc)
  call meta_write("NumTurns",             numl)
  call meta_write("NumSingleElements",    il)
  call meta_write("NumBlockElements",     mblo)
  call meta_write("NumStructureElements", mbloz)

  aperture_napxStart=napx
  iation=abs(ition)
  ib0=0
  dp00=dp1
  if(napx <= 0 .or. imc <= 0) goto 490

  ! A.Mereghetti (CERN, BE-ABP-HSS), 06-03-2018
  ! possible to re-shuffle lattice structure
  call orglat
  call geom_calcDcum

  ! A.Mereghetti (CERN, BE-ABP-HSS), 16-12-2016
  ! initialise aperture of first and last elements of sequence
  if (limifound) then
    write(lout,"(a)") "MAINCR> Check that beginning/end of lattice structure is assigned aperture markers."
    call contour_aperture_markers( iu, 1, .false. )
  end if

#ifdef FLUKA
  if (fluka_enable) then
    call check_coupling_integrity
    call check_coupling_start_point
  end if
#endif

  ! dump aperture model
  if (ldmpaper) then
#ifdef HDF5
    if(h5_useForAPER) then
      call dump_aperture_model_hdf5
    else
      call dump_aperture_model
    end if
#else
    call dump_aperture_model
#endif
  end if
  ! dump x-sections at specific locations
  if (mxsec.gt.0) call dump_aperture_xsecs
  ! map errors, now that the sequence is no longer going to change
  call ord
  if(allocated(zfz)) call fluc_randomReport

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
    if(root_flag .and. root_FLUKA == 1) then
      call root_FLUKA_DumpInsertions
    end if
#endif

  end if
#endif

  do l=1,2
    clo0(l)=clo(l)
    clop0(l)=clop(l)
  end do
  call clorb(zero)

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

  ! beam-beam element
  nlino = nlin
  nlin  = 0
  if(nbeam.ge.1) then
    do i=1,nele
      if(kz(i).eq.20) then
        nlin=nlin+1
        if(nlin.gt.nele) then
          write(lout,"(a)") "MAINCR> ERROR Too many elements for linear optics write-out"
          call prror(-1)
        end if
        bezl(nlin)=bez(i)
      end if
    end do
  end if
  if(isub == 1) call subre(dp1)
  if(ise  == 1) call search(dp1)

  !! Initialize kicks
  izu=0
  do i=1,iu
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
    smiv(i)=sm(ix)+smizf(i) ! Also in initalize_element!
    smi(i)=smiv(i)          ! Also in initalize_element!
    izu=izu+1
    xsiv(i)=xpl(ix)+zfz(izu)*xrms(ix)
    xsi(i)=xsiv(i)
    izu=izu+1
    zsiv(i)=zpl(ix)+zfz(izu)*zrms(ix)
    zsi(i)=zsiv(i)
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
    end if

    ! MULTIPOLE BLOCK
    if(kzz.eq.11) then
      dynk_izuIndex(ix)=izu

      ! Initialize multipoles, combining settings from fort.2 with
      ! coefficients from MULT and random values from FLUC.
      ! Used in program maincr and from initialize_element.

      if(abs(ek(ix)).le.pieni) cycle
      nmz=nmu(ix)
      if(nmz.eq.0) then
        izu=izu+2*mmul
        cycle
      end if
      im=irm(ix)
      do k=1,nmz
        izu=izu+1
        amultip(k,i) = zfz(izu) !To make it easier for Dynk later on
        aaiv(k,i)=(ak0(im,k)+(amultip(k,i)*aka(im,k)))
        izu=izu+1
        bmultip(k,i) = zfz(izu)
        bbiv(k,i)=(bk0(im,k)+(bmultip(k,i)*bka(im,k)))
      end do
      izu=izu+2*mmul-2*nmz
    end if
  end do

#ifdef DEBUG
!     call dumpbin('ado 150',150,150)
!     call abend('ado 150                                           ')
#endif
        dp1=zero
        if(ichrom > 1) then
          itiono=ition
          ition=0
          call chromda
          ition=itiono
          do ncrr=1,iu
            ix=ic(ncrr)
            if(ix.gt.nblo) ix=ix-nblo
            if(ix.eq.is(1).or.iratioe(ix).eq.is(1)) then
              smiv(ncrr)=smi(ncrr)
            else if(ix.eq.is(2).or.iratioe(ix).eq.is(2)) then
              smiv(ncrr)=smi(ncrr)
            endif
          enddo
        else
          itiono = 0 ! -Wmaybe-uninitialized
        endif
        dp1=dp00
        dp0=dp00
        if(imc > 1) then
          ddp1 = (two*dp0)/(real(imc,fPrec)-one)
        else
          ddp1 = zero ! -Wmaybe-uninitialized
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
                  smiv(ncrr)=smi(ncrr)
                else if(ix.eq.iq(2).or.iratioe(ix).eq.iq(2)) then
                  smiv(ncrr)=smi(ncrr)
                endif
              enddo
            endif

            do 190 ib1=1,napx
              ib3=ib1+(ib-1)*napx
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
                  smiv(ncrr)=smi(ncrr)
                else if(ix.eq.iq(2).or.iratioe(ix).eq.iq(2)) then
                  smiv(ncrr)=smi(ncrr)
                endif
              enddo
            endif

            do 170 i=1,napx
              iar=(ib-1)*napx+i
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
          iar=(ib-1)*napx+1

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
            iar=ib1+(ib-1)*napx

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
            nms(ib2)=1

            if(ithick.eq.1) then
              do i=1,nele
                ekv(ib2,i)=ek(i)
              end do
            end if

          end do

          ib0=ib0+napx
  250   continue
#ifdef DEBUG
!     call dumpbin('ado 260',260,260)
!     call abend('ado 260                                           ')
#endif
  napx = napx*imc

#ifdef FLUKA
  ! A.Mereghetti, P. Garcia Ortega, D.Sinuela Pastor, V. Vlachoudis for the FLUKA Team
  ! last modified: 11-06-2014
  ! start connection to FLUKA and initialise max ID
  ! inserted in main code by the 'fluka' compilation flag
  if(fluka_enable) then
    fluka_con = fluka_is_running()
    if(fluka_con == -1) then
      write(lout,"(a)") "FLUKA> ERROR Fluka is expected to run but it is NOT actually the case"
      write(fluka_log_unit,*) "# Fluka is expected to run but it is NOT actually the case"
      call prror(-1)
    end if
    write(lout,"(a)") "FLUKA> Initializing FlukaIO interface ..."
    write(fluka_log_unit,*) "# Initializing FlukaIO interface ..."
    fluka_con = fluka_connect()
    if(fluka_con == -1) then
      write(lout,"(a)") "FLUKA> ERROR Cannot connect to Fluka server"
      write(fluka_log_unit,*) "# Error connecting to Fluka server"
      call prror(-1)
    endif
    write(lout,"(a)") "FLUKA> Successfully connected to Fluka server"
    write(fluka_log_unit,*) "# Successfully connected to Fluka server"
    fluka_connected = .true.
  endif
#endif

#ifdef CR
  write(93,"(a,i0)") "MAINCR> Setting napxo = ",napx
  flush(93)
#endif
  napxo = napx

  if(idp == 0 .or. ition == 0) then
    ! 4D tracking
    if(ithick == 1) then
      call meta_write("TrackingMethod", "Thick 4D")
    else
      call meta_write("TrackingMethod", "Thin 4D")
    end if
    if(iclo6 /= 0) then
      write(lout,"(a,i0)") "MAINCR> ERROR Doing 4D tracking but iclo6 = ",iclo6
      write(lout,"(a)")    "MAINCR>       Expected iclo6 = 0 for 4D tracking."
      call prror(-1)
    end if
  else
    ! 6D tracking
    if(ithick == 1) then
      call meta_write("TrackingMethod", "Thick 6D")
    else
      call meta_write("TrackingMethod", "Thin 6D")
    end if
    if(iclo6 == 0) then
      write(lout,"(a,i0)") "MAINCR> ERROR Doing 6D tracking but iclo6 = ",iclo6
      write(lout,"(a)")    "MAINCR>       Expected iclo6 <> 0 for 6D tracking."
      call prror(-1)
    end if
  end if

  call time_timeStamp(time_afterClosedOrbit)
  call meta_write("4D_ClosedOrbitCorr_x",  clo(1))
  call meta_write("4D_ClosedOrbitCorr_xp", clop(1))
  call meta_write("4D_ClosedOrbitCorr_y",  clo(2))
  call meta_write("4D_ClosedOrbitCorr_yp", clop(2))
  if(iclo6 /= 0) then
    call meta_write("6D_ClosedOrbitCorr_x",     clo6(1))
    call meta_write("6D_ClosedOrbitCorr_xp",    clop6(1))
    call meta_write("6D_ClosedOrbitCorr_y",     clo6(2))
    call meta_write("6D_ClosedOrbitCorr_yp",    clop6(2))
    call meta_write("6D_ClosedOrbitCorr_sigma", clo6(3))
    call meta_write("6D_ClosedOrbitCorr_dp",    clop6(3))
  end if

! ---------------------------------------------------------------------------- !
!  GENERATE THE INITIAL DISTRIBUTION
! ---------------------------------------------------------------------------- !

  do i=1,npart
    pstop(i)  = .false.
    nnumxv(i) = numl
    numxv(i)  = numl
  end do
  rat0 = rat

  ! DIST Block
  if(dist_enable) then
    e0f=sqrt(e0**2-nucm0**2)
    call dist_readDist
    call dist_finaliseDist
    call part_applyClosedOrbit
    if(dist_echo) call dist_echoDist
  end if

  if(idfor /= 2 .and. .not.dist_enable) then
    ! Generated from INIT Distribution Block
    do ia=1,napx,2
      if(st_quiet == 0) write(lout,10050)
      tasia56 = tas(ia,5,6)*c1m3
      bet0x2  = tas(ia,1,3)**2+tas(ia,1,4)**2
      bet0z2  = tas(ia,3,1)**2+tas(ia,3,2)**2
      bet0s1  = tas(ia,5,5)**2+tasia56**2
      dsign   = one
      rat     = rat0
      if(tas(ia,3,3) < (-one*pieni)) rat = -one*rat
      if(rat < (-one*pieni)) dsign = -one*one
      x11    = ampv(ia)/(sqrt(bet0(1))+sqrt(abs(rat)*bet0x2))
      x13    = (x11*dsign)*sqrt(abs(rat))
      amp(2) = (dsign*real(1-iver,fPrec))*(abs(x11)*sqrt(bet0z2)+abs(x13)*sqrt(bet0(2)))
      x1(5)  = zero
      x1(6)  = dpsv(ia)*sqrt(bet0s1)
      chi    = chi0*rad
      dchi   = chid*rad
      do i2=1,2
        i3    = ia+i2-1
        sic   = sin_mb(chi)
        coc   = cos_mb(chi)
        x1(1) = x11*coc
        x1(2) = x11*sic
        x1(3) = x13*coc
        x1(4) = x13*sic
        do ii=1,6
          x2(ii) = zero
          do jj=1,6
            x2(ii) = x2(ii)+tas(ia,ii,jj)*x1(jj)
          end do
        end do
        if(iclo6 == 1 .or. iclo6 == 2) then
          x2(2) = x2(2)/((one+x2(6))+clop6v(3,ia))
          x2(4) = x2(4)/((one+x2(6))+clop6v(3,ia))
        end if
        if(abs(bet0s1) <= pieni) x2(6) = dpsv(ia)
        if(iver == 1) then
          x2(3) = zero
          x2(4) = zero
        end if
        xv1(i3)   = x2(1)+exz(i2,1)
        yv1(i3)   = x2(2)+exz(i2,2)
        xv2(i3)   = x2(3)+exz(i2,3)
        yv2(i3)   = x2(4)+exz(i2,4)
        sigmv(i3) = x2(5)+exz(i2,5)
        dpsv(i3)  = x2(6)
        dpsic     = dpsv(i3)+clop6v(3,ia)
        if(idp == 1 .and. abs(ition) == 1 .and. iclo6 == 0) then
          xv1(i3) = xv1(i3) + di0xs(ia)*dpsic
          xv2(i3) = xv2(i3) + di0zs(ia)*dpsic
          yv1(i3) = yv1(i3) + dip0xs(ia)*dpsic
          yv2(i3) = yv2(i3) + dip0zs(ia)*dpsic
        end if
        chi = chi+dchi
      end do

      epsa(1)    = (ampv(ia)**2/bet0(1))
      epsa(2)    = (amp(2)**2/bet0(2))
      nucm(ia)   = nucm0
      nucm(ia+1) = nucm0

      if(st_quiet == 0) then
        write(lout,10260) ia,nms(ia)*izu0,dpsv(ia)
        write(lout,10060) xv1(ia),yv1(ia),xv2(ia),yv2(ia),sigmv(ia),dpsv(ia), &
          xv1(ia+1),yv1(ia+1),xv2(ia+1),yv2(ia+1),sigmv(ia+1),dpsv(ia+1)
        write(lout,10020) ampv(ia),amp(2),epsa
      end if
    end do
    call part_applyClosedOrbit

  else if(idfor == 2) then
    ! Read from fort.13
    call readFort13
    call part_updatePartEnergy(1)
    ! Note that this effectively overrides the particle delta set in fort.13
  endif

  do ia=1,napx,2
    if(.not.dist_enable .and. st_quiet == 0) then
      write(lout,10090) xv1(ia),yv1(ia),xv2(ia),yv2(ia),sigmv(ia),dpsv(ia),xv1(ia+1),&
        yv1(ia+1),xv2(ia+1),yv2(ia+1),sigmv(ia+1),dpsv(ia+1),e0,ejv(ia),ejv(ia+1)
    end if
    idam  = 3
    icode = 0
    if(abs(xv1(ia)) <= pieni .and. abs(yv1(ia)) <= pieni) then
      idam  = idam-1
    else
      icode = icode+1
    endif
    if(abs(xv2(ia)) <= pieni .and. abs(yv2(ia)) <= pieni) then
      idam  = idam-1
    else
      icode = icode+2
    endif
    if(idp == 0 .or. abs(ition) == 0) then
      idam  = idam-1
    else
      icode = icode+4
    endif
    if(idam  <= 0) idam  = 1
    if(icode <= 0) icode = 1
    ia2 = (ia+1)/2
    if(ntwin /= 2) then
      if(mod(ia+1,2) == 0) then
        xau(1,1) = xv1(ia)
        xau(1,2) = yv1(ia)
        xau(1,3) = xv2(ia)
        xau(1,4) = yv2(ia)
        xau(1,5) = sigmv(ia)
        xau(1,6) = dpsv(ia)
        xau(2,1) = xv1(ia+1)
        xau(2,2) = yv1(ia+1)
        xau(2,3) = xv2(ia+1)
        xau(2,4) = yv2(ia+1)
        xau(2,5) = sigmv(ia+1)
        xau(2,6) = dpsv(ia+1)
        cloau(1) = clo6v(1,ia)
        cloau(2) = clop6v(1,ia)
        cloau(3) = clo6v(2,ia)
        cloau(4) = clop6v(2,ia)
        cloau(5) = clo6v(3,ia)
        cloau(6) = clop6v(3,ia)
        di0au(1) = di0xs(ia)
        di0au(2) = dip0xs(ia)
        di0au(3) = di0zs(ia)
        di0au(4) = dip0zs(ia)

        do ib2=1,6
          do ib3=1,6
            tau(ib2,ib3)=tasau(ia,ib2,ib3)
          end do
        end do

        call distance(xau,cloau,di0au,tau,dam1)
        dam(ia)   = dam1
        dam(ia+1) = dam1
      end if

      ! Write header of track output file(s) used by postprocessing for case ntwin /= 2
#ifndef STF
#ifdef CR
      if(.not.restart) then
#endif
        call writebin_header(ia,ia,91-ia2,ierro,cdate,ctime,progrm)
#ifdef CR
        flush(91-ia2)
        binrecs(ia2)=1
      endif
#endif
#else
#ifdef CR
      if(.not.restart) then
#endif
        call writebin_header(ia,ia,90,ierro,cdate,ctime,progrm)
#ifdef CR
        flush(90)
        binrecs(ia2)=1
      endif
#endif
#endif
    else !ELSE for "if(ntwin.ne.2)"
      ! Write header of track output file(s) used by postprocessing for case ntwin == 2
#ifndef STF
#ifdef CR
      if(.not.restart) then
#endif
        call writebin_header(ia,ia+1,91-ia2,ierro,cdate,ctime,progrm)
#ifdef CR
        flush(91-ia2)
        binrecs(ia2)=1
      endif
#endif
#else
#ifdef CR
      if(.not.restart) then
#endif
        call writebin_header(ia,ia+1,90,ierro,cdate,ctime,progrm)
#ifdef CR
        flush(90)
        binrecs(ia2)=1
      endif
#endif
#endif
    endif !ENDIF (ntwin.ne.2)
    if(ierro /= 0) then
      write(lout,"(a,i0)") "MAINCR> ERROR Problems writing to file #",91-ia2
      write(lout,"(a,i0)") "MAINCR> ERROR Code: ",ierro
      goto 520
    endif
  end do ! napx

#ifdef CR
  if(lhc /= 9) binrec = 1
  ! binrec:  The maximum number of reccords writen for all tracking data files. Thus crbinrecs(:) <= binrec
#endif

  call time_timeStamp(time_afterBeamDist)

! ---------------------------------------------------------------------------- !
!  END GENERATE THE INITIAL DISTRIBUTION
! ---------------------------------------------------------------------------- !

! ---------------------------------------------------------------------------- !
!  PRE-TRACKING INITIALISATION
! ---------------------------------------------------------------------------- !

  if(ithick == 1) then
    ! Compute matrices for linear tracking
    call envarsv(dpsv,moidpsv,rvv,ekv)
    if(idp == 0 .or. ition == 0) then ! Only in case of thck4d
      call blocksv
    end if
  end if

#ifdef FLUKA
  ! P.Garcia Ortega, A.Mereghetti and V.Vlachoudis, for the FLUKA Team
  ! last modified: 26-08-2014
  ! send napx to fluka
  if(fluka_enable) then
    write(lout,"(a,i0)") "FLUKA> Sending napx = ",napx
    write(fluka_log_unit,*) "# Sending napx: ", napx
    fluka_con = fluka_init_max_uid( napx )

    if(fluka_con < 0) then
      write(lout,"(a,i0,a)") "FLUKA> ERROR Failed to send napx ",napx," to fluka "
      write(fluka_log_unit, *) "# failed to send napx to fluka ",napx
      call prror(-1)
    end if

    write(lout,"(a)") "FLUKA> Sending napx successful"
    write(fluka_log_unit,*) "# Sending napx successful;"
    flush(lout)
    flush(fluka_log_unit)
  end if

  ! A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
  ! last modified: 18-01-2016
  ! initialise energy/momentum/rest mass of reference particle in mod_fluka
  !     and synch magnetic rigidity with Fluka (for the time being, consider
  !     only protons);
  if(fluka_enable) then
    write(lout,"(a)") "FLUKA> Updating the reference particle"
    write(fluka_log_unit,*) "# Updating ref particle"
    flush(lout)
    flush(fluka_log_unit)

    fluka_con = fluka_set_synch_part( e0, e0f, nucm0, aa0, zz0)

    if(fluka_con < 0) then
      write(lout,"(a)") "FLUKA> ERROR Failed to update the reference particle"
      write(fluka_log_unit,*) "# failed to update ref particle"
      call prror(-1)
    end if

    write(lout,"(a)") "FLUKA> Updating the reference particle successful"
    write(fluka_log_unit,*) "# Updating ref particle successful;"
    flush(lout)
    flush(fluka_log_unit)
  end if

#endif

  ! Initialise Modules
  call dump_initialise

  call time_timeStamp(time_afterInitialisation)

! ---------------------------------------------------------------------------- !
!  START OF TRACKING
! ---------------------------------------------------------------------------- !
  write(lout,10200)
  call part_setParticleID
  call part_writeState(0)

  time1=0.
  call time_timerCheck(time1)

  ! time1 is now pre-processing CPU
  ! note that this will be reset every restart as we redo pre-processing
  pretime=time1-time0
  part_isTracking = .true.
  if(ithick == 0) call trauthin(nthinerr)
  if(ithick == 1) call trauthck(nthinerr)

  time2=0.
  call time_timerCheck(time2)

  ! trtime is now the tracking time, BUT we must add other time for C/R
  trtime=time2-time1
#ifdef CR
  ! because now crpoint will write tracking time using time3 as a temp and crcheck/crstart will reset crtime3
  trtime=trtime+crtime3
#endif
  if(nthinerr == 3000) goto 520
  if(nthinerr == 3001) goto 460

  ! END OF LOOP OVER TURNS
  460 continue

  ! Set numx=nnuml (for writebin) NOT for LOST particles because all lost set nnuml=numl
  numx = nnuml
  id   = 0

#ifndef FLUKA
  napxto = 0

#ifdef CR
  if(.not.restart) then
    ! If restart is true , we haven't done any tracking and must be running from very last checkpoint
    write(93,"(a)")          "MAINCR> Very last call to WRITEBIN?"
    write(93,"(a,3(1x,i0))") "MAINCR> numlmax, nnuml, numl = ",numlmax,nnuml,numl
    flush(93)
    if(nnuml == numl) then
      ! We REALLY have finished (or all particles lost)
      ! When all lost, nthinerr=3001, we set nnuml=numl
      ! and make sure we do the last WRITEBIN
      write(93,"(a)") "MAINCR> Very last call to WRITEBIN"
      flush(93)
      call writebin(nthinerr)
      if(nthinerr == 3000) goto 520
    else
      ! I assume we are stopping because we have done nnuml turns which should be numlmax and do a writebin only if time
      write(93,"(a)")          "MAINCR> Very last call to WRITEBIN?"
      write(93,"(a,3(1x,i0))") "MAINCR> numlmax, nnuml, numl = ",numlmax,nnuml,numl
      flush(93)
      if(mod(nnuml,nwri) == 0) then
        write(93,"(a)") "MAINCR> Very last call to WRITEBIN"
        flush(93)
        call writebin(nthinerr)
        if(nthinerr == 3000) goto 520
      end if
    end if
    ! do the very last checkpoint
    call callcrp()
  end if
#else
  call writebin(nthinerr)
  if(nthinerr == 3000) goto 520
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
  call time_timeStamp(time_afterTracking)

  if(st_writefort12) then
    call writeFort12
  end if

  if(st_partsum .eqv. .false.) then
    write(lout,"(a)") "MAINCR> NOTE Particle summary report is disabled; either manually, or because npart > 64."
    write(lout,"(a)") "MAINCR>      This is controlled by the PARTICLESUMMARY flag in the SETTINGS block in fort.3."
    write(lout,"(a)") ""
    goto 470
  end if

  write(lout,"(a)") "    PARTICLE SUMMARY:"
  write(lout,"(a)") ""

  do ia=1,napxo,2
    ie=ia+1
    napxto = napxto+numxv(ia)+numxv(ie)

    if(pstop(ia).and.pstop(ie)) then !-- BOTH PARTICLES LOST
      write(lout,10000) ia,nms(ia)*izu0,dp0v(ia),numxv(ia),abs(xv1(ia)),aperv(ia,1),abs(xv2(ia)),aperv(ia,2)
      write(lout,10000) ie,nms(ia)*izu0,dp0v(ia),numxv(ie),abs(xv1(ie)),aperv(ie,1),abs(xv2(ie)),aperv(ie,2)
    end if

    if(.not.pstop(ia).and.pstop(ie)) then !-- SECOND PARTICLE LOST
      if(st_quiet == 0) then
        write(lout,10240) ia,nms(ia)*izu0,dp0v(ia),numxv(ia)
      else if(st_quiet == 1) then
        write(lout,10241) ia,nms(ia)*izu0,dp0v(ia),numxv(ia)
      end if
      write(lout,10000) ie,nms(ia)*izu0,dp0v(ia),numxv(ie),abs(xv1(ie)),aperv(ie,1),abs(xv2(ie)),aperv(ie,2)
    end if

    if(pstop(ia).and..not.pstop(ie)) then !-- FIRST PARTICLE LOST
      write(lout,10000) ia,nms(ia)*izu0,dp0v(ia),numxv(ia),abs(xv1(ia)),aperv(ia,1),abs(xv2(ia)),aperv(ia,2)
      if(st_quiet == 0) then
        write(lout,10240) ie,nms(ia)*izu0,dp0v(ia),numxv(ie)
      else if(st_quiet == 1) then
        write(lout,10241) ie,nms(ia)*izu0,dp0v(ia),numxv(ie)
      end if
    end if

    if(.not.pstop(ia).and..not.pstop(ie)) then !-- BOTH PARTICLES STABLE
      if(st_quiet == 0) then
        write(lout,10270) ia,ie,nms(ia)*izu0,dp0v(ia),numxv(ia)
      else if(st_quiet == 1) then
        write(lout,10271) ia,ie,nms(ia)*izu0,dp0v(ia),numxv(ia)
      end if
    end if

    if(st_quiet == 0) then
      write(lout,"(4x,f47.33)") xv1(ia),yv1(ia),xv2(ia),yv2(ia),sigmv(ia),dpsv(ia), &
        xv1(ie),yv1(ie),xv2(ie),yv2(ie),sigmv(ie),dpsv(ie),e0,ejv(ia),ejv(ie)
    end if
  end do

#else
  ! IFDEF FLUKA
  ! A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
  ! last modified: 17-07-2013
  ! print stable particles only
  write(lout,"(a)") ""
  write(lout,"(a)") str_divLine
  if(napxo > 0) then
    write(lout,"(a)") ""
    write(lout,10350) napxo
    write(lout,"(a)") ""
    write(lout,10360) 'ID', 'GEN', 'WEIGHT', 'X [m]', 'XP []', 'Y [m]', 'YP[]', 'PC [GeV]', 'DE [eV]', 'DT [s]'
    write(lout,"(a)") ""
    do ia=1,napxo
      if(.not.pstop(ia)) then
        write(lout,10370) fluka_uid(ia),fluka_gen(ia),fluka_weight(ia), &
          xv1(ia)*c1m3, yv1(ia)*c1m3, xv2(ia)*c1m3, yv2(ia)*c1m3, &
          ejfv(ia)*c1m3,(ejv(ia)-e0)*c1e6,-c1m3*(sigmv(ia)/clight)*(e0/e0f)
      end if
    end do
  end if
#endif

! ---------------------------------------------------------------------------- !
!  POSTPROCESSING (POSTPR)
! ---------------------------------------------------------------------------- !

470 continue
  ! and we need to open fort.10 unless already opened for BOINC
  call f_open(unit=10, file="fort.10", formatted=.true., mode="rw",err=fErr,recl=8195)
  call f_open(unit=110,file="fort.110",formatted=.false.,mode="w", err=fErr)

  ! Also dump the final state of the particle arrays
  call part_writeState(1)

#ifndef FLUKA
#ifndef STF
  iposc = 0
  if(ipos == 1) then ! Variable IPOS=1 -> postprocessing block present in fort.3
    do ia=1,napxo,2
      ia2=(ia+1)/2
      iposc=iposc+1
#ifndef CR
      call postpr(91-ia2) ! Postprocess file "fort.(91-ia2)"
#else
      write(93,"(a,i0)") "MAINCR> Calling POSTPR nnuml = ",nnuml
      flush(93)
      call postpr(91-ia2,nnuml)
#endif
    end do
    if(iposc >= 1) call sumpos
  end if ! END if(ipos.eq.1)
  goto 520 ! Done postprocessing

490 continue ! GOTO here if(napx.le.0.or.imc.le.0) (skipping tracking)
  if(ipos == 1) then
    ndafi2=ndafi
    do ia=1,ndafi2
      if(ia > ndafi) exit
#ifndef CR
      call postpr(91-ia)
#else
      write(93,"(a,i0)") "MAINCR> Calling POSTPR nnuml = ",nnuml
      flush(93)
      call postpr(91-ia,nnuml)
#endif
    end do
    if(ndafi >= 1) call sumpos
  end if
#else
  ! IFDEF STF
  iposc=0
  if(ipos == 1) then ! Variable IPOS=1 -> postprocessing block present in fort.3
    do ia=1,napxo,2
      iposc=iposc+1
#ifndef CR
      call postpr(ia) ! Postprocess particle ia (and ia+1 if ntwin=2)
#else
      write(93,"(a,i0)") "MAINCR> Calling POSTPR nnuml = ",nnuml
      flush(93)
      call postpr(ia,nnuml)
#endif
    end do
    if(iposc >= 1) call sumpos
  end if
  goto 520 ! Done postprocessing

490 continue ! GOTO here if(napx.le.0.or.imc.le.0) (skipping tracking)
  if(ipos == 1) then
    ndafi2=ndafi
    do ia=1,(2*ndafi2),2
      if(ia > ndafi) exit
#ifndef CR
      call postpr(ia)
#else
      write(93,"(a,i0)") "MAINCR> Calling POSTPR nnuml = ",nnuml
      flush(93)
      call postpr(ia,nnuml)
#endif
    end do
    if(ndafi >= 1) call sumpos
  end if
#endif

! ---------------------------------------------------------------------------- !
!  DONE POSTPROCESSING (POSTPR)
! ---------------------------------------------------------------------------- !

520 continue
  call time_timeStamp(time_afterPostProcessing)
  if(fma_flag) then
    write(lout,"(a)") "MAINCR> Calling FMA_POSTPR"
    call fma_postpr
    call time_timeStamp(time_afterFMA)
  endif
  ! HPLOTTING END
  if(ipos == 1 .and. (idis /= 0 .or. icow /= 0 .or. istw /= 0 .or. iffw /= 0)) then
    call igmeta(999,0)
    call hplend
  endif
#endif

#ifdef FLUKA
  ! A.Mereghetti, for the FLUKA Team
  ! last modified: 28-05-2014
  ! collect a couple of goto statements, sending code flow
  !   to different plotting points, which are not actually
  !   inserted
490 continue
520 continue
  call fluka_close
#endif
#ifdef FFIELD
  ! Modification by B.DALENA and T.PUGNAT
  call ffield_mod_end()
#endif

  time3=0.
  call time_timerCheck(time3)
  ! Note that crpoint no longer destroys time2
  posttime=time3-time2

  ! Make sure all files are flushed before we do stuff with them
  call f_flush

#ifdef HASHLIB
  ! HASH library. Must be before ZIPF
  call hash_fileSums
  call time_timeStamp(time_afterHASH)
#endif

  if(zipf_numfiles > 0) then
    call zipf_dozip
    call time_timeStamp(time_afterZIPF)
  endif

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
  write(lout,"(a,i8)")      "    Particle Turns:           ",meta_nPartTurn
  write(lout,"(a)")         ""
  write(lout,"(a)")         str_divLine

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

  call time_timeStamp(time_beforeExit)
  call time_finalise
  call meta_finalise

! ---------------------------------------------------------------------------- !
!  DONE MAINCR
! ---------------------------------------------------------------------------- !

#ifdef CR
  call abend('                                                  ')
#else
  call closeUnits
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
#ifdef FLUKA
!     A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!     last modified: 17-07-2013
!     print stable particles only (format directives)
10350 format(4X,I8,1X,'SURVIVING PARTICLES:')
10360 format(2(1X,A8),8(1X,A16))
10370 format(2(1X,I8),8(1X,1PE16.9))
#endif
end program maincr
