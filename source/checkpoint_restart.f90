! ================================================================================================ !
!  SixTrack Checkpoint/Restart Module
! ====================================
!  E. McIntosh
!  K.N. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2019-04-29
! ================================================================================================ !
module checkpoint_restart

  use floatPrecision
  use, intrinsic :: iso_fortran_env, only : int16, int32

  implicit none

  ! Checkpoint Files
  integer,           parameter     :: cr_nPoint = 2
  character(len=15), public,  save :: cr_pntFile(cr_nPoint)  = ["crpoint_pri.bin","crpoint_sec.bin"]
  integer,           public,  save :: cr_pntUnit(cr_nPoint)  = -1
  logical,           public,  save :: cr_pntExist(cr_nPoint) = .false.
  logical,           public,  save :: cr_pntRead(cr_nPoint)  = .false.

  ! Logging Files
  character(len=13), parameter     :: cr_errFile  = "cr_stderr.tmp"
  character(len=13), parameter     :: cr_outFile  = "cr_stdout.tmp"
  character(len=13), parameter     :: cr_logFile  = "cr_status.log"
  integer,           parameter     :: cr_errUnit  = 91
  integer,           parameter     :: cr_outUnit  = 92
  integer,           parameter     :: cr_logUnit  = 93

  ! Checkpoint/Restart Flags and Variables
  logical,           public,  save :: cr_rerun    = .false.
  logical,           public,  save :: cr_start    = .true.
  logical,           public,  save :: cr_restart  = .false.

  character(len=21), public,  save :: cr_startMsg = " "
  integer,           public,  save :: cr_numl     = 1
  integer,           public,  save :: binrec      = 0       ! The maximum number of records writen for all tracking data files
  integer,           public,  save :: sixrecs     = 0

  ! Keep length in sync with version.f90
  character(len=8),  private, save :: cr_version  = " "
  character(len=10), private, save :: cr_moddate  = " "

  ! C/R Temp Variables and Arrays
  real(kind=fPrec),  private, save :: cre0
  real(kind=fPrec),  private, save :: crbeta0
  real(kind=fPrec),  private, save :: crbrho
  real(kind=fPrec),  private, save :: crnucmda

  integer,           private, save :: crnpart_old = -1
  integer,           private, save :: crsixrecs
  integer,           public,  save :: crbinrec
  integer,           private, save :: cril
  integer,           private, save :: crnumlcr
  integer,           private, save :: crnuml
  integer,           private, save :: crnapxo
  integer,           private, save :: crnapx

  real(kind=fPrec), allocatable, private, save :: crxv1(:)      ! (npart)
  real(kind=fPrec), allocatable, private, save :: crxv2(:)      ! (npart)
  real(kind=fPrec), allocatable, private, save :: cryv1(:)      ! (npart)
  real(kind=fPrec), allocatable, private, save :: cryv2(:)      ! (npart)
  real(kind=fPrec), allocatable, private, save :: crsigmv(:)    ! (npart)
  real(kind=fPrec), allocatable, private, save :: crdpsv(:)     ! (npart)
  real(kind=fPrec), allocatable, private, save :: crdpsv1(:)    ! (npart)
  real(kind=fPrec), allocatable, private, save :: crejv(:)      ! (npart)
  real(kind=fPrec), allocatable, private, save :: crejfv(:)     ! (npart)
  real(kind=fPrec), allocatable, private, save :: craperv(:,:)  ! (npart,2)
  real(kind=fPrec), allocatable, private, save :: crnucm(:)     ! (npart)
  real(kind=fPrec), allocatable, private, save :: crmtc(:)      ! (npart)

  integer(kind=int16), allocatable, private, save :: crnaa(:)   ! (npart)
  integer(kind=int16), allocatable, private, save :: crnzz(:)   ! (npart)
  integer(kind=int16), allocatable, private, save :: crnqq(:)   ! (npart)
  integer(kind=int32), allocatable, private, save :: crpdgid(:) ! (npart)

  integer(kind=int32), allocatable, private, save :: crpartID(:)   ! (npart)
  integer(kind=int32), allocatable, private, save :: crparentID(:) ! (npart)

  integer,          allocatable, public,  save :: binrecs(:)    ! ((npart+1)/2)
  integer,          allocatable, public,  save :: crbinrecs(:)  ! (npart+1)/2)
  integer,          allocatable, private, save :: crnumxv(:)    ! (npart)
  integer,          allocatable, private, save :: crpairID(:,:) ! (2,npart)

  logical,          allocatable, private, save :: crpstop(:)    ! (npart)
  logical,          allocatable, private, save :: crllostp(:)   ! (npart)

contains

! ================================================================================================ !
!  INIT CRPOINT FILES
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2019-04-29
!  Updated: 2019-04-29
! ================================================================================================ !
subroutine cr_fileInit

  use, intrinsic :: iso_fortran_env, only : output_unit

  use mod_units
  use mod_common, only : fort6
#ifdef ZLIB
  use zipf
#endif

  integer i, zErr
  logical fErr
  character(len=256) fileName

10 continue

  ! First get rid of any previous partial output, and open log file for append (rw+)
  call f_close(cr_errUnit)
  call f_close(cr_outUnit)
  call f_close(cr_logUnit)
  call f_open(unit=cr_errUnit,file=cr_errFile,formatted=.true.,mode="rw", err=fErr,status="replace")
  call f_open(unit=cr_outUnit,file=cr_outFile,formatted=.true.,mode="rw", err=fErr,status="replace")
  call f_open(unit=cr_logUnit,file=cr_logFile,formatted=.true.,mode="rw+",err=fErr,status="unknown")

  ! Now we see if we have a fort.6, which implies that we can perhaps just restart using all exisiting files
  ! including the last checkpoints. If not, we just do a start (with an unzip for BOINC)
  call f_open(unit=output_unit,file=fort6,formatted=.true.,mode="rw",err=fErr,status="old")
  if(fErr) then
#ifdef BOINC
    ! No fort.6 so we do an unzip of Sixin.zip
    ! BUT ONLY IF WE HAVE NOT DONE IT ALREADY
    if(cr_start) then
      cr_start = .false.
      call f_close(cr_errUnit)
      call f_close(cr_outUnit)
      call f_close(cr_logUnit)
      ! Now, if BOINC, after no fort.6, call UNZIP Sixin.zip
      ! Name hardcoded in our boinc_unzip_.
      ! Either it is only the fort.* input data or it is a restart.
      call boincrf("Sixin.zip",fileName)
#ifdef ZLIB
      call minizip_unzip(trim(fileName),".",zErr,len_trim(fileName),1)
      if(zErr /= 0) then
        write(cr_errUnit,"(a)") "SIXTRACR> ERROR Could not extract 'Sixin.zip'"
        call prror
      end if
#else
      write(cr_errUnit,"(a)") "SIXTRACR> ERROR No library available to extract zipfile 'Sixin.zip'"
      call prror
#endif
      goto 10 ! Go to top and check everything again after unzip
    end if
    call f_open(unit=output_unit,file=fort6,formatted=.true.,mode="rw",err=fErr)
#else
    call f_open(unit=output_unit,file=fort6,formatted=.true.,mode="rw",err=fErr,status="new")
#endif
    cr_startMsg = "SIXTRACR> Starts on: "
  else
    cr_startMsg = "SIXTRACR> Reruns on: "
    cr_rerun = .true.
  end if

  ! Open checkpoint files
  do i=1,cr_nPoint
    fErr = .false.
    call f_requestUnit(cr_pntFile(i),cr_pntUnit(i))
    call f_open(unit=cr_pntUnit(i),file=cr_pntFile(i),formatted=.false.,mode="rw",err=fErr,status="old")
    if(fErr) then
      call f_open(unit=cr_pntUnit(i),file=cr_pntFile(i),formatted=.false.,mode="rw",err=fErr,status="new")
    else
      cr_pntExist(i) = .true.
    end if
  end do

end subroutine cr_fileInit

! ================================================================================================ !
!  CR ALLOCATE ARRAYS
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-06-22
! ================================================================================================ !
subroutine cr_expand_arrays(npart_new)

  use mod_alloc
  use numerical_constants, only : zero, one
  use mod_common, only : nucm0, aa0, zz0, qq0, pdgid0

  integer, intent(in) :: npart_new

  integer :: npair_new
  npair_new = npart_new/2 + 1

  call alloc(crxv1,        npart_new, zero,    "crxv1")
  call alloc(crxv2,        npart_new, zero,    "crxv2")
  call alloc(cryv1,        npart_new, zero,    "cryv1")
  call alloc(cryv2,        npart_new, zero,    "cryv2")
  call alloc(crsigmv,      npart_new, zero,    "crsigmv")
  call alloc(crdpsv,       npart_new, zero,    "crdpsv")
  call alloc(crdpsv1,      npart_new, zero,    "crdpsv1")
  call alloc(crejv,        npart_new, zero,    "crejv")
  call alloc(crejfv,       npart_new, zero,    "crejfv")
  call alloc(crnucm,       npart_new, nucm0,   "crnucm")
  call alloc(crmtc,        npart_new, one,     "crmtc")
  call alloc(crnaa,        npart_new, aa0,     "crnaa")
  call alloc(crnzz,        npart_new, zz0,     "crnzz")
  call alloc(crnqq,        npart_new, qq0,     "crnqq")
  call alloc(crpdgid,      npart_new, pdgid0,  "crpdgid")
  call alloc(craperv,   2, npart_new, zero,    "craperv")
  call alloc(binrecs,      npair_new, 0,       "binrecs")
  call alloc(crbinrecs,    npair_new, 0,       "crbinrecs")
  call alloc(crnumxv,      npart_new, 0,       "crnumxv")
  call alloc(crpartID,     npart_new, 0,       "crpartID")
  call alloc(crparentID,   npart_new, 0,       "crparentID")
  call alloc(crpairID,  2, npart_new, 0,       "crpairID")
  call alloc(crpstop,      npart_new, .false., "crpstop")
  call alloc(crllostp,     npart_new, .false., "crllostp")

  crnpart_old = npart_new

end subroutine cr_expand_arrays

! ================================================================================================ !
!  CR KILL SWITCH
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-16
!  This routine will kill SixTrack if the current turn number matches a number in crkillturns
! ================================================================================================ !
subroutine cr_killSwitch(iTurn)

  use crcoall
  use mod_units
  use mod_settings

  integer, intent(in) :: iTurn

  logical killIt, fExist, onKillTurn
  integer pTurn, nKills, i, cUnit, tUnit

  killIt = .false.
  onKillTurn = .false.

  do i=1,size(st_killturns,1)
    if(iTurn == st_killturns(i)) then
      onKillTurn = .true.
    end if
  end do
  if(onKillTurn .eqv. .false.) then
    return
  end if

  call f_requestUnit("crkillswitch.tmp",cUnit)
  call f_requestUnit("crrestartme.tmp", tUnit)

  inquire(file="crkillswitch.tmp",exist=fExist)
  if(fExist .eqv. .false.) then
    call f_open(unit=cUnit,file="crkillswitch.tmp",formatted=.false.,mode="w",access="stream",status="replace")
    write(cUnit) 0,0
    flush(cUnit)
    call f_close(cUnit)
  end if

  call f_open(unit=cUnit,file="crkillswitch.tmp",formatted=.false.,mode="r",access="stream",status="old")
  read(cUnit) pTurn,nKills
  flush(cUnit)
  call f_close(cUnit)
  if(pTurn > 0) then
    write(lout, "(a,i0)") "CRKILLSW> Kill switch previously triggered on turn ",pTurn
    write(crlog,"(a,i0)") "CRKILLSW> Kill switch previously triggered on turn ",pTurn
    flush(lout)
    flush(crlog)
  end if

  do i=1,size(st_killturns,1)
    if(iTurn == st_killturns(i) .and. iTurn > pTurn) then
      killIt = .true.
      exit
    end if
  end do

  if(killIt) then
    nKills = nKills + 1

    write(lout, "(a,i0)") "CRKILLSW> Triggering kill switch on turn ",iTurn
    write(crlog,"(a,i0)") "CRKILLSW> Triggering kill switch on turn ",iTurn
    flush(lout)
    flush(crlog)

    call f_open(unit=tUnit,file="crrestartme.tmp",formatted=.false.,mode="w",access="stream",status="replace")
    write(tUnit) 1
    flush(tUnit)
    call f_close(tUnit)

    call f_open(unit=cUnit,file="crkillswitch.tmp",formatted=.false.,mode="w",access="stream",status="replace")
    write(cUnit) iTurn,nKills
    flush(cUnit)
    call f_close(cUnit)
    stop
  end if

end subroutine cr_killSwitch

! ================================================================================================ !
!  Check CR Files
! ================
!  Last modified: 2018-12-05
!
!  This subroutine checks if the C/R files exist, and if so tries to load them into the cr* variables.
!  This routine also repositions the output files for singletrackfile.dat and various other modules
! ================================================================================================ !
subroutine crcheck

  use, intrinsic :: iso_fortran_env, only : output_unit

  use parpro
  use crcoall
  use mod_common
  use mod_commons
  use mod_common_main
  use mod_version

  use dynk,        only : dynk_enabled,dynk_crcheck_readdata,dynk_crcheck_positionFiles
  use dump,        only : dump_crcheck_readdata,dump_crcheck_positionFiles
  use aperture,    only : limifound, aper_crcheck_readdata, aper_crcheck_positionFiles
  use scatter,     only : scatter_active, scatter_crcheck_readdata, scatter_crcheck_positionFiles
  use elens,       only : nelens, elens_crcheck
  use mod_meta,    only : meta_crcheck
  use mod_time,    only : time_crcheck
  use mod_random,  only : rnd_crcheck
  use collimation, only : coll_crcheck_readdata, coll_crcheck_positionFiles, do_coll
  use coll_db    , only : coll_db_crcheck_readdata

  integer j,k,l,m
  integer nPoint, ioStat

  logical noRestart, rErr
  character(len=1024) inLine

  ! Check the size of CR arrays
  if(npart /= crnpart_old) call cr_expand_arrays(npart)

  cr_restart    = .false.
  cr_pntRead(:) = .false.

  ! Some log entries to fort.93
  write(crlog,"(2(a,l1))") "CR_CHECK> Called with restart = ",cr_restart,", rerun = ",cr_rerun
  flush(crlog)

  noRestart = .false.
  do nPoint=1,cr_nPoint
    noRestart = noRestart .and. .not. cr_pntExist(nPoint)
  end do
  if(noRestart) goto 200

  ! If we do we must have a fort.6 as it was created by CRPOINT
  ! NOT TRUE anymore??? We might be NOT rerun but using a Sixin.zip
!#ifndef BOINC
!  This code is broken
!  if(cr_rerun .eqv. .false.) then
!    write(lerr,"(a)") "CR_CHECK> ERROR Found checkpoint file(s) but no "//trim(fort6)
!    call prror
!  end if
!#endif

  ! Check at least one restart file is readable
  noRestart = .true.
  do nPoint=1,cr_nPoint
    write(crlog,"(a)") "CR_CHECK> Checking file "//cr_pntFile(nPoint)
    flush(crlog)

    if(cr_pntExist(nPoint) .eqv. .false.) then
      write(crlog,"(a)") "CR_CHECK> File "//cr_pntFile(nPoint)//" does not exist"
      cycle
    end if

    rewind(cr_pntUnit(nPoint))

    write(crlog,"(a)") "CR_CHECK>  * Checking header"
    flush(crlog)
    read(cr_pntUnit(nPoint),iostat=ioStat) cr_version,cr_moddate
    if(ioStat /= 0) then
      write(crlog,"(a)") "CR_CHECK> Corrupt or truncated checkpoint file."
      flush(crlog)
      cycle
    end if
    if(cr_version == " " .or. cr_moddate == " ") then
      write(crlog,"(a)") "CR_CHECK> Unknown SixTrack version. Skipping this file."
      flush(crlog)
      cycle
    end if
    if(cr_version /= version .or. cr_moddate /= moddate) then
      write(crlog,"(a)") "CR_CHECK> Checkpoint files "//cr_pntFile(nPoint)//" was written by SixTrack version "//&
        trim(cr_version)//" with release date "//trim(cr_moddate)
      write(crlog,"(a)") "CR_CHECK> This is SixTrack version "//trim(version)//" with release date "//trim(moddate)
      write(crlog,"(a)") "CR_CHECK> Version mismatch. Skipping this file."
      flush(crlog)
      cycle
    end if

    write(crlog,"(a)") "CR_CHECK>  * Tracking variables"
    flush(crlog)
    read(cr_pntUnit(nPoint),iostat=ioStat) crnumlcr,crnuml,crsixrecs,crbinrec,cril,crnapxo, &
      crnapx,cre0,crbeta0,crbrho,crnucmda
    if(ioStat /= 0) cycle

    write(crlog,"(a)") "CR_CHECK>  * Particle arrays"
    flush(crlog)
    read(cr_pntUnit(nPoint),iostat=ioStat) &
      (crbinrecs(j), j=1,(crnapxo+1)/2),   &
      (crnumxv(j),   j=1,crnapxo),         &
      (crpartID(j),  j=1,crnapxo),         &
      (crparentID(j),j=1,crnapxo),         &
      (crpairID(:,j),j=1,crnapxo),         &
      (crpstop(j),   j=1,crnapxo),         &
      (crxv1(j),     j=1,crnapxo),         &
      (cryv1(j),     j=1,crnapxo),         &
      (crxv2(j),     j=1,crnapxo),         &
      (cryv2(j),     j=1,crnapxo),         &
      (crsigmv(j),   j=1,crnapxo),         &
      (crdpsv(j),    j=1,crnapxo),         &
      (crdpsv1(j),   j=1,crnapxo),         &
      (crejv(j),     j=1,crnapxo),         &
      (crejfv(j),    j=1,crnapxo),         &
      (crnucm(j),    j=1,crnapxo),         &
      (crmtc(j),     j=1,crnapxo),         &
      (crnaa(j),     j=1,crnapxo),         &
      (crnzz(j),     j=1,crnapxo),         &
      (crnqq(j),     j=1,crnapxo),         &
      (crpdgid(j),   j=1,crnapxo),         &
      (craperv(:,j), j=1,crnapxo),         &
      (crllostp(j),  j=1,crnapxo)
    if(ioStat /= 0) cycle

    write(crlog,"(a)") "CR_CHECK>  * META variables"
    flush(crlog)
    call meta_crcheck(cr_pntUnit(nPoint),rErr)
    if(rErr) cycle

    write(crlog,"(a)") "CR_CHECK>  * TIME variables"
    flush(crlog)
    call time_crcheck(cr_pntUnit(nPoint),rErr)
    if(rErr) cycle

    write(crlog,"(a)") "CR_CHECK>  * RND variables"
    flush(crlog)
    call rnd_crcheck(cr_pntUnit(nPoint),rErr)
    if(rErr) cycle

    write(crlog,"(a)") "CR_CHECK>  * DUMP variables"
    flush(crlog)
    call dump_crcheck_readdata(cr_pntUnit(nPoint),rErr)
    if(rErr) cycle

    if(dynk_enabled) then
      write(crlog,"(a)") "CR_CHECK>  * DYNK variables"
      flush(crlog)
      call dynk_crcheck_readdata(cr_pntUnit(nPoint),rErr)
      if(rErr) cycle
    end if

    if(scatter_active) then
      write(crlog,"(a)") "CR_CHECK>  * SCATTER variables"
      flush(crlog)
      call scatter_crcheck_readdata(cr_pntUnit(nPoint),rErr)
      if(rErr) cycle
    end if

    if(limifound) then
      write(crlog,"(a)") "CR_CHECK>  * APERTURE variables"
      flush(crlog)
      call aper_crcheck_readdata(cr_pntUnit(nPoint),rErr)
      if(rErr) cycle
    end if

    if(nelens > 0) then
      write(crlog,"(a)") "CR_CHECK>  * ELENS variables"
      flush(crlog)
      call elens_crcheck(cr_pntUnit(nPoint),rErr)
      if(rErr) cycle
    end if

    if(do_coll) then
      write(crlog,"(a)") "CR_CHECK>  * COLLIMATION variables"
      flush(crlog)
      call coll_db_crcheck_readdata(cr_pntUnit(nPoint),rErr)
      call coll_crcheck_readdata(cr_pntUnit(nPoint),rErr)
      if(rErr) cycle
    end if

    write(crlog,"(a)") "CR_CHECK> File "//cr_pntFile(nPoint)//" successfully read"
    flush(crlog)

    cr_pntRead(nPoint) = .true.
    noRestart = .false.
    exit
  end do

  if(noRestart) then
    write(crlog,"(a)") "CR_CHECK> ERROR No complete checkpoint file found"
    flush(crlog)
    goto 200
  end if

  ! If we have successfully read one of the checkpoint files, we need to handle lost particles and ntwin /= 2
  ! Otherwise we just continue with checkpointing as requested
  write(crlog,"(2(a,i8))") "CR_CHECK> Particles  C/R: ",crnapxo,  ", Input:  ",napx*2
  write(crlog,"(2(a,i8))") "CR_CHECK> SixRecords C/R: ",crsixrecs,", Buffer: ",sixrecs
  write(crlog,"(1(a,i8))") "CR_CHECK> BinRecords C/R: ",crbinrec
  flush(crlog)

  !  Position Files
  ! ================

  ! Position fort.6 to last checkpoint
  do j=1,crsixrecs
    read(output_unit,"(a1024)",end=110,err=120,iostat=ioStat) inLine
    sixrecs = sixrecs + 1
  end do
  ! This is not a FLUSH!
  endfile(output_unit,iostat=ierro)
110 continue
  backspace(output_unit,iostat=ierro)
  write(crlog,"(a,i0,a)") "CR_CHECK> Found ",sixrecs," lines in "//trim(fort6)
  flush(crlog)

  call cr_positionTrackFiles

  if(dynk_enabled) then
    write(crlog,"(a)") "CR_CHECK> Repositioning DYNK files"
    flush(crlog)
    call dynk_crcheck_positionFiles
  end if

  write(crlog,"(a)") "CR_CHECK> Repositioning DUMP files"
  flush(crlog)
  call dump_crcheck_positionFiles

  if(scatter_active) then
    write(crlog,"(a)") "CR_CHECK> Repositioning SCATTER files"
    flush(crlog)
    call scatter_crcheck_positionFiles
  end if

  if(limifound) then
    write(crlog,"(a)") "CR_CHECK> Repositioning APERTURE files"
    flush(crlog)
    call aper_crcheck_positionFiles
  end if

  if(do_coll) then
    write(crlog,"(a)") "CR_CHECK> Repositioning COLLIMATION files"
    flush(crlog)
    call coll_crcheck_positionFiles
  end if

  ! Set up flag for tracking routines to call CRSTART
  cr_restart = .true.

  return

  ! We are not checkpointing or we have no checkpoints, or we have no readable checkpoint
  ! If not checkpointing we can just give up on lout and use fort.6. We don't need to count records at all
200 continue
  write(crlog,"(a)") "SIXTRACR> No restart possible"
  flush(crlog)

  return

120 continue
  write(lerr, "(a)") "CR_CHECK> ERROR Cannot read "//trim(fort6)
  write(crlog,"(a)") "CR_CHECK> ERROR Cannot read "//trim(fort6)
  flush(crlog)
  call prror

end subroutine crcheck

! ================================================================================================ !
!  Write Checkpoint
! ==================
!  Last modified: 2019-04-29
!  This subroutine writes the checkpoint data to the binary checpoint files
! ================================================================================================ !
subroutine crpoint

  use crcoall
  use mod_time
  use mod_common
  use mod_commons
  use mod_common_main
  use mod_common_track
  use mod_version
  use mod_settings
  use numerical_constants
  use mod_units,   only : f_flush

  use dynk,        only : dynk_enabled,dynk_getvalue,dynk_fSets_cr,dynk_cSets_unique,dynk_nSets_unique,dynk_crpoint
  use dump,        only : dump_crpoint
  use aperture,    only : aper_crpoint,limifound
  use scatter,     only : scatter_active, scatter_crpoint
  use elens,       only : nelens, elens_crpoint
  use mod_meta,    only : meta_crpoint
  use mod_random,  only : rnd_crpoint
  use collimation, only : coll_crpoint, do_coll
  use coll_db,     only : coll_db_crpoint

  integer j, k, l, m, nPoint
  logical wErr, fErr

! flush everything to be safe
! This is to ensure that the logged file positions match what has been called in
! each case by write()
  call f_flush()

  if(numx >= numl) then
    write(crlog,"(a)") "CR_POINT> Called after last turn"
  else
    write(crlog,"(3(a,i0))") "CR_POINT> Called on turn ",(numx+1)," / ",numl," : interval is ",numlcp
  end if
  flush(crlog)

  if(cr_restart) then
    cr_restart = .false.
    return
  end if

  call time_startClock(time_clockCR)

  ! Copy lout to output_unit
  call cr_copyOut

  crnumlcr = numx+1

  if(dynk_enabled) then ! Store current settings of elements affected by DYNK
    if(st_debug) then
      write(crlog,"(a)") "CR_POINT> Filling DYNK sets"
      flush(crlog)
    end if
    do j=1,dynk_nSets_unique
      dynk_fSets_cr(j) = dynk_getvalue(dynk_cSets_unique(j,1),dynk_cSets_unique(j,2))
    end do
  end if

  !  Write the CR files
  ! ====================

  do nPoint=1,cr_nPoint

    wErr = .false.
    rewind(cr_pntUnit(nPoint))

    if(st_debug) then
      write(crlog,"(a)") "CR_POINT> Writing to checkpoint file "//cr_pntFile(nPoint)
      write(crlog,"(a)") "CR_POINT>  * SixTrack version"
      flush(crlog)
    end if
    write(cr_pntUnit(nPoint),err=100) version, moddate

    if(st_debug) then
      write(crlog,"(a)") "CR_POINT>  * Tracking variables"
      flush(crlog)
    end if
    write(cr_pntUnit(nPoint),err=100) crnumlcr,numl,sixrecs,binrec,il,napxo,napx,e0,beta0,brho,nucmda

    if(st_debug) then
      write(crlog,"(a)") "CR_POINT>  * Particle arrays"
      flush(crlog)
    end if
    write(cr_pntUnit(nPoint),err=100) &
      (binrecs(j), j=1,(napxo+1)/2),  &
      (numxv(j),   j=1,napxo),        &
      (partID(j),  j=1,napxo),        &
      (parentID(j),j=1,napxo),        &
      (pairID(:,j),j=1,napxo),        &
      (pstop(j),   j=1,napxo),        &
      (xv1(j),     j=1,napxo),        &
      (yv1(j),     j=1,napxo),        &
      (xv2(j),     j=1,napxo),        &
      (yv2(j),     j=1,napxo),        &
      (sigmv(j),   j=1,napxo),        &
      (dpsv(j),    j=1,napxo),        &
      (dpsv1(j),   j=1,napxo),        &
      (ejv(j),     j=1,napxo),        &
      (ejfv(j),    j=1,napxo),        &
      (nucm(j),    j=1,napxo),        &
      (mtc(j),     j=1,napxo),        &
      (naa(j),     j=1,napxo),        &
      (nzz(j),     j=1,napxo),        &
      (nqq(j),     j=1,napxo),        &
      (pdgid(j),   j=1,napxo),        &
      (aperv(:,j), j=1,napxo),        &
      (llostp(j),  j=1,napxo)
    flush(cr_pntUnit(nPoint))

    if(st_debug) then
      write(crlog,"(a)") "CR_POINT>  * META variables"
      flush(crlog)
    end if
    call meta_crpoint(cr_pntUnit(nPoint),wErr)
    if(wErr) goto 100

    if(st_debug) then
      write(crlog,"(a)") "CR_POINT>  * TIME variables"
      flush(crlog)
    end if
    call time_crpoint(cr_pntUnit(nPoint),wErr)
    if(wErr) goto 100

    if(st_debug) then
      write(crlog,"(a)") "CR_POINT>  * RND variables"
      flush(crlog)
    end if
    call rnd_crpoint(cr_pntUnit(nPoint),wErr)
    if(wErr) goto 100

    if(st_debug) then
      write(crlog,"(a)") "CR_POINT>  * DUMP variables"
      flush(crlog)
    end if
    call dump_crpoint(cr_pntUnit(nPoint),wErr)
    if(wErr) goto 100

    if(dynk_enabled) then
      if(st_debug) then
        write(crlog,"(a)") "CR_POINT>  * DYNK variables"
        flush(crlog)
      end if
      call dynk_crpoint(cr_pntUnit(nPoint),wErr)
      if(wErr) goto 100
    end if

    if(scatter_active) then
      if(st_debug) then
        write(crlog,"(a)") "CR_POINT>  * SCATTER variables"
        flush(crlog)
      end if
      call scatter_crpoint(cr_pntUnit(nPoint),wErr)
      if(wErr) goto 100
    end if

    if(limifound) then
      if(st_debug) then
        write(crlog,"(a)") "CR_POINT>  * APERTURE variables"
        flush(crlog)
      end if
      call aper_crpoint(cr_pntUnit(nPoint),wErr)
      if(wErr) goto 100
    end if

    if(nelens > 0) then
      if(st_debug) then
        write(crlog,"(a)") "CR_POINT>  * ELENS variables"
        flush(crlog)
      end if
      call elens_crpoint(cr_pntUnit(nPoint),wErr)
      if(wErr) goto 100
    end if

    if(do_coll) then
      if(st_debug) then
        write(crlog,"(a)") "CR_POINT>  * COLLIMATION variables"
        flush(crlog)
      end if
      call coll_db_crpoint(cr_pntUnit(nPoint),wErr)
      call coll_crpoint(cr_pntUnit(nPoint),wErr)
      if(wErr) goto 100
    end if

    flush(crlog)
    flush(cr_pntUnit(nPoint))

  end do ! Loop over nPoint

  call time_stopClock(time_clockCR)

  return

100 continue
  write(lerr ,"(a,i0)") "CR_POINT> ERROR Writing checkpoint file. iostat = ",ierro
  write(crlog,"(a,i0)") "CR_POINT> ERROR Writing checkpoint file. iostat = ",ierro
  flush(crlog)
  call prror

end subroutine crpoint

! ================================================================================================ !
!  Checkpoint Start
! ==================
!  Last modified: 2019-04-29
!  If we are restarting (cr_restart is TRUE), this routine is called in the beginning of the tracking
!  loops. It is used to copy the cr* variables to the normal variables, e.g. crnapx -> napx etc.
! ================================================================================================ !
subroutine crstart

  use parpro
  use crcoall
  use mod_units
  use mod_common
  use mod_commons
  use mod_common_main
  use mod_common_track
  use numerical_constants

  use dynk,        only : dynk_enabled, dynk_crstart
  use scatter,     only : scatter_active, scatter_crstart
  use elens,       only : nelens, elens_crstart
  use mod_meta,    only : meta_crstart
  use mod_time,    only : time_crstart
  use mod_random,  only : rnd_crstart
  use collimation, only : coll_crstart, do_coll

  logical fErr
  integer j, k, l, m, nPoint, ioStat

  write(crlog,"(a,i0)") "CR_START> Starting from checkpoint data from turn ",crnumlcr
  flush(crlog)

  ! We do NOT reset numl so that a run can be extended
  ! for more turns from the last checkpoint
  cr_numl = crnumlcr

  binrec = crbinrec
  napxo  = crnapxo
  napx   = crnapx
  e0     = cre0
  e0f    = sqrt(e0**2-nucm0**2)
  beta0  = crbeta0
  brho   = crbrho
  nucmda = crnucmda

  write(crlog,"(a)") "CR_START> Loading BinRecords"
  do j=1,(napxo+1)/2
    binrecs(j) = crbinrecs(j)
  end do

  write(crlog,"(a)") "CR_START> Loading tracking data"
  flush(crlog)

  partID(1:napxo)   = crpartID(1:napxo)
  parentID(1:napxo) = crparentID(1:napxo)
  pairID(:,1:napxo) = crpairID(:,1:napxo)
  pstop(1:napxo)    = crpstop(1:napxo)
  llostp(1:napxo)   = crllostp(1:napxo)

  xv1(1:napxo)      = crxv1(1:napxo)
  yv1(1:napxo)      = cryv1(1:napxo)
  xv2(1:napxo)      = crxv2(1:napxo)
  yv2(1:napxo)      = cryv2(1:napxo)
  sigmv(1:napxo)    = crsigmv(1:napxo)

  dpsv(1:napxo)     = crdpsv(1:napxo)
  dpsv1(1:napxo)    = crdpsv1(1:napxo)
  ejv(1:napxo)      = crejv(1:napxo)
  ejfv(1:napxo)     = crejfv(1:napxo)

  nucm(1:napxo)     = crnucm(1:napxo)
  mtc(1:napxo)      = crmtc(1:napxo)
  naa(1:napxo)      = crnaa(1:napxo)
  nzz(1:napxo)      = crnzz(1:napxo)
  nqq(1:napxo)      = crnqq(1:napxo)
  pdgid(1:napxo)    = crpdgid(1:napxo)

  numxv(1:napxo)    = crnumxv(1:napxo)
  do j=1,napxo
    if(pstop(j) .eqv. .false.) then
      numxv(j) = numl
    end if
  end do

  ! Recompute from loaded arrays (keep in sync with mod_particles)
  oidpsv(1:napxo)   = one/(one + dpsv(1:napxo))
  moidpsv(1:napxo)  = mtc(1:napxo)/(one + dpsv(1:napxo))
  omoidpsv(1:napxo) = ((one-mtc(1:napxo))*oidpsv(1:napxo))*c1e3
  rvv(1:napxo)      = (ejv(1:napxo)*e0f)/(e0*ejfv(1:napxo))

  ! Recompute the thick arrays
  if(ithick == 1) call synuthck

  ! Recompute the map of particle pairs
  call updatePairMap

  ! Aperture data
  aperv(1:2,1:napxo) = craperv(1:2,1:napxo)

  ! Module data
  call meta_crstart
  call time_crstart
  call rnd_crstart
  if(dynk_enabled) then
    call dynk_crstart
  end if
  if(scatter_active) then
    call scatter_crstart
  end if
  if(nelens > 0) then
    call elens_crstart
  end if
  if(do_coll) then
    call coll_crstart
  end if

  ! Done
  write(crlog,"(3(a,i0))") "CR_START> SixRecords: ",sixrecs,", SixRecords C/R: ",crsixrecs,", BinRecords: ",binrec
  flush(crlog)

  ! Just throw away our fort.92 stuff.
  call f_close(lout)
  call f_open(unit=lout,file=cr_outFile,formatted=.true.,mode="rw",err=fErr,status="replace")

  write(crlog,"(a)") "CR_START> "//repeat("=",80)
  write(crlog,"(a)") "CR_START> SixTrack Restarted"
  write(crlog,"(a)") "CR_START> "//repeat("=",80)
  flush(crlog)

  write(lout, "(a)") "CR_START> "//repeat("=",80)
  write(lout, "(a)") "CR_START> SixTrack Restarted"
  write(lout, "(a)") "CR_START> "//repeat("=",80)
  flush(lout)

end subroutine crstart

! ================================================================================================ !
!  Reposition Track Files
!  Moved from crcheck
!  Last modified: 2019-09-10
! ================================================================================================ !
subroutine cr_positionTrackFiles

  use crcoall
  use mod_units
  use mod_common
  use, intrinsic :: iso_fortran_env, only : int32

  integer j, k, ia, iau
  integer binrecs9x, binrecs94, tUnit

  ! DANGER: If the length of the records in writebin(_header)changes, these arrays must be updated
  integer(kind=int32) hbuff(253),tbuff(35)

  ! We may be re-running with a DIFFERENT number of turns (numl)
  if(numl /= crnuml) then
    if(numl < crnumlcr) then
      write(lerr, "(2(a,i0))") "CR_CHECK> ERROR New numl < crnumlcr : ",numl," < ",crnumlcr
      write(crlog,"(2(a,i0))") "CR_CHECK> ERROR New numl < crnumlcr : ",numl," < ",crnumlcr
      flush(crlog)
      call prror
    end if
    write(crlog,"(2(a,i0))") "CR_CHECK> Resetting numl in binary file headers from ",crnuml," to ",numl
    flush(crlog)

    ! Reposition binary file singletrackfile.dat
    call f_requestUnit("cr_trackfile.tmp",tUnit)
    call f_open(unit=tUnit,file="cr_trackfile.tmp",formatted=.false.,mode="rw")
    ! First, copy crbinrecs(ia)*(crnapx/2) records of data from singletrackfile.dat to temp file
    binrecs9x = 0

    ! Copy headers
    do ia=1,crnapxo/2,1
      read(90,err=105,end=105,iostat=ierro) hbuff
      binrecs9x = binrecs9x + 1
      hbuff(51) = numl ! Reset the number of turns (not very elegant)
      write(tUnit,err=105,iostat=ierro) hbuff
    end do

    ! Copy particle tracking data
    do ia=1,crnapxo/2,1
      do j=2,crbinrecs(ia)
        if(ntwin /= 2) then
          read(90,err=105,end=105,iostat=ierro) (tbuff(k),k=1,17)
          write(tUnit,err=105,iostat=ierro) (tbuff(k),k=1,17)
        else
          read(90,err=105,end=105,iostat=ierro) tbuff
          write(tUnit,err=105,iostat=ierro) tbuff
        end if
        binrecs9x = binrecs9x + 1
      end do
    end do

    ! Second, copy crbinrecs(ia)*(crnapx/2) records of data from temp file to singletrackfile.dat
    rewind(tUnit)
    rewind(90)
    binrecs94 = 0

    ! Copy header
    do ia=1,crnapxo/2,1
      read(tUnit,err=105,end=105,iostat=ierro) hbuff
      binrecs94 = binrecs94 + 1
      write(90,err=105,iostat=ierro) hbuff
    end do

    ! Copy particle tracking data
    do ia=1,crnapxo/2,1
      do j=2,crbinrecs(ia)
        if(ntwin /= 2) then
          read(tUnit,err=105,end=105,iostat=ierro) (tbuff(k),k=1,17)
          write(90,err=105,iostat=ierro) (tbuff(k),k=1,17)
        else
          read(tUnit,err=105,end=105,iostat=ierro) tbuff
          write(90,err=105,iostat=ierro) tbuff
        end if
        binrecs94 = binrecs94 + 1
      end do
    end do

    ! This is not a FLUSH!
    endfile(90,iostat=ierro)
    backspace(90,iostat=ierro)
    call f_freeUnit(tUnit)
  else !ELSE for "if(nnuml.ne.crnuml) then" -> here we treat nnuml.eq.crnuml, i.e. the number of turns have not been changed
    ! Now with the new array crbinrecs we can ignore files which are
    ! basically finished because a particle has been lost.......
    ! Just check crbinrecs against crbinrec
    binrecs9x = 0
    ! Reposition headers
    do ia=1,crnapxo/2,1
      read(90,err=102,end=102,iostat=ierro) hbuff
      binrecs9x=binrecs9x+1
    end do

    ! Reposition track records
    do ia=1,crnapxo/2,1
      do j=2,crbinrecs(ia)
        if(ntwin /= 2) then !ntwin=1
          read(90,err=102,end=102,iostat=ierro) (tbuff(k),k=1,17)
        else                !ntwin=2
          read(90,err=102,end=102,iostat=ierro) tbuff
        end if
        binrecs9x = binrecs9x + 1
      end do
    end do
  end if ! END "if (numl.ne.crnuml) then" and END else
  return

102 continue
  write(lerr,"(2(a,i0))") "CR_CHECK> ERROR Re-reading singletrackfile.dat for ia=",ia," IOSTAT=",ierro
  write(lerr,"(2(a,i0))") "CR_CHECK>       binrecs9x=",binrecs9x," Expected crbinrecs=",crbinrecs(ia)
  call prror
105 continue
  write(lerr,"(2(a,i0))") "CR_CHECK> ERROR Copying particle pair ",ia," IOSTAT=",ierro," from/to singletrackfile.dat"
  write(lerr,"(3(a,i0))") "CR_CHECK>       binrecs9x=",binrecs9x," Expected crbinrecs=",crbinrecs(ia)," binrecs94=",binrecs94
  call prror

end subroutine cr_positionTrackFiles

! ================================================================================================ !
!  Copy lout to output_unit
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2019-04-29
! ================================================================================================ !
subroutine cr_copyOut

  use, intrinsic :: iso_fortran_env, only : output_unit

  use crcoall
  use mod_units
  use mod_common, only : fort6

  character(len=1024) inLine
  integer ioStat, lnSize, nLines

  if(lout == output_unit) return

  flush(lout)
  rewind(lout)
  nLines = 0
10 continue
  read(lout,"(a1024)",end=20,err=20,iostat=ioStat,size=lnSize,advance="no") inLine
  if(ioStat > 0) goto 20 ! Do not use /= 0

  write(output_unit,"(a)",err=30,iostat=ioStat) inLine(1:lnSize)
  if(ioStat /= 0) goto 30

  nLines = nLines + 1
  goto 10

20 continue ! Done copying lout > output_unit
  flush(output_unit)

  call f_close(lout)
  call f_open(unit=lout,file=cr_outFile,formatted=.true.,mode="rw-",status="replace")

  write(crlog,"(2(a,i0))") "COPY_OUT> Copied ",nLines," lines from "//cr_outFile//" to "//trim(fort6)
  flush(crlog)

  sixrecs = sixrecs + nLines

  return

30 continue ! Write error on fort.6
  write(crlog,"(2(a,i0))") "COPY_OUT> Failed to copy "//cr_outFile//" to "//trim(fort6)
  flush(crlog)

end subroutine cr_copyOut

end module checkpoint_restart
