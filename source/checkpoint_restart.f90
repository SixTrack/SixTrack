! ================================================================================================ !
!  SixTrack Checkpoint/Restart Module
! ====================================
!  E. McIntosh
!  K.N. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2019-04-29
! ================================================================================================ !
module checkpoint_restart

  use floatPrecision

  implicit none

  ! Checkpoint Files
  character(len=15), public,  save :: cr_pntFile(2)  = ["crpoint_pri.bin","crpoint_sec.bin"]
  integer,           public,  save :: cr_pntUnit(2)  = -1
  logical,           public,  save :: cr_pntExist(2) = .false.
  logical,           public,  save :: cr_pntRead(2)  = .false.

  ! Logging Files
  character(len=13), public,  save :: cr_errFile  = "cr_stderr.tmp"
  character(len=13), public,  save :: cr_outFile  = "cr_stdout.tmp"
  character(len=13), public,  save :: cr_logFile  = "cr_status.log"
  integer,           parameter     :: cr_errUnit  = 91
  integer,           parameter     :: cr_outUnit  = 92
  integer,           parameter     :: cr_logUnit  = 93

  ! Checkpoint/Restart Flags and Variables
  logical,           public,  save :: cr_rerun    = .false.
  logical,           public,  save :: cr_start    = .true.
  logical,           public,  save :: cr_restart  = .false.
  logical,           public,  save :: cr_checkp   = .false.
  logical,           private, save :: cr_sythck   = .false.

  character(len=21), public,  save :: cr_startMsg = " "
  real,              public,  save :: cr_time     = 0.0
  integer,           public,  save :: cr_numl     = 1
  integer,           public,  save :: binrec      = 0       ! The maximum number of records writen for all tracking data files
  integer,           public,  save :: sixrecs     = 0

  ! Keep length in sync with version.f90
  character(len=8),  private, save :: cr_version  = " "
  character(len=10), private, save :: cr_moddate  = " "

  ! C/R Temp Variables and Arrays
  real(kind=fPrec),  private, save :: cre0
  real(kind=fPrec),  private, save :: crbetrel
  real(kind=fPrec),  private, save :: crbrho

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

  integer,          allocatable, public,  save :: binrecs(:)    ! ((npart+1)/2)
  integer,          allocatable, public,  save :: crbinrecs(:)  ! (npart+1)/2)
  integer,          allocatable, private, save :: crnumxv(:)    ! (npart)
  integer,          allocatable, private, save :: crnnumxv(:)   ! (npart)
  integer,          allocatable, private, save :: crpartID(:)   ! (npart)
  integer,          allocatable, private, save :: crparentID(:) ! (npart)

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

  use mod_units

  integer i
  logical fErr

  do i=1,2
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
  use numerical_constants, only : zero

  integer, intent(in) :: npart_new

  integer :: npair_new
  npair_new = npart_new/2 + 1

  call alloc(crxv1,      npart_new,    zero,    "crxv1")
  call alloc(crxv2,      npart_new,    zero,    "crxv2")
  call alloc(cryv1,      npart_new,    zero,    "cryv1")
  call alloc(cryv2,      npart_new,    zero,    "cryv2")
  call alloc(crsigmv,    npart_new,    zero,    "crsigmv")
  call alloc(crdpsv,     npart_new,    zero,    "crdpsv")
  call alloc(crdpsv1,    npart_new,    zero,    "crdpsv1")
  call alloc(crejv,      npart_new,    zero,    "crejv")
  call alloc(crejfv,     npart_new,    zero,    "crejfv")
  call alloc(craperv,    npart_new, 2, zero,    "craperv")
  call alloc(binrecs,    npair_new,    0,       "binrecs")
  call alloc(crbinrecs,  npair_new,    0,       "crbinrecs")
  call alloc(crnumxv,    npart_new,    0,       "crnumxv")
  call alloc(crnnumxv,   npart_new,    0,       "crnnumxv")
  call alloc(crpartID,   npart_new,    0,       "crpartID")
  call alloc(crparentID, npart_new,    0,       "crparentID")
  call alloc(crpstop,    npart_new,    .false., "crpstop")
  call alloc(crllostp,   npart_new,    .false., "crllostp")

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
  integer pTurn, nKills, i, iUnit

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

  call f_requestUnit("crkillswitch.tmp",iUnit)

  inquire(file="crkillswitch.tmp",exist=fExist)
  if(fExist .eqv. .false.) then
    open(iUnit,file="crkillswitch.tmp",form="unformatted",access="stream",status="replace",action="write")
    write(iUnit) 0,0
    flush(iUnit)
    close(iUnit)
  end if

  open(iUnit,file="crkillswitch.tmp",form="unformatted",access="stream",status="old",action="read")
  read(iUnit) pTurn,nKills
  flush(iUnit)
  close(iUnit)
  if(st_debug .and. pTurn > 0) then
    write(lout,"(a,i0)") "CRKILLSW> Kill switch previously triggered on turn ",pTurn
    write(93,  "(a,i0)") "CRKILLSW> Kill switch previously triggered on turn ",pTurn
    flush(lout)
    flush(93)
  end if

  do i=1,size(st_killturns,1)
    if(iTurn == st_killturns(i) .and. iTurn > pTurn) then
      killIt = .true.
      exit
    end if
  end do

  if(killIt) then
    nKills = nKills + 1

    write(lout,"(a,i0)") "CRKILLSW> Triggering kill switch on turn ",iTurn
    write(93,  "(a,i0)") "CRKILLSW> Triggering kill switch on turn ",iTurn
    flush(lout)
    flush(93)

    open(iUnit,file="crrestartme.tmp",form="unformatted",access="stream",status="replace",action="write")
    write(iUnit) 1
    flush(iUnit)
    close(iUnit)

    open(iUnit,file="crkillswitch.tmp",form="unformatted",access="stream",status="replace",action="write")
    write(iUnit) iTurn,nKills
    flush(iUnit)
    close(iUnit)
    stop
  end if

end subroutine cr_killSwitch

! ================================================================================================ !
!  Check CR Files
! ================
!  Last modified: 2018-12-05
!
!  This subroutine checks if the C/R files exist, and if so tries to load them into the cr* variables.
!  This routine also repositions the output files for fort.90..91-napx/2 and various other modules
! ================================================================================================ !
subroutine crcheck

  use, intrinsic :: iso_fortran_env, only : output_unit

  use parpro
  use crcoall
  use mod_common
  use mod_commons
  use mod_common_main
  use mod_version

  use dynk,      only : dynk_enabled, dynk_noDynkSets,dynk_crcheck_readdata,dynk_crcheck_positionFiles
  use dump,      only : dump_crcheck_readdata,dump_crcheck_positionFiles
  use aperture,  only : limifound, aper_crcheck_readdata, aper_crcheck_positionFiles
  use scatter,   only : scatter_active, scatter_crcheck_readdata, scatter_crcheck_positionFiles
  use mod_hions, only : hions_crcheck_readdata
  use elens,     only : melens, elens_crcheck
  use mod_meta,  only : meta_crcheck

  integer j,k,l,m
  integer nPoint, ioStat

  logical noRestart, rErr
  character(len=1024) inLine

  ! Check the size of CR arrays
  if(npart /= crnpart_old) call cr_expand_arrays(npart)

  cr_restart    = .false.
  cr_pntRead(:) = .false.

  ! Some log entries to fort.93
  write(93,"(3(a,l1))") "CR_CHECK> Called with restart = ",cr_restart,", rerun = ",cr_rerun,", and checkpoint = ",cr_checkp
  flush(93)

  ! We are not checkpoint/restart or we have no restart files
  if(cr_checkp .eqv. .false.) goto 200
  noRestart = .false.
  do nPoint=1,2
    noRestart = noRestart .and. .not. cr_pntExist(nPoint)
  end do
  if(noRestart) goto 200

  ! If we do we must have a fort.6 as it was created by CRPOINT
  ! NOT TRUE anymore??? We might be NOT rerun but using a Sixin.zip
#ifndef BOINC
  if(cr_rerun .eqv. .false.) then
    write(lerr,"(a)") "CR_CHECK> ERROR Found "//cr_pntFile(1)//"/"//cr_pntFile(2)//" but no fort.6"
    call prror
  end if
#endif

  ! Check at least one restart file is readable
  do nPoint=1,2
    write(93,"(a)") "CR_CHECK> Checking file "//cr_pntFile(nPoint)
    flush(93)

    if(cr_pntExist(nPoint) .eqv. .false.) then
      write(93,"(a)") "CR_CHECK> File "//cr_pntFile(nPoint)//" does not exist"
      cycle
    end if

    rewind(cr_pntUnit(nPoint))

    write(93,"(a)") "CR_CHECK>  * SixTrack version"
    flush(93)

    read(cr_pntUnit(nPoint),iostat=ioStat) cr_version,cr_moddate
    if(cr_version == " " .or. cr_moddate == " ") then
      write(93,"(a)") "CR_CHECK> Unknown SixTrack version. Skipping this file."
      flush(93)
      cycle
    end if
    if(cr_version /= version .or. cr_moddate /= moddate) then
      write(93,"(a)") "CR_CHECK> Checkpoint files "//cr_pntFile(nPoint)//" was written by SixTrack version "//&
        trim(cr_version)//" with release date "//trim(cr_moddate)
      write(93,"(a)") "CR_CHECK> This is SixTrack version "//trim(version)//" with release date "//trim(moddate)
      write(93,"(a)") "CR_CHECK> Version mismatch. Skipping this file."
      flush(93)
      cycle
    end if
    if(ioStat /= 0) cycle

    write(93,"(a)") "CR_CHECK>  * Tracking variables"
    flush(93)
    read(cr_pntUnit(nPoint),iostat=ioStat) crnumlcr,crnuml,crsixrecs,crbinrec,cr_sythck,cril, &
      cr_time,crnapxo,crnapx,cre0,crbetrel,crbrho
    if(ioStat /= 0) cycle

    write(93,"(a)") "CR_CHECK>  * Particle arrays"
    flush(93)
    read(cr_pntUnit(nPoint),iostat=ioStat) &
      (crbinrecs(j), j=1,(crnapxo+1)/2), &
      (crnumxv(j),   j=1,crnapxo),       &
      (crnnumxv(j),  j=1,crnapxo),       &
      (crpartID(j),  j=1,crnapxo),       &
      (crparentID(j),j=1,crnapxo),       &
      (crpstop(j),   j=1,crnapxo),       &
      (crxv1(j),     j=1,crnapxo),       &
      (cryv1(j),     j=1,crnapxo),       &
      (crxv2(j),     j=1,crnapxo),       &
      (cryv2(j),     j=1,crnapxo),       &
      (crsigmv(j),   j=1,crnapxo),       &
      (crdpsv(j),    j=1,crnapxo),       &
      (crdpsv1(j),   j=1,crnapxo),       &
      (crejv(j),     j=1,crnapxo),       &
      (crejfv(j),    j=1,crnapxo),       &
      (craperv(j,1), j=1,crnapxo),       &
      (craperv(j,2), j=1,crnapxo),       &
      (crllostp(j),  j=1,crnapxo)
    if(ioStat /= 0) cycle

    write(93,"(a)") "CR_CHECK>  * META variables"
    flush(93)
    call meta_crcheck(cr_pntUnit(nPoint),rErr)
    if(rErr) cycle

    write(93,"(a)") "CR_CHECK>  * DUMP variables"
    flush(93)
    call dump_crcheck_readdata(cr_pntUnit(nPoint),rErr)
    if(rErr) cycle

    write(93,"(a)") "CR_CHECK>  * HION variables"
    flush(93)
    call hions_crcheck_readdata(cr_pntUnit(nPoint),rErr)
    if(rErr) cycle

    if(dynk_enabled) then
      write(93,"(a)") "CR_CHECK>  * DYNK variables"
      flush(93)
      call dynk_crcheck_readdata(cr_pntUnit(nPoint),rErr)
      if(rErr) cycle
    end if

    if(scatter_active) then
      write(93,"(a)") "CR_CHECK>  * SCATTER variables"
      flush(93)
      call scatter_crcheck_readdata(cr_pntUnit(nPoint),rErr)
      if(rErr) cycle
    end if

    if(limifound) then
      write(93,"(a)") "CR_CHECK>  * APERTURE variables"
      flush(93)
      call aper_crcheck_readdata(cr_pntUnit(nPoint),rErr)
      if(rErr) cycle
    end if

    if(melens > 0) then
      write(93,"(a)") "CR_CHECK>  * ELENS variables"
      flush(93)
      call elens_crcheck(cr_pntUnit(nPoint),rErr)
      if(rErr) cycle
    end if

    ! New extended checkpoint for synuthck (ERIC)
    if(cr_sythck) then
      ! and make sure we can read the extended vars before leaving fort.95
      ! We will re-read them in crstart to be sure they are restored correctly
      write(93,"(a)") "CR_CHECK>  * THICK EXTENDED arrays"
      flush(93)
      read(cr_pntUnit(nPoint),iostat=ioStat) &
        ((((al(k,m,j,l),l=1,il),j=1,crnapxo),m=1,2),k=1,6), &
        ((((as(k,m,j,l),l=1,il),j=1,crnapxo),m=1,2),k=1,6), &
        (dpd(j),j=1,crnapxo),(dpsq(j),j=1,crnapxo),(fokqv(j),j=1,crnapxo)
      backspace(cr_pntUnit(nPoint),iostat=ioStat)
      if(ioStat /= 0) cycle
      write(93,"(a)") "CR_CHECK> Read "//cr_pntFile(nPoint)//" EXTENDED OK"
      write(93,"(a)") "CR_CHECK> Leaving "//cr_pntFile(1)//" for CRSTART EXTENDED"
      flush(93)
    end if

    write(93,"(a)") "CR_CHECK> File "//cr_pntFile(nPoint)//" successfully read"
    flush(93)
    cr_pntRead(nPoint) = .true.
    goto 100

  end do

  ! If we're here, we failed to read a checkpoint file
  write(93,"(a)") "CR_CHECK> ERROR Could not read checkpoint files"
  flush(93)
  goto 200

100 continue

  ! If we have successfully read one of the checkpoint files, we need to handle lost particles and ntwin /= 2
  ! Otherwise we just continue with checkpointing as requested
  write(93,"(2(a,i8))") "CR_CHECK> Particles  C/R: ",crnapxo,  ", Input:  ",napx*2
  write(93,"(2(a,i8))") "CR_CHECK> SixRecords C/R: ",crsixrecs,", Buffer: ",sixrecs
  write(93,"(1(a,i8))") "CR_CHECK> BinRecords C/R: ",crbinrec
#ifndef STF
  do j=1,(crnapxo+1)/2
    write(93,"(2(a,i0))") "CR_CHECK>  * Record ",j,": ",crbinrecs(j)
  end do
#endif
  flush(93)

  ! First we position fort.6 to last checkpoint
  do j=1,crsixrecs
    read(output_unit,"(a1024)",end=110,err=120,iostat=ioStat) inLine
    sixrecs = sixrecs + 1
  end do
  ! This is not a FLUSH!
  endfile(output_unit,iostat=ierro)
110 continue
  backspace(output_unit,iostat=ierro)
  write(93,"(a,i0,a)") "CR_CHECK> Found ",sixrecs," records in fort.6"
  flush(93)

  !  Position Files
  ! ================

  call cr_positionTrackFiles

  if(dynk_enabled .and. .not.dynk_noDynkSets) then
    write(93,"(a)") "CR_CHECK> Repositioning dynksets.dat"
    flush(93)
    call dynk_crcheck_positionFiles
  end if

  write(93,"(a)") "CR_CHECK> Repositioning DUMP files"
  flush(93)
  call dump_crcheck_positionFiles

  if(scatter_active) then
    write(93,"(a)") "CR_CHECK> Repositioning SCATTER files"
    flush(93)
    call scatter_crcheck_positionFiles
  end if

  if(limifound) then
    write(93,"(a)") "CR_CHECK> Repositioning APERTURE files"
    flush(93)
    call aper_crcheck_positionFiles
  end if

  ! Set up flag for tracking routines to call CRSTART
  cr_restart = .true.
  write(lout,"(a)") "SIXTRACK> "//repeat("=",80)
  write(lout,"(a)") "SIXTRACK>  Restarted"
  write(lout,"(a)") "SIXTRACK> "//repeat("=",80)
  flush(lout)

  return

  ! We are not checkpointing or we have no checkpoints, or we have no readable checkpoint
  ! If not checkpointing we can just give up on lout and use fort.6. We don't need to count records at all
200 continue
  write(93,"(a,l1)") "SIXTRACK> No restart possible. Checkpoint = ",cr_checkp
  flush(93)
  if(.not.cr_checkp) then
    call cr_copyOut
  end if

  return

120 continue
  write(lerr,"(a)") "CR_CHECK> ERROR Cannot read fort.6"
  write(93,"(a)")   "CR_CHECK> ERROR Cannot read fort.6"
  flush(93)
  call prror

end subroutine crcheck

! ================================================================================================ !
!  Write Checkpoint
! ==================
!  Last modified: 2019-04-29
!  This subroutine writes the checkpoint data to the binary checpoint files
! ================================================================================================ !
subroutine crpoint

  use mod_time
  use mod_common
  use mod_commons
  use mod_common_main
  use mod_common_track
  use mod_version
  use mod_settings
  use numerical_constants

  use dynk,      only : dynk_enabled,dynk_getvalue,dynk_fSets_cr,dynk_cSets_unique,dynk_nSets_unique,dynk_filePos,dynk_crpoint
  use dump,      only : dump_crpoint
  use aperture,  only : aper_crpoint,limifound
  use scatter,   only : scatter_active, scatter_crpoint
  use elens,     only : melens, elens_crpoint
  use mod_meta,  only : meta_crpoint
  use mod_hions, only : hions_crpoint

  integer j, k, l, m, nPoint
  logical wErr, fErr

  write(93,"(3(a,i0))") "CR_POINT> Called on turn ",numx," / ",numl," : frequency is ",numlcp
  flush(93)

  if(cr_restart) then
    cr_restart = .false.
    return
  end if

  ! We need to copy fort.92 (lout) to fort.6 (sixrecs) (if it exists and we are not already using fort.6)
  call cr_copyOut

  ! Hope this is correct
  ! Maybe not!!!! this should be accumulative over multiple C/Rs
  call time_timerCheck(time3)
  time3 = (time3-time1)+cr_time

  crnumlcr = numx+1

  if(dynk_enabled) then ! Store current settings of elements affected by DYNK
    if(st_debug) then
      write(93,"(a)") "CR_POINT> Filling DYNK sets"
      flush(93)
    end if
    do j=1,dynk_nSets_unique
      dynk_fSets_cr(j) = dynk_getvalue(dynk_cSets_unique(j,1),dynk_cSets_unique(j,2))
    end do
  end if

  ! ********************
  !  Write the CR files
  ! ********************

  do nPoint=1,2

    wErr = .false.
    rewind(cr_pntUnit(nPoint))

    if(st_debug) then
      write(93,"(a)") "CR_POINT> Writing to checkpoint file "//cr_pntFile(nPoint)
      write(93,"(a)") "CR_POINT>  * SixTrack version"
      flush(93)
    end if
    write(cr_pntUnit(nPoint),err=100,iostat=ierro) version, moddate

    if(st_debug) then
      write(93,"(a)") "CR_POINT>  * Tracking variables"
      flush(93)
    end if
    write(cr_pntUnit(nPoint),err=100,iostat=ierro) crnumlcr, numl, sixrecs, binrec, sythckcr, il,   &
      time3, napxo, napx, e0, betrel, brho

    if(st_debug) then
      write(93,"(a)") "CR_POINT>  * Particle arrays"
      flush(93)
    end if
    write(cr_pntUnit(nPoint),err=100,iostat=ierro) &
      (binrecs(j), j=1,(napxo+1)/2), &
      (numxv(j),   j=1,napxo),       &
      (nnumxv(j),  j=1,napxo),       &
      (partID(j),  j=1,napxo),       &
      (parentID(j),j=1,napxo),       &
      (pstop(j),   j=1,napxo),       &
      (xv1(j),     j=1,napxo),       &
      (yv1(j),     j=1,napxo),       &
      (xv2(j),     j=1,napxo),       &
      (yv2(j),     j=1,napxo),       &
      (sigmv(j),   j=1,napxo),       &
      (dpsv(j),    j=1,napxo),       &
      (dpsv1(j),   j=1,napxo),       &
      (ejv(j),     j=1,napxo),       &
      (ejfv(j),    j=1,napxo),       &
      (aperv(j,1), j=1,napxo),       &
      (aperv(j,2), j=1,napxo),       &
      (llostp(j),  j=1,napxo)
    flush(cr_pntUnit(nPoint))

    if(st_debug) then
      write(93,"(a)") "CR_POINT>  * META variables"
      flush(93)
    end if
    call meta_crpoint(cr_pntUnit(nPoint),wErr,ierro)
    if(wErr) goto 100

    if(st_debug) then
      write(93,"(a)") "CR_POINT>  * DUMP variables"
      flush(93)
    end if
    call dump_crpoint(cr_pntUnit(nPoint), wErr,ierro)
    if(wErr) goto 100

    if(st_debug) then
      write(93,"(a)") "CR_POINT>  * HION variables"
      flush(93)
    end if
    call hions_crpoint(cr_pntUnit(nPoint),wErr,ierro)
    if(wErr) goto 100

    if(dynk_enabled) then
      if(st_debug) then
        write(93,"(a)") "CR_POINT>  * DYNK variables"
        flush(93)
      end if
      call dynk_crpoint(cr_pntUnit(nPoint),wErr,ierro)
      if(wErr) goto 100
    end if

    if(scatter_active) then
      if(st_debug) then
        write(93,"(a)") "CR_POINT>  * SCATTER variables"
        flush(93)
      end if
      call scatter_crpoint(cr_pntUnit(nPoint),wErr,ierro)
      if(wErr) goto 100
    end if

    if(limifound) then
      if(st_debug) then
        write(93,"(a)") "CR_POINT>  * APERTURE variables"
        flush(93)
      end if
      call aper_crpoint(cr_pntUnit(nPoint),wErr,ierro)
      if(wErr) goto 100
    end if

    if(melens > 0) then
      if(st_debug) then
        write(93,"(a)") "CR_POINT>  * ELENS variables"
        flush(93)
      end if
      call elens_crpoint(cr_pntUnit(nPoint),wErr,ierro)
      if(wErr) goto 100
    end if

    if(sythckcr) then
      if(ithick == 1) then
        if(st_debug) then
          write(93,"(a)") "CR_POINT>  * THICK EXTENDED arrays"
          flush(93)
        end if
        write(cr_pntUnit(nPoint),err=100,iostat=ierro) &
          ((((al(k,m,j,l),l=1,il),j=1,napxo),m=1,2),k=1,6), &
          ((((as(k,m,j,l),l=1,il),j=1,napxo),m=1,2),k=1,6)
        flush(cr_pntUnit(nPoint))
      end if

      if(st_debug) then
        write(93,"(a)") "CR_POINT>  * THICK arrays"
        flush(93)
      end if
      write(cr_pntUnit(nPoint),err=100,iostat=ierro) &
        (dpd(j),j=1,napxo),(dpsq(j),j=1,napxo),(fokqv(j),j=1,napxo)
    end if

    flush(93)
    flush(cr_pntUnit(nPoint))

  end do ! Loop over nPoint

  return

100 continue
  write(93,"(a,i0)") "CR_POINT> ERROR Writing checkpoint file. iostat = ",ierro
  flush(93)
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

  use mod_hions
  use dynk,     only : dynk_enabled, dynk_crstart
  use scatter,  only : scatter_active, scatter_crstart
  use elens,    only : melens, elens_crstart
  use mod_meta, only : meta_crstart

  logical fErr
  integer j,k,l,m

  write(93,"(a,i0)") "CR_START> Starting from turn ",crnumlcr
  flush(93)
  cr_numl = crnumlcr

  ! We do NOT reset numl so that a run can be extended for
  ! for more turns from the last checkpoint
  binrec   = crbinrec
  sythckcr = cr_sythck

  ! the cr_time is required (crtime0/1 removed)
  napxo  = crnapxo
  napx   = crnapx
  e0     = cre0
  e0f    = sqrt(e0**2-nucm0**2)
  betrel = crbetrel
  brho   = crbrho

  write(93,"(a)") "CR_START> Loading BinRecords"
  do j=1,(napxo+1)/2
    binrecs(j) = crbinrecs(j)
  end do

  write(93,"(a)") "CR_START> Loading tracking data"
  flush(93)

  partID(1:napxo)    = crpartID(1:napxo)
  parentID(1:napxo)  = crparentID(1:napxo)
  pstop(1:napxo)     = crpstop(1:napxo)
  llostp(1:napxo)    = crllostp(1:napxo)

  xv1(1:napxo)       = crxv1(1:napxo)
  yv1(1:napxo)       = cryv1(1:napxo)
  xv2(1:napxo)       = crxv2(1:napxo)
  yv2(1:napxo)       = cryv2(1:napxo)
  sigmv(1:napxo)     = crsigmv(1:napxo)

  dpsv(1:napxo)      = crdpsv(1:napxo)
  dpsv1(1:napxo)     = crdpsv1(1:napxo)
  ejv(1:napxo)       = crejv(1:napxo)
  ejfv(1:napxo)      = crejfv(1:napxo)

  numxv(1:napxo)     = crnumxv(1:napxo)
  nnumxv(1:napxo)    = crnnumxv(1:napxo)
  do j=1,napxo
    if(pstop(j) .eqv. .false.) then
      numxv(j)  = numl
      nnumxv(j) = numl
    end if
  end do

  ! Recompute from loaded arrays
  oidpsv(1:napxo) = one/(one + dpsv(1:napxo))
  rvv(1:napxo)    = (ejv(1:napxo)*e0f)/(e0*ejfv(1:napxo))

  ! Aperture data
  aperv(1:napxo,1:2) = craperv(1:napxo,1:2)

  ! Module data
  call hions_crstart
  call meta_crstart
  if(dynk_enabled) then
    call dynk_crstart
  end if
  if(scatter_active) then
    call scatter_crstart
  end if
  if(melens > 0) then
    call elens_crstart
  end if

  ! New extended checkpoint for synuthck (ERIC)
  if(cr_sythck) then
    ! Now read the extended vars from fort.95/96.
    if(cril /= il) then
      write(lout,"(2(a,i0))") "CR_START> ERROR Problem as cril/il are different cril = ",cril,", il = ",il
      write(93,  "(2(a,i0))") "CR_START> ERROR Problem as cril/il are different cril = ",cril,", il = ",il
      flush(93)
      call prror
    end if
    if(cr_pntRead(1)) then
      if(ithick == 1) then
        read(cr_pntUnit(1),end=100,err=100,iostat=ierro) &
          ((((al(k,m,j,l),l=1,il),j=1,napxo),m=1,2),k=1,6),&
          ((((as(k,m,j,l),l=1,il),j=1,napxo),m=1,2),k=1,6)
      end if

      read(cr_pntUnit(1),end=100,err=100,iostat=ierro) &
        (dpd(j),j=1,napxo),(dpsq(j),j=1,napxo),(fokqv(j),j=1,napxo)
      write(93,"(a)") "CR_START> Read "//cr_pntFile(1)//" EXTENDED OK"
      flush(93)
      goto 102
    end if
    if(cr_pntRead(2)) then
      if(ithick == 1) then
        read(cr_pntUnit(2),end=101,err=101,iostat=ierro) &
          ((((al(k,m,j,l),l=1,il),j=1,napxo),m=1,2),k=1,6), &
          ((((as(k,m,j,l),l=1,il),j=1,napxo),m=1,2),k=1,6)
      end if
      read(cr_pntUnit(2),end=101,err=101,iostat=ierro) &
        (dpd(j),j=1,napxo),(dpsq(j),j=1,napxo),(fokqv(j),j=1,napxo)

      write(93,"(a)") "CR_START> Read "//cr_pntFile(2)//" EXTENDED OK"
      flush(93)
      goto 102
    end if
100 continue
    write(93,"(a,i0)") "CR_START> ERROR Could not read checkpoint file "//cr_pntFile(1)//" (extended), iostat = ",ierro
    call prror
101 continue
    write(93,"(a,i0)") "CR_START> ERROR Could not read checkpoint file "//cr_pntFile(2)//" (extended), iostat = ",ierro
    call prror
  end if

102 continue
  write(93,"(3(a,i0))") "CR_START> Sixrecs = ",sixrecs,", crsixrecs = ",crsixrecs,", binrec = ",binrec
  flush(93)

  ! Just throw away our fort.92 stuff.
  call f_close(lout)
  call f_open(unit=lout,file=cr_outFile,formatted=.true.,mode="rw",err=fErr,status="replace")
  write(lout,"(a)") "SIXTRACR> "//repeat("=",80)
  write(lout,"(a)") "SIXTRACR>  Restarted"
  write(lout,"(a)") "SIXTRACR> "//repeat("=",80)
  flush(lout)

end subroutine crstart

! ================================================================================================ !
!  Reposition Track Files
!  Moved from crcheck
!  Last modified: 2019-04-30
! ================================================================================================ !
subroutine cr_positionTrackFiles

  use crcoall
  use mod_units
  use mod_common
  use, intrinsic :: iso_fortran_env, only : int32

  integer j, k, ia, iau
  integer binrecs9x, binrecs94

  ! DANGER: If the length of the records in writebin(_header)changes, these arrays must be updated
  integer(kind=int32) hbuff(253),tbuff(35)

  ! We may be re-running with a DIFFERENT number of turns (numl)
  ! Eric fix this later by reading numl for fort.90
  if(numl /= crnuml) then
    if(numl < crnumlcr) then
      write(lerr,"(2(a,i0))") "CR_CHECK> ERROR New numl < crnumlcr : ",numl," < ",crnumlcr
      write(93,"(2(a,i0))")   "CR_CHECK> ERROR New numl < crnumlcr : ",numl," < ",crnumlcr
      flush(93)
      call prror
    end if
    write(93,"(2(a,i0))") "CR_CHECK> Resetting numl in binary file headers from ",crnuml," to ",numl
    flush(93)

    ! Reposition binary files fort.90 etc. / singletrackfile.dat
    ! fort.94 = temp file where the data from fort.90 etc. is copied to and then back
    call f_open(unit=94,file="fort.94",formatted=.false.,mode="rw")
#ifndef STF
    do ia=1,crnapxo/2,1
      ! First, copy crbinrecs(ia) records of data from fort.91-ia to fort.94
      binrecs9x = 0
      binrecs94 = 0
      iau       = 91-ia

      ! Copy header into integer array hbuff
      read(91-ia,err=105,end=105,iostat=ierro) hbuff
      binrecs9x = binrecs9x + 1
      hbuff(51) = numl ! Reset the number of turns (not very elegant)
      write(94,err=105,iostat=ierro) hbuff

      ! Copy particle tracking data
      do j=2,crbinrecs(ia)
        if(ntwin /= 2) then
          read(91-ia,err=105,end=105,iostat=ierro) (tbuff(k),k=1,17)
          write(94,err=105,iostat=ierro) (tbuff(k),k=1,17)
        else
          read(91-ia,err=105,end=105,iostat=ierro) tbuff
          write(94,err=105,iostat=ierro) tbuff
        end if
        binrecs9x = binrecs9x + 1
      end do

      ! Second, copy crbinrecs(ia) records of data from fort.94 to fort.91-ia
      rewind(94)
      rewind(91-ia)

      ! Copy header
      read(94,err=105,end=105,iostat=ierro) hbuff
      binrecs94 = binrecs94 + 1
      write(91-ia,err=105,iostat=ierro) hbuff

      ! Copy particle tracking data into integer array tbuff
      do j=2,crbinrecs(ia)
        if(ntwin /= 2) then
          read(94,err=105,end=105,iostat=ierro) (tbuff(k),k=1,17)
          write(91-ia,err=105,iostat=ierro) (tbuff(k),k=1,17)
        else
          read(94,err=105,end=105,iostat=ierro) tbuff
          write(91-ia,err=105,iostat=ierro) tbuff
        end if
        binrecs94 = binrecs94 + 1
      end do

      ! This is not a FLUSH!
      endfile(91-ia,iostat=ierro)
      backspace(91-ia,iostat=ierro)
      rewind(94)
    end do
#else
    ! First, copy crbinrecs(ia)*(crnapx/2) records of data from singletrackfile.dat to fort.94
    binrecs9x = 0

    ! Copy headers
    do ia=1,crnapxo/2,1
      read(90,err=105,end=105,iostat=ierro) hbuff
      binrecs9x = binrecs9x + 1
      hbuff(51) = numl ! Reset the number of turns (not very elegant)
      write(94,err=105,iostat=ierro) hbuff
    end do

    ! Copy particle tracking data
    do ia=1,crnapxo/2,1
      do j=2,crbinrecs(ia)
        if(ntwin /= 2) then
          read(90,err=105,end=105,iostat=ierro) (tbuff(k),k=1,17)
          write(94,err=105,iostat=ierro) (tbuff(k),k=1,17)
        else
          read(90,err=105,end=105,iostat=ierro) tbuff
          write(94,err=105,iostat=ierro) tbuff
        end if
        binrecs9x = binrecs9x + 1
      end do
    end do

    ! Second, copy crbinrecs(ia)*(crnapx/2) records of data from fort.94 to singletrackfile.dat
    rewind(94)
    rewind(90)
    binrecs94=0

    ! Copy header
    do ia=1,crnapxo/2,1
      read(94,err=105,end=105,iostat=ierro) hbuff
      binrecs94 = binrecs94 + 1
      write(90,err=105,iostat=ierro) hbuff
    end do

    ! Copy particle tracking data
    do ia=1,crnapxo/2,1
      do j=2,crbinrecs(ia)
        if(ntwin /= 2) then
          read(94,err=105,end=105,iostat=ierro) (tbuff(k),k=1,17)
          write(90,err=105,iostat=ierro) (tbuff(k),k=1,17)
        else
          read(94,err=105,end=105,iostat=ierro) tbuff
          write(90,err=105,iostat=ierro) tbuff
        end if
        binrecs94 = binrecs94 + 1
      end do
    end do

    ! This is not a FLUSH!
    endfile(90,iostat=ierro)
    backspace (90,iostat=ierro)
#endif
    call f_close(94)
  else !ELSE for "if(nnuml.ne.crnuml) then" -> here we treat nnuml.eq.crnuml, i.e. the number of turns have not been changed
    ! Now with the new array crbinrecs we can ignore files which are
    ! basically finished because a particle has been lost.......
    ! Just check crbinrecs against crbinrec
#ifndef STF
    ! Binary files have been rewritten; now re-position
    write(93,"(a)") "CR_CHECK>  * Repositioning binary files"
    do ia=1,crnapxo/2,1
      iau = 91-ia
      if(crbinrecs(ia) >= crbinrec) then
        binrecs9x = 0
        read(91-ia,err=102,end=102,iostat=ierro) hbuff
        do j=2,crbinrecs(ia)
          if(ntwin /= 2) then
            read(91-ia,err=102,end=102,iostat=ierro) (tbuff(k),k=1,17)
          else
            read(91-ia,err=102,end=102,iostat=ierro) tbuff
          end if
          binrecs9x = binrecs9x + 1
        end do

        ! This is not a FLUSH!
        endfile (91-ia,iostat=ierro)
        backspace (91-ia,iostat=ierro)
      else ! Number of ecords written to this file < general number of records written
          ! => Particle has been lost before last checkpoint, no need to reposition.
        write(93,"(2(a,i0))") "CR_CHECK> Ignoring IA ",ia," on unit ",iau
      end if
    end do ! END "do ia=1,crnapxo/2,1"
#else
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
#endif
  end if ! END "if (numl.ne.crnuml) then" and END else
  return

#ifndef STF
102 continue
  write(lerr,"(2(a,i0))") "CR_CHECK> ERROR Re-reading fort.",iau," IOSTAT = ",ierro
  write(lerr,"(3(a,i0))") "CR_CHECK>       Unit ",iau," binrecs9x=",binrecs9x," Expected crbinrecs=",crbinrecs(ia)
  call prror
105 continue
  write(lerr,"(2(a,i0))") "CR_CHECK> ERROR Copying fort.",iau," IOSTAT = ",ierro
  write(lerr,"(4(a,i0))") "CR_CHECK>       Unit ",iau," binrecs9x=",binrecs9x,&
    " Expected crbinrecs=",crbinrecs(ia)," binrecs94=",binrecs94
  call prror
#else
102 continue
  write(lerr,"(2(a,i0))") "CR_CHECK> ERROR Re-reading singletrackfile.dat for ia=",ia," IOSTAT=",ierro
  write(lerr,"(2(a,i0))") "CR_CHECK>       binrecs9x=",binrecs9x," Expected crbinrecs=",crbinrecs(ia)
  call prror
105 continue
  write(lerr,"(2(a,i0))") "CR_CHECK> ERROR Copying particle pair ",ia," IOSTAT=",ierro," from/to singletrackfile.dat"
  write(lerr,"(3(a,i0))") "CR_CHECK>       binrecs9x=",binrecs9x," Expected crbinrecs=",crbinrecs(ia)," binrecs94=",binrecs94
  call prror
#endif
end subroutine cr_positionTrackFiles

! ================================================================================================ !
!  Copy lout to output_unit
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2019-04-29
! ================================================================================================ !
subroutine cr_copyOut

  use crcoall
  use mod_units
  use, intrinsic :: iso_fortran_env, only : output_unit

  character(len=1024) inLine
  integer ioStat, lnSize, nLines

  if(lout == output_unit) return

  flush(lout)
  rewind(lout)
  nLines = 0
10 continue
  read(lout,"(a1024)",end=20,err=20,iostat=ioStat,size=lnSize,advance="no") inLine
  if(ioStat > 0) goto 20 ! End of file (do not use /= 0)

  write(output_unit,"(a)",err=30,iostat=ioStat) inLine(1:lnSize)
  if(ioStat /= 0) goto 30

  nLines = nLines + 1
  goto 10

20 continue ! Done copying lout > output_unit
  flush(output_unit)

  call f_close(lout)
  call f_open(unit=lout,file=cr_outFile,formatted=.true.,mode="rw",status="replace")

  write(93,"(2(a,i0))") "COPY_OUT> Copied ",nLines," from "//cr_outFile//" to fort.6"
  flush(93)

  sixrecs = sixrecs + nLines

  return

30 continue ! Write error on fort.6
  write(93,"(2(a,i0))") "COPY_OUT> Failed to copy "//cr_outFile//" to fort.6"
  flush(93)

end subroutine cr_copyOut

end module checkpoint_restart
