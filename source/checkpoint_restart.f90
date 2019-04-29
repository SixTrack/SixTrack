! ================================================================================================ !
!  SixTrack Checkpoint/Restart Module
! ====================================
!  E. McIntosh
!  K.N. Sjobak, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2019-04-29
! ================================================================================================ !
module checkpoint_restart

  use parpro
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

  real(kind=fPrec), allocatable, private, save :: crxv(:,:)     ! (2,npart)
  real(kind=fPrec), allocatable, private, save :: cryv(:,:)     ! (2,npart)
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

  call alloc(crxv,      2, npart_new,    zero,    "crxv")
  call alloc(cryv,      2, npart_new,    zero,    "cryv")
  call alloc(crsigmv,      npart_new,    zero,    "crsigmv")
  call alloc(crdpsv,       npart_new,    zero,    "crdpsv")
  call alloc(crdpsv1,      npart_new,    zero,    "crdpsv1")
  call alloc(crejv,        npart_new,    zero,    "crejv")
  call alloc(crejfv,       npart_new,    zero,    "crejfv")
  call alloc(craperv,      npart_new, 2, zero,    "craperv")
  call alloc(binrecs,      npair_new,    0,       "binrecs")
  call alloc(crbinrecs,    npair_new,    0,       "crbinrecs")
  call alloc(crnumxv,      npart_new,    0,       "crnumxv")
  call alloc(crnnumxv,     npart_new,    0,       "crnnumxv")
  call alloc(crpartID,     npart_new,    0,       "crpartID")
  call alloc(crparentID,   npart_new,    0,       "crparentID")
  call alloc(crpstop,      npart_new,    .false., "crpstop")
  call alloc(crllostp,     npart_new,    .false., "crllostp")

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
    write(lout,"(a,i0)") "CRKILL> Kill switch previously triggered on turn ",pTurn
    write(93,  "(a,i0)") "SIXTRACR> Kill switch previously triggered on turn ",pTurn
  end if

  do i=1,size(st_killturns,1)
    if(iTurn == st_killturns(i) .and. iTurn > pTurn) then
      killIt = .true.
      exit
    end if
  end do

  if(killIt) then
    nKills = nKills + 1

    write(lout,"(a,i0)") "CRKILL> Triggering kill switch on turn ",iTurn
    write(93,  "(a,i0)") "SIXTRACR> Triggering kill switch on turn ",iTurn

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
!  CRCHECK
!  Last modified: 2018-12-05
!
!  This subroutine checks if the C/R files exist, and if so tries to load them into the cr* variables.
!  This routine also repositions the output files for fort.90..91-napx/2 or STF, DUMP, DYNK and
!     aperture losses
!
!  The file fort.93 is used as a log file for the checkpoint/restarting.
! ================================================================================================ !
subroutine crcheck

  use floatPrecision
  use string_tools
  use numerical_constants
  use dynk,    only : dynk_enabled, dynk_noDynkSets,dynk_crcheck_readdata,dynk_crcheck_positionFiles
  use dump,    only : dump_crcheck_readdata,dump_crcheck_positionFiles
  use aperture,only : aper_crcheck_readdata,aper_crcheck_positionFiles,limifound,losses_filename
  use scatter, only : scatter_active,scatter_crcheck_readdata,scatter_crcheck_positionFiles
  use elens,   only : melens, elens_crcheck
  use, intrinsic :: iso_fortran_env, only : int32
  use crcoall
  use parpro
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da
  use mod_hions
  use mod_version
  use mod_meta

  implicit none

  integer i,j,k,l,m,ia
  integer lstring,myia,mybinrecs,binrecs94

  !DANGER: IF THE LENGTH OF THE RECORDS IN WRITEBIN(_HEADER)CHANGES,
  ! THESE ARRAYS MUST BE UPDATED
  integer(int32) hbuff(253),tbuff(35)

  logical lopen,lerror
  character(len=1024) arecord

#ifdef BOINC
  character(len=256) filename
#endif
  save
  !--------------------------------------------------------------------!
  ! Check the size of CR arrays
  if(npart /= crnpart_old) call cr_expand_arrays(npart)

  cr_restart = .false.
  cr_pntRead(:) = .false.
  ! Some log entries to fort.93
  write(93,"(a,i0,3(a,l1))") "SIXTRACR> CRCHECK CALLED lout=",lout," restart=",cr_restart," rerun=",cr_rerun," checkp=",cr_checkp
  flush(93)
  ! We are not checkpoint/restart or we have no restart files
  if(.not.cr_checkp) goto 605
  if(.not.cr_pntExist(1).and..not.cr_pntExist(2)) goto 605
  ! If we do we must have a fort.6 as they were created by CRPOINT
  ! NOT TRUE anymore??? We might be NOT rerun but using a Sixin.zip
#ifndef BOINC
  if(.not.cr_rerun) then
    write(lerr,"(a)") "SIXTRACR> ERROR CRCHECK Found "//cr_pntFile(1)//"/"//cr_pntFile(2)//" but NO fort.6"
    call prror(-1)
  endif
#endif
  ! Check at least one restart file is readable
  write(93,"(a)") "SIXTRACR> CRCHECK checking "//cr_pntFile(1)//"/96"
  flush(93)
  if(cr_pntExist(1)) then
    write(93,"(a)") "SIXTRACR> CRCHECK reading "//cr_pntFile(1)//" Record 1 VERSION"
    flush(93)

    rewind(cr_pntUnit(1))
    read(cr_pntUnit(1),err=100,end=100) cr_version,cr_moddate
    if ((cr_version /= version) .or. (cr_moddate /= moddate)) then
      write(93,"(a)") "SIXTRACR> CRCHECK "//cr_pntFile(1)//" was written by SixTrack version="//cr_version//" moddate="//cr_moddate
      write(93,"(a)") "          This is SixTrack version="//version//" moddate="//moddate
      write(93,"(a)") "          Version mismatch; giving up on this file."
      flush(93)
      goto 100
    end if

    write(93,"(a)") "SIXTRACR> CRCHECK reading "//cr_pntFile(1)//" Record 2"
    flush(93)
    read(cr_pntUnit(1),err=100,end=100) crnumlcr,crnuml,crsixrecs,crbinrec, &
      cr_sythck,cril,cr_time,crnapxo,crnapx,cre0,crbetrel,crbrho

    write(93,"(a)") "SIXTRACR> CRCHECK reading "//cr_pntFile(1)//" Record 3"
    flush(93)
    read(cr_pntUnit(1),err=100,end=100) &
      (crbinrecs(j),j=1,(crnapxo+1)/2), &
      (crnumxv(j),j=1,crnapxo),         &
      (crnnumxv(j),j=1,crnapxo),        &
      (crpartID(j),j=1,crnapxo),        &
      (crparentID(j),j=1,crnapxo),      &
      (crpstop(j),j=1,crnapxo),         &
      (crxv(1,j),j=1,crnapxo),          &
      (cryv(1,j),j=1,crnapxo),          &
      (crxv(2,j),j=1,crnapxo),          &
      (cryv(2,j),j=1,crnapxo),          &
      (crsigmv(j),j=1,crnapxo),         &
      (crdpsv(j),j=1,crnapxo),          &
      (crdpsv1(j),j=1,crnapxo),         &
      (crejv(j),j=1,crnapxo),           &
      (crejfv(j),j=1,crnapxo),          &
      (craperv(j,1),j=1,crnapxo),       &
      (craperv(j,2),j=1,crnapxo),       &
      (crllostp(j),j=1,crnapxo)

    write(93,"(a)") "SIXTRACR> CRCHECK Reading "//cr_pntFile(1)//" Record META"
    flush(93)
    call meta_crcheck(cr_pntUnit(1),lerror)
    if(lerror) goto 100

    write(93,"(a)") "SIXTRACR> CRCHECK Reading "//cr_pntFile(1)//" Record 5 DUMP"
    flush(93)
    call dump_crcheck_readdata(cr_pntUnit(1),lerror)
    if (lerror) goto 100

    write(93,"(a)") "SIXTRACR> CRCHECK Reading "//cr_pntFile(1)//" Record 5.5 HION"
    flush(93)
    call hions_crcheck_readdata(cr_pntUnit(1),lerror)
    if (lerror) goto 100

    if (dynk_enabled) then
      write(93,"(a)") "SIXTRACR> CRCHECK Reading "//cr_pntFile(1)//" Record 6 DYNK"
      flush(93)
      call dynk_crcheck_readdata(cr_pntUnit(1),lerror)
      if (lerror) goto 100
    end if

    if(scatter_active) then
      write(93,"(a)") "SIXTRACR> CRCHECK Reading "//cr_pntFile(1)//" Record 7 SCATTER"
      flush(93)
      call scatter_crcheck_readdata(cr_pntUnit(1),lerror)
      if (lerror) goto 100
    end if

    if(limifound) then
      write(93,"(a)") "SIXTRACR> CRCHECK Reading "//cr_pntFile(1)//" Record 8 APERTURE LOSSES FILE"
      flush(93)
      call aper_crcheck_readdata(cr_pntUnit(1),lerror)
      if (lerror) goto 100
    end if

    if(melens .gt. 0) then
      write(93,"(a)") "SIXTRACR> CRCHECK Reading "//cr_pntFile(1)//" Record 9 ELENS"
      flush(93)
      call elens_crcheck(cr_pntUnit(1),lerror)
      if (lerror) goto 100
    endif

    !ERIC new extended checkpoint for synuthck
    if(cr_sythck) then
      !ERICVARS
      ! and make sure we can read the extended vars before leaving fort.95
      ! We will re-read them in crstart to be sure they are restored correctly
      write(93,"(a,i0)") "SIXTRACR> CRCHECK verifying Record 10 extended vars "//cr_pntFile(1)//" crnapxo=",crnapxo
      flush(93)
      read(cr_pntUnit(1),end=100,err=100,iostat=ierro) &
        ((((al(k,m,j,l),l=1,il),j=1,crnapxo),m=1,2),k=1,6), &
        ((((as(k,m,j,l),l=1,il),j=1,crnapxo),m=1,2),k=1,6), &
        (dpd(j),j=1,crnapxo),(dpsq(j),j=1,crnapxo),(fokqv(j),j=1,crnapxo)
      backspace(cr_pntUnit(1),iostat=ierro)
      write(93,"(a)") "SIXTRACR> CRCHECK read "//cr_pntFile(1)//" EXTENDED OK"
      flush(93)
      write(93,"(a)") "SIXTRACR> CRCHECK leaving "//cr_pntFile(1)//" for CRSTART EXTENDED"
      flush(93)
    end if
    cr_pntRead(1) = .true.
    goto 103
  end if
100 continue
  if (.not.cr_pntRead(1)) then
    write(93,"(a)") "SIXTRACR> CRCHECK ERROR Could not read checkpoint file.95"
    flush(93)
  end if
  if (cr_pntExist(2)) then
    write(93,"(a)") "SIXTRACR> CRCHECK Trying "//cr_pntFile(2)//" instead"
    flush(93)
    rewind(cr_pntUnit(2))

    write(93,"(a)") "SIXTRACR> CRCHECK Reading "//cr_pntFile(2)//" Record 1 VERSION"
    flush(93)

    read(cr_pntUnit(2),err=101,end=101) cr_version,cr_moddate
    if ((cr_version /= version) .or. (cr_moddate /= moddate)) then
      write(93,"(a)") "SIXTRACR> CRCHECK "//cr_pntFile(2)//" was written by SixTrack version='"//cr_version//&
        "' moddate='"//cr_moddate//"'"
      write(93,"(a)") "          This is SixTrack version='"//version//"' moddate='"//moddate//"'"
      write(93,"(a)") "          Version mismatch; giving up on this file."
      flush(93)
      goto 101
    end if

    write(93,"(a)") "SIXTRACR> CRCHECK Reading "//cr_pntFile(2)//" Record 2"
    flush(93)
    read(cr_pntUnit(2),err=101,end=101,iostat=ierro) crnumlcr,crnuml,crsixrecs,crbinrec,&
      cr_sythck,cril,cr_time,crnapxo,crnapx,cre0,crbetrel,crbrho
    write(93,"(a)") "SIXTRACR> CRCHECK Reading "//cr_pntFile(2)//" Record 3"
    flush(93)
    read(cr_pntUnit(2),err=101,end=101,iostat=ierro) &
      (crbinrecs(j),j=1,(crnapxo+1)/2),  &
      (crnumxv(j),j=1,crnapxo),          &
      (crnnumxv(j),j=1,crnapxo),         &
      (crpartID(j),j=1,crnapxo),         &
      (crparentID(j),j=1,crnapxo),       &
      (crpstop(j),j=1,crnapxo),          &
      (crxv(1,j),j=1,crnapxo),           &
      (cryv(1,j),j=1,crnapxo),           &
      (crxv(2,j),j=1,crnapxo),           &
      (cryv(2,j),j=1,crnapxo),           &
      (crsigmv(j),j=1,crnapxo),          &
      (crdpsv(j),j=1,crnapxo),           &
      (crdpsv1(j),j=1,crnapxo),          &
      (crejv(j),j=1,crnapxo),            &
      (crejfv(j),j=1,crnapxo),           &
      (craperv(j,1),j=1,crnapxo),        &
      (craperv(j,2),j=1,crnapxo),        &
      (crllostp(j),j=1,crnapxo)

    write(93,"(a)") "SIXTRACR> CRCHECK Reading "//cr_pntFile(2)//" Record META"
    flush(93)
    call meta_crcheck(cr_pntUnit(2),lerror)
    if(lerror) goto 101

    write(93,"(a)") "SIXTRACR> CRCHECK Reading "//cr_pntFile(2)//" Record 5 DUMP"
    flush(93)
    call dump_crcheck_readdata(cr_pntUnit(2),lerror)
    if (lerror) goto 101

    if (dynk_enabled) then
      write(93,"(a)") "SIXTRACR> CRCHECK Reading "//cr_pntFile(2)//" Record 6 DYNK"
      flush(93)
      call dynk_crcheck_readdata(cr_pntUnit(2),lerror)
      if (lerror) goto 101
    end if

    if(scatter_active) then
      write(93,"(a)") "SIXTRACR> CRCHECK Reading "//cr_pntFile(2)//" Record 7 SCATTER"
      flush(93)
      call scatter_crcheck_readdata(cr_pntUnit(2),lerror)
      if (lerror) goto 101
    end if

    if(limifound) then
      write(93,"(a)") "SIXTRACR> CRCHECK Reading "//cr_pntFile(2)//" Record 8 APERTURE LOSSES FILE"
      flush(93)
      call aper_crcheck_readdata(cr_pntUnit(2),lerror)
      if (lerror) goto 101
    end if

    if(melens .gt. 0) then
      write(93,"(a)") "SIXTRACR> CRCHECK Reading "//cr_pntFile(2)//" Record 9 ELENS"
      flush(93)
      call elens_crcheck(cr_pntUnit(2),lerror)
      if (lerror) goto 101
    endif

    !ERIC new extended checkpoint for synuthck
    if(cr_sythck) then
      !ERICVARS
      ! and make sure we can read the extended vars before leaving fort.96
      ! We will re-read them in crstart to be sure they are correct
      write(93,"(a,i0)") "SIXTRACR> CRCHECK verifying Record 10 extended vars "//cr_pntFile(2)//", crnapxo=",crnapxo
      flush(93)
      write(93,"(a)") "SIXTRACR> CRCHECK verifying extended vars "//cr_pntFile(2)
      flush(93)
      read(cr_pntUnit(2),end=101,err=101,iostat=ierro)                 &
        ((((al(k,m,j,l),l=1,il),j=1,crnapxo),m=1,2),k=1,6), &
        ((((as(k,m,j,l),l=1,il),j=1,crnapxo),m=1,2),k=1,6), &
        (dpd(j),j=1,crnapxo),(dpsq(j),j=1,crnapxo),(fokqv(j),j=1,crnapxo)
      backspace(cr_pntUnit(2),iostat=ierro)
      write(93,"(a)") "SIXTRACR> CRCHECK Read "//cr_pntFile(2)//" EXTENDED OK"
      flush(93)
      write(93,"(a)") "SIXTRACR> CRCHECK Leaving "//cr_pntFile(2)//" for CRSTART EXTENDED"
      flush(93)
    end if
    cr_pntRead(2) = .true.
    goto 103
  end if
101 continue
  if (.not.cr_pntRead(2)) then
    write(93,"(a)") "SIXTRACR> CRCHECK ERROR Could not read checkpoint file.96"
    flush(93)
  end if
103 continue

  ! If we have successfully read either fort.95 or fort.96
  ! we need to handle lost particles and ntwin .ne. 2
  ! Otherwise we just continue with checkpointing as requested
  if(cr_pntRead(1) .or. cr_pntRead(2)) then
    write(93,"(2(a,l1),7(a,i0))") "SIXTRACR> CRCHECK read95=",cr_pntRead(1),", read96=",cr_pntRead(2),&
      ", crnapxo=",crnapxo,", crbinrec=",crbinrec,", napx=",napx,", sixrecs=",sixrecs,", crsixrecs=",crsixrecs
#ifndef STF
    write(93,"(a)") "SIXTRACR> CRCHECK crbinrecs:"
    do j=1,(crnapxo+1)/2
      write(93,"(2(a,i0))") "SIXTRACR> ",j,": ",crbinrecs(j)
    end do
#endif
    flush(93)

    ! First we position fort.6 to last checkpoint
    do j=1,crsixrecs
      read(6,"(a1024)",end=604,err=106,iostat=ierro) arecord
      sixrecs = sixrecs+1
    end do
    ! This is not a FLUSH!
    endfile (6,iostat=ierro)
604 continue
    backspace (6,iostat=ierro)
    write(93,"(a,i0)") "SIXTRACR> CRCHECK found fort.6 sixrecs=",sixrecs
    flush(93)

    ! We may be re-running with a DIFFERENT number of turns (numl)
    ! Eric fix this later by reading numl for fort.90
    if (numl /= crnuml) then
      if (numl < crnumlcr) then
        write(lerr,"(2(a,i0))") "SIXTRACR> ERROR New numl < crnumlcr : ",numl," < ",crnumlcr
        write(93,"(2(a,i0))")   "SIXTRACR> ERROR New numl < crnumlcr : ",numl," < ",crnumlcr
        flush(93)
        call prror(-1)
      end if
      write(93,"(2(a,i0))") "SIXTRACR> CRCHECK re-sets numl in binary file headers from ",crnuml," to ",numl
      flush(93)

      ! Reposition binary files fort.90 etc. / singletrackfile.dat
      ! fort.94 = temp file where the data from fort.90 etc. is copied to and then back
#ifdef BOINC
      call boincrf('fort.94',filename)
      open(94,file=filename,form='unformatted',status='unknown')
#else
      open(94,file='fort.94',form='unformatted',status='unknown')
#endif
#ifndef STF
      do ia=1,crnapxo/2,1
        ! First, copy crbinrecs(ia) records of data from fort.91-ia to fort.94
        mybinrecs=0
        binrecs94=0
        myia=91-ia
        !Copy header into integer array hbuff
        read(91-ia,err=105,end=105,iostat=ierro) hbuff
        mybinrecs=mybinrecs+1
        hbuff(51)=numl ! Reset the number of turns (not very elegant)
        write(94,err=105,iostat=ierro) hbuff
        ! Copy particle tracking data
        do j=2,crbinrecs(ia)
          if(ntwin.ne.2) then
            read(91-ia,err=105,end=105,iostat=ierro) (tbuff(k),k=1,17)
            write(94,err=105,iostat=ierro) (tbuff(k),k=1,17)
          else
            read(91-ia,err=105,end=105,iostat=ierro) tbuff
            write(94,err=105,iostat=ierro) tbuff
          endif
          mybinrecs=mybinrecs+1
        end do ! END "do j=2,crbinrecs(ia)"

        ! Second, copy crbinrecs(ia) records of data from fort.94 to fort.91-ia
        rewind(94)
        rewind(91-ia)
        !Copy header
        read(94,err=105,end=105,iostat=ierro) hbuff
        binrecs94=binrecs94+1
        write(91-ia,err=105,iostat=ierro) hbuff
        ! Copy particle tracking data into integer array tbuff
        do j=2,crbinrecs(ia)
          if(ntwin.ne.2) then
            read(94,err=105,end=105,iostat=ierro) (tbuff(k),k=1,17)
            write(91-ia,err=105,iostat=ierro) (tbuff(k),k=1,17)
          else
            read(94,err=105,end=105,iostat=ierro) tbuff
            write(91-ia,err=105,iostat=ierro) tbuff
          endif
          binrecs94=binrecs94+1
        end do ! END "j=2,crbinrecs(ia)"
        !This is not a FLUSH!
        endfile(91-ia,iostat=ierro)
        backspace(91-ia,iostat=ierro)
        rewind(94)
      end do ! END "do ia=1,crnapxo/2,1"
#else
      ! First, copy crbinrecs(ia)*(crnapx/2) records of data from singletrackfile.dat to fort.94
      mybinrecs=0
      !Copy headers
      do ia=1,crnapxo/2,1
        read(90,err=105,end=105,iostat=ierro) hbuff
        mybinrecs=mybinrecs+1
        hbuff(51)=numl ! Reset the number of turns (not very elegant)
        write(94,err=105,iostat=ierro) hbuff
      end do
      ! Copy particle tracking data
      do ia=1,crnapxo/2,1
        do j=2,crbinrecs(ia)
          if(ntwin.ne.2) then
            read(90,err=105,end=105,iostat=ierro) (tbuff(k),k=1,17)
            write(94,err=105,iostat=ierro) (tbuff(k),k=1,17)
          else
            read(90,err=105,end=105,iostat=ierro) tbuff
            write(94,err=105,iostat=ierro) tbuff
          endif
          mybinrecs=mybinrecs+1
        end do
      end do

    ! Second, copy crbinrecs(ia)*(crnapx/2) records of data from fort.94 to singletrackfile.dat
      rewind(94)
      rewind(90)
      binrecs94=0
      ! Copy header
      do ia=1,crnapxo/2,1
        read(94,err=105,end=105,iostat=ierro) hbuff
        binrecs94=binrecs94+1
        write(90,err=105,iostat=ierro) hbuff
      end do
      ! Copy particle tracking data
      do ia=1,crnapxo/2,1
        do j=2,crbinrecs(ia)
          if(ntwin.ne.2) then
            read(94,err=105,end=105,iostat=ierro) (tbuff(k),k=1,17)
            write(90,err=105,iostat=ierro) (tbuff(k),k=1,17)
          else
            read(94,err=105,end=105,iostat=ierro) tbuff
            write(90,err=105,iostat=ierro) tbuff
          endif
          binrecs94=binrecs94+1
        end do
      end do
      !This is not a FLUSH!
      endfile   (90,iostat=ierro)
      backspace (90,iostat=ierro)
#endif
      close(94)
    else !ELSE for "if(nnuml.ne.crnuml) then" -> here we treat nnuml.eq.crnuml, i.e. the number of turns have not been changed
      ! Now with the new array crbinrecs we can ignore files which are
      ! basically finished because a particle has been lost.......
      ! Just check crbinrecs against crbinrec
#ifndef STF
      ! Binary files have been rewritten; now re-position
      write(93,"(a)") "SIXTRACR> CRCHECK re-positioning binary files"
      do ia=1,crnapxo/2,1
        myia=91-ia
        if (crbinrecs(ia).ge.crbinrec) then
          mybinrecs=0
          read(91-ia,err=102,end=102,iostat=ierro) hbuff
          do 11 j=2,crbinrecs(ia)
            if(ntwin.ne.2) then
              read(91-ia,err=102,end=102,iostat=ierro) (tbuff(k),k=1,17)
            else
              read(91-ia,err=102,end=102,iostat=ierro) tbuff
            endif
            mybinrecs=mybinrecs+1
11        continue
          !This is not a FLUSH!
          endfile (91-ia,iostat=ierro)
          backspace (91-ia,iostat=ierro)
        else ! Number of ecords written to this file < general number of records written
            ! => Particle has been lost before last checkpoint, no need to reposition.
          write(93,"(2(a,i0))") "SIXTRACR> CRCHECK ignoring IA ",ia," Unit ",myia
        endif
      end do ! END "do ia=1,crnapxo/2,1"
#else
      mybinrecs=0
      ! Reposition headers
      do ia=1,crnapxo/2,1
        read(90,err=102,end=102,iostat=ierro) hbuff
        mybinrecs=mybinrecs+1
      end do
      !Reposition track records
      do ia=1,crnapxo/2,1
        do j=2,crbinrecs(ia)
          if(ntwin.ne.2) then !ntwin=1
            read(90,err=102,end=102,iostat=ierro) (tbuff(k),k=1,17)
          else                !ntwin=2
            read(90,err=102,end=102,iostat=ierro) tbuff
          endif
          mybinrecs=mybinrecs+1
        end do
      end do
#endif
    end if ! END "if (numl.ne.crnuml) then" and END else

    !reposition dynksets.dat
    if (dynk_enabled .and.(.not.dynk_noDynkSets) ) then
      write(93,"(a)") "SIXTRACR> CRCHECK REPOSITIONING dynksets.dat"
        flush(93)
      call dynk_crcheck_positionFiles
    endif !END if (dynk_enabled .and.(.not.dynk_noDynkSets) )

    !Reposition files for DUMP
    write(93,"(a)") "SIXTRACR> CRCHECK REPOSITIONING DUMP files"
    flush(93)
    call dump_crcheck_positionFiles

    if(scatter_active) then
      write(93,"(a)") "SIXTRACR> CRCHECK REPOSITIONING scatter_log.dat and scatter_summary.dat"
      flush(93)
      call scatter_crcheck_positionFiles
    endif

    if(limifound) then
      write(93,"(a)") "SIXTRACR> CRCHECK REPOSITIONING "//trim(losses_filename)
      flush(93)
      call aper_crcheck_positionFiles
    endif

    ! Set up flag for tracking routines to call CRSTART
    cr_restart = .true.
    write(lout,"(a)") "SIXTRACR> "//repeat("=",80)
    write(lout,"(a)") "SIXTRACR>  Restarted"
    write(lout,"(a)") "SIXTRACR> "//repeat("=",80)
    !Flush or truncate?
    endfile(lout,iostat=ierro)
    backspace(lout,iostat=ierro)
    write(93,"(a,i0)") "SIXTRACR> CRCHECK restart=TRUE, crnumlcr=",crnumlcr
    flush(93)
    return
  end if

  goto 605                  !Should not end up here -> checkpoint failed.
                            ! Start simulation over!

  !--   Just abort if we cannot re-position/copy the binary files,
#ifndef STF
102 continue
  write(lerr,"(2(a,i0))") "SIXTRACR> ERROR PROBLEMS RE-READING fort.",myia," IOSTAT=",ierro
  write(lerr,"(3(a,i0))") "          Unit ",myia," mybinrecs=",mybinrecs," Expected crbinrecs=",crbinrecs(ia)
  write(lerr,"(a)")       "SIXTRACR> CRCHECK failure positioning binary files"
  call prror(-1)
105 continue
  write(lerr,"(2(a,i0))") "SIXTRACR> ERROR PROBLEMS COPYING fort.",myia," IOSTAT=",ierro
  write(lerr,"(4(a,i0))") "          Unit ",myia," mybinrecs=",mybinrecs,&
    " Expected crbinrecs=",crbinrecs(ia)," binrecs94=",binrecs94
  write(lerr,"(a)")       "SIXTRACR> CRCHECK failure copying binary files"
  call prror(-1)
#else
102 continue
  write(lerr,"(2(a,i0))") "SIXTRACR> ERROR PROBLEMS RE-READING singletrackfile.dat for ia=",ia," IOSTAT=",ierro
  write(lerr,"(2(a,i0))") "          mybinrecs=",mybinrecs," Expected crbinrecs=",crbinrecs(ia)
  write(lerr,"(a)")       "SIXTRACR> CRCHECK failure positioning binary files"
  call prror(-1)
105 continue
  write(lerr,"(2(a,i0))") "SIXTRACR> ERROR PROBLEMS COPYING particle pair ",ia," IOSTAT=",ierro," from/to singletrackfile.dat"
  write(lerr,"(3(a,i0))") "          mybinrecs=",mybinrecs," Expected crbinrecs=",crbinrecs(ia)," binrecs94=",binrecs94
  write(lerr,"(a)")       "SIXTRACR> CRCHECK failure copying binary files"
  call prror(-1)
#endif
  ! We are not checkpointing or we have no checkpoints
  ! or we have no readable checkpoint
  ! If not checkpointing we can just give up on lout and use
  ! fort.6. We don't need to count records at all
605 continue
  write(93,"(a,l1)") "SIXTRACR> CRCHECK no restart possible checkp=",cr_checkp
  flush(93)
  if(.not.cr_checkp) then
    if(cr_rerun) then
      ! we nevertheless have an existing fort.6
      ! we will just overwrite it for now and delete
      ! 92 to avoid abend copying it again
      write(93,"(a)") "SIXTRACR> CRCHECK overwriting fort.6"
      flush(93)
    end if
    ! and just use fort.6 from now on
    write(93,"(a)") "SIXTRACR> CRCHECK giving up on LOUT"
    flush(93)
    ! Copy the lout to fort.6 (the file, not output_unit)
    ! It seems that FORTRAN will open the file automatically?
    ! There are no open(unit=6) etc. calls anywhere...
    rewind(lout)
3   read(lout,'(a1024)',end=1,err=107,iostat=ierro) arecord
    lstring=1024
    do i=1024,2,-1
      lstring=i
      if (arecord(i:i).ne.' ')goto 2
      lstring=lstring-1
    enddo
2   write(6,'(a)') arecord(1:lstring)
    goto 3
    ! Not a flush?
1   endfile(6,iostat=ierro)
    backspace(6,iostat=ierro)
    ! This is not a FLUSH!
    rewind(lout)
    endfile(lout,iostat=ierro)
    close(lout)
    lout=6
  endif

  return

106 continue
  write(93,"(3(a,i0))") "SIXTRACR> ERROR reading fort.6, iostat = ",ierro,", sixrecs = ",sixrecs,", crsixrecs = ",crsixrecs
  flush(93)
  write(lerr,"(a)") "SIXTRACR> ERROR CRCHECK Failure positioning fort.6"
  call prror

107 continue
  write(93,"(a,i0)") "SIXTRACR> ERROR reading "//cr_outFile//", iostat=",ierro
  flush(93)
  write(lerr,"(a)") "SIXTRACR> ERROR CRCHECK Failure positioning "//cr_outFile
  call prror

end subroutine crcheck

! ================================================================================================ !
!  This subroutine writes the checkpoint data to the binary checpoint files
!  Last modified: 2019-04-29
! ================================================================================================ !
subroutine crpoint

  use floatPrecision
  use numerical_constants

  use dynk,    only : dynk_enabled,dynk_getvalue,dynk_fSets_cr,dynk_cSets_unique,dynk_nSets_unique,dynk_filePos,dynk_crpoint
  use dump,    only : dump_crpoint
  use aperture,only : aper_crpoint,limifound
  use scatter, only : scatter_active, scatter_crpoint
  use elens,   only : melens, elens_crpoint

  use crcoall
  use parpro
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da
  use mod_hions
  use mod_version
  use mod_time
  use mod_meta
  use mod_units
  use mod_settings

  integer j,l,k,m,nPoint
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
      write(93,"(a)") "CR_POINT> Filling dynk_fSets_cr"
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
    if(st_debug) then
      write(93,"(a)") "CR_POINT> Writing to checkpoint file "//cr_pntFile(nPoint)
      write(93,"(a)") "CR_POINT>  * Tracking variables"
      write(93,"(a)") "CR_POINT>  * Particle arrays"
      flush(93)
    end if
    rewind(cr_pntUnit(nPoint))

    write(cr_pntUnit(nPoint),err=100,iostat=ierro) version, moddate
    write(cr_pntUnit(nPoint),err=100,iostat=ierro) crnumlcr, numl, sixrecs, binrec, sythckcr, il,   &
      time3, napxo, napx, e0, betrel, brho
    write(cr_pntUnit(nPoint),err=100,iostat=ierro) &
      (binrecs(j),j=1,(napxo+1)/2), &
      (numxv(j),j=1,napxo),         &
      (nnumxv(j),j=1,napxo),        &
      (partID(j),j=1,napxo),        &
      (parentID(j),j=1,napxo),      &
      (pstop(j),j=1,napxo),         &
      (xv1(j),j=1,napxo),           &
      (yv1(j),j=1,napxo),           &
      (xv2(j),j=1,napxo),           &
      (yv2(j),j=1,napxo),           &
      (sigmv(j),j=1,napxo),         &
      (dpsv(j),j=1,napxo),          &
      (dpsv1(j),j=1,napxo),         &
      (ejv(j),j=1,napxo),           &
      (ejfv(j),j=1,napxo),          &
      (aperv(j,1),j=1,napxo),       &
      (aperv(j,2),j=1,napxo),       &
      (llostp(j),j=1,napxo)
    flush(cr_pntUnit(nPoint))

    if(st_debug) write(93,"(a)") "CR_POINT>  * META variables"
    call meta_crpoint(cr_pntUnit(nPoint),wErr,ierro)
    if(wErr) goto 100

    if(st_debug) write(93,"(a)") "CR_POINT>  * DUMP variables"
    call dump_crpoint(cr_pntUnit(nPoint), wErr,ierro)
    if(wErr) goto 100

    if(st_debug) write(93,"(a)") "CR_POINT>  * HION variables"
    call hions_crpoint(cr_pntUnit(nPoint),wErr,ierro)
    if(wErr) goto 100

    if(dynk_enabled) then
      if(st_debug) write(93,"(a)") "CR_POINT>  * DYNK variables"
      call dynk_crpoint(cr_pntUnit(nPoint),wErr,ierro)
      if(wErr) goto 100
    end if

    if(scatter_active) then
      if(st_debug) write(93,"(a)") "CR_POINT>  * SCATTER variables"
      call scatter_crpoint(cr_pntUnit(nPoint),wErr,ierro)
      if(wErr) goto 100
    end if

    if(limifound) then
      if(st_debug) write(93,"(a)") "CR_POINT>  * APERTURE variables"
      call aper_crpoint(cr_pntUnit(nPoint),wErr,ierro)
      if(wErr) goto 100
    end if

    if(melens > 0) then
      if(st_debug) write(93,"(a)") "CR_POINT>  * ELENS variables"
      call elens_crpoint(cr_pntUnit(nPoint),wErr,ierro)
      if(wErr) goto 100
    end if

    if(sythckcr) then
      if(ithick == 1) then
        if(st_debug) write(93,"(a)") "CR_POINT>  * THICK EXTENDED arrays"
        write(cr_pntUnit(nPoint),err=100,iostat=ierro) &
          ((((al(k,m,j,l),l=1,il),j=1,napxo),m=1,2),k=1,6), &
          ((((as(k,m,j,l),l=1,il),j=1,napxo),m=1,2),k=1,6)
        flush(cr_pntUnit(nPoint))
      end if

      if(st_debug) write(93,"(a)") "CR_POINT>  * THICK arrays"
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
!  If we are restarting (cr_restart is TRUE), this routine is called in the beginning of the tracking
!  loops. It is used to copy the cr* variables to the normal variables, e.g. crnapx -> napx etc.
!  The file fort.93 is used as a log file for the checkpoint/restarting.
! ================================================================================================ !
subroutine crstart

  use floatPrecision
  use numerical_constants
  use dynk, only : dynk_enabled, dynk_crstart
  use scatter, only: scatter_active, scatter_crstart
  use elens,   only : melens, elens_crstart

  use crcoall
  use parpro
  use mod_common
  use mod_common_main
  use mod_commons
  use mod_common_track
  use mod_common_da
  use mod_meta
  use mod_alloc
  use mod_hions
  use mod_units

  implicit none

  real(kind=fPrec) dynk_newValue
  logical fErr
  integer j,l,k,m,i
  character(len=256) filename

  save

  write(93,"(a,i0)") "SIXTRACR> CRSTART called crnumlcr = ",crnumlcr
  flush(93)
  cr_numl = crnumlcr

  ! We do NOT reset numl so that a run can be extended for
  ! for more turns from the last checkpoint
  ! but we need to worry about numxv, nnumxv
  binrec   = crbinrec
  sythckcr = cr_sythck

  ! the cr_time is required (crtime0/1 removed)
  napxo=crnapxo
  napx=crnapx
  e0=cre0
  e0f=sqrt(e0**2-nucm0**2)
  betrel=crbetrel
  brho=crbrho

  write(93,"(a)") "SIXTRACR> CRSTART doing binrecs"
  flush(93)

  do j=1,(napxo+1)/2
    binrecs(j)=crbinrecs(j)
  end do

  write(93,"(a)") "SIXTRACR> CRSTART doing normal NPART vars"
  flush(93)
  do j=1,napxo
    numxv(j)=crnumxv(j)
    nnumxv(j)=crnnumxv(j)
    partID(j)=crpartID(j)
    parentID(j)=crparentID(j)
    pstop(j)=crpstop(j)
    llostp(j)=crllostp(j)
    xv1(j)=crxv(1,j)
    yv1(j)=cryv(1,j)
    xv2(j)=crxv(2,j)
    yv2(j)=cryv(2,j)
    sigmv(j)=crsigmv(j)
    dpsv(j)=crdpsv(j)
    dpsv1(j)=crdpsv1(j)
    ! TEMPORARY? fix for crabamp/multipole problem
    !       oidpsv(j)=croidpsv(j)
    oidpsv(j)=one/(one+dpsv(j))
    moidpsv(j)=mtc(j)/(one+dpsv(j))
    omoidpsv(j)=c1e3*((one-mtc(j))*oidpsv(j))
    ejv(j)=crejv(j)
    ejfv(j)=crejfv(j)
    rvv(j)=(ejv(j)*e0f)/(e0*ejfv(j))
    aperv(j,1)=craperv(j,1)
    aperv(j,2)=craperv(j,2)
    if(pstop(j) .eqv. .false.) then
      numxv(j)=numl
      nnumxv(j)=numl
    endif
  end do

  !ERIC new extended checkpoint for synuthck
  call meta_crstart
  if(dynk_enabled) call dynk_crstart
  if(scatter_active) call scatter_crstart
  call hions_crstart
  if(melens .gt. 0) call elens_crstart

  if(cr_sythck) then
    !ERICVARS now read the extended vars from fort.95/96.
    if(cril /= il) then
      write(lout,"(2(a,i0))") "SIXTRACR> CRSTART Problem as cril/il are different cril = ",cril,", il = ",il
      write(93,  "(2(a,i0))") "SIXTRACR> CRSTART Problem as cril/il are different cril = ",cril,", il = ",il
      flush(93)
      write(lerr,"(a)") "SIXTRACR> ERROR CRSTART Problem wih cril/il extended C/R"
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
      write(93,"(a)") "SIXTRACR> CRSTART Read "//cr_pntFile(1)//" EXTENDED OK"
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

      write(93,"(a)") "SIXTRACR> CRSTART Read "//cr_pntFile(2)//" EXTENDED OK"
      flush(93)
      goto 102
    end if
100 continue
    write(93,"(a,i0)") "SIXTRACR> CRSTART Could not read checkpoint file "//cr_pntFile(1)//" (extended), iostat = ",ierro
    goto 103
101 continue
    write(93,"(a,i0)") "SIXTRACR> CRSTART Could not read checkpoint file "//cr_pntFile(2)//" (extended), iostat = ",ierro
103 continue
    flush(93)
    write(lerr,"(a)") "SIXTRACR> ERROR CRSTART Problem with extended checkpoint"
    call prror
  end if

102 continue
  write(93,"(3(a,i0))") "SIXTRACR> CRSTART sixrecs = ",sixrecs,", crsixrecs = ",crsixrecs,", binrec = ",binrec
  flush(93)

  ! Just throw away our fort.92 stuff.
  rewind(lout)
  endfile(lout,iostat=ierro)
  close(lout)

  call f_open(unit=lout,file=cr_outFile,formatted=.true.,mode="rw",err=fErr)
  ! but also add the rerun message
  write(lout,"(a)") "SIXTRACR> "//repeat("=",80)
  write(lout,"(a)") "SIXTRACR>  Restarted"
  write(lout,"(a)") "SIXTRACR> "//repeat("=",80)
  endfile(lout,iostat=ierro)
  backspace(lout,iostat=ierro)

  return

606 continue
  backspace(6,iostat=ierro)
  write(lout,"(2(a,i0))") "SIXTRACR> CRSTART Problem re-positioning fort.6: sixrecs = ",sixrecs,", crsixrecs = ",crsixrecs
  call prror

end subroutine crstart

subroutine cr_copyOut

  use crcoall
  use mod_units
  use, intrinsic :: iso_fortran_env, only : output_unit

  character(len=1024) inLine
  integer ioStat, lnSize, nLines

  if(lout == output_unit) return

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
