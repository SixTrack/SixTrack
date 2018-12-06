
module checkpoint_restart

  use floatPrecision
  use parpro

  implicit none

  real,                public, save :: crtime3
  real(kind=fPrec),    public, save :: cre0

  character(len=1024), public, save :: arecord
  character(len=20),   public, save :: stxt
  character(len=80),   public, save :: runtim

  integer,             public, save :: crsixrecs
  integer,             public, save :: crbinrec
  integer,             public, save :: crbnlrec
  integer,             public, save :: crbllrec
  integer,             public, save :: cril
  integer,             public, save :: crnumlcr
  integer,             public, save :: crnuml
  integer,             public, save :: crnapxo
  integer,             public, save :: crnapx
  integer,             public, save :: binrec
  integer,             public, save :: bnlrec
  integer,             public, save :: bllrec
  integer,             public, save :: numlcr
  integer,             public, save :: sixrecs

  logical,             public, save :: rerun
  logical,             public, save :: start
  logical,             public, save :: restart
  logical,             public, save :: checkp
  logical,             public, save :: fort95
  logical,             public, save :: fort96
  logical,             public, save :: read95
  logical,             public, save :: read96
  logical,             public, save :: crsythck

  real(kind=fPrec),    allocatable, public, save :: crxv(:,:)    ! (2,npart)
  real(kind=fPrec),    allocatable, public, save :: cryv(:,:)    ! (2,npart)
  real(kind=fPrec),    allocatable, public, save :: crsigmv(:)   ! (npart)
  real(kind=fPrec),    allocatable, public, save :: crdpsv(:)    ! (npart)
  real(kind=fPrec),    allocatable, public, save :: crdpsv1(:)   ! (npart)
  real(kind=fPrec),    allocatable, public, save :: crejv(:)     ! (npart)
  real(kind=fPrec),    allocatable, public, save :: crejfv(:)    ! (npart)
  real(kind=fPrec),    allocatable, public, save :: craperv(:,:) ! (npart,2)
  real(kind=fPrec),    allocatable, public, save :: crxvl(:,:)   ! (2,npart)
  real(kind=fPrec),    allocatable, public, save :: cryvl(:,:)   ! (2,npart)
  real(kind=fPrec),    allocatable, public, save :: crdpsvl(:)   ! (npart)
  real(kind=fPrec),    allocatable, public, save :: crejvl(:)    ! (npart)
  real(kind=fPrec),    allocatable, public, save :: crsigmvl(:)  ! (npart)

  integer,             allocatable, public, save :: binrecs(:)   ! ((npart+1)/2)
  integer,             allocatable, public, save :: crbinrecs(:) ! (npart+1)/2)
  integer,             allocatable, public, save :: crnumxv(:)   ! (npart)
  integer,             allocatable, public, save :: crnnumxv(:)  ! (npart)
  integer,             allocatable, public, save :: crpartID(:)  ! (npart)
  integer,             allocatable, public, save :: crparentID(:) ! (npart)

  logical,             allocatable, public, save :: crpstop(:)   ! (npart)

  integer,                          public, save :: crnpart_old = -1

  ! Others. Keep length in sync with includes/version.f90
  character(len=8),    public, save :: cr_version
  character(len=10),   public, save :: cr_moddate

contains

! ================================================================================================ !
!  CR ALLOCATE ARRAYS
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-06-22
! ================================================================================================ !
subroutine cr_expand_arrays(npart_new)

  use mod_alloc
  use numerical_constants, only : zero

  implicit none

  integer, intent(in) :: npart_new

  integer :: npair_new
  npair_new = npart_new/2 + 1

  call alloc(crxv,      2, npart_new,      zero,    "crxv")
  call alloc(cryv,      2, npart_new,      zero,    "cryv")
  call alloc(crsigmv,      npart_new,      zero,    "crsigmv")
  call alloc(crdpsv,       npart_new,      zero,    "crdpsv")
  call alloc(crdpsv1,      npart_new,      zero,    "crdpsv1")
  call alloc(crejv,        npart_new,      zero,    "crejv")
  call alloc(crejfv,       npart_new,      zero,    "crejfv")
  call alloc(craperv,      npart_new, 2,   zero,    "craperv")
  call alloc(crxvl,     2, npart_new,      zero,    "crxvl")
  call alloc(cryvl,     2, npart_new,      zero,    "cryvl")
  call alloc(crdpsvl,      npart_new,      zero,    "crdpsvl")
  call alloc(crejvl,       npart_new,      zero,    "crejvl")
  call alloc(crsigmvl,     npart_new,      zero,    "crsigmvl")
  call alloc(binrecs,      npair_new,      0,       "binrecs")
  call alloc(crbinrecs,    npair_new,      0,       "crbinrecs")
  call alloc(crnumxv,      npart_new,      0,       "crnumxv")
  call alloc(crnnumxv,     npart_new,      0,       "crnnumxv")
  call alloc(crpartID,     npart_new,      0,       "crpartID")
  call alloc(crparentID,   npart_new,      0,       "crparentID")
  call alloc(crpstop,      npart_new,      .false., "crpstop")

  crnpart_old = npart_new

end subroutine cr_expand_arrays

! ================================================================================================ !
!  CRCHECK
!  Last modified: 2018-06-12
!
!  This subroutine checks if the C/R files fort.95 and fort.96 exists, and if so tries to load
!  them into the cr* variables.
!  This routine also repositions the output files for fort.90..91-napx/2 or STF, DUMP, and DYNK.
!
!  The file fort.93 is used as a log file for the checkpoint/restarting.
! ================================================================================================ !
subroutine crcheck

  use floatPrecision
  use string_tools
  use numerical_constants
  use dynk,    only : dynk_enabled, dynk_noDynkSets,dynk_crcheck_readdata,dynk_crcheck_positionFiles
  use dump,    only : dump_crcheck_readdata,dump_crcheck_positionFiles
  use scatter, only : scatter_active,scatter_crcheck_readdata,scatter_crcheck_positionFiles
  use, intrinsic :: iso_fortran_env, only : int32
  use crcoall
  use parpro
  use mod_common
  use mod_commonmn
  use mod_commons
  use mod_commont
  use mod_commond
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

#ifdef BOINC
  character(len=256) filename
#endif
  save
  !--------------------------------------------------------------------!
  ! Check the size of CR arrays
  if(npart /= crnpart_old) call cr_expand_arrays(npart)

  restart = .false.
  read95  = .false.
  read96  = .false.
  ! Some log entries to fort.93
  write(93,"(a,i0,3(a,l1))") "SIXTRACR> CRCHECK CALLED lout=",lout," restart=",restart," rerun=",rerun," checkp=",checkp
  flush(93)
  ! We are not checkpoint/restart or we have no restart files
  if (.not.checkp) goto 605
  if (.not.fort95.and..not.fort96) goto 605
  ! If we do we must have a fort.6 as they were created by CRPOINT
  ! NOT TRUE anymore??? We might be NOT rerun but using a Sixin.zip
#ifndef BOINC
  if (.not.rerun) then
    write(lout,"(a)") "SIXTRACR> ERROR CRCHECK Found fort.95/fort.96 but NO fort.6"
    call prror(-1)
  endif
#endif
  ! Check at least one restart file is readable
  write(93,"(a)") "SIXTRACR> CRCHECK checking fort.95/96"
  flush(93)
  if (fort95) then
    write(93,"(a)") "SIXTRACR> CRCHECK reading fort.95 Record 1 VERSION"
    flush(93)

    rewind 95
    read(95,err=100,end=100) cr_version,cr_moddate
    if ((cr_version /= version) .or. (cr_moddate /= moddate)) then
      write(93,"(a)") "SIXTRACR> CRCHECK: fort.95 was written by SixTrack version="//cr_version//" moddate="//cr_moddate
      write(93,"(a)") "          This is SixTrack version="//version//" moddate="//moddate
      write(93,"(a)") "          Version mismatch; giving up on this file."
      flush(93)
      goto 100
    end if

    write(93,"(a)") "SIXTRACR> CRCHECK reading fort.95 Record 2"
    flush(93)
    read(95,err=100,end=100) crnumlcr,crnuml,crsixrecs,crbinrec,crbnlrec,crbllrec,crsythck,cril,crtime3,crnapxo,crnapx,cre0

    write(93,"(a)") "SIXTRACR> CRCHECK reading fort.95 Record 3"
    flush(93)
    read(95,err=100,end=100) &
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
      (crxvl(1,j),j=1,crnapxo),         &
      (crxvl(2,j),j=1,crnapxo),         &
      (cryvl(1,j),j=1,crnapxo),         &
      (cryvl(2,j),j=1,crnapxo),         &
      (crdpsvl(j),j=1,crnapxo),         &
      (crejvl(j),j=1,crnapxo),          &
      (crsigmvl(j),j=1,crnapxo)

    write(93,"(a)") "SIXTRACR> CRCHECK reading fort.95 Record META"
    flush(93)
    call meta_crcheck(95,lerror)
    if(lerror) goto 100

    write(93,"(a)") "SIXTRACR> CRCHECK reading fort.95 Record 5 DUMP"
    flush(93)
    call dump_crcheck_readdata(95,lerror)
    if (lerror) goto 100

    write(93,"(a)") "SIXTRACR> CRCHECK reading fort.95 Record 5.5 HION"
    flush(93)
    call hions_crcheck_readdata(95,lerror)
    if (lerror) goto 100

    if (dynk_enabled) then
      write(93,"(a)") "SIXTRACR> CRCHECK reading fort.95 Record 6 DYNK"
      flush(93)
      call dynk_crcheck_readdata(95,lerror)
      if (lerror) goto 100
    end if

    if(scatter_active) then
      write(93,"(a)") "SIXTRACR> CRCHECK reading fort.95 Record 7 SCATTER"
      flush(93)
      call scatter_crcheck_readdata(95,lerror)
      if (lerror) goto 100
    end if

    !ERIC new extended checkpoint for synuthck
    if (crsythck) then
      !ERICVARS
      ! and make sure we can read the extended vars before leaving fort.95
      ! We will re-read them in crstart to be sure they are restored correctly
      write(93,"(a,i0)") "SIXTRACR> CRCHECK verifying Record 8 extended vars fort.95 crnapxo=",crnapxo
      flush(93)
      read(95,end=100,err=100,iostat=ierro) &
        ((((al(k,m,j,l),l=1,il),j=1,crnapxo),m=1,2),k=1,6), &
        ((((as(k,m,j,l),l=1,il),j=1,crnapxo),m=1,2),k=1,6), &
        (aek(j),j=1,crnapxo),                               &
        (afok(j),j=1,crnapxo),                              &
        (as3(j),j=1,crnapxo),                               &
        (as4(j),j=1,crnapxo),                               &
        (as6(j),j=1,crnapxo),                               &
        (co(j),j=1,crnapxo),                                &
        (dpd(j),j=1,crnapxo),                               &
        (dpsq(j),j=1,crnapxo),                              &
        (fi(j),j=1,crnapxo),                                &
        (fok(j),j=1,crnapxo),                               &
        (fok1(j),j=1,crnapxo),                              &
        (fokqv(j),j=1,crnapxo),                             &
        (g(j),j=1,crnapxo),                                 &
        (gl(j),j=1,crnapxo),                                &
        (hc(j),j=1,crnapxo),                                &
        (hi(j),j=1,crnapxo),                                &
        (hi1(j),j=1,crnapxo),                               &
        (hm(j),j=1,crnapxo),                                &
        (hp(j),j=1,crnapxo),                                &
        (hs(j),j=1,crnapxo),                                &
        (rho(j),j=1,crnapxo),                               &
        (rhoc(j),j=1,crnapxo),                              &
        (rhoi(j),j=1,crnapxo),                              &
        (si(j),j=1,crnapxo),                                &
        (siq(j),j=1,crnapxo),                               &
        (sm1(j),j=1,crnapxo),                               &
        (sm12(j),j=1,crnapxo),                              &
        (sm2(j),j=1,crnapxo),                               &
        (sm23(j),j=1,crnapxo),                              &
        (sm3(j),j=1,crnapxo),                               &
        (wf(j),j=1,crnapxo),                                &
        (wfa(j),j=1,crnapxo),                               &
        (wfhi(j),j=1,crnapxo)
      backspace (95,iostat=ierro)
      write(93,"(a)") "SIXTRACR> CRCHECK read fort.95 EXTENDED OK"
      flush(93)
      write(93,"(a)") "SIXTRACR> CRCHECK leaving fort.95 for CRSTART EXTENDED"
      flush(93)
    end if
    read95=.true.
    goto 103
  end if
100 continue
  if (.not.read95) then
    write(93,"(a)") "SIXTRACR> CRCHECK COULD NOT READ CHECKPOINT FILE 95"
    flush(93)
  end if
  if (fort96) then
    write(93,"(a)") "SIXTRACR> CRCHECK trying fort.96 instead"
    flush(93)
    rewind 96

    write(93,"(a)") "SIXTRACR> CRCHECK reading fort.96 Record 1 VERSION"
    flush(93)

    read(96,err=101,end=101) cr_version,cr_moddate
    if ((cr_version /= version) .or. (cr_moddate /= moddate)) then
      write(93,"(a)") "SIXTRACR> CRCHECK: fort.96 was written by SixTrack version='"//cr_version//"' moddate='"//cr_moddate//"'"
      write(93,"(a)") "          This is SixTrack version='"//version//"' moddate='"//moddate//"'"
      write(93,"(a)") "          Version mismatch; giving up on this file."
      flush(93)
      goto 101
    end if

    write(93,"(a)") "SIXTRACR> CRCHECK reading fort.96 Record 2"
    flush(93)
    read(96,err=101,end=101,iostat=ierro) crnumlcr,crnuml,crsixrecs,crbinrec,crbnlrec,crbllrec,&
      crsythck,cril,crtime3,crnapxo,crnapx,cre0
    write(93,"(a)") "SIXTRACR> CRCHECK reading fort.96 Record 3"
    flush(93)
    read(96,err=101,end=101,iostat=ierro) &
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
      (crxvl(1,j),j=1,crnapxo),          &
      (crxvl(2,j),j=1,crnapxo),          &
      (cryvl(1,j),j=1,crnapxo),          &
      (cryvl(2,j),j=1,crnapxo),          &
      (crdpsvl(j),j=1,crnapxo),          &
      (crejvl(j),j=1,crnapxo),           &
      (crsigmvl(j),j=1,crnapxo)

    write(93,"(a)") "SIXTRACR> CRCHECK reading fort.96 Record META"
    flush(93)
    call meta_crcheck(96,lerror)
    if(lerror) goto 101

    write(93,"(a)") "SIXTRACR> CRCHECK reading fort.96 Record 5 DUMP"
    flush(93)
    call dump_crcheck_readdata(96,lerror)
    if (lerror) goto 101

    if (dynk_enabled) then
      write(93,"(a)") "SIXTRACR> CRCHECK reading fort.96 Record 6 DYNK"
      flush(93)
      call dynk_crcheck_readdata(96,lerror)
      if (lerror) goto 101
    end if

    if(scatter_active) then
      write(93,"(a)") "SIXTRACR> CRCHECK reading fort.96 Record 7 SCATTER"
      flush(93)
      call scatter_crcheck_readdata(96,lerror)
      if (lerror) goto 101
    end if

    !ERIC new extended checkpoint for synuthck
    if (crsythck) then
      !ERICVARS
      ! and make sure we can read the extended vars before leaving fort.96
      ! We will re-read them in crstart to be sure they are correct
      write(93,"(a,i0)") "SIXTRACR> CRCHECK verifying Record 8 extended vars fort.96, crnapxo=",crnapxo
      flush(93)
      write(93,"(a)") "SIXTRACR> CRCHECK verifying extended vars fort.96"
      flush(93)
      read(96,end=101,err=101,iostat=ierro)                  &
        ((((al(k,m,j,l),l=1,il),j=1,crnapxo),m=1,2),k=1,6), &
        ((((as(k,m,j,l),l=1,il),j=1,crnapxo),m=1,2),k=1,6), &
        (aek(j),j=1,crnapxo),                               &
        (afok(j),j=1,crnapxo),                              &
        (as3(j),j=1,crnapxo),                               &
        (as4(j),j=1,crnapxo),                               &
        (as6(j),j=1,crnapxo),                               &
        (co(j),j=1,crnapxo),                                &
        (dpd(j),j=1,crnapxo),                               &
        (dpsq(j),j=1,crnapxo),                              &
        (fi(j),j=1,crnapxo),                                &
        (fok(j),j=1,crnapxo),                               &
        (fok1(j),j=1,crnapxo),                              &
        (fokqv(j),j=1,crnapxo),                             &
        (g(j),j=1,crnapxo),                                 &
        (gl(j),j=1,crnapxo),                                &
        (hc(j),j=1,crnapxo),                                &
        (hi(j),j=1,crnapxo),                                &
        (hi1(j),j=1,crnapxo),                               &
        (hm(j),j=1,crnapxo),                                &
        (hp(j),j=1,crnapxo),                                &
        (hs(j),j=1,crnapxo),                                &
        (rho(j),j=1,crnapxo),                               &
        (rhoc(j),j=1,crnapxo),                              &
        (rhoi(j),j=1,crnapxo),                              &
        (si(j),j=1,crnapxo),                                &
        (siq(j),j=1,crnapxo),                               &
        (sm1(j),j=1,crnapxo),                               &
        (sm12(j),j=1,crnapxo),                              &
        (sm2(j),j=1,crnapxo),                               &
        (sm23(j),j=1,crnapxo),                              &
        (sm3(j),j=1,crnapxo),                               &
        (wf(j),j=1,crnapxo),                                &
        (wfa(j),j=1,crnapxo),                               &
        (wfhi(j),j=1,crnapxo)
      backspace (96,iostat=ierro)
      write(93,"(a)") "SIXTRACR> CRCHECK read fort.96 EXTENDED OK"
      flush(93)
      write(93,"(a)") "SIXTRACR> CRCHECK, leaving fort.96 for CRSTART EXTENDED"
      flush(93)
    end if
    read96=.true.
    goto 103
  end if
101 continue
  if (.not.read96) then
    write(93,"(a)") "SIXTRACR> CRCHECK, COULD NOT READ CHECKPOINT FILE 96"
    flush(93)
  end if
103 continue

  ! If we have successfully read either fort.95 or fort.96
  ! we need to handle lost particles and ntwin .ne. 2
  ! Otherwise we just continue with checkpointing as requested
  if (read95.or.read96) then
    write(93,"(2(a,l1),7(a,i0))") "SIXTRACR> CRCHECK read95=",read95," read96=",read96,&
      " crnapxo=",crnapxo," crbinrec=",crbinrec," napx=",napx," sixrecs=",sixrecs,     &
      " crsixrecs=",crsixrecs," crbnlrec=",crbnlrec," crbllrec=",crbllrec
    write(93,"(a)") "SIXTRACR> CRCHECK crbinrecs:"
    do j=1,(crnapxo+1)/2
      write(93,"(2(a,i0))") "SIXTRACR> ",j,": ",crbinrecs(j)
    end do
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
        write(lout,"(2(a,i0))") "SIXTRACR> ERROR New numl < crnumlcr : ",numl," < ",crnumlcr
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
        rewind 94
        rewind 91-ia
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
        endfile (91-ia,iostat=ierro)
        backspace (91-ia,iostat=ierro)
        rewind 94
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
      rewind 94
      rewind 90
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
      write(93,"(a)") "SIXTRACR> CRCHECK REPOSITIONING scatter_log.txt"
      flush(93)
      call scatter_crcheck_positionFiles
    endif

    ! Set up flag for tracking routines to call CRSTART
    restart=.true.
    write(lout,"(a80)") runtim
    !Flush or truncate?
    endfile (lout,iostat=ierro)
    backspace (lout,iostat=ierro)
    write(93,"(a,i0)") "SIXTRACR> CRCHECK restart=TRUE',' crnumlcr=",crnumlcr
    flush(93)
    return
  end if

  goto 605                  !Should not end up here -> checkpoint failed.
                            ! Start simulation over!

  !--   Just abort if we cannot re-position/copy the binary files,
#ifndef STF
102 continue
  write(lout,"(a)") ""
  write(lout,"(2(a,i0))") "SIXTRACR> ERROR PROBLEMS RE-READING fort.",myia," IOSTAT=",ierro
  write(lout,"(3(a,i0))") "          Unit ",myia," mybinrecs=",mybinrecs," Expected crbinrecs=",crbinrecs(ia)
  write(lout,"(a)")       "SIXTRACR> CRCHECK failure positioning binary files"
  call prror(-1)
105 continue
  write(lout,"(a)") ""
  write(lout,"(2(a,i0))") "SIXTRACR> ERROR PROBLEMS COPYING fort.",myia," IOSTAT=",ierro
  write(lout,"(4(a,i0))") "          Unit ",myia," mybinrecs=",mybinrecs,&
    " Expected crbinrecs=",crbinrecs(ia)," binrecs94=",binrecs94
  write(lout,"(a)")       "SIXTRACR> CRCHECK failure copying binary files"
  call prror(-1)
#else
102 continue
  write(lout,"(a)") ""
  write(lout,"(2(a,i0))") "SIXTRACR> ERROR PROBLEMS RE-READING singletrackfile.dat for ia=",ia," IOSTAT=",ierro
  write(lout,"(2(a,i0))") "          mybinrecs=",mybinrecs," Expected crbinrecs=",crbinrecs(ia)
  write(lout,"(a)")       "SIXTRACR> CRCHECK failure positioning binary files"
  call prror(-1)
105 continue
  write(lout,"(a)") ""
  write(lout,"(2(a,i0))") "SIXTRACR> ERROR PROBLEMS COPYING particle pair ",ia," IOSTAT=",ierro," from/to singletrackfile.dat"
  write(lout,"(3(a,i0))") "          mybinrecs=",mybinrecs," Expected crbinrecs=",crbinrecs(ia)," binrecs94=",binrecs94
  write(lout,"(a)")       "SIXTRACR> CRCHECK failure copying binary files"
  call prror(-1)
#endif
  ! We are not checkpointing or we have no checkpoints
  ! or we have no readable checkpoint
  ! If not checkpointing we can just give up on lout and use
  ! fort.6. We don't need to count records at all
605 continue
  write(93,"(a,l1)") "SIXTRACR> CRCHECK no restart possible checkp=",checkp
  flush(93)
  if (.not.checkp) then
    if (rerun) then
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
    rewind lout
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
1   endfile (6,iostat=ierro)
    backspace (6,iostat=ierro)
    ! This is not a FLUSH!
    rewind lout
    endfile (lout,iostat=ierro)
    close(lout)
    lout=6
  endif

  return

106 continue
  write(93,"(a,i0)")    "SIXTRACR> ERROR reading fort.6, iostat=",ierro
  write(93,"(2(a,i0))") "          sixrecs=",sixrecs," crsixrecs=",crsixrecs
  flush(93)
  write(lout,"(a)") "SIXTRACR> CRCHECK failure positioning fort.6"
  call prror(-1)

107 continue
  write(93,"(a,i0)") "SIXTRACR> ERROR reading fort.92, iostat=",ierro
  flush(93)
  write(lout,"(a)") "SIXTRACR> CRCHECK failure positioning fort.92"
  call prror(-1)

end subroutine crcheck

! ================================================================================================ !
!  This subroutine writes the checkpoint data to fort.95/96, and copies the new output from the
!  temporary (lout/fort.92) output file into fort.6.
!  The file fort.93 is used as a log file for the checkpoint/restarting.
! ================================================================================================ !
subroutine crpoint

  use floatPrecision
  use numerical_constants

  use dynk, only : dynk_enabled,dynk_getvalue,dynk_fSets_cr,dynk_cSets_unique,dynk_nSets_unique,dynk_filePos,dynk_crpoint
  use dump, only : dump_crpoint
  use scatter, only : scatter_active, scatter_crpoint

  use crcoall
  use parpro
  use mod_common
  use mod_commonmn
  use mod_commons
  use mod_commont
  use mod_commond
  use mod_hions
  use mod_version
  use mod_time
  use mod_meta
  use mod_units

  implicit none

  integer i,j,l,k,m,lstring,osixrecs,ncalls,maxncalls,crUnit
  logical lerror, fErr
#ifdef BOINC
  character(len=256) filename
#endif
  save

  ncalls = 0
#ifdef DEBUG
  maxncalls = 2000
#else
  maxncalls = 20
#endif

  if(ncalls <= maxncalls .or. numx >= nnuml-maxncalls) then
    write(93,"(2(a,i0))") "SIXTRACR> CRPOINT CALLED numlmax = ",numlmax,", numlcp = ",numlcp
    write(93,"(3(a,i0))") "SIXTRACR> CRPOINT CALLED lout = ",lout,", numx = ",numx,", numl = ",numl
    write(93,"(2(a,i0))") "SIXTRACR> CRPOINT CALLED binrec = ",binrec,", sixrec = ",sixrecs
    endfile(93,iostat=ierro)
    backspace(93,iostat=ierro)
  end if
  ncalls=ncalls+1
  if(restart) then
    restart=.false.
    return
  end if

  ! We need to copy fort.92 (lout) to fort.6 (sixrecs) (if it exists and we are not already using fort.6)
  osixrecs=sixrecs
  rewind lout
3 read(lout,'(a1024)',end=1,err=101,iostat=ierro) arecord
  lstring=1024
  do i=1024,2,-1
    lstring=i
    if (arecord(i:i) /= ' ') goto 2
    lstring=lstring-1
  end do
2 write(6,'(a)',err=102,iostat=ierro) arecord(1:lstring)
  sixrecs=sixrecs+1
  goto 3
1 if(sixrecs /= osixrecs) then
    endfile (6,iostat=ierro)
    backspace (6,iostat=ierro)
    rewind lout
    endfile(lout,iostat=ierro)
    close(lout)
    call units_openUnit(unit=92,fileName="fort.92",formatted=.true.,mode="rw",err=fErr)
#ifndef DEBUG
    if(ncalls <= 5 .or. numx >= numl) then
#endif
      write(93,"(2(a,i0))") "SIXTRACR> CRPOINT copied lout = ",lout,", sixrecs = ",sixrecs
      endfile(93,iostat=ierro)
      backspace(93,iostat=ierro)
#ifndef DEBUG
    end if
#endif
  else
    rewind lout
  end if
  call time_timerCheck(time3)
  ! Hope this is correct
  ! Maybe not!!!! this should be accumulative over multiple C/Rs
  time3=(time3-time1)+crtime3
  crnumlcr=numx+1

  if(dynk_enabled) then ! Store current settings of elements affected by DYNK
    if(ncalls <= maxncalls .or. numx >= nnuml-maxncalls) then
      write(93,"(a)") "SIXTRACR> CRPOINT filling dynk_fSets_cr"
      endfile(93,iostat=ierro)
      backspace(93,iostat=ierro)
    end if
    do j=1,dynk_nSets_unique
      dynk_fSets_cr(j) = dynk_getvalue(dynk_cSets_unique(j,1),dynk_cSets_unique(j,2))
    end do
  end if

  ! ********************
  !  Write the CR files
  ! ********************

  do crUnit=95,96

    lerror = .false.
    if(ncalls <= maxncalls .or. numx >= nnuml-maxncalls) then
      write(93,"(a,i0)") "SIXTRACR> CRPOINT writing fort.",crUnit
      endfile(93,iostat=ierro)
      backspace(93,iostat=ierro)
    end if
    rewind crUnit

    write(crUnit,err=100,iostat=ierro) version, moddate
    write(crUnit,err=100,iostat=ierro) &
      crnumlcr,                     &
      numl,                         &
      sixrecs,                      &
      binrec,                       &
      bnlrec,                       &
      bllrec,                       &
      sythckcr,                     &
      il,                           &
      time3,                        &
      napxo,                        &
      napx,                         &
      e0
    write(crUnit,err=100,iostat=ierro) &
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
      (xvl(1,j),j=1,napxo),         &
      (xvl(2,j),j=1,napxo),         &
      (yvl(1,j),j=1,napxo),         &
      (yvl(2,j),j=1,napxo),         &
      (dpsvl(j),j=1,napxo),         &
      (ejvl(j),j=1,napxo),          &
      (sigmvl(j),j=1,napxo)
    endfile(crUnit,iostat=ierro)
    backspace(crUnit,iostat=ierro)

    if(ncalls <= maxncalls .or. numx >= nnuml-maxncalls) then
      write(93,"(a,i0)") "SIXTRACR> CRPOINT Writing META variables to fort.",crUnit
      endfile(93,iostat=ierro)
      backspace(93,iostat=ierro)
    end if
    call meta_crpoint(crUnit,lerror,ierro)
    if(lerror) goto 100

    if(ncalls <= maxncalls .or. numx >= nnuml-maxncalls) then
      write(93,"(a,i0)") "SIXTRACR> CRPOINT Writing DUMP variables to fort.",crUnit
      endfile(93,iostat=ierro)
      backspace(93,iostat=ierro)
    end if
    call dump_crpoint(crUnit, lerror,ierro)
    if(lerror) goto 100

    if(ncalls <= maxncalls .or. numx >= nnuml-maxncalls) then
      write(93,"(a,i0)") "SIXTRACR> CRPOINT Writing HION variables to fort.",crUnit
      endfile(93,iostat=ierro)
      backspace(93,iostat=ierro)
    end if
    call hions_crpoint(crUnit,lerror,ierro)
    if(lerror) goto 100

    if(dynk_enabled) then
      if(ncalls <= maxncalls .or. numx >= nnuml-maxncalls) then
        write(93,"(a,i0)") "SIXTRACR> CRPOINT Writing DYNK variables to fort.",crUnit
        endfile(93,iostat=ierro)
        backspace(93,iostat=ierro)
      end if
      call dynk_crpoint(crUnit,lerror,ierro)
      if(lerror) goto 100
    end if

    if(scatter_active) then
      if(ncalls <= maxncalls .or. numx >= nnuml-maxncalls) then
        write(93,"(a,i0)") "SIXTRACR> CRPOINT Writing SCATTER variabless to fort.",crUnit
        endfile(93,iostat=ierro)
        backspace(93,iostat=ierro)
      end if
      call scatter_crpoint(crUnit,lerror,ierro)
      if(lerror) goto 100
    end if

    if(sythckcr) then
      if(ncalls <= maxncalls .or. numx >= nnuml-maxncalls) then
        write(93,"(a,i0)") "SIXTRACR> CRPOINT Writing EXTENDED varibless to fort.",crUnit
        endfile(93,iostat=ierro)
        backspace(93,iostat=ierro)
      end if
      if(ithick == 1) then
        if(ncalls <= maxncalls .or. numx >= nnuml-maxncalls) then
          write(93,"(a,i0)") "SIXTRACR> CRPOINT Writing EXTENDED variabless for THICK to fort.",crUnit
          endfile(93,iostat=ierro)
          backspace(93,iostat=ierro)
        end if
        write(crUnit,err=100,iostat=ierro) &
          ((((al(k,m,j,l),l=1,il),j=1,napxo),m=1,2),k=1,6), &
          ((((as(k,m,j,l),l=1,il),j=1,napxo),m=1,2),k=1,6)
        endfile(crUnit,iostat=ierro)
        backspace(crUnit,iostat=ierro)
      end if

      write(crUnit,err=100,iostat=ierro) &
        (aek(j),j=1,napxo),          &
        (afok(j),j=1,napxo),         &
        (as3(j),j=1,napxo),          &
        (as4(j),j=1,napxo),          &
        (as6(j),j=1,napxo),          &
        (co(j),j=1,napxo),           &
        (dpd(j),j=1,napxo),          &
        (dpsq(j),j=1,napxo),         &
        (fi(j),j=1,napxo),           &
        (fok(j),j=1,napxo),          &
        (fok1(j),j=1,napxo),         &
        (fokqv(j),j=1,napxo),        &
        (g(j),j=1,napxo),            &
        (gl(j),j=1,napxo),           &
        (hc(j),j=1,napxo),           &
        (hi(j),j=1,napxo),           &
        (hi1(j),j=1,napxo),          &
        (hm(j),j=1,napxo),           &
        (hp(j),j=1,napxo),           &
        (hs(j),j=1,napxo),           &
        (rho(j),j=1,napxo),          &
        (rhoc(j),j=1,napxo),         &
        (rhoi(j),j=1,napxo),         &
        (si(j),j=1,napxo),           &
        (siq(j),j=1,napxo),          &
        (sm1(j),j=1,napxo),          &
        (sm12(j),j=1,napxo),         &
        (sm2(j),j=1,napxo),          &
        (sm23(j),j=1,napxo),         &
        (sm3(j),j=1,napxo),          &
        (wf(j),j=1,napxo),           &
        (wfa(j),j=1,napxo),          &
        (wfhi(j),j=1,napxo)

      endfile(crUnit,iostat=ierro)
      backspace(crUnit,iostat=ierro)
    end if

  end do ! Loop over crUnit

104 continue
  return

100 continue
  write(93,"(a,i0)") "SIXTRACR> CRPOINT ERROR writing checkpt file, iostat = ",ierro
  goto 103

101 continue
  write(93,"(a,i0)") "SIXTRACR> CRPOINT ERROR reading lout fort.92, iostat = ",ierro
  goto 103

102 continue
  write(93,"(a,i0)") "SIXTRACR> CRPOINT ERROR writing fort.6, iostat = ",ierro

103 continue
  endfile(93,iostat=ierro)
  backspace (93,iostat=ierro)
  write(lout,"(a)") "SIXTRACR> CHECKPOINT I/O Error"
  call prror

end subroutine crpoint

! ================================================================================================ !
!  If we are restarting (restart is TRUE), this routine is called in the beginning of the tracking
!  loops. It is used to copy the cr* variables to the normal variables, e.g. crnapx -> napx etc.
!  The file fort.93 is used as a log file for the checkpoint/restarting.
! ================================================================================================ !
subroutine crstart

  use floatPrecision
  use numerical_constants
  use dynk, only : dynk_enabled, dynk_crstart
  use scatter, only: scatter_active, scatter_crstart

  use crcoall
  use parpro
  use mod_common
  use mod_commonmn
  use mod_commons
  use mod_commont
  use mod_commond
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
  endfile(93,iostat=ierro)
  backspace(93,iostat=ierro)
  numlcr=crnumlcr

  ! We do NOT reset numl so that a run can be extended for
  ! for more turns from the last checkpoint
  ! but we need to worry about numxv, nnumxv
  binrec=crbinrec
  bnlrec=crbnlrec
  bllrec=crbllrec
  sythckcr=crsythck

  ! the crtime3 is required (crtime0/1 removed)
  napxo=crnapxo
  napx=crnapx
  e0=cre0
  e0f=sqrt(e0**2-nucm0**2)

  write(93,"(a)") "SIXTRACR> CRSTART doing binrecs"
  endfile(93,iostat=ierro)
  backspace(93,iostat=ierro)

  do j=1,(napxo+1)/2
    binrecs(j)=crbinrecs(j)
  end do

  write(93,"(a)") "SIXTRACR> CRSTART doing normal NPART vars"
  endfile(93,iostat=ierro)
  backspace (93,iostat=ierro)
  do j=1,napxo
    numxv(j)=crnumxv(j)
    nnumxv(j)=crnnumxv(j)
    partID(j)=crpartID(j)
    parentID(j)=crparentID(j)
    pstop(j)=crpstop(j)
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
    xvl(1,j)=crxvl(1,j)
    xvl(2,j)=crxvl(2,j)
    yvl(1,j)=cryvl(1,j)
    yvl(2,j)=cryvl(2,j)
    dpsvl(j)=crdpsvl(j)
    ejvl(j)=crejvl(j)
    sigmvl(j)=crsigmvl(j)
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

  if (crsythck) then
    !ERICVARS now read the extended vars from fort.95/96.
    if(cril /= il) then
      write(lout,"(2(a,i0))") "SIXTRACR> CRSTART Problem as cril/il are different cril = ",cril,", il = ",il
      write(93,  "(2(a,i0))") "SIXTRACR> CRSTART Problem as cril/il are different cril = ",cril,", il = ",il
      endfile(93,iostat=ierro)
      backspace(93,iostat=ierro)
      write(lout,"(a)") "SIXTRACR> CRSTART Problem wih cril/il extended C/R"
      call prror
    end if

    !ERICVARS now read the extended vars from fort.95/96.
#ifdef DEBUG
  ! Commented out code for multiple records
  ! write(93,"(a)") "SIXTRACR> CRSTART DEBUG DUMP"
  ! call dump('Before xcrstart',0,0)
  ! endfile (93,iostat=ierro)
  ! backspace (93,iostat=ierro)
  ! write(93,"(a)") "SIXTRACR> CRSTART reading EXTENDED vars"
  ! endfile (93,iostat=ierro)
  ! backspace (93,iostat=ierro)
  ! if(read95) then
  !   read(95,end=100,err=100,iostat=ierro) ((((al(k,m,j,l),l=1,il),j=1,napxo),m=1,2),k=1,6)
  !   read(95,end=100,err=100,iostat=ierro) ((((as(k,m,j,l),l=1,il),j=1,napxo),m=1,2),k=1,6)
  !   read(95,end=100,err=100,iostat=ierro) (aek(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (afok(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (as3(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (as4(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (as6(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (co(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (dpd(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (dpsq(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (fi(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (fok(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (fok1(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (fokqv(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (g(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (gl(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (hc(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (hi(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (hi1(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (hm(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (hp(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (hs(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (rho(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (rhoc(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (rhoi(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (si(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (siq(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (sm1(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (sm12(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (sm2(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (sm23(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (sm3(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (wf(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (wfa(j),j=1,napxo)
  !   read(95,end=100,err=100,iostat=ierro) (wfhi(j),j=1,napxo)
  !   go to 102
  ! endif
#endif
    if(read95) then
      if(ithick == 1) then
        read(95,end=100,err=100,iostat=ierro) &
          ((((al(k,m,j,l),l=1,il),j=1,napxo),m=1,2),k=1,6),&
          ((((as(k,m,j,l),l=1,il),j=1,napxo),m=1,2),k=1,6)
      end if

      read(95,end=100,err=100,iostat=ierro) &
        (aek(j),j=1,napxo),   &
        (afok(j),j=1,napxo),  &
        (as3(j),j=1,napxo),   &
        (as4(j),j=1,napxo),   &
        (as6(j),j=1,napxo),   &
        (co(j),j=1,napxo),    &
        (dpd(j),j=1,napxo),   &
        (dpsq(j),j=1,napxo),  &
        (fi(j),j=1,napxo),    &
        (fok(j),j=1,napxo),   &
        (fok1(j),j=1,napxo),  &
        (fokqv(j),j=1,napxo), &
        (g(j),j=1,napxo),     &
        (gl(j),j=1,napxo),    &
        (hc(j),j=1,napxo),    &
        (hi(j),j=1,napxo),    &
        (hi1(j),j=1,napxo),   &
        (hm(j),j=1,napxo),    &
        (hp(j),j=1,napxo),    &
        (hs(j),j=1,napxo),    &
        (rho(j),j=1,napxo),   &
        (rhoc(j),j=1,napxo),  &
        (rhoi(j),j=1,napxo),  &
        (si(j),j=1,napxo),    &
        (siq(j),j=1,napxo),   &
        (sm1(j),j=1,napxo),   &
        (sm12(j),j=1,napxo),  &
        (sm2(j),j=1,napxo),   &
        (sm23(j),j=1,napxo),  &
        (sm3(j),j=1,napxo),   &
        (wf(j),j=1,napxo),    &
        (wfa(j),j=1,napxo),   &
        (wfhi(j),j=1,napxo)
      write(93,"(a)") "SIXTRACR> CRSTART read fort.95 EXTENDED OK"
      endfile(93,iostat=ierro)
      backspace(93,iostat=ierro)
      goto 102
    end if
    if(read96) then
      if(ithick == 1) then
        read(96,end=101,err=101,iostat=ierro) &
          ((((al(k,m,j,l),l=1,il),j=1,napxo),m=1,2),k=1,6), &
          ((((as(k,m,j,l),l=1,il),j=1,napxo),m=1,2),k=1,6)
      end if
      read(96,end=101,err=101,iostat=ierro) &
        (aek(j),j=1,napxo),   &
        (afok(j),j=1,napxo),  &
        (as3(j),j=1,napxo),   &
        (as4(j),j=1,napxo),   &
        (as6(j),j=1,napxo),   &
        (co(j),j=1,napxo),    &
        (dpd(j),j=1,napxo),   &
        (dpsq(j),j=1,napxo),  &
        (fi(j),j=1,napxo),    &
        (fok(j),j=1,napxo),   &
        (fok1(j),j=1,napxo),  &
        (fokqv(j),j=1,napxo), &
        (g(j),j=1,napxo),     &
        (gl(j),j=1,napxo),    &
        (hc(j),j=1,napxo),    &
        (hi(j),j=1,napxo),    &
        (hi1(j),j=1,napxo),   &
        (hm(j),j=1,napxo),    &
        (hp(j),j=1,napxo),    &
        (hs(j),j=1,napxo),    &
        (rho(j),j=1,napxo),   &
        (rhoc(j),j=1,napxo),  &
        (rhoi(j),j=1,napxo),  &
        (si(j),j=1,napxo),    &
        (siq(j),j=1,napxo),   &
        (sm1(j),j=1,napxo),   &
        (sm12(j),j=1,napxo),  &
        (sm2(j),j=1,napxo),   &
        (sm23(j),j=1,napxo),  &
        (sm3(j),j=1,napxo),   &
        (wf(j),j=1,napxo),    &
        (wfa(j),j=1,napxo),   &
        (wfhi(j),j=1,napxo)

      write(93,"(a)") "SIXTRACR> CRSTART read fort.96 EXTENDED OK"
      endfile(93,iostat=ierro)
      backspace(93,iostat=ierro)
      goto 102
    end if
100 continue
    write(93,"(a,i0)") "SIXTRACR> CRSTART Could not read checkpoint file 95 (extended), iostat = ",ierro
    goto 103
101 continue
    write(93,"(a,i0)") "SIXTRACR> CRSTART Could not read checkpoint file 96 (extended), iostat = ",ierro
103 continue
    endfile(93,iostat=ierro)
    backspace(93,iostat=ierro)
    write(lout,"(a)") "SIXTRACR> CRSTART Problem with extended checkpoint"
    call prror
  end if

102 continue
  write(93,"(3(a,i0))") "SIXTRACR> CRSTART sixrecs = ",sixrecs,", crsixrecs = ",crsixrecs,", binrec = ",binrec
  endfile(93,iostat=ierro)
  backspace(93,iostat=ierro)

  ! Just throw away our fort.92 stuff.
  rewind lout
  endfile(lout,iostat=ierro)
  close(lout)

  call units_openUnit(unit=lout,fileName="fort.92",formatted=.true.,mode="rw",err=fErr)
  ! but also add the rerun message
  write(lout,"(a80)") runtim
  runtim(1:20)="SIXTRACR restarted: "
  write(lout,"(a80)") runtim
  endfile(lout,iostat=ierro)
  backspace(lout,iostat=ierro)

  return

606 continue
  backspace(6,iostat=ierro)
  write(lout,"(2(a,i0))") "SIXTRACR> CRSTART Problem re-positioning fort.6: sixrecs = ",sixrecs,", crsixrecs = ",crsixrecs
  call prror

end subroutine crstart

end module checkpoint_restart
