! ================================================================================================ !
! A. Mereghetti, D. Sinuela Pastor and P. Garcia Ortega, for the FLUKA Team
! K. Sjobak, BE-ABP/HSS, BE-ABP/LAT
! V.K. Berglyd Olsen, BE-ABP-HSS
!
! In case the DUMP input block is issued, the beam population is dumped
!   at EACH occurence of the flagged SINGLE ELEMENT(s) in the accelerator structure
! Important remarks:
!   - The chosen SINGLE ELEMENT(s) must be outside a BLOC, and BLOCs cannot be chosen
!   - The special name 'ALL' will trigger dump at all SINGLE ELEMENTs (settings of dump
!     are stored in index 0 of all the usual arrays);
!   - The beam population is ALWAYS dumped at the end of the entry, i.e. AFTER the proper
!     transformation map is applied, and after the aperture check, i.e. AFTER the lost
!     particles are filtered out
!   - A negative or null value of the dump frequency is interpreted as dump at every turn
!   - No check is performed on the logical units, i.e. if the ones selected by the user
!     are used in other places of the code...
!   - The dump format can be changed to the one required by the LHC aperture check
!     post-processing tools, activating the dumpfmt flag (0=off, by default);
! ================================================================================================ !

module dump

  use floatPrecision
  use numerical_constants, only : zero

  implicit none

  logical, private, save :: ldumphighprec = .false. ! High precision printout required at all flagged SINGLE ELEMENTs
  logical, public,  save :: ldumpfront    = .false. ! Dump at the beginning of each element, not at the end.

  logical,          allocatable, public,  save :: ldump(:)      ! flag the SINGLE ELEMENT for dumping
  integer,          allocatable, private, save :: ndumpt(:)     ! Dump every n turns at a flagged SINGLE ELEMENT (dump frequency)
  integer,          allocatable, public,  save :: dumpfirst(:)  ! First turn for DUMP to be active
  integer,          allocatable, public,  save :: dumplast(:)   ! Last turn for this DUMP to be active (-1=all)
  integer,          allocatable, public,  save :: dumpunit(:)   ! Fortran unit for dump at a flagged SINGLE ELEMENT
  integer,          allocatable, public,  save :: dumpfmt(:)    ! Flag the format of the dump
  character(len=:), allocatable, public,  save :: dump_fname(:) ! Filename to write the dump to

  ! tas matrix used for nomalisation of phase space in DUMP and FMA.
  type, public :: dump_normType
    real(kind=fPrec) :: tas(6,6)    = zero
    real(kind=fPrec) :: invtas(6,6) = zero
    real(kind=fPrec) :: orbit(6)    = zero
  end type dump_normType

  type(dump_normType), allocatable, public, save :: dump_normMat(:)
  integer,             allocatable, public, save :: dump_nMatMap(:)
  integer,                          public, save :: dump_nNormMat = 0

#ifdef HDF5
  ! Array to save hdf5 formats for each dump format
  integer, allocatable, private, save :: dump_hdf5Format(:)
  integer, allocatable, private, save :: dump_hdf5DataSet(:)
#endif

#ifdef CR
  ! For resetting file positions
  integer, allocatable, private, save :: dumpfilepos(:), dumpfilepos_cr(:)
#endif

  integer, save :: rootfmt = 1000

  logical, save :: request_root = .false. ! Do we want to dump to a root tree?
! ================================================================================================================================ !
!  THE SUBROUTINES
! ================================================================================================================================ !
contains

! ================================================================================================================================ !
subroutine dump_expand_arrays(nele_new, nblz_new)

  use parpro
  use mod_alloc
  use numerical_constants, only : zero

  integer, intent(in) :: nele_new
  integer, intent(in) :: nblz_new

  call alloc(ldump,                 nele_new, .false., "ldump",       -1)
  call alloc(ndumpt,                nele_new, 0,       "ndumpt",      -1)
  call alloc(dumpfirst,             nele_new, 0,       "dumpfirst",   -1)
  call alloc(dumplast,              nele_new, 0,       "dumplast",    -1)
  call alloc(dumpunit,              nele_new, 0,       "dumpunit",    -1)
  call alloc(dumpfmt,               nele_new, 0,       "dumpfmt",     -1)
  call alloc(dump_fname, mFileName, nele_new, " ",     "dump_fname",  -1)
  call alloc(dump_nMatMap,          nele_new, 0,       "dump_nMatMap",-1)

#ifdef CR
  call alloc(dumpfilepos,           nele_new,-1,       "dumpfilepos",   -1)
  call alloc(dumpfilepos_cr,        nele_new,-1,       "dumpfilepos_cr",-1)
#endif

#ifdef HDF5
  call alloc(dump_hdf5DataSet,      nele_new,0,        "dump_hdf5DataSet",-1)
  call alloc(dump_hdf5Format,       200,     0,        "dump_hdf5Format")
#endif

end subroutine dump_expand_arrays

! ================================================================================================ !
!  V.K. Berglyd Olsen
!  Created: 2019-08-05
!  Updated: 2019-08-05
!  Set the normalisation matrix, compute its inverse, and save it in a compact array of structs.
! ================================================================================================ !
subroutine dump_setTasMatrix(elemID, tasData, cloData)

  use crcoall
  use matrix_inv
  use mod_common, only : bez
  use mod_settings

  integer,          intent(in) :: elemID       ! Single element index
  real(kind=fPrec), intent(in) :: tasData(6,6) ! Tas matrix with consistent units (no mm scaling)
  real(kind=fPrec), intent(in) :: cloData(6)   ! Closed orbit

  type(dump_normType), allocatable :: tmpNorm(:)
  real(kind=fPrec) invData(6,6)
  integer storeID

  if(dump_nMatMap(elemID) == 0) then
    if(allocated(dump_normMat)) then
      allocate(tmpNorm(dump_nNormMat+1))
      tmpNorm(1:dump_nNormMat) = dump_normMat(1:dump_nNormMat)
      call move_alloc(tmpNorm, dump_normMat)
      dump_nNormMat = dump_nNormMat + 1
    else
      allocate(dump_normMat(1))
      dump_nNormMat = 1
    end if
    dump_nMatMap(elemID) = dump_nNormMat
    storeID = dump_nNormMat
    if(st_debug) then
      if(elemID > 0) then
        write(lout,"(a,i0,a)") "DUMP> Saving normalisation matrix for element ",elemID," '"//trim(bez(elemID))//"'"
      else if(elemID == -1) then
        write(lout,"(a)") "DUMP> Saving normalisation matrix for element -1 'StartDUMP'"
      end if
    end if
  else
    storeID = dump_nMatMap(elemID)
    if(st_debug) then
      if(elemID > 0) then
        write(lout,"(a,i0,a)") "DUMP> Updating normalisation matrix for element ",elemID," '"//trim(bez(elemID))//"'"
      else if(elemID == -1) then
        write(lout,"(a)") "DUMP> Updating normalisation matrix for element -1 'StartDUMP'"
      end if
    end if
  end if

  call invert_tas(invData, tasData)

  dump_normMat(storeID)%tas    = tasData
  dump_normMat(storeID)%invtas = invData
  dump_normMat(storeID)%orbit  = cloData

end subroutine dump_setTasMatrix

! ================================================================================================ !
!  V.K. Berglyd Olsen
!  Created: 2019-08-05
!  Updated: 2019-08-05
!  Return the normalisation matrix for a given element. Used by the FMA module.
! ================================================================================================ !
subroutine dump_getTasMatrix(elemID, invData, tasData, cloData)

  use numerical_constants

  integer,          intent(in)  :: elemID
  real(kind=fPrec), intent(out) :: invData(6,6)
  real(kind=fPrec), intent(out) :: tasData(6,6)
  real(kind=fPrec), intent(out) :: cloData(6)

  invData(:,:) = zero
  tasData(:,:) = zero
  cloData(:)   = zero

  if(dump_nMatMap(elemID) == 0) then
    return
  end if

  invData = dump_normMat(dump_nMatMap(elemID))%invtas
  tasData = dump_normMat(dump_nMatMap(elemID))%tas
  cloData = dump_normMat(dump_nMatMap(elemID))%orbit

end subroutine dump_getTasMatrix

subroutine dump_lines(n,i,ix)

  use mod_common_track

  integer, intent(in) :: n,i,ix
  real(kind=fPrec) tmpTas(6,6), tmpInv(6,6), tmpClo(6)

  call dump_getTasMatrix(ix, tmpInv, tmpTas, tmpClo)

  if(ldump(0)) then
    ! Dump at all SINGLE ELEMENTs
    if(ndumpt(0) == 1 .or. mod(n,ndumpt(0)) == 1) then
      if ((n >= dumpfirst(0)) .and. ((n <= dumplast(0)) .or. (dumplast(0) == -1))) then
        call dump_beam_population(n, i, ix, dumpunit(0), dumpfmt(0), ldumphighprec, tmpClo, tmpInv)
      end if
    end if
  end if
  if(ktrack(i) /= 1) then
    ! The next "if" is only safe for SINGLE ELEMENTS, not BLOC where ix<0.
    if(ldump(ix)) then
      ! Dump at this precise SINGLE ELEMENT
      if(ndumpt(ix) == 1 .or. mod(n,ndumpt(ix)) == 1) then
        if((n >= dumpfirst(ix)) .and. ((n <= dumplast(ix)) .or. (dumplast(ix) == -1))) then
          call dump_beam_population(n, i, ix, dumpunit(ix), dumpfmt(ix), ldumphighprec, tmpClo, tmpInv)
        end if
      end if
    end if
  end if

end subroutine dump_lines

! ================================================================================================================================ !
subroutine dump_linesFirst(n)

  integer, intent(in) :: n

  ! StartDUMP - dump on the first element
  if(ldump(-1)) then
    if(ndumpt(-1) == 1 .or. mod(n,ndumpt(-1)) == 1) then
      if((n >= dumpfirst(-1)) .and. ((n <= dumplast(-1)) .or. (dumplast(-1) == -1))) then
        ! dump_nMatMap(-1) is set explicitly in main_cr, so we do not check here
        call dump_beam_population(n, 0, 0, dumpunit(-1), dumpfmt(-1), ldumphighprec, &
          dump_normMat(dump_nMatMap(-1))%orbit, dump_normMat(dump_nMatMap(-1))%invtas)
      end if
    end if
  end if

end subroutine dump_linesFirst

! ================================================================================================================================ !
subroutine dump_parseInputLine(inLine,iErr)

  use crcoall
  use mod_common
  use mod_units
  use string_tools
#ifdef HDF5
  use hdf5_output
#endif

  implicit none

  character(len=*), intent(in)    :: inLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  character(len=mNameLen) elemName
  character(len=mFileName) fileName
  integer i1,i2,i3,i4,i5,kk,j
  integer nSplit
  logical spErr

  ! initialise reading variables, to avoid storing nonsense values
  elemName = " " ! Element Name
  fileName = " " ! File Name
  i1       = 0   ! frequency
  i2       = -1  ! unit
  i3       = 0   ! format
  i4       = 1   ! first turn
  i5       = -1  ! last turn

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lerr,"(a)") "DUMP> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if
  if(nSplit == 0) return

  if(lnSplit(1) == "HIGH") then
    ldumphighprec = .true.
    return
  else if(lnSplit(1) == "FRONT") then
    ldumpfront = .true.
    return
  end if

  if(nSplit < 4 .or. nSplit > 7 .or. nSplit == 6) then
    write(lerr,"(a,i0)") "DUMP> ERROR Expected 4 to 7 (but not 6) arguments, got ",nSplit
    write(lerr,"(a)")   ("DUMP>     * '"//trim(lnSplit(kk))//"' ",kk=1,nSplit)
    iErr = .true.
    return
  end if

  if(len_trim(lnSplit(1)) > mNameLen) then
    write(lerr,"(2(a,i0))") "DUMP> ERROR Element names can be up to ",mNameLen," characters, got ",len_trim(lnSplit(1))
    iErr = .true.
    return
  end if

  elemname = trim(lnSplit(1))
  call chr_cast(lnSplit(2),i1,spErr)
  call chr_cast(lnSplit(3),i2,spErr)
  call chr_cast(lnSplit(4),i3,spErr)
  if(nSplit >= 5) then
    fileName = trim(lnSplit(5))
    if(len_trim(lnSplit(5)) > mFileName) then
      write(lerr,"(2(a,i0))") "DUMP> ERROR File names can be up to ",mFileName," characters, got ",len_trim(lnSplit(5))
      iErr = .true.
      return
    end if
  end if
  if(nSplit == 7) then
    call chr_cast(lnSplit(6),i4,spErr)
    call chr_cast(lnSplit(7),i5,spErr)
  end if

  ! Check that first/last turn is sane
  if(i5 /= -1) then
    if(i5 < i4) then
      write(lerr,"(2(a,i0))") "DUMP> ERROR Expect last turn >= first turn, unless last turn = -1 (infinity), got ", i4,", ",i5
      iErr = .true.
      return
    end if
  end if
  if(i4 < 1) then
    write(lerr,"(a,i0)") "DUMP> ERROR Expect first turn >= 1, got", i4
    iErr = .true.
    return
  end if

  ! Find it in the list of SINGLE ELEMENTs:
  do j=1,il
    if(bez(j) == elemName) then
      if(ldump(j)) then ! Only enable once/element!
        write(lerr,"(a)") "DUMP> ERROR Element '"//trim(elemName)//"' was specified more than once"
        iErr = .true.
        return
      end if

      ! Element was found in SINGLE ELEMENTS list, now do some sanity checks
      if(trim(bez(j)) == "ALL") then
        write(lerr,"(a)") "DUMP> ERROR The element name 'ALL' cannot be used in the SINGLE ELEMENTS list "//&
          "when an 'ALL' special DUMP is active."
        iErr = .true.
        return
      else if(trim(bez(j)) == "StartDUMP") then
        write(lerr,"(a)") "DUMP> ERROR The element name 'StartDUMP' cannot be used in the SINGLE ELEMENTS "// &
          "list when an 'StartDUMP' special DUMP is active."
        iErr = .true.
        return
      end if
      goto 10 ! Element found, store the data
    end if
  end do

  if(elemName == "ALL") then
    j=0
    if(ldump(j)) then
      write(lerr,"(a)") "DUMP> ERROR Element 'ALL' was specified (at least) twice"
      iErr = .true.
      return
      end if
    goto 10 ! Element found, store the data
  end if
  if(elemName == "StartDUMP") then
    j=-1
    if(ldump(j)) then
      write(lerr,"(a)") "DUMP> ERROR Element 'StartDUMP' was specified (at least) twice"
      iErr = .true.
      return
    end if
    goto 10 ! Element found, store the data
  end if

  ! Search failed, fall-through to here:
  write(lerr,"(a)") "DUMP> ERROR Unidentified SINGLE ELEMENT '"//trim(elemName)//"'"
  iErr = .true.
  return

  ! Element found, store the data:
10 continue
  ldump(j)      = .true.
  ndumpt(j)     = i1
  dumpunit(j)   = i2
  dumpfmt(j)    = i3
  dump_fname(j) = fileName
  dumpfirst(j)  = i4
  dumplast(j)   = i5
  if(ndumpt(j) <= 0) ndumpt(j) = 1
#ifdef HDF5
  if(h5_useForDUMP .eqv. .false.) then
#endif
    if(dump_fname(j) == " ") then
      dump_fname(j) = "dump_"//trim(bez(j))
    end if
    call f_requestUnit(trim(dump_fname(j)),dumpunit(j))
#ifdef HDF5
  end if
#endif

end subroutine dump_parseInputLine

! ================================================================================================================================ !
subroutine dump_parseInputDone(iErr)

  use crcoall
  use mod_common
  use string_tools
#ifdef ROOT
  use root_output, only : root_flag
#endif

  implicit none

  logical, intent(inout) :: iErr

  ! Temp variables
  integer ii,jj,kk
  logical have_root

  ! HEADER
! 10460 format(//131('-')//t10,'DATA BLOCK ',a4,' INFOs'/ /t10, 'ELEMENT NAME',8x,'EVERY # TURNs',2x, &
! 'LOGICAL UNIT',2x,'FILENAME',24x,'FORMAT',5x, "FirstTurn",6x,"LastTurn")
! 10470 format(t10,a16,4x,i13,2x,i12,2x,a32,i6,2x,i12,2x,i12)
! 10472 format(t10,a)
  write(lout,"(a)") "DUMP> The last column states the format of the output file (see manual):"
  write(lout,"(a)") "DUMP> Element Name          Every #Turns  Unit#  Filename                        "//&
    "Format  FirstTurn  LastTurn"

  ! ldump(0)=.true. : DUMP all elements found
  if (ldump(0)) then
    write(lout,"(a,a20,2x,i12,2x,i5,2x,a32,i6,2x,i9,2x,i8)") "DUMP> ","All Single Elements ", &
      ndumpt(0),dumpunit(0),chr_rpad(trim(dump_fname(0)),32),dumpfmt(0),dumpfirst(0),dumplast(0)
  end if
  if (ldump(-1)) then
    write(lout,"(a,a20,2x,i12,2x,i5,2x,a32,i6,2x,i9,2x,i8)") "DUMP> ","StartDUMP specified ", &
      ndumpt(0),dumpunit(0),chr_rpad(trim(dump_fname(0)),32),dumpfmt(0),dumpfirst(0),dumplast(0)
  end if

  do ii=1,il
    if(ldump(ii)) then
      write(lout,"(a,a20,2x,i12,2x,i5,2x,a32,i6,2x,i9,2x,i8)") "DUMP> ", &
        bez(ii)(1:20),ndumpt(ii),dumpunit(ii),chr_rpad(trim(dump_fname(ii)),32),dumpfmt(ii),dumpfirst(ii),dumplast(ii)
      ! At which structure indices is this single element found? (Sanity check)
      kk = 0
      do jj=1,mper*mbloz      ! Loop over all structure elements
        if(ic(jj)-nblo == ii) then
          write (lout,"(a,i0)") "DUMP> Found as structure element no. ",jj
          kk = kk + 1
        end if
      end do
      if(kk == 0) then
        write (lerr,"(a)") "DUMP> ERROR No structure elements found for '"//trim(bez(ii))//"'"
        write (lerr,"(a)") "DUMP>       This element is probably only found in a BLOC, or it is not used at all."
        iErr = .true.
        return
      end if
      ! check if we are using a root dump
      if(dumpfmt(ii) .eq. rootfmt) then
        request_root = .true.
      end if
    end if
  end do

  if (ldumphighprec) then
    write(lout,"(a)") "DUMP> Requested high precision dumping"
  end if
  if (ldumpfront) then
    write(lout,"(a)") "DUMP> Requested FRONT dumping"
  end if
  return

! Is there root support compiled in?
  have_root = .false.
#ifdef ROOT
  have_root = .true.
#endif

! If we found a root dump, check we have root support built in.
  if(request_root .eqv. .true.) then
!   Check we have root compiled
    if(have_root .eqv. .true.) then
!     Check we have root enabled
#ifdef ROOT
      if(root_flag .eqv. .false.) then
!     explode and exit
        write(lerr,"(a)") "DUMP> Requested ROOT output in dump, but ROOT output is not enabled!"
        call prror
      end if
#endif
    else
!   No root support compiled in
!   explode and exit
      write(lerr,"(a)") "DUMP> Requested ROOT output in dump, but this executable has no ROOT support!"
      call prror
    end if
  end if
end subroutine dump_parseInputDone

! ================================================================================================================================ !
!  A.Mereghetti, D.Sinuela-Pastor & P.Garcia Ortega, for the FLUKA Team
!  K.Sjobak, A.Santamaria, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-07-13
! ================================================================================================================================ !
subroutine dump_initialise

  use numerical_constants, only : zero
  use crcoall
  use string_tools
  use mod_common
  use mod_units
#ifdef HDF5
  use hdf5_output
#endif
#ifdef ROOT
  use root_output
#endif

  integer i,j,k,l
  logical lOpen, rErr
  character(len=16) tasbuf(6,6)
  real(kind=fPrec) tmpTas(6,6), tmpInv(6,6), tmpClo(6)

#ifdef HDF5
  type(h5_dataField), allocatable :: setFields(:)
  if(h5_useForDUMP) then
    call h5_initForDump()
    goto 3100 ! Skip the normal file initialisation, and do the HDF5 version instead
  end if
#endif

  do i=-1,il
#ifdef CR
    if(dumpfilepos(i) >= 0) then
      ! Expect the file to be opened already, in crcheck
      inquire(unit=dumpunit(i), opened=lopen)
      if(.not.lopen) then
        write(lerr,"(2(a,i0),a)") "DUMP> ERROR The unit ",dumpunit(i)," has dumpfilepos = ", dumpfilepos(i), " >= 0, "//&
          "but the file is NOT open. This is probably a bug."
        call prror
      end if
      cycle ! Everything OK, don't try to open the files again.
    end if
#endif
    if(ldump(i)) then
      ! The same file could be used by more than one SINGLE ELEMENT
      inquire( unit=dumpunit(i), opened=lopen )
      if (.not.lopen) then
        ! Check that the filename is not already taken
        do j=-1,i-1
          if (ldump(j) .and. (dump_fname(j) == dump_fname(i))) then
            write(lerr,"(2(a,i0))") "DUMP> ERROR Output filename '"//trim(dump_fname(i))//&
              "' is used by two DUMPS, but output units differ: ",dumpunit(i)," vs ",dumpunit(j)
            call prror
          end if
        end do
        if (dumpfmt(i) == 3 .or. dumpfmt(i) == 8 .or. dumpfmt(i) == 101) then ! Binary dump
          call f_open(unit=dumpunit(i),file=trim(dump_fname(i)),formatted=.false.,mode="w",status="replace")
#ifdef ROOT
        else if (dumpfmt(i) == rootfmt) then
          !no need to open a file for root
          !We do need to call the init function though
          call root_BunchDumpInit()
#endif
        else ! ASCII dump
          call f_open(unit=dumpunit(i),file=trim(dump_fname(i)),formatted=.true.,mode="w",status="replace")
        end if
#ifdef CR
        dumpfilepos(i) = 0
#endif
      else ! lopen was .TRUE.

        ! Sanity check: If file number i is already open, it should be by another DUMP
        ! (but we can't guarantee that files opened later are handled correctly)
        ! Also, a file should not be shared with element 0 (all) or -1 (StartDUMP)
        ! All dumps writing to the same file (unit) must have the same format and filename.
        ! If everything is OK, add to the header.

        ! Reuse the lopen flag as a temp variable
        lopen = .false.
        do j=-1,i-1 ! Search all possible DUMPs up to but not including the one we're looking at (number i)
          if (ldump(j)) then
            if (dumpunit(j) == dumpunit(i)) then
              if (dumpfmt(j) /= dumpfmt(i)) then
                write(lerr,"(a,i0,a)") "DUMP> ERROR Output unit ",dumpunit(i)," used by two DUMPS, formats are not the same."
                call prror
              else if (j == 0) then
                write(lerr,"(a,i0,a)") "DUMP> ERROR Output unit ",dumpunit(i)," used by two DUMPS, one of which is ALL"
                call prror
              else if (j == -1) then
                write(lerr,"(a,i0,a)") "DUMP> ERROR Output unit ",dumpunit(i)," used by two DUMPS, one of which is StartDUMP"
                call prror
              else if (dump_fname(j) /= dump_fname(i)) then
                write(lerr,"(a,i0,a)") "DUMP> ERROR Output unit ",dumpunit(i)," used by two DUMPS, but filenames differ: '"//&
                  trim(dump_fname(i)),"' vs '",trim(dump_fname(j)),"'"
                call prror
              else
                ! Everything is fine
                lopen = .true.
#ifdef CR
                ! Dumpfilepos is separate for every element, even if they share files.
                dumpfilepos(i) = 0
#endif
              end if
            end if ! IF file unit matches
          end if ! IF ldump(j)
        end do ! DO loop over j

        ! LOPEN not set to true by sanity check in loop above
        ! => File was already open, but not by DUMP.
        if (.not.lopen) then
          write(lerr,"(a,i0,a)") "DUMP> ERROR Unit ",dumpunit(i)," is already open, but not by DUMP. Please pick another unit!"
          write(lerr,"(a)")      "DUMP> Note: This check is not watertight as other parts of the program may later open the "
          write(lerr,"(a)")      "DUMP>       same unit. Althernatively, the unit can be specified as -1 and a unit is assigned."
          call prror
        end if
      end if

      ! Write format-specific headers
      call dump_getTasMatrix(i, tmpInv, tmpTas, tmpClo)
      if(dumpfmt(i) == 1) then
        write(dumpunit(i),'(a)') '# particleID turn s[m] x[mm] xp[mrad] y[mm] yp[mrad] (E-E0)/E0[1] ktrack'
        flush(dumpunit(i))
#ifdef CR
        dumpfilepos(i) = dumpfilepos(i) + 1
#endif
      else if (dumpfmt(i) == 2 .or. dumpfmt(i) == 4 .or. dumpfmt(i) == 5 .or. &
               dumpfmt(i) == 6 .or. dumpfmt(i) == 7 .or. dumpfmt(i) == 9 ) then
        ! Write the general header
        if (i == -1) then  ! STARTdump
          write(dumpunit(i),'(a,i0,a,a16,4(a,i12),2(a,L1))') &
            '# DUMP format #',dumpfmt(i),', START=',bez(1), ', number of particles=',napx, ', dump period=',ndumpt(i), &
            ', first turn=', dumpfirst(i), ', last turn=',dumplast(i), ', HIGH=',ldumphighprec, ', FRONT=',ldumpfront
        else if (i == 0) then ! ALL
          write(dumpunit(i),'(a,i0,a,4(a,i12),2(a,L1))') &
            '# DUMP format #',dumpfmt(i),', ALL ELEMENTS,', ' number of particles=',napx, ', dump period=',ndumpt(i), &
            ', first turn=', dumpfirst(i), ', last turn=',dumplast(i), ', HIGH=',ldumphighprec, ', FRONT=',ldumpfront
        else ! Normal element
          write(dumpunit(i),'(a,i0,a,a16,4(a,i12),2(a,L1))') &
            '# DUMP format #',dumpfmt(i), ', bez=', bez(i), ', number of particles=',napx,', dump period=',ndumpt(i), &
            ', first turn=',dumpfirst(i), ', last turn=',dumplast(i), ', HIGH=',ldumphighprec, ', FRONT=',ldumpfront
        end if

        ! Write the format-specific headers:
        if (dumpfmt(i) == 2) then ! FORMAT 2
          write(dumpunit(i),'(a,a)') '# particleID turn s[m] x[mm] xp[mrad] y[mm] yp[mrad] sigma[mm] (E-E0)/E0[1] ktrack'
        else if (dumpfmt(i) == 200) then ! FORMAT 200
          write(dumpunit(i),'(a,a)') '# particleID turn s[m] x[mm] xp[mrad] y[mm] yp[mrad] sigma[mm] (E-E0)/E0[1] ktrack '// &
&                                    'Z A Q PDGid'
        else if (dumpfmt(i) == 4) then ! FORMAT 4
          write(dumpunit(i),'(a)') '# napx turn s[m] <x>[mm] <xp>[mrad] <y>[mm] <yp>[mrad] <sigma>[mm] <(E-E0)/E0>[1]'
        else if (dumpfmt(i) == 5) then ! FORMAT 5
          write(dumpunit(i),'(a)') '# napx turn s[m] ' //                     &
            '<x>[mm] <xp>[mrad] <y>[mm] <yp>[mrad] <sigma>[mm] <(E-E0)/E0>[1] '// &
            '<x^2> <x*xp> <x*y> <x*yp> <x*sigma> <x*(E-E0)/E0> '//                &
            '<xp^2> <xp*y> <xp*yp> <xp*sigma> <xp*(E-E0)/E0> '//                  &
            '<y^2> <y*yp> <y*sigma> <y*(E-E0)/E0> '//                             &
            '<yp^2> <yp*sigma> <yp*(E-E0)/E0> '//                                 &
            '<sigma^2> <sigma*(E-E0)/E0> '//                                          &
            '<((E-E0)/E0)^2>'
        else if (dumpfmt(i) == 6) then ! FORMAT 6
          write(dumpunit(i),'(a)') '# napx turn s[m] ' //                  &
            '<x>[m] <px>[1] <y>[m] <py>[1] <sigma>[m] <psigma>[1] '//      &
            '<x^2> <x*px> <x*y> <x*py> <x*sigma> <x*psigma> '//            &
            '<px^2> <px*y> <px*py> <px*sigma> <px*psigma> '//              &
            '<y^2> <y*py> <y*sigma> <y*psigma> '//                         &
            '<py^2> <py*sigma> <py*psigma> '//                             &
            '<sigma^2> <sigma*psigma> '//                                  &
            '<psigma^2>'
        else if (dumpfmt(i) == 7 .or. dumpfmt(i) == 9) then
          ! Normalized ASCII dump -> extra headers with matrices and closed orbit
          if (dumpfmt(i) == 7) then ! FORMAT 7
            write(dumpunit(i),'(a)') '# particleID turn s[m] nx[1.e-3 sqrt(m)] npx[1.e-3 sqrt(m)] '// &
              'ny[1.e-3 sqrt(m)] npy[1.e-3 sqrt(m)] nsigma[1.e-3 sqrt(m)] ndp/p[1.e-3 sqrt(m)] ktrack'
          end if
          if (dumpfmt(i) == 9) then ! FORMAT 9
            write(dumpunit(i),'(a)') '# napx turn s[m] ' //                   &
              '<nx>[1.e-3 sqrt(m)] <npx>[1.e-3 sqrt(m)] '//                   &
              '<ny>[1.e-3 sqrt(m)] <npy>[1.e-3 sqrt(m)] '//                   &
              '<nsig>[1.e-3 sqrt(m)] <ndp/p>[1.e-3 sqrt(m)]'//                & ! There should be a space at the end of the string
              '<nx^2> <nx*npx> <nx*ny> <nx*npy> <nx*nsigma> <nx*npsigma> '//  &
              '<npx^2> <npx*ny> <npx*npy> <npx*nsigma> <npx*npsigma> '//      &
              '<ny^2> <ny*npy> <ny*nsigma> <ny*npsigma> '//                   &
              '<npy^2> <npy*nsigma> <npy*npsigma> '//                         &
              '<nsigma^2> <nsigma*npsigma> '//                                &
              '<npsigma^2>'
          end if
          ! closed orbit
          ! units: x,xp,y,yp,sig,dp/p = [mm,mrad,mm,mrad,1]
          ! (note: units are already changed in linopt part)
          do k=1,6
            call chr_fromReal(tmpClo(k),tasbuf(k,1),10,2,rErr)
          end do
          write(dumpunit(i),"(a,1x,6(1x,a16))") "# closed orbit [mm,mrad,mm,mrad,1]", &
            tasbuf(1,1),tasbuf(2,1),tasbuf(3,1),tasbuf(4,1),tasbuf(5,1),tasbuf(6,1)

          do k=1,6
            do l=1,6
              call chr_fromReal(tmpTas(l,k),tasbuf(l,k),10,2,rErr)
            end do
          end do
          write(dumpunit(i),"(a,1x,36(1x,a16))") "# tamatrix [mm,mrad,mm,mrad,1]", &
            tasbuf(1,1),tasbuf(1,2),tasbuf(1,3),tasbuf(1,4),tasbuf(1,5),tasbuf(1,6), &
            tasbuf(2,1),tasbuf(2,2),tasbuf(2,3),tasbuf(2,4),tasbuf(2,5),tasbuf(2,6), &
            tasbuf(3,1),tasbuf(3,2),tasbuf(3,3),tasbuf(3,4),tasbuf(3,5),tasbuf(3,6), &
            tasbuf(4,1),tasbuf(4,2),tasbuf(4,3),tasbuf(4,4),tasbuf(4,5),tasbuf(4,6), &
            tasbuf(5,1),tasbuf(5,2),tasbuf(5,3),tasbuf(5,4),tasbuf(5,5),tasbuf(5,6), &
            tasbuf(6,1),tasbuf(6,2),tasbuf(6,3),tasbuf(6,4),tasbuf(6,5),tasbuf(6,6)

          do k=1,6
            do l=1,6
              call chr_fromReal(tmpInv(l,k),tasbuf(l,k),10,2,rErr)
            end do
          end do
          write(dumpunit(i),"(a,1x,36(1x,a16))") "# inv(tamatrix)", &
            tasbuf(1,1),tasbuf(1,2),tasbuf(1,3),tasbuf(1,4),tasbuf(1,5),tasbuf(1,6), &
            tasbuf(2,1),tasbuf(2,2),tasbuf(2,3),tasbuf(2,4),tasbuf(2,5),tasbuf(2,6), &
            tasbuf(3,1),tasbuf(3,2),tasbuf(3,3),tasbuf(3,4),tasbuf(3,5),tasbuf(3,6), &
            tasbuf(4,1),tasbuf(4,2),tasbuf(4,3),tasbuf(4,4),tasbuf(4,5),tasbuf(4,6), &
            tasbuf(5,1),tasbuf(5,2),tasbuf(5,3),tasbuf(5,4),tasbuf(5,5),tasbuf(5,6), &
            tasbuf(6,1),tasbuf(6,2),tasbuf(6,3),tasbuf(6,4),tasbuf(6,5),tasbuf(6,6)

        end if ! Format-specific headers
        flush(dumpunit(i))
#ifdef CR
        dumpfilepos(i) = dumpfilepos(i) + 2
        ! format 7 also writes clo, tas and tasinv
        if (dumpfmt(i) == 7 .or. dumpfmt(i) == 9) then
          dumpfilepos(i) = dumpfilepos(i) + 3
        end if
#endif
      end if ! END if format 2/4/5/6/7/9 -> General header

      ! Normalized DUMP
      if (dumpfmt(i) == 7 .or. dumpfmt(i) == 8 .or. dumpfmt(i) == 9) then
        ! Have a matrix that's not zero (i.e. did we put a 6d LINE block?)
        if(tmpTas(1,1) == zero .and. tmpTas(1,2) == zero .and. tmpTas(1,3) == zero .and. tmpTas(1,4) == zero) then
          write(lerr,"(a)") "DUMP> ERROR The normalization matrix appears to not be set. Did you forget to put a 6D LINE block?"
          call prror
        end if
        if(idp == 0 .or. ition == 0) then ! We're in the 4D case
          if(i /= -1) then ! Not at StartDUMP
            write(lerr,"(a)") "DUMP> ERROR in normalized DUMP: 4D only supported for StartDUMP!"
            call prror
          end if
        end if
      end if ! END if normalized dump
    end if ! If ldump(i) -> Dump on this element
  end do ! Loop over elements with index i

  return

  ! The following section only applies if HDF5 is available AND enabled
#ifdef HDF5
3100 continue
  do i=-1,il
    if(ldump(i)) then
      select case(dumpfmt(i))

      case(1)
        ! Format 1:
        ! # particleID turn s[m] x[mm] xp[mrad] y[mm] yp[mrad] (E-E0)/E0[1] ktrack
        if(dump_hdf5Format(1) == 0) then
          allocate(setFields(9))
          setFields(1)  = h5_dataField(name="ID",     type=h5_typeInt)
          setFields(2)  = h5_dataField(name="TURN",   type=h5_typeInt)
          setFields(3)  = h5_dataField(name="S",      type=h5_typeReal)
          setFields(4)  = h5_dataField(name="X",      type=h5_typeReal)
          setFields(5)  = h5_dataField(name="XP",     type=h5_typeReal)
          setFields(6)  = h5_dataField(name="Y",      type=h5_typeReal)
          setFields(7)  = h5_dataField(name="YP",     type=h5_typeReal)
          setFields(8)  = h5_dataField(name="dE/E",   type=h5_typeReal)
          setFields(9)  = h5_dataField(name="KTRACK", type=h5_typeInt)
          call h5_createFormat("dumpFormat1", setFields, dump_hdf5Format(1))
        end if
        call h5_createDataSet(dump_fname(i), h5_dumpID, dump_hdf5Format(1), dump_hdf5DataSet(i), napx)
        block
          character(len=:), allocatable :: colNames(:)
          character(len=:), allocatable :: colUnits(:)
          logical spErr
          integer nSplit
          call chr_split("ID turn s x xp y yp (E-E0)/E0 ktrack",colNames,nSplit,spErr)
          call chr_split("1 1 m mm mrad mm mrad 1 1",colUnits,nSplit,spErr)
          call h5_writeDataSetAttr(dump_hdf5DataSet(i),"dumpFormat",1)
          call h5_writeDataSetAttr(dump_hdf5DataSet(i),"colNames",  colNames)
          call h5_writeDataSetAttr(dump_hdf5DataSet(i),"colUnits",  colUnits)
        end block

      case(2)
        ! Format 2:
        ! # particleID turn s[m] x[mm] xp[mrad] y[mm] yp[mrad] sigma[mm] (E-E0)/E0[1] ktrack
        if(dump_hdf5Format(2) == 0) then
          allocate(setFields(10))
          setFields(1)  = h5_dataField(name="ID",     type=h5_typeInt)
          setFields(2)  = h5_dataField(name="TURN",   type=h5_typeInt)
          setFields(3)  = h5_dataField(name="S",      type=h5_typeReal)
          setFields(4)  = h5_dataField(name="X",      type=h5_typeReal)
          setFields(5)  = h5_dataField(name="XP",     type=h5_typeReal)
          setFields(6)  = h5_dataField(name="Y",      type=h5_typeReal)
          setFields(7)  = h5_dataField(name="YP",     type=h5_typeReal)
          setFields(8)  = h5_dataField(name="SIGMA",  type=h5_typeReal)
          setFields(9)  = h5_dataField(name="dE/E",   type=h5_typeReal)
          setFields(10) = h5_dataField(name="KTRACK", type=h5_typeInt)
          call h5_createFormat("dumpFormat2", setFields, dump_hdf5Format(2))
        end if
        call h5_createDataSet(dump_fname(i), h5_dumpID, dump_hdf5Format(2), dump_hdf5DataSet(i), napx)
        block
          character(len=:), allocatable :: colNames(:)
          character(len=:), allocatable :: colUnits(:)
          logical spErr
          integer nSplit
          call chr_split("ID turn s x xp y yp sigma (E-E0)/E0 ktrack",colNames,nSplit,spErr)
          call chr_split("1 1 m mm mrad mm mrad mm 1 1",colUnits,nSplit,spErr)
          call h5_writeDataSetAttr(dump_hdf5DataSet(i),"dumpFormat",2)
          call h5_writeDataSetAttr(dump_hdf5DataSet(i),"colNames",  colNames)
          call h5_writeDataSetAttr(dump_hdf5DataSet(i),"colUnits",  colUnits)
        end block

      case(3)
        ! Format 3:
        ! # particleID turn s[m] x[mm] xp[mrad] y[mm] yp[mrad] sigma[mm] (E-E0)/E0[1] ktrack
        if(dump_hdf5Format(3) == 0) then
          allocate(setFields(10))
          setFields(1)  = h5_dataField(name="ID",     type=h5_typeInt)
          setFields(2)  = h5_dataField(name="TURN",   type=h5_typeInt)
          setFields(3)  = h5_dataField(name="S",      type=h5_typeReal)
          setFields(4)  = h5_dataField(name="X",      type=h5_typeReal)
          setFields(5)  = h5_dataField(name="XP",     type=h5_typeReal)
          setFields(6)  = h5_dataField(name="Y",      type=h5_typeReal)
          setFields(7)  = h5_dataField(name="YP",     type=h5_typeReal)
          setFields(8)  = h5_dataField(name="SIGMA",  type=h5_typeReal)
          setFields(9)  = h5_dataField(name="dE/E",   type=h5_typeReal)
          setFields(10) = h5_dataField(name="KTRACK", type=h5_typeInt)
          call h5_createFormat("dumpFormat3", setFields, dump_hdf5Format(3))
        end if
        call h5_createDataSet(dump_fname(i), h5_dumpID, dump_hdf5Format(3), dump_hdf5DataSet(i), napx)
        block
          character(len=:), allocatable :: colNames(:)
          character(len=:), allocatable :: colUnits(:)
          logical spErr
          integer nSplit
          call chr_split("ID turn s x xp y yp sigma (E-E0)/E0 ktrack",colNames,nSplit,spErr)
          call chr_split("1 1 m mm mrad mm mrad mm 1 1",colUnits,nSplit,spErr)
          call h5_writeDataSetAttr(dump_hdf5DataSet(i),"dumpFormat",3)
          call h5_writeDataSetAttr(dump_hdf5DataSet(i),"colNames",  colNames)
          call h5_writeDataSetAttr(dump_hdf5DataSet(i),"colUnits",  colUnits)
        end block

      case(4)
        ! Format 4:
        ! # napx turn s[m] <x>[mm] <xp>[mrad] <y>[mm] <yp>[mrad] <sigma>[mm] <(E-E0)/E0>[1]
        if(dump_hdf5Format(4) == 0) then
          allocate(setFields(9))
          setFields(1)  = h5_dataField(name="NAPX",   type=h5_typeInt)
          setFields(2)  = h5_dataField(name="TURN",   type=h5_typeInt)
          setFields(3)  = h5_dataField(name="S",      type=h5_typeReal)
          setFields(4)  = h5_dataField(name="X",      type=h5_typeReal)
          setFields(5)  = h5_dataField(name="XP",     type=h5_typeReal)
          setFields(6)  = h5_dataField(name="Y",      type=h5_typeReal)
          setFields(7)  = h5_dataField(name="YP",     type=h5_typeReal)
          setFields(9)  = h5_dataField(name="Z",      type=h5_typeReal)
          setFields(8)  = h5_dataField(name="dE/E",   type=h5_typeReal)
          call h5_createFormat("dumpFormat4", setFields, dump_hdf5Format(4))
        end if
        call h5_createDataSet(dump_fname(i), h5_dumpID, dump_hdf5Format(4), dump_hdf5DataSet(i), numl)
        block
          character(len=:), allocatable :: colNames(:)
          character(len=:), allocatable :: colUnits(:)
          logical spErr
          integer nSplit
          call chr_split("napx turn s <x> <xp> <y> <yp> <sigma> <(E-E0)/E0>",colNames,nSplit,spErr)
          call chr_split("1 1 m mm mrad mm mrad mm 1",colUnits,nSplit,spErr)
          call h5_writeDataSetAttr(dump_hdf5DataSet(i),"dumpFormat",4)
          call h5_writeDataSetAttr(dump_hdf5DataSet(i),"colNames",  colNames)
          call h5_writeDataSetAttr(dump_hdf5DataSet(i),"colUnits",  colUnits)
        end block

      case(5)
        ! Format 5:
        ! # napx turn s[m] <x>[mm] <xp>[mrad] <y>[mm] <yp>[mrad] <sigma>[mm] <(E-E0)/E0>[1]
        ! <x^2> <x*xp> <x*y> <x*yp> <x*sigma> <x*(E-E0)/E0>
        ! <xp^2> <xp*y> <xp*yp> <xp*sigma> <xp*(E-E0)/E0>
        ! <y^2> <y*yp> <y*sigma> <y*(E-E0)/E0>
        ! <yp^2> <yp*sigma> <yp*(E-E0)/E0>
        ! <sigma^2> <sigma*(E-E0)/E0>
        ! <((E-E0)/E0)^2>
        if(dump_hdf5Format(5) == 0) then
          allocate(setFields(30))
          setFields(1)  = h5_dataField(name="NAPX",      type=h5_typeInt)
          setFields(2)  = h5_dataField(name="TURN",      type=h5_typeInt)
          setFields(3)  = h5_dataField(name="S",         type=h5_typeReal)
          setFields(4)  = h5_dataField(name="X",         type=h5_typeReal)
          setFields(5)  = h5_dataField(name="XP",        type=h5_typeReal)
          setFields(6)  = h5_dataField(name="Y",         type=h5_typeReal)
          setFields(7)  = h5_dataField(name="YP",        type=h5_typeReal)
          setFields(8)  = h5_dataField(name="Z",         type=h5_typeReal)
          setFields(9)  = h5_dataField(name="dE/E",      type=h5_typeReal)
          setFields(10) = h5_dataField(name="X^2",       type=h5_typeReal)
          setFields(11) = h5_dataField(name="X*XP",      type=h5_typeReal)
          setFields(12) = h5_dataField(name="X*Y",       type=h5_typeReal)
          setFields(13) = h5_dataField(name="X*YP",      type=h5_typeReal)
          setFields(14) = h5_dataField(name="X*Z",       type=h5_typeReal)
          setFields(15) = h5_dataField(name="X*(dE/E)",  type=h5_typeReal)
          setFields(16) = h5_dataField(name="XP^2",      type=h5_typeReal)
          setFields(17) = h5_dataField(name="XP*Y",      type=h5_typeReal)
          setFields(18) = h5_dataField(name="XP*YP",     type=h5_typeReal)
          setFields(19) = h5_dataField(name="XP*Z",      type=h5_typeReal)
          setFields(20) = h5_dataField(name="XP*(dE/E)", type=h5_typeReal)
          setFields(21) = h5_dataField(name="Y^2",       type=h5_typeReal)
          setFields(22) = h5_dataField(name="Y*YP",      type=h5_typeReal)
          setFields(23) = h5_dataField(name="Y*Z",       type=h5_typeReal)
          setFields(24) = h5_dataField(name="Y*(dE/E)",  type=h5_typeReal)
          setFields(25) = h5_dataField(name="YP^2",      type=h5_typeReal)
          setFields(26) = h5_dataField(name="YP*Z",      type=h5_typeReal)
          setFields(27) = h5_dataField(name="YP*(dE/E)", type=h5_typeReal)
          setFields(28) = h5_dataField(name="Z^2",       type=h5_typeReal)
          setFields(29) = h5_dataField(name="Z*(dE/E)",  type=h5_typeReal)
          setFields(30) = h5_dataField(name="(dE/E)^2",  type=h5_typeReal)
          call h5_createFormat("dumpFormat5", setFields, dump_hdf5Format(5))
        end if
        call h5_createDataSet(dump_fname(i), h5_dumpID, dump_hdf5Format(5), dump_hdf5DataSet(i), numl)
        block
          character(len=:), allocatable :: colNames(:)
          character(len=:), allocatable :: colUnits(:)
          logical spErr
          integer nSplit
          call chr_split("napx turn s <x> <xp> <y> <yp> <sigma> <(E-E0)/E0> <x^2> <x*xp> <x*y> <x*yp> <x*sigma> <x*(E-E0)/E0> "//  &
            "<xp^2> <xp*y> <xp*yp> <xp*sigma> <xp*(E-E0)/E0> <y^2> <y*yp> <y*sigma> <y*(E-E0)/E0> <yp^2> <yp*sigma> "//            &
            "<yp*(E-E0)/E0> <sigma^2> <sigma*(E-E0)/E0> <((E-E0)/E0)^2>",colNames,nSplit,spErr)
          call chr_split("1 1 m mm mrad mm mrad mm 1 mm^2 mm*mrad mm^2 mm*mrad mm*mrad mm^2 mm mm^2 mrad*mm mrad^2 mrad*mm mrad "//&
            "mm^2 mm*mrad mm^2 mm mrad^2 mrad*mm mrad mm^2 mm 1",colUnits,nSplit,spErr)
          call h5_writeDataSetAttr(dump_hdf5DataSet(i),"dumpFormat",5)
          call h5_writeDataSetAttr(dump_hdf5DataSet(i),"colNames",  colNames)
          call h5_writeDataSetAttr(dump_hdf5DataSet(i),"colUnits",  colUnits)
        end block

      case(6)
        ! Format 6:
        ! # napx turn s[m] <x>[m] <px>[1] <y>[m] <py>[1] <sigma>[m] <psigma>[1]
        ! <x^2> <x*px> <x*y> <x*py> <x*sigma> <x*psigma>
        ! <px^2> <px*y> <px*py> <px*sigma> <px*psigma>
        ! <y^2> <y*py> <y*sigma> <y*psigma>
        ! <py^2> <py*sigma> <py*psigma>
        ! <sigma^2> <sigma*psigma>
        ! <psigma^2>
        if(dump_hdf5Format(6) == 0) then
          allocate(setFields(30))
          setFields(1)  = h5_dataField(name="NAPX",         type=h5_typeInt)
          setFields(2)  = h5_dataField(name="TURN",         type=h5_typeInt)
          setFields(3)  = h5_dataField(name="S",            type=h5_typeReal)
          setFields(4)  = h5_dataField(name="X",            type=h5_typeReal)
          setFields(5)  = h5_dataField(name="XP",           type=h5_typeReal)
          setFields(6)  = h5_dataField(name="Y",            type=h5_typeReal)
          setFields(7)  = h5_dataField(name="YP",           type=h5_typeReal)
          setFields(8)  = h5_dataField(name="Z",            type=h5_typeReal)
          setFields(9)  = h5_dataField(name="dE/E",         type=h5_typeReal)
          setFields(10) = h5_dataField(name="X^2",          type=h5_typeReal)
          setFields(11) = h5_dataField(name="X*PX",         type=h5_typeReal)
          setFields(12) = h5_dataField(name="X*Y",          type=h5_typeReal)
          setFields(13) = h5_dataField(name="X*PY",         type=h5_typeReal)
          setFields(14) = h5_dataField(name="X*SIGMA",      type=h5_typeReal)
          setFields(15) = h5_dataField(name="X*PSIGMA",     type=h5_typeReal)
          setFields(16) = h5_dataField(name="PX^2",         type=h5_typeReal)
          setFields(17) = h5_dataField(name="PX*Y",         type=h5_typeReal)
          setFields(18) = h5_dataField(name="PX*PY",        type=h5_typeReal)
          setFields(19) = h5_dataField(name="PX*SIGMA",     type=h5_typeReal)
          setFields(20) = h5_dataField(name="PX*PSIGMA",    type=h5_typeReal)
          setFields(21) = h5_dataField(name="Y^2",          type=h5_typeReal)
          setFields(22) = h5_dataField(name="Y*PY",         type=h5_typeReal)
          setFields(23) = h5_dataField(name="Y*SIGMA",      type=h5_typeReal)
          setFields(24) = h5_dataField(name="Y*PSIGMA",     type=h5_typeReal)
          setFields(25) = h5_dataField(name="PY^2",         type=h5_typeReal)
          setFields(26) = h5_dataField(name="PY*SIGMA",     type=h5_typeReal)
          setFields(27) = h5_dataField(name="PY*PSIGMA",    type=h5_typeReal)
          setFields(28) = h5_dataField(name="SIGMA^2",      type=h5_typeReal)
          setFields(29) = h5_dataField(name="SIGMA*PSIGMA", type=h5_typeReal)
          setFields(30) = h5_dataField(name="PSIGMA^2",     type=h5_typeReal)
          call h5_createFormat("dumpFormat6", setFields, dump_hdf5Format(6))
        end if
        call h5_createDataSet(dump_fname(i), h5_dumpID, dump_hdf5Format(6), dump_hdf5DataSet(i), numl)
        block
          character(len=:), allocatable :: colNames(:)
          character(len=:), allocatable :: colUnits(:)
          logical spErr
          integer nSplit
          call chr_split("napx turn s <x> <px> <y> <py> <sigma> <psigma> <x^2> <x*px> <x*y> <x*py> <x*sigma> <x*psigma> <px^2> "// &
            "<px*y> <px*py> <px*sigma> <px*psigma> <y^2> <y*py> <y*sigma> <y*psigma> <py^2> <py*sigma> <py*psigma> <sigma^2> "//   &
            "<sigma*psigma> <psigma^2>",colNames,nSplit,spErr)
          call chr_split("1 1 m m 1 m 1 m 1 mm^2 mm mm^2 mm mm^2 mm 1 mm 1 mm 1 mm^2 mm mm^2 mm 1 mm 1 mm^2 mm 1",&
            colUnits,nSplit,spErr)
          call h5_writeDataSetAttr(dump_hdf5DataSet(i),"dumpFormat",6)
          call h5_writeDataSetAttr(dump_hdf5DataSet(i),"colNames",  colNames)
          call h5_writeDataSetAttr(dump_hdf5DataSet(i),"colUnits",  colUnits)
        end block

      case(101)
        ! Format 101:
        if(dump_hdf5Format(101) == 0) then
          allocate(setFields(19))
          setFields(1)  = h5_dataField(name="ID",      type=h5_typeInt)
          setFields(2)  = h5_dataField(name="TURN",    type=h5_typeInt)
          setFields(3)  = h5_dataField(name="S",       type=h5_typeReal)
          setFields(4)  = h5_dataField(name="X",       type=h5_typeReal)
          setFields(5)  = h5_dataField(name="XP",      type=h5_typeReal)
          setFields(6)  = h5_dataField(name="Y",       type=h5_typeReal)
          setFields(7)  = h5_dataField(name="YP",      type=h5_typeReal)
          setFields(8)  = h5_dataField(name="dE/E",    type=h5_typeReal)
          setFields(9)  = h5_dataField(name="SIGMA",   type=h5_typeReal)
          setFields(10) = h5_dataField(name="KTRACK",  type=h5_typeInt)
          setFields(11) = h5_dataField(name="E",       type=h5_typeReal)
          setFields(12) = h5_dataField(name="PC",      type=h5_typeReal)
          setFields(13) = h5_dataField(name="P/P0",    type=h5_typeReal)
          setFields(14) = h5_dataField(name="P0/P)",   type=h5_typeReal)
          setFields(15) = h5_dataField(name="RV",      type=h5_typeReal)
          setFields(16) = h5_dataField(name="MASS",    type=h5_typeReal)
          setFields(17) = h5_dataField(name="MTC",     type=h5_typeReal)
          setFields(18) = h5_dataField(name="ENERGY0", type=h5_typeReal)
          setFields(19) = h5_dataField(name="PC0",     type=h5_typeReal)
          call h5_createFormat("dumpFormat101", setFields, dump_hdf5Format(101))
        end if
        call h5_createDataSet(dump_fname(i), h5_dumpID, dump_hdf5Format(101), dump_hdf5DataSet(i), napx)

      case(200)
        ! Format 200:
        ! # particleID turn s[m] x[mm] xp[mrad] y[mm] yp[mrad] sigma[mm] (E-E0)/E0[1] ktrack Z A Q PDGid
        if(dump_hdf5Format(200) == 0) then
          allocate(setFields(14))
          setFields(1)  = h5_dataField(name="ID",     type=h5_typeInt)
          setFields(2)  = h5_dataField(name="TURN",   type=h5_typeInt)
          setFields(3)  = h5_dataField(name="S",      type=h5_typeReal)
          setFields(4)  = h5_dataField(name="X",      type=h5_typeReal)
          setFields(5)  = h5_dataField(name="XP",     type=h5_typeReal)
          setFields(6)  = h5_dataField(name="Y",      type=h5_typeReal)
          setFields(7)  = h5_dataField(name="YP",     type=h5_typeReal)
          setFields(8)  = h5_dataField(name="SIGMA",  type=h5_typeReal)
          setFields(9)  = h5_dataField(name="dE/E",   type=h5_typeReal)
          setFields(10) = h5_dataField(name="KTRACK", type=h5_typeInt)
          setFields(11) = h5_dataField(name="Z",      type=h5_typeInt)
          setFields(12) = h5_dataField(name="A",      type=h5_typeInt)
          setFields(13) = h5_dataField(name="Q",      type=h5_typeInt)
          setFields(14) = h5_dataField(name="PDGid",  type=h5_typeInt)
          call h5_createFormat("dumpFormat200", setFields, dump_hdf5Format(200))
        end if
        call h5_createDataSet(dump_fname(i), h5_dumpID, dump_hdf5Format(200), dump_hdf5DataSet(i), napx)
        block
          character(len=:), allocatable :: colNames(:)
          character(len=:), allocatable :: colUnits(:)
          logical spErr
          integer nSplit
          call chr_split("ID turn s x xp y yp sigma (E-E0)/E0 ktrack Z A Q PDGid",colNames,nSplit,spErr)
          call chr_split("1 1 m mm mrad mm mrad mm 1 1 1 1 1 1",colUnits,nSplit,spErr)
          call h5_writeDataSetAttr(dump_hdf5DataSet(i),"dumpFormat",200)
          call h5_writeDataSetAttr(dump_hdf5DataSet(i),"colNames",  colNames)
          call h5_writeDataSetAttr(dump_hdf5DataSet(i),"colUnits",  colUnits)
        end block

      case(rootfmt)
        ! Format 1000: ROOT
        ! Nothing to do here!
      end select

      if(allocated(setFields) .eqv. .true.) deallocate(setFields)
    end if
  end do
#endif

end subroutine dump_initialise

! ================================================================================================================================ !
!  A.Mereghetti, D.Sinuela-Pastor & P.Garcia Ortega, for the FLUKA Team
!  K.Sjobak, A.Santamaria, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-07-13
!  dump beam particles
!  always in main code
!
!  nturn     : Current turn number
!  i         : Current structure element (0 for StartDUMP)
!  ix        : Corresponding single element (<0 for BLOC, only for ALL; 0 for StartDUMP)
!  unit      : Unit to dump from
!  fmt       : Dump output format (0/1/2/...)
!  lhighprec : High precission output y/n
! ================================================================================================================================ !
subroutine dump_beam_population(nturn, i, ix, unit, fmt, lhighprec, loc_clo, tasinv)

  use floatPrecision
  use physical_constants
  use numerical_constants
  use crcoall
  use string_tools

  use parpro
  use mod_common
  use mod_common_track
  use mod_common_main
  use mod_time

#ifdef HDF5
  use hdf5_output
#endif
#ifdef FLUKA
  use mod_fluka
#endif
#ifdef ROOT
  use root_output
#endif

  implicit none

  ! interface variables:
  integer, intent(in) :: nturn, i, ix, unit, fmt
  logical, intent(in) :: lhighprec
  real(kind=fPrec), intent(in) :: tasinv(6,6) ! normalization matrix in [mm,mrad,mm,mrad,mm,1]
  real(kind=fPrec), intent(in) :: loc_clo(6) ! closed orbit in [mm,mrad,mm,mrad,mm,1]

  ! Temporary variables
  integer j,k,l,m,n
  character(len=mNameLen) localBez
  logical rErr

  real(kind=fPrec) localDcum
  integer localKtrack

  real(kind=fPrec) xyz_particle(6),nxyz_particle(6)
  real(kind=fPrec) xyz(6)
  real(kind=fPrec) xyz2(6,6)

  character(len=16) :: xyz_l(27)
  character(len=25) :: xyz_h(27)

! temp variable for dealing with particle weights (only exist in FLUKA/G4 builds)
  real(kind=fPrec) tmpWeight

#ifdef CR
  ! For accessing dumpfilepos
  integer dumpIdx
  if(unit == dumpunit(0)) then
    ! ALL output must be on separate unit
    dumpIdx = 0
  else if(unit == dumpunit(-1)) then
    ! ALL output must be on separate unit
    dumpIdx = -1
  else
    dumpIdx = ix
  end if
#endif

  call time_startClock(time_clockDUMP)

  ! ------------------------------------------------------------------ !
  !  Format #0
  !  General format
  ! ------------------------------------------------------------------ !
  if(fmt == 0) then
    if(i == 0 .and. ix == 0) then
      localDcum = 0.0
      localBez  = "StartDUMP"
    else
      localDcum = dcum(i)
      if (ktrack(i) /= 1) then
        localBez = bez(ix)
      else                ! BLOCs
        localBez = bezb(ic(i))
      end if
    end if
    if(lhighprec) then
      do j=1,napx
        call chr_fromReal(xv1(j)*c1m3,                    xyz_h(1),19,2,rErr)
        call chr_fromReal(yv1(j)*c1m3,                    xyz_h(2),19,2,rErr)
        call chr_fromReal(xv2(j)*c1m3,                    xyz_h(3),19,2,rErr)
        call chr_fromReal(yv2(j)*c1m3,                    xyz_h(4),19,2,rErr)
        call chr_fromReal(ejfv(j)*c1m3,                    xyz_h(5),19,2,rErr)
        call chr_fromReal((ejv(j)-e0)*c1e6,                xyz_h(6),19,2,rErr)
        call chr_fromReal(-c1m3*(sigmv(j)/clight)*(e0/e0f),xyz_h(7),19,2,rErr)
        write(unit,"(3(1x,i8),1x,a48,1x,f12.5,7(1x,a25))") nturn,i,ix,localBez,localDcum, &
          xyz_h(1),xyz_h(2),xyz_h(3),xyz_h(4),xyz_h(5),xyz_h(6),xyz_h(7)
      end do
    else
      do j=1,napx
        call chr_fromReal(xv1(j)*c1m3,                    xyz_l(1),10,2,rErr)
        call chr_fromReal(yv1(j)*c1m3,                    xyz_l(2),10,2,rErr)
        call chr_fromReal(xv2(j)*c1m3,                    xyz_l(3),10,2,rErr)
        call chr_fromReal(yv2(j)*c1m3,                    xyz_l(4),10,2,rErr)
        call chr_fromReal(ejfv(j)*c1m3,                    xyz_l(5),10,2,rErr)
        call chr_fromReal((ejv(j)-e0)*c1e6,                xyz_l(6),10,2,rErr)
        call chr_fromReal(-c1m3*(sigmv(j)/clight)*(e0/e0f),xyz_l(7),10,2,rErr)
        write(unit,"(3(1x,i8),1x,a48,1x,f12.5,7(1x,a16))") nturn,i,ix,localBez,localDcum, &
          xyz_l(1),xyz_l(2),xyz_l(3),xyz_l(4),xyz_l(5),xyz_l(6),xyz_l(7)
      end do
    end if
    write(unit,"(a)") ""
    write(unit,"(a)") ""
    flush(unit,iostat=ierro)
#ifdef CR
    dumpfilepos(dumpIdx) = dumpfilepos(dumpIdx)+napx+2
#endif

  ! ------------------------------------------------------------------ !
  !  Format #1
  !  Format for aperture check
  ! ------------------------------------------------------------------ !
  else if (fmt == 1) then
    if (i == 0 .and. ix == 0) then
      localDcum = 0.0
      localKtrack = 0
    else
      localDcum = dcum(i)
      localKtrack = ktrack(i)
    end if
#ifdef HDF5
    if(h5_useForDUMP) then
      call h5_prepareWrite(dump_hdf5DataSet(ix), napx)
      call h5_writeData(dump_hdf5DataSet(ix), 1, napx, partID)
      call h5_writeData(dump_hdf5DataSet(ix), 2, napx, nturn)
      call h5_writeData(dump_hdf5DataSet(ix), 3, napx, localDcum)
      call h5_writeData(dump_hdf5DataSet(ix), 4, napx, xv1(:))
      call h5_writeData(dump_hdf5DataSet(ix), 5, napx, yv1(:))
      call h5_writeData(dump_hdf5DataSet(ix), 6, napx, xv2(:))
      call h5_writeData(dump_hdf5DataSet(ix), 7, napx, yv2(:))
      call h5_writeData(dump_hdf5DataSet(ix), 8, napx, (ejv-e0)/e0)
      call h5_writeData(dump_hdf5DataSet(ix), 9, napx, localKtrack)
      call h5_finaliseWrite(dump_hdf5DataSet(ix))
    else
#endif
      if(lhighprec) then
        do j=1,napx
          call chr_fromReal(xv1(j),       xyz_h(1),19,2,rErr)
          call chr_fromReal(yv1(j),       xyz_h(2),19,2,rErr)
          call chr_fromReal(xv2(j),       xyz_h(3),19,2,rErr)
          call chr_fromReal(yv2(j),       xyz_h(4),19,2,rErr)
          call chr_fromReal((ejv(j)-e0)/e0,xyz_h(5),19,2,rErr)
          write(unit,"(2(1x,i8),1x,f12.5,5(1x,a25),1x,i8)") partID(j),nturn,localDcum, &
            xyz_h(1),xyz_h(2),xyz_h(3),xyz_h(4),xyz_h(5),localKtrack
        end do
      else
        do j=1,napx
          call chr_fromReal(xv1(j),       xyz_l(1),10,2,rErr)
          call chr_fromReal(yv1(j),       xyz_l(2),10,2,rErr)
          call chr_fromReal(xv2(j),       xyz_l(3),10,2,rErr)
          call chr_fromReal(yv2(j),       xyz_l(4),10,2,rErr)
          call chr_fromReal((ejv(j)-e0)/e0,xyz_l(5),10,2,rErr)
          write(unit,"(2(1x,i8),1x,f12.5,5(1x,a16),1x,i8)") partID(j),nturn,localDcum, &
            xyz_l(1),xyz_l(2),xyz_l(3),xyz_l(4),xyz_l(5),localKtrack
        end do
      end if
      flush(unit,iostat=ierro)
#ifdef CR
      dumpfilepos(dumpIdx) = dumpfilepos(dumpIdx)+napx
#endif
#ifdef HDF5
    end if
#endif

  ! ------------------------------------------------------------------ !
  !  Format #2
  !  Same as fmt 1, but also include z (for crab cavities etc.)
  ! ------------------------------------------------------------------ !
  else if(fmt == 2) then
    if(i == 0 .and. ix == 0) then
      localDcum   = zero
      localKtrack = 0
    else
      localDcum   = dcum(i)
      localKtrack = ktrack(i)
    end if
#ifdef HDF5
    if(h5_useForDUMP) then
      call h5_prepareWrite(dump_hdf5DataSet(ix), napx)
      call h5_writeData(dump_hdf5DataSet(ix), 1,  napx, partID)
      call h5_writeData(dump_hdf5DataSet(ix), 2,  napx, nturn)
      call h5_writeData(dump_hdf5DataSet(ix), 3,  napx, localDcum)
      call h5_writeData(dump_hdf5DataSet(ix), 4,  napx, xv1(:))
      call h5_writeData(dump_hdf5DataSet(ix), 5,  napx, yv1(:))
      call h5_writeData(dump_hdf5DataSet(ix), 6,  napx, xv2(:))
      call h5_writeData(dump_hdf5DataSet(ix), 7,  napx, yv2(:))
      call h5_writeData(dump_hdf5DataSet(ix), 8,  napx, sigmv)
      call h5_writeData(dump_hdf5DataSet(ix), 9,  napx, (ejv-e0)/e0)
      call h5_writeData(dump_hdf5DataSet(ix), 10, napx, localKtrack)
      call h5_finaliseWrite(dump_hdf5DataSet(ix))
    else
#endif
      if(lhighprec) then
        do j=1,napx
          call chr_fromReal(xv1(j),       xyz_h(1),19,2,rErr)
          call chr_fromReal(yv1(j),       xyz_h(2),19,2,rErr)
          call chr_fromReal(xv2(j),       xyz_h(3),19,2,rErr)
          call chr_fromReal(yv2(j),       xyz_h(4),19,2,rErr)
          call chr_fromReal(sigmv(j),      xyz_h(5),19,2,rErr)
          call chr_fromReal((ejv(j)-e0)/e0,xyz_h(6),19,2,rErr)
          write(unit,"(2(1x,i8),1x,f12.5,6(1x,a25),1x,i8)") partID(j),nturn,localDcum,&
            xyz_h(1),xyz_h(2),xyz_h(3),xyz_h(4),xyz_h(5),xyz_h(6),localKtrack
        end do
      else
        do j=1,napx
          call chr_fromReal(xv1(j),       xyz_l(1),10,2,rErr)
          call chr_fromReal(yv1(j),       xyz_l(2),10,2,rErr)
          call chr_fromReal(xv2(j),       xyz_l(3),10,2,rErr)
          call chr_fromReal(yv2(j),       xyz_l(4),10,2,rErr)
          call chr_fromReal(sigmv(j),      xyz_l(5),10,2,rErr)
          call chr_fromReal((ejv(j)-e0)/e0,xyz_l(6),10,2,rErr)
          write(unit,"(2(1x,i8),1x,f12.5,6(1x,a16),1x,i8)") partID(j),nturn,localDcum,&
            xyz_l(1),xyz_l(2),xyz_l(3),xyz_l(4),xyz_l(5),xyz_l(6),localKtrack
        end do
      end if
      flush(unit,iostat=ierro)
#ifdef CR
      dumpfilepos(dumpIdx) = dumpfilepos(dumpIdx)+napx
#endif
#ifdef HDF5
    end if
#endif

  ! ------------------------------------------------------------------ !
  !  Format #3
  !  Same as fmt 2, but in Fortran binary
  ! ------------------------------------------------------------------ !
  else if(fmt == 3) then
    if(i == 0 .and. ix == 0) then
      localDcum   = zero
      localKtrack = 0
    else
      localDcum   = dcum(i)
      localKtrack = ktrack(i)
    end if
#ifdef HDF5
    if(h5_useForDUMP) then
      call h5_prepareWrite(dump_hdf5DataSet(ix), napx)
      call h5_writeData(dump_hdf5DataSet(ix), 1,  napx, partID)
      call h5_writeData(dump_hdf5DataSet(ix), 2,  napx, nturn)
      call h5_writeData(dump_hdf5DataSet(ix), 3,  napx, localDcum)
      call h5_writeData(dump_hdf5DataSet(ix), 4,  napx, xv1(:))
      call h5_writeData(dump_hdf5DataSet(ix), 5,  napx, yv1(:))
      call h5_writeData(dump_hdf5DataSet(ix), 6,  napx, xv2(:))
      call h5_writeData(dump_hdf5DataSet(ix), 7,  napx, yv2(:))
      call h5_writeData(dump_hdf5DataSet(ix), 8,  napx, sigmv)
      call h5_writeData(dump_hdf5DataSet(ix), 9,  napx, (ejv-e0)/e0)
      call h5_writeData(dump_hdf5DataSet(ix), 10, napx, localKtrack)
      call h5_finaliseWrite(dump_hdf5DataSet(ix))
    else
#endif
      do j=1,napx
        write(unit) partID(j),nturn,localDcum,xv1(j),yv1(j),xv2(j),yv2(j), &
          sigmv(j),(ejv(j)-e0)/e0,localKtrack
      end do
      flush(unit,iostat=ierro)
#ifdef CR
      dumpfilepos(dumpIdx) = dumpfilepos(dumpIdx)+napx
#endif
#ifdef HDF5
    end if
#endif

  ! ------------------------------------------------------------------ !
  !  Format #4
  !  Average bunch position
  ! ------------------------------------------------------------------ !
  else if(fmt == 4) then
    if(i == 0 .and. ix == 0) then
      localDcum = zero
    else
      localDcum = dcum(i)
    end if
    xyz(:) = zero
    do j=1,napx
      xyz(1) = xyz(1) + xv1(j)
      xyz(2) = xyz(2) + yv1(j)
      xyz(3) = xyz(3) + xv2(j)
      xyz(4) = xyz(4) + yv2(j)
      xyz(5) = xyz(5) + sigmv(j)
      xyz(6) = xyz(6) + (ejv(j)-e0)/e0
    end do
    xyz = xyz/napx

#ifdef HDF5
    if(h5_useForDUMP) then
      call h5_prepareWrite(dump_hdf5DataSet(ix), 1)
      call h5_writeData(dump_hdf5DataSet(ix), 1, 1, napx)
      call h5_writeData(dump_hdf5DataSet(ix), 2, 1, nturn)
      call h5_writeData(dump_hdf5DataSet(ix), 3, 1, localDcum)
      call h5_writeData(dump_hdf5DataSet(ix), 4, 1, xyz(1))
      call h5_writeData(dump_hdf5DataSet(ix), 5, 1, xyz(2))
      call h5_writeData(dump_hdf5DataSet(ix), 6, 1, xyz(3))
      call h5_writeData(dump_hdf5DataSet(ix), 7, 1, xyz(4))
      call h5_writeData(dump_hdf5DataSet(ix), 8, 1, xyz(5))
      call h5_writeData(dump_hdf5DataSet(ix), 9, 1, xyz(6))
      call h5_finaliseWrite(dump_hdf5DataSet(ix))
    else
#endif
      if(lhighprec) then
        call chr_fromReal(xyz(1),xyz_h(1),19,2,rErr)
        call chr_fromReal(xyz(2),xyz_h(2),19,2,rErr)
        call chr_fromReal(xyz(3),xyz_h(3),19,2,rErr)
        call chr_fromReal(xyz(4),xyz_h(4),19,2,rErr)
        call chr_fromReal(xyz(5),xyz_h(5),19,2,rErr)
        call chr_fromReal(xyz(6),xyz_h(6),19,2,rErr)
        write(unit,"(2(1x,i8),1x,f12.5,6(1x,a25))") napx,nturn,localDcum,xyz_h(1),xyz_h(2),xyz_h(3),xyz_h(4),xyz_h(5),xyz_h(6)
      else
        call chr_fromReal(xyz(1),xyz_l(1),10,2,rErr)
        call chr_fromReal(xyz(2),xyz_l(2),10,2,rErr)
        call chr_fromReal(xyz(3),xyz_l(3),10,2,rErr)
        call chr_fromReal(xyz(4),xyz_l(4),10,2,rErr)
        call chr_fromReal(xyz(5),xyz_l(5),10,2,rErr)
        call chr_fromReal(xyz(6),xyz_l(6),10,2,rErr)
        write(unit,"(2(1x,i8),1x,f12.5,6(1x,a16))") napx,nturn,localDcum,xyz_l(1),xyz_l(2),xyz_l(3),xyz_l(4),xyz_l(5),xyz_l(6)
      end if
      flush(unit,iostat=ierro)
#ifdef CR
      dumpfilepos(dumpIdx) = dumpfilepos(dumpIdx)+1
#endif
#ifdef HDF5
    end if
#endif

  ! ------------------------------------------------------------------ !
  !  Format #5 or #6
  !  Average beam positon + beam matrix
  ! ------------------------------------------------------------------ !
  else if (fmt == 5 .or. fmt == 6) then
    if (i == 0 .and. ix == 0) then
      localDcum = zero
    else
      localDcum = dcum(i)
    end if

    do l=1,6
      xyz(l) = zero
      do k=1,6
        xyz2(l,k) = zero
      end do
    end do

    if(fmt == 5) then ! Raw
      do j=1,napx
        xyz_particle(6)=(ejv(j)-e0)/e0

        ! Average beam position
        xyz(1) = xyz(1) + xv1(j)
        xyz(2) = xyz(2) + yv1(j)
        xyz(3) = xyz(3) + xv2(j)
        xyz(4) = xyz(4) + yv2(j)
        xyz(5) = xyz(5) + sigmv(j)
        xyz(6) = xyz(6) + xyz_particle(6)

        ! Beam matrix (don't calulate identical elements twice (symmetry))
        xyz2(1,1) = xyz2(1,1) + xv1(j)*xv1(j)
        xyz2(2,1) = xyz2(2,1) + xv1(j)*yv1(j)
        xyz2(3,1) = xyz2(3,1) + xv1(j)*xv2(j)
        xyz2(4,1) = xyz2(4,1) + xv1(j)*yv2(j)
        xyz2(5,1) = xyz2(5,1) + xv1(j)*sigmv(j)
        xyz2(6,1) = xyz2(6,1) + xv1(j)*xyz_particle(6)

        xyz2(2,2) = xyz2(2,2) + yv1(j)*yv1(j)
        xyz2(3,2) = xyz2(3,2) + yv1(j)*xv2(j)
        xyz2(4,2) = xyz2(4,2) + yv1(j)*yv2(j)
        xyz2(5,2) = xyz2(5,2) + yv1(j)*sigmv(j)
        xyz2(6,2) = xyz2(6,2) + yv1(j)*xyz_particle(6)

        xyz2(3,3) = xyz2(3,3) + xv2(j)*xv2(j)
        xyz2(4,3) = xyz2(4,3) + xv2(j)*yv2(j)
        xyz2(5,3) = xyz2(5,3) + xv2(j)*sigmv(j)
        xyz2(6,3) = xyz2(6,3) + xv2(j)*xyz_particle(6)

        xyz2(4,4) = xyz2(4,4) + yv2(j)*yv2(j)
        xyz2(5,4) = xyz2(5,4) + yv2(j)*sigmv(j)
        xyz2(6,4) = xyz2(6,4) + yv2(j)*xyz_particle(6)

        xyz2(5,5) = xyz2(5,5) + sigmv(j)*sigmv(j)
        xyz2(6,5) = xyz2(6,5) + sigmv(j)*xyz_particle(6)

        xyz2(6,6) = xyz2(6,6) + xyz_particle(6)*xyz_particle(6)
      end do

    else if (fmt == 6) then ! Canonical
      do j=1,napx
        xyz_particle(1) = xv1(j)*c1m3                 !x:      [mm]   -> [m]
        xyz_particle(2) = (yv1(j)*c1m3)*(one+dpsv(j)) !px:     [mrad] -> [1]
        xyz_particle(3) = xv2(j)*c1m3                 !y:      [mm]   -> [m]
        xyz_particle(4) = (yv2(j)*c1m3)*(one+dpsv(j)) !py:     [mrad] -> [1]
        xyz_particle(5) = sigmv(j)*c1m3                !sigma:  [mm]   -> [m]
        xyz_particle(6) = (((ejv(j)-e0)*e0)/e0f)/e0f   !psigma: [MeV]  -> [1]

        ! Average beam position
        xyz(1) = xyz(1) + xyz_particle(1)
        xyz(2) = xyz(2) + xyz_particle(2)
        xyz(3) = xyz(3) + xyz_particle(3)
        xyz(4) = xyz(4) + xyz_particle(4)
        xyz(5) = xyz(5) + xyz_particle(5)
        xyz(6) = xyz(6) + xyz_particle(6)

        ! Beam matrix (don't calulate identical elements twice (symmetry))
        xyz2(1,1) = xyz2(1,1) + xyz_particle(1)*xyz_particle(1)
        xyz2(2,1) = xyz2(2,1) + xyz_particle(1)*xyz_particle(2)
        xyz2(3,1) = xyz2(3,1) + xyz_particle(1)*xyz_particle(3)
        xyz2(4,1) = xyz2(4,1) + xyz_particle(1)*xyz_particle(4)
        xyz2(5,1) = xyz2(5,1) + xyz_particle(1)*xyz_particle(5)
        xyz2(6,1) = xyz2(6,1) + xyz_particle(1)*xyz_particle(6)

        xyz2(2,2) = xyz2(2,2) + xyz_particle(2)*xyz_particle(2)
        xyz2(3,2) = xyz2(3,2) + xyz_particle(2)*xyz_particle(3)
        xyz2(4,2) = xyz2(4,2) + xyz_particle(2)*xyz_particle(4)
        xyz2(5,2) = xyz2(5,2) + xyz_particle(2)*xyz_particle(5)
        xyz2(6,2) = xyz2(6,2) + xyz_particle(2)*xyz_particle(6)

        xyz2(3,3) = xyz2(3,3) + xyz_particle(3)*xyz_particle(3)
        xyz2(4,3) = xyz2(4,3) + xyz_particle(3)*xyz_particle(4)
        xyz2(5,3) = xyz2(5,3) + xyz_particle(3)*xyz_particle(5)
        xyz2(6,3) = xyz2(6,3) + xyz_particle(3)*xyz_particle(6)

        xyz2(4,4) = xyz2(4,4) + xyz_particle(4)*xyz_particle(4)
        xyz2(5,4) = xyz2(5,4) + xyz_particle(4)*xyz_particle(5)
        xyz2(6,4) = xyz2(6,4) + xyz_particle(4)*xyz_particle(6)

        xyz2(5,5) = xyz2(5,5) + xyz_particle(5)*xyz_particle(5)
        xyz2(6,5) = xyz2(6,5) + xyz_particle(5)*xyz_particle(6)

        xyz2(6,6) = xyz2(6,6) + xyz_particle(6)*xyz_particle(6)
      end do
    end if

    ! Normalize to get averages
    xyz = xyz/napx

    xyz2(:,1)  = xyz2(:,1) /napx
    xyz2(2:,2) = xyz2(2:,2)/napx
    xyz2(3:,3) = xyz2(3:,3)/napx
    xyz2(4:,4) = xyz2(4:,4)/napx
    xyz2(5:,5) = xyz2(5:,5)/napx
    xyz2(6,6)  = xyz2(6,6) /napx

#ifdef HDF5
    if(h5_useForDUMP) then
      call h5_prepareWrite(dump_hdf5DataSet(ix), 1)
      call h5_writeData(dump_hdf5DataSet(ix), 1,  1, napx)
      call h5_writeData(dump_hdf5DataSet(ix), 2,  1, nturn)
      call h5_writeData(dump_hdf5DataSet(ix), 3,  1, localDcum)
      call h5_writeData(dump_hdf5DataSet(ix), 4,  1, xyz(1))
      call h5_writeData(dump_hdf5DataSet(ix), 5,  1, xyz(2))
      call h5_writeData(dump_hdf5DataSet(ix), 6,  1, xyz(3))
      call h5_writeData(dump_hdf5DataSet(ix), 7,  1, xyz(4))
      call h5_writeData(dump_hdf5DataSet(ix), 8,  1, xyz(5))
      call h5_writeData(dump_hdf5DataSet(ix), 9,  1, xyz(6))
      call h5_writeData(dump_hdf5DataSet(ix), 10, 1, xyz2(1,1))
      call h5_writeData(dump_hdf5DataSet(ix), 11, 1, xyz2(2,1))
      call h5_writeData(dump_hdf5DataSet(ix), 12, 1, xyz2(3,1))
      call h5_writeData(dump_hdf5DataSet(ix), 13, 1, xyz2(4,1))
      call h5_writeData(dump_hdf5DataSet(ix), 14, 1, xyz2(5,1))
      call h5_writeData(dump_hdf5DataSet(ix), 15, 1, xyz2(6,1))
      call h5_writeData(dump_hdf5DataSet(ix), 16, 1, xyz2(2,2))
      call h5_writeData(dump_hdf5DataSet(ix), 17, 1, xyz2(3,2))
      call h5_writeData(dump_hdf5DataSet(ix), 18, 1, xyz2(4,2))
      call h5_writeData(dump_hdf5DataSet(ix), 19, 1, xyz2(5,2))
      call h5_writeData(dump_hdf5DataSet(ix), 20, 1, xyz2(6,2))
      call h5_writeData(dump_hdf5DataSet(ix), 21, 1, xyz2(3,3))
      call h5_writeData(dump_hdf5DataSet(ix), 22, 1, xyz2(4,3))
      call h5_writeData(dump_hdf5DataSet(ix), 23, 1, xyz2(5,3))
      call h5_writeData(dump_hdf5DataSet(ix), 24, 1, xyz2(6,3))
      call h5_writeData(dump_hdf5DataSet(ix), 25, 1, xyz2(4,4))
      call h5_writeData(dump_hdf5DataSet(ix), 26, 1, xyz2(5,4))
      call h5_writeData(dump_hdf5DataSet(ix), 27, 1, xyz2(6,4))
      call h5_writeData(dump_hdf5DataSet(ix), 28, 1, xyz2(5,5))
      call h5_writeData(dump_hdf5DataSet(ix), 29, 1, xyz2(6,5))
      call h5_writeData(dump_hdf5DataSet(ix), 30, 1, xyz2(6,6))
call h5_finaliseWrite(dump_hdf5DataSet(ix))
    else
#endif
      if (lhighprec) then
        call chr_fromReal(xyz(1),xyz_h(1),19,2,rErr)
        call chr_fromReal(xyz(2),xyz_h(2),19,2,rErr)
        call chr_fromReal(xyz(3),xyz_h(3),19,2,rErr)
        call chr_fromReal(xyz(4),xyz_h(4),19,2,rErr)
        call chr_fromReal(xyz(5),xyz_h(5),19,2,rErr)
        call chr_fromReal(xyz(6),xyz_h(6),19,2,rErr)

        call chr_fromReal(xyz2(1,1),xyz_h(7), 19,2,rErr)
        call chr_fromReal(xyz2(2,1),xyz_h(8), 19,2,rErr)
        call chr_fromReal(xyz2(3,1),xyz_h(9), 19,2,rErr)
        call chr_fromReal(xyz2(4,1),xyz_h(10),19,2,rErr)
        call chr_fromReal(xyz2(5,1),xyz_h(11),19,2,rErr)
        call chr_fromReal(xyz2(6,1),xyz_h(12),19,2,rErr)

        call chr_fromReal(xyz2(2,2),xyz_h(13),19,2,rErr)
        call chr_fromReal(xyz2(3,2),xyz_h(14),19,2,rErr)
        call chr_fromReal(xyz2(4,2),xyz_h(15),19,2,rErr)
        call chr_fromReal(xyz2(5,2),xyz_h(16),19,2,rErr)
        call chr_fromReal(xyz2(6,2),xyz_h(17),19,2,rErr)

        call chr_fromReal(xyz2(3,3),xyz_h(18),19,2,rErr)
        call chr_fromReal(xyz2(4,3),xyz_h(19),19,2,rErr)
        call chr_fromReal(xyz2(5,3),xyz_h(20),19,2,rErr)
        call chr_fromReal(xyz2(6,3),xyz_h(21),19,2,rErr)

        call chr_fromReal(xyz2(4,4),xyz_h(22),19,2,rErr)
        call chr_fromReal(xyz2(5,4),xyz_h(23),19,2,rErr)
        call chr_fromReal(xyz2(6,4),xyz_h(24),19,2,rErr)

        call chr_fromReal(xyz2(5,5),xyz_h(25),19,2,rErr)
        call chr_fromReal(xyz2(6,5),xyz_h(26),19,2,rErr)

        call chr_fromReal(xyz2(6,6),xyz_h(27),19,2,rErr)

        write(unit, "(2(1x,i8),1x,f12.5,27(1x,a25))") napx,nturn,localDcum, &
          xyz_h(1), xyz_h(2), xyz_h(3), xyz_h(4), xyz_h(5), xyz_h(6), xyz_h(7), xyz_h(8), xyz_h(9), &
          xyz_h(10),xyz_h(11),xyz_h(12),xyz_h(13),xyz_h(14),xyz_h(15),xyz_h(16),xyz_h(17),xyz_h(18),&
          xyz_h(19),xyz_h(20),xyz_h(21),xyz_h(22),xyz_h(23),xyz_h(24),xyz_h(25),xyz_h(26),xyz_h(27)
      else
        call chr_fromReal(xyz(1),xyz_l(1),10,2,rErr)
        call chr_fromReal(xyz(2),xyz_l(2),10,2,rErr)
        call chr_fromReal(xyz(3),xyz_l(3),10,2,rErr)
        call chr_fromReal(xyz(4),xyz_l(4),10,2,rErr)
        call chr_fromReal(xyz(5),xyz_l(5),10,2,rErr)
        call chr_fromReal(xyz(6),xyz_l(6),10,2,rErr)

        call chr_fromReal(xyz2(1,1),xyz_l(7), 10,2,rErr)
        call chr_fromReal(xyz2(2,1),xyz_l(8), 10,2,rErr)
        call chr_fromReal(xyz2(3,1),xyz_l(9), 10,2,rErr)
        call chr_fromReal(xyz2(4,1),xyz_l(10),10,2,rErr)
        call chr_fromReal(xyz2(5,1),xyz_l(11),10,2,rErr)
        call chr_fromReal(xyz2(6,1),xyz_l(12),10,2,rErr)

        call chr_fromReal(xyz2(2,2),xyz_l(13),10,2,rErr)
        call chr_fromReal(xyz2(3,2),xyz_l(14),10,2,rErr)
        call chr_fromReal(xyz2(4,2),xyz_l(15),10,2,rErr)
        call chr_fromReal(xyz2(5,2),xyz_l(16),10,2,rErr)
        call chr_fromReal(xyz2(6,2),xyz_l(17),10,2,rErr)

        call chr_fromReal(xyz2(3,3),xyz_l(18),10,2,rErr)
        call chr_fromReal(xyz2(4,3),xyz_l(19),10,2,rErr)
        call chr_fromReal(xyz2(5,3),xyz_l(20),10,2,rErr)
        call chr_fromReal(xyz2(6,3),xyz_l(21),10,2,rErr)

        call chr_fromReal(xyz2(4,4),xyz_l(22),10,2,rErr)
        call chr_fromReal(xyz2(5,4),xyz_l(23),10,2,rErr)
        call chr_fromReal(xyz2(6,4),xyz_l(24),10,2,rErr)

        call chr_fromReal(xyz2(5,5),xyz_l(25),10,2,rErr)
        call chr_fromReal(xyz2(6,5),xyz_l(26),10,2,rErr)

        call chr_fromReal(xyz2(6,6),xyz_l(27),10,2,rErr)

        write(unit, "(2(1x,i8),1x,f12.5,27(1x,a16))") napx,nturn,localDcum, &
          xyz_l(1), xyz_l(2), xyz_l(3), xyz_l(4), xyz_l(5), xyz_l(6), xyz_l(7), xyz_l(8), xyz_l(9), &
          xyz_l(10),xyz_l(11),xyz_l(12),xyz_l(13),xyz_l(14),xyz_l(15),xyz_l(16),xyz_l(17),xyz_l(18),&
          xyz_l(19),xyz_l(20),xyz_l(21),xyz_l(22),xyz_l(23),xyz_l(24),xyz_l(25),xyz_l(26),xyz_l(27)
      end if
      flush(unit,iostat=ierro)
#ifdef CR
      dumpfilepos(dumpIdx) = dumpfilepos(dumpIdx)+1
#endif
#ifdef HDF5
    end if
#endif

  ! ------------------------------------------------------------------ !
  !  Format #7, #8 or #9
  !  fmt 7 same as fmt 2,   but in normalized coordinates
  !  fmt 8 same as fmt 3,   but in normalized coordinates
  !  fmt 9 same as fmt 5/6, but in normalized coordinates
  ! ------------------------------------------------------------------ !
  else if(fmt == 7 .or. fmt == 8 .or. fmt == 9) then
    if(i == 0 .and. ix == 0) then
      localDcum   = zero
      localKtrack = 0
    else
      localDcum   = dcum(i)
      localKtrack = ktrack(i)
    end if

    ! initialize parameters for writing of beam moments
    xyz(:)    = zero
    xyz2(:,:) = zero

    ! normalize particle coordinates
    do j=1,napx
      xyz_particle(1) = xv1(j)
      xyz_particle(2) = yv1(j)
      xyz_particle(3) = xv2(j)
      xyz_particle(4) = yv2(j)
      xyz_particle(5) = sigmv(j)
      xyz_particle(6) = (ejv(j)-e0)/e0
      ! Remove closed orbit -> check units used in dumpclo (is x' or px used?)
      do m=1,6
        xyz_particle(m)=xyz_particle(m)-loc_clo(m)
      end do
      ! Convert to canonical variables
      xyz_particle(2)=xyz_particle(2)*((one+xyz_particle(6))+loc_clo(6))
      xyz_particle(4)=xyz_particle(4)*((one+xyz_particle(6))+loc_clo(6))
      ! Normalize nxyz=fma_tas_inv*xyz
      nxyz_particle(:)=zero
      do m=1,6
        do n=1,6
          nxyz_particle(m)=nxyz_particle(m)+tasinv(m,n)*xyz_particle(n)
        end do
        ! a) convert nxyzv(6) to 1.e-3 sqrt(m)
        ! unit: nx,npx,ny,npy,nsig,ndelta all in [1.e-3 sqrt(m)]
        if (m == 6) then
          nxyz_particle(m)=nxyz_particle(m)*c1e3
        end if
      end do

      if(fmt == 7) then
        if(lhighprec) then
          call chr_fromReal(nxyz_particle(1),xyz_h(1),19,2,rErr)
          call chr_fromReal(nxyz_particle(2),xyz_h(2),19,2,rErr)
          call chr_fromReal(nxyz_particle(3),xyz_h(3),19,2,rErr)
          call chr_fromReal(nxyz_particle(4),xyz_h(4),19,2,rErr)
          call chr_fromReal(nxyz_particle(5),xyz_h(5),19,2,rErr)
          call chr_fromReal(nxyz_particle(6),xyz_h(6),19,2,rErr)
          write(unit,"(2(1x,i8),1x,f12.5,6(1x,a25),1x,i8)") partID(j),nturn,localDcum, &
            xyz_h(1),xyz_h(2),xyz_h(3),xyz_h(4),xyz_h(5),xyz_h(6),localKtrack
        else
          call chr_fromReal(nxyz_particle(1),xyz_l(1),10,2,rErr)
          call chr_fromReal(nxyz_particle(2),xyz_l(2),10,2,rErr)
          call chr_fromReal(nxyz_particle(3),xyz_l(3),10,2,rErr)
          call chr_fromReal(nxyz_particle(4),xyz_l(4),10,2,rErr)
          call chr_fromReal(nxyz_particle(5),xyz_l(5),10,2,rErr)
          call chr_fromReal(nxyz_particle(6),xyz_l(6),10,2,rErr)
          write(unit,"(2(1x,i8),1x,f12.5,6(1x,a16),1x,i8)") partID(j),nturn,localDcum, &
            xyz_l(1),xyz_l(2),xyz_l(3),xyz_l(4),xyz_l(5),xyz_l(6),localKtrack
        end if

      else if(fmt == 8) then
        write(unit) partID(j),nturn,localDcum, &
          nxyz_particle(1),nxyz_particle(2),nxyz_particle(3),nxyz_particle(4),nxyz_particle(5),nxyz_particle(6),localKtrack

      else if(fmt == 9) then
        ! Average beam position
        ! here we recycle xyz used also for fmt 5 and 6. These are
        ! all normalized coordinates in units
        ! nx,npx,ny,npy,nsig,ndelta [1.e-3 sqrt(m)]
        xyz(1) = xyz(1) + nxyz_particle(1)
        xyz(2) = xyz(2) + nxyz_particle(2)
        xyz(3) = xyz(3) + nxyz_particle(3)
        xyz(4) = xyz(4) + nxyz_particle(4)
        xyz(5) = xyz(5) + nxyz_particle(5)
        xyz(6) = xyz(6) + nxyz_particle(6)

        ! Beam matrix (don't calulate identical elements twice (symmetry))
        xyz2(1,1) = xyz2(1,1) + nxyz_particle(1)*nxyz_particle(1)
        xyz2(2,1) = xyz2(2,1) + nxyz_particle(1)*nxyz_particle(2)
        xyz2(3,1) = xyz2(3,1) + nxyz_particle(1)*nxyz_particle(3)
        xyz2(4,1) = xyz2(4,1) + nxyz_particle(1)*nxyz_particle(4)
        xyz2(5,1) = xyz2(5,1) + nxyz_particle(1)*nxyz_particle(5)
        xyz2(6,1) = xyz2(6,1) + nxyz_particle(1)*nxyz_particle(6)

        xyz2(2,2) = xyz2(2,2) + nxyz_particle(2)*nxyz_particle(2)
        xyz2(3,2) = xyz2(3,2) + nxyz_particle(2)*nxyz_particle(3)
        xyz2(4,2) = xyz2(4,2) + nxyz_particle(2)*nxyz_particle(4)
        xyz2(5,2) = xyz2(5,2) + nxyz_particle(2)*nxyz_particle(5)
        xyz2(6,2) = xyz2(6,2) + nxyz_particle(2)*nxyz_particle(6)

        xyz2(3,3) = xyz2(3,3) + nxyz_particle(3)*nxyz_particle(3)
        xyz2(4,3) = xyz2(4,3) + nxyz_particle(3)*nxyz_particle(4)
        xyz2(5,3) = xyz2(5,3) + nxyz_particle(3)*nxyz_particle(5)
        xyz2(6,3) = xyz2(6,3) + nxyz_particle(3)*nxyz_particle(6)

        xyz2(4,4) = xyz2(4,4) + nxyz_particle(4)*nxyz_particle(4)
        xyz2(5,4) = xyz2(5,4) + nxyz_particle(4)*nxyz_particle(5)
        xyz2(6,4) = xyz2(6,4) + nxyz_particle(4)*nxyz_particle(6)

        xyz2(5,5) = xyz2(5,5) + nxyz_particle(5)*nxyz_particle(5)
        xyz2(6,5) = xyz2(6,5) + nxyz_particle(5)*nxyz_particle(6)

        xyz2(6,6) = xyz2(6,6) + nxyz_particle(6)*nxyz_particle(6)
      end if
    end do ! END loop over particles (j)

    if(fmt == 7) then
      flush(unit,iostat=ierro)
#ifdef CR
      dumpfilepos(dumpIdx) = dumpfilepos(dumpIdx)+napx
#endif
    else if(fmt == 8) then
      flush(unit,iostat=ierro)
#ifdef CR
      dumpfilepos(dumpIdx) = dumpfilepos(dumpIdx)+napx
#endif
    else if(fmt == 9) then
      ! Normalize to get averages
      xyz = xyz/napx

      xyz2(:,1)  = xyz2(:,1) /napx
      xyz2(2:,2) = xyz2(2:,2)/napx
      xyz2(3:,3) = xyz2(3:,3)/napx
      xyz2(4:,4) = xyz2(4:,4)/napx
      xyz2(5:,5) = xyz2(5:,5)/napx
      xyz2(6,6)  = xyz2(6,6) /napx

      if (lhighprec) then
        call chr_fromReal(xyz(1),xyz_h(1),19,2,rErr)
        call chr_fromReal(xyz(2),xyz_h(2),19,2,rErr)
        call chr_fromReal(xyz(3),xyz_h(3),19,2,rErr)
        call chr_fromReal(xyz(4),xyz_h(4),19,2,rErr)
        call chr_fromReal(xyz(5),xyz_h(5),19,2,rErr)
        call chr_fromReal(xyz(6),xyz_h(6),19,2,rErr)

        call chr_fromReal(xyz2(1,1),xyz_h(7), 19,2,rErr)
        call chr_fromReal(xyz2(2,1),xyz_h(8), 19,2,rErr)
        call chr_fromReal(xyz2(3,1),xyz_h(9), 19,2,rErr)
        call chr_fromReal(xyz2(4,1),xyz_h(10),19,2,rErr)
        call chr_fromReal(xyz2(5,1),xyz_h(11),19,2,rErr)
        call chr_fromReal(xyz2(6,1),xyz_h(12),19,2,rErr)

        call chr_fromReal(xyz2(2,2),xyz_h(13),19,2,rErr)
        call chr_fromReal(xyz2(3,2),xyz_h(14),19,2,rErr)
        call chr_fromReal(xyz2(4,2),xyz_h(15),19,2,rErr)
        call chr_fromReal(xyz2(5,2),xyz_h(16),19,2,rErr)
        call chr_fromReal(xyz2(6,2),xyz_h(17),19,2,rErr)

        call chr_fromReal(xyz2(3,3),xyz_h(18),19,2,rErr)
        call chr_fromReal(xyz2(4,3),xyz_h(19),19,2,rErr)
        call chr_fromReal(xyz2(5,3),xyz_h(20),19,2,rErr)
        call chr_fromReal(xyz2(6,3),xyz_h(21),19,2,rErr)

        call chr_fromReal(xyz2(4,4),xyz_h(22),19,2,rErr)
        call chr_fromReal(xyz2(5,4),xyz_h(23),19,2,rErr)
        call chr_fromReal(xyz2(6,4),xyz_h(24),19,2,rErr)

        call chr_fromReal(xyz2(5,5),xyz_h(25),19,2,rErr)
        call chr_fromReal(xyz2(6,5),xyz_h(26),19,2,rErr)

        call chr_fromReal(xyz2(6,6),xyz_h(27),19,2,rErr)

        write(unit, "(2(1x,i8),1x,f12.5,27(1x,a25))") napx,nturn,localDcum, &
          xyz_h(1), xyz_h(2), xyz_h(3), xyz_h(4), xyz_h(5), xyz_h(6), xyz_h(7), xyz_h(8), xyz_h(9), &
          xyz_h(10),xyz_h(11),xyz_h(12),xyz_h(13),xyz_h(14),xyz_h(15),xyz_h(16),xyz_h(17),xyz_h(18),&
          xyz_h(19),xyz_h(20),xyz_h(21),xyz_h(22),xyz_h(23),xyz_h(24),xyz_h(25),xyz_h(26),xyz_h(27)
      else
        call chr_fromReal(xyz(1),xyz_l(1),10,2,rErr)
        call chr_fromReal(xyz(2),xyz_l(2),10,2,rErr)
        call chr_fromReal(xyz(3),xyz_l(3),10,2,rErr)
        call chr_fromReal(xyz(4),xyz_l(4),10,2,rErr)
        call chr_fromReal(xyz(5),xyz_l(5),10,2,rErr)
        call chr_fromReal(xyz(6),xyz_l(6),10,2,rErr)

        call chr_fromReal(xyz2(1,1),xyz_l(7), 10,2,rErr)
        call chr_fromReal(xyz2(2,1),xyz_l(8), 10,2,rErr)
        call chr_fromReal(xyz2(3,1),xyz_l(9), 10,2,rErr)
        call chr_fromReal(xyz2(4,1),xyz_l(10),10,2,rErr)
        call chr_fromReal(xyz2(5,1),xyz_l(11),10,2,rErr)
        call chr_fromReal(xyz2(6,1),xyz_l(12),10,2,rErr)

        call chr_fromReal(xyz2(2,2),xyz_l(13),10,2,rErr)
        call chr_fromReal(xyz2(3,2),xyz_l(14),10,2,rErr)
        call chr_fromReal(xyz2(4,2),xyz_l(15),10,2,rErr)
        call chr_fromReal(xyz2(5,2),xyz_l(16),10,2,rErr)
        call chr_fromReal(xyz2(6,2),xyz_l(17),10,2,rErr)

        call chr_fromReal(xyz2(3,3),xyz_l(18),10,2,rErr)
        call chr_fromReal(xyz2(4,3),xyz_l(19),10,2,rErr)
        call chr_fromReal(xyz2(5,3),xyz_l(20),10,2,rErr)
        call chr_fromReal(xyz2(6,3),xyz_l(21),10,2,rErr)

        call chr_fromReal(xyz2(4,4),xyz_l(22),10,2,rErr)
        call chr_fromReal(xyz2(5,4),xyz_l(23),10,2,rErr)
        call chr_fromReal(xyz2(6,4),xyz_l(24),10,2,rErr)

        call chr_fromReal(xyz2(5,5),xyz_l(25),10,2,rErr)
        call chr_fromReal(xyz2(6,5),xyz_l(26),10,2,rErr)

        call chr_fromReal(xyz2(6,6),xyz_l(27),10,2,rErr)

        write(unit, "(2(1x,i8),1x,f12.5,27(1x,a16))") napx,nturn,localDcum, &
          xyz_l(1), xyz_l(2), xyz_l(3), xyz_l(4), xyz_l(5), xyz_l(6), xyz_l(7), xyz_l(8), xyz_l(9), &
          xyz_l(10),xyz_l(11),xyz_l(12),xyz_l(13),xyz_l(14),xyz_l(15),xyz_l(16),xyz_l(17),xyz_l(18),&
          xyz_l(19),xyz_l(20),xyz_l(21),xyz_l(22),xyz_l(23),xyz_l(24),xyz_l(25),xyz_l(26),xyz_l(27)
      end if
      flush(unit,iostat=ierro)
#ifdef CR
      dumpfilepos(dumpIdx) = dumpfilepos(dumpIdx)+1
#endif
    end if

  ! ------------------------------------------------------------------
  !  Format #101
  !  Same as fmt 3, but with additional variables
  ! ------------------------------------------------------------------
  else if(fmt == 101) then
    if(i == 0 .and. ix == 0) then
      localDcum   = zero
      localKtrack = 0
    else
      localDcum   = dcum(i)
      localKtrack = ktrack(i)
    end if
#ifdef HDF5
    if(h5_useForDUMP) then
      call h5_prepareWrite(dump_hdf5DataSet(ix), napx)
      call h5_writeData(dump_hdf5DataSet(ix), 1,  napx, partID)
      call h5_writeData(dump_hdf5DataSet(ix), 2,  napx, nturn)
      call h5_writeData(dump_hdf5DataSet(ix), 3,  napx, localDcum)
      call h5_writeData(dump_hdf5DataSet(ix), 4,  napx, xv1(:))
      call h5_writeData(dump_hdf5DataSet(ix), 5,  napx, yv1(:))
      call h5_writeData(dump_hdf5DataSet(ix), 6,  napx, xv2(:))
      call h5_writeData(dump_hdf5DataSet(ix), 7,  napx, yv2(:))
      call h5_writeData(dump_hdf5DataSet(ix), 8,  napx, sigmv)
      call h5_writeData(dump_hdf5DataSet(ix), 9,  napx, (ejv-e0)/e0)
      call h5_writeData(dump_hdf5DataSet(ix), 10, napx, localKtrack)
      call h5_writeData(dump_hdf5DataSet(ix), 11, napx, ejv(j))
      call h5_writeData(dump_hdf5DataSet(ix), 12, napx, ejfv(j))
      call h5_writeData(dump_hdf5DataSet(ix), 13, napx, dpsv(j))
      call h5_writeData(dump_hdf5DataSet(ix), 14, napx, oidpsv(j))
      call h5_writeData(dump_hdf5DataSet(ix), 15, napx, rvv(j))
      call h5_writeData(dump_hdf5DataSet(ix), 16, napx, nucm(j))
      call h5_writeData(dump_hdf5DataSet(ix), 17, napx, mtc(j))
      call h5_writeData(dump_hdf5DataSet(ix), 18, napx, e0)
      call h5_writeData(dump_hdf5DataSet(ix), 19, napx, e0f)
      call h5_finaliseWrite(dump_hdf5DataSet(ix))
    else
#endif
      do j=1,napx
        write(unit) partID(j),nturn,localDcum, &
                    xv1(j),yv1(j),xv2(j),yv2(j), &
                    sigmv(j),(ejv(j)-e0)/e0,localKtrack, &
                    ejv(j), ejfv(j), dpsv(j), oidpsv(j), &
                    rvv(j), nucm(j), mtc(j), e0, e0f
      end do
      flush(unit,iostat=ierro)
#ifdef CR
      dumpfilepos(dumpIdx) = dumpfilepos(dumpIdx)+napx
#endif
#ifdef HDF5
    end if
#endif

#ifdef FLUKA
  ! ------------------------------------------------------------------ !
  !  Format #20
  !  Same as fmt 2, but also include fluka vars.)
  ! ------------------------------------------------------------------ !
  else if(fmt == 20) then
    if(i == 0 .and. ix == 0) then
      localDcum   = zero
      localKtrack = 0
    else
      localDcum   = dcum(i)
      localKtrack = ktrack(i)
    end if
      if(lhighprec) then
        do j=1,napx
          call chr_fromReal(xv1(j),       xyz_h(1),19,2,rErr)
          call chr_fromReal(yv1(j),       xyz_h(2),19,2,rErr)
          call chr_fromReal(xv2(j),       xyz_h(3),19,2,rErr)
          call chr_fromReal(yv2(j),       xyz_h(4),19,2,rErr)
          call chr_fromReal(sigmv(j),      xyz_h(5),19,2,rErr)
          call chr_fromReal((ejv(j)-e0)/e0,xyz_h(6),19,2,rErr)
          write(unit,"(3(1x,i8),1x,2(f12.5),6(1x,a25),4(1x,i8),1x,i12)") partID(j), nturn, parentID(j), partWeight(j), &
            localDcum, xyz_h(1),xyz_h(2),xyz_h(3),xyz_h(4),xyz_h(5),xyz_h(6),localKtrack, naa(j), nzz(j), nqq(j), pdgid(j)
        end do
      else
        do j=1,napx
          call chr_fromReal(xv1(j),       xyz_l(1),10,2,rErr)
          call chr_fromReal(yv1(j),       xyz_l(2),10,2,rErr)
          call chr_fromReal(xv2(j),       xyz_l(3),10,2,rErr)
          call chr_fromReal(yv2(j),       xyz_l(4),10,2,rErr)
          call chr_fromReal(sigmv(j),      xyz_l(5),10,2,rErr)
          call chr_fromReal((ejv(j)-e0)/e0,xyz_l(6),10,2,rErr)
          write(unit,"(3(1x,i8),1x,2(f12.5),6(1x,a16),4(1x,i8),1x,i12)") partID(j), nturn, parentID(j), partWeight(j), &
            localDcum, xyz_l(1),xyz_l(2),xyz_l(3),xyz_l(4),xyz_l(5),xyz_l(6),localKtrack, naa(j), nzz(j), nqq(j), pdgid(j)
        end do
      end if

      ! Flush
      flush(unit,iostat=ierro)
#ifdef CR
      dumpfilepos(dumpIdx) = dumpfilepos(dumpIdx)+napx
#endif
#endif

  ! ------------------------------------------------------------------ !
  !  Format #200
  !  Same as fmt 2, but also include ion variables
  ! ------------------------------------------------------------------ !
  else if(fmt == 200) then
    if(i == 0 .and. ix == 0) then
      localDcum   = zero
      localKtrack = 0
    else
      localDcum   = dcum(i)
      localKtrack = ktrack(i)
    end if
#ifdef HDF5
    if(h5_useForDUMP) then
      call h5_prepareWrite(dump_hdf5DataSet(ix), napx)
      call h5_writeData(dump_hdf5DataSet(ix), 1,  napx, partID)
      call h5_writeData(dump_hdf5DataSet(ix), 2,  napx, nturn)
      call h5_writeData(dump_hdf5DataSet(ix), 3,  napx, localDcum)
      call h5_writeData(dump_hdf5DataSet(ix), 4,  napx, xv1(:))
      call h5_writeData(dump_hdf5DataSet(ix), 5,  napx, yv1(:))
      call h5_writeData(dump_hdf5DataSet(ix), 6,  napx, xv2(:))
      call h5_writeData(dump_hdf5DataSet(ix), 7,  napx, yv2(:))
      call h5_writeData(dump_hdf5DataSet(ix), 8,  napx, sigmv)
      call h5_writeData(dump_hdf5DataSet(ix), 9,  napx, (ejv-e0)/e0)
      call h5_writeData(dump_hdf5DataSet(ix), 10, napx, localKtrack)
      call h5_writeData(dump_hdf5DataSet(ix), 11, napx, int(nzz))
      call h5_writeData(dump_hdf5DataSet(ix), 12, napx, int(naa))
      call h5_writeData(dump_hdf5DataSet(ix), 13, napx, int(nqq))
      call h5_writeData(dump_hdf5DataSet(ix), 14, napx, pdgid)
      call h5_finaliseWrite(dump_hdf5DataSet(ix))
    else
#endif
      if(lhighprec) then
        do j=1,napx
          call chr_fromReal(xv1(j),       xyz_h(1),19,2,rErr)
          call chr_fromReal(yv1(j),       xyz_h(2),19,2,rErr)
          call chr_fromReal(xv2(j),       xyz_h(3),19,2,rErr)
          call chr_fromReal(yv2(j),       xyz_h(4),19,2,rErr)
          call chr_fromReal(sigmv(j),      xyz_h(5),19,2,rErr)
          call chr_fromReal((ejv(j)-e0)/e0,xyz_h(6),19,2,rErr)
          write(unit,"(2(1x,i8),1x,f12.5,6(1x,a25),1x,i8,3(1x,i8),1x,i12)") partID(j),nturn,localDcum,&
            xyz_h(1),xyz_h(2),xyz_h(3),xyz_h(4),xyz_h(5),xyz_h(6),localKtrack, &
            naa(j),nzz(j),nqq(j),pdgid(j)
        end do
      else
        do j=1,napx
          call chr_fromReal(xv1(j),       xyz_l(1),10,2,rErr)
          call chr_fromReal(yv1(j),       xyz_l(2),10,2,rErr)
          call chr_fromReal(xv2(j),       xyz_l(3),10,2,rErr)
          call chr_fromReal(yv2(j),       xyz_l(4),10,2,rErr)
          call chr_fromReal(sigmv(j),      xyz_l(5),10,2,rErr)
          call chr_fromReal((ejv(j)-e0)/e0,xyz_l(6),10,2,rErr)
          write(unit,"(2(1x,i8),1x,f12.5,6(1x,a16),1x,i8,3(1x,i8),1x,i12)") partID(j),nturn,localDcum,&
            xyz_l(1),xyz_l(2),xyz_l(3),xyz_l(4),xyz_l(5),xyz_l(6),localKtrack, &
            naa(j),nzz(j),nqq(j),pdgid(j)
        end do
      end if
      flush(unit,iostat=ierro)
#ifdef CR
      dumpfilepos(dumpIdx) = dumpfilepos(dumpIdx)+napx
#endif
#ifdef HDF5
    end if
#endif

#ifdef ROOT
  !output directly to a root file
  else if(fmt == rootfmt) then
    if(i == 0 .and. ix == 0) then
      localDcum = 0.0
      localBez  = "StartDUMP"
    else
      localDcum = dcum(i)
      if (ktrack(i) /= 1) then
        localBez = bez(ix)
      else                ! BLOCs
        localBez = bezb(ic(i))
      end if
    end if
    tmpWeight = 1.0
    do j=1,napx
#if defined(FLUKA) || defined(G4COLLIMATION)
      tmpWeight = partWeight(j)
#endif
      call root_DumpBunch(localBez, len_trim(localBez), i, ix, nturn, partID(j), parentID(j), pdgid(j), nqq(j), tmpWeight, &
&     localDcum, xv1(j), yv1(j), xv2(j), yv2(j), sigmv(j), (ejv(j)-e0)/e0, spin_x(j), spin_y(j), spin_z(j), nucm(j))
    end do
#endif

  ! Unrecognized format fmt
  ! ------------------------------------------------------------------ !
  else
    write(lerr,"(a,i0,a)") "DUMP> ERROR Format ",fmt," not understood for file '"//trim(dump_fname(ix))//"'"
    call prror
  end if

  call time_stopClock(time_clockDUMP)

  return

end subroutine dump_beam_population

#ifdef CR
! ================================================================================================================================ !
!  Begin Checkpoint Restart
! ================================================================================================================================ !
subroutine dump_crcheck_readdata(fileunit, readerr)

  use parpro
  use crcoall

  integer, intent(in) :: fileunit
  logical, intent(out) :: readerr

  integer j

  read(fileunit,err=100,end=100) (dumpfilepos_cr(j),j=-1,nele)

  readerr = .false.
  return

100 continue
  readerr = .true.
  write(lout, "(a,i0,a)") "CR_CHECK> ERROR Reading C/R file fort.",fileUnit," in DUMP"
  write(crlog,"(a,i0,a)") "CR_CHECK> ERROR Reading C/R file fort.",fileUnit," in DUMP"
  flush(crlog)

end subroutine dump_crcheck_readdata

! ================================================================================================================================ !
subroutine dump_crcheck_positionFiles

  use crcoall
  use mod_units
  use mod_common
  use string_tools

  ! For skipping through binary DUMP files (format 3&8)
  integer tmp_ID, tmp_nturn, tmp_ktrack
  real(kind=fPrec) tmp_dcum, tmp_x, tmp_xp, tmp_y, tmp_yp, tmp_sigma, tmp_dEE

  integer i,j
  logical lerror,lopen
  character(len=1024) arecord

  do i=-1, il
    if(ldump(i)) then
      inquire(unit=dumpunit(i), opened=lopen)
      if(dumpfmt(i) == 3 .or. dumpfmt(i) == 8 .or. dumpfmt(i) == 101) then ! Binary
        if(lopen .eqv. .false.) then
          call f_open(unit=dumpunit(i),file=trim(dump_fname(i)),formatted=.false.,mode="rw-",status="old")
        end if
        dumpfilepos(i) = 0
        do j=1,dumpfilepos_cr(i)
          if(dumpfmt(i) == 3 .or. dumpfmt(i) == 8) then
            read(dumpunit(i),end=111,err=111,iostat=ierro) &
              tmp_ID,tmp_nturn,tmp_dcum,tmp_x,tmp_xp,tmp_y,tmp_yp,tmp_sigma,tmp_dEE,tmp_ktrack
          elseif(dumpfmt(i) == 101) then
            read(dumpunit(i),end=111,err=111,iostat=ierro) &
              tmp_ID,tmp_nturn,tmp_dcum,tmp_x,tmp_xp,tmp_y,tmp_yp,tmp_sigma,tmp_dEE,tmp_ktrack, &
              tmp_x, tmp_x, tmp_x, tmp_x, tmp_x, tmp_x, tmp_x, tmp_x
          else
            write(lout, "(a,i0)") "CR_CHECK> ERROR Failed positioning DUMP file: unknown format ",dumpfmt(i)
            write(crlog,"(a,i0)") "CR_CHECK> ERROR Failed positioning DUMP file: unknown format ",dumpfmt(i)
            flush(crlog)
            call prror
          end if
          dumpfilepos(i) = dumpfilepos(i) + 1
        end do
      else ! Text
        if(lopen .eqv. .false.) then
          call f_open(unit=dumpunit(i),file=trim(dump_fname(i)),formatted=.true.,mode="rw-",status="old")
        end if

        dumpfilepos(i) = 0
        do j=1,dumpfilepos_cr(i)
          read(dumpunit(i),'(a1024)',end=111,err=111,iostat=ierro) arecord
          dumpfilepos(i) = dumpfilepos(i) + 1
        end do
      end if
      write(crlog,"(3(a,i0))") "CR_CHECK> Repositioned DUMP file '"//trim(dump_fname(i))//&
        "' Position: ",dumpfilepos(i),", C/R: ",dumpfilepos_cr(i)
      flush(crlog)
    end if
  end do

  ! Crop DUMP files (if used by multiple DUMPs,
  ! the actual position is the sum of the dumpfileposes
  do i=0,il
    if(ldump(i)) then
      ! This is not a FLUSH!
      endfile(dumpunit(i),iostat=ierro)

      ! Change from 'readwrite' to 'write'
      call f_close(dumpunit(i))
      if(dumpfmt(i) == 3 .or. dumpfmt(i) == 8 .or. dumpfmt(i) == 101) then ! Binary
        call f_open(unit=dumpunit(i),file=trim(dump_fname(i)),formatted=.false.,mode="w+",status="old")
      else ! Text
        call f_open(unit=dumpunit(i),file=trim(dump_fname(i)),formatted=.true.,mode="w+",status="old")
      end if
    end if
  end do

  return

111 continue
  write(crlog,"(2(a,i0))") "CR_CHECK> ERROR Failed positioning DUMP file on unit ",dumpunit(i),", iostat: ",ierro
  write(crlog,"(2(a,i0))") "CR_CHECK>       File position: ",dumpfilepos(i),", C/R position: ",dumpfilepos_cr(i)
  flush(crlog)
  write(lerr,"(a,i0)") "CR_CHECK> ERROR Failed positioning DUMP file on unit ",dumpunit(i)
  call prror

end subroutine dump_crcheck_positionFiles

! ================================================================================================================================ !
subroutine dump_crpoint(fileunit,lerror)

  use parpro
  use crcoall

  integer, intent(in)  :: fileunit
  logical, intent(out) :: lerror

  integer j

  write(fileunit,err=100) (dumpfilepos(j),j=-1,nele)
  flush(fileunit)

  lerror = .false.

  return

100 continue
  lerror = .true.
  write(lout, "(a,i0,a)") "CR_POINT> ERROR Writing C/R file fort.",fileUnit," in DUMP"
  write(crlog,"(a,i0,a)") "CR_POINT> ERROR Writing C/R file fort.",fileUnit," in DUMP"
  flush(crlog)

end subroutine dump_crpoint
! ================================================================================================================================ !
!  End Checkpoint Restart
! ================================================================================================================================ !
#endif
end module dump
