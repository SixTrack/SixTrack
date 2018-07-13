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
  use parpro ! For nele
  use mod_alloc
#ifdef HDF5
  use hdf5_output
#endif

  implicit none

  ! High precision printout required at all flagged SINGLE ELEMENTs
  logical, save :: ldumphighprec
  ! Dump at the beginning of each element, not at the end.
  logical, save :: ldumpfront
  ! flag the SINGLE ELEMENT for dumping
  logical, allocatable, save :: ldump(:)  !(-1:nele)

  ! Dump every n turns at a flagged SINGLE ELEMENT (dump frequency)
  integer, allocatable, save :: ndumpt(:) !(-1:nele)
  ! First turn for DUMP to be active
  integer, allocatable, save :: dumpfirst(:) !(-1:nele)
  ! Last turn for this DUMP to be active (-1=all)
  integer, allocatable, save :: dumplast(:) !(-1:nele)
  ! Fortran unit for dump at a flagged SINGLE ELEMENT
  integer, allocatable, save :: dumpunit(:) !(-1:nele)
  ! Flag the format of the dump
  integer, allocatable, save :: dumpfmt(:) !(-1:nele)
  ! Filename to write the dump to
  character(len=:), allocatable, save :: dump_fname(:) !(mStrLen)(-1:nele)

  ! tas matrix used for nomalisation of phase space in DUMP and FMA.
  ! First index = -1 -> StartDUMP, filled differently than idx > 0; First index = 0  -> Unused.
  real(kind=fPrec), allocatable, save :: dumptas(:,:,:)    ! (-1:nblz,6,6)
  ! inverse matrix of dumptas
  real(kind=fPrec), allocatable, save :: dumptasinv(:,:,:) ! (-1:nblz,6,6)
  ! closed orbit used for normalisation of phase space
  ! TODO: check units used in dumpclo; is x' or px used?
  real(kind=fPrec), allocatable, save :: dumpclo(:,:)      ! (-1:nblz,6)

#ifdef HDF5
  ! Array to save hdf5 formats for each dump format
  integer, allocatable, save :: dump_hdf5Format(:)
  integer, allocatable, save :: dump_hdf5DataSet(:)
#endif

#ifdef CR
  ! For resetting file positions
  integer, allocatable, save :: dumpfilepos(:), dumpfilepos_cr(:) !(-1:nele)
#endif

! ================================================================================================================================ !
!  THE SUBROUTINES
! ================================================================================================================================ !
contains

! ================================================================================================================================ !
subroutine dump_expand_arrays(nele_new, nblz_new)

  use numerical_constants, only : zero

  implicit none

  integer, intent(in) :: nele_new
  integer, intent(in) :: nblz_new

  call alloc(ldump,               nele_new, .false.,    "ldump",      -1)
  call alloc(ndumpt,              nele_new, 0,          "ndumpt",     -1)
  call alloc(dumpfirst,           nele_new, 0,          "dumpfirst",  -1)
  call alloc(dumplast,            nele_new, 0,          "dumplast",   -1)
  call alloc(dumpunit,            nele_new, 0,          "dumpunit",   -1)
  call alloc(dumpfmt,             nele_new, 0,          "dumpfmt",    -1)
  call alloc(dump_fname, mStrLen, nele_new, str_dSpace, "dump_fname", -1)

  call alloc(dumptas,             nblz_new, 6, 6, zero, "dumptas",    -1,1,1)
  call alloc(dumptasinv,          nblz_new, 6, 6, zero, "dumptasinv", -1,1,1)
  call alloc(dumpclo,             nblz_new, 6,    zero, "dumpclo",    -1,1)

#ifdef CR
  call alloc(dumpfilepos,         nele_new,-1,          "dumpfilepos",   -1)
  call alloc(dumpfilepos_cr,      nele_new,-1,          "dumpfilepos_cr",-1)
#endif

#ifdef HDF5
  call alloc(dump_hdf5DataSet,    nele_new,0,           "dump_hdf5DataSet",-1)
  call alloc(dump_hdf5Format,     9,       0,           "dump_hdf5Format")
#endif

end subroutine dump_expand_arrays

subroutine dump_lines(n,i,ix)

  use mod_commont

  implicit none

  integer, intent(in) :: n,i,ix

  if ( ldump(0) ) then
    ! Dump at all SINGLE ELEMENTs
    if (ndumpt(0).eq.1 .or. mod(n,ndumpt(0)).eq.1) then
      if ((n.ge.dumpfirst(0)) .and. ((n.le.dumplast(0)) .or. (dumplast(0).eq.-1))) then
        call dump_beam_population(n, i, ix, dumpunit(0), dumpfmt(0), ldumphighprec, dumpclo(ix,1:6),dumptasinv(ix,1:6,1:6))
      end if
    end if
  end if
  if (ktrack(i) .ne. 1) then
    ! The next "if" is only safe for SINGLE ELEMENTS, not BLOC where ix<0.
    if (ldump(ix)) then
      ! Dump at this precise SINGLE ELEMENT
      if (ndumpt(ix).eq.1 .or. mod(n,ndumpt(ix)).eq.1) then
        if ((n.ge.dumpfirst(ix)) .and. ((n.le.dumplast(ix)) .or. (dumplast(ix).eq.-1))) then
          call dump_beam_population(n, i, ix, dumpunit(ix), dumpfmt(ix), ldumphighprec, dumpclo(ix,1:6),dumptasinv(ix,1:6,1:6))
        end if
      end if
    end if
  end if

end subroutine dump_lines

! ================================================================================================================================ !
subroutine dump_linesFirst(n)

  implicit none
  integer, intent(in) :: n

  ! StartDUMP - dump on the first element
  if (ldump(-1)) then
    if (ndumpt(-1).eq.1 .or. mod(n,ndumpt(-1)).eq.1) then
      if ((n.ge.dumpfirst(-1)) .and. ((n.le.dumplast(-1)) .or. (dumplast(-1).eq.-1))) then
        call dump_beam_population(n, 0, 0, dumpunit(-1), dumpfmt(-1), ldumphighprec, dumpclo(-1,1:6),dumptasinv(-1,1:6,1:6))
      end if
    end if
  end if

end subroutine dump_linesFirst

! ================================================================================================================================ !
!  A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!  Last modified: 01-09-2014
!  Close units for dumping particle population
! ================================================================================================================================ !
subroutine dump_closeUnits

  use mod_common
  implicit none
  logical lopen
  integer i

#ifdef HDF5
  if(.not. h5_useForDUMP) then
#endif
    do i=0,il
      if (ldump(i)) then
        ! The same file could be used by more than one SINGLE ELEMENT
        inquire( unit=dumpunit(i), opened=lopen )
        if (lopen) close(dumpunit(i))
      end if
    end do
#ifdef HDF5
  end if
#endif

end subroutine dump_closeUnits

! ================================================================================================================================ !
subroutine dump_parseInputLine(ch,iErr)

  use crcoall
  use mod_common
  use file_units
  use string_tools

  implicit none

  character(len=*), intent(in)    :: ch
  logical,          intent(inout) :: iErr

  ! Fields split variables
  character gFields(str_maxFields)*(mStrLen) ! Array of fields
  integer   nFields                             ! Number of identified fields
  integer   lFields(str_maxFields)              ! Length of each what:
  logical   fieldsErr                           ! An error flag

  ! Temp variables
  character(len=mNameLen) idat ! Synchronized with daten
  integer i1,i2,i3,i4,i5,kk,j

  character(len=mStrLen) ch1

  ! initialise reading variables, to avoid storing nonsense values
  idat = ' ' ! Name
  i1 = 0     ! frequency
  i2 = -1    ! unit
  i3 = 0     ! format
  i4 = 1     ! first turn
  i5 = -1    ! last turn

  if(ch(:4).eq.'HIGH') then
    ldumphighprec = .true.
    return
  else if(ch(:5).eq.'FRONT') then
    ldumpfront = .true.
    return
  end if

  ! Requested element
  call getfields_split(ch, gFields, lFields, nFields, fieldsErr)
  if (fieldsErr) call prror(-1)

  if ((nFields .lt. 4) .or. (nFields .gt. 7) .or. (nFields .eq. 6)) then
    write(lout,"(a,i0)")    "DUMP> ERROR Expected 4 to 7 (but not 6) arguments, got ",nFields
    write(lout,"(a)")      ("DUMP>       '"//gFields(kk)(1:lFields(kk))//"' ",kk=1,nFields)
    iErr = .true.
    return
  end if

  if (lFields(1) > mNameLen) then
    write(lout,"(a,i0,a)") "DUMP> ERROR Element names are max. ",mNameLen," characters"
    iErr = .true.
    return
  end if

  idat = gFields(1)(1:lFields(1))
  read(gFields(2)(1:lFields(2)),*) i1
  read(gFields(3)(1:lFields(3)),*) i2
  read(gFields(4)(1:lFields(4)),*) i3
  if (nFields .eq. 4) then
    ! Automatic fname
    write(ch1,"(a5,I0)") "fort.", i2
  else if ( (nFields .eq. 5) .or. (nFields .eq. 7)) then
    ! Given fname
    ch1 = gFields(5)(1:lFields(5))
  else
    iErr = .true.
    return
  end if
  if (nFields .eq. 7) then
    read(gFields(6)(1:lFields(6)),*) i4
    read(gFields(7)(1:lFields(7)),*) i5
  end if

  ! Check that first/last turn is sane
  if (i5.ne.-1) then
    if (i5 .lt. i4) then
      write(lout,"(2(a,i0))") "DUMP> ERROR Expect last turn >= first turn, unless last turn = -1 (infinity), got ", i4,", ",i5
      iErr = .true.
      return
      end if
  end if
  if (i4 .lt. 1) then
    write(lout,"(a,i0)") "DUMP> ERROR Expect first turn >= 1, got", i4
    iErr = .true.
    return
  end if

  ! Find it in the list of SINGLE ELEMENTs:
  do j=1,il
    if(bez(j).eq.idat) then
      if (ldump(j)) then ! Only enable once/element!
        write(lout,"(a)") "DUMP> ERROR Element '"//trim(idat)//"' was specified more than once"
        iErr = .true.
        return
      end if

      ! Element was found in SINGLE ELEMENTS list, now do some sanity checks
      if(trim(bez(j)).eq."ALL") then
        write(lout,"(a)") "DUMP> ERROR The element name 'ALL' cannot be used in the SINGLE ELEMENTS list "//&
          "when an 'ALL' special DUMP is active."
        iErr = .true.
        return
      else if(trim(bez(j)).eq."StartDUMP") then
        write(lout,"(a)") "DUMP> ERROR The element name 'StartDUMP' cannot be used in the SINGLE ELEMENTS "// &
          "list when an 'StartDUMP' special DUMP is active."
        iErr = .true.
        return
      end if
      goto 10 ! Element found, store the data
    end if
  end do

  if (idat(:3).eq.'ALL') then
    j=0
    if (ldump(j)) then
      write(lout,"(a)") "DUMP> ERROR Element 'ALL' was specified (at least) twice"
      iErr = .true.
      return
      end if
    goto 10 ! Element found, store the data
  end if
  if (idat(:9).eq.'StartDUMP') then
    j=-1
    if (ldump(j)) then
      write(lout,"(a)") "DUMP> ERROR Element 'StartDUMP' was specified (at least) twice"
      iErr = .true.
      return
      end if
    goto 10 ! Element found, store the data
  end if

  ! Search failed, fall-through to here:
  write(lout,"(a)") "DUMP> ERROR Unidentified SINGLE ELEMENT '"//idat//"'"
  iErr = .true.
  return

  ! Element found, store the data:
10 continue
  ldump(j)     = .true.
  ndumpt(j)    = i1
  dumpunit(j)  = i2
  dumpfmt(j)   = i3
  dump_fname(j)(1:lFields(5)) = ch1(1:lFields(5))
  dumpfirst(j) = i4
  dumplast(j)  = i5
  if(ndumpt(j) <= 0) ndumpt(j) = 1
#ifdef HDF5
  if(h5_useForDUMP .eqv. .false.) then
#endif
    if(dumpunit(j) == -1) then
      call funit_requestUnit(chr_trimZero(dump_fname(j)),dumpunit(j))
    end if
#ifdef HDF5
  end if
#endif

  return

end subroutine dump_parseInputLine

! ================================================================================================================================ !
subroutine dump_parseInputDone(iErr)

  use crcoall
  use mod_common
  use string_tools

  implicit none

  logical, intent(inout) :: iErr

  ! Temp variables
  integer ii,jj,kk
  character(len=mStrLen) ch1

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
        write (lout,"(a)") "DUMP> ERROR No structure elements found for '"//trim(bez(ii))//"'"
        write (lout,"(a)") "DUMP>       This element is probably only found in a BLOC, or it is not used at all."
        iErr = .true.
        return
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

end subroutine dump_parseInputDone

! ================================================================================================================================ !
!  A.Mereghetti, D.Sinuela-Pastor & P.Garcia Ortega, for the FLUKA Team
!  K.Sjobak, A.Santamaria, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-05-08
! ================================================================================================================================ !
subroutine dump_initialise

  use numerical_constants, only : zero
  use crcoall
  use string_tools
  use mod_common

  implicit none

  integer i,j
  logical lopen
  character(len=256) filename

#ifdef HDF5
  type(h5_dataField), allocatable :: setFields(:)
  if(h5_useForDUMP) then
    call h5_initForDump()
    goto 3100 ! Skip the normal file initialisation, and do the HDF5 version instead
  end if
#endif

  do i=-1,il
#ifdef CR
    if (dumpfilepos(i).ge.0) then
      ! Expect the file to be opened already, in crcheck
      inquire( unit=dumpunit(i), opened=lopen )
      if (.not.lopen) then
        write(lout,"(2(a,i0),a)") "DUMP> ERROR The unit ",dumpunit(i)," has dumpfilepos = ", dumpfilepos(i), " >= 0, "//&
          "but the file is NOT open. This is probably a bug."
        call prror(-1)
      end if
      cycle ! Everything OK, don't try to open the files again.
    end if
#endif
    if (ldump(i)) then
      ! The same file could be used by more than one SINGLE ELEMENT
      inquire( unit=dumpunit(i), opened=lopen )
      if (.not.lopen) then
        ! Check that the filename is not already taken
        do j=-1,i-1
          if (ldump(j) .and. (dump_fname(j).eq.dump_fname(i))) then
            write(lout,"(2(a,i0))") "DUMP> ERROR Output filename '"//trim(dump_fname(i))//&
              "' is used by two DUMPS, but output units differ: ",dumpunit(i)," vs ",dumpunit(j)
            call prror(-1)
          end if
        end do
        if (dumpfmt(i).eq.3 .or. dumpfmt(i).eq.8) then ! Binary dump
#ifdef BOINC
          call boincrf(dump_fname(i),filename)
          open(dumpunit(i),file=filename,status='replace',form='unformatted')
#else
          open(dumpunit(i),file=trim(chr_trimZero(dump_fname(i))),status='replace',form='unformatted')
#endif
        else ! ASCII dump
#ifdef BOINC
          call boincrf(dump_fname(i),filename)
          open(dumpunit(i),file=filename,status='replace',form='formatted')
#else
          open(dumpunit(i),file=trim(chr_trimZero(dump_fname(i))),status='replace',form='formatted')
#endif
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
            if (dumpunit(j).eq.dumpunit(i)) then
              if (dumpfmt(j).ne.dumpfmt(i)) then
                write(lout,"(a,i0,a)") "DUMP> ERROR Output unit ",dumpunit(i)," used by two DUMPS, formats are not the same."
                call prror(-1)
              else if (j.eq.0) then
                write(lout,"(a,i0,a)") "DUMP> ERROR Output unit ",dumpunit(i)," used by two DUMPS, one of which is ALL"
                call prror(-1)
              else if (j.eq.-1) then
                write(lout,"(a,i0,a)") "DUMP> ERROR Output unit ",dumpunit(i)," used by two DUMPS, one of which is StartDUMP"
                call prror(-1)
              else if (dump_fname(j).ne.dump_fname(i)) then
                write(lout,"(a,i0,a)") "DUMP> ERROR Output unit ",dumpunit(i)," used by two DUMPS, but filenames differ: '"//&
                  trim(dump_fname(i)),"' vs '",trim(dump_fname(j)),"'"
                call prror(-1)
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
          write(lout,"(a,i0,a)") "DUMP> ERROR Unit",dumpunit(i)," is already open, but not by DUMP. Please pick another unit!"
          write(lout,"(a)")      "DUMP> Note: This check is not watertight as other parts of the program may later open the "
          write(lout,"(a)")      "DUMP>       same unit. Althernatively, the unit can be specified as -1 and a unit is assigned."
          call prror(-1)
        end if
      end if

      ! Write format-specific headers
      if (dumpfmt(i).eq.1) then
        write(dumpunit(i),'(a)') '# ID turn s[m] x[mm] xp[mrad] y[mm] yp[mrad] dE/E[1] ktrack'
        ! Flush file
        endfile   (dumpunit(i))
        backspace (dumpunit(i))
#ifdef CR
        dumpfilepos(i) = dumpfilepos(i) + 1
#endif
      else if (dumpfmt(i).eq.2 .or. dumpfmt(i).eq.4 .or. dumpfmt(i).eq.5 .or. &
               dumpfmt(i).eq.6 .or. dumpfmt(i).eq.7 .or. dumpfmt(i).eq.9 ) then
        ! Write the general header
        if (i.eq.-1) then  ! STARTdump
          write(dumpunit(i),'(a,i0,a,a16,4(a,i12),2(a,L1))') &
            '# DUMP format #',dumpfmt(i),', START=',bez(1), ', number of particles=',napx, ', dump period=',ndumpt(i), &
            ', first turn=', dumpfirst(i), ', last turn=',dumplast(i), ', HIGH=',ldumphighprec, ', FRONT=',ldumpfront
        else if (i.eq.0) then ! ALL
          write(dumpunit(i),'(a,i0,a,4(a,i12),2(a,L1))') &
            '# DUMP format #',dumpfmt(i),', ALL ELEMENTS,', ' number of particles=',napx, ', dump period=',ndumpt(i), &
            ', first turn=', dumpfirst(i), ', last turn=',dumplast(i), ', HIGH=',ldumphighprec, ', FRONT=',ldumpfront
        else ! Normal element
          write(dumpunit(i),'(a,i0,a,a16,4(a,i12),2(a,L1))') &
            '# DUMP format #',dumpfmt(i), ', bez=', bez(i), ', number of particles=',napx,', dump period=',ndumpt(i), &
            ', first turn=',dumpfirst(i), ', last turn=',dumplast(i), ', HIGH=',ldumphighprec, ', FRONT=',ldumpfront
        end if

        ! Write the format-specific headers:
        if (dumpfmt(i).eq.2) then ! FORMAT 2
          write(dumpunit(i),'(a,a)') '# ID turn s[m] x[mm] xp[mrad] y[mm] yp[mrad] z[mm] dE/E[1] ktrack'
        else if (dumpfmt(i).eq.4) then ! FORMAT 4
          write(dumpunit(i),'(a)') '# napx turn s[m] <x>[mm] <xp>[mrad] <y>[mm] <yp>[mrad] <z>[mm] <dE/E>[1]'
        else if (dumpfmt(i).eq.5) then ! FORMAT 5
          write(dumpunit(i),'(a)') '# napx turn s[m] ' //                  &
            '<x>[mm] <xp>[mrad] <y>[mm] <yp>[mrad] <z>[mm] <dE/E>[1] '//   &
            '<x^2> <x*xp> <x*y> <x*yp> <x*z> <x*(dE/E)> '//                &
            '<xp^2> <xp*y> <xp*yp> <xp*z> <xp*(dE/E)> '//                  &
            '<y^2> <y*yp> <y*z> <y*(dE/E)> '//                             &
            '<yp^2> <yp*z> <yp*(dE/E)> '//                                 &
            '<z^2> <z*(dE/E)> '//                                          &
            '<(dE/E)^2>'
        else if (dumpfmt(i).eq.6) then ! FORMAT 6
          write(dumpunit(i),'(a)') '# napx turn s[m] ' //                  &
            '<x>[m] <px>[1] <y>[m] <py>[1] <sigma>[m] <psigma>[1] '//      &
            '<x^2> <x*px> <x*y> <x*py> <x*sigma> <x*psigma> '//            &
            '<px^2> <px*y> <px*py> <px*sigma> <px*psigma> '//              &
            '<y^2> <y*py> <y*sigma> <y*psigma> '//                         &
            '<py^2> <py*sigma> <py*psigma> '//                             &
            '<sigma^2> <sigma*psigma> '//                                  &
            '<psigma^2>'
        else if (dumpfmt(i).eq.7 .or. dumpfmt(i).eq.9) then
          ! Normalized ASCII dump -> extra headers with matrices and closed orbit
          if (dumpfmt(i).eq.7) then ! FORMAT 7
            write(dumpunit(i),'(a)') '# ID turn s[m] nx[1.e-3 sqrt(m)] npx[1.e-3 sqrt(m)] '// &
              'ny[1.e-3 sqrt(m)] npy[1.e-3 sqrt(m)] nsig[1.e-3 sqrt(m)] ndp/p[1.e-3 sqrt(m)] ktrack'
          end if
          if (dumpfmt(i).eq.9) then ! FORMAT 9
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
          write(dumpunit(i),'(a,1x,6(1X,1PE16.9))') '# closed orbit [mm,mrad,mm,mrad,1]', &
            dumpclo(i,1),dumpclo(i,2),dumpclo(i,3),dumpclo(i,4),dumpclo(i,5),dumpclo(i,6)
          write(dumpunit(i),'(a,1x,36(1X,1PE16.9))') '# tamatrix [mm,mrad,mm,mrad,1]', &
            dumptas(i,1,1),dumptas(i,1,2),dumptas(i,1,3),dumptas(i,1,4),dumptas(i,1,5),dumptas(i,1,6), &
            dumptas(i,2,1),dumptas(i,2,2),dumptas(i,2,3),dumptas(i,2,4),dumptas(i,2,5),dumptas(i,2,6), &
            dumptas(i,3,1),dumptas(i,3,2),dumptas(i,3,3),dumptas(i,3,4),dumptas(i,3,5),dumptas(i,3,6), &
            dumptas(i,4,1),dumptas(i,4,2),dumptas(i,4,3),dumptas(i,4,4),dumptas(i,4,5),dumptas(i,4,6), &
            dumptas(i,5,1),dumptas(i,5,2),dumptas(i,5,3),dumptas(i,5,4),dumptas(i,5,5),dumptas(i,5,6), &
            dumptas(i,6,1),dumptas(i,6,2),dumptas(i,6,3),dumptas(i,6,4),dumptas(i,6,5),dumptas(i,6,6)
          write(dumpunit(i),'(a,1x,36(1X,1PE16.9))') '# inv(tamatrix)', &
            dumptasinv(i,1,1),dumptasinv(i,1,2),dumptasinv(i,1,3),dumptasinv(i,1,4),dumptasinv(i,1,5),dumptasinv(i,1,6), &
            dumptasinv(i,2,1),dumptasinv(i,2,2),dumptasinv(i,2,3),dumptasinv(i,2,4),dumptasinv(i,2,5),dumptasinv(i,2,6), &
            dumptasinv(i,3,1),dumptasinv(i,3,2),dumptasinv(i,3,3),dumptasinv(i,3,4),dumptasinv(i,3,5),dumptasinv(i,3,6), &
            dumptasinv(i,4,1),dumptasinv(i,4,2),dumptasinv(i,4,3),dumptasinv(i,4,4),dumptasinv(i,4,5),dumptasinv(i,4,6), &
            dumptasinv(i,5,1),dumptasinv(i,5,2),dumptasinv(i,5,3),dumptasinv(i,5,4),dumptasinv(i,5,5),dumptasinv(i,5,6), &
            dumptasinv(i,6,1),dumptasinv(i,6,2),dumptasinv(i,6,3),dumptasinv(i,6,4),dumptasinv(i,6,5),dumptasinv(i,6,6)
        end if ! Format-specific headers

        ! Flush file
        endfile   (dumpunit(i))
        backspace (dumpunit(i))
#ifdef CR
        dumpfilepos(i) = dumpfilepos(i) + 2
        ! format 7 also writes clo, tas and tasinv
        if (dumpfmt(i).eq.7 .or. dumpfmt(i).eq.9) then
          dumpfilepos(i) = dumpfilepos(i) + 3
        end if
#endif
      end if ! END if format 2/4/5/6/7/9 -> General header

      ! Normalized DUMP
      if (dumpfmt(i).eq.7 .or. dumpfmt(i).eq.8 .or. dumpfmt(i).eq.9) then
        ! Have a matrix that's not zero (i.e. did we put a 6d LINE block?)
        if (dumptas(i,1,1).eq.zero .and. dumptas(i,1,2).eq.zero .and. &
            dumptas(i,1,3).eq.zero .and. dumptas(i,1,4).eq.zero) then
          write(lout,"(a)") "DUMP> ERROR The normalization matrix appears to not be set. Did you forget to put a 6D LINE block?"
          call prror(-1)
        end if
        if(idp.eq.0 .or. ition.eq.0) then ! We're in the 4D case
          if(imc.ne.1) then ! Energy scan
            write(lout,"(a)") "DUMP> ERROR in normalized DUMP: Energy scan (imc != 1) not supported!"
            call prror(-1)
          end if
          if(i.ne.-1) then ! Not at StartDUMP
            write(lout,"(a)") "DUMP> ERROR in normalized DUMP: 4D only supported for StartDUMP!"
            call prror(-1)
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
        ! # ID turn s[m] x[mm] xp[mrad] y[mm] yp[mrad] dE/E[1] ktrack
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

      case(2)
        ! Format 2:
        ! # ID turn s[m] x[mm] xp[mrad] y[mm] yp[mrad] z[mm] dE/E[1] ktrack
        if(dump_hdf5Format(2) == 0) then
          allocate(setFields(10))
          setFields(1)  = h5_dataField(name="ID",     type=h5_typeInt)
          setFields(2)  = h5_dataField(name="TURN",   type=h5_typeInt)
          setFields(3)  = h5_dataField(name="S",      type=h5_typeReal)
          setFields(4)  = h5_dataField(name="X",      type=h5_typeReal)
          setFields(5)  = h5_dataField(name="XP",     type=h5_typeReal)
          setFields(6)  = h5_dataField(name="Y",      type=h5_typeReal)
          setFields(7)  = h5_dataField(name="YP",     type=h5_typeReal)
          setFields(8)  = h5_dataField(name="dE/E",   type=h5_typeReal)
          setFields(9)  = h5_dataField(name="Z",      type=h5_typeReal)
          setFields(10) = h5_dataField(name="KTRACK", type=h5_typeInt)
          call h5_createFormat("dumpFormat2", setFields, dump_hdf5Format(2))
        end if
        call h5_createDataSet(dump_fname(i), h5_dumpID, dump_hdf5Format(2), dump_hdf5DataSet(i), napx)

      case(3)
        ! Format 3:
        ! # ID turn s[m] x[mm] xp[mrad] y[mm] yp[mrad] z[mm] dE/E[1] ktrack
        if(dump_hdf5Format(3) == 0) then
          allocate(setFields(10))
          setFields(1)  = h5_dataField(name="ID",     type=h5_typeInt)
          setFields(2)  = h5_dataField(name="TURN",   type=h5_typeInt)
          setFields(3)  = h5_dataField(name="S",      type=h5_typeReal)
          setFields(4)  = h5_dataField(name="X",      type=h5_typeReal)
          setFields(5)  = h5_dataField(name="XP",     type=h5_typeReal)
          setFields(6)  = h5_dataField(name="Y",      type=h5_typeReal)
          setFields(7)  = h5_dataField(name="YP",     type=h5_typeReal)
          setFields(8)  = h5_dataField(name="dE/E",   type=h5_typeReal)
          setFields(9)  = h5_dataField(name="Z",      type=h5_typeReal)
          setFields(10) = h5_dataField(name="KTRACK", type=h5_typeInt)
          call h5_createFormat("dumpFormat3", setFields, dump_hdf5Format(3))
        end if
        call h5_createDataSet(dump_fname(i), h5_dumpID, dump_hdf5Format(3), dump_hdf5DataSet(i), napx)

      case(4)
        ! Format 4:
        ! # napx turn s[m] <x>[mm] <xp>[mrad] <y>[mm] <yp>[mrad] <z>[mm] <dE/E>[1]
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

      case(5)
        ! Format 5:
        ! # napx turn s[m] <x>[mm] <xp>[mrad] <y>[mm] <yp>[mrad] <z>[mm] <dE/E>[1]
        ! <x^2> <x*xp> <x*y> <x*yp> <x*z> <x*(dE/E)>
        ! <xp^2> <xp*y> <xp*yp> <xp*z> <xp*(dE/E)>
        ! <y^2> <y*yp> <y*z> <y*(dE/E)>
        ! <yp^2> <yp*z> <yp*(dE/E)>
        ! <z^2> <z*(dE/E)>
        ! <(dE/E)^2>
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

      end select

      if(allocated(setFields) .eqv. .true.) deallocate(setFields)
    end if
  end do
#endif

end subroutine dump_initialise

! ================================================================================================================================ !
!  A.Mereghetti, D.Sinuela-Pastor & P.Garcia Ortega, for the FLUKA Team
!  K.Sjobak, A.Santamaria, V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-05-08
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
  use mod_commont
  use mod_commonmn

#ifdef COLLIMAT
  use collimation
#endif

  implicit none

  ! interface variables:
  integer, intent(in) :: nturn, i, ix, unit, fmt
  logical, intent(in) ::  lhighprec
  real(kind=fPrec), intent(in) :: tasinv(6,6) ! normalization matrix in [mm,mrad,mm,mrad,mm,1]
  real(kind=fPrec), intent(in) :: loc_clo(6) ! closed orbit in [mm,mrad,mm,mrad,mm,1]

  ! Temporary variables
  integer j,k,l,m,n
  character(len=mNameLen) localBez

  real(kind=fPrec) localDcum
  integer localKtrack

  real(kind=fPrec) xyz_particle(6),nxyz_particle(6)
  real(kind=fPrec) xyz(6)
  real(kind=fPrec) xyz2(6,6)

#ifdef CR
  ! For accessing dumpfilepos
  integer dumpIdx
  if (unit .eq. dumpunit(0)) then
    ! ALL output must be on separate unit
    dumpIdx = 0
  elseif (unit .eq. dumpunit(-1)) then
    ! ALL output must be on separate unit
    dumpIdx = -1
  else
    dumpIdx = ix
  end if
#endif

  ! ------------------------------------------------------------------ !
  !  Format #0
  !  General format
  ! ------------------------------------------------------------------ !
  if (fmt == 0) then
    if (i == 0 .and. ix == 0) then
      localDcum = 0.0
      localBez = "StartDUMP"
    else
      localDcum = dcum(i)
      if (ktrack(i) .ne. 1) then
        localBez = bez(ix)
      else                ! BLOCs
        localBez = bezb(ic(i))
      end if
    end if
    if (lhighprec) then
      do j=1,napx
        write(unit,1981) nturn, i, ix, localBez, localDcum, xv(1,j)*c1m3, yv(1,j)*c1m3, xv(2,j)*c1m3, &
          yv(2,j)*c1m3, ejfv(j)*c1m3, (ejv(j)-e0)*c1e6, -c1m3*(sigmv(j)/clight)*(e0/e0f)
      end do
    else
      do j=1,napx
        write(unit,1982) nturn, i, ix, localBez, localDcum, xv(1,j)*c1m3, yv(1,j)*c1m3, xv(2,j)*c1m3, &
          yv(2,j)*c1m3, ejfv(j)*c1m3, (ejv(j)-e0)*c1e6, -c1m3*(sigmv(j)/clight)*(e0/e0f)
      end do
    end if
    write(unit,*) ''
    write(unit,*) ''

    ! Flush
    endfile (unit,iostat=ierro)
    backspace (unit,iostat=ierro)
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
      call h5_writeData(dump_hdf5DataSet(ix), 1, napx, nlostp)
      call h5_writeData(dump_hdf5DataSet(ix), 2, napx, nturn)
      call h5_writeData(dump_hdf5DataSet(ix), 3, napx, localDcum)
      call h5_writeData(dump_hdf5DataSet(ix), 4, napx, xv(1,:))
      call h5_writeData(dump_hdf5DataSet(ix), 5, napx, yv(1,:))
      call h5_writeData(dump_hdf5DataSet(ix), 6, napx, xv(2,:))
      call h5_writeData(dump_hdf5DataSet(ix), 7, napx, yv(2,:))
      call h5_writeData(dump_hdf5DataSet(ix), 8, napx, (ejv-e0)/e0)
      call h5_writeData(dump_hdf5DataSet(ix), 9, napx, localKtrack)
      call h5_finaliseWrite(dump_hdf5DataSet(ix))
    else
#endif
      if (lhighprec) then
        do j=1,napx
          write(unit,1983) nlostp(j),nturn,localDcum, &
            xv(1,j),yv(1,j),xv(2,j),yv(2,j),(ejv(j)-e0)/e0,localKtrack
        end do
      else
        do j=1,napx
          write(unit,1984) nlostp(j),nturn,localDcum, &
            xv(1,j),yv(1,j),xv(2,j),yv(2,j),(ejv(j)-e0)/e0,localKtrack
        end do
      end if

      ! Flush
      endfile (unit,iostat=ierro)
      backspace (unit,iostat=ierro)
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
  else if (fmt .eq. 2) then
    if (i.eq.0 .and. ix.eq.0) then
      localDcum = zero
      localKtrack = 0
    else
      localDcum = dcum(i)
      localKtrack = ktrack(i)
    end if
#ifdef HDF5
    if(h5_useForDUMP) then
      call h5_prepareWrite(dump_hdf5DataSet(ix), napx)
      call h5_writeData(dump_hdf5DataSet(ix), 1,  napx, nlostp)
      call h5_writeData(dump_hdf5DataSet(ix), 2,  napx, nturn)
      call h5_writeData(dump_hdf5DataSet(ix), 3,  napx, localDcum)
      call h5_writeData(dump_hdf5DataSet(ix), 4,  napx, xv(1,:))
      call h5_writeData(dump_hdf5DataSet(ix), 5,  napx, yv(1,:))
      call h5_writeData(dump_hdf5DataSet(ix), 6,  napx, xv(2,:))
      call h5_writeData(dump_hdf5DataSet(ix), 7,  napx, yv(2,:))
      call h5_writeData(dump_hdf5DataSet(ix), 8,  napx, sigmv)
      call h5_writeData(dump_hdf5DataSet(ix), 9,  napx, (ejv-e0)/e0)
      call h5_writeData(dump_hdf5DataSet(ix), 10, napx, localKtrack)
      call h5_finaliseWrite(dump_hdf5DataSet(ix))
    else
#endif
      if (lhighprec) then
        do j=1,napx
          write(unit,1985) nlostp(j),nturn,localDcum,xv(1,j),yv(1,j),xv(2,j),yv(2,j), &
            sigmv(j),(ejv(j)-e0)/e0,localKtrack
        end do
      else
        do j=1,napx
          write(unit,1986) nlostp(j),nturn,localDcum,xv(1,j),yv(1,j),xv(2,j),yv(2,j), &
            sigmv(j),(ejv(j)-e0)/e0,localKtrack
        end do
      end if

      ! Flush
      endfile (unit,iostat=ierro)
      backspace (unit,iostat=ierro)
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
  else if (fmt == 3) then
    if (i == 0 .and. ix == 0) then
      localDcum = zero
      localKtrack = 0
    else
      localDcum = dcum(i)
      localKtrack = ktrack(i)
    end if
#ifdef HDF5
    if(h5_useForDUMP) then
      call h5_prepareWrite(dump_hdf5DataSet(ix), napx)
      call h5_writeData(dump_hdf5DataSet(ix), 1,  napx, nlostp)
      call h5_writeData(dump_hdf5DataSet(ix), 2,  napx, nturn)
      call h5_writeData(dump_hdf5DataSet(ix), 3,  napx, localDcum)
      call h5_writeData(dump_hdf5DataSet(ix), 4,  napx, xv(1,:))
      call h5_writeData(dump_hdf5DataSet(ix), 5,  napx, yv(1,:))
      call h5_writeData(dump_hdf5DataSet(ix), 6,  napx, xv(2,:))
      call h5_writeData(dump_hdf5DataSet(ix), 7,  napx, yv(2,:))
      call h5_writeData(dump_hdf5DataSet(ix), 8,  napx, sigmv)
      call h5_writeData(dump_hdf5DataSet(ix), 9,  napx, (ejv-e0)/e0)
      call h5_writeData(dump_hdf5DataSet(ix), 10, napx, localKtrack)
      call h5_finaliseWrite(dump_hdf5DataSet(ix))
    else
#endif
      do j=1,napx
        write(unit) nlostp(j),nturn,localDcum,xv(1,j),yv(1,j),xv(2,j),yv(2,j), &
          sigmv(j),(ejv(j)-e0)/e0,localKtrack
      end do

      ! Flush
      endfile (unit,iostat=ierro)
      backspace (unit,iostat=ierro)
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
  else if (fmt .eq. 4) then
    if (i.eq.0 .and. ix.eq.0) then
      localDcum = zero
    else
      localDcum = dcum(i)
    end if
    do l=1,6
      xyz(l) = zero
    end do
    do j=1,napx
      xyz(1) = xyz(1) + xv(1,j)
      xyz(2) = xyz(2) + yv(1,j)
      xyz(3) = xyz(3) + xv(2,j)
      xyz(4) = xyz(4) + yv(2,j)
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
      if (lhighprec) then
        write(unit,1989) napx, nturn, localDcum, xyz(1),xyz(2),xyz(3),xyz(4),xyz(5),xyz(6)
      else
        write(unit,1990) napx, nturn, localDcum, xyz(1),xyz(2),xyz(3),xyz(4),xyz(5),xyz(6)
      end if

      ! Flush
      endfile (unit,iostat=ierro)
      backspace (unit,iostat=ierro)
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
        xyz(1) = xyz(1) + xv(1,j)
        xyz(2) = xyz(2) + yv(1,j)
        xyz(3) = xyz(3) + xv(2,j)
        xyz(4) = xyz(4) + yv(2,j)
        xyz(5) = xyz(5) + sigmv(j)
        xyz(6) = xyz(6) + xyz_particle(6)

        ! Beam matrix (don't calulate identical elements twice (symmetry))
        xyz2(1,1) = xyz2(1,1) + xv(1,j)*xv(1,j)
        xyz2(2,1) = xyz2(2,1) + xv(1,j)*yv(1,j)
        xyz2(3,1) = xyz2(3,1) + xv(1,j)*xv(2,j)
        xyz2(4,1) = xyz2(4,1) + xv(1,j)*yv(2,j)
        xyz2(5,1) = xyz2(5,1) + xv(1,j)*sigmv(j)
        xyz2(6,1) = xyz2(6,1) + xv(1,j)*xyz_particle(6)

        xyz2(2,2) = xyz2(2,2) + yv(1,j)*yv(1,j)
        xyz2(3,2) = xyz2(3,2) + yv(1,j)*xv(2,j)
        xyz2(4,2) = xyz2(4,2) + yv(1,j)*yv(2,j)
        xyz2(5,2) = xyz2(5,2) + yv(1,j)*sigmv(j)
        xyz2(6,2) = xyz2(6,2) + yv(1,j)*xyz_particle(6)

        xyz2(3,3) = xyz2(3,3) + xv(2,j)*xv(2,j)
        xyz2(4,3) = xyz2(4,3) + xv(2,j)*yv(2,j)
        xyz2(5,3) = xyz2(5,3) + xv(2,j)*sigmv(j)
        xyz2(6,3) = xyz2(6,3) + xv(2,j)*xyz_particle(6)

        xyz2(4,4) = xyz2(4,4) + yv(2,j)*yv(2,j)
        xyz2(5,4) = xyz2(5,4) + yv(2,j)*sigmv(j)
        xyz2(6,4) = xyz2(6,4) + yv(2,j)*xyz_particle(6)

        xyz2(5,5) = xyz2(5,5) + sigmv(j)*sigmv(j)
        xyz2(6,5) = xyz2(6,5) + sigmv(j)*xyz_particle(6)

        xyz2(6,6) = xyz2(6,6) + xyz_particle(6)*xyz_particle(6)
      end do

    else if (fmt == 6) then ! Canonical
      do j=1,napx
        xyz_particle(1) = xv(1,j)*c1m3                 !x:      [mm]   -> [m]
        xyz_particle(2) = (yv(1,j)*c1m3)*(one+dpsv(j)) !px:     [mrad] -> [1]
        xyz_particle(3) = xv(2,j)*c1m3                 !y:      [mm]   -> [m]
        xyz_particle(4) = (yv(2,j)*c1m3)*(one+dpsv(j)) !py:     [mrad] -> [1]
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
        write(unit,1991) napx, nturn, localDcum, xyz(1),xyz(2),xyz(3),xyz(4),xyz(5),xyz(6), &
          xyz2(1,1),xyz2(2,1),xyz2(3,1),xyz2(4,1),xyz2(5,1),xyz2(6,1),&
                    xyz2(2,2),xyz2(3,2),xyz2(4,2),xyz2(5,2),xyz2(6,2),&
                              xyz2(3,3),xyz2(4,3),xyz2(5,3),xyz2(6,3),&
                                        xyz2(4,4),xyz2(5,4),xyz2(6,4),&
                                                  xyz2(5,5),xyz2(6,5),&
                                                            xyz2(6,6)
      else
        write(unit,1992) napx, nturn, localDcum, xyz(1),xyz(2),xyz(3),xyz(4),xyz(5),xyz(6), &
          xyz2(1,1),xyz2(2,1),xyz2(3,1),xyz2(4,1),xyz2(5,1),xyz2(6,1),&
                    xyz2(2,2),xyz2(3,2),xyz2(4,2),xyz2(5,2),xyz2(6,2),&
                              xyz2(3,3),xyz2(4,3),xyz2(5,3),xyz2(6,3),&
                                        xyz2(4,4),xyz2(5,4),xyz2(6,4),&
                                                  xyz2(5,5),xyz2(6,5),&
                                                            xyz2(6,6)
      end if

      ! Flush
      endfile (unit,iostat=ierro)
      backspace (unit,iostat=ierro)
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
  else if (fmt .eq. 7 .or. fmt .eq. 8 .or. fmt .eq. 9) then
    if (i.eq.0 .and. ix.eq.0) then
      localDcum = zero
      localKtrack = 0
    else
      localDcum = dcum(i)
      localKtrack = ktrack(i)
    end if

    ! initialize parameters for writing of beam moments
    do l=1,6
      xyz(l) = zero
      do k=1,6
        xyz2(l,k) = zero
      end do
    end do

    ! normalize particle coordinates
    do j=1,napx
      xyz_particle(1) = xv(1,j)
      xyz_particle(2) = yv(1,j)
      xyz_particle(3) = xv(2,j)
      xyz_particle(4) = yv(2,j)
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
      ! Initialize nxyz
      do m=1,6
        nxyz_particle(m)=zero
      end do
      do m=1,6
        do n=1,6
          nxyz_particle(m)=nxyz_particle(m)+tasinv(m,n)*xyz_particle(n)
        end do
        ! a) convert nxyzv(6) to 1.e-3 sqrt(m)
        ! unit: nx,npx,ny,npy,nsig,ndelta all in [1.e-3 sqrt(m)]
        if (m.eq.6) then
          nxyz_particle(m)=nxyz_particle(m)*c1e3
        end if
      end do

      if (fmt .eq. 7) then
        if (lhighprec) then
          write(unit,1985) nlostp(j),nturn,localDcum, &
            nxyz_particle(1),nxyz_particle(2),nxyz_particle(3),nxyz_particle(4),nxyz_particle(5),nxyz_particle(6),localKtrack
        else
          write(unit,1986) nlostp(j),nturn,localDcum, &
            nxyz_particle(1),nxyz_particle(2),nxyz_particle(3),nxyz_particle(4),nxyz_particle(5),nxyz_particle(6),localKtrack
        end if

      else if (fmt .eq. 8) then
        write(unit) nlostp(j),nturn,localDcum, &
          nxyz_particle(1),nxyz_particle(2),nxyz_particle(3),nxyz_particle(4),nxyz_particle(5),nxyz_particle(6),localKtrack

      else if (fmt .eq. 9) then
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

    if (fmt .eq. 7) then
      ! Flush
      endfile (unit,iostat=ierro)
      backspace (unit,iostat=ierro)
#ifdef CR
      dumpfilepos(dumpIdx) = dumpfilepos(dumpIdx)+napx
#endif
    else if (fmt .eq. 8) then
      ! Flush
      endfile (unit,iostat=ierro)
      backspace (unit,iostat=ierro)
#ifdef CR
      dumpfilepos(dumpIdx) = dumpfilepos(dumpIdx)+napx
#endif
    else if (fmt .eq. 9) then
      ! Normalize to get averages
      xyz = xyz/napx

      xyz2(:,1)  = xyz2(:,1) /napx
      xyz2(2:,2) = xyz2(2:,2)/napx
      xyz2(3:,3) = xyz2(3:,3)/napx
      xyz2(4:,4) = xyz2(4:,4)/napx
      xyz2(5:,5) = xyz2(5:,5)/napx
      xyz2(6,6)  = xyz2(6,6) /napx

      if (lhighprec) then
        write(unit,1991) napx,nturn,localDcum,xyz(1),xyz(2),xyz(3),xyz(4),xyz(5),xyz(6), &
          xyz2(1,1),xyz2(2,1),xyz2(3,1),xyz2(4,1),xyz2(5,1),xyz2(6,1),&
                    xyz2(2,2),xyz2(3,2),xyz2(4,2),xyz2(5,2),xyz2(6,2),&
                              xyz2(3,3),xyz2(4,3),xyz2(5,3),xyz2(6,3),&
                                        xyz2(4,4),xyz2(5,4),xyz2(6,4),&
                                                  xyz2(5,5),xyz2(6,5),&
                                                            xyz2(6,6)
      else
        write(unit,1992) napx,nturn,localDcum,xyz(1),xyz(2),xyz(3),xyz(4),xyz(5),xyz(6), &
          xyz2(1,1),xyz2(2,1),xyz2(3,1),xyz2(4,1),xyz2(5,1),xyz2(6,1),&
                    xyz2(2,2),xyz2(3,2),xyz2(4,2),xyz2(5,2),xyz2(6,2),&
                              xyz2(3,3),xyz2(4,3),xyz2(5,3),xyz2(6,3),&
                                        xyz2(4,4),xyz2(5,4),xyz2(6,4),&
                                                  xyz2(5,5),xyz2(6,5),&
                                                            xyz2(6,6)
      end if

      ! Flush
      endfile (unit,iostat=ierro)
      backspace (unit,iostat=ierro)
#ifdef CR
      dumpfilepos(dumpIdx) = dumpfilepos(dumpIdx)+1
#endif
    end if

  ! Unrecognized format fmt
  else
    write (lout,"(2(a,i0))") "DUMP> ERROR Format ",fmt," not understood for unit ",unit
    call prror(-1)
  end if

  return

1981 format (3(1X,I8),1X,A16,1X,F12.5,7(1X,1PE25.18)) !fmt 0 / hiprec
1982 format (3(1X,I8),1X,A16,1X,F12.5,7(1X,1PE16.9))  !fmt 0 / not hiprec

1983 format (2(1x,I8),1X,F12.5,5(1X,1PE25.18),1X,I8)  !fmt 1 / hiprec
1984 format (2(1x,I8),1X,F12.5,5(1X,1PE16.9),1X,I8)   !fmt 1 / not hiprec

1985 format (2(1x,I8),1X,F12.5,6(1X,1PE25.18),1X,I8)  !fmt 2&7 / hiprec
1986 format (2(1x,I8),1X,F12.5,6(1X,1PE16.9),1X,I8)   !fmt 2&7 / not hiprec

1989 format (2(1x,I8),1X,F12.5,6(1X,1PE25.18))        !fmt 4 / hiprec
1990 format (2(1x,I8),1X,F12.5,6(1X,1PE16.9))         !fmt 4 / not hiprec

1991 format (2(1x,I8),1X,F12.5,27(1X,1PE25.18))       !fmt 5&6 / hiprec
1992 format (2(1x,I8),1X,F12.5,27(1X,1PE16.9))        !fmt 5&6 / not hiprec

end subroutine dump_beam_population

! ================================================================================================================================ !
subroutine dump_comnul

  implicit none

  integer i1,i2,i3,i

  ldumphighprec = .false.
  ldumpfront    = .false.

  do i1=-1,nblz
    do i2=1,6
      dumpclo(i1,i2)=0
      do i3=1,6
        dumptas(i1,i2,i3)=0
      end do
    end do
  end do

  do i=-1,nele
    ldump(i)    = .false.
    ndumpt(i)   = 0
    dumpfirst(i) = 0
    dumplast(i)  = 0
    dumpunit(i) = 0
    dumpfmt(i)  = 0
#ifdef CR
    dumpfilepos(i) = -1
#endif
  end do

end subroutine dump_comnul

! ================================================================================================================================ !
!  Begin Checkpoint Restart
! ================================================================================================================================ !
#ifdef CR

! ================================================================================================================================ !
subroutine dump_crcheck_readdata(fileunit, readerr)

  implicit none

  integer, intent(in) :: fileunit
  logical, intent(out) :: readerr

  integer j

  read(fileunit,err=100,end=100) (dumpfilepos_cr(j),j=-1,nele)

  readerr = .false.
  return

100 continue
  readerr = .true.

end subroutine dump_crcheck_readdata

! ================================================================================================================================ !
subroutine dump_crcheck_positionFiles

  use crcoall
  use string_tools
  use mod_common

  implicit none

  ! For skipping through binary DUMP files (format 3&8)
  integer tmp_ID, tmp_nturn, tmp_ktrack
  real(kind=fPrec) tmp_dcum, tmp_x,tmp_xp,tmp_y,tmp_yp,tmp_sigma,tmp_dEE

  integer i,j
  logical lerror,lopen
  character(len=256) filename
  character(len=1024) arecord

  do i=-1, il
    if (ldump(i)) then
      write(93,*) "SIXTRACR CRCHECK REPOSITIONING DUMP file"
      if (i .gt. 0) then
        write(93,*) "element=",bez(i), "unit=",dumpunit(i)," filename='"//trim(chr_trimZero(dump_fname(i)))// &
                    "' format=",dumpfmt(i)
      else if (i.eq.0) then
        write(93,*) "element=","ALL" , "unit=",dumpunit(i)," filename='"//trim(chr_trimZero(dump_fname(i)))// &
                    "' format=",dumpfmt(i)
      else if(i .eq. -1) then
        write(93,*) "element=","StartDump" , "unit=",dumpunit(i)," filename='"//trim(chr_trimZero(dump_fname(i)))// &
                    "' format=",dumpfmt(i)
      else
        write(93,*) "Error - index=",i,"is unknown"
        goto 111
      end if
      flush(93)

      inquire( unit=dumpunit(i), opened=lopen )
      if (dumpfmt(i).ne.3 .and. dumpfmt(i).ne.8) then ! ASCII
        if (.not. lopen) then
#ifdef BOINC
          call boincrf(dump_fname(i),filename)
          open(dumpunit(i),file=filename, status='old',form='formatted',action='readwrite')
#else
          open(dumpunit(i),file=dump_fname(i), status='old',form='formatted',action='readwrite')
#endif
        end if

        dumpfilepos(i) = 0
        do j=1,dumpfilepos_cr(i)
702       read(dumpunit(i),'(a1024)',end=111,err=111,iostat=ierro) arecord
          dumpfilepos(i) = dumpfilepos(i) + 1
        end do

      else                         ! BINARY (format = 3 & 8)
        if (.not. lopen) then
#ifdef BOINC
          call boincrf(dump_fname(i),filename)
          open(dumpunit(i),file=filename,status='old',form='unformatted',action='readwrite')
#else
          open(dumpunit(i),file=dump_fname(i),status='old',form='unformatted',action='readwrite')
#endif
        end if
        dumpfilepos(i) = 0
        do j=1,dumpfilepos_cr(i)
703       read(dumpunit(i),end=111,err=111,iostat=ierro) &
            tmp_ID,tmp_nturn,tmp_dcum,tmp_x,tmp_xp,tmp_y,tmp_yp,tmp_sigma,tmp_dEE,tmp_ktrack
          dumpfilepos(i) = dumpfilepos(i) + 1
        end do
      end if
    end if
  end do

  ! Crop DUMP files (if used by multiple DUMPs,
  ! the actual position is the sum of the dumpfileposes
  do i=0,il
    if (ldump(i)) then
      ! This is not a FLUSH!
      endfile (dumpunit(i),iostat=ierro)

      ! Change from 'readwrite' to 'write'
      close(dumpunit(i))
      if (dumpfmt(i).ne.3 .and. dumpfmt(i).ne.8) then ! ASCII
#ifdef BOINC
        call boincrf(dump_fname(i),filename)
        open(dumpunit(i),file=filename, status='old',position='append',form='formatted',action='write')
#else
        open(dumpunit(i),file=dump_fname(i), status='old',position='append',form='formatted',action='write')
#endif
      else ! Binary (format = 3)
#ifdef BOINC
        call boincrf(dump_fname(i),filename)
        open(dumpunit(i),file=filename, status='old',position='append',form='unformatted',action='write')
#else
        open(dumpunit(i),file=dump_fname(i), status='old',position='append',form='unformatted',action='write')
#endif
      end if
    end if
  end do

  return

111 continue
  write(93,*) 'SIXTRACR DUMP_CRCHECK_POSITIONFILES *** ERROR *** reading DUMP file#', dumpunit(i),' iostat=',ierro
  write(93,*) 'dumpfilepos=',dumpfilepos(i),' dumpfilepos_cr=',dumpfilepos_cr(i)
  flush(93)
  write(lout,"(a)") "SIXTRACR> ERROR DUMP_CRCHECK_POSITIONFILES failure positioning DUMP file"
  call prror(-1)

end subroutine dump_crcheck_positionFiles

! ================================================================================================================================ !
subroutine dump_crpoint(fileunit,lerror,ierro)

  use parpro !nele
  implicit none

  integer, intent(in) :: fileunit
  logical, intent(out) :: lerror
  integer, intent(inout) :: ierro
  integer j

  write(95,err=100,iostat=ierro) (dumpfilepos(j),j=-1,nele)
  return

100 continue
  lerror = .true.
  return

end subroutine dump_crpoint
! ================================================================================================================================ !

#endif
! ================================================================================================================================ !
!  End Checkpoint Restart
! ================================================================================================================================ !
end module dump
