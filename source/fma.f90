module fma

  implicit none

  integer, parameter     :: fma_max           = 200                ! Max. number of FMAs
  integer, parameter     :: fma_nturn_max     = 10000              ! Max. number of turns used for fft
  integer, private, save :: fma_numfiles      = 0                  ! Number of FMAs
  logical, public,  save :: fma_flag          = .false.            ! FMA input block exists
  logical, private, save :: fma_writeNormDUMP = .true.             ! Writing out the normalised DUMP files
  character(len=:), allocatable, private, save :: fma_fname(:)     ! Name of input file from dump
  character(len=:), allocatable, private, save :: fma_method(:)    ! Method used to find the tunes
  integer,          allocatable, private, save :: fma_first(:)     ! First turn used for FMA
  integer,          allocatable, private, save :: fma_last(:)      ! Last turn used for FMA
  integer,          allocatable, private, save :: fma_norm_flag(:) ! normalise phase space before FFT

contains

subroutine fma_allocate
  use mod_alloc
  use parpro
  call alloc(fma_fname,  mFileName, fma_max, " ", "fma_fname")
  call alloc(fma_method, mStrLen,   fma_max, " ", "fma_method")
  call alloc(fma_first,             fma_max, 0,   "fma_first")
  call alloc(fma_last,              fma_max, 0,   "fma_last")
  call alloc(fma_norm_flag,         fma_max, 1,   "fma_norm_flag")
end subroutine fma_allocate

subroutine fma_parseInputline(inLine,iErr)

  use crcoall
  use string_tools

  implicit none

  character(len=*), intent(in)    :: inLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit
  logical spErr, cErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lerr,"(a)") "FMA> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if
  if(nSplit == 0) return

  if(lnSplit(1) == "NoNormDUMP") then
    fma_writeNormDUMP = .false.
    return
  endif

  if(fma_numfiles >= fma_max) then
    write(lerr,"(a,i0,a)") "FMA> ERROR You can only do ",fma_max," number of FMAs"
    iErr = .true.
    return
  end if

  fma_numfiles=fma_numfiles+1

  if(nSplit == 1 .or. nSplit == 4 .or. nSplit >= 6) then
    write(lerr,"(a,i0)") "FMA> ERROR Wrong number of input parameters. Expected 2, 3 or 5, got ",nSplit
    iErr = .true.
    return
  end if

  fma_fname(fma_numfiles)  = trim(lnSplit(1))
  fma_method(fma_numfiles) = trim(lnSplit(2))
#ifndef NAFF
  if(fma_method(fma_numfiles) == "NAFF") then
    write(lerr,"(a)") "FMA> ERROR NAFF requested, but SixTrack was not built with the NAFF flag."
    iErr = .true.
    return
  end if
#endif
  if(nSplit == 2) then
    fma_norm_flag(fma_numfiles) = 1 ! default: normalise phase space
  else if(nSplit == 3) then
    call chr_cast(lnSplit(3),fma_norm_flag(fma_numfiles),cErr)
  else if(nSplit == 5) then
    call chr_cast(lnSplit(3),fma_norm_flag(fma_numfiles),cErr)
    call chr_cast(lnSplit(4),fma_first(fma_numfiles),    cErr)
    call chr_cast(lnSplit(5),fma_last(fma_numfiles),     cErr)
  end if

  ! Input sanity checks
  if(.not.(&
          fma_method(fma_numfiles) == "TUNELASK"  &
     .or. fma_method(fma_numfiles) == "TUNEFFTI"  &
     .or. fma_method(fma_numfiles) == "TUNEFFT"   &
     .or. fma_method(fma_numfiles) == "TUNEAPA"   &
     .or. fma_method(fma_numfiles) == "TUNEFIT"   &
     .or. fma_method(fma_numfiles) == "TUNENEWT"  &
     .or. fma_method(fma_numfiles) == "TUNEABT2"  &
     .or. fma_method(fma_numfiles) == "TUNEABT"   &
     .or. fma_method(fma_numfiles) == "TUNENEWT1" &
     .or. fma_method(fma_numfiles) == "NAFF"      &
    )) then
    write(lerr,"(a,i0)") "FMA> ERROR The method '"//trim(fma_method(fma_numfiles))//"' is unknown. FMA index = ",fma_numfiles
    write(lerr,"(a)")    "FMA>       Please use one of TUNELASK, TUNEFFTI, TUNEFFT, "// &
      "TUNEAPA, TUNEFIT, TUNENEWT, TUNEABT, TUNEABT2, TUNENEWT1, NAFF."
    write(lout,"(a)")    "FMA>       Note that it is case-sensitive, so use uppercase only."
    iErr = .true.
    return
  end if

  if(.not.(fma_norm_flag(fma_numfiles).eq.0 .or. fma_norm_flag(fma_numfiles).eq.1)) then
    write(lerr,"(2(a,i0))") "FMA> ERROR Expected fma_norm_flag = 1 or 0, Got: ",fma_norm_flag(fma_numfiles),&
      ", FMA index =",fma_numfiles
    iErr = .true.
    return
  end if

  fma_flag = .true.

end subroutine fma_parseInputline

! ================================================================================================ !
!  FMA POSTPROCESSING
!  M.Fitterer, R. De Maria, K.Sjobak, BE-ABP-HSS
!  Updated by: V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-07-23
!  purpose: return files used for fma analysis
!           -> calculate particle amplitudes and tunes using the
!              normalised coordinates for input files
!              fma_fname(fma_numfiles)
!  output format: q1,q2,q3,eps1_min,eps2_min,eps3_min,eps1_max,
!                 eps2_max,eps3_max,eps1_avg, eps2_avg,eps3_avg,
!                 eps1_0,eps2_0,eps3_0,phi1_0,phi2_0,phi3_0
! ================================================================================================ !
subroutine fma_postpr

  use floatPrecision
  use string_tools
  use mathlib_bouncer
  use platofma
  use dump, only : dump_fname, dumpfmt, dumpunit, dumpfirst, dumplast, dump_getTasMatrix
  use numerical_constants, only : zero, one, c1e3

  use crcoall
  use parpro
  use mod_common
  use mod_common_track
  use mod_alloc
  use mod_units

  implicit none

  character(len=:), allocatable :: lnSplit(:)
  character(len=mInputLn)       :: rLine
  character(len=mStrLen)        :: ch, ch1
  integer nSplit, fmaUnit, tmpUnit, dumpLastTurn, numModes
  logical spErr, fErr, cErr, isOpen, fExist
  integer i,j,k,l,m,n

  ! Current turn no (particle, rel. turn no)
  integer, allocatable :: turn(:,:)
  ! Number of turns to analyze for this particle
  integer, allocatable :: nturns(:)
  ! Have we written a normDump file for this element before?
  logical, allocatable :: hasNormDumped(:)
  ! Number of turns used for fft for this FMA
  integer, allocatable :: fma_nturn(:)
  ! Phase space variables (x,x',y,y',z,dE/E) [mm,mrad,mm,mrad,mm,1.e-3]
  real(kind=fPrec), allocatable :: xyzv(:,:,:)
  ! Normalised phase space variables [sqrt(m) 1.e-3]
  real(kind=fPrec), allocatable :: nxyzv(:,:,:)
  ! Normalised emittances
  real(kind=fPrec), allocatable :: epsnxyzv(:,:,:)
  real(kind=fPrec) dumptas(6,6), dumpclo(6), dumptasinv(6,6)

#ifdef NAFF
  interface
    real(c_double) function tunenaff(x,xp,maxn,plane_idx,norm_flag,fft_naff) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), intent(in), dimension(1) :: x,xp
      integer(c_int), intent(in), value :: maxn, plane_idx, norm_flag
      real(c_double), intent(in), value :: fft_naff
    end function tunenaff
  end interface
#endif

  ! need to pass a single dimension array to naff,
  !  since the stride in the xyzv/nxyzv arrays are difficult to pass correctly to c++.
  ! (we can't interpret the struct that fortran is passing us;
  !  see the naff_interface.cpp for more info                 )
  real(kind=fPrec), allocatable :: naff_xyzv1(:)
  real(kind=fPrec), allocatable :: naff_xyzv2(:)

  ! dummy variables for readin + normalisation + loops
  integer :: id,kt,counter,thisturn
  real(kind=fPrec) :: pos, fft_naff
  real(kind=fPrec), dimension(6) :: xyzvdummy,xyzvdummy2,nxyzvdummy !phase space variables x,x',y,y',sig,delta
  real(kind=fPrec), dimension(3) :: q123 !tune q1,q2,q3
  real(kind=fPrec), dimension(3) :: eps123_0,eps123_min,eps123_max,eps123_avg !initial,minimum,maximum,average emittance
  real(kind=fPrec), dimension(3) :: phi123_0  !initial phase

  call alloc(turn,         napx,fma_nturn_max,   0,      "turn")
  call alloc(nturns,       napx,                 0,      "nturns")
  call alloc(hasNormDumped,nele,                 .false.,"hasNormDumped",-1)
  call alloc(fma_nturn,    fma_max,              0,      "fma_nturn")
  call alloc(xyzv,         napx,fma_nturn_max,6, zero,   "xyzv")
  call alloc(nxyzv,        napx,fma_nturn_max,6, zero,   "nxyzv")
  call alloc(epsnxyzv,     napx,fma_nturn_max,3, zero,   "epsnxyzv")
  call alloc(naff_xyzv1,   fma_nturn_max,        zero,   "naff_xyzv1")
  call alloc(naff_xyzv2,   fma_nturn_max,        zero,   "naff_xyzv2")

  ! fma_six = data file for storing the results of the FMA analysis
  call f_requestUnit("fma_sixtrack",fmaUnit)
  call f_open(unit=fmaUnit,file="fma_sixtrack",formatted=.true.,mode="w",err=fErr,status="replace")
  if(fErr) then
    write(lerr, "(a)") "FMA> ERROR Cannot open file 'fma_sixtrack' for writing."
    call prror
  end if

  if(idp == 0 .or. ition == 0) then
    numModes = 2 ! 4D tracking
  else
    numModes = 3 ! 6D tracking
  end if

  ! write the header of fma_sixtrack
  write(fmaUnit,"(a)") "# eps1*,eps2*,eps3* all in 1.e-6*m, phi* [rad]"
  write(fmaUnit,"(a)") "# inputfile method id q1 q2 q3 eps1_min eps2_min eps3_min eps1_max eps2_max eps3_max "//&
    "eps1_avg eps2_avg eps3_avg eps1_0 eps2_0 eps3_0 phi1_0 phi2_0 phi3_0 norm_flag first_turn last_turn"

  ! Start FMA analysis: loop over all files, calculate tunes, write output file
  do i=1,fma_numfiles
    fExist = .false.
    do j=-1,nele ! START: loop over dump files = loop over single elements
      if(fma_fname(i) == dump_fname(j)) then
        call dump_getTasMatrix(j,dumptasinv,dumptas, dumpclo)
        fExist = .true.
        write(lout,"(3(a,i0))") "FMA> Start FMA analysis using file '"//trim(fma_fname(i))//&
          "': Number of particles = ",napx,", first turn = ",fma_first(i),", last turn = ",fma_last(i)

        ! Check the format, if dumpfmt != 2,3 (physical) or 7,8 (normalised) then abort
        if(.not.(dumpfmt(j) == 2 .or. dumpfmt(j) == 3 .or. dumpfmt(j) == 7 .or. dumpfmt(j) == 8)) then
          write(lerr,"(a)") "FMA> ERROR Input file has wrong format. Choose format 2, 3, 7 or 8 in DUMP block."
          call prror
        end if

        ! Open dump file for reading, resume to original position before exiting the subroutine
        inquire(unit=dumpunit(j),opened=isOpen)
        if(isOpen) then
          call f_close(dumpunit(j))
        else ! File has to be open if nothing went wrong
          write(lerr,"(a)") "FMA> ERROR Expected file '"//trim(dump_fname(j))//"' to be open."
          call prror
        end if

        if(dumpfmt(j) == 2 .or. dumpfmt(j) == 7) then
          call f_open(unit=dumpunit(j),file=dump_fname(j),formatted=.true.,mode="r",err=fErr,status="old")
          if(fErr) then
            write(lerr,"(a,i0,a)") "FMA> ERROR Opening file 'NORM_"//trim(dump_fname(j))//"' (dumpfmt=",dumpfmt(j),")"
            call prror
          end if
        else if(dumpfmt(j) == 3 .or. dumpfmt(j) == 8) then
          call f_open(unit=dumpunit(j),file=dump_fname(j),formatted=.false.,mode="r",err=fErr,status="old")
          if(fErr) then
            write(lerr,"(a,i0,a)") "FMA> ERROR Opening file 'NORM_"//trim(dump_fname(j))//"' (dumpfmt=",dumpfmt(j),")"
            call prror
          end if
        else
          write(lerr,"(a,i0,a)") "FMA> ERROR Got dumpfmt = ",dumpfmt(j),", but expected 2,3,7 or 8."
          call prror
        end if

        ! Define first/last turn for FMA
        ! - If first and last turn are not defined in FMA block, take all turns saved in DUMP file
        if(fma_first(i) == 0 .and. fma_last(i) == 0) then
          fma_first(i) = dumpfirst(j)
          fma_last(i)  = dumplast(j)
        end if

        ! If -1, take the last turn of the dump file or the maximum number of turns if dumplast = -1
        if(fma_last(i) == -1) then
          if(dumplast(j) == -1) then
            fma_last(i) = numl
          else
            fma_last(i) = dumplast(j)
          end if
        end if

        ! Now check that first turn are compatible with turns saved in dump file
        if(fma_first(i) < dumpfirst(j)) then
          write(lerr,"(2(a,i0))") "FMA> ERROR First turn in FMA block is smaller than first turn in DUMP block: "//&
            "fma_first = ",fma_first(i)," < dumpfirst = ",dumpfirst(j)
          call prror
        end if

        ! Now check last turn
        ! - If fma_last = -1, we already have fma_last = numl check if fma_last < 0 and !=-1
        if(fma_last(i) <= 0) then
          write(lerr,"(a,i0)") "FMA> ERROR Last turn in FMA block must be -1 or a positive integer, "//&
            "but fma_last = ",fma_last(i)
          call prror
        end if

        ! If fma_last >0 check that fma_last < dump_last
        if(dumplast(j) == -1) then
          if(fma_last(i) > numl) then
            write(lerr,"(2(a,i0))") "FMA> ERROR Last turn in FMA block is larger than number of turns tracked. "//&
              " fma_last = ",fma_last(i)," > turns tracked = ",numl
            call prror
          end if
        else
          if(fma_last(i) > dumplast(j)) then
            write(lerr,"(2(a,i0))") "FMA> ERROR Last turn in FMA block is larger than number of turns tracked in DUMP block. "//&
              "fma_last = ",fma_last(i)," > dumplast = ",dumplast(j)
          end if
        end if

        ! Now we can set the number of turns used for the FMA required for the PLATO routines
        fma_nturn(i) = fma_last(i)-fma_first(i)+1
        do m=1,napx
          nturns(m) = fma_nturn(i)
        end do
        if(fma_nturn(i) > fma_nturn_max) then
          write(lerr,"(a,i0,a,i0,a)") "FMA> ERROR Only ",fma_nturn_max," turns allowed for fma, but ",fma_nturn(i)," used."
          write(lerr,"(a,i0)")        "FMA>       -> reset fma_nturn_max > ",fma_nturn_max
          call prror
        end if

        ! Now we can start reading in the file
        if(dumpfmt(j) == 2 .or. dumpfmt(j) == 7) then ! ASCII -> skip the header
          counter = 1
          do
            read(dumpunit(j),"(a)",iostat=ierro) ch
            if(ierro /= 0) then
              write(lerr,"(a)") "FMA> ERROR Reading file '"//trim(dump_fname(j))//"'"
              call prror
            end if
            ch1=adjustl(trim(ch))
            if(ch1(1:1) /= "#") exit
            if(counter > 500) then
              write(lerr,"(a)") "FMA> ERROR Something is wrong with your dumpfile '"//trim(dump_fname(j))//&
                "'> Found more than 500 header lines."
              call prror
            end if
            counter = counter+1
          end do
          backspace(dumpunit(j),iostat=ierro)
        end if

        ! Format 7 and 8 use normalised coordinates -> set fma_norm_flag =1
        if(dumpfmt(j) == 7 .or. dumpfmt(j) == 8) then
          if(fma_norm_flag(i) /= 1 ) then
            ! For format 7 and 8, the particles are already normalised by the DUMP block
            write(lerr,"(a,i0)") "FMA> ERROR For FMA #",i
            write(lerr,"(a)")    "FMA>       Cannot do FMA on physical coordinates if normalised DUMP is used (format 7 or 8)"
            call prror
          end if
        else ! Reading physical coordinates
          if(fma_norm_flag(i) == 1 ) then
            ! Have a matrix that's not zero (i.e. did we put a 6d LINE block?)
            if(dumptas(1,1) == zero .and. dumptas(1,2) == zero .and. &
               dumptas(1,3) == zero .and. dumptas(1,4) == zero) then
              write(lerr,"(a)") "FMA> ERROR The normalisation matrix appears to not be set? Did you forget to put a 6D LINE block?"
              call prror
            end if
            if(idp == 0 .or. ition == 0) then ! We're in the 4D case
              if(j /= -1) then ! Not at StartDUMP
                write(lerr,"(a)") "FMA> ERROR normalised coordinates: 4D only supported for StartDUMP."
                call prror
              end if
            end if
          end if
        end if
        ! Now we have done all checks

        !Normalized copy of the dump
        if(fma_writeNormDUMP .and. .not.(dumpfmt(j) == 7 .or. dumpfmt(j) == 8) .and. .not.hasNormDumped(j)) then
          ! Get a file unit, if needed
          call f_requestUnit("NORM_"//dump_fname(j),tmpUnit)
          call f_open(unit=tmpUnit,file="NORM_"//dump_fname(j),formatted=.true.,mode="w",err=fErr,status="replace")
          if(fErr) then
            write(lerr,"(a)") "FMA> ERROR Opening file 'NORM_"//trim(dump_fname(j))//"'"
            call prror
          end if

          ! Write the file headers
          write(lout,"(a)") "FMA> Writing normalised DUMP for '"//trim(dump_fname(j))// "' ..."
          ! Dump normalised particle amplitudes for debugging (tmpUnit)
          !  units: dumptas, dumptasinv, dumpclo [mm,mrad,mm,mrad,1]
          !  note: closed orbit dumpclo already converted in linopt part to [mm,mrad,mm,mrad,1]
          !        tas matrix in linopt part in [mm,mrad,mm,mrad,1.e-3]

          ! - write closed orbit in header of file with normalised phase space coordinates (tmpUnit)
          !   units: x,xp,y,yp,sig,dp/p = [mm,mrad,mm,mrad,1] (note: units are already changed in linopt part)
          write(tmpUnit,"(a,1x,6(1x,1pe16.9))") "# closorb", &
            dumpclo(1),dumpclo(2),dumpclo(3), &
            dumpclo(4),dumpclo(5),dumpclo(6)

          ! - write tas-matrix and its inverse in header of file with normalised phase space coordinates (tmpUnit)
          !   units: x,px,y,py,sig,dp/p [mm,mrad,mm,mrad,1]
          write(tmpUnit,"(a,1x,36(1x,1pe16.9))") "# tamatrix [mm,mrad,mm,mrad,1]", &
            dumptas(1,1),dumptas(1,2),dumptas(1,3),dumptas(1,4), &
            dumptas(1,5),dumptas(1,6),dumptas(2,1),dumptas(2,2), &
            dumptas(2,3),dumptas(2,4),dumptas(2,5),dumptas(2,6), &
            dumptas(3,1),dumptas(3,2),dumptas(3,3),dumptas(3,4), &
            dumptas(3,5),dumptas(3,6),dumptas(4,1),dumptas(4,2), &
            dumptas(4,3),dumptas(4,4),dumptas(4,5),dumptas(4,6), &
            dumptas(5,1),dumptas(5,2),dumptas(5,3),dumptas(5,4), &
            dumptas(5,5),dumptas(5,6),dumptas(6,1),dumptas(6,2), &
            dumptas(6,3),dumptas(6,4),dumptas(6,5),dumptas(6,6)

          write(tmpUnit,"(a,1x,36(1x,1pe16.9))") "# inv(tamatrix)", &
            dumptasinv(1,1),dumptasinv(1,2),dumptasinv(1,3), &
            dumptasinv(1,4),dumptasinv(1,5),dumptasinv(1,6), &
            dumptasinv(2,1),dumptasinv(2,2),dumptasinv(2,3), &
            dumptasinv(2,4),dumptasinv(2,5),dumptasinv(2,6), &
            dumptasinv(3,1),dumptasinv(3,2),dumptasinv(3,3), &
            dumptasinv(3,4),dumptasinv(3,5),dumptasinv(3,6), &
            dumptasinv(4,1),dumptasinv(4,2),dumptasinv(4,3), &
            dumptasinv(4,4),dumptasinv(4,5),dumptasinv(4,6), &
            dumptasinv(5,1),dumptasinv(5,2),dumptasinv(5,3), &
            dumptasinv(5,4),dumptasinv(5,5),dumptasinv(5,6), &
            dumptasinv(6,1),dumptasinv(6,2),dumptasinv(6,3), &
            dumptasinv(6,4),dumptasinv(6,5),dumptasinv(6,6)

          write(tmpUnit,"(a)") "# id turn pos[m] nx[1.e-3 sqrt(m)] npx[1.e-3 sqrt(m)] ny[1.e-3 sqrt(m)] "//&
            "npy[1.e-3 sqrt(m)] nsig[1.e-3 sqrt(m)] ndp/p[1.e-3 sqrt(m)] kt"
        end if ! END IF fma_writeNormDUMP

        ! Read in particle amplitudes a(part,turn), x,xp,y,yp,sigma,dE/E [mm,mrad,mm,mrad,mm,1]
        ! TODO: This logic breaks apart if there are particle losses;
        !  it is checked for, but it only triggers a "call prror".

        ! If normalization within FMA, we now have to always write the full NORM_* file
        ! Otherwise  one would overwrite the NORM_* file constantly if different FMAs are done
        ! on the same DUMP file

        if(dumplast(j) == -1) then
          dumpLastTurn = numl
        else
          dumpLastTurn = dumplast(j)
        end if

        ! Loop over all turns in the DUMP file; this is neccessary since we're writing normalised DUMP files.
        do k=dumpfirst(j),dumpLastTurn ! Loop over turns, use the dump files
          ! Loop over particles
          do l=1,napx
            if(dumpfmt(j) == 2 .or. dumpfmt(j) == 7) then ! Read an ASCII dump
              read(dumpunit(j),"(a)",iostat=ierro) rLine
              call chr_split(rLine, lnSplit, nSplit, spErr)
              if(spErr) then
                write(lerr,"(a,i0,a)") "FMA> ERROR Failed to parse line from file '"//trim(dump_fname(j))//&
                  "' (dumpfmt = ",dumpfmt(j),")"
                call prror
              end if
              if(nSplit > 0) call chr_cast(lnSplit(1), id,          cErr)
              if(nSplit > 1) call chr_cast(lnSplit(2), thisturn,    cErr)
              if(nSplit > 2) call chr_cast(lnSplit(3), pos,         cErr)
              if(nSplit > 3) call chr_cast(lnSplit(4), xyzvdummy(1),cErr)
              if(nSplit > 4) call chr_cast(lnSplit(5), xyzvdummy(2),cErr)
              if(nSplit > 5) call chr_cast(lnSplit(6), xyzvdummy(3),cErr)
              if(nSplit > 6) call chr_cast(lnSplit(7), xyzvdummy(4),cErr)
              if(nSplit > 7) call chr_cast(lnSplit(8), xyzvdummy(5),cErr)
              if(nSplit > 8) call chr_cast(lnSplit(9), xyzvdummy(6),cErr)
              if(nSplit > 9) call chr_cast(lnSplit(10),kt,          cErr)
            else if(dumpfmt(j) == 3 .or. dumpfmt(j) == 8) then ! Read a binary dump
              read(dumpunit(j),iostat=ierro) &
                id,thisturn,pos,xyzvdummy(1),xyzvdummy(2),xyzvdummy(3),xyzvdummy(4),xyzvdummy(5),xyzvdummy(6),kt
              if(ierro /= 0) then
                write(lerr,"(a,i0,a)") "FMA> ERROR Failed to parse line from file '"//trim(dump_fname(j))//&
                  "' (dumpfmt = ",dumpfmt(j),")"
                call prror
              end if
            end if

            ! Check for losses
            if(l /= id .or. k /= thisturn) then
              if(k < nturns(l)+fma_first(l)-1) then
                nturns(l) = k-fma_first(l)
              end if

              !TODO: Actually handle those losses.
              write(lerr,"(2(a,i0))") "FMA> ERROR Reading DUMP file #",j, " for FMA #",i
              write(lerr,"(2(a,i0))") "FMA>       Expected turn ",k," and particle ID ",l
              write(lerr,"(2(a,i0))") "FMA>       Got turn ",thisturn," and particle ID ",id
              write(lerr,"(a)")       "FMA>       Reading probably got unsynchronized because of particle losses,"//&
                " which is currently not handled in FMA."
              call prror
            end if

            ! Normalization
            if(dumpfmt(j) == 2 .or. dumpfmt(j) == 3) then
              ! Case: The file isn't pre-normalised -> Compute normalization
              !
              ! At this point fma_norm_flag doesn't matter;
              ! we anyway compute the normalised coordinates.
              !
              ! units: dumptas, dumptasinv, dumpclo [mm,mrad,mm,mrad,1]

              ! remove closed orbit -> check units used in dumpclo (is x' or px used?)
              do m=1,6
                xyzvdummy2(m)=xyzvdummy(m)-dumpclo(m)
              end do

              ! For use in with normalised coordinates: convert to canonical variables
              xyzvdummy2(2)=xyzvdummy2(2) * ((one+xyzvdummy2(6))+dumpclo(6))
              xyzvdummy2(4)=xyzvdummy2(4) * ((one+xyzvdummy2(6))+dumpclo(6))

              ! Normalise nxyz=dumptasinv*xyz2
              do m=1,6
                nxyzvdummy(m)=zero
                do n=1,6
                  nxyzvdummy(m)=nxyzvdummy(m) + dumptasinv(m,n)*xyzvdummy2(n)
                end do
                ! convert nxyzvdummy(6) to 1.e-3 sqrt(m)
                ! unit: nx,npx,ny,npy,nsig,ndelta all in [1.e-3 sqrt(m)]
                if(m == 6) then
                  nxyzvdummy(m)=nxyzvdummy(m)*c1e3
                end if
              end do

              ! Write normalised particle amplitudes
              ! (only when reading physical coordinates)
              if(fma_writeNormDUMP .and. .not.(dumpfmt(j) == 7 .or. dumpfmt(j) == 8) .and. .not.hasNormDumped(j) ) then
                write(tmpUnit,"(2(1x,i8),1x,f12.5,6(1x,1pe16.9),1x,i8)") &
                  id,thisturn,pos,nxyzvdummy(1),nxyzvdummy(2),nxyzvdummy(3),nxyzvdummy(4),nxyzvdummy(5),nxyzvdummy(6),kt
              end if

            else if(dumpfmt(j) == 7 .or. dumpfmt(j) == 8) then
              ! Case: we are already normalised;
              ! just copy the data into the relevant array
              do m=1,6
                nxyzvdummy(m) = xyzvdummy(m)
              end do
            end if ! END IF already normalised or not

            ! Copy the data into the final arrays
            if(thisturn >= fma_first(i) .and. thisturn <= fma_last(i) ) then
              turn(l,k-fma_first(i)+1) = thisturn
              do m=1,6
                ! For FMA in physical coordinates, convert units to [mm,mrad,mm,mrad,mm,1.e-3]
                if(m == 6) then
                  xyzv(l,k-fma_first(i)+1,m)=xyzvdummy(m)*c1e3
                else
                  xyzv(l,k-fma_first(i)+1,m)=xyzvdummy(m)
                end if

                nxyzv(l,k-fma_first(i)+1,m) = nxyzvdummy(m)
                ! calculate emittance of mode 1,2,3
                if(mod(m,2) == 0) then
                  epsnxyzv(l,k-fma_first(i)+1,m/2) = nxyzvdummy((m-1))**2+nxyzvdummy(m)**2
                end if
              end do
            end if ! END if fma_first <= thisturn <= fma_last

          end do ! END loop over particles l
        end do ! END loop over turns k

        ! Calculate tunes of particles using the methods in plato_seq.f and NAFF
        !  for fma_norm_flag == 0: use physical coordinates x,x',y,y',sig,dp/p
        !  for fma_norm_flag == 1: use normalised coordinates
        do l=1,napx ! loop over particles
          !TODO particle losses - detect if nturns(l) is too small & skip that particle.
          ! (probably just write a line of mostly zeros to the file)

          do m=1,numModes ! Loop over modes (hor.,vert.,long.)
            select case(trim(fma_method(i)))
            case("TUNELASK")
              if(fma_norm_flag(i) == 0) then
                q123(m) =  tunelask( xyzv(l,1:nturns(l),2*(m-1)+1), xyzv(l,1:nturns(l),2*m), nturns(l) )
              else
                q123(m) =  tunelask(nxyzv(l,1:nturns(l),2*(m-1)+1),nxyzv(l,1:nturns(l),2*m), nturns(l) )
              end if

            case("TUNEFFTI")
              if(fma_norm_flag(i) == 0) then
                q123(m) =  tuneffti( xyzv(l,1:nturns(l),2*(m-1)+1), xyzv(l,1:nturns(l),2*m), nturns(l) )
              else
                q123(m) =  tuneffti(nxyzv(l,1:nturns(l),2*(m-1)+1),nxyzv(l,1:nturns(l),2*m), nturns(l) )
              end if

            case("TUNEFFT")
              if(fma_norm_flag(i) == 0) then
                q123(m) =   tunefft( xyzv(l,1:nturns(l),2*(m-1)+1), xyzv(l,1:nturns(l),2*m), nturns(l) )
              else
                q123(m) =   tunefft(nxyzv(l,1:nturns(l),2*(m-1)+1),nxyzv(l,1:nturns(l),2*m), nturns(l) )
              end if

            case("TUNEAPA")
              if(fma_norm_flag(i) == 0) then
                q123(m) =   tuneapa( xyzv(l,1:nturns(l),2*(m-1)+1), xyzv(l,1:nturns(l),2*m), nturns(l) )
              else
                q123(m) =   tuneapa(nxyzv(l,1:nturns(l),2*(m-1)+1),nxyzv(l,1:nturns(l),2*m), nturns(l) )
              end if

            case("TUNEFIT")
              if(fma_norm_flag(i) == 0) then
                q123(m) =   tunefit( xyzv(l,1:nturns(l),2*(m-1)+1), xyzv(l,1:nturns(l),2*m), nturns(l) )
              else
                q123(m) =   tunefit(nxyzv(l,1:nturns(l),2*(m-1)+1),nxyzv(l,1:nturns(l),2*m), nturns(l) )
              end if

            case("TUNENEWT")
              if(fma_norm_flag(i) == 0) then
                q123(m) =  tunenewt( xyzv(l,1:nturns(l),2*(m-1)+1), xyzv(l,1:nturns(l),2*m), nturns(l) )
              else
                q123(m) =  tunenewt(nxyzv(l,1:nturns(l),2*(m-1)+1),nxyzv(l,1:nturns(l),2*m), nturns(l) )
              end if

            case("TUNEABT2")
              if(fma_norm_flag(i) == 0) then
                q123(m) =  tuneabt2( xyzv(l,1:nturns(l),2*(m-1)+1), xyzv(l,1:nturns(l),2*m), nturns(l) )
              else
                q123(m) =  tuneabt2(nxyzv(l,1:nturns(l),2*(m-1)+1),nxyzv(l,1:nturns(l),2*m), nturns(l) )
              end if

            case("TUNEABT")
              if(fma_norm_flag(i) == 0) then
                q123(m) =   tuneabt( xyzv(l,1:nturns(l),2*(m-1)+1), xyzv(l,1:nturns(l),2*m), nturns(l) )
              else
                q123(m) =   tuneabt(nxyzv(l,1:nturns(l),2*(m-1)+1),nxyzv(l,1:nturns(l),2*m), nturns(l) )
              end if

            case("TUNENEWT1")
              if(fma_norm_flag(i) == 0) then
                q123(m) = tunenewt1( xyzv(l,1:nturns(l),2*(m-1)+1), xyzv(l,1:nturns(l),2*m), nturns(l) )
              else
                q123(m) = tunenewt1(nxyzv(l,1:nturns(l),2*(m-1)+1),nxyzv(l,1:nturns(l),2*m), nturns(l) )
              endif

#ifdef NAFF
            case("NAFF")
              flush(lout)
              ! Copy the relevant contents of the arrays
              ! into a new temporary array with stride=1
              ! for passing to C++.
              if(fma_norm_flag(i) == 0) then
                naff_xyzv1 = xyzv (l, 1:nturns(l), 2*(m-1)+1)
                naff_xyzv2 = xyzv (l, 1:nturns(l), 2*m)
              else
                naff_xyzv1 = nxyzv(l,1:nturns(l),  2*(m-1)+1)
                naff_xyzv2 = nxyzv(l,1:nturns(l),  2*m)
              endif
              fft_naff = tunefft(naff_xyzv1, naff_xyzv2, nturns(l))
#ifndef DOUBLE_MATH
              q123(m)  = real(tunenaff(real(naff_xyzv1,kind=real64), real(naff_xyzv2,kind=real64), &
                nturns(l), m, fma_norm_flag(i), real(fft_naff,kind=real64)),kind=fPrec)
#else
              q123(m)  = tunenaff(naff_xyzv1, naff_xyzv2, nturns(l), m, fma_norm_flag(i), fft_naff)
#endif

              flush(lout)
#endif

            case default
              write(lerr,"(a)") "FMA> ERROR Method '"//trim(fma_method(i))//&
                "' not known. Note that the method name must be in capital letters."
              call prror
            end select

            ! mode 3 rotates anticlockwise, mode 1 and 2 rotate clockwise -> synchroton tune is negative,
            ! but define it as convention positive
            if(m == 3) q123(m)=one-q123(m)

            ! Some general calculations
            eps123_0(m)   = epsnxyzv(l,1,m)                               ! initial amplitude
            phi123_0(m)   = atan2_mb(nxyzv(l,1,2*m),nxyzv(l,1,2*(m-1)+1)) ! inital phase
            eps123_min(m) = minval( epsnxyzv(l,1:nturns(l),m) )           ! minimum emittance
            eps123_max(m) = maxval( epsnxyzv(l,1:nturns(l),m) )           ! maximum emittance
            eps123_avg(m) = sum(epsnxyzv(l,1:nturns(l),m))/nturns(l)      ! average emittance

          end do ! Loop over modes (hor.,vert.,long.)

          if(numModes == 2) then
            q123(3)       = zero
            eps123_min(3) = zero
            eps123_max(3) = zero
            eps123_avg(3) = zero
            eps123_0(3)   = zero
            phi123_0(3)   = zero
          end if

          ! Write the FMA output file "fma_sixtrack"
          ! TODO losses: fma_first and fma_last may not be the right start/stop variables...
          write(fmaUnit,"(2(1x,a20),1x,i8,18(1x,1pe16.9),3(1x,i8))") &
            trim(fma_fname(i)), trim(fma_method(i)),l,q123(1),q123(2),q123(3), &
            eps123_min(1),eps123_min(2),eps123_min(3),eps123_max(1), &
            eps123_max(2),eps123_max(3),eps123_avg(1),eps123_avg(2), &
            eps123_avg(3),eps123_0(1),eps123_0(2),eps123_0(3), &
            phi123_0(1),phi123_0(2),phi123_0(3),fma_norm_flag(i), &
            fma_first(i),fma_last(i)

        end do ! END loop over particles l

        if(fma_writeNormDUMP .and. .not.(dumpfmt(j) == 7 .or. dumpfmt(j) == 8) .and. .not.hasNormDumped(j)) then
          ! filename NORM_* (normalised particle amplitudes)
          call f_close(tmpUnit)
          call f_freeUnit(tmpUnit)
          hasNormDumped(j) = .true.
        end if

        ! resume initial position of dumpfile = end of file
        call f_close(dumpunit(j))
        if(dumpfmt(j) == 2 .or. dumpfmt(j) == 7) then ! ASCII
          call f_open(unit=dumpunit(j),file=dump_fname(j),formatted=.true.,mode="rw+",err=fErr)
          if(fErr) then
            write(lerr,"(a,i0,a)") "FMA> ERROR Resuming file '"//trim(dump_fname(j))//"' (dumpfmt = ",dumpfmt(j),")"
            call prror
          end if
        elseif (dumpfmt(j).eq.3 .or. dumpfmt(j).eq.8) then !BINARY
          call f_open(unit=dumpunit(j),file=dump_fname(j),formatted=.false.,mode="rw+",err=fErr)
          if(fErr) then
            write(lerr,"(a,i0,a)") "FMA> ERROR Resuming file '"//trim(dump_fname(j))//"' (dumpfmt = ",dumpfmt(j),")"
            call prror
          end if
        end if
      end if ! END: if fma_fname(i) matches dump_fname(j)

      ! If file has been already found, jump to next file fma_fname(i)
      if(fExist) exit
    end do !END: loop over dump files

    if(.not. fExist) then ! if no dumpfile has been found, raise error and abort
      write(lerr,"(a)") "FMA> ERROR DUMP file '"//trim(fma_fname(i))//&
        "' does not exist. Please check that filenames in FMA block agree with the ones in the DUMP block."
      call prror
    end if

  end do ! END: loop over fma files

  call f_close(fmaUnit) !filename: fma_sixtrack

  call dealloc(turn,         "turn")
  call dealloc(nturns,       "nturns")
  call dealloc(hasNormDumped,"hasNormDumped")
  call dealloc(fma_nturn,    "fma_nturn")
  call dealloc(xyzv,         "xyzv")
  call dealloc(nxyzv,        "nxyzv")
  call dealloc(epsnxyzv,     "epsnxyzv")
  call dealloc(naff_xyzv1,   "naff_xyzv1")
  call dealloc(naff_xyzv2,   "naff_xyzv2")

end subroutine fma_postpr

end module fma
