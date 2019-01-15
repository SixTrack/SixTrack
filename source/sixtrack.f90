! ================================================================================================ !
!  DATEN - INPUT PARSING
! ~~~~~~~~~~~~~~~~~~~~~~~
!  Rewritten by V.K. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-06-25
!  Reads input data from files fort.2 and fort.3
! ================================================================================================ !
subroutine daten

  use crcoall
  use floatPrecision
  use sixtrack_input
  use parpro
  use parbeam, only : beam_expflag
  use mod_settings
  use mod_common
  use mod_commons
  use mod_commont
  use mod_commond
  use physical_constants
  use numerical_constants
  use string_tools
  use mod_alloc
  use mod_units

  use mod_dist,  only : dist_enable, dist_parseInputLine
  use scatter,   only : scatter_active,scatter_debug,scatter_dumpdata,scatter_parseInputLine,scatter_allocate
  use dynk,      only : dynk_enabled,dynk_debug,dynk_dumpdata,dynk_inputsanitycheck,dynk_allocate,dynk_parseInputLine
  use fma,       only : fma_parseInputLine, fma_allocate
  use dump,      only : dump_parseInputLine,dump_parseInputDone
  use zipf,      only : zipf_parseInputLine,zipf_parseInputDone
  use bdex,      only : bdex_parseInputLine,bdex_parseInputDone
  use mod_fluc,  only : fluc_parseInputLine,fluc_readInputs
  use wire,      only : wire_parseInputLine,wire_parseInputDone
  use elens,     only : elens_parseInputLine,elens_parseInputDone,elens_postInput
  use aperture
  use mod_hions
#ifdef HASHLIB
  use mod_hash
#endif
#ifdef FLUKA
  use mod_fluka, only : fluka_parsingDone,fluka_parseInputLine,fluka_enable
#endif
#ifdef HDF5
  use hdf5_output
#endif
#ifdef ROOT
  use root_output
#endif
  use collimation

  implicit none

  character(len=mInputLn) inLine, pLines(5)
  character(len=mNameLen) ic0(10)
  character(len=60)       iHead
  character(len=8)        cPad
  character(len=4)        currBlock, cCheck

  integer nUnit,lineNo2,lineNo3,nGeom
  integer blockLine,blockCount

  logical blockOpened,blockClosed,blockReopen,openBlock,closeBlock
  logical inErr,fErr,parseFort2

  integer icc,il1,ilin0,iMod,i,j,k,k10,k11,kk,l,ll,l1,l2,l3,l4,mblozz,nac,nfb,nft

! ================================================================================================ !
!  SET DEFAULT VALUES
! ================================================================================================ !

  ! Main Variables
  iHead        = " "
  sixtit       = " "
  ic0(:)       = " "
  cCheck       = " "
  kanf         = 1

  ! TRACKING PARAMETERS
  numl         = 1
  napx         = 0
  amp0         = zero
  amp(1)       = c1m3
  ird          = 0
  imc          = 0
  numlcp       = 1000
  numlmax      = 1000000000
  iclo6        = 0
  iclo6r       = 0

  idz(:)       = 1
  idfor        = 0
  irew         = 0
  iclo         = 0

  nde(:)       = 0
  nwr(:)       = 1
  nwr(4)       = 10000
  ntwin        = 1

  ! INITIAL COORDINATES
  itra         = 0
  chi0         = zero
  chid         = zero
  rat          = zero

  ! CHROMATICITY ADJUSTMENTS
  ichrom       = 0

  ! TUNE ADJUSTMENTS
  iqmod        = 0

  ! LINEAR OPTICS CALCULATION
  ilin         = 0
  sixin_ilin0  = 1
  nlin         = 0

  ! DIFFERENTIAL ALGEBRA
  preda        = c1m38

  ! SYNCHROTRON OSCILLATIONS
  sixin_alc    = c1m3
  sixin_harm   = one
  sixin_phag   = zero
  idp          = 0
  ncy          = 0
  tlen         = one
  pma          = pmap
  ition        = 0
  dpscor       = one
  sigcor       = one
  qs           = zero

  ! MULTIPOLE COEFFICIENTS
  sixin_im     = 0

  ! RF - MULTIPOLE
  sixin_rfm    = 0

  ! RANDOM FLUCTUATIONS
  izu0         = 0
  mcut         = 0
  mout2        = 0

  ! SUB-RESONANCE CALCULATION
  isub         = 0

  ! ORGANISATION OF RANDOM NUMBERS
  iorg         = 0

  ! RESONANCE COMPENSATION
  irmod2       = 0

  ! DECOUPLING OF MOTION
  iskew        = 0

  ! NORMAL FORMS
  idial        = 0

  ! SEARCH
  ise          = 0

  ! ITERATION ERRORS
  itco         = 500
  dma          = c1m12
  dmap         = c1m15
  itcro        = 10
  dech         = c1m10
  de0          = c1m9
  ded          = c1m9
  dsi          = c1m9
  dsm0         = c1m10
  itqv         = 10
  dkq          = c1m10
  dqq          = c1m10

  ! BEAM-BEAM ELEMENT
  sixin_emitNX = zero
  sixin_emitNY = zero

  ! PHASE TROMBONE
  sixin_imtr0  = 0
  ntr          = 1
  call alloc(cotr,1,6,  zero,"cotr")
  call alloc(rrtr,1,6,6,zero,"rrtr")

  ! POST PROCESSING
  ipos         = 0
  iconv        = 0
  imad         = 0
  iskip        = 1
  cma1         = one
  cma2         = one
  iav          = 1
  iwg          = 1
  dphix        = zero
  dphiz        = zero
  qx0          = zero
  qz0          = zero
  ivox         = 1
  ivoz         = 1
  ires         = 1
  dres         = one
  ifh          = 0
  dfft         = one
  idis         = 0
  icow         = 0
  istw         = 0
  iffw         = 0
  nprint       = 1
  ndafi        = 1

  ! HIONS MODULE
  zz0          = 1
  aa0          = 1
  nucm0        = pmap
  has_hion     = .false.

! ================================================================================================ !
!  SET DEFAULT INPUT PARSING VALUES
! ================================================================================================ !

  ! SIXTRACK INPUT MODULE
  inErr       = .false.
  sixin_ncy2  = 0
  sixin_icy   = 0

  call alloc(sixin_bez0,mNameLen,nele,str_nmSpace,"sixin_bez0")

  ! DATEN INTERNAL
  nGeom       = 0
  lineNo2     = 0
  lineNo3     = 0
  pLines(:)   = " "

! ================================================================================================ !
!  READ FORT.3 HEADER
! ================================================================================================ !

  call f_open(unit=3,file="fort.3",formatted=.true.,mode="r",err=fErr)
  if(fErr) then
    write(lout,"(a)") "INPUT> ERROR Could not open fort.3"
    call prror
  end if

90 continue
  read(3,"(a4,a8,a60)",end=9997,iostat=ierro) cCheck,cPad,iHead
  if(ierro > 0) then
    write(lout,"(a)") "INPUT> ERROR Could not read from fort.3"
    call prror
  end if
  pLines(5) = cCheck//cPad//iHead
  lineNo3 = lineNo3+1
  if(cCheck(1:1) == "/") goto 90
  if(cCheck(1:1) == "!") goto 90

  select case(cCheck)
  case("FREE") ! Mode FREE. Elements in fort.3
    iMod       = 1
    parseFort2 = .false.
  case("GEOM") ! Mode GEOM. Elements in fort.2
    iMod       = 2
    parseFort2 = .true.
    call f_open(unit=2,file="fort.2",formatted=.true.,mode="r",err=fErr)
    if(fErr) then
      write(lout,"(a)") "INPUT> ERROR Could not open fort.2"
      call prror
    end if
  case default
    write(lout,"(a)") "INPUT> ERROR Unknown mode '"//cCheck//"'"
    goto 9999
  end select

  write(lout,"(a)") ""
  write(lout,"(a)") "    OOOOOOOOOOOOOOOOOOOOOO"
  write(lout,"(a)") "    OO                  OO"
  write(lout,"(a)") "    OO  SixTrack Input  OO"
  write(lout,"(a)") "    OO                  OO"
  write(lout,"(a)") "    OOOOOOOOOOOOOOOOOOOOOO"
  write(lout,"(a)") ""
  if(ihead /= " ") write(lout,"(a)") "    TITLE: "//trim(iHead)
  if(imod  == 1)   write(lout,"(a)") "    MODE:  Free Format Input (fort.3)"
  if(imod  == 2)   write(lout,"(a)") "    MODE:  Geometry Strength File (fort.2)"
  write(lout,"(a)") ""
  write(lout,"(a)") str_divLine

  sixtit(1:60) = iHead

! ================================================================================================ !
!  BEGIN PARSING FORT.2 AND FORT.3
! ================================================================================================ !

  currBlock   = "NONE"  ! The current block being parsed
  openBlock   = .false. ! Whether the current block was opened on this pass
  blockOpened = .false. ! Whether the current block is open
  closeBlock  = .false. ! Whether the current block should be closed after the pass
  blockClosed = .false. ! Whether the current block is now closed, and should not be opened again

110 continue ! fort.3 loop

  ! We have our three geometry blocks, stop parsing fort.2
  if(nGeom >= 3) then
    parseFort2 = .false.
  end if

  ! Select unit, and increment line number for error output
  if(parseFort2) then
    nUnit   = 2
    lineNo2 = lineNo2 + 1
  else
    nUnit   = 3
    lineNo3 = lineNo3 + 1
  end if

  read(nUnit,"(a)",end=9998,iostat=iErro) inLine
  if(iErro > 0) then
    write(lout,"(a,i0)") "INPUT> ERROR Could not read from fort.",nUnit
    call prror
  end if

  ! Keep the last few lines for error output, but only fort fort.3
  if(nUnit == 3) then
    do i=1,4
      pLines(i) = pLines(i+1)
    end do
    pLines(5) = inLine
  end if

  if(len_trim(inLine) == 0) goto 110 ! Empty line, ignore
  if(inLine(1:1) == "/")    goto 110 ! Comment line, ignore
  if(inLine(1:1) == "!")    goto 110 ! Comment line, ignore
  read(inLine,"(a4)") cCheck

  ! Check for end of block flag
  if(cCheck == "NEXT") then
    if(currBlock == "NONE") then
      ! Catch orphaned NEXT blocks here.
      write(lout,"(a)") "INPUT> ERROR Unexpected NEXT block encountered. There is no open block to close."
      goto 9999
    else
      ! Actual close check is done after a last pass so
      ! each block can finalise any necessary initialisation
      closeBlock = .true.
    end if
  end if

  ! Check for end of fort.3 input
  if(cCheck == "ENDE") then
    if(iMod == 2) call f_close(2)
    call f_close(3)
    goto 9000
  end if

  ! Check if no block is active. If so, there should be a new one if input is sane.
  if(currBlock == "NONE") then
    currBlock = cCheck
    openBlock = .true.
  end if

  ! By default blocks can only be opened once. The exceptions are set here.
  blockReopen = .false.
  if(currBlock == "LINE") blockReopen = .true.
  if(currBlock == "MULT") blockReopen = .true.
  if(currBlock == "TROM") blockReopen = .true.
  if(currBlock == "RFMU") blockReopen = .true.

  ! Check the status of the current block
  call sixin_checkBlock(currBlock, nUnit, blockOpened, blockClosed, blockLine, blockCount)

  ! Check if the current block has already been seen and closed.
  ! If so, the block exists more than once in the input files. It shouldn't unless intended to.
  if(blockCount > 1 .and. .not. blockReopen) then
    write(lout,"(a)") "INPUT> ERROR Block '"//currBlock//"' encountered more than once."
    goto 9999
  end if

  ! Block-Wise Parsing Code
  !=========================
  ! These are parsed one line at a time, including the block declaration and the NEXT.
  ! When the block is first encountered, the 'openBlock' flag is .true. for one pass.
  ! When the block is closed, the 'closeBlock' flag is .true. for the final pass.
  ! Otherwise, we have an input line that should be parsed by a single subroutine call.
  ! The in-block line number is available as the variable 'blockLine' if needed.

  select case(currBlock)

  case("PRIN") ! Enable the PRINT flag
    if(openBlock) then
      st_print = .true.
      write(lout,"(a)") "INPUT> Printout of input parameters ENABLED"
      write(lout,"(a)") "INPUT> WARNING The PRINT block is deprectaed and will be removed in a future release. Please use:"
      write(lout,"(a)") ""
      write(lout,"(a)") "SETTINGS"
      write(lout,"(a)") "  PRINT"
      write(lout,"(a)") "NEXT"
      write(lout,"(a)") ""
    elseif(closeBlock) then
      continue
    else
      write(lout,"(a)") "INPUT> ERROR PRINT block does not take any parameters. Did you forget to close it with a NEXT?"
      goto 9999
    end if

  case("SETT") ! Global Settings Block
    if(openBlock) then
      continue
    elseif(closeBlock) then
      continue
    else
      call sixin_parseInputLineSETT(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("COMM") ! Comment Block
    if(openBlock) then
      continue
    elseif(closeBlock) then
      continue
    else
      if(len(inLine) > 80) then
        commen = inLine(1:80)
      else
        commen = inLine
      end if
    end if

  case("SING") ! Single Elements Block
    if(openBlock) then
      sixin_nSing = 1
    elseif(closeBlock) then
      nGeom = nGeom + 1
    else
      call sixin_parseInputLineSING(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("BLOC") ! Block Definitions
    if(openBlock) then
      sixin_nBloc = 0
    elseif(closeBlock) then
      nGeom = nGeom + 1
    else
      call sixin_parseInputLineBLOC(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("STRU") ! Structure Input
    if(openBlock) then
      sixin_nStru = 0
    elseif(closeBlock) then
      nGeom = nGeom + 1
    else
      call sixin_parseInputLineSTRU(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("DISP") ! Displacement of Elements Block
    if(openBlock) then
      continue
    elseif(closeBlock) then
      continue
    else
      call sixin_parseInputLineDISP(inLine,inErr)
      if(inErr) goto 9999
    end if

  case("INIT") ! Initial Coordinates
    if(openBlock) then
      continue
    elseif(closeBlock) then
      dp1 = exz(1,6)
    else
      call sixin_parseInputLineINIT(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("TRAC") ! Tracking Parameters
    if(openBlock) then
      continue
    elseif(closeBlock) then
      continue
    else
      call sixin_parseInputLineTRAC(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("DIFF") ! Differential Algebra
    if(openBlock) then
      continue
    elseif(closeBlock) then
      continue
    else
      call sixin_parseInputLineDIFF(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("CHRO") ! Chromaticity Adjustment
    if(openBlock) then
      ichrom = 1
    elseif(closeBlock) then
      continue
    else
      call sixin_parseInputLineCHRO(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("TUNE") ! Tune Adjustment
    if(openBlock) then
      iqmod = 1
    elseif(closeBlock) then
      call sixin_parseInputLineTUNE(inLine,-1,inErr)
    else
      call sixin_parseInputLineTUNE(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("LINE") ! Linear Optics Calculation Block
    if(openBlock) then
      ilin0 = 1
    elseif(closeBlock) then
      continue
    else
      call sixin_parseInputLineLINE(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("SYNC") ! Synchrotron Oscillations
    if(openBlock) then
      continue
    elseif(closeBlock) then
      continue
    else
      call sixin_parseInputLineSYNC(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("MULT") ! Multipole Coefficients
    if(openBlock) then
      continue
    elseif(closeBlock) then
      continue
    else
      call sixin_parseInputLineMULT(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("RFMU") ! RF - Multipole
    if(openBlock) then
      continue
    elseif(closeBlock) then
      continue
    else
      call sixin_parseInputLineRFMU(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("FLUC") ! Fluctuation Random Starting Number
    if(openBlock) then
      continue
    elseif(closeBlock) then
      call fluc_readInputs
    else
      call fluc_parseInputLine(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("SUBR") ! Sub-Resonance Calculation
    if(openBlock) then
      write(lout,"(a)") "SUBR> WARNING This block is inhertited from older versions of SixTrack and is not covered by tests."
      write(lout,"(a)") "SUBR>         It therefore may not produce the rosults expected."
      write(lout,"(a)") "SUBR>         Please report any bugs to the dev team."
    elseif(closeBlock) then
      isub = 1
    else
      call sixin_parseInputLineSUBR(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("ORGA") ! Organisation of Random Numbers
    if(openBlock) then
      write(lout,"(a)") "ORGA> WARNING This block is inhertited from older versions of SixTrack and is not covered by tests."
      write(lout,"(a)") "ORGA>         It therefore may not produce the rosults expected."
      write(lout,"(a)") "ORGA>         Please report any bugs to the dev team."
    elseif(closeBlock) then
      continue
    else
      call sixin_parseInputLineORGA(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("ITER") ! Iteration Errors
    if(openBlock) then
      continue
    elseif(closeBlock) then
      continue
    else
      call sixin_parseInputLineITER(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("ORBI") ! Orbit Correction
    if(openBlock) then
      write(lout,"(a)") "ORBI> WARNING This block is inhertited from older versions of SixTrack and is not covered by tests."
      write(lout,"(a)") "ORBI>         It therefore may not produce the rosults expected."
      write(lout,"(a)") "ORBI>         Please report any bugs to the dev team."
    elseif(closeBlock) then
      continue
    else
      call sixin_parseInputLineORBI(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("COMB") ! Combination of Elements
    if(openBlock) then
      write(lout,"(a)") "COMB> WARNING This block is inhertited from older versions of SixTrack and is not covered by tests."
      write(lout,"(a)") "COMB>         It therefore may not produce the rosults expected."
      write(lout,"(a)") "COMB>         Please report any bugs to the dev team."
    elseif(closeBlock) then
      continue
    else
      call sixin_parseInputLineCOMB(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("RESO") ! Resonance Compensation
    if(openBlock) then
      write(lout,"(a)") "RESO> WARNING This block is inhertited from older versions of SixTrack and is not covered by tests."
      write(lout,"(a)") "RESO>         It therefore may not produce the rosults expected."
      write(lout,"(a)") "RESO>         Please report any bugs to the dev team."
    elseif(closeBlock) then
      continue
    else
      call sixin_parseInputLineRESO(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("SEAR") ! Search for Optimum Places to Compensate Resonances
    if(openBlock) then
      write(lout,"(a)") "SEAR> WARNING This block is inhertited from older versions of SixTrack and is not covered by tests."
      write(lout,"(a)") "SEAR>         It therefore may not produce the rosults expected."
      write(lout,"(a)") "SEAR>         Please report any bugs to the dev team."
    elseif(closeBlock) then
      continue
    else
      call sixin_parseInputLineSEAR(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("DECO") ! Decoupling of Motion in the Transverse Planes
    if(openBlock) then
      write(lout,"(a)") "DECO> WARNING This block is inhertited from older versions of SixTrack and is not covered by tests."
      write(lout,"(a)") "DECO>         It therefore may not produce the rosults expected."
      write(lout,"(a)") "DECO>         Please report any bugs to the dev team."
    elseif(closeBlock) then
      continue
    else
      call sixin_parseInputLineDECO(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("NORM") ! Normal Forms
    if(openBlock) then
      write(lout,"(a)") "NORM> WARNING This block is inhertited from older versions of SixTrack and is not covered by tests."
      write(lout,"(a)") "NORM>         It therefore may not produce the rosults expected."
      write(lout,"(a)") "NORM>         Please report any bugs to the dev team."
    elseif(closeBlock) then
      continue
    else
      call sixin_parseInputLineNORM(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("POST") ! Post-Processing
    if(openBlock) then
      continue
    elseif(closeBlock) then
      continue
    else
      call sixin_parseInputLinePOST(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("BEAM") ! Beam-Beam Element
    if(openBlock) then
      continue
    elseif(closeBlock) then
      continue
    else
      if(beam_expflag == 1) then
        call sixin_parseInputLineBEAM_EXP(inLine,blockLine,inErr)
      else
        call sixin_parseInputLineBEAM(inLine,blockLine,inErr)
      end if
      if(inErr) goto 9999
    end if

  case("TROM") ! “Phase Trombone” Element
    if(openBlock) then
      continue
    elseif(closeBlock) then
      call sixin_parseInputLineTROM(inLine,-1,inErr)
    else
      call sixin_parseInputLineTROM(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("FLUK") ! Fluka Coupling
#ifndef FLUKA
    if(openBlock) then
      write(lout,"(a)") "INPUT> ERROR SixTrack was not compiled with the FLUKA flag."
      goto 9999
    else
      continue
    end if
#else
    if(openBlock) then
      continue
    elseif(closeBlock) then
      call fluka_parsingDone
    else
      call fluka_parseInputLine(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if
#endif

  case("BDEX") ! Beam Distribution EXchange
    if(openBlock) then
      continue
    elseif(closeBlock) then
      call bdex_parseInputDone
    else
      call bdex_parseInputLine(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("WIRE") ! Wire
    if(openBlock) then
      continue
    elseif(closeBlock) then
      call wire_parseInputDone(inErr)
      if(inErr) goto 9999
    else
      call wire_parseInputLine(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("ELEN") ! Electron Lens
    if(openBlock) then
      continue
    elseif(closeBlock) then
      call elens_parseInputDone(inErr)
      if(inErr) goto 9999
    else
      call elens_parseInputLine(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("DIST") ! Beam Distribution
    if(openBlock) then
      dist_enable = .true.
    elseif(closeBlock) then
      continue
    else
      call dist_parseInputLine(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("CORR") ! Tuneshift Corrections
    write(lout,"(a)") "INPUT> ERROR CORR module is deprecated."
    goto 9999

  case("RIPP") ! Power Supply Ripple Block
    write(lout,"(a)") "INPUT> ERROR RIPP module is deprecated and replaced by DYNK."
    write(lout,"(a)") "INPUT>       The script rippconvert.py in the pytools folder can be used to convert the fort.3 file."
    goto 9999

  case("LIMI") ! Aperture Limitations
    if(openBlock) then
      lbacktracking = .true.
    elseif(closeBlock) then
      call aper_inputParsingDone
    else
      call aper_parseInputLine(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("COLL") ! Collimation Block
#ifdef CR
    write(lout,'(a)') "INPUT> ERROR Collimation incompatible with checkpoint/restart (CR)"
    goto 9999
#endif
    if(openBlock) then
      ! If a collimation block is present, even disabled, allocate the storage.
      ! This mimmics the old compiler flag.
      if(ilin /= 1) then
        write(lout,"(a)") "INPUT> ERROR Incompatible flag with collimation version detected in the LINEAR OPTICS block."
        write(lout,"(a)") "INPUT>       You have not chosen ilin=1 (4D mode), which is required for the collimation version."
        write(lout,"(a)") "INPUT>       Note that the ilin=2 (6D mode) is not compatible with the collimation version."
        goto 9999
      end if
    elseif(closeBlock) then
      continue
    else
      call collimate_parseInputLine(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("FMA") ! FMA Block
    if(openBlock) then
      call fma_allocate
    elseif(closeBlock) then
      continue
    else
      call fma_parseInputLine(inLine,inErr)
      if(inErr) goto 9999
    end if

  case("DUMP") ! DUMP Block
    if(openBlock) then
      continue
    elseif(closeBlock) then
      call dump_parseInputDone(inErr)
    else
      call dump_parseInputLine(inLine,inErr)
      if(inErr) goto 9999
    end if

  case("DYNK") ! Dynamic Kick Block
    if(openBlock) then
      call dynk_allocate
    elseif(closeBlock) then
      if(dynk_debug) call dynk_dumpdata
      dynk_enabled = .true.
      call dynk_inputsanitycheck
    else
      call dynk_parseInputLine(inLine,inErr)
      if(inErr) goto 9999
    end if

  case("HION") ! Heavy Ion Input Block
    if(openBlock) then
      continue
    elseif(closeBlock) then
      has_hion = .true.
    else
      call hions_parseInputLine(inLine,blockLine,inErr)
      if(inErr) goto 9999
    end if

  case("ZIPF") ! ZIPF Input Block
    if(openBlock) then
      continue
    elseif(closeBlock) then
      call zipf_parseInputDone
    else
      call zipf_parseInputline(inLine,inErr)
      if(inErr) goto 9999
    end if

  case("SCAT") ! SCATTER Input Block
    if(openBlock) then
      scatter_active = .true.
      call scatter_allocate
    elseif(closeBlock) then
      if(scatter_debug) call scatter_dumpData
    else
      call scatter_parseInputLine(string(inLine))
    end if

  case("HDF5") ! HDF5 Input Block
#ifndef HDF5
    write(lout,"(a)") "INPUT> ERROR SixTrack was not compiled with the HDF5 flag."
    goto 9999
#else
    if(openBlock) then
      h5_isActive = .true.
    elseif(closeBlock) then
      continue
    else
      call h5_parseInputLine(string(inLine))
    end if
#endif

  case("HASH") ! HASH Library
#ifndef HASHLIB
    write(lout,"(a)") "INPUT> ERROR SixTrack was not compiled with the HASHLIB flag."
    goto 9999
#else
    if(openBlock) then
      continue
    elseif(closeBlock) then
      continue
    else
      call hash_parseInputLine(inLine,inErr)
      if(inErr) goto 9999
    end if
#endif

  case("ROOT") ! ROOT Input Block
#ifndef ROOT
    write(lout,"(a)") "INPUT> ERROR SixTrack was not compiled with the ROOT flag."
    goto 9999
#else
  if(openBlock) then
    continue
  elseif(closeBlock) then
    call root_parseInputDone
  else
    call root_daten(inLine,inErr)
    if(inErr) goto 9999
  end if
#endif

  case default ! Unknown Block, Time to Panic
    write(lout,"(a)") "INPUT> ERROR Unknown block '"//currBlock//"' encountered. Check your input file."
    goto 9999

  end select

  ! Some cleanup and record keeping

  if(openBlock) then
    ! Nothing to do other than switch off the flag
    openBlock = .false.
  end if

  if(closeBlock) then
    ! Close the block and check that everything is normal.
    if(blockClosed) then
      write(lout,"(a)") "INPUT> Trying to close block '"//currBlock//"', which is already closed."
      goto 9999
    end if
    call sixin_closeBlock(currBlock)
    currBlock  = "NONE"  ! No active blocks
    closeBlock = .false. ! Switch off the flag
  end if

  goto 110 ! Get next line

! ================================================================================================ !
!  DONE PARSING FORT.2 AND FORT.3
! ================================================================================================ !

! ================================================================================================ !
!  ENDE WAS REACHED
! ================================================================================================ !
9000 continue

  if(napx >= 1) then
    if(e0 < pieni .or. e0 < pma) then
      write(lout,"(a)") "ENDE> ERROR Kinetic energy of the particle is less or equal to zero."
      call prror
    end if

    call hions_postInput
    gammar = nucm0/e0
    betrel = sqrt((one+gammar)*(one-gammar))

    if(nbeam >= 1) then
      parbe14 = (((((-one*crad)*partnum)/four)/pi)/sixin_emitNX)*c1e6
    end if
    crad   = (((two*crad)*partnum)*gammar)*c1e6
    emitx  = sixin_emitNX*gammar
    emity  = sixin_emitNY*gammar

    if (do_coll) then
      call collimate_postInput(gammar)
    endif

    !Check for incompatible flags
    if (ipos == 1) then
      if (do_coll) then
        write(lout,'(a)') "ENDE> ERROR COLLimation block and POSTprocessing block are not compatible."
        call prror(-1)
      endif

      if (scatter_active) then
        write(lout,'(a)') "ENDE> ERROR SCATTER block and POSTprocessing block are not compatible."
        call prror(-1)
      endif
#ifdef FLUKA
      if (fluka_enable) then
        write(lout,'(a)') "ENDE> ERROR FLUKA block and POSTprocessing block are not compatible."
        call prror(-1)
      endif
#endif
    endif

  end if

  call elens_postInput

  if(idp == 0 .or. ition == 0 .or. nbeam < 1) then
    do j=1,il
      ! Converting 6D lenses to 4D
      if(beam_expflag == 1) then
        if(parbe(j,2) > 0) then
          parbe(j,2) = zero
          parbe(j,1) = parbe(j,7)
          parbe(j,3) = parbe(j,10)
        end if
      else
        parbe(j,2) = zero
      end if
    end do
  else
    do j=1,il
      if(parbe(j,2) > real(mbea,fPrec)) then
        write(lout,"(3(a,i5))") "ENDE> ERROR Requested ",int(parbe(j,2))," slices for 6D beam-beam element"//&
          " #",j," named '"//trim(bez(j))//"', maximum is mbea = ",mbea
        parbe(j,2) = real(mbea,fPrec)
        call prror(-1) ! Treat this warning as an error
      end if
    end do
  end if

  ! Done with checks. Write the report
  call sixin_blockReport

! ================================================================================================ !

  ! This is where the PRINT spam happens
  if(.not.st_print) goto 9500 ! Skip it

  write(lout,"(a)") ""
  write(lout,"(a)") "  *** RING PARAMETERS ***"
  write(lout,"(a)") ""

  ! Print Single Elements
  write(lout,"(a)") "  SINGLE ELEMENTS:"
  write(lout,"(a)") ""
  write(lout,"(a)") "   NO NAME                TYP  1/RHO         STRENGTH      LENGTH        X-POS"//&
    "         X-RMS         Y-POS         Y-RMS"
  write(lout,"(a)") str_divLine
  il1=il
  if(sixin_ncy2 == 0) il1 = il-1
  do k=1,il1
    if(abs(kz(k)) == 12) then
      write(lout,"(i5,1x,a20,1x,i2,7(1x,e13.6))") k,bez(k)(1:20),kz(k),ed(k),ek(k),phasc(k),xpl(k),xrms(k),zpl(k),zrms(k)
      kz(k)=abs(kz(k))
      phasc(k)=phasc(k)*rad
    else
      write(lout,"(i5,1x,a20,1x,i2,7(1x,e13.6))") k,bez(k)(1:20),kz(k),ed(k),ek(k),el(k),xpl(k),xrms(k),zpl(k),zrms(k)
    end if
  end do
  write(lout,"(a)") str_divLine

  ! Print Ring Structure
  write(lout,"(a)")    ""
  write(lout,"(a)")    "  RING STRUCTURE:"
  write(lout,"(a)")    ""
  write(lout,"(a,i8)") "  Superperiods:     ",mper
  do k=1,mper
    write(lout,"(a,i8)") "  Symmetry:         ",msym(k)
  end do
  write(lout,"(a,i8)") "  Unique Blocks:    ",mblo
  write(lout,"(a,i8)") "  Block per Period: ",mbloz
  write(lout,"(a)")    ""
  write(lout,"(a)")    str_divLine

  ! Print Block Structure
  write(lout,"(a)") ""
  write(lout,"(a)") "  BLOCK STRUCTURE:"
  write(lout,"(a)") ""
  write(lout,"(a)") "   NO NAME                NUM  SINGLE ELEMENTS"
  write(lout,"(a)") str_divLine
  do l=1,mblo
    kk = mel(l)
    ll = kk/6
    if(ll /= 0) then
      do l1=1,ll
        l2 = (l1-1)*6+1
        l3 = l2+5
        if(l2 == 1) then
          write(lout,"(i5,1x,a20,1x,i3,6(1x,a20))") l,bezb(l),kk,(sixin_beze(l,k),k=1,6)
        else
          write(lout,"(30x,6(1x,a20))") (sixin_beze(l,k),k=l2,l3)
        end if
      end do
      if(mod(kk,6) /= 0) then
        l4 = ll*6+1
        write(lout,"(30x,6(1x,a20))") (sixin_beze(l,k),k=l4,kk)
      end if
    else
      write(lout,"(i5,1x,a20,1x,i3,6(1x,a20))") l,bezb(l),kk,(sixin_beze(l,k),k=1,kk)
    end if
  end do
  write(lout,"(a)") str_divLine

  ! Print Block Structure of Super Periods
  write(lout,"(a)") ""
  write(lout,"(a)") "  BLOCK STRUCTURE OF SUPER PERIODS:"
  write(lout,"(a)") ""
  write(lout,"(a)") "   NO NAME                NUM  SINGLE ELEMENTS"
  write(lout,"(a)") str_divLine
  mblozz=mbloz/5+1
  do k=1,mblozz
    k10 = (k-1)*5
    if((mbloz-k10) == 0) cycle
    do l=1,5
      if((k10+l) > mbloz) ic0(l) = " "
      if((k10+l) > mbloz) cycle
      icc = ic(k10+l)
      if(icc <= nblo) then
        ic0(l) = bezb(icc)
      else
        ic0(l) = sixin_bez0(icc-nblo)
      end if
    end do
    k11 = k10+1
    write(lout,"(i5,5(1x,a20))") k11,(ic0(l),l=1,5)
  end do
  write(lout,"(a)") str_divLine

  if(idp == 0) goto 8000
  ! Write out with BB parameters
  write(lout,"(a)")           ""
  write(lout,"(a)")           "  SYNCHROTRON OSCILLATIONS AND BEAM-BEAM:"
  write(lout,"(a)")           ""
  write(lout,"(a,i20)")       "  Number of cavities:                    ",ncy
  write(lout,"(a,f30.9)")     "  Momentum amplitude dP/P:               ",dp1
  write(lout,"(a,f30.9)")     "  Offset momentum amplitude dP/P:        ",dppoff
  write(lout,"(a,f30.9)")     "  Machine length in (m):                 ",tlen
  write(lout,"(a,f30.9)")     "  Particle mass (MeV):                   ",pma
  if(nbeam >= 1) then
    write(lout,"(a,f30.9)")   "  Particle number (1e9):                 ",abs(partnum/c1e9)
    if(partnum > zero) then
      write(lout,"(a,a20)")   "  Beams have same charge:                ","YES"
    else
      write(lout,"(a,a20)")   "  Beams have opposite charge:            ","YES"
    end if
    write(lout,"(a,f30.9)")   "  Beam-beam parameter:                   ",parbe14
    if(ibeco == 0) then
      write(lout,"(a,a20)")   "  Closed orbit due to beam-beam kick:    ","LEFT"
    else
      write(lout,"(a,a20)")   "  Closed orbit due to beam-beam kick:    ","SUBTRACTED"
    end if
    if(ibtyp == 0) then
      write(lout,"(a,a20)")   "  Fast beam-beam kick switch:            ","OFF"
    else
      write(lout,"(a,a20)")   "  Fast beam-beam kick switch:            ","ON"
    end if
    if(beam_expflag == 0) then
      if(ibb6d == 0) then
        write(lout,"(a,a20)") "  Hirata 6D:                             ","OFF"
      else
        write(lout,"(a,a20)") "  Hirata 6D:                             ","ON"
      end if
    end if
    if(ibbc == 0) then
      write(lout,"(a,a20)")   "  Consider linear coupling for BB:       ","OFF"
    else
      write(lout,"(a,a20)")   "  Consider linear coupling for BB:       ","ON"
    end if
    write(lout,"(a,f30.9)")   "  Bunch length:                          ",sigz
    write(lout,"(a,f30.9)")   "  Energy spread:                         ",sige
    write(lout,"(a,f30.9)")   "  Normalized horizontal emmittance (um): ",sixin_emitNX
    write(lout,"(a,f30.9)")   "  Normalized vertical emmittance (um):   ",sixin_emitNY
  end if
  write(lout,"(a,f30.9)")     "  Energy in (MeV):                       ",e0
  if(sixin_ncy2.eq.0) then
    write(lout,"(a,f30.9)")   "  Harmonic number:                       ",sixin_harm
    write(lout,"(a,f30.9)")   "  Circumf. voltage (MV):                 ",sixin_u0
    write(lout,"(a,f30.9)")   "  Equilibrium phase (deg):               ",sixin_phag
    write(lout,"(a,f30.9)")   "  Frequency (units of rev. freq.):       ",qs
    write(lout,"(a,f30.9)")   "  Momentum compaction:                   ",sixin_alc
  end if
  if(beam_expflag == 0) then
    if(ibb6d == 1) then
      write(lout,"(a)") ""
      write(lout,"(a)") "  HIRATA's 6D BEAM-BEAM ELEMENTS"
      write(lout,"(a)") ""
      write(lout,"(a)") "ELEMENT           #_OF_SLICES    CROSSING_ANGLE     CROSSING_PLANE     COUPLING_ANGLE"
      write(lout,"(a)") str_divLine
      do j=1,il
        if(parbe(j,2) > zero) then
          write(lout,"(t10,a16,5x,i4,7x,d17.10,2x,d17.10)") bez(j),int(parbe(j,2)),parbe(j,1),parbe(j,3)
        end if
      end do
    end if
    write(lout,"(a)") str_divLine
  elseif(beam_expflag == 1) then
    write(lout,"(a)") ""
    write(lout,"(a)") "  HIRATA's 6D BEAM-BEAM ELEMENTS"
    write(lout,"(a)") ""
    write(lout,"(a)") "ELEMENT           #_OF_SLICES    XING_ANGLE  XING_PLANE   HOR_SEP     VER_SEP"//&
      "        S11        S12        S22         S33         S34         S44         S13         S14         S23         S24"
    write(lout,"(a)") repeat("-",200)
    do j=1,il
      if(kz(j) == 20 .and. parbe(j,17) == 1) then
        write(lout,"(t10,a16,5x,i4,7x,1pe10.3,2x,1pe10.3,2x,1pe10.3,2x,1pe10.3,2x,1pe10.3,2x,1pe10.3&
        &,2x,1pe10.3,2x,1pe10.3,2x,1pe10.3,2x,1pe10.3,2x,1pe10.3,2x,1pe10.3,2x,1pe10.3,2x,1pe10.3)")&
        &bez(j),int(parbe(j,2)),parbe(j,1),parbe(j,3),parbe(j,5),parbe(j,6),parbe(j,7),parbe(j,8),  &
        &parbe(j,9),parbe(j,10),parbe(j,11),parbe(j,12),parbe(j,13),parbe(j,14),parbe(j,15),parbe(j,16)
      end if
    end do
    write(lout,"(a)") repeat("-",200)
    write(lout,"(a)") ""
    write(lout,"(a)") "  4D BEAM-BEAM ELEMENTS"
    write(lout,"(a)") ""
    write(lout,"(a)") "ELEMENT           #_OF_SLICES        S11        S22       HOR_SEP     VER_SEP"
    write(lout,"(a)") str_divLine
    do j=1,il
      if (kz(j) == 20 .and. parbe(j,17) == 0) then
        write(lout,"(t10,a16,5x,i4,7x,1pe10.3,2x,1pe10.3,2x,1pe10.3,2x,1pe10.3)") &
          bez(j),int(parbe(j,2)),parbe(j,1),parbe(j,3),parbe(j,5),parbe(j,6)
      end if
    end do
    write(lout,"(a)") str_divLine
  end if

8000 continue
  write(lout,"(a)") ""
  write(lout,"(a)") "  *** TRACKING PARAMETERS ***"
  write(lout,"(a)") ""

  nfb = nde(1)
  nac = nde(2)
  nft = numl-nde(2)
  if(numl <= nde(2)) nft = 0
  if(numl <= nde(2)) nac = numl
  if(numl <= nde(1)) nac = 0
  if(numl <= nde(1)) nfb = numl

  write(lout,"(a,i20)")       "  Number of revolutions:                 ",numl
  write(lout,"(a,i20)")       "  Number of reverse-revolutions:         ",numlr
  write(lout,"(a,i20)")       "  Turns per coor.-printout:              ",nwr(4)
  write(lout,"(a,i20)")       "  Flat bottom up to turn:                ",nfb
  write(lout,"(a,i20)")       "  Turns per print on dataset:            ",nwr(1)
  write(lout,"(a,i20)")       "  Acceleration up to turn:               ",nac
  write(lout,"(a,i20)")       "  Turns per print on dataset:            ",nwr(2)
  write(lout,"(a,i20)")       "  Flat top number of turns:              ",nft
  write(lout,"(a,i20)")       "  Turns per print on dataset:            ",nwr(3)
  write(lout,"(a,i20)")       "  Tracking start at element no.:         ",kanf
  write(lout,"(a,f30.9)")     "  Initial amplitude-h in (mm):           ",amp(1)
  write(lout,"(a,f30.9)")     "  Coupling  eps-y/eps-x:                 ",rat
  write(lout,"(a,i20)")       "  Number of C.O. iterations:             ",itco
  write(lout,"(a,e34.9)")     "  Precision of C.O. deviation:           ",dma
  write(lout,"(a,e34.9)")     "  Precision of C.O. slope:               ",dmap
  write(lout,"(a,i20)")       "  Number of q-adj. iterations:           ",itqv
  write(lout,"(a,e34.9)")     "  Change in k-strength by:               ",dkq
  write(lout,"(a,e34.9)")     "  Precision of q-adjustement:            ",dqq
  write(lout,"(a,i20)")       "  Number of chromat.-adj. iter.:         ",itcro
  write(lout,"(a,e34.9)")     "  Change in sex.-strength by:            ",dsm0
  write(lout,"(a,e34.9)")     "  Precision of chromat.-adj.:            ",dech
  write(lout,"(a,e34.9)")     "  DP-interval f. cromat.-adj.:           ",de0
  write(lout,"(a,e34.9)")     "  DP-interval for dispersion:            ",ded
  write(lout,"(a,e34.9)")     "  Precision for C.O. RMS:                ",dsi
  write(lout,"(a)") ""
  write(lout,"(a)") str_divLine
  write(lout,"(a)") ""
  write(lout,"(a)") "    OOOOOOOOOOOOOOOOOOOOO"
  write(lout,"(a)") "    OO                 OO"
  write(lout,"(a)") "    OO  PREPROCESSING  OO"
  write(lout,"(a)") "    OO                 OO"
  write(lout,"(a)") "    OOOOOOOOOOOOOOOOOOOOO"
  write(lout,"(a)") ""
  write(lout,"(a)") str_divLine

9500 continue

  call dealloc(sixin_bez0,mNameLen,"sixin_bez0")

  return

! ================================================================================================ !
!  END OF INPUT PARSING
! ================================================================================================ !

9997 continue
  write(lout,"(a)") "INPUT> ERROR Header could not be read from fort.3"
  call prror
  return

9998 continue
  write(lout,"(a,i0,a)") "INPUT> ERROR fort.",nUnit," is missing or empty, or end was reached without an ENDE flag."
  call prror
  return

9999 continue
  if(nUnit == 2) then
    write(lout,"(a)")      "INPUT> ERROR in fort.2"
    write(lout,"(a,i0,a)") "INPUT> Line ",lineNo2,": '"//trim(inLine)//"'"
  else
    write(lout,"(a)")      "INPUT> ERROR in fort.3"
    write(lout,"(a,i0,a)") "INPUT> Line ",lineNo3,": '"//trim(inLine)//"'"
    write(lout,"(a)")      "INPUT> Previous Lines:"
    write(lout,"(a)")      ""
    do i=1,5
      if(lineNo3-5+i <= 0) cycle
      write(lout,"(i5,a)") lineNo3-5+i," | "//trim(pLines(i))
    end do
  end if
  call prror
  return

end subroutine daten

! ================================================================================================ !
! purpose:                                                             *
!   modification of wwerf, real(kind=fPrec) complex error function,    *
!   written at cern by k. koelbig.                                     *
!   taken from mad8                                                    *
! input:                                                               *
!   xx, yy    (real)    argument to cerf.                              *
! output:                                                              *
!   wx, wy    (real)    function result.                               *
! ================================================================================================ !
subroutine errf(xx,yy,wx,wy)
  ! real(kind=fPrec) version.
  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  use mod_common, only : cc,xlim,ylim
  implicit none

  integer n,nc,nu
  real(kind=fPrec) h,q,rx,ry,saux,sx,sy,tn,tx,ty,wx,wy,x,xh,xl,xx,y,yh,yy
  dimension rx(33),ry(33)
  save
!-----------------------------------------------------------------------
  x=abs(xx)
  y=abs(yy)
  if(y.lt.ylim.and.x.lt.xlim) then
    q=(one-y/ylim)*sqrt(one-(x/xlim)**2)
    h=one/(3.2_fPrec*q)
    nc=7+int(23.0_fPrec*q)                                               !hr05
!       xl=h**(1-nc)
    xl=exp_mb((1-nc)*log_mb(h))                                      !yil11
#ifdef DEBUG
!       call wda('errfq',q,nc,0,0,0)
!       call wda('errfh',h,nc,0,0,0)
!       call wda('errfxl',xl,nc,0,0,0)
#endif
#ifdef DEBUG
!       call wda('errfxlrn',xl,nc,0,0,0)
#endif
    xh=y+half/h
    yh=x
    nu=10+int(21.0_fPrec*q)
#ifdef DEBUG
!       call wda('errfxh',xh,nu,0,0,0)
!       call wda('errfyh',yh,nu,0,0,0)
#endif
    rx(nu+1)=zero
    ry(nu+1)=zero
    do 10 n=nu,1,-1
      tx=xh+real(n,fPrec)*rx(n+1)                                          !hr05
      ty=yh-real(n,fPrec)*ry(n+1)                                          !hr05
      tn=tx**2+ty**2                                                 !hr05
      rx(n)=(half*tx)/tn                                            !hr05
      ry(n)=(half*ty)/tn                                            !hr05
10   continue
    sx=zero
    sy=zero
    do 20 n=nc,1,-1
      saux=sx+xl
      sx=rx(n)*saux-ry(n)*sy
      sy=rx(n)*sy+ry(n)*saux
      xl=h*xl
20   continue
    wx=cc*sx
    wy=cc*sy
  else
    xh=y
    yh=x
    rx(1)=zero
    ry(1)=zero
    do 30 n=9,1,-1
      tx=xh+real(n,fPrec)*rx(1)                                            !hr05
      ty=yh-real(n,fPrec)*ry(1)                                            !hr05
      tn=tx**2+ty**2                                                 !hr05
      rx(1)=(half*tx)/tn                                            !hr05
      ry(1)=(half*ty)/tn                                            !hr05
30   continue
    wx=cc*rx(1)
    wy=cc*ry(1)
  endif
!      if(y.eq.0.) wx=exp(-x**2)
  if(yy.lt.zero) then
    wx=(two*exp_mb(y**2-x**2))*cos_mb((two*x)*y)-wx                  !hr05
    wy=((-one*two)*exp_mb(y**2-x**2))*sin_mb((two*x)*y)-wy           !hr05
    if(xx.gt.zero) wy=-one*wy                                        !hr05
  else
    if(xx.lt.zero) wy=-one*wy
  endif
end subroutine errf

! ================================================================================================ !
!  subroutine wzsubv
!
!  This subroutine sets u=real(w(z)) and v=imag(w(z)), where z=x+i*y and
!  where w(z) is the complex error function defined by formula 7.1.3 in
!  "Handbook of Mathematical functions [eds. M.Abramowitz & I.A.Stegun,
!  Washington, 1966].  The absolute error of the computed value is less
!  than 1E-8.
!
!  *** Note.  Subroutine WZSET must have been called before this sub-
!  routine can be used.
!
!  For (x,y) inside the rectangle with opposite corners (xcut,0) and
!  (0,ycut), where xcut and ycut have been set by WZSET, an interpo-
!  lation formula is used.  For (x,y) outside this rectangle, a two-
!  term rational approximation is used.
!
!  (G.A.Erskine, 29.09.1997)
!
!  Vectorised for up to 64 argument values by E.McIntosh, 30.10.1997.
!  Much impoved using short vector buffers Eric 1st May, 2014.
!
!  Third-order divided-difference interpolation over the corners of a
!  square [e.g. formula (2.5.1) in "Introduction to Numerical Analysis"
!  (F.B.Hildebrand New York, 1957), but with complex nodes and
!  function values].
!
!  In the interpolation formula the corners of the grid square contain-
!  ing (x,y) are numbered (0,0)=3, (h,0)=4, (h,h)=1, (0,h)=2.
!  Identifiers d, dd and ddd denote divided-differences of orders 1, 2
!  and 3 respectively, and a preceding 't' indicates twice the value.
!
!
!  Two-term rational approximation to w(z) [Footnote to Table 7.9
!  in "Handbook of Mathematical Functions (eds. M.Abramowitz &
!  I.A.Stegun, Washington, 1966), but with additional digits in
!  the constants]:
!              u+i*v = i*z*( a1/(z**2-b1) + a2/(z**2-b2) ).
!  Maximum absolute error:
!        <1.E-6  for  x>=4.9  or  y>=4.4
!        <1.E-7  for  x>=6.1  or  y>=5.7
!        <1.E-8  for  x>=7.8  or  y>=7.5
!
! ================================================================================================ !
subroutine wzsubv(n,vx,vy,vu,vv)

  use parpro, only : npart
  use floatPrecision
  use numerical_constants
  implicit none

  dimension vx(*),vy(*),vu(*),vv(*)
  integer i,j,k,n,vmu,vnu
  real(kind=fPrec) a1,a2,b1,b2,vd12i,vd12r,vd23i,vd23r,vd34i,vd34r,vp,vq,vqsq,vr,vsimag,vsreal,vt,  &
    vtdd13i,vtdd13r,vtdd24i,vtdd24r,vtdddi,vtdddr,vti,vtr,vu,vusum,vusum3,vv,vvsum,vvsum3,vw1i,vw1r,&
    vw2i,vw2r,vw3i,vw3r,vw4i,vw4r,vx,vxh,vxhrel,vy,vyh,vyhrel
  integer idim,kstep,nx,ny
  real(kind=fPrec) h,hrecip,wtimag,wtreal,xcut,ycut
  parameter ( xcut = 7.77_fPrec, ycut = 7.46_fPrec )
  parameter ( h = one/63.0_fPrec )
  parameter ( nx = 490, ny = 470 )
  parameter ( idim = (nx+2)*(ny+2) )
  common /wzcom1/ hrecip, kstep
  common /wzcom2/ wtreal(idim), wtimag(idim)
  parameter ( a1 = 0.5124242248_fPrec, a2 = 0.0517653588_fPrec )
  parameter ( b1 = 0.2752551286_fPrec, b2 = 2.7247448714_fPrec )
  real(kind=fPrec) xm,xx,yy
  parameter (xm=1e16_fPrec)
!     temporary arrays to facilitate vectorisation
  integer in,out,ins,outs
  dimension ins(npart),outs(npart)
!-----------------------------------------------------------------------
  save
  in=0
  out=0
  do i=1,n
    if (vx(i).ge.xcut.or.vy(i).ge.ycut) then
      out=out+1
      outs(out)=i
      if (out.eq.npart) then
!     everything outside the rectangle so approximate
!     write (*,*) 'ALL outside'
!     write (*,*) 'i=',i
        do j=1,out
          xx=vx(outs(j))
          yy=vy(outs(j))
          if (xx.ge.xm) xx=xm
          if (yy.ge.xm) yy=xm
          vp=xx**2-yy**2
          vq=(two*xx)*yy
          vqsq=vq**2
          !  First term.
          vt=vp-b1
          vr=a1/(vt**2+vqsq)
          vsreal=vr*vt
          vsimag=-vr*vq
          !  Second term
          vt=vp-b2
          vr=a2/(vt**2+vqsq)
          vsreal=vsreal+vr*vt
          vsimag=vsimag-vr*vq
          !  Multiply by i*z.
          vu(outs(j))=-(yy*vsreal+xx*vsimag)
          vv(outs(j))=xx*vsreal-yy*vsimag
        enddo
        out=0
      endif
    else
      in=in+1
      ins(in)=i
      if (in.eq.npart) then
!     everything inside the square, so interpolate
!     write (*,*) 'ALL inside'
        do j=1,in
          vxh = hrecip*vx(ins(j))
          vyh = hrecip*vy(ins(j))
          vmu = int(vxh)
          vnu = int(vyh)
!  Compute divided differences.
          k = 2 + vmu + vnu*kstep
          vw4r = wtreal(k)
          vw4i = wtimag(k)
          k = k - 1
          vw3r = wtreal(k)
          vw3i = wtimag(k)
          vd34r = vw4r - vw3r
          vd34i = vw4i - vw3i
          k = k + kstep
          vw2r = wtreal(k)
          vw2i = wtimag(k)
          vd23r = vw2i - vw3i
          vd23i = vw3r - vw2r
          vtr = vd23r - vd34r
          vti = vd23i - vd34i
          vtdd24r = vti - vtr
          vtdd24i = -one* ( vtr + vti )                             !hr05
          k = k + 1
          vw1r = wtreal(k)
          vw1i = wtimag(k)
          vd12r = vw1r - vw2r
          vd12i = vw1i - vw2i
          vtr = vd12r - vd23r
          vti = vd12i - vd23i
          vtdd13r = vtr + vti
          vtdd13i = vti - vtr
          vtdddr = vtdd13i - vtdd24i
          vtdddi = vtdd24r - vtdd13r
!  Evaluate polynomial.
          vxhrel = vxh - real(vmu,fPrec)
          vyhrel = vyh - real(vnu,fPrec)
          vusum3=half*(vtdd13r+(vxhrel*vtdddr-vyhrel*vtdddi))
          vvsum3=half*(vtdd13i+(vxhrel*vtdddi+vyhrel*vtdddr))
          vyhrel = vyhrel - one
          vusum=vd12r+(vxhrel*vusum3-vyhrel*vvsum3)
          vvsum=vd12i+(vxhrel*vvsum3+vyhrel*vusum3)
          vxhrel = vxhrel - one
          vu(ins(j))=vw1r+(vxhrel*vusum-vyhrel*vvsum)
          vv(ins(j))=vw1i+(vxhrel*vvsum+vyhrel*vusum)
        enddo
        in=0
      endif
    endif
  enddo
!     everything outside the rectangle so approximate
!     write (*,*) 'ALL outside'
!     write (*,*) 'i=',i
  do j=1,out
    xx=vx(outs(j))
    yy=vy(outs(j))
    if (xx.ge.xm) xx=xm
    if (yy.ge.xm) yy=xm
    vp=xx**2-yy**2
    vq=(two*xx)*yy
    vqsq=vq**2
!  First term.
    vt=vp-b1
    vr=a1/(vt**2+vqsq)
    vsreal=vr*vt
    vsimag=-vr*vq
!  Second term
    vt=vp-b2
    vr=a2/(vt**2+vqsq)
    vsreal=vsreal+vr*vt
    vsimag=vsimag-vr*vq
!  Multiply by i*z.
    vu(outs(j))=-(yy*vsreal+xx*vsimag)
    vv(outs(j))=xx*vsreal-yy*vsimag
  enddo
!     everything inside the square, so interpolate
!     write (*,*) 'ALL inside'
  do j=1,in
    vxh = hrecip*vx(ins(j))
    vyh = hrecip*vy(ins(j))
    vmu = int(vxh)
    vnu = int(vyh)
!  Compute divided differences.
    k = 2 + vmu + vnu*kstep
    vw4r = wtreal(k)
    vw4i = wtimag(k)
    k = k - 1
    vw3r = wtreal(k)
    vw3i = wtimag(k)
    vd34r = vw4r - vw3r
    vd34i = vw4i - vw3i
    k = k + kstep
    vw2r = wtreal(k)
    vw2i = wtimag(k)
    vd23r = vw2i - vw3i
    vd23i = vw3r - vw2r
    vtr = vd23r - vd34r
    vti = vd23i - vd34i
    vtdd24r = vti - vtr
    vtdd24i = -one* ( vtr + vti )                             !hr05
    k = k + 1
    vw1r = wtreal(k)
    vw1i = wtimag(k)
    vd12r = vw1r - vw2r
    vd12i = vw1i - vw2i
    vtr = vd12r - vd23r
    vti = vd12i - vd23i
    vtdd13r = vtr + vti
    vtdd13i = vti - vtr
    vtdddr = vtdd13i - vtdd24i
    vtdddi = vtdd24r - vtdd13r
!  Evaluate polynomial.
    vxhrel = vxh - real(vmu,fPrec)
    vyhrel = vyh - real(vnu,fPrec)
    vusum3=half*(vtdd13r+(vxhrel*vtdddr-vyhrel*vtdddi))
    vvsum3=half*(vtdd13i+(vxhrel*vtdddi+vyhrel*vtdddr))
    vyhrel = vyhrel - one
    vusum=vd12r+(vxhrel*vusum3-vyhrel*vvsum3)
    vvsum=vd12i+(vxhrel*vvsum3+vyhrel*vusum3)
    vxhrel = vxhrel - one
    vu(ins(j))=vw1r+(vxhrel*vusum-vyhrel*vvsum)
    vv(ins(j))=vw1i+(vxhrel*vvsum+vyhrel*vusum)
  enddo
  return
end subroutine wzsubv

! ================================================================================================ !
!  subroutine wzsub
!
!  This subroutine sets u=real(w(z)) and v=imag(w(z)), where z=x+i*y and
!  where w(z) is the complex error function defined by formula 7.1.3 in
!  "Handbook of Mathematical functions [eds. M.Abramowitz & I.A.Stegun,
!  Washington, 1966].  The absolute error of the computed value is less
!  than 1E-8.
!
!  *** Note.  Subroutine WZSET must have been called before this sub-
!  routine can be used.
!
!  For (x,y) inside the rectangle with opposite corners (xcut,0) and
!  (0,ycut), where xcut and ycut have been set by WZSET, an interpo-
!  lation formula is used.  For (x,y) outside this rectangle, a two-
!  term rational approximation is used.
!
!  (G.A.Erskine, 29.09.1997)
!
!
!  Third-order divided-difference interpolation over the corners of a
!  square [e.g. formula (2.5.1) in "Introduction to Numerical Analysis"
!  (F.B.Hildebrand New York, 1957), but with complex nodes and
!  function values].
!
!  In the interpolation formula the corners of the grid square contain-
!  ing (x,y) are numbered (0,0)=3, (h,0)=4, (h,h)=1, (0,h)=2.
!  Identifiers d, dd and ddd denote divided-differences of orders 1, 2
!  and 3 respectively, and a preceding 't' indicates twice the value.
!
! ================================================================================================ !
subroutine wzsub(x,y,u,v)

  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  use parpro
  use parbeam
  implicit none
  integer k,mu,nu
  real(kind=fPrec) a1,a2,b1,b2,d12i,d12r,d23i,d23r,d34i,d34r,p,q,qsq,r,simag,sreal,t,tdd13i,tdd13r, &
    tdd24i,tdd24r,tdddi,tdddr,ti,tr,u,usum,usum3,v,vsum,vsum3,w1i,w1r,w2i,w2r,w3i,w3r,w4i,w4r,x,xh, &
    xhrel,y,yh,yhrel
  parameter ( a1 = 0.5124242248_fPrec, a2 = 0.0517653588_fPrec )
  parameter ( b1 = 0.2752551286_fPrec, b2 = 2.7247448714_fPrec )
  save
!-----------------------------------------------------------------------
  if ( x.ge.xcut .or. y.ge.ycut ) goto 1000
  xh = hrecip*x
  yh = hrecip*y
  mu = int(xh)
  nu = int(yh)
!  Compute divided differences.
  k = 2 + mu + nu*kstep
  w4r = wtreal(k)
  w4i = wtimag(k)
  k = k - 1
  w3r = wtreal(k)
  w3i = wtimag(k)
  d34r = w4r - w3r
  d34i = w4i - w3i
  k = k + kstep
  w2r = wtreal(k)
  w2i = wtimag(k)
  d23r = w2i - w3i
  d23i = w3r - w2r
  tr = d23r - d34r
  ti = d23i - d34i
  tdd24r = ti - tr
  tdd24i = -one* ( tr + ti )                                         !hr05
  k = k + 1
  w1r = wtreal(k)
  w1i = wtimag(k)
  d12r = w1r - w2r
  d12i = w1i - w2i
  tr = d12r - d23r
  ti = d12i - d23i
  tdd13r = tr + ti
  tdd13i = ti - tr
  tdddr = tdd13i - tdd24i
  tdddi = tdd24r - tdd13r
!  Evaluate polynomial.
  xhrel = xh - real(mu,fPrec)
  yhrel = yh - real(nu,fPrec)
  usum3 = half*( tdd13r + ( xhrel*tdddr - yhrel*tdddi ) )
  vsum3 = half*( tdd13i + ( xhrel*tdddi + yhrel*tdddr ) )
  yhrel = yhrel - one
  usum = d12r + ( xhrel*usum3 - yhrel*vsum3 )
  vsum = d12i + ( xhrel*vsum3 + yhrel*usum3 )
  xhrel = xhrel - one
  u = w1r + ( xhrel*usum - yhrel*vsum )
  v = w1i + ( xhrel*vsum + yhrel*usum )
  return
!
!  Two-term rational approximation to w(z) [Footnote to Table 7.9
!  in "Handbook of Mathematical Functions (eds. M.Abramowitz &
!  I.A.Stegun, Washington, 1966), but with additional digits in
!  the constants]:
!              u+i*v = i*z*( a1/(z**2-b1) + a2/(z**2-b2) ).
!  Maximum absolute error:
!        <1.E-6  for  x>=4.9  or  y>=4.4
!        <1.E-7  for  x>=6.1  or  y>=5.7
!        <1.E-8  for  x>=7.8  or  y>=7.5
!
1000 p=x**2-y**2
  q=(2.d0*x)*y                                                       !hr05
  qsq=q**2
!  First term.
  t=p-b1
  r=a1/(t**2+qsq)
  sreal=r*t
  simag=(-one*r)*q                                                   !hr05
!  Second term
  t=p-b2
  r=a2/(t**2+qsq)
  sreal=sreal+r*t
  simag=simag-r*q
!  Multiply by i*z.
  u=-one*(y*sreal+x*simag)                                           !hr05
  v=x*sreal-y*simag
  return
!
end subroutine wzsub

subroutine initialize_element(ix,lfirst)
!-----------------------------------------------------------------------
!     K.Sjobak & A.Santamaria, BE-ABP/HSS
!     last modified: 23-12-2016
!     Initialize a lattice element with index elIdx,
!     such as done when reading fort.2 (GEOM) and in DYNK.
!
!     Never delete an element from the lattice, even if it is not making a kick.
!     If the element is not recognized, do nothing (for now).
!     If trying to initialize an element (not lfirst) which is disabled,
!     print an error and exit.
!-----------------------------------------------------------------------

      use floatPrecision
      use dynk, only : dynk_elemData, dynk_izuIndex
      use numerical_constants
      use crcoall
      use string_tools
      use parpro
      use parbeam, only : beam_expflag,beam_expfile_open
      use mod_common
      use mod_commont
      use mod_commonmn
      use mod_hions
      use elens
      use wire
      use mathlib_bouncer
      implicit none

      integer, intent(in) :: ix
      logical, intent(in) :: lfirst

      !Temp variables
      integer i, m, k, im, nmz, izu, ibb, ii,j
      real(kind=fPrec) r0, r0a, bkitemp,sfac1,sfac2,sfac2s,sfac3,sfac4,sfac5,  crkveb_d, cikveb_d, &
      rho2b_d,tkb_d,r2b_d,rb_d,rkb_d,xrb_d,zrb_d,cbxb_d,cbzb_d,crxb_d,crzb_d,xbb_d,zbb_d, napx0
      real(kind=fPrec) crkveb(npart),cikveb(npart),rho2b(npart),tkb(npart),r2b(npart),rb(npart),        &
      rkb(npart),xrb(npart),zrb(npart),xbb(npart),zbb(npart),crxb(npart),crzb(npart),cbxb(npart),     &
      cbzb(npart)

      integer :: nbeaux(nbb)

!--Nonlinear Elements
! TODO: Merge these cases into 1 + subcases?
      if(abs(kz(ix)).eq.1) then
         if(.not.lfirst) then
            do i=1,iu
               if ( ic(i)-nblo.eq.ix ) then
                 if(ktrack(i).eq.31) goto 100 !ERROR
                 sm(ix)=ed(ix)             ! Also done in envar() which is called from clorb()
                 smiv(i)=sm(ix)+smizf(i) ! Also done in program maincr
                 smi(i)=smiv(i)          ! Also done in program maincr
#include "include/stra01.f90"
               endif
            enddo
         endif

      elseif(abs(kz(ix)).eq.2) then
         if(.not.lfirst) then
            do i=1,iu
               if ( ic(i)-nblo.eq.ix ) then
                 if(ktrack(i).eq.31) goto 100 !ERROR
                 sm(ix)=ed(ix)
                 smiv(i)=sm(ix)+smizf(i)
                 smi(i)=smiv(i)
#include "include/stra02.f90"
               endif
            enddo
         endif
      elseif(abs(kz(ix)).eq.3) then
         if(.not.lfirst) then
            do i=1,iu
               if ( ic(i)-nblo.eq.ix ) then
                 if(ktrack(i).eq.31) goto 100 !ERROR
                 sm(ix)=ed(ix)
                 smiv(i)=sm(ix)+smizf(i)
                 smi(i)=smiv(i)
#include "include/stra03.f90"
               endif
            enddo
         endif

      elseif(abs(kz(ix)).eq.4) then
         if(.not.lfirst) then
            do i=1,iu
               if ( ic(i)-nblo.eq.ix ) then
                 if(ktrack(i).eq.31) goto 100 !ERROR
                 sm(ix)=ed(ix)
                 smiv(i)=sm(ix)+smizf(i)
                 smi(i)=smiv(i)
#include "include/stra04.f90"
               endif
            enddo
         endif

      elseif(abs(kz(ix)).eq.5) then
         if(.not.lfirst) then
            do i=1,iu
               if ( ic(i)-nblo.eq.ix ) then
                 if(ktrack(i).eq.31) goto 100 !ERROR
                 sm(ix)=ed(ix)
                 smiv(i)=sm(ix)+smizf(i)
                 smi(i)=smiv(i)
#include "include/stra05.f90"
               endif
            enddo
         endif

      elseif(abs(kz(ix)).eq.6) then
         if(.not.lfirst) then
            do i=1,iu
               if ( ic(i)-nblo.eq.ix ) then
                 if(ktrack(i).eq.31) goto 100 !ERROR
                 sm(ix)=ed(ix)
                 smiv(i)=sm(ix)+smizf(i)
                 smi(i)=smiv(i)
#include "include/stra06.f90"
               endif
            enddo
         endif

      elseif(abs(kz(ix)).eq.7) then
         if(.not.lfirst) then
            do i=1,iu
               if ( ic(i)-nblo.eq.ix ) then
                 if(ktrack(i).eq.31) goto 100 !ERROR
                 sm(ix)=ed(ix)
                 smiv(i)=sm(ix)+smizf(i)
                 smi(i)=smiv(i)
#include "include/stra07.f90"
               endif
            enddo
         endif

      elseif(abs(kz(ix)).eq.8) then
         if(.not.lfirst) then
            do i=1,iu
               if ( ic(i)-nblo.eq.ix ) then
                 if(ktrack(i).eq.31) goto 100 !ERROR
                 sm(ix)=ed(ix)
                 smiv(i)=sm(ix)+smizf(i)
                 smi(i)=smiv(i)
#include "include/stra08.f90"
               endif
            enddo
         endif

      elseif(abs(kz(ix)).eq.9) then
         if(.not.lfirst) then
            do i=1,iu
               if ( ic(i)-nblo.eq.ix ) then
                 if(ktrack(i).eq.31) goto 100 !ERROR
                 sm(ix)=ed(ix)
                 smiv(i)=sm(ix)+smizf(i)
                 smi(i)=smiv(i)
#include "include/stra09.f90"
               endif
            enddo
         endif

      elseif(abs(kz(ix)).eq.10) then
         if(.not.lfirst) then
            do i=1,iu
               if ( ic(i)-nblo.eq.ix ) then
                 if(ktrack(i).eq.31) goto 100 !ERROR
                 sm(ix)=ed(ix)
                 smiv(i)=sm(ix)+smizf(i)
                 smi(i)=smiv(i)
#include "include/stra10.f90"
               endif
            enddo
         endif

!--Multipoles
      elseif(kz(ix).eq.11) then
        if (lfirst) then
           if (abs(el(ix)+one).le.pieni) then
              dki(ix,1) = ed(ix)
              dki(ix,3) = ek(ix)
              ed(ix) = one
              ek(ix) = one
              el(ix) = zero
           else if(abs(el(ix)+two).le.pieni) then
              dki(ix,2) = ed(ix)
              dki(ix,3) = ek(ix)
              ed(ix) = one
              ek(ix) = one
              el(ix) = zero
           endif
        else
          do i=1,iu
            if ( ic(i)-nblo.eq.ix ) then
              nmz=nmu(ix)
              im=irm(ix)
              do k=1,nmz
                aaiv(k,i)=scalemu(im)*(ak0(im,k)+amultip(k,i)*aka(im,k))
                bbiv(k,i)=scalemu(im)*(bk0(im,k)+bmultip(k,i)*bka(im,k))
              end do
            endif
          enddo
        endif


!--Cavities (ktrack = 2 for thin)
      elseif(abs(kz(ix)).eq.12) then
         !Moved from daten
         phasc(ix) = el(ix)
         el(ix) = zero
         dynk_elemData(ix,3) = phasc(ix)
         if (.not.lfirst) then

            ! Doesn't work, as i is not initialized here.
            !if (.not.ktrack(i).eq.2) goto 100 !ERROR

            phasc(ix) = phasc(ix)*rad

            hsyc(ix) = ((two*pi)*ek(ix))/tlen         ! daten SYNC block
            hsyc(ix)=(c1m3*hsyc(ix))*real(itionc(ix),fPrec) ! trauthin/trauthck
         endif
!--BEAM-BEAM
      elseif(kz(ix).eq.20) then

        if (lfirst) then
          ptnfac(ix)=el(ix)
          el(ix)=zero
          parbe(ix,5) = ed(ix)
          ed(ix)=zero
          parbe(ix,6) = ek(ix)
          ek(ix)=zero
          endif
! This is to inialize all the beam-beam element before the tracking (or to update it for DYNK).
        if (.not.lfirst) then
          do i=1,iu
            if ( ic(i)-nblo.eq.ix ) then
              ibb=imbb(i)
              if(parbe(ix,2).gt.zero) then
                if(beam_expflag.eq.1) then
                 bbcu(ibb,1)=parbe(ix,7)
                 bbcu(ibb,4)=parbe(ix,8)
                 bbcu(ibb,6)=parbe(ix,9)
                 bbcu(ibb,2)=parbe(ix,10)
                 bbcu(ibb,9)=parbe(ix,11)
                 bbcu(ibb,10)=parbe(ix,12)
                 bbcu(ibb,3)=parbe(ix,13)
                 bbcu(ibb,5)=parbe(ix,14)
                 bbcu(ibb,7)=parbe(ix,15)
                 bbcu(ibb,8)=parbe(ix,16)
                  do ii=1,10
                    bbcu(ibb,ii)=bbcu(ibb,ii)*c1m6
                  enddo
                endif
                ktrack(i)=44
                parbe(ix,4)=(((-one*crad)*ptnfac(ix))*half)*c1m6
                if(ibeco.eq.1) then
                  track6d(1,1)=parbe(ix,5)*c1m3
                  track6d(2,1)=zero
                  track6d(3,1)=parbe(ix,6)*c1m3
                  track6d(4,1)=zero
                  track6d(5,1)=zero
                  track6d(6,1)=zero
                  napx0=napx
                  napx=1
                  call beamint(napx,track6d,parbe,sigz,bbcu,imbb(i),ix,ibtyp,ibbc, mtc)
                  beamoff(1,imbb(i))=track6d(1,1)*c1e3
                  beamoff(2,imbb(i))=track6d(3,1)*c1e3
                  beamoff(3,imbb(i))=track6d(5,1)*c1e3
                  beamoff(4,imbb(i))=track6d(2,1)*c1e3
                  beamoff(5,imbb(i))=track6d(4,1)*c1e3
                  beamoff(6,imbb(i))=track6d(6,1)
                  napx=napx0

                endif

              else if(parbe(ix,2).eq.zero) then
                if(beam_expflag.eq.1) then
                   bbcu(ibb,1)=parbe(ix,1)
                   bbcu(ibb,2)=parbe(ix,3)
                   bbcu(ibb,3)=parbe(ix,13)
                endif
                if(ibbc.eq.1) then
                  sfac1=bbcu(ibb,1)+bbcu(ibb,2)
                  sfac2=bbcu(ibb,1)-bbcu(ibb,2)
                  sfac2s=one
                  if(sfac2.lt.zero) sfac2s=-one
                  sfac3=sqrt(sfac2**2+(four*bbcu(ibb,3))*bbcu(ibb,3))
                  if(sfac3 > sfac1) then
                    write(lout,"(a)") "BEAMBEAM> ERROR 6D beam-beam with tilt not possible."
                    call prror(-1)
                  end if
                  sfac4=(sfac2s*sfac2)/sfac3
                  sfac5=(((-one*sfac2s)*two)*bbcu(ibb,3))/sfac3
                  sigman(1,ibb)=sqrt(((sfac1+sfac2*sfac4)+(two*bbcu(ibb,3))*sfac5)*half)
                  sigman(2,ibb)=sqrt(((sfac1-sfac2*sfac4)-(two*bbcu(ibb,3))*sfac5)*half)
                  bbcu(ibb,11)=sqrt(half*(one+sfac4))
                  bbcu(ibb,12)=(-one*sfac2s)*sqrt(half*(one-sfac4))
                  if(bbcu(ibb,3).lt.zero) bbcu(ibb,12)=-one*bbcu(ibb,12)
                else
                  bbcu(ibb,11)=one
                  sigman(1,ibb)=sqrt(bbcu(ibb,1))
                  sigman(2,ibb)=sqrt(bbcu(ibb,2))
                endif

!--round beam
                nbeaux(imbb(i))=0
                if(sigman(1,imbb(i)).eq.sigman(2,imbb(i))) then
                  if(nbeaux(imbb(i)).eq.2.or.nbeaux(imbb(i)).eq.3) then
                    write(lout,"(a)") "BEAMBEAM> ERROR At each interaction point the beam must be either "//&
                    "round or elliptical for all particles"
                    call prror(-1)
                  else
                    nbeaux(imbb(i))=1
                    sigman2(1,imbb(i))=sigman(1,imbb(i))**2
                  endif
                endif
  !--elliptic beam x>z
                if(sigman(1,imbb(i)).gt.sigman(2,imbb(i))) then
                  if(nbeaux(imbb(i)).eq.1.or.nbeaux(imbb(i)).eq.3) then
                    write(lout,"(a)") "BEAMBEAM> ERROR At each interaction point the beam must be either "//&
                    "round or elliptical for all particles"
                    call prror(-1)
                  else
                    nbeaux(imbb(i))=2
                    ktrack(i)=42
                    sigman2(1,imbb(i))=sigman(1,imbb(i))**2
                    sigman2(2,imbb(i))=sigman(2,imbb(i))**2
                    sigmanq(1,imbb(i))=sigman(1,imbb(i))/sigman(2,imbb(i))
                    sigmanq(2,imbb(i))=sigman(2,imbb(i))/sigman(1,imbb(i))
                  endif
                endif
  !--elliptic beam z>x
                if(sigman(1,imbb(i)).lt.sigman(2,imbb(i))) then
                  if(nbeaux(imbb(i)).eq.1.or.nbeaux(imbb(i)).eq.2) then
                    write(lout,"(a)") "BEAMBEAM> ERROR At each interaction point the beam must be either "//&
                    "round or elliptical for all particles"
                    call prror(-1)
                  else
                    nbeaux(imbb(i))=3
                    ktrack(i)=43
                    sigman2(1,imbb(i))=sigman(1,imbb(i))**2
                    sigman2(2,imbb(i))=sigman(2,imbb(i))**2
                    sigmanq(1,imbb(i))=sigman(1,imbb(i))/sigman(2,imbb(i))
                    sigmanq(2,imbb(i))=sigman(2,imbb(i))/sigman(1,imbb(i))
                  endif
                endif


                strack(i)=crad*ptnfac(ix)
                if(ibbc.eq.0) then
                  crkveb_d=parbe(ix,5)
                  cikveb_d=parbe(ix,6)
                else
                  crkveb_d=parbe(ix,5)*bbcu(imbb(i),11)+parbe(ix,6)*bbcu(imbb(i),12)
                  cikveb_d=parbe(ix,6)*bbcu(imbb(i),11)-parbe(ix,5)*bbcu(imbb(i),12)
                endif

                if(nbeaux(imbb(i)).eq.1) then
                  ktrack(i)=41
                  if(ibeco.eq.1) then
                  rho2b_d=crkveb_d**2+cikveb_d**2
                  tkb_d=rho2b_d/(two*sigman2(1,imbb(i)))
                  beamoff(4,imbb(i))=((strack(i)*crkveb_d)/rho2b_d)*(one-exp_mb(-one*tkb_d))
                  beamoff(5,imbb(i))=((strack(i)*cikveb_d)/rho2b_d)*(one-exp_mb(-one*tkb_d))
                  endif
                endif

                if(ktrack(i) == 42) then
                  if(ibeco == 1) then
                    r2b_d=two*(sigman2(1,imbb(i))-sigman2(2,imbb(i)))
                    rb_d=sqrt(r2b_d)
                    rkb_d=(strack(i)*pisqrt)/rb_d
                    xrb_d=abs(crkveb_d)/rb_d
                    zrb_d=abs(cikveb_d)/rb_d
                    if(ibtyp == 0) then
                      call errf(xrb_d,zrb_d,crxb_d,crzb_d)
                      tkb_d=(crkveb_d**2/sigman2(1,imbb(i))+cikveb_d**2/sigman2(2,imbb(i)))*half
                      xbb_d=sigmanq(2,imbb(i))*xrb_d
                      zbb_d=sigmanq(1,imbb(i))*zrb_d
                      call errf(xbb_d,zbb_d,cbxb_d,cbzb_d)
                    else if(ibtyp == 1) then
                      tkb_d=(crkveb_d**2/sigman2(1,imbb(i))+cikveb_d**2/sigman2(2,imbb(i)))*half
                      xbb_d=sigmanq(2,imbb(i))*xrb_d
                      zbb_d=sigmanq(1,imbb(i))*zrb_d
                    else
                      tkb_d = zero ! -Wmaybe-uninitialized
                    endif
                  else
                    rkb_d = zero ! -Wmaybe-uninitialized
                    tkb_d = zero ! -Wmaybe-uninitialized
                  end if
                  beamoff(4,imbb(i))=(rkb_d*(crzb_d-exp_mb(-one*tkb_d)*cbzb_d))*sign(one,crkveb_d)
                  beamoff(5,imbb(i))=(rkb_d*(crxb_d-exp_mb(-one*tkb_d)*cbxb_d))*sign(one,cikveb_d)
                endif

                if(ktrack(i) == 43) then
                  if(ibeco == 1) then
                    r2b_d=two*(sigman2(2,imbb(i))-sigman2(1,imbb(i)))
                    rb_d=sqrt(r2b_d)
                    rkb_d=(strack(i)*pisqrt)/rb_d
                    xrb_d=abs(crkveb_d)/rb_d
                    zrb_d=abs(cikveb_d)/rb_d
                    if(ibtyp == 0) then
                      call errf(zrb_d,xrb_d,crzb_d,crxb_d)
                      tkb_d=(crkveb_d**2/sigman2(1,imbb(i))+cikveb_d**2/sigman2(2,imbb(i)))*half
                      xbb_d=sigmanq(2,imbb(i))*xrb_d
                      zbb_d=sigmanq(1,imbb(i))*zrb_d
                      call errf(zbb_d,xbb_d,cbzb_d,cbxb_d)
                    else if(ibtyp == 1) then
                      tkb_d=(crkveb_d**2/sigman2(1,imbb(i))+cikveb_d**2/sigman2(2,imbb(i)))*half
                      xbb_d=sigmanq(2,imbb(i))*xrb_d
                      zbb_d=sigmanq(1,imbb(i))*zrb_d
                    else
                      tkb_d = zero ! -Wmaybe-uninitialized
                    endif
                  else
                    rkb_d = zero ! -Wmaybe-uninitialized
                    tkb_d = zero ! -Wmaybe-uninitialized
                  end if
                  beamoff(4,imbb(i))=(rkb_d*(crzb_d-exp_mb(-one*tkb_d)*cbzb_d))*sign(one,crkveb_d)
                  beamoff(5,imbb(i))=(rkb_d*(crxb_d-exp_mb(-one*tkb_d)*cbxb_d))*sign(one,cikveb_d)
                endif
              endif
            endif
          enddo
        endif

!--Crab Cavities
!   Note: If setting something else than el(),
!   DON'T call initialize_element on a crab, it will reset the phase to 0.
      elseif(abs(kz(ix)).eq.23) then
         !Moved from daten()
         crabph(ix)=el(ix)
         el(ix)=zero
!--CC Mult kick order 2
      elseif(abs(kz(ix)).eq.26) then
         !Moved from daten()
         crabph2(ix)=el(ix)
         el(ix)=zero
!--CC Mult kick order 3
      elseif(abs(kz(ix)).eq.27) then
         !Moved from daten()
         crabph3(ix)=el(ix)
         el(ix)=zero
!--CC Mult kick order 4
      else if(abs(kz(ix)).eq.28) then
         !Moved from daten()
         crabph4(ix)=el(ix)
         el(ix)=zero
!--Wire
      else if(kz(ix).eq.15) then
         ed(ix)=zero
         ek(ix)=zero
         el(ix)=zero
!--e-lens
      else if(kz(ix).eq.29) then
         ed(ix)=zero
         ek(ix)=zero
         el(ix)=zero
      endif

      return

      !Error handlers
 100  continue
      write (lout,*) "ERROR in initialize_element, tried to set"
      write (lout,*) "the strength of an element which is disabled."
      write (lout,*) "bez = ", bez(ix)
      call prror(-1)

end subroutine initialize_element

subroutine rounderr(errno,fields,f,value)

  use floatPrecision
  use crcoall

  implicit none

  integer errno,f,l
  character(len=*) fields(*)
  character(len=999) localstr
  real(kind=fPrec) value

  l=len(fields(f))
  localstr=fields(f)(1:l)
  write (lout,"(a,i0)")     "ROUND> ERROR Data Input Error Overfow/Underflow in round_near. Errno: ",errno
  write (lout,"(a,i0,a)")   "ROUND>       f:fieldf:",f,":"//localstr
  write (lout,"(a,e24.16)") "ROUND>       Function fround (rounderr) returning: ",value

  call prror(-1)

end subroutine rounderr

      subroutine wzset
!  *********************************************************************
!
!  This subroutine must be called before subroutine WZSUB can be used to
!  compute values of the complex error function w(z).
!
!  Parameters xcut and ycut specify the opposite corners (xcut,0) and
!  (0,ycut) of the rectangle inside which interpolation is to be used
!  by subroutine WZSUB.
!
!  Parameter h is the side of the squares of the interpolation grid.
!
!  Parameters nx and ny must be set to the nearest integers to xcut/h
!  and ycut/h respectively (or to larger values).
!
!  Calls MYWWERF new version of (CERN library) WWERF (C335)
!
!  (G.A.Erskine, 29.09.1995)
!
!  *********************************************************************
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      use parpro
      use parbeam
      implicit none
      integer i,j,k
      real(kind=fPrec) wi,wr,x,y
      save
!-----------------------------------------------------------------------
      hrecip = one/h
      kstep = nx+2
      k = 0
      do 2 j=0,ny+1
         do 1 i=0,nx+1
            k = k+1
            x=real(i,fPrec)*h                                                  !hr05
            y=real(j,fPrec)*h                                                  !hr05
            call mywwerf(x,y,wr,wi)
            wtreal(k)=wr
            wtimag(k)=wi
 1       continue
 2    continue
      end subroutine wzset

      subroutine mywwerf(x,y,wr,wi)
      use floatPrecision
      use mathlib_bouncer
      use numerical_constants
      implicit none
      integer n
      real(kind=fPrec) c,c1,c2,c3,c4,hf,p,rr,ri,sr0,sr,si,tr,ti,vi,vr,  &
     &wi,wr,x,xa,xl,y,ya,zhi,zhr,z1,z10
      parameter (z1=one,hf=z1/two,z10=c1e1)
      parameter (c1=74.0_fPrec/z10,c2=83.0_fPrec/z10)
      parameter (c3=z10/32.0_fPrec,c4=16.0_fPrec/z10)
!     parameter (c=1.12837916709551257d0,p=(2d0*c4)**33)
      parameter (c=1.12837916709551257_fPrec)
      parameter (p=46768052394588893.3825_fPrec)
      dimension rr(37),ri(37)
      save
!-----------------------------------------------------------------------
      xa=abs(x)
      ya=abs(y)
      if(ya.lt.c1.and.xa.lt.c2) then
!        zh=dcmplx(ya+c4,xa)
        zhr=ya+c4
        zhi=xa
        rr(37)=zero
        ri(37)=zero
        do n=36,1,-1
!          t=zh+n*dconjg(r(n+1))
          tr=zhr+real(n,fPrec)*rr(n+1)                                         !hr05
          ti=zhi-real(n,fPrec)*ri(n+1)                                         !hr05
!          r(n)=hf*t/(dreal(t)**2+dimag(t)**2)
          rr(n)=(hf*tr)/(tr**2+ti**2)                                    !hr05
          ri(n)=(hf*ti)/(tr**2+ti**2)                                    !hr05
        enddo
        xl=p
        sr=zero
        si=zero
        do n=33,1,-1
          xl=c3*xl
!          s=r(n)*(s+xl)
          sr0=rr(n)*(sr+xl)-ri(n)*si
          si=rr(n)*si+ri(n)*(sr+xl)
          sr=sr0
        enddo
!        v=c*s
        vr=c*sr
        vi=c*si
      else
        zhr=ya
        zhi=xa
        rr(1)=zero
        ri(1)=zero
        do n=9,1,-1
!          t=zh+n*dconjg(r(1))
          tr=zhr+real(n,fPrec)*rr(1)                                           !hr05
          ti=zhi-real(n,fPrec)*ri(1)                                           !hr05
!          r(1)=hf*t/(dreal(t)**2+dimag(t)**2)
          rr(1)=(hf*tr)/(tr**2+ti**2)                                    !hr05
          ri(1)=(hf*ti)/(tr**2+ti**2)                                    !hr05
        enddo
!        v=c*r(1)
        vr=c*rr(1)
        vi=c*ri(1)
      endif
      if(ya.eq.zero) then                                                 !hr05
!        v=dcmplx(exp(-xa**2),dimag(v))
        vr=exp_mb(-one*xa**2)                                            !hr05
      endif
      if(y.lt.zero) then
!        v=2*exp(-dcmplx(xa,ya)**2)-v
        vr=(two*exp_mb(ya**2-xa**2))*cos_mb((two*xa)*ya)-vr              !hr05
        vi=(-two*exp_mb(ya**2-xa**2))*sin_mb((two*xa)*ya)-vi             !hr05
        if(x.gt.zero) vi=-one*vi                                          !hr05
      else
        if(x.lt.zero) vi=-one*vi                                          !hr05
      endif
      wr=vr
      wi=vi
      return
end

!-----------------------------------------------------------------------
!  CALCULATION OF : MOMENTUM-DEPENDING ELEMENT-MATRICES AND
!                   CHANGE OF PATH LENGTHS FOR EACH PARTICLE.
!-----------------------------------------------------------------------
subroutine envars(j,dpp,rv)

  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use parpro
  use mod_common
  use mod_commons
  use mod_commont
  use mod_commond
  implicit none
  integer i,ih,j,kz1,l,ll
  real(kind=fPrec) aek,afok,as3,as4,as6,co,dpd,dpp,dpsq,fi,fok,fok1,fokq,g,gl,hc,hi,hi1,hm,&
                   hp,hs,rho,rhoc,rhoi,rv,si,siq,sm1,sm12,sm2,sm23,sm3,sm5,sm6,wf,wfa,wfhi
  integer jm,k,m,na,ne
  save

  dpd  = one+dpp
  dpsq = sqrt(dpd)
  do i=1,il
    do ll=1,6
      do l=1,2
        al(ll,l,j,i) = zero
        as(ll,l,j,i) = zero
        at(ll,l,j,i) = zero
        a2(ll,l,j,i) = zero
      end do
    end do

    if(abs(el(i)).le.pieni) cycle
    kz1 = kz(i)+1
    select case (kz1)
    case (1) ! DRIFTLENGTH
      do l=1,2
        al(1,l,j,i) = one
        al(2,l,j,i) = el(i)
        al(3,l,j,i) = zero
        al(4,l,j,i) = one
        as(6,l,j,i) = ((-one*rv)*al(2,l,j,i))/c2e3 ! hr05
      end do
      as(1,1,j,i) = (el(i)*(one-rv))*c1e3          ! hr05

    case (2,5)
      ! 2: RECTANGULAR MAGNET HORIZONTAL
      ! 5: RECTANGULAR MAGNET VERTICAL
      if (kz1.eq.2) then
        ih = 1
      else
        ih = 2
      end if
      fok=(el(i)*ed(i))/dpsq                                           !hr05
      if(abs(fok).le.pieni) then
        do l=1,2
          al(1,l,j,i) = one
          al(2,l,j,i) = el(i)
          al(3,l,j,i) = zero
          al(4,l,j,i) = one
          as(6,l,j,i) = ((-one*rv)*al(2,l,j,i))/c2e3 ! hr05
        end do
        as(1,1,j,i) = (el(i)*(one-rv))*c1e3          ! hr05
        cycle
      end if
      rho=(one/ed(i))*dpsq
      fok1=(tan_mb(fok*half))/rho
      si=sin_mb(fok)
      co=cos_mb(fok)
      al(1,ih,j,i)=one
      al(2,ih,j,i)=rho*si
      al(3,ih,j,i)=zero
      al(4,ih,j,i)=one
      al(5,ih,j,i)=((-one*dpp)*((rho*(one-co))/dpsq))*c1e3             !hr05
      al(6,ih,j,i)=((-one*dpp)*((two*tan_mb(fok*half))/dpsq))*c1e3     !hr05
      sm1=cos_mb(fok)
      sm2=sin_mb(fok)*rho
      sm3=(-one*sin_mb(fok))/rho                                       !hr05
      sm5=((-one*rho)*dpsq)*(one-sm1)                                  !hr05
      sm6=((-one*sm2)*dpsq)/rho                                        !hr05
      sm12=el(i)-sm1*sm2
      sm23=sm2*sm3
      as3=(-one*rv)*(((dpp*rho)/(two*dpsq))*sm23+sm5)                  !hr05
      as4=((-one*rv)*sm23)/c2e3                                        !hr05
      as6=((-one*rv)*(el(i)+sm1*sm2))/c4e3                             !hr05
      as(1,ih,j,i)=(el(i)*(one-rv)-rv*((dpp**2/(four*dpd))*sm12+dpp*(el(i)-sm2)))*c1e3 !hr05
      as(2,ih,j,i)=fok1*as3-rv*((dpp/((two*rho)*dpsq))*sm12+sm6)       !hr05
      as(3,ih,j,i)=as3
      as(4,ih,j,i)=as4+(two*as6)*fok1                                  !hr05
      as(5,ih,j,i)=(as6*fok1*2+fok1*as4)-(rv*sm12)/(c4e3*rho**2)       !hr05
      as(6,ih,j,i)=as6
      ! VERTIKAL
      ih=ih+1
      if(ih.gt.2) ih=1
      g=tan_mb(fok*half)/rho
      gl=el(i)*g
      al(1,ih,j,i)=one-gl
      al(2,ih,j,i)=el(i)
      al(3,ih,j,i)=(-one*g)*(two-gl)                                   !hr05
      al(4,ih,j,i)=al(1,ih,j,i)
      as6=((-one*rv)*al(2,ih,j,i))/c2e3                                !hr05
      as(4,ih,j,i)=((-one*two)*as6)*fok1                               !hr05
      as(5,ih,j,i)=as6*fok1**2                                         !hr05
      as(6,ih,j,i)=as6

    case (3)
      ! QUADRUPOLE
      ! FOCUSSING
      fok=ek(i)/(one+dpp)
      aek=abs(fok)
      if(abs(fok).le.pieni) then
        do l=1,2
          al(1,l,j,i) = one
          al(2,l,j,i) = el(i)
          al(3,l,j,i) = zero
          al(4,l,j,i) = one
          as(6,l,j,i) = ((-one*rv)*al(2,l,j,i))/c2e3 ! hr05
        end do
        as(1,1,j,i) = (el(i)*(one-rv))*c1e3          ! hr05
        cycle
      end if
      ih=0
      hi=sqrt(aek)
      fi=el(i)*hi
      if(fok.gt.zero) goto 110
100   ih=ih+1
      al(1,ih,j,i)=cos_mb(fi)
      hi1=sin_mb(fi)
      al(2,ih,j,i)=hi1/hi
      al(3,ih,j,i)=(-one*hi1)*hi                                       !hr05
      al(4,ih,j,i)=al(1,ih,j,i)
      as(1,ih,j,i)=(el(i)*(one-rv))*c1e3                               !hr05
      as(4,ih,j,i)=(((-one*rv)*al(2,ih,j,i))*al(3,ih,j,i))/c2e3        !hr05
      as(5,ih,j,i)=(((-one*rv)*(el(i)-al(1,ih,j,i)*al(2,ih,j,i)))*aek)/c4e3 !hr05
      as(6,ih,j,i)=((-one*rv)*(el(i)+al(1,ih,j,i)*al(2,ih,j,i)))/c4e3  !hr05
      if(ih.eq.2) cycle
      !--DEFOCUSSING
110   ih=ih+1
      hp=exp_mb(fi)
      hm=one/hp
      hc=(hp+hm)*half
      hs=(hp-hm)*half
      al(1,ih,j,i)=hc
      al(2,ih,j,i)=hs/hi
      al(3,ih,j,i)=hs*hi
      al(4,ih,j,i)=hc
      as(4,ih,j,i)=(((-one*rv)*al(2,ih,j,i))*al(3,ih,j,i))/c2e3        !hr05
      as(5,ih,j,i)=((rv*(el(i)-al(1,ih,j,i)*al(2,ih,j,i)))*aek)/c4e3   !hr05
      as(6,ih,j,i)=((-one*rv)*(el(i)+al(1,ih,j,i)*al(2,ih,j,i)))/c4e3  !hr05
      if(ih.eq.1) goto 100

    case (4)
      ! 4: SEKTORMAGNET HORIZONTAL
      ! 6: SEKTORMAGNET VERTICAL
      if (kz1.eq.4) then
        ih = 1
      else
        ih = 2
      end if
      fok=(el(i)*ed(i))/dpsq                                           !hr05
      if(abs(fok).le.pieni) then
        do l=1,2
          al(1,l,j,i) = one
          al(2,l,j,i) = el(i)
          al(3,l,j,i) = zero
          al(4,l,j,i) = one
          as(6,l,j,i) = ((-one*rv)*al(2,l,j,i))/c2e3 ! hr05
        end do
        as(1,1,j,i) = (el(i)*(one-rv))*c1e3          ! hr05
        cycle
      end if
      rho=(one/ed(i))*dpsq
      si=sin_mb(fok)
      co=cos_mb(fok)
      rhoc=(rho*(one-co))/dpsq                                         !hr05
      siq=si/dpsq
      al(1,ih,j,i)=co
      al(2,ih,j,i)=rho*si
      al(3,ih,j,i)=(-one*si)/rho                                       !hr05
      al(4,ih,j,i)=co
      al(5,ih,j,i)=((-one*dpp)*rhoc)*c1e3                              !hr05
      al(6,ih,j,i)=((-one*dpp)*siq)*c1e3                               !hr05
      sm12=el(i)-al(1,ih,j,i)*al(2,ih,j,i)
      sm23=al(2,ih,j,i)*al(3,ih,j,i)
      as(1,ih,j,i)=(el(i)*(one-rv)-rv*((dpp**2/(four*dpd))*sm12+dpp*(el(i)-al(2,ih,j,i))))*c1e3 !hr05
      as(2,ih,j,i)=(-one*rv)*((dpp/((two*rho)*dpsq))*sm12-dpd*siq)     !hr05
      as(3,ih,j,i)=(-one*rv)*(((dpp*rho)/(two*dpsq))*sm23-dpd*rhoc)    !hr05
      as(4,ih,j,i)=((-one*rv)*sm23)/c2e3                               !hr05
      as(5,ih,j,i)=((-one*rv)*sm12)/(c4e3*rho**2)                      !hr05
      as(6,ih,j,i)=((-one*rv)*(el(i)+al(1,ih,j,i)*al(2,ih,j,i)))/c4e3  !hr05
      ! VERTIKAL
      ih=ih+1
      if(ih.gt.2) ih=1
      al(1,ih,j,i)=one
      al(2,ih,j,i)=el(i)
      al(3,ih,j,i)=zero
      al(4,ih,j,i)=one
      as(6,ih,j,i)=((-one*rv)*al(2,ih,j,i))/c2e3                       !hr05

    case (7,8)
      ! 7: COMBINED FUNCTION MAGNET HORIZONTAL
      ! 8: COMBINED FUNCTION MAGNET VERTICAL
      if (kz1.eq.7) then
        ih   = 0
        fokq = ek(i)
      else
        ih   = 1
        fokq = -one*ek(i) ! hr05
      end if
      ! FOCUSSING
      wf=ed(i)/dpsq
      fok=fokq/(dpd)-wf**2                                             !hr05
      if(abs(fok).le.pieni) then
        do l=1,2
          al(1,l,j,i) = one
          al(2,l,j,i) = el(i)
          al(3,l,j,i) = zero
          al(4,l,j,i) = one
          as(6,l,j,i) = ((-one*rv)*al(2,l,j,i))/c2e3 ! hr05
        end do
        as(1,1,j,i) = (el(i)*(one-rv))*c1e3          ! hr05
        cycle
      end if
      afok=abs(fok)
      hi=sqrt(afok)
      fi=hi*el(i)
      if(fok.gt.zero) goto 160
140   ih=ih+1
      si=sin_mb(fi)
      co=cos_mb(fi)
      wfa=((wf/afok)*(one-co))/dpsq                                    !hr05
      wfhi=((wf/hi)*si)/dpsq                                           !hr05
      al(1,ih,j,i)=co
      al(2,ih,j,i)=si/hi
      al(3,ih,j,i)=(-one*si)*hi                                        !hr05
      al(4,ih,j,i)=co
      al(5,ih,j,i)=((-one*wfa)*dpp)*c1e3                               !hr05
      al(6,ih,j,i)=((-one*wfhi)*dpp)*c1e3                              !hr05
      sm12=el(i)-al(1,ih,j,i)*al(2,ih,j,i)
      sm23=al(2,ih,j,i)*al(3,ih,j,i)
      as(1,ih,j,i)=(el(i)*(one-rv)-((rv*((dpp**2/(four*dpd))*sm12+dpp*(el(i)-al(2,ih,j,i))))/afok)*wf**2)*c1e3 !hr05
      as(2,ih,j,i)=(-one*rv)*(((dpp*wf)/(two*dpsq))*sm12-dpd*wfhi)     !hr05
      as(3,ih,j,i)=(-one*rv)*(((((dpp*half)/afok)/dpd)*ed(i))*sm23-dpd*wfa) !hr05
      as(4,ih,j,i)=((-one*rv)*sm23)/c2e3                               !hr05
      as(5,ih,j,i)=(((-one*rv)*sm12)*afok)/c4e3                        !hr05
      as(6,ih,j,i)=((-one*rv)*(el(i)+al(1,ih,j,i)*al(2,ih,j,i)))/c4e3  !hr05
      ih=ih+1
      if(ih.gt.2) ih=1
      aek=abs(ek(i)/dpd)
      hi=sqrt(aek)
      fi=hi*el(i)
      hp=exp_mb(fi)
      hm=one/hp
      hc=(hp+hm)*half
      hs=(hp-hm)*half
      al(1,ih,j,i)=hc
      al(2,ih,j,i)=el(i)
      if(abs(hi).le.pieni) goto 150
      al(2,ih,j,i)=hs/hi
150   al(3,ih,j,i)=hs*hi
      al(4,ih,j,i)=hc
      as(4,ih,j,i)=(((-one*rv)*al(2,ih,j,i))*al(3,ih,j,i))/c2e3        !hr05
      as(5,ih,j,i)=((rv*(el(i)-al(1,ih,j,i)*al(2,ih,j,i)))*aek)/c4e3   !hr05
      as(6,ih,j,i)=((-one*rv)*(el(i)+al(1,ih,j,i)*al(2,ih,j,i)))/c4e3  !hr05
      cycle
      ! DEFOCUSSING
160   ih=ih+1
      hp=exp_mb(fi)
      hm=one/hp
      hc=(hp+hm)*half
      hs=(hp-hm)*half
      al(1,ih,j,i)=hc
      al(2,ih,j,i)=hs/hi
      al(3,ih,j,i)=hs*hi
      al(4,ih,j,i)=hc
      wfa=((wf/afok)*(one-hc))/dpsq                                    !hr05
      wfhi=((wf/hi)*hs)/dpsq                                           !hr05
      al(5,ih,j,i)= (wfa*dpp)*c1e3                                     !hr05
      al(6,ih,j,i)=((-one*wfhi)*dpp)*c1e3                              !hr05
      sm12=el(i)-al(1,ih,j,i)*al(2,ih,j,i)
      sm23=al(2,ih,j,i)*al(3,ih,j,i)
      as(1,ih,j,i)=(((rv*((dpp**2/(four*dpd))*sm12+dpp*(el(i)-al(2,ih,j,i))))/afok)*wf**2+el(i)*(one-rv))*c1e3 !hr05
      as(2,ih,j,i)=(-one*rv)*(((dpp*wf)/(two*dpsq))*sm12-dpd*wfhi)     !hr05
      as(3,ih,j,i)=rv*(((((dpp*half)/afok)/dpd)*ed(i))*sm23-dpd*wfa)   !hr05
      as(4,ih,j,i)=((-one*rv)*sm23)/c2e3                               !hr05
      as(5,ih,j,i)=((rv*sm12)*afok)/c4e3                               !hr05
      as(6,ih,j,i)=((-one*rv)*(el(i)+al(1,ih,j,i)*al(2,ih,j,i)))/c4e3  !hr05
      ih=ih+1
      if(ih.gt.2) ih=1
      aek=abs(ek(i)/dpd)
      hi=sqrt(aek)
      fi=hi*el(i)
      si=sin_mb(fi)
      co=cos_mb(fi)
      al(1,ih,j,i)=co
      al(2,ih,j,i)=si/hi
      al(3,ih,j,i)=(-one*si)*hi                                        !hr05
      al(4,ih,j,i)=co
      as(4,ih,j,i)=(((-one*rv)*al(2,ih,j,i))*al(3,ih,j,i))/c2e3        !hr05
      as(5,ih,j,i)=(((-one*rv)*(el(i)-al(1,ih,j,i)*al(2,ih,j,i)))*aek)/c4e3 !hr05
      as(6,ih,j,i)=((-one*rv)*(el(i)+al(1,ih,j,i)*al(2,ih,j,i)))/c4e3  !hr05

    case (9)
      ! EDGE FOCUSSING
      rhoi=ed(i)/dpsq
      fok=rhoi*tan_mb((el(i)*rhoi)*half)                               !hr05
      al(1,1,j,i)=one
      al(2,1,j,i)=zero
      al(3,1,j,i)=fok
      al(4,1,j,i)=one
      al(1,2,j,i)=one
      al(2,2,j,i)=zero
      al(3,2,j,i)=-one*fok                                             !hr05
      al(4,2,j,i)=one

    end select
  end do

  do k=1,mblo
    jm=mel(k)
    do m=1,jm
      na=mtyp(k,m)
      ne=mtyp(k,jm-m+1)
      do l=1,2
        at(1,l,j,na)=as(1,l,j,ne)
        at(2,l,j,na)=as(2,l,j,ne)
        at(3,l,j,na)=as(3,l,j,ne)
        at(4,l,j,na)=as(4,l,j,ne)
        at(5,l,j,na)=as(5,l,j,ne)
        at(6,l,j,na)=as(6,l,j,ne)
        a2(1,l,j,na)=al(1,l,j,ne)
        a2(2,l,j,na)=al(2,l,j,ne)
        a2(3,l,j,na)=al(3,l,j,ne)
        a2(4,l,j,na)=al(4,l,j,ne)
        a2(5,l,j,na)=al(5,l,j,ne)
        a2(6,l,j,na)=al(6,l,j,ne)
      end do
    end do
  end do
  return

end subroutine envars

!-----------------------------------------------------------------------
!  SUBROUTINE TO SET THE ALL COMMON VARIABLES TO ZERO
!-----------------------------------------------------------------------
subroutine comnul

  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  use parpro
  use parbeam, only : beam_expflag,beam_expfile_open
  use mod_common
  use mod_commonmn
  use mod_commons
  use mod_commont
  use mod_commond

  use aperture
  use elens
  use wire
  use scatter,     only : scatter_comnul
  use dump,        only : dump_comnul
  use bdex,        only : bdex_comnul
  use collimation, only : collimation_comnul
#ifdef HDF5
  use hdf5_output, only : h5_comnul
#endif

  implicit none

  integer i,i1,i2,i3,i4,j

  save

  ! CHROMATICITY ADJUSTMENT BLOCK
  ichrom  = 0     ! mod_common

  ! SUBRESONANCE CALCULATION BLOCK
  isub    = 0     ! mod_common
  nta     = 0     ! mod_common
  nte     = 0     ! mod_common
  ipt     = 0     ! mod_common

  ! RESONANCE COMPENSATION BLOCK
  irmod2  = 0     ! mod_common
  nre     = 0     ! mod_common

  ! TODO
      nur=0
      nch=0
      nqc=0
      npp=0
      ipos=0
      iconv=0
      imad=0
      nstart=0
      nstop=0
      iskip=1
      iav=0
      iwg=0
      ivox=0
      ivoz=0
      ires=0
      ifh=0
      idis=0
      icow=0
      istw=0
      iffw=0
      idial=0
      nord=0
      nvar=0
      nvar2=0
      ndimf=0
      nordf=0
      nvarf=0
      nord1=1
      nsix=0
      nvar2=0
      ncor=0
      idptr=0
      ibb6d=0
!-----------------------------------------------------------------------
      inorm=0
      imod1=0
      imod2=0
!-----------------------------------------------------------------------
      dphix=zero
      dphiz=zero
      qx0=zero
      qz0=zero
      dres=zero
      dfft=zero
      preda=zero
      damp=zero
      ampt=zero
!-----------------------------------------------------------------------
      do 10 i=1,2
        is(i)=0
        bet0(i)=zero
        alf0(i)=zero
        clo(i)=zero
        clop(i)=zero
        cro(i)=zero
        qwsk(i)=zero
        betx(i)=zero
        betz(i)=zero
        alfx(i)=zero
        alfz(i)=zero
   10 continue

      do 20 i=1,3
        clon(i)=zero
        wxys(i)=zero
        do i1=1,3
          corr(i,i1)=zero
        enddo
   20 continue

      corr(1,1)=zero
      corr(1,2)=zero
      chromc(1)=9.999999e23_fPrec
      chromc(2)=9.999999e23_fPrec

      do 40 i=1,5
        ipr(i)=0
        nrr(i)=0
        nu(i)=0
   40 continue

      do 50 i=1,6
        nskew(i)=0
   50 continue

      do 60 i=1,10
        dtr(i)=zero
   60 continue

      do 70 i=1,12
        ire(i)=0
   70 continue

      do i1=1,9
        do i2=1,18
          do i3=1,10
            do i4=1,5
              rtc(i1,i2,i3,i4)=zero
              rts(i1,i2,i3,i4)=zero
            end do
          end do
        end do
      end do

!--NUMBER OF ELEMENTS---------------------------------------------------
      do i=1,nele
        dki(i,1)=zero
        dki(i,2)=zero
        dki(i,3)=zero
        acdipph(i)=zero
      end do

!     From the FLUKA version
      do i=1,nele
         call SELNUL(i)
      end do

!-- BEAM-EXP------------------------------------------------------------
      beam_expflag = 0
      beam_expfile_open = .false.

!--# OF TRAJECTORIES----------------------------------------------------
      do 220 i=1,mpa
        do 210 i1=1,2
          x(i,i1)=zero
          y(i,i1)=zero
  210   continue
  220 continue

!--PAW------------------------------------------------------------------
      do i=1,nplo
        hmal(i)=0.0
      end do

!--TROMBONES------------------------------------------------------------
      do i=1,ntr
        do i1=1,6
          cotr(i,i1)=zero
          do i2=1,6
            rrtr(i,i1,i2)=zero
          end do
        end do
      end do

!--DUMP BEAM POPULATION-------------------------------------------------
!     A.Mereghetti, D.Sinuela Pastor and P.Garcia Ortega, for the FLUKA Team
!     K.Sjobak, BE-ABP/HSS
!     last modified: 03-09-2015
!     initialise common
!     always in main code
      call dump_comnul

!--ELEN - ELECTRON LENS---------------------------------------------------------
!     M. Fitterer (FNAL), A. Mereghetti
!     last modified: 09-02-2018
!     always in main code
!     elensparam - used for tracking (parameters of single element)
      do i=1,nele
         ielens(i) = 0
      end do
      melens=0
      do i=1,nelens
        elens_type(i)          = 0
        elens_theta_r2(i)      = zero
        elens_r2(i)            = zero
        elens_r1(i)            = zero
        elens_offset_x(i)      = zero
        elens_offset_y(i)      = zero
        elens_sig(i)           = zero
        elens_geo_norm(i)      = zero
        elens_len(i)           = zero
        elens_I(i)             = zero
        elens_Ek(i)            = zero
        elens_lThetaR2(i)      = .false.
        elens_lAllowUpdate(i)  = .true.
        elens_iCheby(i)        = 0
        elens_cheby_angle(i)   = zero
      end do
!     table with coefficients of chebyshev polynominals
      do i=1,nelens_cheby_tables
         do j=1,16
            elens_cheby_filename(i)(j:j)=' '
         end do
         do i1=0,elens_cheby_order
            do i2=0,elens_cheby_order
               elens_cheby_coeffs(i1,i2,i)=zero
            end do
         end do
         elens_cheby_refCurr(i)=zero
         elens_cheby_refBeta(i)=zero
         elens_cheby_refRadius(i)=zero
      end do
!--WIRE - WIRE ELEMENT---------------------------------------------------------
!     M. Fitterer (FNAL), A. Patapenka (NIU)
!     last modified: 22-12-2016
! 1)  wireparam - used for tracking (parameters of single element)
      do i=1,nele
        wire_flagco(i)  = 0
        wire_current(i) = 0
        wire_lint(i)    = 0
        wire_lphys(i)   = 0
        wire_dispx(i)   = 0
        wire_dispy(i)   = 0
        wire_tiltx(i)   = 0
        wire_tilty(i)   = 0
      end do

! 2) loop over structure elements
      do i=1,nblz
        wire_num(i)=0
      end do

! 3) loop over number of wires
      do i=1,wire_max
        do j=1,6
          wire_clo(j,i)=zero
        end do
      end do

!--APERTURE-------------------------------------------------------------
!     P.G.Ortega and A.Mereghetti, for the FLUKA Team
!     last modified: 02-03-2018
!     initialise common
!     always in main code
      call aperture_comnul

!--SCATTER--------------------------------------------------------------
      call scatter_comnul
!--HDF5-----------------------------------------------------------------
#ifdef HDF5
      call h5_comnul
#endif
!--COLLIMATION----------------------------------------------------------
      call collimation_comnul
!--BDEX-----------------------------------------------------------------
      call bdex_comnul
!-----------------------------------------------------------------------
      return
end subroutine comnul

subroutine SELNUL( iel )
!-----------------------------------------------------------------------
!     A.Mereghetti, 2016-03-14
!     initialise a single element to empty
!-----------------------------------------------------------------------
      use floatPrecision
      use numerical_constants
      use parpro
      use mod_common
      use mod_commons
      use aperture
      implicit none

!     local variables
      integer iel, i1, i3, i4

      kz(iel)=0
      kp(iel)=0
      irm(iel)=0
      imtr(iel)=0
      nmu(iel)=0
      kpa(iel)=0
      isea(iel)=0
      ncororb(iel)=0
      iratioe(iel)=0
      itionc(iel)=0
      dki(iel,1)=zero
      dki(iel,2)=zero
      dki(iel,3)=zero
      ed(iel)=zero
      el(iel)=zero
      ek(iel)=zero
      sm(iel)=zero
      xpl(iel)=zero
      xrms(iel)=zero
      zpl(iel)=zero
      zrms(iel)=zero
      benkc(iel)=zero
      r00(iel)=zero

      ratioe(iel)=one
      hsyc(iel)=zero
      phasc(iel)=zero
      ptnfac(iel)=zero
!      wirel(iel)=zero
      acdipph(iel)=zero
      crabph(iel)=zero
      crabph2(iel)=zero
      crabph3(iel)=zero
      crabph4(iel)=zero
      bez(iel)=' '
      bezl(iel)=' '

      do i3=1,2
        do i4=1,6
          a(iel,i3,i4)=zero
        end do
      end do

      do i1=1,mmul
         bk0(iel,i1)=zero
         ak0(iel,i1)=zero
         bka(iel,i1)=zero
         aka(iel,i1)=zero
      end do

      do i1=1,3
         bezr(i1,iel)=' '
      end do

      ! JBG increasing parbe to dimension 5
      do i1=1,5
         parbe(iel,i1)=zero
      end do

      call aperture_nul(iel)

      return
end subroutine SELNUL

subroutine STRNUL( iel )
!-----------------------------------------------------------------------
!     A.Mereghetti, 2016-03-14
!     initialise an element in lattice structure to empty
!-----------------------------------------------------------------------
      use floatPrecision
      use numerical_constants
      use parpro
      use mod_common
      use mod_commonmn
      implicit none

!     local variables
      integer iel, i1, i2, i3

      ic(iel)=0
      mzu(iel)=0
      icext(iel)=0
      icextal(iel)=0
      ! extalign(iel,1)=zero
      ! extalign(iel,2)=zero
      ! extalign(iel,3)=zero
      sigmoff(iel)=zero
      tiltc(iel)=one
      tilts(iel)=zero

!--Beam-Beam------------------------------------------------------------
      imbb(iel)=0
      ! do j=1,40
      !    exterr(iel,j)=zero
      ! enddo
      xsi(iel)=zero
      zsi(iel)=zero
      smi(iel)=zero
      smizf(iel)=zero
      do i3=1,mmul
          aaiv(i3,iel)=zero
          bbiv(i3,iel)=zero
      enddo
      return
end subroutine STRNUL

! ================================================================================================ !
!     by A.Mereghetti
!     last modified: 01-12-2016
!     Insert a New Empty Element in Lattice Structure
!     interface variables:
!     - iEl: index in lattice structure where to insert the element
!     always in main code
! ================================================================================================ !
integer function INEELS( iEl )

  use floatPrecision
  use numerical_constants
  use crcoall
  use parpro
  use mod_common
  use mod_commont
  use mod_commonmn
  implicit none

!     interface variables
  integer iEl

!     temporary variables
  integer i,iInsert

  if ( iu.gt.nblz-3) then
    call expand_arrays(nele, npart, nblz+100, nblo)
  end if
  iu=iu+1
  if ( iEl.eq.0 ) then
!        append
    iInsert=iu
  elseif ( iEl .lt. 0 ) then
    iInsert=iu+iEl
  else
    iInsert=iEl
  end if
  if ( iInsert.le.iu ) then
!     shift by one all lattice elements, to make room for the new
!        starting marker
    do i=iu,iInsert+1,-1
      ic(i)=ic(i-1)
      icext(i)=icext(i-1)
      icextal(i)=icextal(i-1)
      dcum(i)=dcum(i-1)
    enddo
  endif

!     initialise element to empty
  call STRNUL(iInsert)
!     update dcum of added element
  dcum(iInsert)=dcum(iInsert-1)
!     return iu
  INEELS=iu
end function INEELS

integer function INEESE()
!-----------------------------------------------------------------------
!     by A.Mereghetti
!     last modified: 01-12-2016
!     Insert a New Empty Element (empty) in SINGLE ELEMENTS
!     for the moment, it only appends the new single element
!     always in main code
!-----------------------------------------------------------------------
  use floatPrecision
  use numerical_constants
  use crcoall

  use parpro
  use mod_common
  use mod_commont
  use mod_commonmn
  implicit none

  il=il+1
  if( il.gt.nele-2 ) then
    call expand_arrays(nele+50, npart, nblz, nblo )
    if (ithick.eq.1) then
      call expand_thickarrays(nele, npart, nblz, nblo )
    end if
  end if

! initialise element to empty
  call SELNUL(il)

! returned variable
  INEESE=il

end function INEESE

integer function check_SE_unique( iEl, ixEl )
!-----------------------------------------------------------------------
!     by A.Mereghetti
!     last modified: 01-12-2016
!     check that a given entry in the sequence is unique
!     interface variables:
!     - iEl:  index in lattice structure to be checked
!     - ixEl: index in array of SINGLE ELEMENTs of the element to be checked
!     always in main code
!-----------------------------------------------------------------------
  use floatPrecision
  use numerical_constants

  use parpro
  use mod_common
  use mod_commont
  use mod_commonmn
  implicit none

! interface variables
  integer iEl, ixEl

! temporary variables
  integer i,ix

  check_SE_unique=-1

  do i=1,iu
    ix=ic(i)-nblo
    if(ix.gt.0) then
!     SINGLE ELEMENT
      if( i.ne.iEl .and. ix.eq.ixEl ) then
        check_SE_unique=i
        exit
      end if
    end if
  end do

  return

end function check_SE_unique

subroutine find_entry_at_s( sLoc, llast, iEl, ixEl, lfound )
!-----------------------------------------------------------------------
!     by A.Mereghetti (CERN, BE/ABP-HSS), 2018-03-22
!     find the element at location sLoc
!     interface variables:
!     - sLoc: s-coordinate where element should be
!     - iEl:  index in lattice structure of found element
!     - ixEl: index in array of SINGLE ELEMENTs of found element
!     - llast: if true, return last lens at sLoc
!     always in main code
!-----------------------------------------------------------------------
  use floatPrecision
  use mod_common     ! for iu, tlen, ic
  use parpro
  use crcoall
  use numerical_constants, only : zero

  implicit none

! interface variables
  integer iEl, ixEl
  real(kind=fPrec) sLoc
  logical llast, lfound

! temporary variables
  integer i, iDelta, iCheck, iMax, iStep
  logical lSlide

  iEl=-1
  ixEl=-1

  if (sLoc.gt.tlen.or.sLoc.lt.zero) then
     write(lout,'(a,f11.4,a,f11.4,a)')&
 &         "requested s-location: ",sLoc," is not inside accelerator range: [0:",tlen,"]"
     call prror(-1)
  endif

  ! fast search
  iCheck=iu
  iDelta=iu
  do while(iDelta.gt.1.or.iCheck.gt.0.or.iCheck.lt.iu)
     if (dcum(iCheck).eq.sLoc) exit
     iDelta=nint(real(iDelta/2))
     if (dcum(iCheck).lt.sLoc) then
        iCheck=iCheck+iDelta
     else
        iCheck=iCheck-iDelta
     endif
  end do

  ! finalise search
  if (dcum(iCheck).lt.sLoc.or.(dcum(iCheck).eq.sLoc.and.llast)) then
     iMax=iu
     iStep=1
     lslide=llast
  else
     iMax=1
     iStep=-1
     lSlide=.not.llast
  endif
  do i=iCheck,iMax,iStep
     if (dcum(i).lt.sLoc) continue
     if (dcum(i).gt.sLoc) then
        if (iEl.eq.-1) iEl=i
     else
        iEl=i
        if (lSlide) continue
     endif
     exit
  enddo

  if (lfound) then
     ixEl=ic(iEl)-nblo
     if (ixEl.lt.0) ixEl=ic(iEl) ! drift
  else
     write(lout,'(a,f11.4,a,f11.4,a)')&
 &         "s-location: ",sLoc," was not found in acclerator range: [0:",tlen,"]"
     write(lout,*)"this is actually embarassing..."
  endif
  return

end subroutine find_entry_at_s

subroutine distance(x,clo,di0,t,dam)
!-----------------------------------------------------------------------
!  CALCULATION OF DISTANCE IN PHASE SPACE FOR POST-PROCESSING
!-----------------------------------------------------------------------
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      use parpro
      use mod_common, only : dpscor,sigcor,icode,idam,its6d
      implicit none
      integer i,ii,iq,j,jq
      real(kind=fPrec) clo,cx,dam,di0,phi,sx,t,x,x1
      dimension x(2,6),x1(2,6),clo(6),di0(4),t(6,6),phi(3)
      save
!-----------------------------------------------------------------------
      if(icode.ge.4.and.its6d.eq.0) then
        do i=1,2
          do j=1,4
            x(i,j)=x(i,j)-di0(j)*x(i,6)
          end do
        end do
      endif

      do 60 i=1,2
        do 20 j=1,6
          x(i,j)=x(i,j)-clo(j)
   20   continue

        if(its6d.eq.1) then
          x(i,2)=x(i,2)/((one+x(i,6))+clo(6))                            !hr06
          x(i,4)=x(i,4)/((one+x(i,6))+clo(6))                            !hr06
        endif

        do 40 iq=1,6
          x1(i,iq)=zero
          do 30 jq=1,6
            x1(i,iq)=x1(i,iq)+t(jq,iq)*x(i,jq)
   30     continue
   40   continue

        do 50 j=1,6
          x(i,j)=x1(i,j)
   50   continue

   60 continue

      do 70 i=1,2
        x(i,5)=x(i,5)*sigcor
        x(i,6)=x(i,6)*dpscor
   70 continue

      do 80 i=1,3
        ii=2*i
        sx=x(2,ii-1)*x(1,ii)-x(1,ii-1)*x(2,ii)
        cx=x(1,ii-1)*x(2,ii-1)+x(1,ii)*x(2,ii)
        if(abs(sx).gt.c1m15.or.abs(cx).gt.c1m15) then
          phi(i)=atan2_mb(sx,cx)
        else
          phi(i)=zero
        endif
   80 continue
      dam=sqrt((phi(1)**2+phi(2)**2+phi(3)**2)/real(idam,fPrec))/pi            !hr06
!-----------------------------------------------------------------------
      return
end subroutine distance

subroutine betalf(dpp,qw)
!-----------------------------------------------------------------------
!  CALCULATION OF : OPT. PARAMETERS AT THE STARTING POSITION:
!                   BETA-, ALFA-FUNCTIONS, Q-VALUES
!-----------------------------------------------------------------------
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      use parpro
      use mod_common
      use mod_commons
      use mod_commont
      implicit none
      integer i,j
      real(kind=fPrec) am,det,detb,detc,dpp,egwg1,egwg2,f0,f1,f2,fak1,  &
     &fak2,qw,rca1,rca2,rclam1,rclam2,rcw1(4),rcw2(4),rn1,rn2,spa,spd,  &
     &sqrn,yca1,yca2,yclam1,yclam2,ycw1(4),ycw2(4)
      dimension am(4,4)
      dimension qw(2)
      save
!-----------------------------------------------------------------------
      ierro=0
      call matrix(dpp,am)
!--CALCULATION OF EIGENVALUES
   10 spa=am(1,1)+am(2,2)
      spd=am(3,3)+am(4,4)
      det=(am(1,3)+am(4,2))*(am(2,4)+am(3,1))                           &
     &-(am(1,4)-am(3,2))*(am(2,3)-am(4,1))
      f0=spa-spd
      f1=spa+spd
      f2=f0**2+four*det                                                  !hr06
      if(f2 .lt. zero) goto 160
      f2=sqrt(f2)
      if(f0.lt.zero) goto 30                                              !hr06
      if(f0.ge.zero) goto 20                                              !hr06
   20 egwg1=(f1+f2)*half
      egwg2=(f1-f2)*half
      goto 40
   30 egwg1=(f1-f2)*half
      egwg2=(f1+f2)*half
   40 continue
      f1=egwg1**2-four                                                   !hr06
      f2=egwg2**2-four                                                   !hr06
      rca1=f1
      yca1=zero
      rca2=f2
      yca2=zero
      if (rca1.ge.0) then
        rca1=sqrt(rca1)
      else
        yca1=sqrt(-one*rca1)                                             !hr06
        rca1=zero
      endif
      if (rca2.ge.0) then
        rca2=sqrt(rca2)
      else
        yca2=sqrt(-one*rca2)                                             !hr06
        rca2=zero
      endif
      rclam1=(egwg1+rca1)*half
      yclam1=yca1*half
      rclam2=(egwg2+rca2)*half
      yclam2=yca2*half
      if(egwg1**2 .ge. four) goto 160                                    !hr06
      if(egwg2**2 .ge. four) goto 160                                    !hr06
   50 continue
      detb=am(1,3)*am(2,4)-am(1,4)*am(2,3)
      detc=am(3,1)*am(4,2)-am(3,2)*am(4,1)
      fak1=spd-egwg1
      if(abs(fak1).gt.pieni) then
        rcw1(1)=am(1,2)-(am(1,3)*am(3,2)+am(1,4)*am(4,2))/fak1           !hr06
        ycw1(1)=zero
        rcw1(2)=((am(1,3)*am(3,1)+am(1,4)*am(4,1))+detb)/fak1-(am(1,1)  &!hr06
     &-rclam1)                                                           !hr06
        ycw1(2)=yclam1
      rcw1(3)=-one*((am(3,1)+am(2,4))*rcw1(1)+(am(3,2)-am(1,4))*rcw1(2))&!hr06
     &/fak1                                                              !hr06
      ycw1(3)=-one*((am(3,1)+am(2,4))*ycw1(1)+(am(3,2)-am(1,4))*ycw1(2))&!hr06
     &/fak1                                                              !hr06
      rcw1(4)=-one*((am(4,1)-am(2,3))*rcw1(1)+(am(4,2)+am(1,3))*rcw1(2))&!hr06
     &/fak1                                                              !hr06
      ycw1(4)=-one*((am(4,1)-am(2,3))*ycw1(1)+(am(4,2)+am(1,3))*ycw1(2))&!hr06
     &/fak1                                                              !hr06
      else
        rcw1(1)=am(1,2)
        ycw1(1)=zero
        rcw1(2)=rclam1-am(1,1)                                           !hr06
        ycw1(2)=yclam1
        rcw1(3)=zero
        ycw1(3)=zero
        rcw1(4)=zero
        ycw1(4)=zero
      endif
      fak2=spa-egwg2
      if(abs(fak2).gt.pieni) then
        rcw2(3)=am(3,4)-(am(3,1)*am(1,4)+am(3,2)*am(2,4))/fak2           !hr06
        ycw2(3)=zero
        rcw2(4)=((am(3,1)*am(1,3)+am(3,2)*am(2,3))+detc)/fak2-(am(3,3)  &!hr06
     &-rclam2)                                                           !hr06
        ycw2(4)=yclam2
      rcw2(1)=-one*((am(1,3)+am(4,2))*rcw2(3)+(am(1,4)-am(3,2))*rcw2(4))&!hr06
     &/fak2                                                              !hr06
      ycw2(1)=-one*((am(1,3)+am(4,2))*ycw2(3)+(am(1,4)-am(3,2))*ycw2(4))&!hr06
     &/fak2                                                              !hr06
      rcw2(2)=-one*((am(2,3)-am(4,1))*rcw2(3)+(am(2,4)+am(3,1))*rcw2(4))&!hr06
     &/fak2                                                              !hr06
      ycw2(2)=-one*((am(2,3)-am(4,1))*ycw2(3)+(am(2,4)+am(3,1))*ycw2(4))&!hr06
     &/fak2                                                              !hr06
      else
        rcw2(3)=am(3,4)
        ycw2(3)=zero
        rcw2(4)=rclam2-am(3,3)                                           !hr06
        ycw2(4)=yclam2
        rcw2(1)=zero
        ycw2(1)=zero
        rcw2(2)=zero
        ycw2(2)=zero
      endif

!--LEAVING COMPLEX NUMBERS
      do 60 i=1,4
        ta(i,1)=rcw1(i)
        ta(i,3)=rcw2(i)
        ta(i,2)=ycw1(i)
        ta(i,4)=ycw2(i)
   60 continue

!--NORMALISATION OF EIGENVALUES
      rn1=((ta(1,1)*ta(2,2)-ta(2,1)*ta(1,2))                            &!hr06
     &+ta(3,1)*ta(4,2))-ta(4,1)*ta(3,2)                                  !hr06
      if(rn1.lt.zero) goto 70                                             !hr06
      if(rn1.eq.zero) goto 160                                            !hr06
      if(rn1.gt.zero) goto 90                                             !hr06
   70 yclam1=-one*yclam1                                                 !hr06

      do i=1,4
        ta(i,2)=-one*ta(i,2)                                               !hr06
      end do

   90 sqrn=sqrt(abs(rn1))

      do i=1,4
        ta(i,1)=ta(i,1)/sqrn
        ta(i,2)=ta(i,2)/sqrn
      end do

      rn2=((ta(1,3)*ta(2,4)-ta(2,3)*ta(1,4))                            &!hr06
     &+ta(3,3)*ta(4,4))-ta(4,3)*ta(3,4)                                  !hr06
      if(rn2.lt.zero) goto 110                                           !hr06
      if(rn2.eq.zero) goto 160                                           !hr06
      if(rn2.gt.zero) goto 130                                           !hr06
  110 yclam2=-one*yclam2                                                 !hr06

      do i=1,4
        ta(i,4)=-one*ta(i,4)                                               !hr06
      end do

  130 sqrn=sqrt(abs(rn2))

      do i=1,4
        ta(i,3)=ta(i,3)/sqrn
        ta(i,4)=ta(i,4)/sqrn
      end do

      qw(1)= atan_mb(yclam1/(one+rclam1))/pi
      qw(2)= atan_mb(yclam2/(one+rclam2))/pi

!-----------------------------------------------------------------------
!  OPTICAL PARAMETERS AT THE STARTING POINT
!-----------------------------------------------------------------------
      betx(1)=ta(1,1)**2+ta(1,2)**2                                      !hr06
      alfx(1)=-one*(ta(1,1)*ta(2,1)+ta(1,2)*ta(2,2))                     !hr06
      betx(2)=ta(1,3)**2+ta(1,4)**2                                      !hr06
      alfx(2)=-one*(ta(1,3)*ta(2,3)+ta(1,4)*ta(2,4))                     !hr06
      betz(1)=ta(3,1)**2+ta(3,2)**2                                      !hr06
      alfz(1)=-one*(ta(3,1)*ta(4,1)+ta(3,2)*ta(4,2))                     !hr06
      betz(2)=ta(3,3)**2+ta(3,4)**2                                      !hr06
      alfz(2)=-one*(ta(3,3)*ta(4,3)+ta(3,4)*ta(4,4))                     !hr06
      bet0(1)=betx(1)
      alf0(1)=alfx(1)
      bet0(2)=betz(2)
      alf0(2)=alfz(2)

      if(ta(1,1).lt.-pieni) then
        do i=1,4
          do j=1,4
            ta(i,j)=-one*ta(i,j)
          end do
        end do
      endif

      return
!-----------------------------------------------------------------------
  160 ierro=1
      return
end

subroutine block
!-----------------------------------------------------------------------
!  COMBINATION OF LINEAR ELEMENTS TO ONE MATRIX
!-----------------------------------------------------------------------
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      use parpro
      use mod_common
      use mod_commons
      use mod_commont
      implicit none
      integer i,j,jm,k,l,m,n
      real(kind=fPrec) g,h
      dimension h(nblo,2,6),g(nblo,2,6)
      save
!-----------------------------------------------------------------------
      do k=1,mblo
        jm=mel(k)
        i=mtyp(k,1)
        n=mtyp(k,jm)

        do l=1,2
          do m=1,6
            h(1,l,m)=a(i,l,m)
            g(1,l,m)=a(n,l,m)
          end do
        end do

        if(jm.eq.1) goto 40

        do j=2,jm
          i=mtyp(k,j)
          n=mtyp(k,jm-j+1)
          do l=1,2
            h(j,l,1)=h(j-1,l,1)*a(i,l,1)+h(j-1,l,3)*a(i,l,2)
            h(j,l,2)=h(j-1,l,2)*a(i,l,1)+h(j-1,l,4)*a(i,l,2)
            h(j,l,3)=h(j-1,l,1)*a(i,l,3)+h(j-1,l,3)*a(i,l,4)
            h(j,l,4)=h(j-1,l,2)*a(i,l,3)+h(j-1,l,4)*a(i,l,4)
            g(j,l,1)=g(j-1,l,1)*a(n,l,1)+g(j-1,l,3)*a(n,l,2)
            g(j,l,2)=g(j-1,l,2)*a(n,l,1)+g(j-1,l,4)*a(n,l,2)
            g(j,l,3)=g(j-1,l,1)*a(n,l,3)+g(j-1,l,3)*a(n,l,4)
            g(j,l,4)=g(j-1,l,2)*a(n,l,3)+g(j-1,l,4)*a(n,l,4)
            h(j,l,5)=(h(j-1,l,5)*a(i,l,1)+h(j-1,l,6)*a(i,l,2))+a(i,l,5)  !hr06
            h(j,l,6)=(h(j-1,l,5)*a(i,l,3)+h(j-1,l,6)*a(i,l,4))+a(i,l,6)  !hr06
            g(j,l,5)=(g(j-1,l,5)*a(n,l,1)+g(j-1,l,6)*a(n,l,2))+a(n,l,5)  !hr06
            g(j,l,6)=(g(j-1,l,5)*a(n,l,3)+g(j-1,l,6)*a(n,l,4))+a(n,l,6)  !hr06
          end do
        end do

   40   do l=1,2
          do m=1,6
            bl1(k,l,m)=h(jm,l,m)
            bl2(k,l,m)=g(jm,l,m)
          end do
        end do

      end do

      return
end

subroutine blockdis(aeg,bl1eg,bl2eg)
!-----------------------------------------------------------------------
!  COMBINATION OF LINEAR ELEMENTS TO ONE MATRIX, USED FOR DISPERSION
!-----------------------------------------------------------------------
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      use parpro
      use mod_common
      use mod_commons
      use mod_commont
      implicit none
      integer i,j,jm,k,l,m,n
      real(kind=fPrec) aeg,bl1eg,bl2eg,g,h
      dimension h(nblo,2,6),g(nblo,2,6)
      dimension aeg(nele,2,6),bl1eg(nblo,2,6),bl2eg(nblo,2,6)
      save
!-----------------------------------------------------------------------
      do k=1,mblo
        jm=mel(k)
        i=mtyp(k,1)
        n=mtyp(k,jm)

        do l=1,2
          do m=1,6
            h(1,l,m)=aeg(i,l,m)
            g(1,l,m)=aeg(n,l,m)
          end do
         end do

        if(jm.eq.1) goto 40

        do j=2,jm
          i=mtyp(k,j)
          n=mtyp(k,jm-j+1)
          do l=1,2
            h(j,l,1)=h(j-1,l,1)*aeg(i,l,1)+h(j-1,l,3)*aeg(i,l,2)
            h(j,l,2)=h(j-1,l,2)*aeg(i,l,1)+h(j-1,l,4)*aeg(i,l,2)
            h(j,l,3)=h(j-1,l,1)*aeg(i,l,3)+h(j-1,l,3)*aeg(i,l,4)
            h(j,l,4)=h(j-1,l,2)*aeg(i,l,3)+h(j-1,l,4)*aeg(i,l,4)
            g(j,l,1)=g(j-1,l,1)*aeg(n,l,1)+g(j-1,l,3)*aeg(n,l,2)
            g(j,l,2)=g(j-1,l,2)*aeg(n,l,1)+g(j-1,l,4)*aeg(n,l,2)
            g(j,l,3)=g(j-1,l,1)*aeg(n,l,3)+g(j-1,l,3)*aeg(n,l,4)
            g(j,l,4)=g(j-1,l,2)*aeg(n,l,3)+g(j-1,l,4)*aeg(n,l,4)
            h(j,l,5)=(h(j-1,l,5)*aeg(i,l,1)+h(j-1,l,6)*aeg(i,l,2))+aeg  &!hr06
     &(i,l,5)                                                            !hr06
            h(j,l,6)=(h(j-1,l,5)*aeg(i,l,3)+h(j-1,l,6)*aeg(i,l,4))+aeg  &!hr06
     &(i,l,6)                                                            !hr06
            g(j,l,5)=(g(j-1,l,5)*aeg(n,l,1)+g(j-1,l,6)*aeg(n,l,2))+aeg  &!hr06
     &(n,l,5)                                                            !hr06
            g(j,l,6)=(g(j-1,l,5)*aeg(n,l,3)+g(j-1,l,6)*aeg(n,l,4))+aeg  &!hr06
     &(n,l,6)                                                            !hr06
          end do
        end do

   40   do l=1,2
          do m=1,6
            bl1eg(k,l,m)=h(jm,l,m)
            bl2eg(k,l,m)=g(jm,l,m)
          end do
        end do

      end do

      return
end

subroutine chroma
!-----------------------------------------------------------------------
!  CALCULATION OF CHROMATICITY FROM 5 ENERGIE-VALUES
!-----------------------------------------------------------------------
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      use crcoall
      use parpro
      use mod_common
      use mod_commons
      use mod_commont
      implicit none
      integer i,ii,isl,j,jj,l,n
      real(kind=fPrec) cor,coro,cro0,de2,det,dm,dpp,dsm,ox,oz,qwc,sens,sm0,su2,suxy,suzy,xi,zi
      dimension dsm(2,4),sens(2,4),xi(2),zi(2),dm(2),sm0(2)
      dimension qwc(3),cro0(2)
      save
!-----------------------------------------------------------------------
      cor=zero
      coro=1.0e38_fPrec

      do i=1,2
        do j=1,4
          dsm(i,j)=zero
          sens(i,j)=zero
        end do
      end do

      do i=1,2
        xi(i)=zero
        zi(i)=zero
        dm(i)=zero
        sm0(i)=zero
        qwc(i)=zero
        cro0(i)=zero
      end do

      qwc(3)=zero
      write(lout,10010)
      dsm(1,2)=dsm0
      dsm(2,3)=dsm0
      de2=de0*half
      do 90 jj=1,itcro
        do 80 ii=1,4
          su2=zero
          suxy=zero
          suzy=zero
          do 30 l=1,2
            isl=is(l)
            if(kz(isl).ne.3) then
              write(lout,"(a)") "CHROMA> ERROR Element specified for chromaticity correction is not a sextupole."
              call prror(-1)
            end if
            ed(isl)=ed(isl)+dsm(l,ii)
            if(kp(isl).eq.5) call combel(isl)
   30     continue
          do 40 n=1,5
            dpp=de2*real(3-n,fPrec)                                            !hr06
            call clorb(dpp)
            if(ierro.gt.0) then
              write(lout,"(a)") "CHROMA> ERROR Unstable closed orbit during chromaticity correction."
              call prror(-1)
            end if
            call phasad(dpp,qwc)
            if(ierro.gt.0) then
              write(lout,"(a)") "CHROMA> ERROR No optical solution during chromaticity correction."
              call prror(-1)
            end if
            ox=qwc(1)
            oz=qwc(2)
            su2=su2+dpp**2                                               !hr06
            suxy=suxy+ox*dpp
            suzy=suzy+oz*dpp
   40     continue
          do 50 l=1,2
            isl=is(l)
            ed(isl)=ed(isl)-dsm(l,ii)
            if(kp(isl).eq.5) call combel(isl)
   50     continue
          sens(1,ii)=suxy/su2
          sens(2,ii)=suzy/su2
          if(ii.ne.3) goto 80

!--COMPENSATION OF CHROMATICITY
          do l=1,2
            cro0(l)=sens(l,1)-cro(l)
            xi(l)=(sens(1,l+1)-sens(1,1))/dsm0
            zi(l)=(sens(2,l+1)-sens(2,1))/dsm0
          end do

          cor=sqrt(cro0(1)**2+cro0(2)**2)                                !hr06
          if(jj.eq.1.or.cor.lt.coro) then
            coro=cor
            det=xi(1)*zi(2)-zi(1)*xi(2)
            dm(1)=(cro0(2)*xi(2)-cro0(1)*zi(2))/det                      !hr06
            dm(2)=(cro0(1)*zi(1)-cro0(2)*xi(1))/det                      !hr06

            do l=1,2
              sm0(l)=ed(is(l))
              isl=is(l)
              ed(isl)=ed(isl)+dm(l)
              if(kp(isl).eq.5) call combel(isl)
            end do
          else
            write(lout,10035)
            return
          endif
   80   continue
        write(lout,10020) sens(1,1),sens(1,4),sens(2,1),sens(2,4)
        chromc(1)=sens(1,4)*c1m3
        chromc(2)=sens(2,4)*c1m3
        write(lout,10030) sm0(1),ed(is(1)),bez(is(1)), sm0(2),ed(is(2)),&
     &bez(is(2))
        write(lout,10040) xi,zi
        write(lout,10010)
        if(abs(sens(1,4)-cro(1)).lt.dech.and.abs(sens(2,4)-cro(2))      &
     &.lt.dech) return
   90 continue
      write(lout,10000) itcro
!-----------------------------------------------------------------------
      return
10000 format(/131('-')//t10,'CHROMATICITY CORRECTION'/t10,              &
     &'MAXIMUM NUMBER OF ITERATIONS ACHIEVED--->',2x,i4/ t10,           &
     &'PROCEDURE MAY NOT HAVE CONVERGED')
10010 format(/131('-'))
10020 format(/131('-')//t10,'DATA BLOCK CHROMATICITY CORRECTION'/t10,   &
     &'CHROMATICITIES         BEFORE           AFTER CORRECTION'/t10,   &
     &'HORIZONTAL       ',d17.10,7x,d17.10/ t10,'VERTICAL         ',d17.&
     &10,7x,d17.10/)
10040 format(t10,'SEXTUPOLE SENSITIVITIES    XI/M1 XI/M2 YI/M1 YI/M2  ',&
     &4d15.8)
10030 format(t10,'SEXTUP.STRENGTHS ',g17.10,7x,g17.10,'   INDEX   ',a16/&
     &t10,'                 ',g17.10,7x,g17.10,'           ',a16)
10035 format(/t5,'---- NO Improvement in last Step ----'/)
end subroutine chroma

      subroutine chromda
!-----------------------------------------------------------------------
!  CHROMATICITY CORRECTION VIA DA
!-----------------------------------------------------------------------
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      use crcoall
      use parpro
      use mod_common
      use mod_commons
      use mod_commont
      use mod_commond
      implicit none
      integer icht,iq1,iq2,ix,ncorr,ncorruo,nd,nd2
      real(kind=fPrec) cor,coro,dps0,dq1,dq2,edcor1,edcor2,qw,qwc
      dimension qw(2),qwc(3)
      save
!-----------------------------------------------------------------------
      write(lout,10000)
      nd=2
      nd2=4
      dps(1)=dp1+dppoff
      ncorruo=ncorru
      ncorru=1
      call clorb(dp1)
      call betalf(dp1,qw)
      call phasad(dp1,qwc)
      if(nbeam.ge.1) then
#include "include/beamcou.f90"
#ifdef DEBUG
!     call dumpbin('abeamcou2',3,33)
!     call abend('after beam coupling                               ')
#endif
      endif
      ncorru=ncorruo
      iq1=is(1)
      iq2=is(2)
      edcor(1)=ed(iq1)
      edcor(2)=ed(iq2)
      edcor1=edcor(1)
      edcor2=edcor(2)
      coro=1e38_fPrec
      cor=0
      ncorr=0
      do ncorr=1,itcro+1
        ichromc=2
        call mydaini(1,1,nd2,nd,nd2,1)
        ichromc=1
        call mydaini(2,4,7,2,5,1)
        dq1=corr(1,1)-cro(1)*c1m3
        dq2=corr(1,2)-cro(2)*c1m3
        if(ncorr.eq.1) cor=c1e3*sqrt(dq1**2+dq2**2)                      !hr06
        if(cor.gt.dech) then
          cor=c1e3*sqrt(dq1**2+dq2**2)                                   !hr06
          if(ncorr.eq.1.or.cor.lt.coro) then
            coro=cor
            ed(iq1)=(ed(iq1)-corr(2,1)*dq1)-corr(2,2)*dq2                !hr06
            ed(iq2)=(ed(iq2)-corr(3,1)*dq1)-corr(3,2)*dq2                !hr06
            do icht=1,iu
              ix=ic(icht)
              if(ix.gt.nblo) then
                ix=ix-nblo
                if(ix.eq.iq1.or.iratioe(ix).eq.iq1) then
                  smi(icht)=ed(iq1)*ratioe(ix)+smizf(icht)
                else if(ix.eq.iq2.or.iratioe(ix).eq.iq2) then
                  smi(icht)=ed(iq2)*ratioe(ix)+smizf(icht)
                endif
              endif
            enddo
            edcor(1)=ed(iq1)
            edcor(2)=ed(iq2)
            if(ncorr.eq.1) then
              write(lout,10010) cro(1),corr(1,1)*c1e3,cro(2),           &
     &corr(1,2)*c1e3,ncorr-1,cor
              write(lout,10030) edcor1,ed(iq1),bez(iq1),edcor2,ed(iq2), &
     &bez(iq2)
            else
              write(lout,10020) cro(1),corr(1,1)*c1e3,cro(2),           &
     &corr(1,2)*c1e3,ncorr-1,cor
              write(lout,10030) edcor1,ed(iq1),bez(iq1),edcor2,ed(iq2), &
     &bez(iq2)
            endif
          else
            write(lout,10040) ncorr-1
            goto 1
          endif
        else
          write(lout,10050) ncorr-1
          goto 1
        endif
      enddo
 1    continue
      chromc(1)=corr(1,1)
      chromc(2)=corr(1,2)
      if(ncorr.eq.itcro+1) write(lout,10060) itcro
      if(ncorr.eq.1) then
        write(lout,10010) cro(1),corr(1,1)*c1e3,cro(2),                 &
     &corr(1,2)*c1e3,ncorr-1,cor
      else
        write(lout,10020) cro(1),corr(1,1)*c1e3,cro(2),corr(1,2)*c1e3,  &
     &ncorr-1,cor
      endif
      write(lout,10030) edcor1,ed(iq1),bez(iq1),edcor2,ed(iq2),bez(iq2)
!-----------------------------------------------------------------------
10000 format(/131('-')/t10,'ENTERING DA CHROMATICITY CORRECTION'/)
10010 format(/131('-')/t10,                                             &
     &'CHROMATICITY'   ,18x,'THEORET.        BEFORE CORRECTION'/ t10,   &
     &'HORIZONTAL'     ,15x,G21.14,1x,G21.14/ t10,                      &
     &'VERTICAL'       ,17x,G21.14,1x,G21.14// t10,                     &
     &'ITERATION:'     ,21x,i3/ t10,                                    &
     &'ACCURACY:'      ,17x,g17.10/)
10020 format(/131('-')/t10,                                             &
     &'CHROMATICITY'   ,18x,'THEORET.        AFTER CORRECTION'/ t10,    &
     &'HORIZONTAL'     ,15x,G21.14,1x,G21.14/ t10,                      &
     &'VERTICAL'       ,17x,G21.14,1x,G21.14// t10,                     &
     &'ITERATION:'     ,21x,i3/ t10,                                    &
     &'ACCURACY:'      ,17x,g17.10/)
10030 format(t10,'SEXTUPOLE STRENGTH',5x,g17.10,2x,g17.10,'   TYP     ',&
     &a16/t10,                  23x,g17.10,2x,g17.10,'           ',     &
     &a16)
10040 format(/t5,'---- NO IMPROVEMENT OF DA CHROMATICITY CORRECTION ',  &
     &'IN ITERATION: ',i4/)
10050 format(t5/t10,'DA CHROMATICITY CORRECTION SUCCESSFUL IN ',        &
     &'ITERATION: ',i4/)
10060 format(/t10,'DA CHROMATICITY CORRECTION'/ t10,                    &
     &'MAXIMUM NUMBER OF ITERATIONS ACHIEVED--->',2x,i4/ t10,           &
     &'PROCEDURE MAY NOT HAVE CONVERGED')
end

subroutine clorb(dpp)
!-----------------------------------------------------------------------
!  CALCULATION OF THE CLOSED ORBIT   'CLO(2),CLOP(2)'
!-----------------------------------------------------------------------
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      use crcoall
      use parpro
      use mod_settings
      use mod_common
      use mod_commons
      use mod_commont

      implicit none
      integer ierr,ii,l,ll
      real(kind=fPrec) am,cor,dclo,dclop,dcx,dcxp,dcz,dczp,det,dpp,dx,  &
     &dy,x0,x1,y0,y1
      dimension x1(2),y1(2),x0(2),y0(2)
      dimension dclo(2),dclop(2)
      dimension dx(2),dy(2),am(4,4)
      save ! Saving DPP?
!-----------------------------------------------------------------------
      ierro=0
      do 10 l=1,2
        clo(l)=dpp*di0(l)
        clop(l)=dpp*dip0(l)
        dx(l)=c1e6
        dy(l)=c1e6
   10 continue
      call envar(dpp)
      call umlauf(dpp,1,ierr)
      ierro=ierr
      if(ierro.ne.0) return
      do 40 ii=1,itco
        dcx=abs(dx(1))
        dcxp=abs(dy(1))
        dcz=abs(dx(2))
        dczp=abs(dy(2))
        if(dcx.le.dma.and.dcz.le.dma.and.dcxp.le.dmap.and.dczp.le.dmap) &
     &goto 50

        do l=1,2
          x(1,l)=clo(l)
          y(1,l)=clop(l)
          x0(l)=x(1,l)
         y0(l)=y(1,l)
        end do

        call matrix(dpp,am)
        if(ierro.ne.0) return
        do 30 l=1,2
          ll=2*l
          x1(l)=x(1,l)
          y1(l)=y(1,l)
          det=(two-am(ll-1,ll-1))-am(ll,ll)                              !hr06
          dx(l)=x0(l)-x1(l)
          dy(l)=y0(l)-y1(l)
          dclo(l)=(dx(l)*(am(ll,ll)-one)-dy(l)*am(ll-1,ll))/det
          dclop(l)=(dy(l)*(am(ll-1,ll-1)-one)-dx(l)*am(ll,ll-1))/det
          clo(l)=clo(l)+dclo(l)
          clop(l)=clop(l)+dclop(l)
   30   continue
   40 continue
      if(ncorru.ne.1) write(lout,10000) itco
   50 cor=c1e3*sqrt(dcx**2+dcz**2)                                       !hr06
      if(st_print .and. ncorru /= 1) then
        write(lout,10010) dpp,clo(1),clop(1),clo(2),clop(2),ii,cor
#ifdef DEBUG
!     call warr('dpp',dpp,0,0,0,0)
!     call warr('dpp',dpp,0,0,0,0)
!     call warr('clo(1)',clo(1),0,0,0,0)
!     call warr('clop(1)',clop(1),0,0,0,0)
!     call warr('clo(2)',clo(2),0,0,0,0)
!     call warr('clop(2)',clop(2),0,0,0,0)
!     call warr('ii',0d0,ii,0,0,0)
!     call warr('cor',cor,0,0,0,0)
#endif
      endif
!-----------------------------------------------------------------------
      return
10000 format(t5/t10,'CLOSED ORBIT CALCULATION'/ t10,                    &
     &'MAXIMUM NUMBER OF ITERATIONS ACHIEVED--->',2x,i4/ t10,           &
     &'PROCEDURE MAY NOT HAVE CONVERGED')
10010 format(t5,'---- ENTRY CLORB ----/DPP=',f8.5,' /CLOX/', 2f10.5,    &
     &' /CLOY/',2f10.5,' /ITERAT.=',i3,'/ ACCURACY=',d13.6)
end

subroutine clorb2(dpp)
!-----------------------------------------------------------------------
!  CALCULATION OF THE CLOSED ORBIT - NO WRITEOUT
!-----------------------------------------------------------------------
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      use parpro
      use crcoall
      use mod_common
      use mod_commons
      use mod_commont
      implicit none
      integer ierr,ii,l,ll
      real(kind=fPrec) am,dclo,dclop,dcx,dcxp,dcz,dczp,det,dpp,dx,dy,x0,x1,y0,y1
      dimension x1(2),y1(2),x0(2),y0(2)
      dimension dclo(2),dclop(2)
      dimension dx(2),dy(2),am(4,4)
      save
!-----------------------------------------------------------------------
      ierro=0
      do 10 l=1,2
        clo(l)=dpp*di0(l)
        clop(l)=dpp*dip0(l)
        dx(l)=c1e6                                                        !hr06
        dy(l)=c1e6                                                        !hr06
   10 continue

      call envar(dpp)
      call umlauf(dpp,1,ierr)
      ierro=ierr
      if(ierro /= 0) then
        write(lout,"(a)") "CLORB> ERROR No convergence in rmod."
        call prror(-1)
      end if

      do 40 ii=1,itco
        dcx=abs(dx(1))
        dcxp=abs(dy(1))
        dcz=abs(dx(2))
        dczp=abs(dy(2))
        if(dcx.le.dma.and.dcz.le.dma.and.dcxp.le.dmap.and.dczp.le.dmap) &
     &return

        do l=1,2
          x(1,l)=clo(l)
          y(1,l)=clop(l)
          x0(l)=x(1,l)
          y0(l)=y(1,l)
        end do

        call matrix(dpp,am)

        if(ierro /= 0) then
          write(lout,"(a)") "CLORB> ERROR No convergence in rmod."
          call prror(-1)
        end if

        do 30 l=1,2
          ll=2*l
          x1(l)=x(1,l)
          y1(l)=y(1,l)
          det=two-am(ll-1,ll-1)-am(ll,ll)
          dx(l)=x0(l)-x1(l)
          dy(l)=y0(l)-y1(l)
          dclo(l)=(dx(l)*(am(ll,ll)-one)-dy(l)*am(ll-1,ll))/det
          dclop(l)=(dy(l)*(am(ll-1,ll-1)-one)-dx(l)*am(ll,ll-1))/det
          clo(l)=clo(l)+dclo(l)
          clop(l)=clop(l)+dclop(l)
   30   continue

   40 continue
!-----------------------------------------------------------------------
      return
end subroutine clorb2

subroutine combel(iql)
!-----------------------------------------------------------------------
!  COMBINATION OF ELEMENTS
!-----------------------------------------------------------------------
      use floatPrecision
      use crcoall
      use numerical_constants
      use mathlib_bouncer
      use parpro
      use mod_common
      use mod_commons
      use mod_commont
      implicit none
      integer ico,ico0,iql,j,m
      save
!-----------------------------------------------------------------------
      do 20 j=1,icoe
        ico0=icomb0(j)
        if(iql.ne.ico0) goto 20
        do 10 m=1,20
          ico=icomb(j,m)
          if(ico.eq.0) goto 10
          if(kz(ico0).ne.kz(ico)) then
            write(lout,"(a)") "COMBEL> ERROR Elements of different types are combined in data block combination of elements."
            call prror(-1)
          end if
          if(abs(el(ico0)).gt.pieni) then
            if(abs(el(ico)).gt.pieni) then
              ek(ico)=ek(ico0)*ratio(j,m)
            else
              ed(ico)=ek(ico0)*ratio(j,m)
            endif
          endif
          if(abs(el(ico0)).le.pieni) then
            if(abs(el(ico)).le.pieni) then
              ed(ico)=ed(ico0)*ratio(j,m)
            else
              ek(ico)=ed(ico0)*ratio(j,m)
            endif
          endif
   10   continue
   20 continue
!-----------------------------------------------------------------------
      return
end subroutine combel

!-----------------------------------------------------------------------
!  CALCULATION OF ELEMENT MATRICES
!  Rewritten from computed goto to select case 16/11/2017, VKBO
!-----------------------------------------------------------------------
subroutine envar(dpp)

    use floatPrecision
    use numerical_constants
    use mathlib_bouncer

    use parpro
    use mod_common
    use mod_commons
    use mod_commont
    implicit none
    integer i,ih,kz1,l,ll
    real(kind=fPrec) afok,co,dpd,dpp,dpsq,fi,fok,fokq,g,gl,hc,hi,hi1,hm,hp,hs,rho,rhoi,si,wf

    save

    dpd  = one+dpp
    dpsq = sqrt(dpd)
    do i=1,il
        do ll=1,6
            do l=1,2
                a(i,l,ll)=zero
            end do
        end do
        if (abs(el(i)).le.pieni) then ! NONLINEAR INSERTION
            sm(i)=ed(i)
            cycle
        end if

        kz1 = kz(i)+1
        select case(kz1)
        case (1) ! DRIFTLENGTH
            do l=1,2
                a(i,l,1) = one
                a(i,l,2) = el(i)
                a(i,l,3) = zero
                a(i,l,4) = one
            end do

        case (2, 4, 5, 6)
            ! 2: RECTANGULAR MAGNET
            ! 4: SEKTORMAGNET
            ! 5: RECTANGULAR MAGNET VERTICAL
            ! 6: SEKTORMAGNET VERTICAL
            fok = el(i)*ed(i)/dpsq
            if(abs(fok).le.pieni) then
                do l=1,2
                    a(i,l,1) = one
                    a(i,l,2) = el(i)
                    a(i,l,3) = zero
                    a(i,l,4) = one
                end do
                cycle
            end if
            rho = (one/ed(i))*dpsq
            si  = sin_mb(fok)
            co  = cos_mb(fok)
            g   = tan_mb(fok*half)/rho
            gl  = el(i)*g
            select case (kz1)
            case (2)
                ! HORIZONTAL
                a(i,1,1) = one
                a(i,1,2) = rho*si
                a(i,1,3) = zero
                a(i,1,4) = one
                a(i,1,5) = ((-one*rho)*(one-co))/dpsq         ! hr06
                a(i,1,6) = ((-one*two)*tan_mb(fok*half))/dpsq ! hr06
                ! VERTICAL
                a(i,2,1) = one-gl
                a(i,2,2) = el(i)
                a(i,2,3) = (-one*g)*(two-gl) ! hr06
                a(i,2,4) = a(i,2,1)
            case (4)
                ! HORIZONTAL
                a(i,1,1) = co
                a(i,1,2) = rho*si
                a(i,1,3) = (-one*si)/rho              ! hr06
                a(i,1,4) = co
                a(i,1,5) = ((-one*rho)*(one-co))/dpsq ! hr06
                a(i,1,6) = (-one*si)/dpsq             ! hr06
                ! VERTICAL
                a(i,2,1) = one
                a(i,2,2) = el(i)
                a(i,2,3) = zero
                a(i,2,4) = one
            case (5)
                ! HORIZONTAL
                a(i,2,1) = one
                a(i,2,2) = rho*si
                a(i,2,3) = zero
                a(i,2,4) = one
                a(i,2,5) = ((-one*rho)*(one-co))/dpsq         ! hr06
                a(i,2,6) = ((-one*two)*tan_mb(fok*half))/dpsq ! hr06
                ! VERTIKAL
                a(i,1,1) = one-gl
                a(i,1,2) = el(i)
                a(i,1,3) = (-one*g)*(two-gl) ! hr06
                a(i,1,4) = a(i,1,1)
            case (6)
                ! HORIZONTAL
                a(i,2,1) = co
                a(i,2,2) = rho*si
                a(i,2,3) = (-one*si)/rho              ! hr06
                a(i,2,4) = co
                a(i,2,5) = ((-one*rho)*(one-co))/dpsq ! hr06
                a(i,2,6) = (-one*si)/dpsq             ! hr06
                ! VERTIKAL
                a(i,1,1) = one
                a(i,1,2) = el(i)
                a(i,1,3) = zero
                a(i,1,4) = one
            end select

        case (3) ! QUADRUPOLE
            fok=ek(i)/(one+dpp)
            if(abs(fok).le.pieni) then
                do l=1,2
                    a(i,l,1) = one
                    a(i,l,2) = el(i)
                    a(i,l,3) = zero
                    a(i,l,4) = one
                end do
                cycle
            end if
            ih = 0
            hi = sqrt(abs(fok))
            fi = el(i)*hi
            if(fok.gt.zero) goto 110
    100     ih = ih+1
            a(i,ih,1) = cos_mb(fi)
            hi1 = sin_mb(fi)
            a(i,ih,2) = hi1/hi
            a(i,ih,3) = (-one*hi1)*hi  ! hr06
            a(i,ih,4) = a(i,ih,1)
            if(ih.eq.2) cycle
            !--DEFOCUSSING
    110     ih = ih+1
            hp = exp_mb(fi)
            hm = one/hp
            hc = (hp+hm)*half
            hs = (hp-hm)*half
            a(i,ih,1) = hc
            a(i,ih,2) = hs/hi
            a(i,ih,3) = hs*hi
            a(i,ih,4) = hc
            if(ih.eq.1) goto 100

        case (7, 8)
            ! 7: COMBINED FUNCTION MAGNET HORIZONTAL
            ! 8: COMBINED FUNCTION MAGNET VERTICAL
            if (kz1.eq.7) then
                ih   = 0
                fokq = ek(i)
            else
                ih   = 1
                fokq = -one*ek(i) ! hr06
            end if
            wf  = ed(i)/dpsq
            fok = fokq/dpd-wf**2  ! hr06
            if(abs(fok).le.pieni) then
                do l=1,2
                    a(i,l,1) = one
                    a(i,l,2) = el(i)
                    a(i,l,3) = zero
                    a(i,l,4) = one
                end do
                cycle
            end if
            afok = abs(fok)
            hi   = sqrt(afok)
            fi   = hi*el(i)
            if(fok.gt.zero) goto 160
    140     ih = ih+1
            si = sin_mb(fi)
            co = cos_mb(fi)
            a(i,ih,1) = co
            a(i,ih,2) = si/hi
            a(i,ih,3) = (-one*si)*hi                     ! hr06
            a(i,ih,4) = co
            a(i,ih,5) = (((-one*wf)/afok)*(one-co))/dpsq ! hr06
            a(i,ih,6) = (((-one*wf)/hi)*si)/dpsq         ! hr06
            ih = ih+1
            if(ih.gt.2) ih = 1
            hi = sqrt(abs(ek(i)/dpd))
            fi = hi*el(i)
            hp = exp_mb(fi)
            hm = one/hp
            hc = (hp+hm)*half
            hs = (hp-hm)*half
            a(i,ih,1) = hc
            a(i,ih,2) = el(i)
            if(abs(hi).le.pieni) goto 150
            a(i,ih,2) = hs/hi
    150     a(i,ih,3) = hs*hi
            a(i,ih,4) = hc
            cycle
            ! DEFOCUSSING
    160     ih = ih+1
            hp = exp_mb(fi)
            hm = one/hp
            hc = (hp+hm)*half
            hs = (hp-hm)*half
            a(i,ih,1) = hc
            a(i,ih,2) = hs/hi
            a(i,ih,3) = hs*hi
            a(i,ih,4) = hc
            a(i,ih,5) = ((wf/afok)*(one-hc))/dpsq ! hr06
            a(i,ih,6) = (((-one*wf)/hi)*hs)/dpsq  ! hr06
            ih = ih+1
            if(ih.gt.2) ih = 1
            hi = sqrt(abs(ek(i)/dpd))
            fi = hi*el(i)
            si = sin_mb(fi)
            co = cos_mb(fi)
            a(i,ih,1) = co
            a(i,ih,2) = si/hi
            a(i,ih,3) = (-one*si)*hi ! hr06
            a(i,ih,4) = co

        case (9) ! EDGE FOCUSSING
            rhoi = ed(i)/dpsq
            fok  = rhoi*tan_mb((el(i)*rhoi)*half) ! hr06
            a(i,1,1) = one
            a(i,1,2) = zero
            a(i,1,3) = fok
            a(i,1,4) = one
            a(i,2,1) = one
            a(i,2,2) = zero
            a(i,2,3) = -one*fok ! hr06
            a(i,2,4) = one

        end select
    end do
    call block
    return

end subroutine envar

!-----------------------------------------------------------------------
!  CALCULATION OF ELEMENT MATRICES
!  Rewritten from computed goto to select case 16/11/2017, VKBO
!-----------------------------------------------------------------------
subroutine envardis(dpp,aeg,bl1eg,bl2eg)

    use floatPrecision
    use numerical_constants
    use mathlib_bouncer

    use parpro
    use mod_common
    use mod_commons
    use mod_commont
    implicit none
    integer i,ih,kz1,l,ll
    real(kind=fPrec) aeg,afok,bl1eg,bl2eg,co,dpd,dpp,dpsq,fi,fok,fokq,g,gl,hc,hi,hi1,hm,hp,hs,rho,rhoi,si,wf
    dimension aeg(nele,2,6),bl1eg(nblo,2,6),bl2eg(nblo,2,6)
    save

    dpd  = one+dpp
    dpsq = sqrt(dpd)
    do i=1,il
        do ll=1,6
            do l=1,2
              aeg(i,l,ll) = zero
            end do
        end do
        if(abs(el(i)).le.pieni) cycle
        kz1 = kz(i)+1
        select case (kz1)
        case (1)
            do l=1,2
                aeg(i,l,1) = one
                aeg(i,l,2) = el(i)
                aeg(i,l,3) = zero
                aeg(i,l,4) = one
            end do

        case (2, 4, 5, 6)
            ! 2: RECTANGULAR MAGNET
            ! 4: SEKTORMAGNET
            ! 5: RECTANGULAR MAGNET VERTICAL
            ! 6: SEKTORMAGNET VERTICAL
            fok = el(i)*ed(i)/dpsq
            if(abs(fok).le.pieni) then
                do l=1,2
                    aeg(i,l,1) = one
                    aeg(i,l,2) = el(i)
                    aeg(i,l,3) = zero
                    aeg(i,l,4) = one
                end do
            end if
            rho = (one/ed(i))*dpsq
            si  = sin_mb(fok)
            co  = cos_mb(fok)
            select case (kz1)
            case (2)
                ! HORIZONTAL
                aeg(i,1,1) = one
                aeg(i,1,2) = rho*si
                aeg(i,1,3) = zero
                aeg(i,1,4) = one
                aeg(i,1,5) = ((-one*rho)*(one-co))/dpsq         ! hr06
                aeg(i,1,6) = ((-one*two)*tan_mb(fok*half))/dpsq ! hr06
                ! VERTICAL
                g  = tan_mb(fok*half)/rho
                gl = el(i)*g
                aeg(i,2,1) = one-gl
                aeg(i,2,2) = el(i)
                aeg(i,2,3) = (-one*g)*(two-gl)                  ! hr06
                aeg(i,2,4) = aeg(i,2,1)
            case (4)
                ! HORIZONTAL
                aeg(i,1,1)=co
                aeg(i,1,2)=rho*si
                aeg(i,1,3)=(-one*si)/rho                        ! hr06
                aeg(i,1,4)=co
                aeg(i,1,5)=((-one*rho)*(one-co))/dpsq           ! hr06
                aeg(i,1,6)=(-one*si)/dpsq                       ! hr06
                ! VERTICAL
                aeg(i,2,1)=one
                aeg(i,2,2)=el(i)
                aeg(i,2,3)=zero
                aeg(i,2,4)=one
            case (5)
                ! HORIZONTAL
                aeg(i,2,1) = one
                aeg(i,2,2) = rho*si
                aeg(i,2,3) = zero
                aeg(i,2,4) = one
                aeg(i,2,5) = ((-one*rho)*(one-co))/dpsq         ! hr06
                aeg(i,2,6) = ((-one*two)*tan_mb(fok*half))/dpsq ! hr06
                ! VERTICAL
                g  = tan_mb(fok*half)/rho
                gl = el(i)*g
                aeg(i,1,1) = one-gl
                aeg(i,1,2) = el(i)
                aeg(i,1,3) = (-one*g)*(two-gl)                  ! hr06
                aeg(i,1,4) = aeg(i,1,1)
            case (6)
                ! HORIZONTAL
                aeg(i,2,1)=co
                aeg(i,2,2)=rho*si
                aeg(i,2,3)=(-one*si)/rho                        ! hr06
                aeg(i,2,4)=co
                aeg(i,2,5)=((-one*rho)*(one-co))/dpsq           ! hr06
                aeg(i,2,6)=(-one*si)/dpsq                       ! hr06
                ! VERTICAL
                aeg(i,1,1)=one
                aeg(i,1,2)=el(i)
                aeg(i,1,3)=zero
                aeg(i,1,4)=one
            end select

        case (3) ! QUADRUPOLE
            ! FOCUSSING
            fok = ek(i)/(one+dpp)
            if(abs(fok).le.pieni) then
                do l=1,2
                    aeg(i,l,1) = one
                    aeg(i,l,2) = el(i)
                    aeg(i,l,3) = zero
                    aeg(i,l,4) = one
                end do
                cycle
            end if
            ih = 0
            hi = sqrt(abs(fok))
            fi = el(i)*hi
            if(fok.gt.zero) goto 110
    100     ih = ih+1
            aeg(i,ih,1) = cos_mb(fi)
            hi1 = sin_mb(fi)
            aeg(i,ih,2) = hi1/hi
            aeg(i,ih,3) = (-one*hi1)*hi ! hr06
            aeg(i,ih,4) = aeg(i,ih,1)
            if(ih.eq.2) cycle
            ! DEFOCUSSING
    110     ih = ih+1
            hp = exp_mb(fi)
            hm = one/hp
            hc = (hp+hm)*half
            hs = (hp-hm)*half
            aeg(i,ih,1) = hc
            aeg(i,ih,2) = hs/hi
            aeg(i,ih,3) = hs*hi
            aeg(i,ih,4) = hc
            if(ih.eq.1) goto 100

        case (7, 8)
            ! 7: COMBINED FUNCTION MAGNET HORIZONTAL
            ! 8: COMBINED FUNCTION MAGNET VERTICAL
            if (kz1.eq.7) then
                ih   = 0
                fokq = ek(i)
            else
                ih   = 1
                fokq = -one*ek(i)      ! hr06
            end if
            wf  = ed(i)/dpsq
            fok = fokq/dpd-wf**2       ! hr06
            if(abs(fok).le.pieni) then
                do l=1,2
                    aeg(i,l,1) = one
                    aeg(i,l,2) = el(i)
                    aeg(i,l,3) = zero
                    aeg(i,l,4) = one
                end do
                cycle
            end if
            afok = abs(fok)
            hi   = sqrt(afok)
            fi   = hi*el(i)
            if(fok.gt.zero) goto 160
    140     ih = ih+1
            si = sin_mb(fi)
            co = cos_mb(fi)
            aeg(i,ih,1) = co
            aeg(i,ih,2) = si/hi
            aeg(i,ih,3) = (-one*si)*hi                      ! hr06
            aeg(i,ih,4) = co
            aeg(i,ih,5) = (((-one*wf)/afok)*(one-co))/dpsq  ! hr06
            aeg(i,ih,6) = (((-one*wf)/hi)*si)/dpsq          ! hr06
            ih = ih+1
            if(ih.gt.2) ih=1
            hi = sqrt(abs(ek(i)/dpd))
            fi = hi*el(i)
            hp = exp_mb(fi)
            hm = one/hp
            hc = (hp+hm)*half
            hs = (hp-hm)*half
            aeg(i,ih,1) = hc
            aeg(i,ih,2) = el(i)
            if(abs(hi).le.pieni) goto 150
            aeg(i,ih,2) = hs/hi
    150     aeg(i,ih,3) = hs*hi
            aeg(i,ih,4) = hc
            ! DEFOCUSSING
    160     ih = ih+1
            hp = exp_mb(fi)
            hm = one/hp
            hc = (hp+hm)*half
            hs = (hp-hm)*half
            aeg(i,ih,1) = hc
            aeg(i,ih,2) = hs/hi
            aeg(i,ih,3) = hs*hi
            aeg(i,ih,4) = hc
            aeg(i,ih,5) = ((wf/afok)*(one-hc))/dpsq ! hr06
            aeg(i,ih,6) = (((-one*wf)/hi)*hs)/dpsq  ! hr06
            ih = ih+1
            if(ih.gt.2) ih = 1
            hi = sqrt(abs(ek(i)/dpd))
            fi = hi*el(i)
            si = sin_mb(fi)
            co = cos_mb(fi)
            aeg(i,ih,1) = co
            aeg(i,ih,2) = si/hi
            aeg(i,ih,3) = (-one*si)*hi              ! hr06
            aeg(i,ih,4) = co

        case (9) ! EDGE FOCUSSING
            rhoi = ed(i)/dpsq
            fok  = rhoi*tan_mb((el(i)*rhoi)*half) ! hr06
            aeg(i,1,1) = one
            aeg(i,1,2) = zero
            aeg(i,1,3) = fok
            aeg(i,1,4) = one
            aeg(i,2,1) = one
            aeg(i,2,2) = zero
            aeg(i,2,3) = -one*fok
            aeg(i,2,4) = one

        end select
    end do
    call blockdis(aeg,bl1eg,bl2eg)
    return

end subroutine envardis

!-----------------------------------------------------------------------
!  LINEAR PARAMETERS AT THE POSITION OF EVERY ELEMENT OR BLOCK
!-----------------------------------------------------------------------
subroutine linopt(dpp)

  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  use crcoall
  use parpro
  use mod_common
  use mod_commons
  use mod_commont

#ifdef ROOT
  use root_output
#endif

#ifdef HDF5
  use hdf5_output
  use hdf5_linopt
#endif

  use collimation

  implicit none

  integer i,iiii,im,ium,ix,izu,j,jj,jk,jm,k,kpz,kzz,l,l1,ll,nmz,nr,dj
  real(kind=fPrec) aa,aeg,alfa,bb,benkr,beta,bexi,bezii,bl1eg,bl2eg,ci,cikve,clo0,clop0,cr,crkve, &
    crkveuk,di00,dip00,dphi,dpp,dpp1,dppi,dpr,dyy1,dyy2,ekk,etl,phi,phibf,pie,puf,qu,qv,qw,qwc,r0,&
    r0a,t,xl,xs,zl,zs,quz,qvz
#ifdef TILT
  real(kind=fPrec) dyy11,qu1,tiltck,tiltsk
#endif
  character(len=mNameLen) idum

  dimension t(6,4)
  dimension beta(2),alfa(2),phibf(2),phi(2)
  dimension clo0(2),clop0(2),di00(2),dip00(2),qw(2),qwc(3)
  dimension aa(mmul),bb(mmul),dpr(6)
  dimension cr(mmul),ci(mmul)
  dimension aeg(nele,2,6),bl1eg(nblo,2,6),bl2eg(nblo,2,6)
  data dpr/6*zero/
  save
!-----------------------------------------------------------------------

      nhmoni=0
      nvmoni=0
      nhcorr=0
      nvcorr=0
      ium=6
      pie=two*pi

      if(ncorru.eq.0) then
        write(lout,10010)
        write(lout,10000)
      endif

      do i=1,ium
        dpr(i)=zero
      end do

      do i=1,ium
        do j=1,4
          t(i,j)=zero
        end do
      end do

      do i=1,2
        beta(i)=zero
        alfa(i)=zero
        phibf(i)=zero
        phi(i)=zero
        clo0(i)=zero
        clop0(i)=zero
        di00(i)=zero
        dip00(i)=zero
        qw(i)=zero
        qwc(i)=zero
      end do

      qwc(3)=zero

      do i=1,mmul
        aa(i)=zero
        bb(i)=zero
        cr(i)=zero
        ci(i)=zero
      end do

      etl=zero
      dpr(1)=dpp*c1e3
      dpr(6)=one
      dpp1=dpp+ded
      call clorb(dpp1)

      do l=1,2
        clo0(l)=clo(l)
        clop0(l)=clop(l)
      end do

      call clorb(dpp)

      do l=1,2
        ll=2*l
        di0(l)=(clo0(l)-clo(l))/ded
        dip0(l)=(clop0(l)-clop(l))/ded
        t(6,ll-1)=di0(l)
        t(6,ll)=dip0(l)
      end do

      if(ncorru.eq.0) then
        write(lout,10010)
        write(lout,10050) (di0(l),dip0(l),l=1,2)
      endif

      call betalf(dpp,qw)
      call phasad(dpp,qwc)

      if(ierro /= 0) then
        write(lout,"(a)") "LINOPT> ERROR No optical solution."
        call prror(-1)
      end if
      if(ncorru.eq.0) write(lout,10040) dpp,qwc(1),qwc(2)

      call envar(dpp)

      if(ithick.eq.1) call envardis(dpp1,aeg,bl1eg,bl2eg)

!--STARTVALUES OF THE TRAJECTORIES
      do l=1,2
        ll=2*l
        t(1,ll-1)=clo(l)
        t(1,ll)=clop(l)
      end do

      do i=1,4
        do j=1,4
          t(i+1,j)=ta(j,i)
          t(i+1,j)=ta(j,i)
        end do
      end do

      if(ncorru.eq.0 .and. iprint.eq.1) then
        write(lout,10010)
        write(lout,10030)
        write(lout,10020)
        write(lout,10010)
      endif

!--START OF THE MACHINE
      idum='START'
      nr=0
      call writelin(nr,idum,etl,phi,t,1,.false.,0)
      if(ntco.ne.0) then
        if(mod(nr,ntco).eq.0) call cpltwis(idum,t,etl,phi)
      endif

#ifdef ROOT
      if(root_flag .and. root_Optics.eq.1) then
        call OpticsRootWrite()
      end if
#endif

!--STRUCTURE ELEMENT LOOP
      if(nt.le.0.or.nt.gt.iu) nt=iu
      izu=0
#ifdef HDF5
      if(h5_writeOptics) call h5lin_init
#endif

      STRUCTLOOP: do k=1,nt
        ix=ic(k)
        if(ix.gt.nblo) goto 220 !Not a BLOCK
        if(ithick.eq.1.and.iprint.eq.1) goto 160

        jj=0 !initial idx
        dj=1 !step

        if (ix.le.0) then
           ix=-1*ix             !hr13
           jj=mel(ix)+1         !initial idx
           dj=-1                !step
        endif
        jm=mel(ix)
!-- Loop over elements inside the block
        do 150 j=1,jm
          jj=jj+dj       ! Subelement index of current sub=element
          jk=mtyp(ix,jj) ! Single-element index of the current sub-element
          if(ithick.eq.1.and.kz(jk).ne.0) goto 120
          if(ithick.eq.0.and.kz(jk).ne.0) then
            etl=etl+el(jk)

!c$$$            nr=nr+1
!c$$$            call writelin(nr,bez(jk),etl,phi,t,ix,.true.,k)
!c$$$            if(ntco.ne.0) then
!c$$$              if(mod(nr,ntco).eq.0) call cpltwis(bez(jk),t,etl,phi)
!c$$$            endif

            write(lout,"(a)") "LINOPT> ERROR In block '"//trim(bezb(ix))//"': found a thick non-drift element '"//&
              trim(bez(jk))//"' while ithick=1. This should not be possible!"
            call prror(-1)
            cycle STRUCTLOOP
          endif

!--IN BLOCK: PURE DRIFTLENGTH (above: If ITHICK=1 and kz!=0, goto 120->MAGNETELEMENT)
          etl=etl+el(jk)

          do l=1,2
            ll=2*l
            if(abs(t(ll,ll-1)).gt.pieni) then
              phibf(l)=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
            else
              phibf(l)=pi2
            endif
            do i=1,ium
              t(i,ll-1)=t(i,ll-1)+t(i,ll)*(el(jk))
            end do
          end do

          do l=1,2
            ll=2*l
            if(abs(t(ll,ll-1)).gt.pieni) then
              dphi=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
            else
              dphi=pi2-phibf(l)
            endif
            if((-one*dphi).gt.pieni) dphi=dphi+pi                        !hr06
            phi(l)=phi(l)+dphi/pie
          end do

          nr=nr+1
          call writelin(nr,bez(jk),etl,phi,t,ix,.true.,k)
          if(ntco.ne.0) then
            if(mod(nr,ntco).eq.0) call cpltwis(bez(jk),t,etl, phi)
          endif
#ifdef ROOT
      if(root_flag .and. root_Optics.eq.1) then
        call OpticsRootWrite()
      end if
#endif

          goto 150

!--IN BLOCK: MAGNETELEMENT
  120     continue
          if(kz(jk).ne.8) etl=etl+el(jk)
          do l=1,2
            ll=2*l

            if(abs(t(ll,ll-1)).gt.pieni) then
              phibf(l)=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
            else
              phibf(l)=zero
            endif

            puf=t(6,ll-1)
         t(6,ll-1)=(((((aeg(jk,l,1)*(t(1,ll-1)+puf*ded)+ aeg(jk,l,2)*(t &!hr06
     &(1,ll)+t(6,ll)*ded))+aeg(jk,l,5)*dpp1*c1e3)- a(jk,l,1)*t          &!hr06
     &(1,ll-1))-a(jk,l,2)*t(1,ll))- a(jk,l,5)*dpr(1))/ded                !hr06
           t(6,ll)=(((((aeg(jk,l,3)*(t(1,ll-1)+puf*ded)+ aeg(jk,l,4)*(t &!hr06
     &(1,ll)+t(6,ll)*ded))+aeg(jk,l,6)*dpp1*c1e3)- a(jk,l,3)*t          &!hr06
     &(1,ll-1))-a(jk,l,4)*t(1,ll))- a(jk,l,6)*dpr(1))/ded                !hr06

            do i=1,ium-1
              puf=t(i,ll-1)
            t(i,ll-1)=(puf*a(jk,l,1)+t(i,ll)*a(jk,l,2))+dpr(i)*a(jk,l,5) !hr06
            t(i,ll)=(puf*a(jk,l,3)+t(i,ll)*a(jk,l,4))+dpr(i)*a(jk,l,6)   !hr06
            enddo
          enddo

          do l=1,2
            ll=2*l

            if(abs(t(ll,ll-1)).gt.pieni) then
              dphi=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
            else
              dphi=-one*phibf(l)                                         !hr06
            endif

            if(kz(jk).ne.8.and.-one*dphi.gt.pieni) dphi=dphi+pi          !hr06
            phi(l)=phi(l)+dphi/pie
          enddo

          nr=nr+1
          call writelin(nr,bez(jk),etl,phi,t,ix,.true.,k)
          if(ntco.ne.0) then
            if(mod(nr,ntco).eq.0) call cpltwis(bez(jk),t,etl, phi)
          endif
#ifdef ROOT
          if(root_flag .and. root_Optics.eq.1) then
            call OpticsRootWrite()
          end if
#endif

  150   continue !End of loop over elements inside block

        nr=nr+1
        call writelin(nr,bezb(ix),etl,phi,t,ix,.true.,k)
        if(ntco.ne.0) then
          if(mod(nr,ntco).eq.0) call cpltwis(bezb(ix),t,etl,phi)
        endif
#ifdef ROOT
        if(root_flag .and. root_Optics.eq.1) then
          call OpticsRootWrite()
        end if
#endif

        cycle STRUCTLOOP

!--BETACALCULATION FOR SERIES OF BLOCKS (ix.ge.nblo.and.ithick.eq.1.and.iprint.eq.1)
  160   continue !if ithick=1 and iprint=1:
        if(ix.le.0) goto 190
!--REGULAR RUN THROUGH BLOCKS
        etl=etl+elbe(ix)

        do l=1,2
          ll=2*l

          if(abs(t(ll,ll-1)).gt.pieni) then
            phibf(l)=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
          else
            phibf(l)=zero
          endif

          puf=t(6,ll-1)
      t(6,ll-1)=(((((bl1eg(ix,l,1)*(t(1,ll-1)+puf*ded)+ bl1eg(ix,l,2)*(t&!hr06
     &(1,ll)+t(6,ll)*ded))+ bl1eg(ix,l,5)*dpp1*c1e3)- bl1(ix,l,1)*t     &!hr06
     &(1,ll-1))-bl1(ix,l,2)*t(1,ll))- bl1(ix,l,5)*dpr(1))/ded            !hr06
      t(6,ll)=(((((bl1eg(ix,l,3)*(t(1,ll-1)+puf*ded)+ bl1eg(ix,l,4)*(t  &!hr06
     &(1,ll)+t(6,ll)*ded))+ bl1eg(ix,l,6)*dpp1*c1e3)- bl1(ix,l,3)*t     &!hr06
     &(1,ll-1))-bl1(ix,l,4)*t(1,ll))- bl1(ix,l,6)*dpr(1))/ded            !hr06

          do i=1,ium-1
            puf=t(i,ll-1)
            t(i,ll-1)=(bl1(ix,l,1)*puf+bl1(ix,l,2)*t(i,ll))+dpr(i)*bl1  &!hr06
     &(ix,l,5)                                                           !hr06
        t(i,ll)=(bl1(ix,l,3)*puf+bl1(ix,l,4)*t(i,ll))+dpr(i)*bl1(ix,l,6) !hr06
          end do
        end do

        do l=1,2
          ll=2*l
          if(abs(t(ll,ll-1)).gt.pieni) then
            dphi=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
          else
            dphi=-one*phibf(l)                                           !hr06
          endif
          if(-one*dphi.gt.pieni) dphi=dphi+pi                            !hr06
          phi(l)=phi(l)+dphi/pie
        end do

        nr=nr+1
        call writelin(nr,bezb(ix),etl,phi,t,ix,.true.,k)
        if(ntco.ne.0) then
          if(mod(nr,ntco).eq.0) call cpltwis(bezb(ix),t,etl,phi)
        endif
#ifdef ROOT
        if(root_flag .and. root_Optics.eq.1) then
          call OpticsRootWrite()
        end if
#endif

        cycle STRUCTLOOP

!--REVERSE RUN THROUGH BLOCKS (ix.le.0)
  190   ix=-ix
        etl=etl+elbe(ix)
        do l=1,2
          ll=2*l

          if(abs(t(ll,ll-1)).gt.pieni) then
            phibf(l)=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
          else
            phibf(l)=zero
          endif

          puf=t(6,ll-1)
      t(6,ll-1)=(((((bl2eg(ix,l,1)*(t(1,ll-1)+puf*ded)+ bl2eg(ix,l,2)*(t&!hr06
     &(1,ll)+t(6,ll)*ded))+ bl2eg(ix,l,5)*dpp1*c1e3)- bl2(ix,l,1)*t     &!hr06
     &(1,ll-1))-bl2(ix,l,2)*t(1,ll))- bl2(ix,l,5)*dpr(1))/ded            !hr06
      t(6,ll)=(((((bl2eg(ix,l,3)*(t(1,ll-1)+puf*ded)+ bl2eg(ix,l,4)*(t  &!hr06
     &(1,ll)+t(6,ll)*ded))+ bl2eg(ix,l,6)*dpp1*c1e3)- bl2(ix,l,3)*t     &!hr06
     &(1,ll-1))-bl2(ix,l,4)*t(1,ll))- bl2(ix,l,6)*dpr(1))/ded            !hr06

          do i=1,ium-1
            puf=t(i,ll-1)
            t(i,ll-1)=(bl2(ix,l,1)*puf+bl2(ix,l,2)*t(i,ll))+dpr(i)*bl2  &!hr06
     &(ix,l,5)                                                           !hr06
        t(i,ll)=(bl2(ix,l,3)*puf+bl2(ix,l,4)*t(i,ll))+dpr(i)*bl2(ix,l,6) !hr06
          end do
        end do

        do l=1,2
          ll=2*l

          if(abs(t(ll,ll-1)).gt.pieni) then
            dphi=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
          else
            dphi=-phibf(l)
          endif

          if(-one*dphi.gt.pieni) dphi=dphi+pi                            !hr06
          phi(l)=phi(l)+dphi/pie
        end do

        nr=nr+1
        call writelin(nr,bezb(ix),etl,phi,t,ix,.true.,k)
        if(ntco.ne.0) then
          if(mod(nr,ntco).eq.0) call cpltwis(bezb(ix),t,etl,phi)
        endif
#ifdef ROOT
        if(root_flag .and. root_Optics.eq.1) then
          call OpticsRootWrite()
        end if
#endif

        cycle STRUCTLOOP

!--NOT A BLOCK / Nonlinear insertion
  220   ix=ix-nblo
        qu=zero
        qv=zero
        dyy1=zero
        dyy2=zero
        kpz=kp(ix)
        kzz=kz(ix)

 ! Cavity
        if( kpz.eq.6 .or.  abs(kzz).eq.12 ) then
          nr=nr+1
          call writelin(nr,bez(ix),etl,phi,t,ix,.false.,k)
          if(ntco.ne.0) then
            if(mod(nr,ntco).eq.0) call cpltwis(bez(ix),t,etl,phi)
          endif
#ifdef ROOT
          if(root_flag .and. root_Optics.eq.1) then
            call OpticsRootWrite()
          end if
#endif

          cycle STRUCTLOOP
        endif

        !Beam Beam element .and. fort.3 has BB block
        if(kzz.eq.20.and.nbeam.ge.1) then
          nbeam=k
          nr=nr+1
          call writelin(nr,bez(ix),etl,phi,t,ix,.false.,k)
          if(ntco.ne.0) then
            if(mod(nr,ntco).eq.0) call cpltwis(bez(ix),t,etl,phi)
          endif
#ifdef ROOT
          if(root_flag .and. root_Optics.eq.1) then
            call OpticsRootWrite()
          end if
#endif
          cycle STRUCTLOOP
        endif
        ! if kzz==22, starts a do over l; Update t matrix
        if(kzz == 22) then
          do l=1,2
            ll=2*l
            if(abs(t(ll,ll-1)).gt.pieni) then
              phibf(l)=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
            else
              phibf(l)=zero
            end if
            do i=1,ium
              puf=t(i,ll-1)
              t(i,ll-1)=(puf*rrtr(imtr(ix),ll-1,ll-1)+t(i,ll)*rrtr(imtr(ix),ll-1,ll))+dpr(i)*rrtr(imtr(ix),ll-1,6)
              t(i,ll)=(puf*rrtr(imtr(ix),ll,ll-1)+t(i,ll)*rrtr(imtr(ix),ll,ll))+dpr(i)*rrtr(imtr(ix),ll,6)
            end do
            t(1,ll-1)=t(1,ll-1)+cotr(imtr(ix),ll-1)
            t(1,ll)=t(1,ll)+cotr(imtr(ix),ll)
            if(abs(t(ll,ll-1)) > pieni) then
              dphi=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
            else
              dphi=-one*phibf(l)
            end if
            if(-one*dphi.gt.pieni) dphi=dphi+pi
            phi(l)=phi(l)+dphi/pie
          enddo
        endif

!+if collimat.or.bnlelens
        ! Marker, beam-beam, phase-trombone, crab cavity (incl. multipole), or wire
        if(kzz.eq.0.or.kzz.eq.20.or.kzz.eq.22                           &
     &     .or.abs(kzz).eq.23.or.abs(kzz).eq.26                         &
     &     .or.abs(kzz).eq.27.or.abs(kzz).eq.28                         &
     &     .or.abs(kzz).eq.15) then

          nr=nr+1
          call writelin(nr,bez(ix),etl,phi,t,ix,.false.,k)
          if(ntco.ne.0) then
            if(mod(nr,ntco).eq.0) call cpltwis(bez(ix),t,etl,phi)
          endif
#ifdef ROOT
          if(root_flag .and. root_Optics.eq.1) then
            call OpticsRootWrite()
          end if
#endif
          cycle STRUCTLOOP
        endif
!+ei
!+if .not.collimat.and..not.bnlelens
!        ! Marker, beam-beam or phase-trombone -> next element
!        if(kzz.eq.0.or.kzz.eq.20.or.kzz.eq.22) then
!           cycle STRUCTLOOP
!        endif
!        ! Wire -> next element
!        if(abs(kzz).eq.15) then
!           cycle STRUCTLOOP
!        endif
!        ! RF CC Multipoles -> next element
!        if (abs(kzz).eq.23.or.abs(kzz).eq.26.or.                        &
!     &       abs(kzz).eq.27.or.abs(kzz).eq.28) then
!           cycle STRUCTLOOP
!        endif
!+ei

        ! Update the matrix etc. for supported blocks
        dyy1=zero
        dyy2=zero
        if(iorg.lt.0) mzu(k)=izu
        izu=mzu(k)+1
        ekk=(sm(ix)+zfz(izu)*ek(ix))/(one+dpp)
        izu=izu+1
        xs=xpl(ix)+zfz(izu)*xrms(ix)
        izu=izu+1
        zs=zpl(ix)+zfz(izu)*zrms(ix)
#include "include/alignl.f90"

      if (kzz .ge. 0) then
        select case(kzz)

        case (1)
!--HORIZONTAL DIPOLE
           ekk=ekk*c1e3
#include "include/kickl01h.f90"
#include "include/kickq01h.f90"
!--NORMAL QUADRUPOLE
        case(2)
#include "include/kicklxxh.f90"
#include "include/kickq02h.f90"
!--   NORMAL SEXTUPOLE
        case(3)
           ekk=ekk*c1m3
#include "include/kickq03h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
!--NORMAL OCTUPOLE
        case(4)
           ekk=ekk*c1m6
#include "include/kicksho.f90"
#include "include/kickq04h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
!--NORMAL DECAPOLE
        case(5)
           ekk=ekk*c1m9
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq05h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
!--NORMAL DODECAPOLE
        case(6)
           ekk=ekk*c1m12
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq06h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
!--NORMAL 14-POLE
        case(7)
           ekk=ekk*c1m15
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq07h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
!--NORMAL 16-POLE
        case(8)
           ekk=ekk*c1m18
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq08h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
!--NORMAL 18-POLE
        case(9)
           ekk=ekk*c1m21
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq09h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
!--NORMAL 20-POLE
        case(10)
           ekk=ekk*c1m24
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq10h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
!--Multipole block
        case(11)
        r0=ek(ix)
        if(abs(dki(ix,1)).gt.pieni) then
          if(abs(dki(ix,3)).gt.pieni) then
#include "include/multl01.f90"
#include "include/multl08.f90"
            do 340 i=2,ium
#include "include/multl02.f90"
  340       continue
          else
#include "include/multl03.f90"
#include "include/multl09.f90"
          endif
        endif
        if(abs(dki(ix,2)).gt.pieni) then
          if(abs(dki(ix,3)).gt.pieni) then
#include "include/multl04.f90"
#include "include/multl10.f90"
            do 350 i=2,ium
#include "include/multl05.f90"
  350       continue
          else
#include "include/multl06.f90"
#include "include/multl11.f90"
          endif
        endif
        if(abs(r0).le.pieni) then
           cycle STRUCTLOOP
        endif
        nmz=nmu(ix)
        if(nmz.eq.0) then
          izu=izu+2*mmul

          nr=nr+1
          call writelin(nr,bez(ix),etl,phi,t,ix,.false.,k)
          if(ntco.ne.0) then
            if(mod(nr,ntco).eq.0) call cpltwis(bez(ix),t,etl,phi)
          endif
#ifdef ROOT
          if(root_flag .and. root_Optics.eq.1) then
            call OpticsRootWrite()
          end if
#endif

          cycle STRUCTLOOP
        endif
        im=irm(ix)
        r0a=one
        benkr=ed(ix)/(one+dpp)
        do 360 l=1,nmz
#include "include/multl07a.f90"
  360   continue
        if(nmz.ge.2) then
#include "include/multl07b.f90"
          do 365 l=3,nmz
#include "include/multl07c.f90"
  365     continue
        else
#include "include/multl07d.f90"
        endif
#ifdef TILT
#include "include/multl07e.f90"
#endif
        izu=izu+2*mmul-2*nmz

!--Skipped elements
        case(12,13,14,15,16,17,18,19,20,21,22,23)
           cycle STRUCTLOOP

!--DIPEDGE ELEMENT
        case(24)
#include "include/kickldpe.f90"
#include "include/kickqdpe.f90"
!--solenoid
        case(25)
#include "include/kicklso1.f90"
#include "include/kickqso1.f90"

!--Skipped elements
        case(26,27,28)
           cycle STRUCTLOOP

!--Unrecognized element (incl. cav with kp.ne.6 for non-collimat/bnlelens)
        case default
           nr=nr+1
           call writelin(nr,bez(ix),etl,phi,t,ix,.false.,k)
           if(ntco.ne.0) then
              if(mod(nr,ntco).eq.0) call cpltwis(bez(ix),t,etl,phi)
           endif
#ifdef ROOT
           if(root_flag .and. root_Optics.eq.1) then
             call OpticsRootWrite()
           end if
#endif
           cycle STRUCTLOOP
        end select


!--SKEW ELEMENTS
        else if(kzz .lt. 0) then
           kzz=-kzz             !Make it positive
           select case(kzz)
           case(1)
!--VERTICAL DIPOLE
              ekk=ekk*c1e3
#include "include/kickl01v.f90"
#include "include/kickq01v.f90"
!--SKEW QUADRUPOLE
           case(2)
#include "include/kicklxxv.f90"
#include "include/kickq02v.f90"
!--SKEW SEXTUPOLE
           case(3)
              ekk=ekk*c1m3
#include "include/kickq03v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
!--SKEW OCTUPOLE
           case(4)
              ekk=ekk*c1m6
#include "include/kicksho.f90"
#include "include/kickq04v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
!--SKEW DECAPOLE
           case(5)
              ekk=ekk*c1m9
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq05v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
!--SKEW DODECAPOLE
           case(6)
              ekk=ekk*c1m12
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq06v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
!--SKEW 14-POLE
           case(7)
              ekk=ekk*c1m15
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq07v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
!--SKEW 16-POLE
           case(8)
              ekk=ekk*c1m18
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq08v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
!--SKEW 18-POLE
           case(9)
              ekk=ekk*c1m21
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq09v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
!--SKEW 20-POLE
           case(10)
              ekk=ekk*c1m24
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq10v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"

!     Unrecognized skew element (including kzz=-12,kp.ne.6 for non-collimat/bnlelens)
           case default
              nr=nr+1
              call writelin(nr,bez(ix),etl,phi,t,ix,.false.,k)
              if(ntco.ne.0) then
                 if(mod(nr,ntco).eq.0) call cpltwis(bez(ix),t,etl,phi)
              endif
#ifdef ROOT
              if(root_flag .and. root_Optics.eq.1) then
                call OpticsRootWrite()
              end if
#endif
              cycle STRUCTLOOP
           end select
        endif

        !Done processing an element: go here!
        t(6,2)=t(6,2)-dyy1/(one+dpp)
        t(6,4)=t(6,4)-dyy2/(one+dpp)
        t(1,2)=t(1,2)+dyy1
        t(1,4)=t(1,4)+dyy2
        do i=2,ium
          if(kzz.eq.24) then
            t(i,2)=(t(i,2)+t(i,1)*qu)-qv*t(i,3)                          !hr06
            t(i,4)=(t(i,4)-t(i,3)*quz)-qvz*t(i,1)                        !hr06
!Contains elseif statements
#include "include/phas1so1.f90"
#include "include/phas2so1.f90"
#include "include/phas3so1.f90"
          else
            t(i,4)=(t(i,4)-t(i,3)*qu)-qv*t(i,1)                          !hr06
            t(i,2)=(t(i,2)+t(i,1)*qu)-qv*t(i,3)                          !hr06
          endif
        end do
        bexi=t(2,1)**2+t(3,1)**2                                         !hr06
        bezii=t(4,3)**2+t(5,3)**2                                        !hr06
        if(ncorru.eq.0) then
           if(kz(ix).eq.11) then
              if(abs(aa(2)).gt.pieni.and.nmz.gt.1) then
                 write(34,10070) etl,bez(ix),-2,aa(2),bexi,bezii,phi
              endif
              do iiii=3,nmz
                 if(abs(bb(iiii)).gt.pieni) then
                    write(34,10070) etl,bez(ix),iiii,bb(iiii),bexi,bezii,phi
                 endif
                 if(abs(aa(iiii)).gt.pieni) then
                    write(34,10070) etl,bez(ix),-iiii,aa(iiii),bexi,bezii,phi
                 endif
              enddo
           elseif(abs(ekk).gt.pieni.and.abs(kz(ix)).ge.3) then
              write(34,10070) etl,bez(ix),kz(ix),ekk,bexi,bezii,phi
           elseif(abs(ekk).gt.pieni.and.kz(ix).eq.-2) then
              write(34,10070) etl,bez(ix),kz(ix),ekk,bexi,bezii,phi
           endif
        endif

        nr=nr+1
        call writelin(nr,bez(ix),etl,phi,t,ix,.false.,k)
        if(ntco.ne.0) then
          if(mod(nr,ntco).eq.0) call cpltwis(bez(ix),t,etl,phi)
        endif
#ifdef ROOT
        if(root_flag .and. root_Optics.eq.1) then
          call OpticsRootWrite()
        end if
#endif

      end do STRUCTLOOP ! END LOOP OVER ELEMENTS

#ifdef HDF5
      if(h5_writeOptics) call h5lin_saveData
#endif

      call clorb(ded)
      do 510 l=1,2
        clo0(l)=clo(l)
        clop0(l)=clop(l)
  510 continue
      call clorb(zero)
      do 520 l=1,2
        ll=2*l
        di0(l)=(clo0(l)-clo(l))/ded
        dip0(l)=(clop0(l)-clop(l))/ded
  520 continue
      iiii=100
      idum='END'
      bexi=t(2,1)**2+t(3,1)**2                                           !hr06
      bezii=t(4,3)**2+t(5,3)**2                                          !hr06
      if(ncorru.eq.0) write(34,10070) etl,idum,iiii,zero,bexi,bezii,phi
      if(ncorru.eq.0) write(lout,10060)
!-----------------------------------------------------------------------
      return
10000 format(t5 ,'---- ENTRY LINOPT ----')
10010 format(132('-'))
10020 format('  NR     TYP      L-TOTAL    P     PHI          ',        &
     &'BETA         ALFA         GAMMA        DIS        DISP         ',&
     &'CLO        CLOP'/ 1x,                                            &
     &'                    (M)           (2*PI)        ',               &
     &'(M)          (RAD)         (M)         (M)        (RAD)        ',&
     &'(MM)       (MRAD)')
10030 format('  LINEAR OPTICS CALCULATION WITH PRINTOUT ',              &
     &'AFTER EACH BLOCK'/                                               &
     &'   A T T E N T I O N : BETATRON PHASE CALCULATION MIGHT BE WRONG'&
     &,' BY A MULTIPLE OF 0.5 FOR EACH LARGE BLOCK'/)
10040 format(/10x,'RELATIVE ENERGY DEVIATION  ',t40,f23.16/ 10x,         &
     &'TUNES -HORIZONTAL',t40,f23.16/ 10x,'      -VERTICAL  ',t40,f23.16/)
10050 format(t8,'  PLANE          DISP(MM)                 DISP(MRAD)'/ &
     &t6,'      X  ',2(f20.12,6x)/t10,'  Y  ',2(f20.12,6x)/)
10060 format(//131('-')//)
10070 format(1x,1pg21.14,1x,a,1x,i4,5(1x,1pg21.14))
end subroutine linopt

!-----------------------------------------------------------------------
!  WRITE OUT LINEAR OPTICS PARAMETERS AND IF COLLIMATION, SAVE STUFF.
!-----------------------------------------------------------------------
subroutine writelin(nr,typ,tl,p1,t,ixwl,isBLOC,ielem)
  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  use crcoall
  use parpro
  use mod_settings
  use mod_common
  use mod_commons
  use mod_commont

#ifdef ROOT
  use iso_c_binding, only: C_NULL_CHAR
  use root_output
#endif

#ifdef HDF5
  use hdf5_output
  use hdf5_linopt
#endif

  use collimation

  implicit none

  integer i,iwrite,ixwl,l,ll,nr
  real(kind=fPrec) al1,al2,b1,b2,c,cp,d,dp,g1,g2,p1,t,tl
  character(len=mNameLen) typ
  ! isBLOC.eq.TRUE if ixwl currently refers to a BLOC index, FALSE if it is a SINGLE ELEMENT index
  logical isBLOC
  dimension p1(2),t(6,4),b1(2),b2(2),al1(2),al2(2),g1(2),g2(2)
  dimension d(2),dp(2),c(2),cp(2)
  integer ielem

#ifdef HDF5
    real(kind=fPrec) hdf5Data(17)
#endif

  save
!-----------------------------------------------------------------------
  iwrite=0
  if(nlin.eq.0) then
    iwrite=1
  else
    do i=1,nlin
      if(typ.eq.bezl(i)) iwrite=1
    end do
  end if
  if(iwrite.eq.1) then
    do l=1,2
      ll=2*l
      b1(l)=t(ll,ll-1)**2+t(ll+1,ll-1)**2                            !hr06
      b2(l)=t(6-ll,ll-1)**2+t(7-ll,ll-1)**2                          !hr06
      al1(l)=-one*(t(ll,ll-1)*t(ll,ll)+t(ll+1,ll-1)*t(ll+1,ll))      !hr06
      al2(l)=-one*(t(6-ll,ll-1)*t(6-ll,ll)+t(7-ll,ll-1)*t(7-ll,ll))  !hr06
      g1(l)=t(ll,ll)**2+t(ll+1,ll)**2                                !hr06
      g2(l)=t(6-ll,ll)**2+t(7-ll,ll)**2                              !hr06
      d(l)=t(6,ll-1)*c1m3
      dp(l)=t(6,ll)*c1m3
      c(l)=t(1,ll-1)
      cp(l)=t(1,ll)
    end do

#ifdef ROOT
    if(root_flag .and. root_Optics.eq.1) then
      call OpticsRootWriteLin(nr, typ // C_NULL_CHAR,len(typ),tl,c(1),cp(1),c(2),cp(2),&
        b1(1),b1(2),al1(1),al1(2),d(1),d(2),dp(1),dp(2))
    end if
#endif
#ifdef HDF5
  if(h5_writeOptics) then
      hdf5Data(:) = (/tl,&
        p1(1),b1(1),al1(1),g1(1),d(1),dp(1),c(1),cp(1),&
        p1(2),b1(2),al1(2),g1(2),d(2),dp(2),c(2),cp(2)/)
      call h5lin_writeLine(nr, typ, hdf5Data)
    end if
#endif

    if (do_coll) then
      tbetax(max(ielem,1))  = b1(1)
      tbetay(max(ielem,1))  = b1(2)
      talphax(max(ielem,1)) = al1(1)
      talphay(max(ielem,1)) = al1(2)
      torbx(max(ielem,1))   = c(1)
      torbxp(max(ielem,1))  = cp(1)
      torby(max(ielem,1))   = c(2)
      torbyp(max(ielem,1))  = cp(2)
      tdispx(max(ielem,1))  = d(1)
      tdispy(max(ielem,1))  = d(2)
    endif

    if(ncorru == 0) then
      if(st_quiet == 0) then
        write(lout,10000) nr,typ(:8),tl,p1(1),b1(1),al1(1),g1(1),d(1),dp(1),c(1),cp(1)
        write(lout,10010) b2(1),al2(1),g2(1)
        write(lout,10030) typ(9:16)
        write(lout,10020) p1(2),b1(2),al1(2),g1(2),d(2),dp(2),c(2),cp(2)
        write(lout,10010) b2(2),al2(2),g2(2)
        write(lout,10040)
      end if
    else
      if(.not.isBLOC) then
        if(kp(ixwl).eq.3) then
          nhmoni=nhmoni+1
          betam(nhmoni,1)=b1(1)
          pam(nhmoni,1)=(p1(1)*two)*pi
          bclorb(nhmoni,1)=c(1)
        else if(kp(ixwl).eq.4) then
          nhcorr=nhcorr+1
          betac(nhcorr,1)=b1(1)
          pac(nhcorr,1)=(p1(1)*two)*pi
        else if(kp(ixwl).eq.-3) then
          nvmoni=nvmoni+1
          betam(nvmoni,2)=b1(2)
          pam(nvmoni,2)=(p1(2)*two)*pi
          bclorb(nvmoni,2)=c(2)
        else if(kp(ixwl).eq.-4) then
          nvcorr=nvcorr+1
          betac(nvcorr,2)=b1(2)
          pac(nvcorr,2)=(p1(2)*two)*pi
        end if
      end if
    end if
  end if
!-----------------------------------------------------------------------
  return
10010 format('|',6x,'|',8x,'|',12x,'|',1x,'|',12x,'|',f12.6,'|', f13.7,'|',f11.6,'|',11x,'|',11x,'|',11x,'|',11x,'|')
10020 format('|',6x,'|',8x,'|',12x,'|','Y','|',f12.7,'|',f12.6,'|', f13.7,'|',f11.6,'|',f11.7,'|',f11.7,'|',f11.7,'|',f11.7,'|')
10040 format(132('-'))
10000 format('|',i6,'|',a8,'|',f12.5,'|','X','|',f12.7,'|',f12.6,'|',f13.7,'|',f11.6,'|',f11.7,'|',f11.7,'|',f11.7,'|',f11.7,'|')
10030 format('|',6x,'|',a8,'|',12x,'|',102('-'))
end subroutine writelin

subroutine cpltwis(typ,t,etl,phi)
!-----------------------------------------------------------------------
!  CALCULATES COUPLED TWISS PARAMETERS AROUND THE RING AND ALSO THE
!  ANGLE OF THE MAJOR AXIS OF A ELLIPSE IN THE X-Y PROJECTION WITH
!  THE X-AXIS. THE 4-D ELLIPSOID IS GIVEN BY THE BOUNDARY OF A
!  DISTRIBUTION OF PARTICLES WITH MAXIMUM EMITANCE OF MODE I AND II,
!  EUI AND EUII RESPECTIVELY.
!  BINARY PRINT ON FILE 11 OF 22 VALUES :
!  POSITION [M],
!  BET(1-4), ALF(1-4), GAM(1-4), COOR-PHI(1-4), COOR-PRIME-PHI(1-4),
!  COUUANGL
!-----------------------------------------------------------------------
      use floatPrecision
#ifdef ROOT
      use root_output
#endif
      use numerical_constants
      use mathlib_bouncer
      use parpro
      use mod_common
      use mod_commons
      use mod_commont
      implicit none
      integer i,iwrite
      real(kind=fPrec) alxi,alxii,alzi,alzii,bexi,bexii,bezi,bezii,     &
     &couuang,etl,gaxi,gaxii,gazi,gazii,phi,phxi,phxii,phxpi,phxpii,    &
     &phzi,phzii,phzpi,phzpii,t
      character(len=mNameLen) typ
      dimension t(6,4),phi(2)
      save
!-----------------------------------------------------------------------
      iwrite=0
      if(nlin.eq.0) then
        iwrite=1
      else
        do 10 i=1,nlin
          if(typ.eq.bezl(i)) iwrite=1
   10   continue
      endif
      if(iwrite.eq.1) then
        bexi=t(2,1)**2+t(3,1)**2                                         !hr06
        bexii=t(4,1)**2+t(5,1)**2                                        !hr06
        bezi=t(2,3)**2+t(3,3)**2                                         !hr06
        bezii=t(4,3)**2+t(5,3)**2                                        !hr06
        alxi=-one*(t(2,1)*t(2,2)+t(3,1)*t(3,2))                          !hr06
        alxii=-one*(t(4,1)*t(4,2)+t(5,1)*t(5,2))                         !hr06
        alzi=-one*(t(2,3)*t(2,4)+t(3,3)*t(3,4))                          !hr06
        alzii=-one*(t(4,3)*t(4,4)+t(5,3)*t(5,4))                         !hr06
        gaxi=t(2,2)**2+t(3,2)**2                                         !hr06
        gaxii=t(4,2)**2+t(5,2)**2                                        !hr06
        gazi=t(2,4)**2+t(3,4)**2                                         !hr06
        gazii=t(4,4)**2+t(5,4)**2                                        !hr06
        if(abs(t(2,1)).gt.pieni) phxi=atan2_mb(t(3,1),t(2,1))
        if(abs(t(4,1)).gt.pieni) phxii=atan2_mb(t(5,1),t(4,1))
        if(abs(t(4,1)).gt.pieni) phxii=atan2_mb(t(5,1),t(4,1))
        if(abs(t(2,3)).gt.pieni) phzi=atan2_mb(t(3,3),t(2,3))
        if(abs(t(4,3)).gt.pieni) phzii=atan2_mb(t(5,3),t(4,3))
        if(abs(t(2,2)).gt.pieni) phxpi=atan2_mb(t(3,2),t(2,2))
        if(abs(t(4,2)).gt.pieni) phxpii=atan2_mb(t(5,2),t(4,2))
        if(abs(t(2,4)).gt.pieni) phzpi=atan2_mb(t(3,4),t(2,4))
        if(abs(t(4,4)).gt.pieni) phzpii=atan2_mb(t(5,4),t(4,4))
        if(abs(t(2,1)).le.pieni) phxi=pi*half
        if(abs(t(4,1)).le.pieni) then
          if(bexii.gt.pieni) phxii=pi*half
          if(bexii.le.pieni) phxii=zero
        endif
        if(abs(t(2,3)).le.pieni) then
          if(bezi.gt.pieni) phzi=pi*half
          if(bezi.le.pieni) phzi=zero
        endif
        if(abs(t(4,3)).le.pieni) phzii=pi*half
        if(abs(t(2,2)).le.pieni) phxpi=pi*half
        if(abs(t(4,2)).le.pieni) then
          if(gaxii.gt.pieni) phxpii=pi*half
          if(gaxii.le.pieni) phxpii=zero
        endif
        if(abs(t(2,4)).le.pieni) then
          if(gazi.gt.pieni) phzpi=pi*half
          if(gazi.le.pieni) phzpi=zero
        endif
        if(abs(t(4,4)).le.pieni) phzpii=pi*half
        if(abs(eui*(bexi-bezi)+euii*(bexii-bezii)).gt.pieni) then
          couuang=half*atan_mb((two*((eui*sqrt(bexi*bezi))*             &!hr06
     &cos_mb(phxi-phzi)+                                                &!hr06
     &(euii*sqrt(bexii*bezii))*cos_mb(phxii-phzii)))/ (eui*(bexi-bezi)  &!hr06
     &+euii*(bexii-bezii)))                                              !hr06
        else
          couuang=zero
        endif
        write(11,*) typ,etl,phi,bexi,bexii,bezi,bezii, alxi,alxii,alzi, &
     &alzii, gaxi,gaxii,gazi,gazii,phxi,phxii,phzi,phzii, phxpi,        &
     &phxpii,phzpi,phzpii,couuang,t(6,1),t(6,2),t(6,3),t(6,4),t(1,1),   &
     &t(1,2),t(1,3),t(1,4)

#ifdef ROOT
      if(root_flag .and. root_Optics.eq.1) then
        call OpticsRootWriteCpl(phi(1), phi(2),bexi,bexii,bezi,bezii,       &
 &                                  alxi,alxii,alzi,alzii,       &
 &                                  gaxi,gaxii,gazi,gazii,       &
 &                                  phxi,phxii,phzi,phzii,       &
 &                                  phxpi,phxpii,phzpi,phzpii,   &
 &                                  couuang,                     &
 &                                  t(6,1),t(6,2),t(6,3),t(6,4), &
 &                                  t(1,1),t(1,2),t(1,3),t(1,4))
      end if
#endif

      endif
      return
end subroutine cpltwis

subroutine loesd (rmat, vec,dimakt,dimtot,kod)
!-----------------------------------------------------------------------
!  SOLUTION OF A SYSTEM OF LINEAR EQUATIONS
!  VEC1 = VEC2 * RMAT , WITH VEC2 AS RESULT
!-----------------------------------------------------------------------
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      use crcoall
      use parpro
      implicit none
      integer ik,indi,j,jk,jy,k,kk,kod,l,n,n1,dimtot,dimakt
      real(kind=fPrec) emax,eps,r,rmat,vec
      dimension rmat(dimtot,dimakt),vec(dimakt)
      data eps /1e-20_fPrec/
      save
!-----------------------------------------------------------------------
      kod=1
      do j=1,dimakt
        emax=zero
        do 10 ik=j,dimakt
          if(abs(emax).gt.abs(rmat(j,ik)) .or.emax.ne.emax) goto 10
          emax=rmat(j,ik)
          indi=ik
   10   continue
        if(abs(emax).lt.eps) then
          write(lout,*) '  ****   ERROR IN LOESD   **** '
          return
        endif

   20   do l=j,dimakt
          r=rmat(l,j)
          rmat(l,j)=rmat(l,indi)
          rmat(l,indi)=r
          rmat(l,j)=rmat(l,j)/emax
        end do

        r=vec(indi)
        vec(indi)=vec(j)
        vec(j)=r/emax
        if(j.eq.dimakt) goto 60
        jy=j+1

        do jk=jy,dimakt
          r=rmat(j,jk)

          do kk=jy,dimakt
            rmat(kk,jk)= rmat(kk,jk)-r*rmat(kk,j)
          end do
          vec(jk)=vec(jk)-vec(j)*r
        end do
      end do

   60 n=dimakt
      n1=dimakt-1

      do j=1,n1
        do k=1,j
          vec(n-j)=vec(n-j)-rmat(n-k+1,n-j)*vec(n-k+1)
        end do
      end do

      kod = 0
      return
end subroutine loesd

subroutine matrix(dpp,am)
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      use parpro
      use mod_common
      use mod_commons
      use mod_commont
      implicit none
      integer i,ierr,l
      real(kind=fPrec) am,dpp
      dimension am(4,4)
      save
!-----------------------------------------------------------------------
      do i=2,5
        do l=1,2
           x(i,l)=zero
           y(i,l)=zero
        end do
      end do

      x(2,1)=one
      y(3,1)=one
      x(4,2)=one
      y(5,2)=one

      do l=1,2
        x(1,l)=clo(l)
        y(1,l)=clop(l)
      end do

      call umlauf(dpp,5,ierr)
      ierro=ierr

      do i=1,4
        am(1,i)=x(i+1,1)
        am(2,i)=y(i+1,1)
        am(3,i)=x(i+1,2)
        am(4,i)=y(i+1,2)
      end do
!-----------------------------------------------------------------------
      return
end subroutine matrix

subroutine corrorb
!-----------------------------------------------------------------------
!  CORRECTION OF CLOSED ORBIT FIRST (MOST EFFECTIV CORRECTOR STRATEGY
!  USING MICADO), THEN
!  SCALING OF DIPOLE-ERRORS FOR RMS-VALUES OF THE CLOSED ORBIT
!-----------------------------------------------------------------------
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      use crcoall
      use parpro
      use mod_common
      use mod_commons
      use mod_commont
      implicit none
      integer i,icflag,ihflag,ii,ij,im,iprinto,ivflag,j,k,kpz,kzz,l,nlino,ntcoo,nto,nx
      real(kind=fPrec) ar(nmon1,ncor1)
      real(kind=fPrec) b(nmon1),orbr(nmon1),xinc(ncor1)
      real(kind=fPrec) rmsx,ptpx,rmsz,ptpz,rzero,rzero1
      real(kind=fPrec) clo0,clop0,hfac,qwc1,vfac
      character(len=mNameLen) bezlo(nele)
      dimension clo0(2),clop0(2)
      dimension qwc1(3),nx(ncor1)
      save
!-----------------------------------------------------------------------
      rzero=zero
      rzero1=zero
      do l=1,2
        clo0(l)=zero
        clop0(l)=zero
        di0(l)=zero
       dip0(l)=zero
      end do

      call clorb(ded)
      if(ierro.gt.0) then
        write(lout,"(a)") "CLORB> ERROR Unstable closed orbit during initial dispersion calculation."
        write(lout,"(a)") "CLORB>       Instability occurred for small relative energy deviation."
        call prror(-1)
      end if

      do l=1,2
        clo0(l)=clo(l)
        clop0(l)=clop(l)
      end do

      call clorb(zero)
      if(ierro.gt.0) then
        write(lout,"(a)") "CLORB> ERROR Unstable closed orbit for zero energy deviation."
        call prror(-1)
      end if

      do l=1,2
        di0(l)=(clo0(l)-clo(l))/ded
        dip0(l)=(clop0(l)-clop(l))/ded
      end do

      do l=1,ncor1
        xinc(l)=zero
        nx(l)=0
      end do

      if(iclo.eq.0) return

!-- ORBIT CORRECTION
      ihflag=0
      ivflag=0
      icflag=0

      write(lout,*)
      write(lout,10000)

      if(ncorru == 0) then
        write(lout,"(a)") "CLORB> ERROR Number of orbit correctors is zero."
        call prror(-1)
      else
        if(ncorrep.le.0) then
          write(lout,10010) ncorru,sigma0(1),sigma0(2)
        else
          write(lout,10020) ncorru,ncorrep
        endif
      endif

      write(lout,*)

!-- SAVE OLD 'LINOPT' SETTINGS
      iprinto=iprint
      nto=nt
      ntcoo=ntco

      do i=1,nlin
        bezlo(i)=bezl(i)
      end do

      nlino=nlin

!-- PUT MONITORS AND CORRECTORS INTO LINOPT SETTINGS
!-- GET TWISS PARAMETERS AND DISTORTED ORBIT BACK
      iprint=0
      ntco=0
      nlin=0

      do i=1,il
        if(kp(i).eq.3.or.kp(i).eq.4.or. kp(i).eq.-3.or.kp(i).eq.-4) bezl&
     &(i)=bez(i)
        nlin=nlin+1
      end do

      call linopt(zero)
      call phasad(zero,qwc1)

!-- CHECK SOME CONDITIONS
      write(lout,10100) nhmoni,nhcorr,nvmoni,nvcorr
      if(nhmoni.gt.nmon1) then
        write(lout,10070) nhmoni,nmon1
        return
      endif
      if(nvmoni.gt.nmon1) then
        write(lout,10070) nvmoni,nmon1
        return
      endif
      if(nhcorr.gt.ncor1) then
        write(lout,10080) nhcorr,ncor1
        return
      endif
      if(nvcorr.gt.ncor1) then
        write(lout,10080) nvcorr,ncor1
        return
      endif
      if(nhmoni.lt.nhcorr.or.nvmoni.lt.nvcorr) write(lout,10090)

      write(lout,*)
      call orbinit

!-- CORRECT BOTH PLANES
      if(ncorrep.eq.0) then
        icflag=1
        ncorrep=itco
      endif
      do 110 ii=1,ncorrep

!-- HORIZONTAL PLANE FIRST
        do i=1,nhmoni
          b(i)=bclorb(i,1)
          do j=1,nhcorr
      ar(i,j)=((sqrt(betam(i,1)*betac(j,1))*cos_mb(abs(pam(i,1)-pac     &!hr06
     &(j,1))-qwc1(1)*pi))*c1e3)/(two*sin_mb(qwc1(1)*pi))                 !hr06
          end do
        end do

        call calrms(b,nhmoni,rmsx,ptpx)

!-- MICADO WITH HOUSEHOLDER TRANSFORMATION
        call htls(ar,b,nhmoni,nhcorr,xinc,nx,orbr,ncorru,rzero,rzero1)

!-- VERTICAL PLANE HERE
        do i=1,nvmoni
          b(i)=bclorb(i,2)                                               !hr06
          do j=1,nvcorr
      ar(i,j)=((sqrt(betam(i,2)*betac(j,2))*cos_mb(abs(pam(i,2)-pac     &!hr06
     &(j,2))-qwc1(2)*pi))*c1e3)/(two*sin_mb(qwc1(2)*pi))                 !hr06
          end do
        end do

        call calrms(b,nvmoni,rmsz,ptpz)
        write(lout,10030) ii-1,rmsx,rmsz
        write(lout,10040) ii-1,ptpx,ptpz
        if(icflag.eq.1.and.sigma0(1).gt.rmsx.and.ihflag.eq.0) then
          write(lout,10110)
          ihflag=1
        endif
        if(icflag.eq.1.and.sigma0(2).gt.rmsz.and.ivflag.eq.0) then
          write(lout,10120)
          ivflag=1
        endif

        if(ihflag.eq.0) then
          write(lout,*)

          do ij=1,ncorru/10
            write(lout,10050) (nx(10*(ij-1)+k), k=1,10)
          end do

          if(mod(ncorru,10).gt.0) then
            write(lout,10050) (nx(10*(ij-1)+k), k=1,mod(ncorru,10))
          endif
          call putorb(xinc,nx,1)
        endif

!-- MICADO WITH HOUSEHOLDER TRANSFORMATION
        call htls(ar,b,nvmoni,nvcorr,xinc,nx,orbr,ncorru,rzero,rzero1)

        if(ivflag.eq.0) then
          write(lout,*)
          do 100 ij=1,ncorru/10
            write(lout,10060) (nx(10*(ij-1)+k), k=1,10)
  100     continue

          if(mod(ncorru,10).gt.0) then
            write(lout,10060) (nx(10*(ij-1)+k), k=1,mod(ncorru,10))
          endif
          call putorb(xinc,nx,2)
        endif

        if(ihflag.eq.1.and.ivflag.eq.1) goto 140
        call linopt(zero)
        call phasad(zero,qwc1)
  110 continue

!-- GET LAST VALUES AFTER CORRECTION
      do 120 i=1,nhmoni
        b(i)=bclorb(i,1)                                                 !hr06
  120 continue

      call calrms(b,nhmoni,rmsx,ptpx)

      do 130 i=1,nvmoni
        b(i)=bclorb(i,2)                                                 !hr06
  130 continue

      call calrms(b,nvmoni,rmsz,ptpz)
      write(lout,10030) ncorrep,rmsx,rmsz
      write(lout,10040) ncorrep,ptpx,ptpz
      write(lout,*)

  140 continue
      if((ii-1).eq.itco) write(lout,10130) itco

!-- SCALE TO DESIRED RMS VALUE IF IT IS GREATER THAN ZERO
      if(sigma0(1).gt.pieni.or.sigma0(2).gt.pieni) then
        do 180 ii=1,itco
          write(lout,10140)
          hfac=sigma0(1)/rmsx                                            !hr06
          vfac=sigma0(2)/rmsz                                            !hr06
          do 150 i=1,il
            kzz=kz(i)
            kpz=kp(i)
            if(kzz.eq.1.and.el(i).lt.pieni) then
              ed(i)=ed(i)*hfac
              ek(i)=ek(i)*hfac
            endif
            if(kzz.eq.-1.and.el(i).lt.pieni) then
              ed(i)=ed(i)*vfac
              ek(i)=ek(i)*vfac
            endif
            if(kzz.eq.11) then
              im=irm(i)
              ak0(im,1)=ak0(im,1)*vfac
              aka(im,1)=aka(im,1)*vfac
              bk0(im,1)=bk0(im,1)*hfac
              bka(im,1)=bka(im,1)*hfac
            endif
  150     continue
          call linopt(zero)

          do 160 i=1,nhmoni
            b(i)=bclorb(i,1)                                             !hr06
  160     continue

          call calrms(b,nhmoni,rmsx,ptpx)

          do 170 i=1,nvmoni
            b(i)=bclorb(i,2)                                             !hr06
  170     continue

          call calrms(b,nvmoni,rmsz,ptpz)

          write(lout,10150) ii,rmsx,rmsz
          write(lout,10160) ii,ptpx,ptpz
          write(lout,*)
          if(abs(real(rmsx,fPrec)-sigma0(1)).lt.dsi.and.                      &!hr06
     &       abs(real(rmsz,fPrec)-sigma0(2)).lt.dsi)                          &!hr06
     &goto 190
  180   continue
      endif

      if((ii-1).eq.itco) write(lout,10130) itco
  190 continue

!-- WRITE OUT ADJUSTED CLOSED ORBIT
      do 200 i=1,nhmoni
        write(28,*) i,bclorb(i,1)
  200 continue

      do 210 i=1,nhmoni
        write(29,*) i,bclorb(i,2)
  210 continue

!-- CHANGE BACK TO OLD 'LINOPT' SETTINGS
      iprint=iprinto
      nt=nto
      ntco=ntcoo
      nlin=nlino

      do 220 i=1,nlin
        bezl(i)=bezlo(i)
  220 continue

      ncorru=0
!-----------------------------------------------------------------------
      return
10000 format(t5,'---- ORBIT CORRECTION WITH MOST EFFCTIVE CORRECTOR ',  &
     &'STRATEGY ----')
10010 format(t5,'     ORBIT CORRECTION WITH ',i4,' CORRECTORS UNTIL',/, &
     &t5,'       HOR. RMS SMALLER THAN ',f6.3,' MM',/, t5,              &
     &'       VER. RMS SMALLER THAN ',f6.3,' MM')
10020 format(t5,'     ORBIT CORRECTION WITH ',i4,' CORRECTORS AND ',i4, &
     &' ITERATIONS.')
10030 format(t5,'---- CORRECTION ITERATION NO. ',i4,' HOR.-RMS: ',f6.3, &
     &' VER.-RMS: ',f6.3)
10040 format(t5,'---- CORRECTION ITERATION NO. ',i4,' HOR.-PTP: ',f6.3, &
     &' VER.-PTP: ',f6.3)
10050 format(t5,'     HORIZONTAL CORRECTORS USED:', i4,i4,i4,i4,i4,i4,  &
     &i4,i4,i4,i4)
10060 format(t5,'     VERTICAL   CORRECTORS USED:', i4,i4,i4,i4,i4,i4,  &
     &i4,i4,i4,i4)
10070 format(/,t5,'ERROR: NUMBER OF MONITORS TOO BIG.',/                &
     &'    THERE ARE ',i4,' MONITORS SET, BUT ONLY ',i4, ' ALLOWED.',/  &
     &'    NO CORRECTION DONE.',/)
10080 format(/,t5,'ERROR: NUMBER OF CORRECTORS TOO BIG.',/              &
     &'    THERE ARE ',i4,' MONITORS SET, BUT ONLY ',i4, ' ALLOWED.',/  &
     &'    NO CORRECTION DONE.',/)
10090 format(/,t5,'WARNING: NUMBER OF MONITORS IS SMALLER THAN NUMBER', &
     &' OF CORRECTORS.',/ '    NUMERICAL PROBLEMS MIGHT BE ENCOUNTERED.'&
     &)
10100 format(/,t5,'NUMBER OF HOR. MONITORS: ',i4,                       &
     &'  NUMBER OF HOR. CORRECTORS: ',i4,/, t5,                         &
     &'NUMBER OF VER. MONITORS: ',i4, '  NUMBER OF VER. CORRECTORS: ',  &
     &i4)
10110 format(t10,'HORIZONTAL RMS GOAL REACHED')
10120 format(t10,'VERTICAL RMS GOAL REACHED')
10130 format(t10,'MAXIMUM NUMBER OF ITERATIONS ACHIVED: ',i4,/ ,t10,    &
     &'INCREASE ITCO TO INCREASE THE NUMBER OF ' ,                      &
     &'CLOSED ORBIT ITERATIONS',/)
10140 format(t5,'---- ORBIT SCALING USING ALL POSSIBLE ELEMENTS ')
10150 format(t5,'---- SCALING ITERATION NO. ',i4,' HOR.-RMS: ',f6.3,    &
     &' VER.-RMS: ',f6.3)
10160 format(t5,'---- SCALING ITERATION NO. ',i4,' HOR.-PTP: ',f6.3,    &
     &' VER.-PTP: ',f6.3)
end subroutine corrorb

      subroutine putorb(xinc,nx,npflag)
!-----------------------------------------------------------------------
!  PUT ORBIT CHANGES FROM MICADO TO THE GIVEN ORBIT CORRECTORS
!-----------------------------------------------------------------------
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      use crcoall
      use parpro
      use mod_common
      use mod_commons
      use mod_commont
      implicit none

      integer i,im,ix,izu,j,k,kcorr,kcorru,kpz,kzz,nmz,npflag,nx
      real(kind=fPrec) xinc(ncor1)
      real(kind=fPrec) ckicknew,ckickold,r0,r0a
      dimension nx(ncor1)
      save
!-----------------------------------------------------------------------
      kcorru=0
      kcorr=0
      izu=0

      do 60 i=1,iu
        ix=ic(i)
        if(ix.le.nblo) goto 60
        ix=ix-nblo
        kpz=kp(ix)
        kzz=kz(ix)
        if(kpz.eq.6.or.kzz.eq.0.or.kzz.eq.20.or.kzz.eq.22) goto 60
        if(kzz.eq.15) goto 60
        if(iorg.lt.0) mzu(i)=izu
        izu=mzu(i)+1
        if(kpz.eq.4.and.kzz.eq.1.and.npflag.eq.1.or.                    &
     &kpz.eq.-4.and.kzz.eq.-1.and.npflag.eq.2) then
          kcorr=kcorr+1
          do 10 j=1,ncorru
            if(nx(j).eq.kcorr) then
              kcorru=kcorru+1
              ckickold=sm(ix)+zfz(izu)*ek(ix)
              zfz(izu)=zfz(izu)+real(xinc(j),fPrec)/ek(ix)                     !hr06
              ckicknew=sm(ix)+zfz(izu)*ek(ix)
              write(lout,10000) kcorru,kcorr,bez(ix), ckickold*c1e3,    &
     &ckicknew*c1e3
            endif
   10     continue
        endif
        izu=izu+2

        if(kzz.eq.11) then
          r0=ek(ix)
          if(abs(r0).le.pieni) goto 60
          nmz=nmu(ix)
          if(nmz.eq.0) then
            izu=izu+2*mmul
            goto 60
          endif
          im=irm(ix)
          r0a=one
          do 50 k=1,nmz
            izu=izu+1
            if(kpz.eq.-4.and.npflag.eq.2.and.k.eq.1) then
              kcorr=kcorr+1
              do 30, j=1,ncorru
                if(nx(j).eq.kcorr) then
                  kcorru=kcorru+1
                  ckickold=ed(ix)*(ak0(im,k)+zfz(izu)* aka(im,k))/r0a
           zfz(izu)=zfz(izu)+(c1e3*(real(xinc(j),fPrec)/(r0a*ed(ix))-ak0&!hr06
     &(im,k)))/aka(im,k)                                                 !hr06
                  ckicknew=(ed(ix)*(ak0(im,k)+zfz(izu)* aka(im,k)))/r0a  !hr06
                  write(lout,10000) kcorru,kcorr,bez(ix), ckickold,     &
     &ckicknew
                endif
   30         continue
            endif
            izu=izu+1
            if(kpz.eq.4.and.npflag.eq.1.and.k.eq.1) then
              kcorr=kcorr+1
              do 40, j=1,ncorru
                if(nx(j).eq.kcorr) then
                  kcorru=kcorru+1
                  ckickold=(ed(ix)*(bk0(im,k)+zfz(izu)* bka(im,k)))/r0a  !hr06
           zfz(izu)=zfz(izu)+(c1e3*(real(xinc(j),fPrec)/(r0a*ed(ix))-bk0&!hr06
     &(im,k)))/bka(im,k)                                                 !hr06
                  ckicknew=(ed(ix)*(bk0(im,k)+zfz(izu)* bka(im,k)))/r0a  !hr06
                  write(lout,10000) kcorru,kcorr,bez(ix), ckickold,     &
     &ckicknew
                endif
   40         continue
            endif
   50     continue
          izu=izu+2*mmul-2*nmz
        endif
   60 continue
!-----------------------------------------------------------------------
      return
10000 format(t5,i4,i4,' ',a16,'  OLD: ',d14.7,' MRAD   NEW: ' ,d14.7,   &
     &' MRAD')
      end

      subroutine orbinit
!-----------------------------------------------------------------------
!  INITIALIZES THE RANDOM NUMBER OF NOT SET CORRCTORS
!-----------------------------------------------------------------------
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      use parpro
      use mod_common
      use mod_commons
      use mod_commont
      implicit none
      integer i,im,ix,izu,kpz,kzz,nmz
      real(kind=fPrec) r0
      save
!-----------------------------------------------------------------------
      izu=0
      do 10 i=1,iu
        ix=ic(i)
        if(ix.le.nblo) goto 10
        ix=ix-nblo
        kpz=kp(ix)
        kzz=kz(ix)
        if(kpz.eq.6.or.kzz.eq.0.or.kzz.eq.20.or.kzz.eq.22) goto 10
        if(kzz.eq.15) goto 10
        if(iorg.lt.0) mzu(i)=izu
        izu=mzu(i)+1
        if((kpz.eq.4.and.kzz.eq.1).or.(kpz.eq.-4.and.kzz.eq.-1)) then
          zfz(izu)=zero
          ek(ix)=one
          ncororb(ix)=1
        endif
        izu=izu+2

        if(kzz.eq.11) then
          r0=ek(ix)
          if(abs(r0).le.pieni) goto 10
          nmz=nmu(ix)
          if(nmz.eq.0) then
            izu=izu+2*mmul
            goto 10
          endif
          im=irm(ix)

          izu=izu+1
          if(kpz.eq.-4) then
            zfz(izu)=zero
            aka(im,1)=one
          endif
          izu=izu+1
          if(kpz.eq.4) then
            zfz(izu)=zero
            bka(im,1)=one
          endif
          izu=izu+2*mmul-2
        endif
   10 continue
      return
      end

subroutine htls(a,b,m,n,x,ipiv,r,iter,rms,ptp)
!*********************************************************************
!     Subroutine HTLS to make Householder transform                  *
!                                                                    *
!     Authors:     many                Date:  17.09.1989             *
!                                                                    *
!     DIMENSION OF ARRAY RHO SHOULD BE 3*NCOR1                       *
!     M    - NUMBER OF AVAILABLE MONITORS                            *
!     N    - NUMBERR OF AVAILABLE INDEPENDENT CORRECTORS             *
!     ITER - NUMBER OF CORRECTORS TO BE USED                         *
!     RMS  - RMS VALUE TO CORRECT FOR                                *
!     PTP  - PEAK TO PEAK VALUE TO CORRECT FOR                       *
!*********************************************************************
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      use crcoall
      implicit none
      integer i,iii,ij1,ip,ipiv,iter,j,j1,k,k2,k3,ki,kk,kpiv,m,n,ncor1, &
     &nmon1
      real(kind=fPrec) a,b,piv,pivt,ptop,r,rho,rmss,x,xiter,xptp,xrms
      real(kind=fPrec) rms,ptp
      real(kind=fPrec) g,h,sig,beta
      parameter (nmon1 = 600)
      parameter (ncor1 = 600)
      dimension a(nmon1,ncor1),b(nmon1),x(ncor1),ipiv(ncor1),r(nmon1)
      dimension rho(3*ncor1),xiter(ncor1),xrms(ncor1),xptp(ncor1)
      dimension rmss(ncor1),ptop(ncor1)
      save
!-----------------------------------------------------------------------

! --- calcul du premier pivot

!============================
      beta=zero

      do ij1=1,500
        rho(ij1)=zero
      end do

      k2=n + 1
      piv=zero

      do k=1,n
        ipiv(k)=k
        h=zero                                                           !hr06
        g=zero                                                           !hr06

        do i=1,m
          h=h+a(i,k)*a(i,k)
          g=g+a(i,k)*b(i)
        end do

        rho(k)=h
        rho(k2) = g
        pivt = g**2/h                                                    !hr06
        if(pivt-piv.le.0) goto 40
        if(pivt-piv.gt.0) goto 30
   30   piv = pivt

        kpiv=k
        k2 = k2 + 1
   40   continue
      end do

! --- boucle pour chaque iteration

      do 150 k=1,iter
        if(kpiv.eq.k)goto 60

! --- on echange les K et KPIV si KPIV plus grand que K
        h=rho(k)
        rho(k)=rho(kpiv)
        rho(kpiv)=h
        k2=n+k
        k3=n+kpiv
        g = rho(k2)
        rho(k2) = rho(k3)
        rho(k3) = g
        do i=1,m
          h=a(i,k)
          a(i,k)=a(i,kpiv)
          a(i,kpiv)=h
        end do

! --- calcul de beta,sigma et uk dans htul
   60   continue
        call htul(a,m,k,sig,beta)

! --- on garde SIGMA dans RHO(N+K)
        j=n+k
        rho(j)=-one*sig                                                  !hr06
        ip=ipiv(kpiv)
        ipiv(kpiv)=ipiv(k)
        ipiv(k)=ip
        if(k.eq.n) goto 70

! --- transformation de A dans HTAL
        call htal(a,m,n,k,beta)

! --- transformation de B dans HTBL
   70   continue
        call htbl(a,b,m,k,beta)

! --- recherche du pivot (K+1)
!=============================

        rho(k)=sqrt(piv)
        if(k.eq.n) goto 90
        piv=zero                                                          !hr06
        kpiv = k + 1
        j1 = kpiv
        k2=n + j1

        do j=j1,n
          h=rho(j)-(a(k,j))*(a(k,j))

          if(h.lt.c1m7) then
            write(lout,*)
            write(lout,*) 'CORRECTION PROCESS ABORTED.'
            write(lout,*) 'DIVISION BY ZERO EXPECTED.'
            write(lout,*) 'PROBABLY TWO CORRECTORS TOO CLOSE.'
            write(lout,10000) ' SUSPECTED CORRECTOR: ',j
            write(lout,'(a)') "Error '777' in subroutine htls"
            call prror(-1)
          endif

          rho(j)=h
          g=rho(k2)-(a(k,j))*(b(k))
          rho(k2) = g
          pivt = g**2/h                                                  !hr06
          if(pivt.lt.piv)goto 80
          kpiv=j
          piv=pivt
          k2 = k2 + 1
   80     continue
        end do

! --- calcul des X
   90   x(k)=b(k)/rho(n+k)
        if(k.eq.1)goto 120
        do i=2,k
          kk=k-i+1
          x(kk)=b(kk)
          ki=kk+1
          do j=ki,k
            x(kk)=x(kk)-a(kk,j)*x(j)
          end do
          x(kk)=x(kk)/rho(n+kk)
        end do
  120   continue

! --- save residual orbit and inverse sign of corrections (convention!)
        do iii= 1,m
          r(iii) = b(iii)
        end do
        do iii= 1,k
          x(iii) =-one*x(iii)                                           !hr06
        end do

! --- calcul du vecteur residuel dans HTRL
!=========================================

!     transform orbit R back to "normal space"
        call htrl(a,r,m,n,k,rho)
        call calrms(r,m,rmss(k),ptop(k))
        xiter(k+1) = k
        xrms(k+1) = rmss(k)
        xptp(k+1) = ptop(k)

        if(ptop(k).le.ptp)goto 160
        if(rmss(k).le.rms)goto 160
  150 continue
      return

! --- correction is already good enough:
!=======================================

  160 ptp=ptop(k)
      rms=rmss(k)
10000 format(a,i4)
end subroutine htls

!*********************************************************************
!     Subroutine HTAL to make Householder transform                  *
!                                                                    *
!     Authors:     many                Date:  17.09.1989             *
!                                                                    *
!     Householder transform of matrix A
!*********************************************************************
subroutine htal(a,m,n,k,beta)
  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  implicit none
  integer j,k,k1,m,n,nc,ncor1,nmon1
  real(kind=fPrec) a,beta,h
  parameter (nmon1 = 600)
  parameter (ncor1 = 600)
  dimension a(nmon1,ncor1)
  save
  nc=n-k
  do j=1,nc
    h=zero
    do k1=k,m
      h=h+a(k1,k)*a(k1,k+j)
    end do
    h=beta*h
    do k1=k,m
      a(k1,k+j)=a(k1,k+j)-a(k1,k)*h
    end do
  end do
end subroutine htal

!*********************************************************************
!     Subroutine HTBL to make Householder transform                  *
!                                                                    *
!     Authors:     many                Date:  17.09.1989             *
!                                                                    *
!     Householder transform of vector B
!*********************************************************************
subroutine htbl(a,b,m,k,beta)
  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  implicit none
  integer k,k1,m,ncor1,nmon1
  real(kind=fPrec) a,b,beta,h
  parameter (nmon1 = 600)
  parameter (ncor1 = 600)
  dimension a(nmon1,ncor1),b(nmon1)
  save
  h=zero
  do k1=k,m
    h=h+a(k1,k)*b(k1)
  end do
  h=beta*h
  do k1=k,m
    b(k1)=b(k1)-a(k1,k)*h
  end do
end subroutine htbl

!*********************************************************************
!     Subroutine HTRL to make Householder transform                  *
!                                                                    *
!     Authors:     many                Date:  17.09.1989             *
!                                                                    *
!     calculate residual orbit vector
!*********************************************************************
subroutine htrl(a,b,m,n,k,rho)
  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  implicit none
  integer i,k,kk,kl,kn,lv,m,n,ncor1,nmon1
  real(kind=fPrec) a,b,beta,rho
  parameter (nmon1 = 600)
  parameter (ncor1 = 600)
  dimension a(nmon1,ncor1),b(nmon1),rho(3*ncor1)
  save
  do i= 1,k,1
    b(i)= zero
  end do
  do kk=1,k
    lv=m-k+kk
    kn=n+k-kk+1
    kl=k-kk+1
    beta=-one/(rho(kn)*a(kl,kl))
    call htbl(a,b,m,kl,beta)
  end do
end subroutine htrl

!*********************************************************************
!     Subroutine HTUL to make Householder transform                  *
!                                                                    *
!     Authors:     many                Date:  17.09.1989             *
!                                                                    *
!     calculate vector U
!*********************************************************************
subroutine htul(a,m,k,sig,beta)
  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  implicit none
  integer i,k,m,ncor1,nmon1
  real(kind=fPrec) a,beta,h,sig
  parameter (nmon1 = 600)
  parameter (ncor1 = 600)
  dimension a(nmon1,ncor1)
  save
  sig=zero
  do i=k,m
    sig=sig+a(i,k)* a(i,k)
  end do
  sig=sqrt(sig)
  ! on choisit le signe correct pour SIG:
  h=a(k,k)
  if(h.lt.zero)sig=-one*sig
  beta=h + sig
  a(k,k)=beta
  beta=one/(sig*beta)
end subroutine htul

!*********************************************************************
!     Subroutine CALRMS to calculate rms                             *
!                                                                    *
!     Authors:     many                Date:  17.09.1989             *
!                                                                    *
!     calculates rms and p.to.p value of R(1) .... R(M)
!*********************************************************************
subroutine calrms(r,m,rms,ptp)
  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  implicit none
  integer i,imax,imin,m,maxmin
  real(kind=fPrec) ave,ptp,r,rms,xave,xrms
  dimension r(m)
  save
  xave = zero
  xrms = zero
  do i=1,m
    xave = xave + r(i)
    xrms = xrms + r(i)**2
  end do
  ave = xave / real(m,fPrec)
  rms = xrms / real(m,fPrec)
  imax=maxmin(r(1),m,1)
  imin=maxmin(r(1),m,0)
  ptp=r(imax)-r(imin)
  rms=sqrt(rms)
  return
end subroutine calrms

!-----------------------------------------------------------------------
!     if M=0, MAXMIN=lowest index of minimum element in A
!     if M=1, MAXMIN=lowest index of maximun element in A
!     if N<1, MAXMIN=1
!-----------------------------------------------------------------------
function maxmin (a,n,m)
  use floatPrecision
  use mathlib_bouncer
  implicit none
  integer i,m,maxmin,n
  real(kind=fPrec) a,curent
  dimension a(n)
  save
  maxmin=1
  if (n.lt.1) return
  curent=a(1)
  do i=2,n
    if ((m.eq.0).and.(a(i).ge.curent)) goto 10
    if ((m.eq.1).and.(a(i).le.curent)) goto 10
    curent=a(i)
    maxmin=i
10   continue
  end do
  return
end function maxmin

! ================================================================================================ !
!  ORGANISATION OF NONLINEAR ELEMENTS AND RANDOM NUMBERS
!  reserving places in zfz:
!  - 1+2 for misalignment (h/v - DISP block / fort.8 / fort.30);
!  - 2xmmul for multipole errors (fort.16);
!  mapping also position of errors for a given element in lattice
! ================================================================================================ !
subroutine ord

  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  use crcoall
  use parpro
  use mod_common
  use mod_commons
  use mod_commont
  use mod_fluc

  implicit none
  integer i,inz,iran,ix,izu,j,jra,jra3,kpz,kzz,kzz1,kzz2,nra1
  dimension jra(nele,5),iran(nele),inz(nele)
  save
!-----------------------------------------------------------------------
  ! initialisation
  do i=1,nele
    iran(i)=0
    inz(i)=0
    do j=1,5
      jra(i,j)=0
    end do
  end do

  ! ORGANISATION OF RANDOM NUMBERS
  if(niu(1).lt.0) niu(1)=iabs(niu(1))
  if(niu(2).lt.0) niu(2)=iabs(niu(2))
  if(niu(1).eq.0) niu(1)=1
  if(niu(2).eq.0) niu(2)=iu
  if(niu(1).gt.iu) niu(1)=1
  if(niu(2).gt.iu) niu(2)=iu
  izu=0
  nra1=nran
  iorg=iorg-1
  if(iorg.ge.0) then
    if(iorg.eq.0) then !iorg == 0
      do i=1,iu
        ix=ic(i)
        ! skip blocks:
        if(ix.le.nblo) cycle
        ix=ix-nblo
        kpz=kp(ix)
        kzz=kz(ix)
        ! skip RF cavity, inactive non-linear elements, BB lenses, phase-trombones, wires
        if(kpz.eq.6.or.kzz.eq.0.or.kzz.eq.20.or.kzz.eq.22.or.kzz.eq.15) cycle
        !  map position of errors for present element in lattice structure
        mzu(i)=izu
        izu=izu+3
        if(kzz.eq.11.and.abs(ek(ix)).gt.pieni) izu=izu+2*mmul
        if(izu > nran) then
          write(lout,"(a,i0,a)") "ORD> ERROR The random number: ",nran," for the initial structure is too small."
          call prror(-1)
        end if
        if(izu > nzfz) then
          call fluc_moreRandomness
        endif
      end do
    else ! iorg.gt.0
      do i=1,iorg
        do j=1,il
          if(bez(j).eq.bezr(1,i)) then
            jra(i,1)=j
            if(kz(j) == 0 .or. kz(j) == 20 .or. kz(j) == 22) then
              write(lout,"(a)") "ORD> ERROR Elements that need random numbers have a kz not equal to 0, 20 or 22."
              call prror(-1)
            end if
            jra(i,2)=kz(j)
          endif
          if(bez(j).eq.bezr(2,i)) then
            jra(i,3)=j
            if(kz(j) == 0 .or. kz(j) == 20 .or. kz(j) == 22) then
              write(lout,"(a)") "ORD> ERROR Elements that need random numbers have a kz not equal to 0, 20 or 22."
              call prror(-1)
            end if
            jra(i,4)=kz(j)
          endif
        end do
        kzz1=jra(i,2)
        kzz2=jra(i,4)
        if(kzz1.ne.0.and.kzz2.eq.0) then
          jra(i,5)=nra1
          nra1=nra1+mran*3
          if(kzz1.eq.11.and.abs(ek(jra(i,1))).gt.pieni) nra1=nra1+mran*2*mmul
          if(nra1 > nzfz) then
            call fluc_moreRandomness
          end if
        endif
        if(kzz1 == 11 .and. (kzz2 /= 11 .and. kzz2 /= 0)) then
          write(lout,"(a)") "ORD> ERROR To use the same random numbers for 2 elements, the inserted element "//&
            "must not need more of such numbers than the reference element."
          call prror(-1)
        end if
      end do
      do i=1,iu
        ix=ic(i)
        ! skip blocks:
        if(ix.le.nblo) cycle
        ix=ix-nblo
        kpz=kp(ix)
        kzz=kz(ix)
        ! skip RF cavity, inactive non-linear elements, BB lenses, phase-trombones, wires
        if(kpz.eq.6.or.kzz.eq.0.or.kzz.eq.20.or.kzz.eq.22.or.kzz.eq.15) cycle
        do j=1,iorg
          if(bez(ix).eq.bezr(1,j)) goto 90
        end do
        goto 100
90       jra3=jra(j,3)
        if(jra3.ne.0) then
          ! map position of errors for present element in lattice structure
          mzu(i)=iran(jra3)
          iran(ix)=mzu(i)
        else
          inz(j)=inz(j)+1
          if(inz(j) > mran) then
            write(lout,"(a,i0,a)") "ORD> ERROR Not more than ",mran," of each type of inserted elements can be used."
            call prror(-1)
          end if
          ! map position of errors for present element in lattice structure
          mzu(i)=jra(j,5)
          iran(ix)=mzu(i)
          jra(j,5)=jra(j,5)+3
          if(jra(j,2).eq.11) jra(j,5)=jra(j,5)+2*mmul
        endif
        cycle
        ! map position of errors for present element in lattice structure
100       mzu(i)=izu
        iran(ix)=izu
        izu=izu+3
        if(kzz.eq.11.and.abs(ek(ix)).gt.pieni) izu=izu+2*mmul
        if(izu > nran) then
          write(lout,"(a,i0,a)") "ORD> ERROR The random number: ",nran," for the initial structure is too small."
          call prror(-1)
        end if
      end do
    endif
  else !iorg < 0 (in case of no ORGA block in fort.3)
    do i=1,iu
      ix=ic(i)
      ! skip blocks:
      if(ix.le.nblo) cycle
      ix=ix-nblo
      kpz=kp(ix)
      kzz=kz(ix)
      ! skip RF cavity, inactive non-linear elements, BB lenses, phase-trombones, wires
      if(kpz.eq.6.or.kzz.eq.0.or.kzz.eq.20.or.kzz.eq.22.or.kzz.eq.15) cycle
      izu=izu+3
      if(kzz.eq.11.and.abs(ek(ix)).gt.pieni) izu=izu+2*mmul
      ! why just checking? shouldn't we map on mzu(i)?
      if(izu > nran) then
        write(lout,"(a,i0,a)") "ORD> ERROR The random number: ",nran," for the initial structure is too small."
        call prror(-1)
      end if
      if(izu > nzfz) then
        call fluc_moreRandomness
      end if
    end do
  end if

  ! Misalignments
  izu = 0
  do i=1,iu
    ix = ic(i)
    if(ix.le.nblo) cycle ! Skip blocks
    ix  = ix-nblo
    kpz = kp(ix)
    kzz = kz(ix)
    ! Skip RF cavity, inactive non-linear elements, BB lenses, phase-trombones, wires
    if(kpz == 6 .or. kzz == 0 .or. kzz == 20 .or. kzz == 22 .or. kzz == 15) cycle
    if(icextal(i) > 0) then
      ! Use values from fort.8
      izu        = izu+3
      xrms(ix)   = one
      zrms(ix)   = one
      zfz(izu-1) = fluc_errAlign(1,icextal(i))
      zfz(izu)   = fluc_errAlign(2,icextal(i))
      tiltc(i)   = cos_mb(fluc_errAlign(3,icextal(i))*c1m3)
      tilts(i)   = sin_mb(fluc_errAlign(3,icextal(i))*c1m3)
    else if(icextal(i) < 0) then
      ! Use values from fort.30
      izu        = izu+3
      xrms(ix)   = one
      zrms(ix)   = one
      zfz(izu-2) = fluc_errZFZ(1,-icextal(i))
      zfz(izu-1) = fluc_errZFZ(2,-icextal(i))
      zfz(izu)   = fluc_errZFZ(3,-icextal(i))
      tiltc(i)   = cos_mb(fluc_errZFZ(4,-icextal(i))*c1m3)
      tilts(i)   = sin_mb(fluc_errZFZ(4,-icextal(i))*c1m3)
    else
      izu = izu+3
    end if
    if(kzz == 11 .and. abs(ek(ix)) > pieni .and. icext(i) /= 0) then
      do j=1,mmul
        izu     = izu+1
        zfz(izu)= fluc_errExt(20+j,icext(i))
        izu     = izu+1
        zfz(izu)= fluc_errExt(j,icext(i))
      end do
    else if(kzz == 11 .and. abs(ek(ix)) > pieni .and. icext(i) == 0) then
      izu = izu+2*mmul
    end if
  end do

end subroutine ord

! ================================================================================================ !
!  ORGANISATION OF LATTICE SEQUENCE
!  A.Mereghetti, last modified 2016-03-14
!  extract pieces of code about lattice structure from original ord subroutine
! ================================================================================================ !
subroutine orglat

  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  use crcoall
  use parpro
  use mod_common
  use mod_commons
  use mod_commont
  implicit none
  integer i,icext1,icextal1,ihi,ii,ilf,ilfr,j,kanf1
  real(kind=fPrec) extalig1 !,exterr1
  dimension ilf(nblz),ilfr(nblz)
  dimension extalig1(nblz,3),icext1(nblz),icextal1(nblz) !,exterr1(nblz,40)
  save
!-----------------------------------------------------------------------
  if(mper.gt.1) then
    do i=2,mper
      do j=1,mbloz
        ii=(i-1)*mbloz+j
        ihi=j
        if(msym(i).lt.0) ihi=mbloz-j+1
        ic(ii)=msym(i)*ic(ihi)
        if(ic(ii).lt.-nblo) ic(ii)=-ic(ii)
      end do
    end do
  end if

  ! set iu
  iu=mper*mbloz

  ! "GO" was not the first structure element -> Reshuffle the structure
  if(kanf.ne.1) then
    write(lout,"(a)") "ORGLAT> Reshuffling lattice structure following existence "//&
      "of GO keyword not in first position of lattice definition!"
    !        initialise some temporary variables
    do i=1,nblz
      ilf(i)=0
      ilfr(i)=0
      icext1(i)=0
      icextal1(i)=0
      do ii=1,3
          extalig1(i,ii)=zero
      enddo
    enddo

    !--Re-saving of the starting point (UMSPEICHERUNG AUF DEN STARTPUNKT)
    kanf1=kanf-1
    do i=1,kanf1
      if(iorg.ge.0) ilfr(i)=mzu(i)
      ilf(i)=ic(i)
      icext1(i)=icext(i)
      icextal1(i)=icextal(i)
    end do
    do i=kanf,iu
      if(iorg.ge.0) mzu(i-kanf1)=mzu(i)
      ic(i-kanf1)=ic(i)
      icext(i-kanf1)=icext(i)
      icextal(i-kanf1)=icextal(i)
    end do
    do i=1,kanf1
      if(iorg.ge.0) mzu(iu-kanf1+i)=ilfr(i)
      ic(iu-kanf1+i)=ilf(i)
      icext(iu-kanf1+i)=icext1(i)
      icextal(iu-kanf1+i)=icextal1(i)
    end do
  endif
  return
end subroutine orglat

!-----------------------------------------------------------------------
!  ADDITIONAL ADJUSTMENT OF THE X-PHASEADVANCE BETWEEN 2 POSITIONS
!-----------------------------------------------------------------------
subroutine phasad(dpp,qwc)
  ! Rewritten to remove computed gotos by V.K.B.Olsen on 20/11/2017
  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  use crcoall
  use parpro
  use mod_common
  use mod_commons
  use mod_commont
  implicit none
  integer i,ikpv,im,ium,ix,izu,j,jj,jk,jm,k,kpv,kpz,kzz,l,l1,ll,nmz,dj
  real(kind=fPrec) aa,alfa,bb,benkr,beta,ci,cikve,cr,crkve,crkveuk,dphi,dpp,dppi,dpr,&
                   dyy1,dyy2,ekk,phi,phibf,pie,puf,qu,qv,qw,qwc,qxsa,qxse,r0,r0a,t,xl,xs,zl,zs,quz,qvz
#ifdef TILT
  real(kind=fPrec) dyy11,qu1,tiltck,tiltsk
#endif
  dimension t(5,4)
  dimension beta(2),alfa(2),phi(2),phibf(2)
  dimension qw(2),qwc(3)
  dimension aa(mmul),bb(mmul),dpr(5)
  dimension cr(mmul),ci(mmul)
      save
!-----------------------------------------------------------------------
      ium=5
!GRD
      qxsa = zero
      qxse = zero
!GRD
      do i=1,ium
        dpr(i)=zero
      end do

      do i=1,ium
        do j=1,4
          t(i,j)=zero
        end do
      end do

      do i=1,2
        beta(i)=zero
        alfa(i)=zero
        phi(i)=zero
        phibf(i)=zero
        qw(i)=zero
        qwc(i)=zero
      end do

      qwc(3)=zero

      do i=1,mmul
        aa(i)=zero
        bb(i)=zero
        cr(i)=zero
        ci(i)=zero
      end do

      pie=two*pi
      ikpv=0
      dpr(1)=dpp*c1e3
      call clorb(dpp)
      call betalf(dpp,qw)
      if(ierro /= 0) then
        write(lout,"(a)") "PHASAD> ERROR No optical solution."
        call prror(-1)
      end if
      call envar(dpp)
#ifdef DEBUG
!     call warr('qw',qw(1),1,0,0,0)
!     call warr('qw',qw(2),2,0,0,0)
!     call warr('qwc',qwc(1),1,0,0,0)
!     call warr('qwc',qwc(2),2,0,0,0)
!     call warr('qwc',qwc(3),3,0,0,0)
!     call dumpbin('aenvarqmod',88,R88
!     call abend('aenvarqmod                                        ')
#endif

!--STARTVALUES OF THE TRAJECTORIES
      do l=1,2
        ll=2*l
        alfa(l)=alf0(l)
        beta(l)=bet0(l)
        t(1,ll-1)=clo(l)
        t(1,ll)=clop(l)
      end do

      do i=1,4
        do j=1,4
          t(i+1,j)=ta(j,i)
          t(i+1,j)=ta(j,i)
        end do
      end do

!--SINGLE TURN BLOCKLOOP
      izu=0
      do 450 k=1,iu
        ix=ic(k)
        if(ix.gt.nblo) goto 140
        jj=0
        dj=1
        if(ix.gt.0) goto 70
        ix=-ix
        jj=mel(ix)+1
        dj=-1
   70   jm=mel(ix)

!--BLOCKELEMENTLOOP
        do 130 j=1,jm
          jj=jj+dj
          jk=mtyp(ix,jj)
          if(ithick.eq.1.and.kz(jk).ne.0) goto 100
          if(ithick.eq.0.and.kz(jk).ne.0) goto 450

!--PURE DRIFTLENGTH
          do l=1,2
            ll=2*l

            if(abs(t(ll,ll-1)).gt.pieni) then
              phibf(l)=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
            else
              phibf(l)=pi2
            endif

            do i=1,ium
              t(i,ll-1)=t(i,ll-1)+t(i,ll)*(el(jk))
            end do
          end do

          do l=1,2
            ll=2*l
            beta(l)=t(ll,ll-1)**2+t(ll+1,ll-1)**2                        !hr06
            alfa(l)=-one*(t(ll,ll-1)*t(ll,ll)+t(ll+1,ll-1)*t(ll+1,ll))   !hr06

            if(abs(t(ll,ll-1)).gt.pieni) then
              dphi=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
            else
              dphi=pi2-phibf(l)
            endif

            if(-one*dphi.gt.pieni) dphi=dphi+pi                          !hr06
            phi(l)=phi(l)+dphi/pie
          end do

          goto 130
!--MAGNETELEMENT
  100     continue
          do l=1,2
            ll=2*l
            if(abs(t(ll,ll-1)).gt.pieni) then
              phibf(l)=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
            else
              phibf(l)=zero
            endif
            do i=1,ium
              puf=t(i,ll-1)
            t(i,ll-1)=(puf*a(jk,l,1)+t(i,ll)*a(jk,l,2))+dpr(i)*a(jk,l,5) !hr06
            t(i,ll)=(puf*a(jk,l,3)+t(i,ll)*a(jk,l,4))+dpr(i)*a(jk,l,6)   !hr06
            enddo
          enddo
          do l=1,2
            ll=2*l
            beta(l)=t(ll,ll-1)**2+t(ll+1,ll-1)**2                        !hr06
            alfa(l)=-one*(t(ll,ll-1)*t(ll,ll)+t(ll+1,ll-1)*t(ll+1,ll))   !hr06
            if(abs(t(ll,ll-1)).gt.pieni) then
              dphi=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
            else
              dphi=-one*phibf(l)                                         !hr06
            endif
            if(kz(jk).ne.8.and.-one*dphi.gt.pieni) dphi=dphi+pi          !hr06
            phi(l)=phi(l)+dphi/pie
          enddo
  130   continue
        goto 450
!--NL-INSERTION
  140   ix=ix-nblo
        qu=zero
        qv=zero
        dyy1=zero
        dyy2=zero
        kpz=kp(ix)
        if(kpz.eq.6) goto 450
        kzz=kz(ix)
        kpv=kpa(ix)
        if(kpv.ne.1) goto 150
        qxsa=phi(1)
  150   if(kpv.ne.2.or.ikpv.eq.1) goto 160
        qxse=phi(1)
        ikpv=1
  160   continue
        if(kzz == 22) then
          do l=1,2
            ll=2*l
            if(abs(t(ll,ll-1)).gt.pieni) then
              phibf(l)=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
            else
              phibf(l)=zero
            end if
            do i=1,ium
              puf=t(i,ll-1)
              t(i,ll-1)=(puf*rrtr(imtr(ix),ll-1,ll-1)+t(i,ll)*rrtr(imtr(ix),ll-1,ll))+dpr(i)*rrtr(imtr(ix),ll-1,6)
              t(i,ll)=(puf*rrtr(imtr(ix),ll,ll-1)+t(i,ll)*rrtr(imtr(ix),ll,ll))+dpr(i)*rrtr(imtr(ix),ll,6)
            end do
            t(1,ll-1)=t(1,ll-1)+cotr(imtr(ix),ll-1)
            t(1,ll)=t(1,ll)+cotr(imtr(ix),ll)
            beta(l)=t(ll,ll-1)**2+t(ll+1,ll-1)**2
            alfa(l)=-one*(t(ll,ll-1)*t(ll,ll)+t(ll+1,ll-1)*t(ll+1,ll))
            if(abs(t(ll,ll-1)) > pieni) then
              dphi=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
            else
              dphi=-one*phibf(l)
            end if
            if(-one*dphi.gt.pieni) dphi=dphi+pi
            phi(l)=phi(l)+dphi/pie
          end do
        end if
        if(kzz.eq.0.or.kzz.eq.20.or.kzz.eq.22) goto 450
        if(kzz.eq.15) goto 450
! JBG RF CC Multipoles to 450
!        if(kzz.eq.26.or.kzz.eq.27.or.kzz.eq.28) write(*,*)'out'
!        if(kzz.eq.26.or.kzz.eq.27.or.kzz.eq.28) goto 450
        dyy1=zero
        dyy2=zero
        if(iorg.lt.0) mzu(k)=izu
        izu=mzu(k)+1
        ekk=(sm(ix)+zfz(izu)*ek(ix))/(one+dpp)
        izu=izu+1
        xs=xpl(ix)+zfz(izu)*xrms(ix)
        izu=izu+1
        zs=zpl(ix)+zfz(izu)*zrms(ix)
#include "include/alignl.f90"
#ifdef DEBUG
!     call warr('qw',qw(1),1,0,0,0)
!     call warr('qw',qw(2),2,0,0,0)
!     call warr('qwc',qwc(1),1,0,0,0)
!     call warr('qwc',qwc(2),2,0,0,0)
!     call warr('qwc',qwc(3),3,0,0,0)
!     call warr('kzz',0d0,kzz,0,0,0)
!     call dumpbin('bkzz      ',77 777)
!     call abend('bkzz                                              ')
#endif
      select case (kzz)
      case (1) ! HORIZONTAL DIPOLE
        ekk=ekk*c1e3
#include "include/kickl01h.f90"
#include "include/kickq01h.f90"
        goto 420
      case (2) ! NORMAL QUADRUPOLE
#include "include/kicklxxh.f90"
#include "include/kickq02h.f90"
        goto 420
      case (3) ! NORMAL SEXTUPOLE
        ekk=ekk*c1m3
#include "include/kickq03h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
        goto 420
      case (4) ! NORMAL OCTUPOLE
        ekk=ekk*c1m6
#include "include/kicksho.f90"
#include "include/kickq04h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
        goto 420
      case (5) ! NORMAL DECAPOLE
        ekk=ekk*c1m9
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq05h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
        goto 420
      case (6) ! NORMAL DODECAPOLE
        ekk=ekk*c1m12
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq06h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
        goto 420
      case (7) ! NORMAL 14-POLE
        ekk=ekk*c1m15
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq07h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
        goto 420
      case (8) ! NORMAL 16-POLE
        ekk=ekk*c1m18
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq08h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
        goto 420
      case (9) ! NORMAL 18-POLE
        ekk=ekk*c1m21
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq09h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
        goto 420
      case (10) ! NORMAL 20-POLE
        ekk=ekk*c1m24
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq10h.f90"
#include "include/kicksho.f90"
#include "include/kicklxxh.f90"
        goto 420
      case (11)
        r0=ek(ix)
        if(abs(dki(ix,1)).gt.pieni) then
          if(abs(dki(ix,3)).gt.pieni) then
#include "include/multl01.f90"
            do i=2,ium
#include "include/multl02.f90"
            end do
          else
#include "include/multl03.f90"
          end if
        end if
        if(abs(dki(ix,2)).gt.pieni) then
          if(abs(dki(ix,3)).gt.pieni) then
#include "include/multl04.f90"
            do i=2,ium
#include "include/multl05.f90"
            end do
          else
#include "include/multl06.f90"
          end if
        end if
        if(abs(r0).le.pieni) goto 450
        nmz=nmu(ix)
        if(nmz.eq.0) then
          izu=izu+2*mmul
          goto 450
        end if
        im=irm(ix)
        r0a=one
        benkr=ed(ix)/(one+dpp)
        do l=1,nmz
#include "include/multl07a.f90"
        end do
        if(nmz.ge.2) then
#include "include/multl07b.f90"
          do l=3,nmz
#include "include/multl07c.f90"
          end do
        else
#include "include/multl07d.f90"
        end if
#ifdef TILT
#include "include/multl07e.f90"
#endif
        izu=izu+2*mmul-2*nmz
        goto 420
      case (12,13,14,15,16,17,18,19,20,21,22,23)
        goto 450
      case (24) ! DIPEDGE ELEMENT
#include "include/kickldpe.f90"
#include "include/kickqdpe.f90"
        goto 420
      case (25) ! Solenoid
#include "include/kicklso1.f90"
#include "include/kickqso1.f90"
        goto 420
      case (26,27,28)
        goto 450

      !-----------------
      !--SKEW ELEMENTS--
      !------------------
      case (-1)  ! VERTICAL DIPOLE
        ekk=ekk*c1e3
#include "include/kickl01v.f90"
#include "include/kickq01v.f90"
        goto 420
      case (-2)  ! SKEW QUADRUPOLE
#include "include/kicklxxv.f90"
#include "include/kickq02v.f90"
        goto 420
      case (-3)  ! SKEW SEXTUPOLE
        ekk=ekk*c1m3
#include "include/kickq03v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
        goto 420
      case (-4)  ! SKEW OCTUPOLE
        ekk=ekk*c1m6
#include "include/kicksho.f90"
#include "include/kickq04v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
        goto 420
      case (-5)  ! SKEW DECAPOLE
        ekk=ekk*c1m9
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq05v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
        goto 420
      case (-6)  ! SKEW DODECAPOLE
        ekk=ekk*c1m12
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq06v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
        goto 420
      case (-7)  ! SKEW 14-POLE
        ekk=ekk*c1m15
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq07v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
        goto 420
      case (-8)  ! SKEW 16-POLE
        ekk=ekk*c1m18
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq08v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
        goto 420
      case (-9)  ! SKEW 18-POLE
        ekk=ekk*c1m21
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq09v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
        goto 420
      case (-10) ! SKEW 20-POLE
        ekk=ekk*c1m24
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kickq10v.f90"
#include "include/kicksho.f90"
#include "include/kicklxxv.f90"
      case default
        goto 450
      end select
      goto 450

  420 continue
      t(1,2)=t(1,2)+dyy1
      t(1,4)=t(1,4)+dyy2
      do i=2,ium
        if(kzz.eq.24) then
          t(i,2)=(t(i,2)+t(i,1)*qu)-qv*t(i,3)                          !hr06
          t(i,4)=(t(i,4)-t(i,3)*quz)-qvz*t(i,1)                        !hr06
#include "include/phas1so1.f90"
#include "include/phas2so1.f90"
#include "include/phas3so1.f90"
        else
          t(i,2)=(t(i,2)+t(i,1)*qu)-qv*t(i,3)                          !hr06
          t(i,4)=(t(i,4)-t(i,3)*qu)-qv*t(i,1)                          !hr06
        end if
      end do

      do l=1,2
        ll=2*l
        alfa(l)=-one*(t(ll,ll-1)*t(ll,ll)+t(ll+1,ll-1)*t(ll+1,ll))       !hr06
      end do

  450 continue
      qwc(1)=phi(1)
      qwc(2)=phi(2)
      if(qxse.ge.qxsa) then
        qwc(3)=qxse-qxsa
      else
        qwc(3)=(phi(1)+qxse)-qxsa                                        !hr06
      endif
#ifdef DEBUG
!     call warr('qw',qw(1),1,0,0,0)
!     call warr('qw',qw(2),2,0,0,0)
!     call warr('qwc',qwc(1),1,0,0,0)
!     call warr('qwc',qwc(2),2,0,0,0)
!     call warr('qwc',qwc(3),3,0,0,0)
!     call dumpbin('aphasad',97,997)
!     call abend('aphasad                                           ')
#endif
!-----------------------------------------------------------------------
  return
end subroutine phasad

subroutine qmod0
!-----------------------------------------------------------------------
!  ADJUSTMENT OF THE Q-VALUES PLUS AN ADDITIONAL ADJUSTMENT OF A
!  X-PHASEADVANCE BETWEEN 2 POSITIONS IN THE MACHINE
!-----------------------------------------------------------------------
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      use crcoall
      use parpro
      use mod_common
      use mod_commons
      use mod_commont
      implicit none
      integer i,ierr,ii,iq1,iq2,iq3,iql,j,l,n,nite
      real(kind=fPrec) a11,a12,a13,a21,a22,a23,a31,a32,a33,aa,aa1,bb,   &
     &dpp,dq1,dq2,dq3,qwc,qx,qz,sens,sm0,sqx,sqxh,sqz
      dimension sens(3,5),aa(3,3),bb(3),qx(3),qz(3),sm0(3),qwc(3)
      dimension aa1(2,2)
      save
!-----------------------------------------------------------------------
      do i=1,3
        bb(i)=zero
        qx(i)=zero
        qz(i)=zero
        sm0(i)=zero
        qwc(i)=zero
        do j=1,3
          aa(i,j)=zero
        end do
      end do

      do i=1,3
        do j=1,5
          sens(i,j)=zero
        end do
      end do

      do i=1,2
        do j=1,2
          aa1(i,j)=zero
        end do
      end do

      write(lout,10010)
      sqx=zero
      sqz=zero
      sqxh=zero
      dpp=zero
      iq1=iq(1)
      iq2=iq(2)
      if(kz(iq1).ne.2.or.kz(iq2).ne.2) then
        write(lout,"(a)") "QMOD> ERROR Element is not a quadrupole."
        call prror(-1)
      end if

      if (abs(el(iq1)).le.pieni.or.abs(el(iq2)).le.pieni) then
        sm0(1)=ed(iq1)
        sm0(2)=ed(iq2)
      else
        sm0(1)=ek(iq1)
        sm0(2)=ek(iq2)
      endif

      if(kp(iq1).eq.5) call combel(iq1)
      if(kp(iq2).eq.5) call combel(iq2)
      sens(1,1)=qw0(1)
      sens(2,1)=qw0(2)

      if(abs(qw0(3)).gt.pieni) then
        iq3=iq(3)
        if(kz(iq3).ne.2) then
          write(lout,"(a)") "QMOD> ERROR Element is not a quadrupole."
          call prror(-1)
        end if
        if (abs(el(iq3)).le.pieni) then
          sm0(3)=ed(iq3)
        else
          sm0(3)=ek(iq3)
        endif
        if(kp(iq3).eq.5) call combel(iq3)
        nite=3
      else
        nite=2
      endif

      call clorb(dpp)
      if(ierro.gt.0) then
        write(lout,"(a)") "QMOD> ERROR Unstable closed orbit during tune variation."
        call prror(-1)
      end if
      call phasad(dpp,qwc)
      sens(1,5)=qwc(1)
      sens(2,5)=qwc(2)
      if(nite.eq.3) then
        sens(3,1)=qw0(3)
        sens(3,5)=qwc(3)
        write(lout,10100)
        write(lout,10120) qwc,qw0
      else
        write(lout,10110)
        write(lout,10130) qwc(1),qwc(2),qw0(1),qw0(2)
      endif
      do 60 ii=1,itqv
        do 40 n=1,nite
          iql=iq(n)
          if (abs(el(iql)).le.pieni) then
            ed(iql)=ed(iql)+dkq
          else
            ek(iql)=ek(iql)+dkq
          endif
          if(kp(iql).eq.5) call combel(iql)
          call clorb(dpp)
          if(ierro.gt.0) then
            write(lout,"(a)") "QMOD> ERROR Unstable closed orbit during tune variation."
            call prror(-1)
          end if
          call phasad(dpp,qwc)
          sens(1,n+1)=qwc(1)
          sens(2,n+1)=qwc(2)
          if(nite.eq.3) then
            sens(3,n+1)=qwc(3)
            write(lout,10140) ii,n,qwc
          else
            write(lout,10150) ii,n,qwc(1),qwc(2)
          endif
          if (abs(el(iql)).le.pieni) then
            ed(iql)=ed(iql)-dkq
          else
            ek(iql)=ek(iql)-dkq
          endif
          if(kp(iql).eq.5) call combel(iql)
   40   continue
!--Q-VALUE ADJUSTMENT
        aa1(1,1)=(sens(1,2)-sens(1,5))/dkq
        aa1(1,2)=(sens(2,2)-sens(2,5))/dkq
        aa1(2,1)=(sens(1,3)-sens(1,5))/dkq
        aa1(2,2)=(sens(2,3)-sens(2,5))/dkq
        a11=aa1(1,1)
        a12=aa1(1,2)
        a21=aa1(2,1)
        a22=aa1(2,2)
        bb(1)=sens(1,5)-sens(1,1)
        bb(2)=sens(2,5)-sens(2,1)
        sqx=sqx+abs(bb(1))
        sqz=sqz+abs(bb(2))
        if(nite.eq.3) then
          aa(1,1)=a11
          aa(1,2)=a12
          aa(1,3)=(sens(3,2)-sens(3,5))/dkq
          aa(2,1)=a21
          aa(2,2)=a22
          aa(2,3)=(sens(3,3)-sens(3,5))/dkq
          aa(3,1)=(sens(1,4)-sens(1,5))/dkq
          aa(3,2)=(sens(2,4)-sens(2,5))/dkq
          aa(3,3)=(sens(3,4)-sens(3,5))/dkq
          a13=aa(1,3)
          a23=aa(2,3)
          a31=aa(3,1)
          a32=aa(3,2)
          a33=aa(3,3)
          bb(3)=sens(3,5)-sens(3,1)
          sqxh=sqxh+abs(bb(3))
          call loesd(aa,bb,nite,nite,ierr)
        else
          call loesd(aa1,bb,nite,nite,ierr)
        endif
        if(ierr == 1) then
          write(lout,"(a)") "QMOD> ERROR Problems during matrix-inversion."
          call prror(-1)
        end if
        do 50 l=1,nite
          iql=iq(l)
          if (abs(el(iql)).le.pieni) then
            ed(iql)=ed(iql)-bb(l)
          else
            ek(iql)=ek(iql)-bb(l)
          endif
          if(kp(iql).eq.5) call combel(iql)
   50   continue
        call clorb(dpp)
        if(ierro.gt.0) then
          write(lout,"(a)") "QMOD> ERROR Unstable closed orbit during tune variation."
          call prror(-1)
        end if
        call phasad(dpp,qwc)
        sens(1,5)=qwc(1)
        sens(2,5)=qwc(2)
        if(nite.eq.3) then
          sens(3,5)=qwc(3)
          write(lout,10020) qw0(1),qwc(1),qw0(2),qwc(2),qw0(3),qwc(3)
          if (abs(el(iq1)).le.pieni) then
            write(lout,10040) sm0(1),ed(iq1),bez(iq1),sm0(2),ed(iq2),bez&
     &(iq2),sm0(3),ed(iq3),bez(iq3)
          else
            write(lout,10040) sm0(1),ek(iq1),bez(iq1),sm0(2),ek(iq2),bez&
     &(iq2),sm0(3),ek(iq3),bez(iq3)
          endif
          write(lout,10080) sqx,sqz,sqxh
          write(lout,10060) a11,a12,a13,a21,a22,a23,a31,a32,a33
        else
          write(lout,10030) qw0(1),qwc(1),qw0(2),qwc(2)
          if (abs(el(iq1)).le.pieni) then
            write(lout,10050) sm0(1),ed(iq1),bez(iq1),sm0(2),ed(iq2),bez&
     &(iq2)
          else
            write(lout,10050) sm0(1),ek(iq1),bez(iq1),sm0(2),ek(iq2),bez&
     &(iq2)
          endif
          write(lout,10090) sqx,sqz
          write(lout,10070) a11,a12,a21,a22
        endif
        if (abs(el(iq(1))).le.pieni) then
          sm0(1)=ed(iq(1))
          sm0(2)=ed(iq(2))
        else
          sm0(1)=ek(iq(1))
          sm0(2)=ek(iq(2))
        endif
        dq1=abs(qwc(1)-qw0(1))
        dq2=abs(qwc(2)-qw0(2))
        if(nite.eq.3) then
          if (abs(el(iq(3))).le.pieni) then
            sm0(3)=ed(iq(3))
          else
            sm0(3)=ek(iq(3))
          endif
          dq3=abs(qwc(3)-qw0(3))
          if(dq1.lt.dqq.and.dq2.lt.dqq.and.dq3.lt.dqq) return
        else
          if(dq1.lt.dqq.and.dq2.lt.dqq) return
        endif
   60 continue
      write(lout,10000) itqv
!-----------------------------------------------------------------------
      return
10000 format(t5/t10,'TUNE ADJUSTMENT'/ t10,                             &
     &'MAXIMUM NUMBER OF ITERATIONS ACHIEVED--->',2x,i4/ t10,           &
     &'PROCEDURE MAY NOT HAVE CONVERGED')
10010 format(/131('-'))
10020 format(//131('-')//t10,'DATA BLOCK TUNE-VARIATION' / /t10,        &
     &'TUNE'           ,26x,'THEORET.     AFTER CORRECTION'/ t10,       &
     &'HORIZONTAL'     ,17x,g17.10,2x,g17.10/ t10,                      &
     &'VERTICAL'       ,19x,g17.10,2x,g17.10/ t10,                      &
     &'PART-HORIZONTAL',12x,g17.10,2x,g17.10/)
10030 format(//131('-')//t10,'DATA BLOCK TUNE-VARIATION' / /t10,        &
     &'TUNE'           ,26x,'THEORET.      AFTER CORRECTION'/ t10,      &
     &'HORIZONTAL'     ,17x,g17.10,2x,g17.10/ t10,                      &
     &'VERTICAL'       ,19x,g17.10,2x,g17.10/)
10060 format(t10,'QUADRUPOLE SENSITIVITIES',6x,'D-QX',14x,'D-QY',14x,   &
     &'D-QXH'/29x,'QF   ',d15.8,3x,d15.8,3x,d15.8/29x,                  &
     &'QD   ',d15.8,3x,d15.8,3x,d15.8/29x,                              &
     &'QF2  ',d15.8,3x,d15.8,3x,d15.8//131('-')//)
10070 format(t10,'QUADRUPOLE SENSITIVITIES',6x,'D-QX',14x,'D-QY', /29x, &
     &'QF   ',d15.8,3x,d15.8/29x,'QD   ',d15.8,3x,d15.8 //131('-')//)
10080 format(t10,'TOTAL TUNE SHIFT',10x,'QX =',f10.7,'    QY =',f10.7,  &
     &'   QXH =',f10.7)
10090 format(t10,'TOTAL TUNE SHIFT',10x,'QX =',f10.7,'    QY =',f10.7)
10100 format(t5,'---- QMOD FOR SPLIT-Q-VALUES ENTRY ---- ',             &
     &'(ZERO MOMENTUM-DEVIATION)')
10110 format(t5,'---- QMOD ENTRY ---- (ZERO MOMENTUM-DEVIATION)')
10120 format(t10,'START-QX-QY-QXH',3f12.7,' END-QX-QY-QXH',3f12.7)
10130 format(t10,'START-QX-QY',2f12.7,' END-QX-QY',2f12.7)
10140 format(t10,'ITER=',i3,'/QUAD=',i3,'/QX-QY-QXH',3f12.7)
10150 format(t10,'ITER=',i3,'/QUAD=',i3,'/QX-QY',2f12.7)
10040 format(t10,'QUADRU.STRENGTHS',7x,g17.10,2x,g17.10,'   TYP     ',  &
     &a16/t10,                  23x,g17.10,2x,g17.10,'           ',     &
     &a16)
10050 format(t10,'QUADRU.STRENGTHS',7x,g17.10,2x,g17.10,'   TYP     ',  &
     &a16/t10,                  23x,g17.10,2x,g17.10,'           ',     &
     &a16)
end subroutine qmod0

      subroutine qmodda(mm,qwc)
!-----------------------------------------------------------------------
!  ADJUSTMENT OF THE Q-VALUES VIA DA
!-----------------------------------------------------------------------
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      use crcoall
      use parpro
      use mod_common
      use mod_commons
      use mod_commont
      use mod_commond
      implicit none
      integer i,intwq,ix,mm,ncorr,ncorruo,ncrr,nd,nd2,ndh
      real(kind=fPrec) cor,coro,dq1,dq2,dps0,edcor1,edcor2,qwc
      dimension intwq(3),qwc(3)
      save
!-----------------------------------------------------------------------
#ifdef DEBUG
!     call warr('qwc',qwc(1),1,0,0,0)
!     call warr('qwc',qwc(2),2,0,0,0)
!     call warr('qwc',qwc(3),3,0,0,0)
#endif
      ncorruo=ncorru
      ncorru=1
      nd2=2*mm
      ndh=nd2+2
      intwq(1)=int(qwc(1))
      intwq(2)=int(qwc(2))
      intwq(3)=0
#ifdef DEBUG
!     call warr('intwq',0d0,intwq(1),1,0,0)
!     call warr('intwq',0d0,intwq(2),2,0,0)
!     call warr('intwq',0d0,intwq(3),3,0,0)
!     call warr('clo6(1)',clo6(1),1,0,0,0)
!     call warr('clo6(2)',clo6(2),2,0,0,0)
!     call warr('clo6(3)',clo6(3),3,0,0,0)
!     call warr('clop6(1)',clop6(1),1,0,0,0)
!     call warr('clop6(2)',clop6(2),2,0,0,0)
!     call warr('clop6(3)',clop6(3),3,0,0,0)
!sqmodda
!     write(*,*) 'qmodda called!'
!     call dumpbin('sqmodda',80,800)
!     call abend('sqmodda                                           ')
!     write(*,*) 'mm=',mm
#endif
      dq1=zero
      dq2=zero
      if(iqmod6.eq.1) then
        if(el(iq(1)).le.pieni) then
          edcor(1)=ed(iq(1))
        else
          edcor(1)=ek(iq(1))
        endif
        if(el(iq(2)).le.pieni) then
          edcor(2)=ed(iq(2))
        else
          edcor(2)=ek(iq(2))
        endif
        edcor1=edcor(1)
        edcor2=edcor(2)
#ifdef DEBUG
!       call warr('edcor1',edcor1,1,0,0,0)
!       call warr('edcor2',edcor2,2,0,0,0)
#endif
        cor=zero
        coro=1.0e38_fPrec
      endif
      do ncorr=1,itqv+1
        if(nbeam.ge.1) then
          nd=mm
#include "include/beamcou.f90"
        endif
        if(iqmod6.eq.1) write(lout,10080) nd2
        if(iqmod6.ne.1) write(lout,10090) nd2
        if(mm.eq.2) then
          write(lout,10010) clo(1),clop(1)
          write(lout,10010) clo(2),clop(2)
        elseif(mm.eq.3) then
          write(lout,10010) clo6(1),clop6(1)
          write(lout,10010) clo6(2),clop6(2)
          write(lout,10010) clo6(3),clop6(3)
        endif
        iqmodc=2
#ifdef DEBUG
!     call warr('qwc',qwc(1),1,0,0,0)
!     call warr('qwc',qwc(2),2,0,0,0)
!     call warr('qwc',qwc(3),3,0,0,0)
!     call warr('intwq',0d0,intwq(1),1,0,0)
!     call warr('intwq',0d0,intwq(2),2,0,0)
!     call warr('intwq',0d0,intwq(3),3,0,0)
!     call warr('clo6(1)',clo6(1),1,0,0,0)
!     call warr('clo6(2)',clo6(2),2,0,0,0)
!     call warr('clo6(3)',clo6(3),3,0,0,0)
!     call warr('clop6(1)',clop6(1),1,0,0,0)
!     call warr('clop6(2)',clop6(2),2,0,0,0)
!     call warr('clop6(3)',clop6(3),3,0,0,0)
!     call dumpbin('bdaini',96,996)
!     call abend('before daini                                      ')
#endif
        call mydaini(1,1,nd2,mm,nd2,1)
#ifdef DEBUG
!     call warr('qwc',qwc(1),1,0,0,0)
!     call warr('qwc',qwc(2),2,0,0,0)
!     call warr('qwc',qwc(3),3,0,0,0)
!     call warr('intwq',0d0,intwq(1),1,0,0)
!     call warr('intwq',0d0,intwq(2),2,0,0)
!     call warr('intwq',0d0,intwq(3),3,0,0)
!     call warr('clo6(1)',clo6(1),1,0,0,0)
!     call warr('clo6(2)',clo6(2),2,0,0,0)
!     call warr('clo6(3)',clo6(3),3,0,0,0)
!     call warr('clop6(1)',clop6(1),1,0,0,0)
!     call warr('clop6(2)',clop6(2),2,0,0,0)
!     call warr('clop6(3)',clop6(3),3,0,0,0)
!     call dumpbin('adaini',96,996)
!     call abend('after  daini                                      ')
#endif
        if(iqmod6.eq.1) then
          write(lout,10000) nd2
          iqmodc=1
          call mydaini(2,3,ndh,mm,nd2,1)
#ifdef DEBUG
!     call warr('qwc',qwc(1),1,0,0,0)
!     call warr('qwc',qwc(2),2,0,0,0)
!     call warr('qwc',qwc(3),3,0,0,0)
!     call warr('intwq',0d0,intwq(1),1,0,0)
!     call warr('intwq',0d0,intwq(2),2,0,0)
!     call warr('intwq',0d0,intwq(3),3,0,0)
!     call dumpbin('adaini',99,999)
!     call abend('after  daini                                      ')
#endif
          do i=1,mm
            qwc(i)=real(intwq(i),fPrec)+corr(1,i)                              !hr06
          enddo
          dq1=qwc(1)-qw0(1)
          dq2=qwc(2)-qw0(2)
          if(ncorr.eq.1) cor=sqrt(dq1**2+dq2**2)                         !hr06
          if(abs(dq1).gt.dqq.or.abs(dq2).gt.dqq) then
            cor=sqrt(dq1**2+dq2**2)                                      !hr06
            if(ncorr.eq.1.or.cor.lt.coro) then
              coro=cor
              if(el(iq(1)).le.pieni) then
                ed(iq(1))=(ed(iq(1))-corr(2,1)*dq1)-corr(2,2)*dq2        !hr06
              else
                ek(iq(1))=(ek(iq(1))-corr(2,1)*dq1)-corr(2,2)*dq2        !hr06
              endif
              if(el(iq(2)).le.pieni) then
                ed(iq(2))=(ed(iq(2))-corr(3,1)*dq1)-corr(3,2)*dq2        !hr06
              else
                ek(iq(2))=(ek(iq(2))-corr(3,1)*dq1)-corr(3,2)*dq2        !hr06
              endif
              do ncrr=1,iu
                ix=ic(ncrr)
                if(ix.gt.nblo) then
                  ix=ix-nblo
                  if(ix.eq.iq(1).or.iratioe(ix).eq.iq(1)) then
                    smi(ncrr)=ed(iq(1))*ratioe(ix)+smizf(ncrr)
                  else if(ix.eq.iq(2).or.iratioe(ix).eq.iq(2)) then
                    smi(ncrr)=ed(iq(2))*ratioe(ix)+smizf(ncrr)
                  endif
                endif
              enddo
              if(el(iq(1)).le.pieni) then
                edcor(1)=ed(iq(1))
              else
                edcor(1)=ek(iq(1))
              endif
              if(el(iq(2)).le.pieni) then
                edcor(2)=ed(iq(2))
              else
                edcor(2)=ek(iq(2))
              endif
              if(ncorr.eq.1) then
                write(lout,10020) nd2,qw0(1),qwc(1),qw0(2),qwc(2),      &
     &ncorr-1,                                                          &
     &cor
              else
                write(lout,10030) nd2,qw0(1),qwc(1),qw0(2),qwc(2),      &
     &ncorr-1,                                                          &
     &cor
              endif
              if(el(iq(1)).le.pieni.and.el(iq(2)).le.pieni) then
                write(lout,10040) edcor1,ed(iq(1)),bez(iq(1)),edcor2,   &
     &ed(iq(2)),bez(iq(2))
              elseif(el(iq(1)).le.pieni.and.el(iq(2)).gt.pieni) then
                write(lout,10040) edcor1,ed(iq(1)),bez(iq(1)),edcor2,   &
     &ek(iq(2)),bez(iq(2))
              elseif(el(iq(1)).gt.pieni.and.el(iq(2)).le.pieni) then
                write(lout,10040) edcor1,ek(iq(1)),bez(iq(1)),edcor2,   &
     &ed(iq(2)),bez(iq(2))
              else
                write(lout,10040) edcor1,ek(iq(1)),bez(iq(1)),edcor2,   &
     &ek(iq(2)),bez(iq(2))
              endif
            else
              write(lout,10050) nd2,ncorr-1
              goto 1
            endif
          else
            write(lout,10060) nd2,ncorr-1
            goto 1
          endif
        else
          iqmodc=3
          call mydaini(2,2,nd2,mm,nd2,1)
          do i=1,mm
            qwc(i)=real(intwq(i),fPrec)+wxys(i)                                !hr06
          enddo
#ifdef DEBUG
!     call warr('qwc',qwc(1),1,0,0,0)
!     call warr('qwc',qwc(2),2,0,0,0)
!     call warr('qwc',qwc(3),3,0,0,0)
!     call warr('intwq',0d0,intwq(1),1,0,0)
!     call warr('intwq',0d0,intwq(2),2,0,0)
!     call warr('intwq',0d0,intwq(3),3,0,0)
!     call dumpbin('adaini',98,998)
!     call abend('after  daini 98                                   ')
#endif
          goto 1
        endif
      enddo
 1    continue
      if(iqmod6.eq.1) then
        do ncrr=1,iu
          ix=ic(ncrr)
          if(ix.le.nblo) then
            if(iratioe(ix).eq.iq(1)) ek(ix)=ek(iq(1))*ratioe(ix)
            if(iratioe(ix).eq.iq(2)) ek(ix)=ek(iq(2))*ratioe(ix)
          endif
        enddo
        iqmodc=3
        call mydaini(2,2,nd2,mm,nd2,1)
        do i=1,mm
          qwc(i)=real(intwq(i),fPrec)+wxys(i)                                  !hr06
        enddo
        if(ncorr.eq.itqv+1) write(lout,10070) nd2,itqv
        if(ncorr.eq.1) then
          write(lout,10020) nd2,qw0(1),qwc(1),qw0(2),qwc(2),ncorr-1,cor
        else
          write(lout,10030) nd2,qw0(1),qwc(1),qw0(2),qwc(2),ncorr-1,cor
        endif
        if(el(iq(1)).le.pieni.and.el(iq(2)).le.pieni) then
          write(lout,10040)edcor1,ed(iq(1)),bez(iq(1)),edcor2,ed(iq(2)),&
     &bez(iq(2))
        elseif(el(iq(1)).le.pieni.and.el(iq(2)).gt.pieni) then
          write(lout,10040)edcor1,ed(iq(1)),bez(iq(1)),edcor2,ek(iq(2)),&
     &bez(iq(2))
        elseif(el(iq(1)).gt.pieni.and.el(iq(2)).le.pieni) then
          write(lout,10040)edcor1,ek(iq(1)),bez(iq(1)),edcor2,ed(iq(2)),&
     &bez(iq(2))
        else
          write(lout,10040)edcor1,ek(iq(1)),bez(iq(1)),edcor2,ek(iq(2)),&
     &bez(iq(2))
        endif
      endif
      ncorru=ncorruo
#ifdef DEBUG
!     call dumpbin('end qmodda',7,999)
!     call abend('end qmodda 7 999                                  ')
#endif
!-----------------------------------------------------------------------
10000 format(/131('-')/t10,'ENTERING ',i1,'D DA TUNE-VARIATION')
10010 format(1x,f47.33/1x,f47.33)
10020 format(/131('-')/t10,i1,'D DA TUNE-VARIATION'/t10,                &
     &'TUNE'           ,26x,'THEORET.       BEFORE CORRECTION'/ t10,    &
     &'HORIZONTAL'     ,15x,G21.14,G21.14/ t10,                         &
     &'VERTICAL'       ,17x,G21.14,G21.14// t10,                        &
     &'ITERATION:'     ,21x,i3/ t10,                                    &
     &'ACCURACY:'      ,17x,g17.10/)
10030 format(/131('-')/t10,i1,'D DA TUNE-VARIATION'/t10,                &
     &'TUNE'           ,26x,'THEORET.       AFTER CORRECTION'/ t10,     &
     &'HORIZONTAL'     ,15x,G21.14,G21.14/ t10,                         &
     &'VERTICAL'       ,17x,G21.14,G21.14// t10,                        &
     &'ITERATION:'     ,21x,i3/ t10,                                    &
     &'ACCURACY:'      ,17x,g17.10/)
10040 format(t10,'QUADRUPOLE STRENGTH',6x,g17.10,4x,g17.10,'   TYP     '&
     &,a16/t10,                  25x,g17.10,4x,g17.10,'           ',    &
     &a16)
10050 format(/t5,'---- NO IMPROVEMENT OF ',i1,'D DA TUNE-VARIATION ',   &
     &'IN ITERATION: ',i4/)
10060 format(/t10,i1,'D DA TUNE-VARIATION SUCCESSFUL IN ITERATION: ',   &
     &i4/)
10070 format(/t10,i1,'D DA TUNE-VARIATION'/ t10,                        &
     &'MAXIMUM NUMBER OF ITERATIONS ACHIEVED--->',2x,i4/ t10,           &
     &'PROCEDURE MAY NOT HAVE CONVERGED')
10080 format(/t10,'Initial ',i1,'-D DA CLOSED ORBIT IN QMODDA')
10090 format(/t10,'Initial ',i1,'-D DA CLOSED ORBIT IN QMODDA (NO ',    &
     &'TUNE ADJUSTEMENT)')
end

!-----------------------------------------------------------------------
!     ONE TURN-TRANSFORMATION (INCLUDING QUADRUPOLE CONTRIBUTIONS)
!-----------------------------------------------------------------------
subroutine umlauf(dpp,ium,ierr)
  ! Rewritten to remove computed gotos by V.K.B.Olsen on 23/11/2017
  use floatPrecision
  use mathlib_bouncer
  use numerical_constants
  use parpro
  use crcoall
  use mod_common
  use mod_commons
  use mod_commont
  implicit none
  integer i,ierr,im,ium,ix,izu,j,k,kpz,kx,kzz,l,ll,l1,nmz
  real(kind=fPrec) aa,bb,benkr,ci,cikve,cr,crkve,crkveuk,dpp,dpr,dyy1,dyy2,ekk,puf,qu,qv,quz,qvz,r0,r0a,xl,xs,zl,zs
#ifdef TILT
  real(kind=fPrec) dyy11,qu1,tiltck,tiltsk
#endif
  dimension aa(mmul),bb(mmul),dpr(5)
  dimension cr(mmul),ci(mmul)
  save

  do i=1,mmul
    aa(i)=zero
    bb(i)=zero
    cr(i)=zero
    ci(i)=zero
  end do
  do i=1,5
    dpr(i)=zero
  end do
  ierr=0
  dpr(1)=dpp*c1e3
  izu=0
  do 350 k=1,iu
    ix=ic(k)
    if(ix.gt.nblo) goto 60
    if(ix.le.0) goto 40

    do j=1,ium
      do kx=1,2
        if(ithick.eq.1) then
          puf=x(j,kx)
          x(j,kx)=(bl1(ix,kx,1)*puf+bl1(ix,kx,2)*y(j,kx))+dpr(j)*bl1(ix,kx,5) !hr06
          y(j,kx)=(bl1(ix,kx,3)*puf+bl1(ix,kx,4)*y(j,kx))+dpr(j)*bl1(ix,kx,6) !hr06
        else
          x(j,kx)=x(j,kx)+bl1(ix,kx,2)*y(j,kx)
        end if
      end do
    end do
    goto 350

40  ix=-ix
    do j=1,ium
      do kx=1,2
        if(ithick.eq.1) then
          puf=x(j,kx)
          x(j,kx)=(bl2(ix,kx,1)*puf+bl2(ix,kx,2)*y(j,kx))+dpr(j)*bl2(ix,kx,5) !hr06
          y(j,kx)=(bl2(ix,kx,3)*puf+bl2(ix,kx,4)*y(j,kx))+dpr(j)*bl2(ix,kx,6) !hr06
        else
          x(j,kx)=x(j,kx)+bl2(ix,kx,2)*y(j,kx)
        end if
      end do
    end do

    goto 350

60  ix=ix-nblo
    qu=zero
    qv=zero
    dyy1=zero
    dyy2=zero
    kpz=kp(ix)
    if(kpz.eq.6) goto 350
    kzz=kz(ix)
    if(abs(x(1,1)).lt.aper(1).and.abs(x(1,2)).lt.aper(2)) goto 70
    ierr=1
    write(lout,"(a)") "UMLAUF> Error amplitudes exceed the maximum values."
    call prror(-1)
    return

70  continue
    if(kzz.eq.22) then
      do j=1,ium
        do kx=1,2
          ll=kx*2
          puf=x(j,kx)
          x(j,kx)=((cotr(imtr(ix),ll-1)+rrtr(imtr(ix),ll-1,ll-1)*puf)+rrtr(imtr(ix),ll-1,ll)*y(j,kx))+dpr(j)*rrtr(imtr(ix),ll-1,6)
          y(j,kx)=((cotr(imtr(ix),ll)+rrtr(imtr(ix),ll,ll-1)*puf)+rrtr(imtr(ix),ll,ll)*y(j,kx))+dpr(j)*rrtr(imtr(ix),ll,6)
        end do
      end do
    end if
    if(kzz.eq.0.or.kzz.eq.20.or.kzz.eq.22) goto 350
    if(kzz.eq.15) goto 350
! JBG RF CC Multipoles to 350
!        if(kzz.eq.26.or.kzz.eq.27.or.kzz.eq.28) write(*,*)'out'
!        if(kzz.eq.26.or.kzz.eq.27.or.kzz.eq.28) goto 350
    if(iorg.lt.0) mzu(k)=izu
    izu=mzu(k)+1
    ekk=(sm(ix)+zfz(izu)*ek(ix))/(one+dpp)
    izu=izu+1
    xs=xpl(ix)+zfz(izu)*xrms(ix)
    izu=izu+1
    zs=zpl(ix)+zfz(izu)*zrms(ix)
#include "include/alignu.f90"

    select case (kzz)
    case (1) ! HORIZONTAL DIPOLE
      ekk=ekk*c1e3
#include "include/kicku01h.f90"
      goto 350
    case (2) ! NORMAL QUADRUPOLE
#include "include/kickuxxh.f90"
      if(ium.eq.1) goto 350
#include "include/kickq02h.f90"
      goto 330
    case (3) ! NORMAL SEXTUPOLE
      ekk=ekk*c1m3
      if(ium.ne.1) then
#include "include/kickq03h.f90"
      end if
#include "include/kicksho.f90"
#include "include/kickuxxh.f90"
      if(ium.eq.1) goto 350
      goto 330
    case (4) ! NORMAL OCTUPOLE
      ekk=ekk*c1m6
#include "include/kicksho.f90"
      if(ium.ne.1) then
#include "include/kickq04h.f90"
      end if
#include "include/kicksho.f90"
#include "include/kickuxxh.f90"
      if(ium.eq.1) goto 350
      goto 330
    case (5) ! NORMAL DECAPOLE
      ekk=ekk*c1m9
#include "include/kicksho.f90"
#include "include/kicksho.f90"
      if(ium.ne.1) then
#include "include/kickq05h.f90"
      end if
#include "include/kicksho.f90"
#include "include/kickuxxh.f90"
      if(ium.eq.1) goto 350
      goto 330
    case (6) ! NORMAL DODECAPOLE
      ekk=ekk*c1m12
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
      if(ium.ne.1) then
#include "include/kickq06h.f90"
      end if
#include "include/kicksho.f90"
#include "include/kickuxxh.f90"
      if(ium.eq.1) goto 350
      goto 330
    case (7) ! NORMAL 14-POLE
      ekk=ekk*c1m15
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
      if(ium.ne.1) then
#include "include/kickq07h.f90"
      end if
#include "include/kicksho.f90"
#include "include/kickuxxh.f90"
      if(ium.eq.1) goto 350
      goto 330
    case (8) ! NORMAL 16-POLE
      ekk=ekk*c1m18
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
      if(ium.ne.1) then
#include "include/kickq08h.f90"
      end if
#include "include/kicksho.f90"
#include "include/kickuxxh.f90"
      if(ium.eq.1) goto 350
      goto 330
    case (9) ! NORMAL 18-POLE
      ekk=ekk*c1m21
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
      if(ium.ne.1) then
#include "include/kickq09h.f90"
      end if
#include "include/kicksho.f90"
#include "include/kickuxxh.f90"
      if(ium.eq.1) goto 350
      goto 330
    case (10) ! NORMAL 20-POLE
      ekk=ekk*c1m24
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
      if(ium.ne.1) then
#include "include/kickq10h.f90"
      end if
#include "include/kicksho.f90"
#include "include/kickuxxh.f90"
      if(ium.eq.1) goto 350
      goto 330
    case (11)
      r0=ek(ix)
      if(abs(dki(ix,1)).gt.pieni) then
        if(abs(dki(ix,3)).gt.pieni) then
#include "include/multu01.f90"
          do j=2,ium
#include "include/multu02.f90"
          end do
        else
#include "include/multu03.f90"
        end if
      end if
      if(abs(dki(ix,2)).gt.pieni) then
        if(abs(dki(ix,3)).gt.pieni) then
#include "include/multu04.f90"
          do j=2,ium
#include "include/multu05.f90"
          end do
        else
#include "include/multu06.f90"
        end if
      end if
      if(abs(r0).le.pieni) goto 350
      nmz=nmu(ix)
      if(nmz.eq.0) then
        izu=izu+2*mmul
        goto 350
      end if
      im=irm(ix)
      r0a=one
      benkr=ed(ix)/(one+dpp)
      do l=1,nmz
#include "include/multl07a.f90"
      end do
      if(nmz.ge.2) then
#include "include/multl07b.f90"
        do l=3,nmz
#include "include/multl07c.f90"
        end do
      else
#include "include/multl07d.f90"
      end if
#ifdef TILT
#include "include/multl07e.f90"
#endif
      izu=izu+2*mmul-2*nmz
      y(1,1)=y(1,1)+dyy1
      y(1,2)=y(1,2)+dyy2
      if(ium.eq.1) goto 350
      goto 330
    case (12,13,14,15,16,17,18,19,20,21,22,23)
      goto 350
    case (24) ! DIPEDGE ELEMENT
#include "include/kickudpe.f90"
      if(ium.eq.1) goto 350
#include "include/kickqdpe.f90"
      goto 330
    case (25) ! Solenoid
#include "include/kickuso1.f90"
      if(ium.eq.1) goto 350
#include "include/kickqso1.f90"
      goto 330
    case (26,27,28)
      goto 350

    !-----------------
    !--SKEW ELEMENTS--
    !------------------
    case (-1) ! VERTICAL DIPOLE
      ekk=ekk*c1e3
#include "include/kicku01v.f90"
      goto 350
    case (-2) ! SKEW QUADRUPOLE
#include "include/kickuxxv.f90"
      if(ium.eq.1) goto 350
#include "include/kickq02v.f90"
      goto 330
    case (-3) ! SKEW SEXTUPOLE
      ekk=ekk*c1m3
      if(ium.ne.1) then
#include "include/kickq03v.f90"
      end if
#include "include/kicksho.f90"
#include "include/kickuxxv.f90"
      if(ium.eq.1) goto 350
      goto 330
    case (-4) ! SKEW OCTUPOLE
      ekk=ekk*c1m6
#include "include/kicksho.f90"
      if(ium.ne.1) then
#include "include/kickq04v.f90"
      end if
#include "include/kicksho.f90"
#include "include/kickuxxv.f90"
      if(ium.eq.1) goto 350
      goto 330
    case (-5) ! SKEW DECAPOLE
      ekk=ekk*c1m9
#include "include/kicksho.f90"
#include "include/kicksho.f90"
      if(ium.ne.1) then
#include "include/kickq05v.f90"
      end if
#include "include/kicksho.f90"
#include "include/kickuxxv.f90"
      if(ium.eq.1) goto 350
      goto 330
    case (-6) ! SKEW DODECAPOLE
      ekk=ekk*c1m12
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
      if(ium.ne.1) then
#include "include/kickq06v.f90"
      end if
#include "include/kicksho.f90"
#include "include/kickuxxv.f90"
      if(ium.eq.1) goto 350
      goto 330
    case (-7) ! SKEW 14-POLE
      ekk=ekk*c1m15
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
      if(ium.ne.1) then
#include "include/kickq07v.f90"
      end if
#include "include/kicksho.f90"
#include "include/kickuxxv.f90"
      if(ium.eq.1) goto 350
      goto 330
    case (-8) ! SKEW 16-POLE
      ekk=ekk*c1m18
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
      if(ium.ne.1) then
#include "include/kickq08v.f90"
      end if
#include "include/kicksho.f90"
#include "include/kickuxxv.f90"
      if(ium.eq.1) goto 350
      goto 330
    case (-9) ! SKEW 18-POLE
      ekk=ekk*c1m21
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
      if(ium.ne.1) then
#include "include/kickq09v.f90"
      end if
#include "include/kicksho.f90"
#include "include/kickuxxv.f90"
      if(ium.eq.1) goto 350
      goto 330
    case (-10) ! SKEW 20-POLE
      ekk=ekk*c1m24
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
#include "include/kicksho.f90"
      if(ium.ne.1) then
#include "include/kickq10v.f90"
      end if
#include "include/kicksho.f90"
#include "include/kickuxxv.f90"
      if(ium.eq.1) goto 350

    case default
      goto 350
    end select
    goto 350
330 continue
    do j=2,ium
      if(kzz.eq.24) then
        y(j,1)=(y(j,1)+x(j,1)*qu)-qv*x(j,2)                          !hr06
        y(j,2)=(y(j,2)-x(j,2)*quz)-qvz*x(j,1)                        !hr06
      elseif(kzz.eq.25) then
        crkve=y(j,1)-(x(j,1)*qu)*qv                                  !hr06
        cikve=y(j,2)-(x(j,2)*qu)*qv                                  !hr06
        y(j,1)=crkve*cos_mb(qv)+cikve*sin_mb(qv)                     !hr09
        y(j,2)=cikve*cos_mb(qv)-crkve*sin_mb(qv)                     !hr09
        crkve=x(j,1)*cos_mb(qv)+x(j,2)*sin_mb(qv)                    !hr09
        cikve=x(j,2)*cos_mb(qv)-x(j,1)*sin_mb(qv)                    !hr09
        x(j,1)=crkve
        x(j,2)=cikve
      else
        y(j,1)=(y(j,1)+x(j,1)*qu)-qv*x(j,2)                          !hr06
        y(j,2)=(y(j,2)-x(j,2)*qu)-qv*x(j,1)                          !hr06
      endif
    end do
350 continue

  return

end subroutine umlauf

!-----------------------------------------------------------------------
!  CALCULATION OF DRIVINGTERMS OF RESONANCES INCLUDING SUBRESONANCE
!  USED FOR RMOD
!-----------------------------------------------------------------------
subroutine resex(dpp)
  ! Modified for Fortran 2015 by V.K.B. Olsen on 19/11/2017
  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  use crcoall
  use parpro
  use mod_common
  use mod_commons
  use mod_commont
  implicit none
  integer i,i1,i2,ii,ik,im,ip,ium,ix,izu,j,jj,jk,jl,jm,k,k1,kpz,kzz,l,l1,l2,ll,lmin,m2,m4,m6,min,&
          mm,mpe,mx,n,n2,n2e,nf1,nf3,nf4,nkk,nmz,nn1,nn2,nnf,np,np2,ns,nv,nv1,nv11,nv2,nv21,nz2,dj
  real(kind=fPrec) aa,ab1,ab2,alfa,b,b1,b2,bb,benkr,beta,btc,bts,chy,ci,cikve,cr,crkve,cxzi,cxzr,&
          cxzyi,cxzyr,cxzyrr,del,dphi,dpp,dppi,dpr,dt,dyy1,dyy2,e,ea,eb,ekk,ep,etl,gerad,phi,phibf,&
          phy,pie,puf,qu,qv,qw,r0,r0a,radi,re,re1,res,rn2,sb1,sb2,sea,seb,shy,t,vdt1,vdt2,vdt3,xl,&
          xs,zl,zs,quz,qvz
#ifdef TILT
  real(kind=fPrec) dyy11,qu1,tiltck,tiltck1,tiltck2,tiltck3,tiltck4,tiltck5,tiltckuk,tiltsk,&
          tiltsk1,tiltsk2,tiltsk3,tiltsk4,tiltsk5
#endif
  dimension t(5,4)
  dimension beta(2),alfa(2),phi(2),phibf(2)
  dimension qw(2)
  dimension aa(mmul),bb(mmul),dpr(5)
  dimension nnf(10),ep(2)
  dimension ab1(10),ab2(10),re(10,18),ip(10,18)
  dimension b(10,10),nz2(9),e(10,10)
  dimension chy(9,18),shy(9,18),min(5)
  dimension cr(mmul),ci(mmul)
      save

      ium=5
      do i=1,ium
        dpr(i)=zero
      end do

      do i=1,ium
        do j=1,4
          t(i,j)=zero
        end do
      end do

      do i=1,2
        beta(i)=zero
        alfa(i)=zero
        phi(i)=zero
        phibf(i)=zero
        qw(i)=zero
        ep(i)=zero
      end do

      do i=1,10
        nnf(i)=0
        do j=1,18
          ip(i,j)=0
          re(i,j)=zero
        end do
      end do

      do i=1,mmul
        aa(i)=zero
        bb(i)=zero
        cr(i)=zero
        ci(i)=zero
      end do

      do i=1,9
        nz2(i)=0
        do j=1,18
          chy(i,j)=zero
          shy(i,j)=zero
          do k=1,10
            do ii=1,10
              e(k,ii)=zero
              b(k,ii)=zero
            end do
            do l=1,5
              rtc(i,j,k,l)=zero
              rts(i,j,k,l)=zero
              min(l)=0
            end do
          end do
        end do
      end do

      btc=zero
      bts=zero
      phy=zero
      dt=zero
      del=zero
      ns=0
      ik=0
      pie=two*pi
      etl=zero
      radi=totl/pie
      dpr(1)=dpp*c1e3

      call clorb(dpp)
      call betalf(dpp,qw)

      if(ierro /= 0) then
        write(lout,"(a)") "RESEX> ERROR No optical solution."
        call prror(-1)
      end if
      call envar(dpp)

!--STARTVALUES OF THE TRAJECTORIES
      do l=1,2
        ll=2*l
        alfa(l)=alf0(l)
        beta(l)=bet0(l)
        t(1,ll-1)=clo(l)
         t(1,ll)=clop(l)
      end do

      do i=1,4
        do j=1,4
          t(i+1,j)=ta(j,i)
          t(i+1,j)=ta(j,i)
         end do
      end do

!--EP=EMITTANCE IN PI*MM*MRAD
      ep(1)=tam1**2/beta(1)                                              !hr06
      ep(2)=tam2**2/beta(2)                                              !hr06

!--SINGLE TURN BLOCKLOOP
      izu=0
      do 770 k=1,iu

        do k1=1,10
          ab1(k1)=zero
          ab2(k1)=zero
        end do

        ix=ic(k)
        if(ix.gt.nblo) goto 210
        jj=0
        dj=1
        if(ix.gt.0) goto 140
        ix=-ix
        jj=mel(ix)+1
        dj=-1
  140   jm=mel(ix)
!--BLOCKELEMENTLOOP
        do 200 j=1,jm
          jj=jj+dj
          jk=mtyp(ix,jj)
          if(ithick.eq.1.and.kz(jk).ne.0) goto 170
          if(ithick.eq.0.and.kz(jk).ne.0) goto 770
!--PURE DRIFTLENGTH
          etl=etl+el(jk)

          do l=1,2
            ll=2*l
            if(abs(t(ll,ll-1)).gt.pieni) then
              phibf(l)=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
            else
              phibf(l)=pi2
            endif
            do i=1,ium
              t(i,ll-1)=t(i,ll-1)+t(i,ll)*(el(jk))
            end do
          end do

          do l=1,2
            ll=2*l
            beta(l)=t(ll,ll-1)**2+t(ll+1,ll-1)**2                        !hr06
            alfa(l)=-one*(t(ll,ll-1)*t(ll,ll)+t(ll+1,ll-1)*t(ll+1,ll))   !hr06
            if(abs(t(ll,ll-1)).gt.pieni) then
              dphi=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
            else
              dphi=pi2-phibf(l)
            endif
            if(-one*dphi.gt.pieni) dphi=dphi+pi
            phi(l)=phi(l)+dphi
          end do

          goto 200
!--MAGNETELEMENT
  170     continue
          if(kz(jk).ne.8) etl=etl+el(jk)
          do l=1,2
            ll=2*l
            if(abs(t(ll,ll-1)).gt.pieni) then
              phibf(l)=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
            else
              phibf(l)=zero
            endif
            do i=1,ium
              puf=t(i,ll-1)
            t(i,ll-1)=(puf*a(jk,l,1)+t(i,ll)*a(jk,l,2))+dpr(i)*a(jk,l,5) !hr06
            t(i,ll)=(puf*a(jk,l,3)+t(i,ll)*a(jk,l,4))+dpr(i)*a(jk,l,6)   !hr06
            enddo
          enddo
          do l=1,2
            ll=2*l
            beta(l)=t(ll,ll-1)**2+t(ll+1,ll-1)**2                        !hr06
            alfa(l)=-one*(t(ll,ll-1)*t(ll,ll)+t(ll+1,ll-1)*t(ll+1,ll))   !hr06
            if(abs(t(ll,ll-1)).gt.pieni) then
              dphi=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
            else
              dphi=-one*phibf(l)
            endif
            if(kz(jk).ne.8.and.-dphi.gt.pieni) dphi=dphi+pi
            phi(l)=phi(l)+dphi
          enddo
  200   continue
        goto 770
!--NL-INSERTION
  210   ix=ix-nblo
        qu=zero
        qv=zero
        kpz=kp(ix)
        if(kpz.eq.6) goto 770
        kzz=kz(ix)
        if(kzz == 22) then
          do l=1,2
            ll=2*l
            if(abs(t(ll,ll-1)).gt.pieni) then
              phibf(l)=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
            else
              phibf(l)=zero
            end if
            do i=1,ium
              puf=t(i,ll-1)
              t(i,ll-1)=(puf*rrtr(imtr(ix),ll-1,ll-1)+t(i,ll)*rrtr(imtr(ix),ll-1,ll))+dpr(i)*rrtr(imtr(ix),ll-1,6)
              t(i,ll)=(puf*rrtr(imtr(ix),ll,ll-1)+t(i,ll)*rrtr(imtr(ix),ll,ll))+dpr(i)*rrtr(imtr(ix),ll,6)
            end do
            t(1,ll-1)=t(1,ll-1)+cotr(imtr(ix),ll-1)
            t(1,ll)=t(1,ll)+cotr(imtr(ix),ll)
            beta(l)=t(ll,ll-1)**2+t(ll+1,ll-1)**2
            alfa(l)=-one*(t(ll,ll-1)*t(ll,ll)+t(ll+1,ll-1)*t(ll+1,ll))
            if(abs(t(ll,ll-1)) > pieni) then
              dphi=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
            else
              dphi=-one*phibf(l)
            end if
            if(-one*dphi.gt.pieni) dphi=dphi+pi
                        phi(l)=phi(l)+dphi
          enddo
        endif
        if(kzz.eq.0.or.kzz.eq.20.or.kzz.eq.22) goto 770
        if(kzz.eq.15) goto 770
! JBG RF CC Multipoles to 770
        if(kzz.eq.26.or.kzz.eq.27.or.kzz.eq.28) goto 770
        if(kzz.eq.-26.or.kzz.eq.-27.or.kzz.eq.-28) goto 770
        dyy1=zero
        dyy2=zero
        if(iorg.lt.0) mzu(k)=izu
        izu=mzu(k)+1
        ekk=(sm(ix)+zfz(izu)*ek(ix))/(one+dpp)
        izu=izu+1
        xs=xpl(ix)+zfz(izu)*xrms(ix)
        izu=izu+1
        zs=zpl(ix)+zfz(izu)*zrms(ix)
#include "include/alignl.f90"

        select case (kzz)
        case (1)  ! HORIZONTAL DIPOLE
          ekk=ekk*c1e3
#include "include/kicka01h.f90"
        case (2)  ! NORMAL QUADRUPOLE
#include "include/kicka02h.f90"
        case (3)  ! NORMAL SEXTUPOLE
          ekk=ekk*c1m3
#include "include/kicka03h.f90"
        case (4)  ! NORMAL OCTUPOLE
          ekk=ekk*c1m6
#include "include/kicka04h.f90"
        case (5)  ! NORMAL DECAPOLE
          ekk=ekk*c1m9
#include "include/kicka05h.f90"
        case (6)  ! NORMAL DODECAPOLE
          ekk=ekk*c1m12
#include "include/kicka06h.f90"
        case (7)  ! NORMAL 14-POLE
          ekk=ekk*c1m15
#include "include/kicka07h.f90"
        case (8)  ! NORMAL 16-POLE
          ekk=ekk*c1m18
#include "include/kicka08h.f90"
        case (9)  ! NORMAL 18-POLE
          ekk=ekk*c1m21
#include "include/kicka09h.f90"
        case (10) ! NORMAL 20-POLE
          ekk=ekk*c1m24
#include "include/kicka10h.f90"
        case (11)
          r0 = ek(ix)
          if(abs(dki(ix,1)).gt.pieni) then
            if(abs(dki(ix,3)).gt.pieni) then
#include "include/multl01.f90"
              do i=2,ium
#include "include/multl02.f90"
              end do
            else
#include "include/multl03.f90"
            end if
          end if
          if(abs(dki(ix,2)).gt.pieni) then
            if(abs(dki(ix,3)).gt.pieni) then
#include "include/multl04.f90"
              do i=2,ium
#include "include/multl05.f90"
              end do
            else
#include "include/multl06.f90"
            end if
          end if
          mpe = 9
          mx  = 0
          if(abs(r0).le.pieni) goto 770
          nmz = nmu(ix)
          if(nmz.eq.0) then
            izu = izu+2*mmul
            goto 770
          end if
          im    = irm(ix)
          r0a   = one
          benkr = ed(ix)/(one+dpp)
          cr(1) = one
          cr(2) = xl
          ci(2) = zl
          cxzyr = xl
          cxzyi = zl
          cxzr  = cxzyr
          cxzi  = cxzyi
          dyy1  = zero
          dyy2  = zero
          qu    = zero
          qv    = zero
          lmin  = 3
          if(nmz.eq.1) lmin=2
          do l=lmin,mmul
            cr(l)=zero
            ci(l)=zero
          end do
          do l=1,nmz
#include "include/multl13.f90"
          end do
#ifdef TILT
#include "include/multl07e.f90"
#endif
          izu = izu+2*mmul-2*nmz
          goto 480
        case (12,13,14,15,16,17,18,19,20,21,22,23)
          goto 480
        case (24) ! DIPEDGE ELEMENT
#include "include/kickadpe.f90"
        case (25) ! Solenoid
#include "include/kickaso1.f90"
        case (26,27,28)
          goto 480

        !-----------------
        !--SKEW ELEMENTS--
        !-----------------
        case (-1)  ! VERTICAL DIPOLE
          ekk=ekk*c1e3
#include "include/kicka01v.f90"
        case (-2)  ! SKEW QUADRUPOLE
#include "include/kicka02v.f90"
        case (-3)  ! SKEW SEXTUPOLE
          ekk=ekk*c1m3
#include "include/kicka03v.f90"
        case (-4)  ! SKEW OCTUPOLE
          ekk=ekk*c1m6
#include "include/kicka04v.f90"
        case (-5)  ! SKEW DECAPOLE
          ekk=ekk*c1m9
#include "include/kicka05v.f90"
        case (-6)  ! SKEW DODECAPOLE
          ekk=ekk*c1m12
#include "include/kicka06v.f90"
        case (-7)  ! SKEW 14-POLE
          ekk=ekk*c1m15
#include "include/kicka07v.f90"
        case (-8)  ! SKEW 16-POLE
          ekk=ekk*c1m18
#include "include/kicka08v.f90"
        case (-9)  ! SKEW 18-POLE
          ekk=ekk*c1m21
#include "include/kicka09v.f90"
        case (-10) ! SKEW 20-POLE
          ekk=ekk*c1m24
#include "include/kicka10v.f90"
        end select
        goto 770
  480   continue
        t(1,2)=t(1,2)+dyy1
        t(1,4)=t(1,4)+dyy2
        do 490 i=2,ium
          if(kzz.eq.24) then
            t(i,2)=(t(i,2)+t(i,1)*qu)-qv*t(i,3)                          !hr06
            t(i,4)=(t(i,4)-t(i,3)*quz)-qvz*t(i,1)                        !hr06
#include "include/phas1so1.f90"
#include "include/phas2so1.f90"
#include "include/phas3so1.f90"
          else
            t(i,2)=(t(i,2)+t(i,1)*qu)-qv*t(i,3)                          !hr06
            t(i,4)=(t(i,4)-t(i,3)*qu)-qv*t(i,1)                          !hr06
          endif
  490   continue
        do l=1,2
          ll=2*l
          alfa(l)=-one*(t(ll,ll-1)*t(ll,ll)+t(ll+1,ll-1)*t(ll+1,ll))     !hr06
        end do
        if(mpe.gt.9.or.(mpe.eq.9.and.nmz.le.1)) goto 770
        if(mpe.lt.nta) goto 770
        if(mpe.gt.nte) mpe=nte
        if(nta.gt.2) goto 520
        if(mx.eq.-1.or.mx.eq.1.or.mx.eq.2.or.mx.eq.3.or.mx.eq.4 .or.mx.eq.5.or.mx.eq.6.or.mx.eq.7) goto 520

!-----------------------------------------------------------------------
!  SKEW-QUADRUPOLE;MULTIPOLES UP TO 9-TH ORDER
!-----------------------------------------------------------------------
        do l=2,nmz
          l1=l-1
          ab2(2)=ab2(2)+real(l1,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))               !hr06
        end do

  520   b1=beta(1)
        b2=beta(2)
        sb1=sqrt(b1)
        sb2=sqrt(b2)
        b(3,1)=b1
        b(1,3)=b2
        b(2,2)=sb1*sb2
        if(nta.gt.3) goto 540
        if(mpe.eq.2.or.(mpe.eq.9.and.nmz.le.2)) goto 670
        if(mx.eq.1.or.mx.eq.2.or.mx.eq.3.or.mx.eq.4 .or.mx.eq.5.or.mx.eq.6.or.mx.eq.7) goto 540

!-----------------------------------------------------------------------
!  REGULAR-SEXTUPOLE;MULTIPOLES UP TO 9-TH ORDER
!-----------------------------------------------------------------------
        l2=1
        do l=3,nmz
          l1=l-2
          ab1(3)=ab1(3)+real(l2,fPrec)*(aa(l)*ci(l1)+bb(l)*cr(l1))             !hr06
          ab2(3)=ab2(3)+real(l2,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))             !hr06
          l2=l2*l/l1
        end do

  540   b(4,1)=b1*sb1
        b(1,4)=b2*sb2
        b(3,2)=b1*sb2
        b(2,3)=b2*sb1
        if(nta.gt.4) goto 560
        if(mpe.eq.3.or.(mpe.eq.9.and.nmz.le.3)) goto 670
        if(mx.eq.2.or.mx.eq.3.or.mx.eq.4 .or.mx.eq.5.or.mx.eq.6.or.mx.eq.7) goto 560

!-----------------------------------------------------------------------
!  REGULAR-OCTUPOLE;MULTIPOLES UP TO 9-TH ORDER
!-----------------------------------------------------------------------
        l2=1
        do l=4,nmz
          l1=l-3
          ab1(4)=ab1(4)+real(l2,fPrec)*(aa(l)*ci(l1)+bb(l)*cr(l1))             !hr06
          ab2(4)=ab2(4)+real(l2,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))             !hr06
          l2=l2*l/l1
        end do

  560   b(5,1)=b1**2                                                     !hr06
        b(1,5)=b2**2                                                     !hr06
        b(4,2)=b(3,2)*sb1
        b(2,4)=b(2,3)*sb2
        b(3,3)=b1*b2
        if(nta.gt.5) goto 580
        if(mpe.eq.4.or.(mpe.eq.9.and.nmz.le.4)) goto 670
        if(mx.eq.3.or.mx.eq.4 .or.mx.eq.5.or.mx.eq.6.or.mx.eq.7) goto 580

!-----------------------------------------------------------------------
!  REGULAR-DEKAPOLE;MULTIPOLES UP TO 9-TH ORDER
!-----------------------------------------------------------------------
        l2=1
        do l=5,nmz
          l1=l-4
          ab1(5)=ab1(5)+real(l2,fPrec)*(aa(l)*ci(l1)+bb(l)*cr(l1))             !hr06
          ab2(5)=ab2(5)+real(l2,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))             !hr06
          l2=l2*l/l1
        end do

  580   b(6,1)=b(5,1)*sb1
        b(1,6)=b(1,5)*sb2
        b(5,2)=b(4,2)*sb1
        b(2,5)=b(2,4)*sb2
        b(4,3)=b(4,2)*sb2
        b(3,4)=b(2,4)*sb1
        if(nta.gt.6) goto 600
        if(mpe.eq.5.or.(mpe.eq.9.and.nmz.le.5)) goto 670
        if(mx.eq.4 .or.mx.eq.5.or.mx.eq.6.or.mx.eq.7) goto 600

!-----------------------------------------------------------------------
!  REGULAR-12-POLE;MULTIPOLES UP TO 9-TH ORDER
!-----------------------------------------------------------------------
        l2=1
        do l=6,nmz
          l1=l-5
          ab1(6)=ab1(6)+real(l2,fPrec)*(aa(l)*ci(l1)+bb(l)*cr(l1))             !hr06
          ab2(6)=ab2(6)+real(l2,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))             !hr06
          l2=l2*l/l1
        end do

  600   b(7,1)=b(6,1)*sb1
        b(1,7)=b(1,6)*sb2
        b(6,2)=b(5,2)*sb1
        b(2,6)=b(2,5)*sb2
        b(5,3)=b(5,2)*sb2
        b(3,5)=b(2,5)*sb1
        b(4,4)=b(3,4)*sb1
        if(nta.gt.7) goto 620
        if(mpe.eq.6.or.(mpe.eq.9.and.nmz.le.6)) goto 670
        if(mx.eq.5.or.mx.eq.6.or.mx.eq.7) goto 620

!-----------------------------------------------------------------------
!  REGULAR-14-POLE;MULTIPOLES UP TO 9-TH ORDER
!-----------------------------------------------------------------------
        l2=1
        do l=7,nmz
          l1=l-6
          ab1(7)=ab1(7)+real(l2,fPrec)*(aa(l)*ci(l1)+bb(l)*cr(l1))             !hr06
          ab2(7)=ab2(7)+real(l2,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))             !hr06
          l2=l2*l/l1
        end do

  620   b(8,1)=b(7,1)*sb1
        b(1,8)=b(1,7)*sb2
        b(7,2)=b(7,1)*sb2
        b(2,7)=b(1,7)*sb1
        b(6,3)=b(5,3)*sb1
        b(3,6)=b(3,5)*sb2
        b(5,4)=b(4,4)*sb1
        b(4,5)=b(4,4)*sb2
        if(nta.gt.8) goto 640
        if(mpe.eq.7.or.(mpe.eq.9.and.nmz.le.7)) goto 670
        if(mx.eq.6.or.mx.eq.7) goto 640
!-----------------------------------------------------------------------
!  REGULAR-16-POLE;MULTIPOLES UP TO 9-TH ORDER
!-----------------------------------------------------------------------
        l2=1
        do l=8,nmz
          l1=l-7
          ab1(8)=ab1(8)+real(l2,fPrec)*(aa(l)*ci(l1)+bb(l)*cr(l1))             !hr06
          ab2(8)=ab2(8)+real(l2,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))             !hr06
          l2=l2*l/l1
        end do

  640   b(9,1)=b(8,1)*sb1
        b(1,9)=b(1,8)*sb2
        b(8,2)=b(8,1)*sb2
        b(2,8)=b(1,8)*sb1
        b(7,3)=b(7,2)*sb2
        b(3,7)=b(2,7)*sb1
        b(6,4)=b(6,3)*sb2
        b(4,6)=b(3,6)*sb1
        b(5,5)=b(4,5)*sb1
        if(mpe.eq.8.or.(mpe.eq.9.and.nmz.le.8)) goto 670
        if(mx.eq.7) goto 660
!-----------------------------------------------------------------------
!  REGULAR-18-POLE
!-----------------------------------------------------------------------
        l2=1
        do l=9,nmz
          l1=l-8
          ab1(9)=ab1(9)+real(l2,fPrec)*(aa(l)*ci(l1)+bb(l)*cr(l1))             !hr06
          ab2(9)=ab2(9)+real(l2,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))             !hr06
          l2=l2*l/l1
        end do

  660   b(10,1)=b(9,1)*sb1
        b(1,10)=b(1,9)*sb2
        b(9,2)=b(9,1)*sb2
        b(2,9)=b(1,9)*sb1
        b(8,3)=b(8,2)*sb2
        b(3,8)=b(2,8)*sb1
        b(4,7)=b(3,7)*sb1
        b(7,4)=b(7,3)*sb2
        b(5,6)=b(4,6)*sb1
        b(6,5)=b(6,4)*sb2
!-----------------------------------------------------------------------
  670   do 700 np=1,mpe
          n2e=2*np
          do 690 nv=1,n2e
            n2=nv-np
            nn2=abs(n2)
            nn1=np-nn2
            re1=real(nn1,fPrec)*qxt+real(n2,fPrec)*qzt                   !hr06
            ipt=0

            do ii=1,nre
              if(n2.eq.nrr(ii)) ipt=ipr(ii)
            end do

            ip(np,nv)=int(re1+half)+ipt
            if(-one*re1.gt.pieni) ip(np,nv)=-int(abs(re1)+half)-ipt
!--RE=DISTANCE FROM THE RESONANCE
            re(np,nv)=re1-real(ip(np,nv),fPrec)                          !hr06
            res=re(np,nv)/radi
           chy(np,nv)=cos_mb((real(nn1,fPrec)*phi(1)+real(n2,fPrec)*phi(2))-res*etl) !hr06
           shy(np,nv)=sin_mb((real(nn1,fPrec)*phi(1)+real(n2,fPrec)*phi(2))-res*etl) !hr06
  690     continue
  700   continue
        do 760 np=nta,mpe
          np2=np
          nkk=0
  710     nkk=nkk+1
          n2e=2*np2
          do 750 i=1,nkk
            do 740 nv=1,n2e
              nn2=abs(nv-np2)
              nv1=np2-nn2+(i-1)*2+1
              nv2=np-nv1+2
              rn2=real(nn2,fPrec)*half                                    !hr06
!--EVENESS OF N2
              mm=0
              gerad=rn2-aint(rn2)
              if(abs(gerad).le.pieni) mm=1
!--MM=0 =>N2 UNEVEN, MM=1 => N2 EVEN
              if (mm.eq.0) goto 720
              btc=(ab1(np)*b(nv1,nv2))*chy(np2,nv)                       !hr06
              bts=(ab1(np)*b(nv1,nv2))*shy(np2,nv)                       !hr06
              goto 730
  720         btc=(ab2(np)*b(nv1,nv2))*chy(np2,nv)                       !hr06
              bts=(ab2(np)*b(nv1,nv2))*shy(np2,nv)                       !hr06
  730         rtc(np2,nv,np,i)=rtc(np2,nv,np,i)+btc
              rts(np2,nv,np,i)=rts(np2,nv,np,i)+bts
  740       continue
  750     continue
          np2=np2-2
          if(np2.ge.1) goto 710
  760   continue
  770 continue
      nnf(1)=1
      nnf(2)=1
      nnf(3)=2
      nz2(2)=2
      sea=sqrt(ep(1))
      seb=sqrt(ep(2))
      ea=ep(1)
      eb=ep(2)
      e(3,1)=one/eb
      e(1,3)=one/ea
      e(2,2)=(one/seb)/sea                                               !hr06
      nnf(4)=6
      nz2(3)=4
      e(4,1)=sea/eb
      e(1,4)=seb/ea
      e(3,2)=one/seb
      e(2,3)=one/sea
      nnf(5)=24
      nz2(4)=8
      e(5,1)=ea/eb
      e(1,5)=eb/ea
      e(4,2)=sea/seb
      e(2,4)=seb/sea
      e(3,3)=one
      nnf(6)=120
      nz2(5)=16
      e(6,1)=e(5,1)*sea
      e(1,6)=e(1,5)*seb
      e(5,2)=ea/seb
      e(2,5)=eb/sea
      e(4,3)=sea
      e(3,4)=seb
      nnf(7)=720
      nz2(6)=32
      e(7,1)=e(6,1)*sea
      e(1,7)=e(1,6)*seb
      e(6,2)=e(5,2)*sea
      e(2,6)=e(2,5)*seb
      e(5,3)=ea
      e(3,5)=eb
      e(4,4)=sea*seb
      nnf(8)=5040
      nz2(7)=64
      e(8,1)=e(7,1)*sea
      e(1,8)=e(1,7)*seb
      e(7,2)=e(6,2)*sea
      e(2,7)=e(2,6)*seb
      e(6,3)=ea*sea
      e(3,6)=eb*seb
      e(5,4)=ea*seb
      e(4,5)=sea*eb
      nnf(9)=40320
      nz2(8)=128
      e(9,1)=e(8,1)*sea
      e(1,9)=e(1,8)*seb
      e(8,2)=e(7,2)*sea
      e(2,8)=e(2,7)*seb
      e(7,3)=ea**2                                                       !hr06
      e(3,7)=eb**2                                                       !hr06
      e(6,4)=e(5,4)*sea
      e(4,6)=e(4,5)*seb
      e(5,5)=ea*eb
      nnf(10)=362880
      nz2(9)=256
      e(10,1)=e(9,1)*sea
      e(1,10)=e(1,9)*seb
      e(9,2)=e(8,2)*sea
      e(2,9)=e(2,8)*seb
      e(8,3)=e(7,3)*sea
      e(3,8)=e(3,7)*seb
      e(7,4)=e(6,4)*sea
      e(4,7)=e(4,6)*seb
      e(6,5)=e(5,5)*sea
      e(5,6)=e(5,5)*seb
      do 810 np=nta,nte
        vdt1=real(nnf(np),fPrec)/(real(nz2(np),fPrec)*pi)                            !hr06
        np2=np
        nkk=0
  780   nkk=nkk+1
        n2e=2*np2
        do 800 i=1,nkk
          do 790 nv=1,n2e
            n2=nv-np2
            nn2=abs(n2)
            nn1=np2-nn2
            nv1=nn1+(i-1)*2+1
            nv2=np-nv1+2
            nv11=nv1-1
            nv21=nv2-1
            nf1=nn1+i
            nf3=nkk-i+1
            nf4=nf3+nn2
      vdt2=(vdt1*e(nv1,nv2))/real(((nnf(nf1)*nnf(i))*nnf(nf3))*nnf(nf4),fPrec) !hr06
            vdt3=real(nn2,fPrec)*ea+real(nn1,fPrec)*eb                         !hr06
            if(n2.ge.0) vdt3=real(n2*nv21,fPrec)*ea+real(nn1*nv11,fPrec)*eb    !hr06
            rtc(np2,nv,np,i)=rtc(np2,nv,np,i)*vdt2*vdt3
            rts(np2,nv,np,i)=rts(np2,nv,np,i)*vdt2*vdt3
  790     continue
  800   continue
        np2=np2-2
        if(np2.ge.1) goto 780
  810 continue
      if(nur.eq.0) goto 840

      do 830 j=1,nur
        jk=j*2
        do i=1,nur
          jl=nu(i)-npp-jk
          if(jl.eq.0) min(j)=1
          if(jl.eq.0) goto 830
        end do
  830 continue

  840 m2=npp+2
      m4=npp+4
      m6=npp+6
      do 850 i=1,nre
        i2=2*i
        i1=i2-1
        n=nrr(i)+npp
        dtr(i1)=rtc(npp,n,npp,1)+(min(1)*(rtc(npp,n,m2,2)-              &!hr06
     &rtc(npp,n,m2,1))+min(2)*((rtc(npp,n,m4,1)-rtc(npp,n,m4,2))+rtc    &!hr06
     &(npp,n,m4,3)))+min(3)*(((rtc(npp,n,m6,2)-rtc(npp,n,m6,1))-rtc     &!hr06
     &(npp,n,m6,3))+ rtc(npp,n,m6,4))                                    !hr06
        dtr(i2)=rts(npp,n,npp,1)+(min(1)*(rts(npp,n,m2,2)-              &!hr06
     &rts(npp,n,m2,1))+min(2)*((rts(npp,n,m4,1)-rts(npp,n,m4,2))+rts    &!hr06
     &(npp,n,m4,3)))+min(3)*(((rts(npp,n,m6,2)-rts(npp,n,m6,1))-rts     &!hr06
     &(npp,n,m6,3))+rts(npp,n,m6,4))                                     !hr06
  850 continue
  return
end subroutine resex

subroutine rmod(dppr)
!-----------------------------------------------------------------------
!  CALCULATION OF THE STRENGTH OF CORRECTION-ELEMENTS
!-----------------------------------------------------------------------
  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  use crcoall
  use parpro
  use mod_common
  use mod_commons
  use mod_commont
  implicit none
  integer i,i1,i2,ierr,irr,j,j1,j2,j3,j4,jj1,jj2,jjr,k,n,no,ntao,nteo
  real(kind=fPrec) aa,bb,d1,de2,dpp,dppr,dsm,ox,oz,qwc,se11,se12,se2,sen,sen15,sen16,sen17,sen18,sn,ss
      dimension aa(10,10),bb(10),dsm(10),sn(10),sen(10),ss(10)
      dimension qwc(3),d1(10),irr(12)
      save
!-----------------------------------------------------------------------
      ntao=nta
      nteo=nte
      nta=npp
      nte=npp
      dpp=dppr

      do i=1,10
        bb(i)=zero
        dsm(i)=zero
        sn(i)=zero
        sen(i)=zero
        ss(i)=zero
        d1(i)=zero

        do j=1,10
          aa(i,j)=zero
        end do
      end do

      do i=1,12
        irr(i)=0
      end do

      do i=1,3
        qwc(i)=zero
      end do

      k=1
      jj1=0
      jj2=0
      jjr=2*nre
      de2=de0*half
      if(nre.eq.0) goto 50
      write(lout,10000)
      write(lout,10010) npp,totl,qxt,qzt,tam1
      call resex(dpp)
      do 40 i=1,nre
        i2=2*i
        i1=i2-1
        irr(i1)=ire(i1)
        irr(i2)=ire(i2)
        sn(i1)=ed(irr(i1))
        sn(i2)=ed(irr(i2))
        dsm(i1)=dsm0
        dsm(i2)=dsm0
        write(lout,10020) i,nrr(i),ipr(i)
        sen(i1)=dtr(i1)
        bb(i1)=sen(i1)
        sen(i2)=dtr(i2)
        bb(i2)=sen(i2)
        ss(i1)=sen(i1)
        ss(i2)=sen(i2)
   40 continue
      j2=jjr
   50 if(nur.eq.0) goto 70
      write(lout,10030) nur
      do 60 i=1,nur
        write(lout,10040) nu(i),i
   60 continue
   70 if(nch.eq.0) goto 90
      write(lout,10050)
      j1=j2+1
      j2=j2+2
      irr(j1)=ire(7)
      irr(j2)=ire(8)
      sn(j1)=ed(irr(j1))
      sn(j2)=ed(irr(j2))
      dsm(j1)=dsm0
      dsm(j2)=dsm0
      se2=zero
      se11=zero
      se12=zero
      do 80 n=1,5
        dpp=de2*real(3-n,fPrec)                                                !hr06
        call clorb2(dpp)
        call phasad(dpp,qwc)
        ox=qwc(1)
        oz=qwc(2)
        se2=se2+dpp*dpp
        se11=se11+ox*dpp
        se12=se12+oz*dpp
   80 continue
      sen(j1)=se11/se2
      sen(j2)=se12/se2
      bb(j1)=sen(j1)
      bb(j2)=sen(j2)
      ss(j1)=sen(j1)
      ss(j2)=sen(j2)
   90 if(nqc.eq.0) goto 100
      write(lout,10060)
      j1=j2+1
      j2=j2+2
      jj1=j1
      jj2=j2
      irr(j1)=ire(9)
      irr(j2)=ire(10)
      if (abs(el(irr(j1))).le.pieni.or.abs(el(irr(j2))).le.pieni) then
        sn(j1)=ed(irr(j1))
        sn(j2)=ed(irr(j2))
      else
        sn(j1)=ek(irr(j1))
        sn(j2)=ek(irr(j2))
      endif
      dsm(j1)=dkq
      dsm(j2)=dkq
      dpp=zero
      call clorb2(dpp)
      call phasad(dpp,qwc)
      sen(j1)=qwc(1)
      sen(j2)=qwc(2)
      bb(j1)=sen(j1)-qw0(1)
      bb(j2)=sen(j2)-qw0(2)
      ss(j1)=sen(j1)
      ss(j2)=sen(j2)
  100 do 330 no=1,itcro
        do 160 i=1,j2
          if(i.ne.jj1.and.i.ne.jj2) ed(irr(i))=ed(irr(i))+dsm(i)
          if(i.eq.jj1.or.i.eq.jj2) then
            if (abs(el(irr(i))).le.pieni) then
              ed(irr(i))=ed(irr(i))+dsm(i)
            else
              ek(irr(i))=ek(irr(i))+dsm(i)
            endif
          endif
          if(kp(irr(i)).eq.5) call combel(irr(i))
          if(nre.eq.0) goto 120
          call resex(dpp)
          do 110 j=1,jjr
            aa(i,j)=(dtr(j)-ss(j))/dsm(i)
  110     continue
  120     if(nch.eq.0) goto 140
          j3=jjr+1
          j4=jjr+2
          se2=zero
          se11=zero
          se12=zero
          do 130 n=1,5
            dpp=de2*real(3-n,fPrec)                                            !hr06
            call clorb2(dpp)
            call phasad(dpp,qwc)
            ox=qwc(1)
            oz=qwc(2)
            se2=se2+dpp*dpp
            se11=se11+ox*dpp
            se12=se12+oz*dpp
  130     continue
          sen15=se11/se2
          sen16=se12/se2
          aa(i,j3)=(sen15-ss(j3))/dsm(i)
          aa(i,j4)=(sen16-ss(j4))/dsm(i)
  140     if(nqc.eq.0) goto 150
          dpp=zero
          call clorb2(dpp)
          call phasad(dpp,qwc)
          sen17=qwc(1)
          sen18=qwc(2)
          aa(i,j1)=(sen17-ss(j1))/dsm(i)
          aa(i,j2)=(sen18-ss(j2))/dsm(i)
  150     continue
          if(i.eq.jj1.or.i.eq.jj2) then
            if (abs(el(irr(i))).le.pieni) then
              ed(irr(i))=ed(irr(i))-dsm(i)
            else
              ek(irr(i))=ek(irr(i))-dsm(i)
            endif
          endif
          if(i.ne.jj1.and.i.ne.jj2)ed(irr(i))=ed(irr(i))-dsm(i)
          if(kp(irr(i)).eq.5) call combel(irr(i))
  160   continue
        call loesd(aa,bb,j2,10,ierr)
        if(ierr == 1) then
          write(lout,"(a)") "RMOD> ERROR Problems during matrix-inversion."
          call prror(-1)
        end if
        do 170 i=1,j2
          if(i.eq.jj1.or.i.eq.jj2) then
            if (abs(el(irr(i))).le.pieni) then
              ed(irr(i))=ed(irr(i))-bb(i)
            else
              ek(irr(i))=ek(irr(i))-bb(i)
            endif
          endif
          if(i.ne.jj1.and.i.ne.jj2)ed(irr(i))=ed(irr(i))-bb(i)
          if(kp(irr(i)).eq.5) call combel(irr(i))
  170   continue
        if(nre.eq.0) goto 190
        call resex(dpp)

        do i=1,jjr
          ss(i)=dtr(i)
          d1(i)=abs(ss(i))
        end do

  190   if(nch.eq.0) goto 210
        se2=zero
        se11=zero
        se12=zero
        do 200 n=1,5
          dpp=de2*real(3-n,fPrec)                                              !hr06
          call clorb2(dpp)
          call phasad(dpp,qwc)
          ox=qwc(1)
          oz=qwc(2)
          se2=se2+dpp*dpp
          se11=se11+ox*dpp
          se12=se12+oz*dpp
  200   continue
        ss(j3)=se11/se2
        ss(j4)=se12/se2
        d1(j3)=abs(ss(j3))
        d1(j4)=abs(ss(j4))
  210   if(nqc.eq.0) goto 220
        dpp=zero
        call clorb2(dpp)
        call phasad(dpp,qwc)
        ss(j1)=qwc(1)
        ss(j2)=qwc(2)
        d1(j1)=abs(qwc(1)-qw0(1))
        d1(j2)=abs(qwc(2)-qw0(2))
  220   write(lout,10070)
        if(nre.eq.0) goto 270
        write(lout,10080) no,nrr(1),sen(1),ss(1),sen(2),ss(2)
        if(nre.eq.1) goto 240

        do i=2,nre
          i2=2*i
          i1=i2-1
          write(lout,10090) nrr(i),sen(i1),ss(i1),sen(i2),ss(i2)
        end do

  240   write(lout,10100)
        write(lout,10110)bez(irr(1)),sn(1),ed(irr(1)),bez(irr(2)),sn(2),&
     &ed(irr(2))
        if(nre.eq.1) goto 260

        do i=2,nre
          i2=2*i
          i1=i2-1
        write(lout,10110)bez(irr(i1)),sn(i1),ed(irr(i1)),bez(irr(i2)),sn&
     &(i2), ed(irr(i2))
        end do

  260   write(lout,10070)
  270   if(nch.eq.0) goto 280
        write(lout,10120) sen(j3),ss(j3),sen(j4),ss(j4)
        write(lout,10110)bez(irr(j3)),sn(j3),ed(irr(j3)),bez(irr(j4)),sn&
     &(j4), ed(irr(j4))
        write(lout,10070)
  280   if(nqc.eq.0) goto 290
        write(lout,10130) qw0(1),qwc(1),qw0(2),qwc(2)
        if (abs(el(irr(j1))).le.pieni) then
          write(lout,10140) sn(j1),ed(irr(j1)),irr(j1),sn(j2),          &
     &ed(irr(j2)),                                                      &
     &irr(j2)
        else
          write(lout,10140) sn(j1),ek(irr(j1)),irr(j1),sn(j2),          &
     &ek(irr(j2)),                                                      &
     &irr(j2)
        endif

  290   do i=1,j2
          if(d1(i).gt.dsi) goto 310
        end do

        nta=ntao
        nte=nteo
        return

  310   do i=1,j2
          bb(i)=ss(i)
        end do

        if(nqc.eq.1) bb(j1)=bb(j1)-qw0(1)
        if(nqc.eq.1) bb(j2)=bb(j2)-qw0(2)
  330 continue
      nta=ntao
      nte=nteo
!-----------------------------------------------------------------------
      return
10000 format(t5,'---- ENTRY RMOD ----')
10010 format(/10x,'N=',i1,' IS THE ORDER OF RESONACE, THAT WILL BE',    &
     &' COMPENSATED'// 10x,'L=',f15.6,'; QX=',f10.5,'; QY=',f10.5,      &
     &'; AMAX=',f10.5)
10020 format(/10x,i1,' RESONANCE; NY=',i2,';CHANGE OF P=',i2)
10030 format(/10x,'NUMBER OF SUBRESONANCES THAT ARE CONSIDERED IS ',i2)
10040 format(/10x,'NU=',i2,' IS THE ',i1,' SUBRESONANCE-MULTIPOLE-ORDER'&
     &,i2)
10050 format(/10x,'CHROMATICITY IS COMPENSATED')
10060 format(/10x,'Q-VALUES ARE ADJUSTED')
10070 format(131('-'))
10080 format(/10x,'RESONANCE-CORRECTION     ITERATION #',i2// 15x,      &
     &'DRIVING-TERM',13x,'BEFORE         AFTER     COMPENSATION'// 10x, &
     &'NY=',i2,'  COS-COMPONENT  ',2g15.5/ 17x,'SIN-COMPONENT  ',2g15.5/&
     &)
10090 format(10x,'NY=',i2,'  COS-COMPONENT  ',2g15.5/ 17x,              &
     &'SIN-COMPONENT  ',2g15.5/)
10100 format(10x,'  ELEMENT NAME'/)
10130 format(10x,'Q-VARIATION' / 10x,                                   &
     &'Q-VALUE            THEORET.        AFTER     COMPENSATION'/ 10x, &
     &'HORIZONTAL     ',2g15.7/ 10x,'VERTICAL       ',2g15.7/)
10140 format(10x,'QUADRU.STRENGTH',2g15.8,'   INDEX ',i3/ 10x,          &
     &'               ',2g15.8,'         ',i3)
10120 format(10x,'CHROMATICITY-CORRECTION'/ 15x,'CHROMATICITY',13x,     &
     &'BEFORE         AFTER     COMPENSATION'// 19x,'HORIZONTAL   ',2g15&
     &.5/ 19x,'VERTICAL     ',2g15.5/ 10x,'   SEXTUPOLE'/)
10110 format(14x,a16,2x,g17.10,1x,g17.10/14x,a16,2x,g17.10,1x,g17.10)
end subroutine rmod

subroutine search(dpp)
!-----------------------------------------------------------------------
!  FINDING THE BEST POSITIONS FOR CORRECTION-ELEMENTS
!-----------------------------------------------------------------------
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      use crcoall
      use parpro
      use mod_common
      use mod_commons
      use mod_commont
      implicit none
      integer i,id,n21,n22,n23,ntao,nteo
      real(kind=fPrec) b,c,c1,c2,c3,d,dpp,e,f,g,s1,s2,s3
      character(len=mNameLen) ref
      save
!-----------------------------------------------------------------------
      ntao=nta
      nteo=nte
      nta=mp
      nte=mp
      ref='REFERENCE       '
      id=0
      write(lout,10010)
      write(lout,10000)
      write(lout,10010)
      write(lout,10020) mp
      write(lout,10010)
      write(lout,10030) m21,ise1,m22,ise2,m23,ise3
      write(lout,10010)
      write(lout,10040)
      write(lout,10010)
      n21=m21+mp
      n22=m22+mp
      n23=m23+mp
      ipt=ise1
      call subsea(dpp)
      c1=rtc(mp,n21,mp,1)
      s1=rts(mp,n21,mp,1)
      ipt=ise2
      call subsea(dpp)
      c2=rtc(mp,n22,mp,1)
      s2=rts(mp,n22,mp,1)
      ipt=ise3
      call subsea(dpp)
      c3=rtc(mp,n23,mp,1)
      s3=rts(mp,n23,mp,1)
      write(lout,10050) ref,id,c1,s1,c2,s2,c3,s3
      do 10 i=1,mesa
        ed(isea(i))=ed(isea(i))+dsm0
        if(kp(isea(i)).eq.5) call combel(isea(i))
        ipt=ise1
        call subsea(dpp)
        b=rtc(mp,n21,mp,1)-c1
        c=rts(mp,n21,mp,1)-s1
        ipt=ise2
        call subsea(dpp)
        d=rtc(mp,n22,mp,1)-c2
        e=rts(mp,n22,mp,1)-s2
        ipt=ise3
        call subsea(dpp)
        f=rtc(mp,n23,mp,1)-c3
        g=rts(mp,n23,mp,1)-s3
        write(lout,10050) bez(isea(i)),i,b,c,d,e,f,g
        ed(isea(i))=ed(isea(i))-dsm0
        if(kp(isea(i)).eq.5) call combel(isea(i))
   10 continue
      nta=ntao
      nte=nteo
!-----------------------------------------------------------------------
      return
10000 format(t5,'---- ENTRY SEARCH ----')
10010 format(1x ,131('-'))
10020 format(10x,///'RESONANCES OF ORDER',i4,'  ARE CONSIDERED'//)
10030 format(24x ,'|',6x,'NY =',i4,';D-P= ',i4,7x, '|',6x,'NY =',i4,    &
     &';D-P= ',i4,7x,'|',6x,'NY =',i4,';D-P= ',i4,7x, '|')
10040 format(1x,'ELEMENT          | POS |',6x,'COS',13x,'SIN',6x,'|',   &
     &6x,'COS',13x,'SIN',6x,'|', 6x,'COS',13x,'SIN',6x,'|')
10050 format(1x,a16,1x,'|',i3,'  |',g15.5,'|',g15.5,'|',g15.5,'|',      &
     &g15.5,'|',g15.5,'|',g15.5,'|')
end

!-----------------------------------------------------------------------
!  CALCULATION OF RESONANCE- AND SUBRESONANCE-DRIVINGTERMS
!-----------------------------------------------------------------------
subroutine subre(dpp)
  ! Rewritten to remove computed gotos by V.K.B.Olsen on 23/11/2017
  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  use crcoall
  use parpro
  use mod_common
  use mod_commons
  use mod_commont
  implicit none

  integer i,ii,ik,im,ip,ipc,ipcc,ipl,ium,iv,ix,izu,j,jj,jk,jm,k,k1,kpz,kzz,l,l1,l2,ll,lmin,min1,min2,&
          mis,mm,mpe,mx,n2,n22,n2e,nf1,nf3,nf4,nkk,nmz,nn1,nn2,nnf,np,np2,nph,nr,ns,ntx,nv,nv1,nv11,&
          nv2,nv21,nz2,dj
  real(kind=fPrec) aa,ab1,ab2,alfa,b,b1,b2,bb,benkr,beta,btc,bts,cc_r,chy,ci,cikve,clo0,clop0,cr,crkve,&
          cxzi,cxzr,cxzyi,cxzyr,cxzyrr,del,dfac,dphi,dpp,dpp1,dppi,dpr,dt,dtu,dtup,dyy1,dyy2,e,ea,eb,&
          ekk,ekko,ep,etl,gerad,gtu1,gtu2,phi,phibf,phy,pie,puf,qu,qv,qw,qwc,r0,r0a,radi,rc,re,re1,res,&
          rn2,rs,sb1,sb2,sdel,sdel2,sea,seb,shy,ss,t,vdt1,vdt2,vdt3,vdt4,xl,xs,zl,zs,quz,qvz
#ifdef TILT
  real(kind=fPrec) dyy11,qu1,tiltck,tiltck1,tiltck2,tiltck3,tiltck4,tiltck5,tiltck6,tiltck8,tiltck10,&
          tiltckuk,tiltsk,tiltsk1,tiltsk2,tiltsk3,tiltsk4,tiltsk5,tiltsk6,tiltsk8,tiltsk10
#endif
      dimension t(6,4)
      dimension beta(2),alfa(2),phi(2),phibf(2)
      dimension clo0(2),clop0(2)
      dimension aa(mmul),bb(mmul)
      dimension qw(2),qwc(3),dpr(6)
      dimension nnf(10),ep(2)
      dimension ab1(10),ab2(10),re(10,18),ip(10,18)
      dimension b(10,10),nz2(9),e(10,10)
      dimension chy(9,18),shy(9,18)
      dimension dfac(10),dtu(2,5),dtup(2,5,0:4,0:4)
      dimension cr(mmul),ci(mmul)
      save
!-----------------------------------------------------------------------
      ium=5
      ipl=1
      gtu1=zero
      gtu2=zero
      dfac(1)=one
      dfac(2)=one
      dfac(3)=two
      dfac(4)=six                                                        !hr13
      dfac(5)=24.0                                                       !hr13
      dfac(6)=120.0_fPrec                                                !hr13
      dfac(7)=720.0_fPrec                                                !hr13
      dfac(8)=5040.0_fPrec                                               !hr13
      dfac(9)=40320.0_fPrec                                              !hr13
      dfac(10)=362880.0_fPrec                                            !hr13
      if(ipt.eq.1) ipl=3

      do 940 ipcc=1,ipl
        ipc=ipcc-ipl+1
        if(ipt.eq.0) ipc=0
        btc=zero
        bts=zero
        phy=zero
        dt=zero
        del=zero
        ns=0
        ik=0

        do i=1,ium
          dpr(i)=zero
        end do

        do i=1,ium
          do j=1,4
            t(i,j)=zero
          end do
        end do

        do i=1,2
          beta(i)=zero
          alfa(i)=zero
          phi(i)=zero
          phibf(i)=zero
          qw(i)=zero
          qwc(i)=zero
          clo0(i)=zero
          clop0(i)=zero
          ep(i)=zero
        end do

        qwc(3)=zero
        do i=1,10
          nnf(i)=0
          do j=1,18
            re(i,j)=zero
            ip(i,j)=0
          end do
        end do

        do i=1,mmul
          aa(i)=zero
          bb(i)=zero
          cr(i)=zero
          ci(i)=zero
        end do

        do i=1,2
          do j=1,5
            dtu(i,j)=zero
          end do
        end do

        do i=1,5
          do j=0,4
            do k=0,4
              dtup(1,i,j,k)=zero
              dtup(2,i,j,k)=zero
            end do
          end do
        end do

        do 120 i=1,9
          nz2(i)=0
          do 110 j=1,18
            chy(i,j)=zero
            shy(i,j)=zero
            do 100 k=1,10
              do 80 ii=1,10
                e(k,ii)=zero
                b(k,ii)=zero
   80         continue
              do 90 l=1,5
                rtc(i,j,k,l)=zero
                rts(i,j,k,l)=zero
   90         continue
  100       continue
  110     continue
  120   continue

        write(lout,10030)
        write(lout,10020)
        pie=two*pi
        etl=zero
        radi=totl/pie
        nr=0
        dpr(1)=dpp*c1e3
        dpr(6)=c1e3
        dpp1=dpp+ded

        call clorb(dpp1)

        do l=1,2
          clo0(l)=clo(l)
          clop0(l)=clop(l)
        end do

        call clorb(dpp)

        do l=1,2
          di0(l)=(clo0(l)-clo(l))/ded
          dip0(l)=(clop0(l)-clop(l))/ded
        end do

        write(lout,10030)
        write(lout,10120) (di0(l),dip0(l),l=1,2)
        call betalf(dpp,qw)
        call phasad(dpp,qwc)
        if(ierro /= 0) then
          write(lout,"(a)") "SUBRE> ERROR No optical solution."
          call prror(-1)
        end if
        write(lout,10070) dpp,qwc(1),qwc(2)
        call envar(dpp)

!--STARTVALUES OF THE TRAJECTORIES
        do l=1,2
          ll=2*l
          alfa(l)=alf0(l)
          beta(l)=bet0(l)
          t(1,ll-1)=clo(l)
          t(1,ll)=clop(l)
          clo0(l)=clo(l)
          clop0(l)=clop(l)
        end do

        do i=1,4
          do j=1,4
            t(i+1,j)=ta(j,i)
            t(i+1,j)=ta(j,i)
          end do
        end do

        write(lout,10030)
        write(lout,10040)
        write(lout,10030)
        write(lout,10010) nr,'START   ',zero,zero,(beta(l),alfa(l),phi(l),&
                          di0(l),dip0(l),clo0(l),clop0(l),l=1,2)

!--EP=EMITTANCE IN PI*MM*MRAD
        ep(1)=tam1**2/beta(1)                                            !hr06
        ep(2)=tam2**2/beta(2)                                            !hr06
        write(lout,10050) tam1,ep(1),tam2,ep(2)
        write(lout,10030)

!--SINGLE TURN BLOCKLOOP
        izu=0
        do 790 k=1,iu
          do k1=1,10
            ab1(k1)=zero
            ab2(k1)=zero
          end do

          ix=ic(k)
          if(ix.gt.nblo) goto 250
          jj=0
          dj=1
          if(ix.gt.0) goto 180
          ix=-1*ix
          jj=mel(ix)+1
          dj=-1
  180     jm=mel(ix)

!--SINGLE TURN BLOCKLOOP
          do 240 j=1,jm
            jj=jj+dj
            jk=mtyp(ix,jj)
            if(ithick.eq.1.and.kz(jk).ne.0) goto 210
            if(ithick.eq.0.and.kz(jk).ne.0) goto 790

!--PURE DRIFTLENGTH
            etl=etl+el(jk)
            do l=1,2
              ll=2*l
              if(abs(t(ll,ll-1)).gt.pieni) then
                phibf(l)=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
              else
                phibf(l)=pi2
              endif
              do i=1,ium
                t(i,ll-1)=t(i,ll-1)+t(i,ll)*(el(jk))
              end do
            end do

            do l=1,2
              ll=2*l
              beta(l)=t(ll,ll-1)**2+t(ll+1,ll-1)**2                      !hr06
              alfa(l)=-one*(t(ll,ll-1)*t(ll,ll)+t(ll+1,ll-1)*t(ll+1,ll)) !hr06
              clo0(l)=t(1,ll-1)
              clop0(l)=t(1,ll)

              if(abs(t(ll,ll-1)).gt.pieni) then
                dphi=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
              else
                dphi=pi2-phibf(l)
              endif

              if(-one*dphi.gt.pieni) dphi=dphi+pi                        !hr06
              phi(l)=phi(l)+dphi/pie
            end do

            nr=nr+1
            goto 240

!--MAGNETELEMENT
  210       continue
            if(kz(jk).ne.8) etl=etl+el(jk)
            do l=1,2
              ll=2*l
              if(abs(t(ll,ll-1)).gt.pieni) then
                phibf(l)=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
              else
                phibf(l)=zero
              endif
              do i=1,ium
                puf=t(i,ll-1)
                t(i,ll-1)=(puf*a(jk,l,1)+t(i,ll)*a(jk,l,2))+dpr(i)*a(jk,l,5)
                t(i,ll)=(puf*a(jk,l,3)+t(i,ll)*a(jk,l,4))+dpr(i)*a(jk,l,6) !hr06
              enddo
            enddo
            do l=1,2
              ll=2*l
              beta(l)=t(ll,ll-1)**2+t(ll+1,ll-1)**2                      !hr06
              alfa(l)=-one*(t(ll,ll-1)*t(ll,ll)+t(ll+1,ll-1)*t(ll+1,ll)) !hr06
              clo0(l)=t(1,ll-1)
              clop0(l)=t(1,ll)
              if(abs(t(ll,ll-1)).gt.pieni) then
                dphi=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
              else
                dphi=-phibf(l)
              endif
              if(kz(jk).ne.8.and.-dphi.gt.pieni) dphi=dphi+pi
              phi(l)=phi(l)+dphi/pie
            enddo
            nr=nr+1
  240     continue
          goto 790

!--NL-INSERTION
  250     ix=ix-nblo
          qu=zero
          qv=zero
          kpz=kp(ix)
          if(kpz.eq.6) goto 790
          kzz=kz(ix)
          if(kzz == 22) then
            do l=1,2
              ll=2*l
              if(abs(t(ll,ll-1)).gt.pieni) then
                phibf(l)=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
              else
                phibf(l)=zero
              end if
              do i=1,ium
                puf=t(i,ll-1)
                t(i,ll-1)=(puf*rrtr(imtr(ix),ll-1,ll-1)+t(i,ll)*rrtr(imtr(ix),ll-1,ll))+dpr(i)*rrtr(imtr(ix),ll-1,6)
                t(i,ll)=(puf*rrtr(imtr(ix),ll,ll-1)+t(i,ll)*rrtr(imtr(ix),ll,ll))+dpr(i)*rrtr(imtr(ix),ll,6)
              end do
              t(1,ll-1)=t(1,ll-1)+cotr(imtr(ix),ll-1)
              t(1,ll)=t(1,ll)+cotr(imtr(ix),ll)
              beta(l)=t(ll,ll-1)**2+t(ll+1,ll-1)**2
              alfa(l)=-one*(t(ll,ll-1)*t(ll,ll)+t(ll+1,ll-1)*t(ll+1,ll))
              if(abs(t(ll,ll-1)) > pieni) then
                dphi=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
              else
                dphi=-one*phibf(l)
              end if
              if(-one*dphi.gt.pieni) dphi=dphi+pi
                          phi(l)=phi(l)+dphi
          enddo
        endif
          clo0(1)=t(1,1)
          clop0(1)=t(1,2)
          clo0(2)=t(2,3)
          clop0(2)=t(2,4)
          if(kzz.eq.0.or.kzz.eq.20.or.kzz.eq.22) goto 790
          if(kzz.eq.15) goto 790
! JBG RF CC Multipoles to 790
          if(kzz.eq.26.or.kzz.eq.27.or.kzz.eq.28) goto 790
          if(kzz.eq.-26.or.kzz.eq.-27.or.kzz.eq.-28) goto 790
          dyy1=zero
          dyy2=zero
          if(iorg.lt.0) mzu(k)=izu
          izu=mzu(k)+1
          ekk=(sm(ix)+zfz(izu)*ek(ix))/(one+dpp)
          izu=izu+1
          xs=xpl(ix)+zfz(izu)*xrms(ix)
          izu=izu+1
          zs=zpl(ix)+zfz(izu)*zrms(ix)
#include "include/alignl.f90"

        select case (kzz)
        case (1) ! HORIZONTAL DIPOLE
          ekk=ekk*c1e3
#include "include/kicka01h.f90"
        case (2) ! NORMAL QUADRUPOLE
#include "include/kicka02h.f90"
        case (3) ! NORMAL SEXTUPOLE
          ekk=ekk*c1m3
#include "include/kicka03h.f90"
        case (4) ! NORMAL OCTUPOLE
          ekk=ekk*c1m6
#include "include/kicka04h.f90"
          call detune(2,ekk,ep,beta,dtu,dtup,dfac)
        case (5) ! NORMAL DECAPOLE
          ekk=ekk*c1m9
#include "include/kicka05h.f90"
          call detune(2,ab1(4),ep,beta,dtu,dtup,dfac)
        case (6) ! NORMAL DODECAPOLE
          ekk=ekk*c1m12
#include "include/kicka06h.f90"
          call detune(2,ab1(4),ep,beta,dtu,dtup,dfac)
          call detune(3,ekk,ep,beta,dtu,dtup,dfac)
        case (7) ! NORMAL 14-POLE
          ekk=ekk*c1m15
#include "include/kicka07h.f90"
          call detune(2,ab1(4),ep,beta,dtu,dtup,dfac)
          call detune(3,ab1(6),ep,beta,dtu,dtup,dfac)
        case (8) ! NORMAL 16-POLE
          ekk=ekk*c1m18
#include "include/kicka08h.f90"
          call detune(2,ab1(4),ep,beta,dtu,dtup,dfac)
          call detune(3,ab1(6),ep,beta,dtu,dtup,dfac)
          call detune(4,ekk,ep,beta,dtu,dtup,dfac)
        case (9) ! NORMAL 18-POLE
          ekk=ekk*c1m21
#include "include/kicka09h.f90"
          call detune(2,ab1(4),ep,beta,dtu,dtup,dfac)
          call detune(3,ab1(6),ep,beta,dtu,dtup,dfac)
          call detune(4,ab1(8),ep,beta,dtu,dtup,dfac)
        case (10) ! NORMAL 20-POLE
          ekk=ekk*c1m24
#include "include/kicka10h.f90"
#include "include/kispa10h.f90"
        case (11)
          r0=ek(ix)
          if(abs(dki(ix,1)).gt.pieni) then
            if(abs(dki(ix,3)).gt.pieni) then
#include "include/multl01.f90"
              do i=2,ium
#include "include/multl02.f90"
              end do
            else
#include "include/multl03.f90"
            end if
          end if
          if(abs(dki(ix,2)).gt.pieni) then
            if(abs(dki(ix,3)).gt.pieni) then
#include "include/multl04.f90"
              do i=2,ium
#include "include/multl05.f90"
              end do
            else
#include "include/multl06.f90"
            end if
          end if
          mpe=9
          mx=0
          if(abs(r0).le.pieni) goto 790
          nmz=nmu(ix)
          if(nmz.eq.0) then
            izu=izu+2*mmul
            goto 790
          end if
          im=irm(ix)
          r0a=one
          benkr=ed(ix)/(one+dpp)
          cr(1)=one
          cr(2)=xl
          ci(2)=zl
          cxzyr=xl
          cxzyi=zl
          cxzr=cxzyr
          cxzi=cxzyi
          dyy1=zero
          dyy2=zero
          qu=zero
          qv=zero
          lmin=3
          if(nmz.eq.1) lmin=2
          do l=lmin,mmul
            aa(l)=zero
            bb(l)=zero
            cr(l)=zero
            ci(l)=zero
          end do
          do l=1,nmz
#include "include/multl13.f90"
          end do
#ifdef TILT
#include "include/multl07e.f90"
#endif
          izu=izu+2*mmul-2*nmz
          do iv=2,5
#include "include/multl12.f90"
          end do
        case (12,13,14,15,16,17,18,19,20,21,22,23)
          goto 790
        case (24) ! DIPEDGE ELEMENT
#include "include/kickadpe.f90"
        case (25) ! Solenoid
#include "include/kickaso1.f90"
        case (26,27,28)
          goto 790

        !-----------------
        !--SKEW ELEMENTS--
        !-----------------
        case (-1) ! VERTICAL DIPOLE
          ekk=ekk*c1e3
#include "include/kicka01v.f90"
        case (-2) ! SKEW QUADRUPOLE
#include "include/kicka02v.f90"
        case (-3) ! SKEW SEXTUPOLE
          ekk=ekk*c1m3
#include "include/kicka03v.f90"
        case (-4) ! SKEW OCTUPOLE
          ekk=ekk*c1m6
#include "include/kicka04v.f90"
        case (-5) ! SKEW DECAPOLE
          ekk=ekk*c1m9
#include "include/kicka05v.f90"
          call detune(2,ab1(4),ep,beta,dtu,dtup,dfac)
        case (-6) ! SKEW DODECAPOLE
          ekk=ekk*c1m12
#include "include/kicka06v.f90"
          call detune(2,ab1(4),ep,beta,dtu,dtup,dfac)
        case (-7) ! SKEW 14-POLE
          ekk=ekk*c1m15
#include "include/kicka07v.f90"
          call detune(2,ab1(4),ep,beta,dtu,dtup,dfac)
          call detune(3,ab1(6),ep,beta,dtu,dtup,dfac)
        case (-8) ! SKEW 16-POLE
          ekk=ekk*c1m18
#include "include/kicka08v.f90"
          call detune(2,ab1(4),ep,beta,dtu,dtup,dfac)
          call detune(3,ab1(6),ep,beta,dtu,dtup,dfac)
        case (-9) ! SKEW 18-POLE
          ekk=ekk*c1m21
#include "include/kicka09v.f90"
          call detune(2,ab1(4),ep,beta,dtu,dtup,dfac)
          call detune(3,ab1(6),ep,beta,dtu,dtup,dfac)
          call detune(4,ab1(8),ep,beta,dtu,dtup,dfac)
        case (-10) ! SKEW 20-POLE
          ekk=ekk*c1m24
#include "include/kicka10v.f90"
#include "include/kispa10v.f90"
        case default
          goto 790
        end select

          t(1,2)=t(1,2)+dyy1
          t(1,4)=t(1,4)+dyy2
          do i=2,ium
            if(kzz.eq.24) then
              t(i,2)=(t(i,2)+t(i,1)*qu)-qv*t(i,3)                          !hr06
              t(i,4)=(t(i,4)-t(i,3)*quz)-qvz*t(i,1)                        !hr06
#include "include/phas1so1.f90"
#include "include/phas2so1.f90"
#include "include/phas3so1.f90"
            else
              t(i,2)=(t(i,2)+t(i,1)*qu)-qv*t(i,3)                          !hr06
              t(i,4)=(t(i,4)-t(i,3)*qu)-qv*t(i,1)                          !hr06
            end if
          end do

          do l=1,2
            ll=2*l
            alfa(l)=-one*(t(ll,ll-1)*t(ll,ll)+t(ll+1,ll-1)*t(ll+1,ll))   !hr06
            clop0(l)=t(1,ll)
          end do

          if(mpe.gt.9.or.(mpe.eq.9.and.nmz.le.1)) goto 790
          if(mpe.lt.nta) goto 790
          if(mpe.gt.nte) mpe=nte
          if(nta.gt.2) goto 550
          if(mx.eq.-1.or.mx.eq.1.or.mx.eq.2.or.mx.eq.3.or.mx.eq.4 .or. mx.eq.5.or.mx.eq.6.or.mx.eq.7) goto 550

!-----------------------------------------------------------------------
!  SKEW-QUADRUPOLE;MULTIPOLES UP TO 9-TH ORDER
!-----------------------------------------------------------------------
          do l=2,nmz
            l1=l-1
            ab2(2)=ab2(2)+real(l1,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))             !hr06
          end do

  550     b1=beta(1)
          b2=beta(2)
          sb1=sqrt(b1)
          sb2=sqrt(b2)
          b(3,1)=b1
          b(1,3)=b2
          b(2,2)=sb1*sb2
          if(nta.gt.3) goto 570
          if(mpe.eq.2.or.(mpe.eq.9.and.nmz.le.2)) goto 700
          if(mx.eq.1.or.mx.eq.2.or.mx.eq.3.or.mx.eq.4 .or.mx.eq.5.or.mx .eq.6.or.mx.eq.7) goto 570

!-----------------------------------------------------------------------
!  REGULAR-SEXTUPOLE;MULTIPOLES UP TO 9-TH ORDER
!-----------------------------------------------------------------------
          l2=1
          do l=3,nmz
            l1=l-2
            ab1(3)=ab1(3)+real(l2,fPrec)*(aa(l)*ci(l1)+bb(l)*cr(l1))           !hr06
            ab2(3)=ab2(3)+real(l2,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))           !hr06
            l2=l2*l/l1
          end do

  570     b(4,1)=b1*sb1
          b(1,4)=b2*sb2
          b(3,2)=b1*sb2
          b(2,3)=b2*sb1
          if(nta.gt.4) goto 590
          if(mpe.eq.3.or.(mpe.eq.9.and.nmz.le.3)) goto 700
          if(mx.eq.2.or.mx.eq.3.or.mx.eq.4 .or.mx.eq.5.or.mx.eq.6.or.mx .eq.7) goto 590

!-----------------------------------------------------------------------
!  REGULAR-OCTUPOLE;MULTIPOLES UP TO 9-TH ORDER
!-----------------------------------------------------------------------
          l2=1
          do l=4,nmz
            l1=l-3
            ab1(4)=ab1(4)+real(l2,fPrec)*(aa(l)*ci(l1)+bb(l)*cr(l1))     !hr06
            ab2(4)=ab2(4)+real(l2,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))     !hr06
            l2=l2*l/l1
          end do

  590     b(5,1)=b1**2                                                   !hr06
          b(1,5)=b2**2                                                   !hr06
          b(4,2)=b(3,2)*sb1
          b(2,4)=b(2,3)*sb2
          b(3,3)=b1*b2
          if(nta.gt.5) goto 610
          if(mpe.eq.4.or.(mpe.eq.9.and.nmz.le.4)) goto 700
          if(mx.eq.3.or.mx.eq.4 .or.mx.eq.5.or.mx.eq.6.or.mx.eq.7) goto 610

!-----------------------------------------------------------------------
!  REGULAR-DEKAPOLE;MULTIPOLES UP TO 9-TH ORDER
!-----------------------------------------------------------------------
          l2=1
          do l=5,nmz
            l1=l-4
            ab1(5)=ab1(5)+real(l2,fPrec)*(aa(l)*ci(l1)+bb(l)*cr(l1))    !hr06
            ab2(5)=ab2(5)+real(l2,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))    !hr06
            l2=l2*l/l1
          end do

  610     b(6,1)=b(5,1)*sb1
          b(1,6)=b(1,5)*sb2
          b(5,2)=b(4,2)*sb1
          b(2,5)=b(2,4)*sb2
          b(4,3)=b(4,2)*sb2
          b(3,4)=b(2,4)*sb1
          if(nta.gt.6) goto 630
          if(mpe.eq.5.or.(mpe.eq.9.and.nmz.le.5)) goto 700
          if(mx.eq.4 .or.mx.eq.5.or.mx.eq.6.or.mx.eq.7) goto 630

!-----------------------------------------------------------------------
!  REGULAR-12-POLE;MULTIPOLES UP TO 9-TH ORDER
!-----------------------------------------------------------------------
          l2=1
          do l=6,nmz
            l1=l-5
            ab1(6)=ab1(6)+real(l2,fPrec)*(aa(l)*ci(l1)+bb(l)*cr(l1))    !hr06
            ab2(6)=ab2(6)+real(l2,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))    !hr06
            l2=l2*l/l1
          end do

  630     b(7,1)=b(6,1)*sb1
          b(1,7)=b(1,6)*sb2
          b(6,2)=b(5,2)*sb1
          b(2,6)=b(2,5)*sb2
          b(5,3)=b(5,2)*sb2
          b(3,5)=b(2,5)*sb1
          b(4,4)=b(3,4)*sb1
          if(nta.gt.7) goto 650
          if(mpe.eq.6.or.(mpe.eq.9.and.nmz.le.6)) goto 700
          if(mx.eq.5.or.mx.eq.6.or.mx.eq.7) goto 650

!-----------------------------------------------------------------------
!  REGULAR-14-POLE;MULTIPOLES UP TO 9-TH ORDER
!-----------------------------------------------------------------------
          l2=1
          do l=7,nmz
            l1=l-6
            ab1(7)=ab1(7)+real(l2,fPrec)*(aa(l)*ci(l1)+bb(l)*cr(l1))    !hr06
            ab2(7)=ab2(7)+real(l2,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))    !hr06
            l2=l2*l/l1
          end do

  650     b(8,1)=b(7,1)*sb1
          b(1,8)=b(1,7)*sb2
          b(7,2)=b(7,1)*sb2
          b(2,7)=b(1,7)*sb1
          b(6,3)=b(5,3)*sb1
          b(3,6)=b(3,5)*sb2
          b(5,4)=b(4,4)*sb1
          b(4,5)=b(4,4)*sb2
          if(nta.gt.8) goto 670
          if(mpe.eq.7.or.(mpe.eq.9.and.nmz.le.7)) goto 700
          if(mx.eq.6.or.mx.eq.7) goto 670

!-----------------------------------------------------------------------
!  REGULAR-16-POLE;MULTIPOLES UP TO 9-TH ORDER
!-----------------------------------------------------------------------
          l2=1
          do l=8,nmz
            l1=l-7
            ab1(8)=ab1(8)+real(l2,fPrec)*(aa(l)*ci(l1)+bb(l)*cr(l1))    !hr06
            ab2(8)=ab2(8)+real(l2,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))    !hr06
            l2=l2*l/l1
          end do

  670     b(9,1)=b(8,1)*sb1
          b(1,9)=b(1,8)*sb2
          b(8,2)=b(8,1)*sb2
          b(2,8)=b(1,8)*sb1
          b(7,3)=b(7,2)*sb2
          b(3,7)=b(2,7)*sb1
          b(6,4)=b(6,3)*sb2
          b(4,6)=b(3,6)*sb1
          b(5,5)=b(4,5)*sb1
          if(mpe.eq.8.or.(mpe.eq.9.and.nmz.le.8)) goto 700
          if(mx.eq.7) goto 690

!-----------------------------------------------------------------------
!  REGULAR-18-POLE
!-----------------------------------------------------------------------
          l2=1
          do l=9,nmz
            l1=l-8
            ab1(9)=ab1(9)+real(l2,fPrec)*(aa(l)*ci(l1)+bb(l)*cr(l1))     !hr06
            ab2(9)=ab2(9)+real(l2,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))     !hr06
            l2=l2*l/l1
          end do

  690     b(10,1)=b(9,1)*sb1
          b(1,10)=b(1,9)*sb2
          b(9,2)=b(9,1)*sb2
          b(2,9)=b(1,9)*sb1
          b(8,3)=b(8,2)*sb2
          b(3,8)=b(2,8)*sb1
          b(4,7)=b(3,7)*sb1
          b(7,4)=b(7,3)*sb2
          b(5,6)=b(4,6)*sb1
          b(6,5)=b(6,4)*sb2
!-----------------------------------------------------------------------
  700     do 720 np=1,mpe
            n2e=2*np
            do 710 nv=1,n2e
              n2=nv-np
              nn2=abs(n2)
              nn1=np-nn2
              re1=real(nn1,fPrec)*qxt+real(n2,fPrec)*qzt                             !hr06
              ip(np,nv)=int(re1+half)+ipc
              if(-one*re1.gt.pieni) ip(np,nv)=-int(abs(re1)+half)-ipc
!--RE=DISTANCE FROM THE RESONANCE
              re(np,nv)=re1-real(ip(np,nv),fPrec)
              res=re(np,nv)/radi
          chy(np,nv)=cos_mb((real(nn1,fPrec)*pie*phi(1)+real(n2,fPrec)*pie*phi(2))-res*etl) !hr06
          shy(np,nv)=sin_mb((real(nn1,fPrec)*pie*phi(1)+real(n2,fPrec)*pie*phi(2))-res*etl) !hr06
  710       continue
  720     continue
          do 780 np=nta,mpe
            np2=np
            nkk=0
  730       nkk=nkk+1
            n2e=2*np2
            do 770 i=1,nkk
              do 760 nv=1,n2e
                nn2=abs(nv-np2)
                nv1=np2-nn2+(i-1)*2+1
                nv2=np-nv1+2
                rn2=real(nn2,fPrec)*half                                            !hr06
!--EVENESS OF N2
                mm=0
                gerad=rn2-aint(rn2)
                if(abs(gerad).le.pieni) mm=1
!--MM=0 =>N2 UNEVEN, MM=1 => N2 EVEN
                if (mm.eq.0) goto 740
                btc=ab1(np)*b(nv1,nv2)*chy(np2,nv)
                bts=ab1(np)*b(nv1,nv2)*shy(np2,nv)
                goto 750
  740           btc=ab2(np)*b(nv1,nv2)*chy(np2,nv)
                bts=ab2(np)*b(nv1,nv2)*shy(np2,nv)
  750           rtc(np2,nv,np,i)=rtc(np2,nv,np,i)+btc
                rts(np2,nv,np,i)=rts(np2,nv,np,i)+bts
  760         continue
  770       continue
            np2=np2-2
            if(np2.ge.1) goto 730
  780     continue
          nr=nr+1
  790   continue
        nnf(1)=1
        nnf(2)=1
        nnf(3)=2
        nz2(2)=2
        sea=sqrt(ep(1))
        seb=sqrt(ep(2))
        ea=ep(1)
        eb=ep(2)
        e(3,1)=one/eb
        e(1,3)=one/ea
        e(2,2)=(one/seb)/sea                                             !hr06
        nnf(4)=6
        nz2(3)=4
        e(4,1)=sea/eb
        e(1,4)=seb/ea
        e(3,2)=one/seb
        e(2,3)=one/sea
        nnf(5)=24
        nz2(4)=8
        e(5,1)=ea/eb
        e(1,5)=eb/ea
        e(4,2)=sea/seb
        e(2,4)=seb/sea
        e(3,3)=one
        nnf(6)=120
        nz2(5)=16
        e(6,1)=e(5,1)*sea
        e(1,6)=e(1,5)*seb
        e(5,2)=ea/seb
        e(2,5)=eb/sea
        e(4,3)=sea
        e(3,4)=seb
        nnf(7)=720
        nz2(6)=32
        e(7,1)=e(6,1)*sea
        e(1,7)=e(1,6)*seb
        e(6,2)=e(5,2)*sea
        e(2,6)=e(2,5)*seb
        e(5,3)=ea
        e(3,5)=eb
        e(4,4)=sea*seb
        nnf(8)=5040
        nz2(7)=64
        e(8,1)=e(7,1)*sea
        e(1,8)=e(1,7)*seb
        e(7,2)=e(6,2)*sea
        e(2,7)=e(2,6)*seb
        e(6,3)=ea*sea
        e(3,6)=eb*seb
        e(5,4)=ea*seb
        e(4,5)=sea*eb
        nnf(9)=40320
        nz2(8)=128
        e(9,1)=e(8,1)*sea
        e(1,9)=e(1,8)*seb
        e(8,2)=e(7,2)*sea
        e(2,8)=e(2,7)*seb
        e(7,3)=ea**2                                                     !hr06
        e(3,7)=eb**2                                                     !hr06
        e(6,4)=e(5,4)*sea
        e(4,6)=e(4,5)*seb
        e(5,5)=ea*eb
        nnf(10)=362880
        nz2(9)=256
        e(10,1)=e(9,1)*sea
        e(1,10)=e(1,9)*seb
        e(9,2)=e(8,2)*sea
        e(2,9)=e(2,8)*seb
        e(8,3)=e(7,3)*sea
        e(3,8)=e(3,7)*seb
        e(7,4)=e(6,4)*sea
        e(4,7)=e(4,6)*seb
        e(6,5)=e(5,5)*sea
        e(5,6)=e(5,5)*seb
        write(lout,10000)
        write(lout,10030)
        write(lout,10010)nr,'END     ',etl,zero,(beta(l),alfa(l),phi(l),di0(l),dip0(l),clo0(l),clop0(l),l=1,2)
        write(lout,10030)
        write(lout,10110) etl,qwc(1),qwc(2)
        write(lout,10030)
        do 800 iv=2,5
          gtu1=gtu1+dtu(1,iv)
          gtu2=gtu2+dtu(2,iv)
  800   continue
        write(lout,10150) dtu(1,2),dtu(1,3),dtu(1,4),dtu(1,5),gtu1, dtu (2,2),dtu(2,3),dtu(2,4),dtu(2,5),gtu2

        do i=1,2
          do j=1,5
            do l=0,4
              do k=0,4
                if(i.eq.2.and.j.eq.1.and.k.eq.1.and.l.eq.1) write (lout,10160)
                if(abs(dtup(i,j,k,l)).gt.pieni) write(lout,'(10X,G17.10,3X,I2,2X,I2)') dtup(i,j,k,l),k,l
              end do
            end do
          end do
        end do

        write(lout,10060)
        write(lout,10030)
        do 880 np=nta,nte
          write(lout,10080) np
          write(lout,10030)
          vdt1=real(nnf(np),fPrec)/(real(nz2(np),fPrec)*pi)               !hr06
          np2=np
          nkk=0
          write(lout,10090) np
          goto 830
  820     write(lout,10100) np,np2
  830     nkk=nkk+1
          n2e=2*np2
          do 850 i=1,nkk
            do 840 nv=1,n2e
              n2=nv-np2
              nn2=abs(n2)
              nn1=np2-nn2
              nv1=(nn1+(i-1)*2)+1                                        !hr06
              nv2=(np-nv1)+2                                             !hr06
              nv11=nv1-1
              nv21=nv2-1
              nf1=nn1+i
              nf3=nkk-i+1
              nf4=nf3+nn2
              vdt2=vdt1*e(nv1,nv2)/real(nnf(nf1)*nnf(i)*nnf(nf3)*nnf(nf4),fPrec) !hr06
              vdt3=real(nn2,fPrec)*ea+real(nn1,fPrec)*eb                 !hr06
              vdt4=vdt3
              if(n2.ge.0) vdt3=real(n2*nv21,fPrec)*ea + real(nn1*nv11,fPrec)*eb  !hr06
              rtc(np2,nv,np,i)=rtc(np2,nv,np,i)*vdt2*vdt3
              rts(np2,nv,np,i)=rts(np2,nv,np,i)*vdt2*vdt3
  840       continue
  850     continue
          do 870 nv=1,n2e
            mis=1
            rc=zero
            rs=zero
            do 860 i=1,nkk
              rc=rc+real(mis,fPrec)*rtc(np2,nv,np,i)                     !hr06
              rs=rs+real(mis,fPrec)*rts(np2,nv,np,i)                     !hr06
              mis=-mis
  860       continue
            sdel2=sqrt(rc**2+rs**2)                                      !hr06
            n22=nv-np2
            write(lout,10140) n22,ip(np2,nv),ipc,rc,rs,re(np2,nv),sdel2
  870     continue
          np2=np2-2
          if(np2.ge.1) goto 820
  880   continue
        ntx=nte-2
        write(lout,10130)
        do 930 np=1,nte
          write(lout,10090) np
          n2e=2*np
          do 920 nv=1,n2e
            n2=nv-np
            nkk=2
            nph=np+2
            min1=-1
  890       min2=min1
            do 900 i=1,nkk
             rtc(np,nv,np,1)=rtc(np,nv,np,1)+real(min2,fPrec)*rtc(np,nv,nph,i) !hr06
             rts(np,nv,np,1)=rts(np,nv,np,1)+real(min2,fPrec)*rts(np,nv,nph,i) !hr06
              min2=-min2
  900       continue
            nph=nph+2
            if(nph.gt.nte) goto 910
            nkk=nkk+1
            min1=-min1
            goto 890
  910       cc_r=rtc(np,nv,np,1)
            ss=rts(np,nv,np,1)
            sdel=sqrt(cc_r**2+ss**2)                                       !hr06
            write(lout,10140) n2,ip(np,nv),ipc,cc_r,ss,re(np,nv),sdel
  920     continue
  930   continue
  940 continue
      call clorb(ded)
      do 950 l=1,2
        clo0(l)=clo(l)
        clop0(l)=clop(l)
  950 continue
      call clorb(zero)
      do 960 l=1,2
        ll=2*l
        di0(l)=(clo0(l)-clo(l))/ded
        dip0(l)=(clop0(l)-clop(l))/ded
  960 continue
!-----------------------------------------------------------------------
      return
10000 format(1x,i4,27x,f7.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f6.3,1x,f6.2,1x, &
     &f6.3,1x,f7.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f6.3,1x,f6.2,1x,f6.3)
10020 format(t5,'---- ENTRY SUBRES ----')
10030 format(131('-'))
10040 format('  NR  TYP      L-TOTAL  LENGTH   BETAH  ALFAH  ',         &
     &' PHIH   DISH  DISPH   CLOH  CLOPH',                              &
     &'   BETAV  ALFAV   PHIV   DISV  DISPV   CLOV  CLOPV'/ 1x,         &
     &'                 (M)      (M)     (M)           ',               &
     &'(QE)   (M)   (RAD)   (MM) (MRAD)',                               &
     &'    (M)           (QE)   (M)   (RAD)   (MM) (MRAD)')
10050 format(//7x,'INIT. X-AMPLITUDE=',g15.8,'X-EMITTANCE=',g15.8,/40x, &
     &/7x,'INIT. Y-AMPLITUDE=',g15.8,'Y-EMITTANCE=',g15.8,              &
     &'UNITS ARE (PI X MM X MRAD)'//)
10060 format(//10x,'E=NX*QX+NY*QY-P',//10x,'CLOSESET P-VALUE CHANGED ', &
     &'BY D-P',//10x,'DELTA-E STANDS FOR THE RESONANCE-WIDTH' //10x)
10070 format(/10x,'RELATIVE ENERGY DEVIATION  ',t40,f23.16/ 10x,        &
     &'TUNES -HORIZONTAL',t40,f23.16/ 10x,'      -VERTICAL  ',t40,f23.16)
10080 format(/10x,'RESONANCE EXCITING MULTIPOLE-ORDER = ',i2)
10090 format(//20x,'RESONANCE-ORDER =',i2/20x,100('-')/ 20x,'| NY |',   &
     &'   P  | D-P |',2x,'DRIVING-COS ',3x,'|', 2x,'DRIVING-SIN ',3x,'|'&
     &, 8x,'E',8x,'|',5x,'DELTA-E',5x,'|')
10100 format(//20x,'RESONANCE-ORDER =',i2,5x,'SUBRESONANCE-ORDER = ',i2,&
     &/20x,100('-')/ 20x,'| NY |','   P  | D-P |',2x,'DRIVING-COS ',3x, &
     &'|', 2x,'DRIVING-SIN ',3x,'|', 8x,'E',8x,'|',5x,'DELTA-E',5x,'|')
10110 format(/10x,'PRECISE LENGTH OF THE MACHINE : ',f43.33/ /10x,      &
     &'   PRECISE HORIZONTAL Q-VALUE : ',f43.33/ /10x,                  &
     &'     PRECISE VERTICAL Q-VALUE : ',f43.33/)
10120 format(t8,'  PLANE     DISP(MM)     DISP(MRAD)   '/ t6,'      X  '&
     &,2(f12.3,3x)/t10,'  Y  ',2(f12.3,3x)/)
10130 format(//10x,'E=NX*QX+NY*QY-P',//10x,'CLOSESET P-VALUE CHANGED ', &
     &'BY D-P',//10x,'DELTA-E STANDS FOR THE RESONANCE-WIDTH' //10x,    &
     &'!!!! ALL SUBRESONANCES ARE INCLUDED !!!! ')
10140 format(20x,'| ',i2,' | ',i4,' | ',i3,' |', g16.8,' |',g16.8,' |', &
     &g16.8,' |',g16.8,' |')
10150 format(/10x,'NONLINEAR DETUNING  '// 10x,'CHANGE IN QX'/ 10x,     &
     &' 4. ORDER ',f15.12/ 10x,' 6. ORDER ',f15.12/ 10x,' 8. ORDER ',f15&
     &.12/ 10x,'10. ORDER ',f15.12/ 10x,'   TOTAL  ',f15.12/ 10x,       &
     &'CHANGE IN QY'/ 10x,' 4. ORDER ',f15.12/ 10x,' 6. ORDER ',f15.12/ &
     &10x,' 8. ORDER ',f15.12/ 10x,'10. ORDER ',f15.12/ 10x,'   TOTAL  '&
     &,f15.12// 10x,'DETUNING ORDER BY ORDER'// 10x,                    &
     &'Qx - COEFFICIENT   Ex  EY'/ 10x,'-------------------------')
10160 format(/ 10x,'Qy - COEFFICIENT   Ex  Ey'/ 10x,                    &
     &'-------------------------')
10010 format(1x,i4,1x,a8,1x,f8.2,1x,f7.3,1x, f7.2,1x,f6.2,1x,f6.2,1x,f6.&
     &2,1x,f6.3,1x,f6.2,1x,f6.3,1x, f7.2,1x,f6.2,1x,f6.2,1x,f6.2,1x,f6. &
     &3,1x,f6.2,1x,f6.3)
end subroutine subre

      subroutine detune(iv,ekk,ep,beta,dtu,dtup,dfac)
!-----------------------------------------------------------------------
!  USED FOR SUBRE - CALCULATES DETUNING
!-----------------------------------------------------------------------
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      use crcoall
      use parpro
      implicit none
      integer iv,iv2,iv3,iv4,iv5,iv6
      real(kind=fPrec) beta,dfac,dtu,dtu1,dtu2,dtup,ekk,ep,vor,vtu1,vtu2
      dimension dfac(10),dtu(2,5),ep(2),beta(2),dtup(2,5,0:4,0:4)
      save
!-----------------------------------------------------------------------
      if(iv.lt.2) then
        write(lout,*)
        write(lout,*) '       ***** ERROR IN DETUNE *****'
        write(lout,*)
        write(lout,*) '       IV LESS THAN 2, NO DETUNING POSSIBLE'
        write(lout,*)
        return
      endif
      iv2=2*iv
      iv3=iv+1
!      vtu1=(((-one*ekk)*(half**iv2))*dfac(iv2))/pi                       !hr06
      vtu1=(((-one*ekk)*exp_mb(real(iv2,fPrec)*log_mb(half)))*          &
     &dfac(iv2))/pi    !hr13
      dtu1=zero
      dtu2=zero
      do 10 iv4=1,iv3
        iv5=iv4-1
        iv6=iv-iv5
        vor=one
        if(mod(iv6,2).ne.0) vor=-one                                 !hr06
!        vtu2=vor/(dfac(iv5+1)**2)/(dfac(iv6+1)**2)*(beta(1)**iv5)* (beta&
!     &(2)**iv6)
        vtu2=(((vor/(dfac(iv5+1)**2))/(dfac(iv6+1)**2))*                &!hr13
     &exp_mb(real(iv5,fPrec)*log_mb(beta(1))))*                         &!hr13
     &exp_mb(real(iv6,fPrec)*log_mb(beta(2)))                            !hr13
        if(iv5.ne.0) then
!          dtu1=dtu1+((vtu2*dble(iv5))*(ep(1)**(iv5-1)))*(ep(2)**iv6)    !hr06
         dtu1=dtu1+((vtu2*real(iv5,fPrec))*exp_mb(real(iv5-1,fPrec)*    &
     &log_mb(ep(1))))*exp_mb(real(iv6,fPrec)*log_mb(ep(2)))                           !hr13
         dtup(1,iv,iv5-1,iv6)=dtup(1,iv,iv5-1,iv6)+                     &
     &(vtu2*real(iv5,fPrec))*vtu1 !hr06
        endif
        if(iv6.ne.0) then
!          dtu2=dtu2+((vtu2*dble(iv6))*(ep(1)**iv5))*(ep(2)**(iv6-1))     !hr06
          dtu2=dtu2+((vtu2*real(iv6,fPrec))*exp_mb(real(iv5,fPrec)*     &
     &log_mb(ep(1))))*exp_mb(real(iv6-1,fPrec)*log_mb(ep(2)))                  !hr13
         dtup(2,iv,iv5,iv6-1)=dtup(2,iv,iv5,iv6-1)+                     &
     &(vtu2*real(iv6,fPrec))*vtu1 !hr06
        endif
   10 continue
      dtu(1,iv)=dtu(1,iv)+vtu1*dtu1
      dtu(2,iv)=dtu(2,iv)+vtu1*dtu2
      return
end

!-----------------------------------------------------------------------
!  CALCULATION OF DRIVINGTERMS OF RESONANCES INCLUDING SUBRESONANCE
!  USED FOR SEARCH
!-----------------------------------------------------------------------
subroutine subsea(dpp)
  ! Rewritten to remove computed gotos by V.K.B.Olsen on 23/11/2017
  use floatPrecision
  use numerical_constants
  use mathlib_bouncer
  use parpro
  use crcoall
  use mod_common
  use mod_commons
  use mod_commont
  implicit none
  integer i,ii,ik,im,ip,ium,ix,izu,j,jj,jk,jm,k,k1,kpz,kzz,l,l1,l2,ll,lmin,mm,mpe,mx,n2,n2e,nf1,nf3,&
          nf4,nkk,nmz,nn1,nn2,nnf,np,np2,ns,nv,nv1,nv11,nv2,nv21,nz2,dj
  real(kind=fPrec) aa,ab1,ab2,alfa,b,b1,b2,bb,benkr,beta,btc,bts,chy,ci,cikve,cr,crkve,cxzi,cxzr,&
          cxzyi,cxzyr,cxzyrr,del,dphi,dpp,dppi,dpr,dt,dyy1,dyy2,e,ea,eb,ekk,ep,etl,gerad,phi,phibf,&
          phy,pie,puf,qu,qv,qw,r0,r0a,radi,re,re1,res,rn2,sb1,sb2,sea,seb,shy,t,vdt1,vdt2,vdt3,xl,xs,zl,zs,quz,qvz
#ifdef TILT
  real(kind=fPrec) dyy11,qu1,tiltck,tiltck1,tiltck2,tiltck3,tiltck4,tiltck5,tiltckuk,tiltsk,tiltsk1,tiltsk2,tiltsk3,tiltsk4,tiltsk5
#endif
      dimension t(5,4)
      dimension beta(2),alfa(2),phi(2),phibf(2)
      dimension aa(mmul),bb(mmul)
      dimension qw(2),dpr(5)
      dimension nnf(10),ep(2)
      dimension ab1(10),ab2(10),re(10,18),ip(10,18)
      dimension b(10,10),nz2(9),e(10,10)
      dimension chy(9,18),shy(9,18)
      dimension cr(mmul),ci(mmul)
      save
!-----------------------------------------------------------------------
      ium=5
      do 10 i=1,ium
        dpr(i)=zero
   10 continue
      do i=1,ium
        do j=1,4
          t(i,j)=zero
        end do
      end do

      do 30 i=1,2
        beta(i)=zero
        alfa(i)=zero
        phi(i)=zero
        phibf(i)=zero
        qw(i)=zero
        ep(i)=zero
   30 continue

      do i=1,10
        nnf(i)=0
        do j=1,18
          ip(i,j)=0
          re(i,j)=zero
        end do
      end do

      do 50 i=1,mmul
        aa(i)=zero
        bb(i)=zero
        cr(i)=zero
        ci(i)=zero
   50 continue

      do 100 i=1,9
        nz2(i)=0
        do 90 j=1,18
          chy(i,j)=zero
          shy(i,j)=zero
          do 80 k=1,10
            do 60 ii=1,10
              e(k,ii)=zero
              b(k,ii)=zero
   60       continue
            do 70 l=1,5
              rtc(i,j,k,l)=zero
              rts(i,j,k,l)=zero
   70       continue
   80     continue
   90   continue
  100 continue

      btc=zero
      bts=zero
      phy=zero
      dt=zero
      del=zero
      ns=0
      ik=0
      pie=two*pi
      etl=zero
      radi=totl/pie
      dpr(1)=dpp*c1e3
      call clorb2(dpp)
      call betalf(dpp,qw)
      if(ierro /= 0) then
        write(lout,"(a)") "SUBSEA> ERROR No optical solution."
        call prror(-1)
      end if
      call envar(dpp)

!--STARTVALUES OF THE TRAJECTORIES
      do l=1,2
        ll=2*l
        alfa(l)=alf0(l)
        beta(l)=bet0(l)
        t(1,ll-1)=clo(l)
        t(1,ll)=clop(l)
      end do

      do i=1,4
        do j=1,4
          t(i+1,j)=ta(j,i)
          t(i+1,j)=ta(j,i)
        end do
      end do

!--EP=EMITTANCE IN PI*MM*MRAD
      ep(1)=tam1**2/beta(1)                                              !hr06
      ep(2)=tam2**2/beta(2)                                              !hr06

!--SINGLE TURN BLOCKLOOP
      izu=0
      do 740 k=1,iu
        do k1=1,10
          ab1(k1)=zero
          ab2(k1)=zero
        end do

        ix=ic(k)
        if(ix.gt.nblo) goto 210
        jj=0
        dj=1
        if(ix.gt.0) goto 140
        ix=-ix
        jj=mel(ix)+1
        dj=-1
  140   jm=mel(ix)
!--BLOCKELEMENTLOOP
        do 200 j=1,jm
          jj=jj+dj
          jk=mtyp(ix,jj)
          if(ithick.eq.1.and.kz(jk).ne.0) goto 170
          if(ithick.eq.0.and.kz(jk).ne.0) goto 740

!--PURE DRIFTLENGTH
          etl=etl+el(jk)
          do l=1,2
            ll=2*l

            if(abs(t(ll,ll-1)).gt.pieni) then
              phibf(l)=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
            else
              phibf(l)=pi2
            endif

            do i=1,ium
              t(i,ll-1)=t(i,ll-1)+t(i,ll)*(el(jk))
            end do
          end do

          do l=1,2
            ll=2*l
            beta(l)=t(ll,ll-1)**2+t(ll+1,ll-1)**2                        !hr06
            alfa(l)=-one*(t(ll,ll-1)*t(ll,ll)+t(ll+1,ll-1)*t(ll+1,ll))   !hr06

            if(abs(t(ll,ll-1)).gt.pieni) then
              dphi=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
            else
              dphi=pi2-phibf(l)
            endif

            if(-one*dphi.gt.pieni) dphi=dphi+pi                          !hr06
            phi(l)=phi(l)+dphi
          end do

          goto 200

!--MAGNETELEMENT
  170     continue
          if(kz(jk).ne.8) etl=etl+el(jk)
          do l=1,2
            ll=2*l
            if(abs(t(ll,ll-1)).gt.pieni) then
              phibf(l)=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
            else
              phibf(l)=zero
            endif
            do i=1,ium
              puf=t(i,ll-1)
            t(i,ll-1)=(puf*a(jk,l,1)+t(i,ll)*a(jk,l,2))+dpr(i)*a(jk,l,5) !hr06
            t(i,ll)=(puf*a(jk,l,3)+t(i,ll)*a(jk,l,4))+dpr(i)*a(jk,l,6)   !hr06
            enddo
          enddo
          do l=1,2
            ll=2*l
            beta(l)=t(ll,ll-1)**2+t(ll+1,ll-1)**2                        !hr06
            alfa(l)=-one*(t(ll,ll-1)*t(ll,ll)+t(ll+1,ll-1)*t(ll+1,ll))   !hr06
            if(abs(t(ll,ll-1)).gt.pieni) then
              dphi=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
            else
              dphi=-phibf(l)
            endif
            if(kz(jk).ne.8.and.-one*dphi.gt.pieni) dphi=dphi+pi          !hr06
            phi(l)=phi(l)+dphi
          enddo
  200   continue
        goto 740
!--NL-INSERTION
  210   ix=ix-nblo
        qu=zero
        qv=zero
        kpz=kp(ix)
        if(kpz.eq.6) goto 740
        kzz=kz(ix)
        if(kzz == 22) then
          do l=1,2
            ll=2*l
            if(abs(t(ll,ll-1)).gt.pieni) then
              phibf(l)=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))
            else
              phibf(l)=zero
            end if
            do i=1,ium
              puf=t(i,ll-1)
              t(i,ll-1)=(puf*rrtr(imtr(ix),ll-1,ll-1)+t(i,ll)*rrtr(imtr(ix),ll-1,ll))+dpr(i)*rrtr(imtr(ix),ll-1,6)
              t(i,ll)=(puf*rrtr(imtr(ix),ll,ll-1)+t(i,ll)*rrtr(imtr(ix),ll,ll))+dpr(i)*rrtr(imtr(ix),ll,6)
            end do
            t(1,ll-1)=t(1,ll-1)+cotr(imtr(ix),ll-1)
            t(1,ll)=t(1,ll)+cotr(imtr(ix),ll)
            beta(l)=t(ll,ll-1)**2+t(ll+1,ll-1)**2
            alfa(l)=-one*(t(ll,ll-1)*t(ll,ll)+t(ll+1,ll-1)*t(ll+1,ll))
            if(abs(t(ll,ll-1)) > pieni) then
              dphi=atan_mb(t(ll+1,ll-1)/t(ll,ll-1))-phibf(l)
            else
              dphi=-one*phibf(l)
            end if
            if(-one*dphi.gt.pieni) dphi=dphi+pi
                        phi(l)=phi(l)+dphi
          enddo
        endif
        if(kzz.eq.0.or.kzz.eq.20.or.kzz.eq.22) goto 740
        if(kzz.eq.15) goto 740
! JBG RF CC Multipoles to 740
        if(kzz.eq.26.or.kzz.eq.27.or.kzz.eq.28) goto 740
        if(kzz.eq.-26.or.kzz.eq.-27.or.kzz.eq.-28) goto 740
        dyy1=zero
        dyy2=zero
        if(iorg.lt.0) mzu(k)=izu
        izu=mzu(k)+1
        ekk=(sm(ix)+zfz(izu)*ek(ix))/(one+dpp)
        izu=izu+1
        xs=xpl(ix)+zfz(izu)*xrms(ix)
        izu=izu+1
        zs=zpl(ix)+zfz(izu)*zrms(ix)
#include "include/alignl.f90"

      select case (kzz)
      case (1) ! HORIZONTAL DIPOLE
        ekk=ekk*c1e3
#include "include/kicka01h.f90"
      case (2) ! NORMAL QUADRUPOLE
#include "include/kicka02h.f90"
      case (3) ! NORMAL SEXTUPOLE
        ekk=ekk*c1m3
#include "include/kicka03h.f90"
      case (4) ! NORMAL OCTUPOLE
        ekk=ekk*c1m6
#include "include/kicka04h.f90"
      case (5) ! NORMAL DECAPOLE
        ekk=ekk*c1m9
#include "include/kicka05h.f90"
      case (6) ! NORMAL DODECAPOLE
        ekk=ekk*c1m12
#include "include/kicka06h.f90"
      case (7) ! NORMAL 14-POLE
        ekk=ekk*c1m15
#include "include/kicka07h.f90"
      case (8) ! NORMAL 16-POLE
        ekk=ekk*c1m18
#include "include/kicka08h.f90"
      case (9) ! NORMAL 18-POLE
        ekk=ekk*c1m21
#include "include/kicka09h.f90"
      case (10) ! NORMAL 20-POLE
        ekk=ekk*c1m24
#include "include/kicka10h.f90"
      case (11)
        r0=ek(ix)
        if(abs(dki(ix,1)).gt.pieni) then
          if(abs(dki(ix,3)).gt.pieni) then
#include "include/multl01.f90"
            do i=2,ium
#include "include/multl02.f90"
            end do
          else
#include "include/multl03.f90"
          end if
        end if
        if(abs(dki(ix,2)).gt.pieni) then
          if(abs(dki(ix,3)).gt.pieni) then
#include "include/multl04.f90"
            do i=2,ium
#include "include/multl05.f90"
            end do
          else
#include "include/multl06.f90"
          end if
        end if
        mpe=9
        mx=0
        if(abs(r0).le.pieni) goto 740
        nmz=nmu(ix)
        if(nmz.eq.0) then
          izu=izu+2*mmul
          goto 740
        end if
        im=irm(ix)
        r0a=one
        benkr=ed(ix)/(one+dpp)
        cr(1)=one
        cr(2)=xl
        ci(2)=zl
        cxzyr=xl
        cxzyi=zl
        cxzr=cxzyr
        cxzi=cxzyi
        dyy1=zero
        dyy2=zero
        qu=zero
        qv=zero
        lmin=3
        if(nmz.eq.1) lmin=2

        do l=lmin,mmul
          cr(l)=zero
          ci(l)=zero
        end do

        do l=1,nmz
#include "include/multl13.f90"
        end do
#ifdef TILT
#include "include/multl07e.f90"
#endif
        izu=(izu+2*mmul)-2*nmz                                           !hr06
      case (12,13,14,15,16,17,18,19,20,21,22,23)
        goto 740
      case (24) ! DIPEDGE ELEMENT
#include "include/kickadpe.f90"
      case (25) ! Solenoid
#include "include/kickaso1.f90"
      case (26,27,28)
        goto 740

        !-----------------
        !--SKEW ELEMENTS--
        !------------------
      case (-1) ! VERTICAL DIPOLE
        ekk=ekk*c1e3
#include "include/kicka01v.f90"
      case (-2) ! SKEW QUADRUPOLE
#include "include/kicka02v.f90"
      case (-3) ! SKEW SEXTUPOLE
        ekk=ekk*c1m3
#include "include/kicka03v.f90"
      case (-4) ! SKEW OCTUPOLE
        ekk=ekk*c1m6
#include "include/kicka04v.f90"
      case (-5) ! SKEW DECAPOLE
        ekk=ekk*c1m9
#include "include/kicka05v.f90"
      case (-6) ! SKEW DODECAPOLE
        ekk=ekk*c1m12
#include "include/kicka06v.f90"
      case (-7) ! SKEW 14-POLE
        ekk=ekk*c1m15
#include "include/kicka07v.f90"
      case (-8) ! SKEW 16-POLE
        ekk=ekk*c1m18
#include "include/kicka08v.f90"
      case (-9) ! SKEW 18-POLE
        ekk=ekk*c1m21
#include "include/kicka09v.f90"
      case (-10) ! SKEW 20-POLE
        ekk=ekk*c1m24
#include "include/kicka10v.f90"

      case default
        goto 740
      end select

      t(1,2)=t(1,2)+dyy1
      t(1,4)=t(1,4)+dyy2
      do i=2,ium
        if(kzz.eq.24) then
          t(i,2)=(t(i,2)+t(i,1)*qu)-qv*t(i,3)                          !hr06
          t(i,4)=(t(i,4)-t(i,3)*quz)-qvz*t(i,1)                        !hr06
#include "include/phas1so1.f90"
#include "include/phas2so1.f90"
#include "include/phas3so1.f90"
        else
          t(i,2)=(t(i,2)+t(i,1)*qu)-qv*t(i,3)                          !hr06
          t(i,4)=(t(i,4)-t(i,3)*qu)-qv*t(i,1)                          !hr06
        end if
      end do
      do l=1,2
        ll=2*l
        alfa(l)=-one*(t(ll,ll-1)*t(ll,ll)+t(ll+1,ll-1)*t(ll+1,ll))     !hr06
      end do
      if(mpe.gt.9.or.(mpe.eq.9.and.nmz.le.1)) goto 740
      if(mpe.lt.nta) goto 740
      if(mpe.gt.nte) mpe=nte
      if(nta.gt.2) goto 500
      if(mx.eq.-1.or.mx.eq.1.or.mx.eq.2.or.mx.eq.3.or.mx.eq.4 .or.mx.eq.5.or.mx.eq.6.or.mx.eq.7) goto 500

!-----------------------------------------------------------------------
!  SKEW-QUADRUPOLE;MULTIPOLES UP TO 9-TH ORDER
!-----------------------------------------------------------------------
        do l=2,nmz
          l1=l-1
          ab2(2)=ab2(2)+real(l1,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))               !hr06
        end do

  500   b1=beta(1)
        b2=beta(2)
        sb1=sqrt(b1)
        sb2=sqrt(b2)
        b(3,1)=b1
        b(1,3)=b2
        b(2,2)=sb1*sb2
        if(nta.gt.3) goto 520
        if(mpe.eq.2.or.(mpe.eq.9.and.nmz.le.2)) goto 650
        if(mx.eq.1.or.mx.eq.2.or.mx.eq.3.or.mx.eq.4 .or.mx.eq.5.or.mx.eq&
     &.6.or.mx.eq.7) goto 520

!-----------------------------------------------------------------------
!  REGULAR-SEXTUPOLE;MULTIPOLES UP TO 9-TH ORDER
!-----------------------------------------------------------------------
        l2=1
        do l=3,nmz
          l1=l-2
          ab1(3)=ab1(3)+real(l2,fPrec)*(aa(l)*ci(l1)+bb(l)*cr(l1))             !hr06
          ab2(3)=ab2(3)+real(l2,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))             !hr06
          l2=l2*l/l1
        end do

  520   b(4,1)=b1*sb1
        b(1,4)=b2*sb2
        b(3,2)=b1*sb2
        b(2,3)=b2*sb1
        if(nta.gt.4) goto 540
        if(mpe.eq.3.or.(mpe.eq.9.and.nmz.le.3)) goto 650
        if(mx.eq.2.or.mx.eq.3.or.mx.eq.4 .or.mx.eq.5.or.mx.eq.6.or.mx.eq&
     &.7) goto 540

!-----------------------------------------------------------------------
!  REGULAR-OCTUPOLE;MULTIPOLES UP TO 9-TH ORDER
!-----------------------------------------------------------------------
        l2=1
        do l=4,nmz
          l1=l-3
          ab1(4)=ab1(4)+real(l2,fPrec)*(aa(l)*ci(l1)+bb(l)*cr(l1))             !hr06
          ab2(4)=ab2(4)+real(l2,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))             !hr06
          l2=l2*l/l1
        end do

  540   b(5,1)=b1**2                                                     !hr06
        b(1,5)=b2**2                                                     !hr06
        b(4,2)=b(3,2)*sb1
        b(2,4)=b(2,3)*sb2
        b(3,3)=b1*b2
        if(nta.gt.5) goto 560
        if(mpe.eq.4.or.(mpe.eq.9.and.nmz.le.4)) goto 650
        if(mx.eq.3.or.mx.eq.4 .or.mx.eq.5.or.mx.eq.6.or.mx.eq.7)        &
     &goto 560

!-----------------------------------------------------------------------
!  REGULAR-DEKAPOLE;MULTIPOLES UP TO 9-TH ORDER
!-----------------------------------------------------------------------
        l2=1
        do l=5,nmz
          l1=l-4
          ab1(5)=ab1(5)+real(l2,fPrec)*(aa(l)*ci(l1)+bb(l)*cr(l1))             !hr06
          ab2(5)=ab2(5)+real(l2,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))             !hr06
          l2=l2*l/l1
        end do

  560   b(6,1)=b(5,1)*sb1
        b(1,6)=b(1,5)*sb2
        b(5,2)=b(4,2)*sb1
        b(2,5)=b(2,4)*sb2
        b(4,3)=b(4,2)*sb2
        b(3,4)=b(2,4)*sb1
        if(nta.gt.6) goto 580
        if(mpe.eq.5.or.(mpe.eq.9.and.nmz.le.5)) goto 650
        if(mx.eq.4 .or.mx.eq.5.or.mx.eq.6.or.mx.eq.7) goto 580

!-----------------------------------------------------------------------
!  REGULAR-12-POLE;MULTIPOLES UP TO 9-TH ORDER
!-----------------------------------------------------------------------
        l2=1
        do l=6,nmz
          l1=l-5
          ab1(6)=ab1(6)+real(l2,fPrec)*(aa(l)*ci(l1)+bb(l)*cr(l1))             !hr06
          ab2(6)=ab2(6)+real(l2,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))             !hr06
          l2=l2*l/l1
        end do

  580   b(7,1)=b(6,1)*sb1
        b(1,7)=b(1,6)*sb2
        b(6,2)=b(5,2)*sb1
        b(2,6)=b(2,5)*sb2
        b(5,3)=b(5,2)*sb2
        b(3,5)=b(2,5)*sb1
        b(4,4)=b(3,4)*sb1
        if(nta.gt.7) goto 600
        if(mpe.eq.6.or.(mpe.eq.9.and.nmz.le.6)) goto 650
        if(mx.eq.5.or.mx.eq.6.or.mx.eq.7) goto 600

!-----------------------------------------------------------------------
!  REGULAR-14-POLE;MULTIPOLES UP TO 9-TH ORDER
!-----------------------------------------------------------------------
        l2=1
        do l=7,nmz
          l1=l-6
          ab1(7)=ab1(7)+real(l2,fPrec)*(aa(l)*ci(l1)+bb(l)*cr(l1))             !hr06
          ab2(7)=ab2(7)+real(l2,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))             !hr06
          l2=l2*l/l1
        end do

  600   b(8,1)=b(7,1)*sb1
        b(1,8)=b(1,7)*sb2
        b(7,2)=b(7,1)*sb2
        b(2,7)=b(1,7)*sb1
        b(6,3)=b(5,3)*sb1
        b(3,6)=b(3,5)*sb2
        b(5,4)=b(4,4)*sb1
        b(4,5)=b(4,4)*sb2
        if(nta.gt.8) goto 620
        if(mpe.eq.7.or.(mpe.eq.9.and.nmz.le.7)) goto 650
        if(mx.eq.6.or.mx.eq.7) goto 620

!-----------------------------------------------------------------------
!  REGULAR-16-POLE;MULTIPOLES UP TO 9-TH ORDER
!-----------------------------------------------------------------------
        l2=1
        do l=8,nmz
          l1=l-7
          ab1(8)=ab1(8)+real(l2,fPrec)*(aa(l)*ci(l1)+bb(l)*cr(l1))             !hr06
          ab2(8)=ab2(8)+real(l2,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))             !hr06
          l2=l2*l/l1
        end do

  620   b(9,1)=b(8,1)*sb1
        b(1,9)=b(1,8)*sb2
        b(8,2)=b(8,1)*sb2
        b(2,8)=b(1,8)*sb1
        b(7,3)=b(7,2)*sb2
        b(3,7)=b(2,7)*sb1
        b(6,4)=b(6,3)*sb2
        b(4,6)=b(3,6)*sb1
        b(5,5)=b(4,5)*sb1
        if(mpe.eq.8.or.(mpe.eq.9.and.nmz.le.8)) goto 650
        if(mx.eq.7) goto 640

!-----------------------------------------------------------------------
!  REGULAR-18-POLE
!-----------------------------------------------------------------------
        l2=1
        do l=9,nmz
          l1=l-8
          ab1(9)=ab1(9)+real(l2,fPrec)*(aa(l)*ci(l1)+bb(l)*cr(l1))             !hr06
          ab2(9)=ab2(9)+real(l2,fPrec)*(aa(l)*cr(l1)-bb(l)*ci(l1))             !hr06
         l2=l2*l/l1
        end do

  640   b(10,1)=b(9,1)*sb1
        b(1,10)=b(1,9)*sb2
        b(9,2)=b(9,1)*sb2
        b(2,9)=b(1,9)*sb1
        b(8,3)=b(8,2)*sb2
        b(3,8)=b(2,8)*sb1
        b(4,7)=b(3,7)*sb1
        b(7,4)=b(7,3)*sb2
        b(5,6)=b(4,6)*sb1
        b(6,5)=b(6,4)*sb2
!-----------------------------------------------------------------------
  650   do 670 np=1,mpe
          n2e=2*np
          do 660 nv=1,n2e
            n2=nv-np
            nn2=abs(n2)
            nn1=np-nn2
            re1=real(nn1,fPrec)*qxt+real(n2,fPrec)*qzt                               !hr06
            ip(np,nv)=int(re1+half)+ipt
            if(-one*re1.gt.pieni) ip(np,nv)=-int(abs(re1)+half)-ipt      !hr06
!--RE=DISTANCE FROM THE RESONANCE
            re(np,nv)=re1-real(ip(np,nv),fPrec)                          !hr06
            res=re(np,nv)/radi
           chy(np,nv)=cos_mb((real(nn1,fPrec)*phi(1)+real(n2,fPrec)*    &
     &phi(2))-res*etl) !hr06
           shy(np,nv)=sin_mb((real(nn1,fPrec)*phi(1)+real(n2,fPrec)*    &
     &phi(2))-res*etl) !hr06
  660     continue
  670   continue
        do 730 np=nta,mpe
          np2=np
          nkk=0
  680     nkk=nkk+1
          n2e=2*np2
          do 720 i=1,nkk
            do 710 nv=1,n2e
              nn2=abs(nv-np2)
              nv1=np2-nn2+(i-1)*2+1
              nv2=np-nv1+2
              rn2=real(nn2,fPrec)*half                                   !hr06
!--EVENESS OF N2
              mm=0
              gerad=rn2-aint(rn2)
              if(abs(gerad).le.pieni) mm=1
!--MM=0 =>N2 UNEVEN, MM=1 => N2 EVEN
              if (mm.eq.0) goto 690
              btc=(ab1(np)*b(nv1,nv2))*chy(np2,nv)                       !hr06
              bts=(ab1(np)*b(nv1,nv2))*shy(np2,nv)                       !hr06
              goto 700
  690         btc=(ab2(np)*b(nv1,nv2))*chy(np2,nv)                       !hr06
              bts=(ab2(np)*b(nv1,nv2))*shy(np2,nv)                       !hr06
  700         rtc(np2,nv,np,i)=rtc(np2,nv,np,i)+btc
              rts(np2,nv,np,i)=rts(np2,nv,np,i)+bts
  710       continue
  720     continue
          np2=np2-2
          if(np2.ge.1) goto 680
  730   continue
  740 continue
      nnf(1)=1
      nnf(2)=1
      nnf(3)=2
      nz2(2)=2
      sea=sqrt(ep(1))
      seb=sqrt(ep(2))
      ea=ep(1)
      eb=ep(2)
      e(3,1)=one/eb
      e(1,3)=one/ea
      e(2,2)=(one/seb)/sea                                               !hr06
      nnf(4)=6
      nz2(3)=4
      e(4,1)=sea/eb
      e(1,4)=seb/ea
      e(3,2)=one/seb
      e(2,3)=one/sea
      nnf(5)=24
      nz2(4)=8
      e(5,1)=ea/eb
      e(1,5)=eb/ea
      e(4,2)=sea/seb
      e(2,4)=seb/sea
      e(3,3)=one
      nnf(6)=120
      nz2(5)=16
      e(6,1)=e(5,1)*sea
      e(1,6)=e(1,5)*seb
      e(5,2)=ea/seb
      e(2,5)=eb/sea
      e(4,3)=sea
      e(3,4)=seb
      nnf(7)=720
      nz2(6)=32
      e(7,1)=e(6,1)*sea
      e(1,7)=e(1,6)*seb
      e(6,2)=e(5,2)*sea
      e(2,6)=e(2,5)*seb
      e(5,3)=ea
      e(3,5)=eb
      e(4,4)=sea*seb
      nnf(8)=5040
      nz2(7)=64
      e(8,1)=e(7,1)*sea
      e(1,8)=e(1,7)*seb
      e(7,2)=e(6,2)*sea
      e(2,7)=e(2,6)*seb
      e(6,3)=ea*sea
      e(3,6)=eb*seb
      e(5,4)=ea*seb
      e(4,5)=sea*eb
      nnf(9)=40320
      nz2(8)=128
      e(9,1)=e(8,1)*sea
      e(1,9)=e(1,8)*seb
      e(8,2)=e(7,2)*sea
      e(2,8)=e(2,7)*seb
      e(7,3)=ea**2                                                       !hr06
      e(3,7)=eb**2                                                       !hr06
      e(6,4)=e(5,4)*sea
      e(4,6)=e(4,5)*seb
      e(5,5)=ea*eb
      nnf(10)=362880
      nz2(9)=256
      e(10,1)=e(9,1)*sea
      e(1,10)=e(1,9)*seb
      e(9,2)=e(8,2)*sea
      e(2,9)=e(2,8)*seb
      e(8,3)=e(7,3)*sea
      e(3,8)=e(3,7)*seb
      e(7,4)=e(6,4)*sea
      e(4,7)=e(4,6)*seb
      e(6,5)=e(5,5)*sea
      e(5,6)=e(5,5)*seb
      do 780 np=nta,nte
        vdt1=real(nnf(np),fPrec)/(real(nz2(np),fPrec)*pi)                            !hr06
        np2=np
        nkk=0
  750   nkk=nkk+1
        n2e=2*np2
        do 770 i=1,nkk
          do 760 nv=1,n2e
            n2=nv-np2
            nn2=abs(n2)
            nn1=np2-nn2
            nv1=nn1+(i-1)*2+1
            nv2=np-nv1+2
            nv11=nv1-1
            nv21=nv2-1
            nf1=nn1+i
            nf3=nkk-i+1
            nf4=nf3+nn2
            vdt2=vdt1*e(nv1,nv2)/                                       &
     &real(nnf(nf1)*nnf(i)*nnf(nf3)*nnf(nf4),fPrec) !hr06
            vdt3=real(nn2,fPrec)*ea+real(nn1,fPrec)*eb                               !hr06
            if(n2.ge.0) then
              vdt3=real(n2*nv21,fPrec)*ea+real(nn1*nv11,fPrec)*eb          !hr06
            end if
            rtc(np2,nv,np,i)=rtc(np2,nv,np,i)*vdt2*vdt3
            rts(np2,nv,np,i)=rts(np2,nv,np,i)*vdt2*vdt3
  760     continue
  770   continue
        np2=np2-2
        if(np2.ge.1) goto 750
  780 continue

  return

end subroutine subsea

subroutine decoup
!-----------------------------------------------------------------------
!  DECOUPLING USING MATRIX ELEMENTS
!
!-----------------------------------------------------------------------
      use floatPrecision
      use numerical_constants
      use mathlib_bouncer
      use crcoall
      use parpro
      use mod_common
      use mod_commons
      use mod_commont
      implicit none
      integer i,ierr,j,no
      real(kind=fPrec) aa,bb,d1,dpp,dsm,qw,qwc,sen,sn,ss
      dimension aa(6,6),bb(6),dsm(6),sn(6),sen(6),ss(6)
      dimension qwc(3),qw(2),d1(6)
      save
!-----------------------------------------------------------------------
      do i=1,6
        bb(i)=zero
        dsm(i)=zero
        sn(i)=zero
        sen(i)=zero
        ss(i)=zero
        d1(i)=zero

        do j=1,6
          aa(i,j)=zero
        end do
      end do

      do 20 i=1,3
        qwc(i)=zero
   20 continue
      dpp=zero
      write(lout,10000)
      call betalf(dpp,qw)
      call phasad(dpp,qwc)
      sen(1)=ta(3,1)
      sen(2)=ta(3,2)
      sen(3)=ta(4,1)
      sen(4)=ta(4,2)
      if(iskew.eq.1) then
        sen(5)=qwc(1)
        sen(6)=qwc(2)
      endif
      do 30 i=1,6
        if(iskew.eq.2.and.i.gt.4) goto 30
        if(i.le.4) then
          sn(i)=ed(nskew(i))
          dsm(i)=dsm0
          bb(i)=sen(i)
        else
          if (abs(el(nskew(i))).le.pieni) then
            sn(i)=ed(nskew(i))
          else
            sn(i)=ek(nskew(i))
          endif
          dsm(i)=dkq
          bb(i)=sen(i)-qwsk(i-4)
        endif
        ss(i)=sen(i)
   30 continue
      do 100 no=1,itcro
        do 40 i=1,6
          if(iskew.eq.2.and.i.gt.4) goto 40
          if(i.le.4) then
            ed(nskew(i))=ed(nskew(i))+dsm(i)
          else
            if (abs(el(nskew(i))).le.pieni) then
              ed(nskew(i))=ed(nskew(i))+dsm(i)
            else
              ek(nskew(i))=ek(nskew(i))+dsm(i)
            endif
          endif
          if(kp(nskew(i)).eq.5) call combel(nskew(i))
          call betalf(dpp,qw)
          call phasad(dpp,qwc)
          aa(i,1)=(ta(3,1)-ss(1))/dsm(i)
          aa(i,2)=(ta(3,2)-ss(2))/dsm(i)
          aa(i,3)=(ta(4,1)-ss(3))/dsm(i)
          aa(i,4)=(ta(4,2)-ss(4))/dsm(i)
          if(iskew.eq.1) then
            aa(i,5)=(qwc(1)-ss(5))/dsm(i)
            aa(i,6)=(qwc(2)-ss(6))/dsm(i)
          endif
          if(i.le.4) then
            ed(nskew(i))=ed(nskew(i))-dsm(i)
          else
            if (abs(el(nskew(i))).le.pieni) then
              ed(nskew(i))=ed(nskew(i))-dsm(i)
            else
              ek(nskew(i))=ek(nskew(i))-dsm(i)
            endif
          endif
          if(kp(nskew(i)).eq.5) call combel(nskew(i))
   40   continue
        if(iskew.eq.1) then
          call loesd(aa,bb,6,6,ierr)
        else if(iskew.eq.2) then
          call loesd(aa,bb,4,4,ierr)
        endif
        if(ierr == 1) then
          write(lout,"(a)") "DECOUP> ERROR Problems during matrix-inversion."
          call prror(-1)
        end if
        do 50 i=1,6
          if(iskew.eq.2.and.i.gt.4) goto 50
          if(i.le.4) then
            ed(nskew(i))=ed(nskew(i))-bb(i)
          else
            if (abs(el(nskew(i))).le.pieni) then
              ed(nskew(i))=ed(nskew(i))-bb(i)
            else
              ek(nskew(i))=ek(nskew(i))-bb(i)
            endif
          endif
          if(kp(nskew(i)).eq.5) call combel(nskew(i))
   50   continue
        call betalf(dpp,qw)
        call phasad(dpp,qwc)
        ss(1)=ta(3,1)
        ss(2)=ta(3,2)
        ss(3)=ta(4,1)
        ss(4)=ta(4,2)
        if(iskew.eq.1) then
          ss(5)=qwc(1)
          ss(6)=qwc(2)
        endif
        write(lout,10010)
        write(lout,10020) no,sen(1),ss(1),sen(2),ss(2),sen(3),ss(3), sen&
     &(4),ss(4)
        write(lout,10030) bez(nskew(1)),sn(1),ed(nskew(1)),             &
     &bez(nskew(2)),sn                                                  &
     &(2),ed(nskew(2)),bez(nskew(3)),sn(3),ed(nskew(3)), bez            &
     &(nskew(4)),sn(4),ed(nskew(4))
        if(iskew.eq.1) then
          write(lout,10010)
          write(lout,10040) qwsk(1),qwc(1),qwsk(2),qwc(2)
          if (abs(el(nskew(5))).le.pieni) then
            write(lout,10060) sn(5),ed(nskew(5)),nskew(5),sn(6),ed      &
     &(nskew(6)), nskew(6)
          else
            write(lout,10060) sn(5),ek(nskew(5)),nskew(5),sn(6),ek      &
     &(nskew(6)), nskew(6)
          endif
        else if(iskew.eq.2) then
          write(lout,10010)
          write(lout,10050) qwc(1),qwc(2)
        endif
        do 60 i=1,6
          if(iskew.eq.2.and.i.gt.4) goto 60
          if(i.le.4) then
            d1(i)=abs(ss(i))
          else
            d1(i)=abs(ss(i)-qwsk(i-4))
          endif
   60   continue
        do 70 i=1,6
          if(iskew.eq.2.and.i.gt.4) goto 70
          if(d1(i).gt.dsi) goto 80
   70   continue
        return
   80   do 90 i=1,6
          if(iskew.eq.2.and.i.gt.4) goto 90
          if(i.le.4) then
            bb(i)=ss(i)
          else
            bb(i)=ss(i)-qwsk(i-4)
          endif
   90   continue
  100 continue
!-----------------------------------------------------------------------
      return
10000 format(t5,'---- ENTRY DECOUP ----')
10010 format(131('-'))
10020 format(/10x,'DECOUPLING ROUTINE  ITERATION #',i2// 30x,           &
     &'BEFORE         AFTER     DECOUPLING'// 17x,'   M(3,1)      ',2g15&
     &.5/ 17x,'   M(3,2)      ',2g15.5/ 17x,'   M(4,1)      ',2g15.5/ 17&
     &x,'   M(4,2)      ',2g15.5// 5x,'SKEW QUDRUPOLE STRENGTHS')
10040 format(10x,'Q-VARIATION' / 10x,                                   &
     &'Q-VALUE            THEORET.        AFTER     COMPENSATION'/ 10x, &
     &'HORIZONTAL     ',2g15.7/ 10x,'VERTICAL       ',2g15.7/)
10050 format(10x,'CURRENT TUNE' / 10x,'Q-VALUE'/ 10x,'HORIZONTAL     ', &
     &g15.7/ 10x,'VERTICAL       ',g15.7/)
10060 format(10x,'QUADRU.STRENGTH',2g15.8,'   INDEX ',i3/ 10x,          &
     &'               ',2g15.8,'         ',i3)
10030 format(14x,a16,2x,g17.10,1x,g17.10/14x,a16,2x,g17.10,1x,          &
     &g17.10/14x,a16,2x,g17.10,1x,g17.10/14x,a16,2x,g17.10,1x,g17.10)
end subroutine decoup

subroutine datime(nd,nt)
      implicit none
! Fill common slate for usage by hmachi call as per z007 writeup.        !hr08
      common /slate/ isl(40)                                             !hr08

      integer isl                                                        !hr08
!
!-    call datime (nd,nt)   returns integer date   nd = yymmdd
!-                                  integer time   nt =   hhmm
!     integer nd,nt,mm(3),nn(3)
!     call idate (mm(1),mm(2),mm(3))
!     call itime (nn)
      character(len=8) date
      character(len=10) time
      character(len=5) zone
      integer values(8),mm(3),nd,nt
      save
      call date_and_time(date,time,zone,values)
      mm(3)=mod(values(1),100)
!     mm(3) = mod (mm(3),100)
      mm(2)=values(3)
      mm(1)=values(2)
      isl(1)= mm(3)                                                      !hr08
      isl(2)= mm(2)                                                      !hr08
      isl(3)= mm(1)                                                      !hr08
      isl(4)= values(5)                                                  !hr08
      isl(5)= values(6)                                                  !hr08
      isl(6)= 0                                                          !hr08
      nd = (mm(3)*100+mm(1))*100 + mm(2)
!     nt =            nn(1) *100 + nn(2)
      nt=values(5)*100+values(6)
      return
end subroutine datime
