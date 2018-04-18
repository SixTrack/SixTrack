!     This file contains variables that are shared between all
!     parts of the code, such as lout and crlibm functions
      
+cd crlibco
      double precision sin_rn,cos_rn,tan_rn,sinh_rn,cosh_rn,asin_rn,    &
     &acos_rn,atan_rn,atan2_rn,exp_rn,log_rn,log10_rn
+cd crcoall
!     Standard output unit
!     For CR version, this is the "buffer file" fort.92;
!     Otherwise write directly to "*" aka iso_fortran_env::output_unit (usually unit 6)
      integer lout
      common /crflags/lout
+cd ffieldcommon
!-----------------------------------------------------------------------
!    Fringe Field interface parametter with SixTrack (ASIMONA, BDALENA, TPUGNAT)
!-----------------------------------------------------------------------
      integer :: ffNLn                                           ! Number of line in the file
      integer :: FFindex
      character(len=100), pointer, dimension(:) :: ffQNames      ! Names Quad
      common/FringeField/ FFindex(nele),ffQNames, ffNLn
+cd ffieldcommon1
      logical ffReady                                            ! Check data read with success
      integer :: ffNLFile                                        ! Number of line in the file
      double precision, pointer, dimension(:,:) :: ffParam       ! Kin, Lin, Corin, Kex, Lex, Corex
      character(len=300), pointer, dimension(:) :: ffFNames
      integer, pointer, dimension(:,:) :: ffQ2File
      common/FringeFieldTab/ ffNLFile, ffReady, ffParam, ffQ2File, ffFNames
