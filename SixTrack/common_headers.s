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
