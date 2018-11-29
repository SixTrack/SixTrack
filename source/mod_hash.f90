! ================================================================================================ !
!  HASH MODULE
!  V.k. Berglyd Olsen, BE-ABP-HSS
!  Last modified: 2018-11-29
!
!  This module provides an interface to the MD5 implementation written by Ronald L. Rivest (MIT).
!  The source code is available under source/md5. Note that this implementation is not compliant
!  with RFC1321, as it produces different hash values than the standard.
! ================================================================================================ !

module mod_hash

  use floatPrecision

  implicit none

  ! C Interface
  interface

    subroutine hash_digestFloatArray(inArr, arrLen) bind(C, name="md5wrapper_digestFloatArray")
      use, intrinsic :: iso_c_binding
      real(kind=C_DOUBLE),        intent(in) :: inArr
      integer(kind=C_INT), value, intent(in) :: arrLen
    end subroutine hash_digestFloatArray

    subroutine hash_md5Init(nInst) bind(C, name="md5wrapper_md5Init")
      use, intrinsic :: iso_c_binding
      integer(kind=C_INT), value, intent(in) :: nInst
    end subroutine hash_md5Init

    subroutine hash_md5Update(ctxID, inStr, strLen) bind(C, name="md5wrapper_md5Update")
      use, intrinsic :: iso_c_binding
      integer(kind=C_INT), value,   intent(in) :: ctxID
      character(kind=C_CHAR,len=1), intent(in) :: inStr
      integer(kind=C_INT), value,   intent(in) :: strLen
    end subroutine hash_md5Update

    subroutine hash_md5Final(ctxID) bind(C, name="md5wrapper_md5Final")
      use, intrinsic :: iso_c_binding
      integer(kind=C_INT), value, intent(in) :: ctxID
    end subroutine hash_md5Final

  end interface

contains



end module mod_hash
