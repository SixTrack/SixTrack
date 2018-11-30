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

    subroutine hash_md5FinalPtr(ctxID, md5Vals, md5Size) bind(C, name="md5wrapper_md5Final")
      use, intrinsic :: iso_c_binding
      integer(kind=C_INT), value, intent(in)    :: ctxID
      integer(kind=C_INT),        intent(inout) :: md5Vals(*)
      integer(kind=C_INT), value, intent(in)    :: md5Size
    end subroutine hash_md5FinalPtr

  end interface

contains

subroutine hash_md5Final(ctxID, md5Digest)

  use, intrinsic :: iso_c_binding
  use string_tools

  integer,           intent(in)  :: ctxID
  character(len=32), intent(out) :: md5Digest

  integer(kind=C_INT) tmpVals(16)
  integer     i


  tmpVals(:) = 0
  call hash_md5FinalPtr(ctxID, tmpVals, 16)
  write(md5Digest,"(16(z2.2))") tmpVals
  md5Digest = chr_toLower(md5Digest)
  ! do i=1,16
  !   write(6,"(a,i2,a,i3,a,z2.2)") "F> ",i," = ",tmpVals(i)," = 0x",tmpVals(i)
  ! end do

end subroutine hash_md5Final

end module mod_hash
