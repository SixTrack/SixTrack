#ifdef HDF5

! ================================================================================================ !
!  Special HDF5 Module for Writing Linear Optics Parameters
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-09
! ================================================================================================ !
module hdf5_linopt

  use hdf5_output
  use mod_alloc
  use mod_commons

  implicit none

  integer,                       private, save :: h5lin_nElem(2)
  integer,          allocatable, private, save :: h5lin_iData(:)
  character(len=:), allocatable, private, save :: h5lin_cData(:)
  real(kind=fPrec), allocatable, private, save :: h5lin_rData(:,:)

  ! The HDF5 chunk size for writing linear optics
  ! This should be less than 1Mb/188b = 5577
  integer, parameter :: h5lin_chunkSize = 2000

contains

subroutine h5lin_init

  h5lin_nElem(1) = 0
  h5lin_nElem(2) = h5lin_chunkSize

  call alloc(h5lin_iData,               h5lin_nElem(2), 0,                            "h5lin_iData")
  call alloc(h5lin_cData, mNameLen, h5lin_nElem(2), repeat(char(0),mNameLen), "h5lin_cData")
  call alloc(h5lin_rData, 17,           h5lin_nElem(2), 0.0_fPrec,                    "h5lin_rData")

end subroutine h5lin_init

subroutine h5lin_writeLine(lineNo, elemType, arrData)

  integer,          intent(in) :: lineNo
  character(len=*), intent(in) :: elemType
  real(kind=fPrec), intent(in) :: arrData(17)

  h5lin_nElem(1) = h5lin_nElem(1) + 1

  if(h5lin_nElem(1) > h5lin_nElem(2)) then
    h5lin_nElem(2) = h5lin_nElem(2) + h5lin_chunkSize
    call alloc(h5lin_iData,               h5lin_nElem(2), 0,                            "h5lin_iData")
    call alloc(h5lin_cData, mNameLen, h5lin_nElem(2), repeat(char(0),mNameLen), "h5lin_cData")
    call alloc(h5lin_rData, 17,           h5lin_nElem(2), 0.0_fPrec,                    "h5lin_rData")
  end if

  h5lin_iData(h5lin_nElem(1))   = lineNo
  h5lin_cData(h5lin_nElem(1))   = elemType
  h5lin_rData(:,h5lin_nElem(1)) = arrData

end subroutine h5lin_writeLine

subroutine h5lin_saveData

  type(h5_dataField), allocatable :: setFields(:)
  integer optFmt, optSet, rIdx

  allocate(setFields(19))

  setFields(1)  = h5_dataField(name="NR",     type=h5_typeInt)
  setFields(2)  = h5_dataField(name="TYP",    type=h5_typeChar, size=mNameLen)
  setFields(3)  = h5_dataField(name="LTOT",   type=h5_typeReal)
  setFields(4)  = h5_dataField(name="PHIX",   type=h5_typeReal)
  setFields(5)  = h5_dataField(name="BETAX",  type=h5_typeReal)
  setFields(6)  = h5_dataField(name="ALPHAX", type=h5_typeReal)
  setFields(7)  = h5_dataField(name="GAMMAX", type=h5_typeReal)
  setFields(8)  = h5_dataField(name="DISX",   type=h5_typeReal)
  setFields(9)  = h5_dataField(name="DISXP",  type=h5_typeReal)
  setFields(10) = h5_dataField(name="CLOX",   type=h5_typeReal)
  setFields(11) = h5_dataField(name="CLOXP",  type=h5_typeReal)
  setFields(12) = h5_dataField(name="PHIY",   type=h5_typeReal)
  setFields(13) = h5_dataField(name="BETAY",  type=h5_typeReal)
  setFields(14) = h5_dataField(name="ALPHAY", type=h5_typeReal)
  setFields(15) = h5_dataField(name="GAMMAY", type=h5_typeReal)
  setFields(16) = h5_dataField(name="DISY",   type=h5_typeReal)
  setFields(17) = h5_dataField(name="DISYP",  type=h5_typeReal)
  setFields(18) = h5_dataField(name="CLOY",   type=h5_typeReal)
  setFields(19) = h5_dataField(name="CLOYP",  type=h5_typeReal)

  call h5_createFormat("linearOptics", setFields, optFmt)
  call h5_createDataSet("linopt", h5_rootID, optFmt, optSet, h5lin_chunkSize)

  call h5_prepareWrite(optSet, h5lin_nElem(1))
  call h5_writeData(optSet, 1, h5lin_nElem(1), h5lin_iData)
  call h5_writeData(optSet, 2, h5lin_nElem(1), h5lin_cData)
  do rIdx=1,17
    call h5_writeData(optSet, rIdx+2, h5lin_nElem(1), h5lin_rData(rIdx,:))
  end do
  call h5_finaliseWrite(optSet)

  deallocate(setFields)

  call dealloc(h5lin_iData,               "h5lin_iData")
  call dealloc(h5lin_cData, mNameLen, "h5lin_cData")
  call dealloc(h5lin_rData,               "h5lin_rData")

end subroutine h5lin_saveData

end module hdf5_linopt

! ================================================================================================ !
!  Special HDF5 Module for Writing Tracks2 to HDF5
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-09
! ================================================================================================ !
module hdf5_tracks2

  use hdf5_output
  use mod_alloc
  use mod_commons
  use crcoall

  implicit none

  integer,          allocatable, private, save :: h5tr2_iData(:,:)
  real(kind=fPrec), allocatable, private, save :: h5tr2_rData(:,:)

  integer,                       private, save :: h5tr2_nElem
  integer,                       private, save :: h5tr2_fmtID
  integer,                       private, save :: h5tr2_setID

  ! The HDF5 chunk size for writing tracks2
  ! This should be less than 1Mb/60b = 17476
  integer, parameter :: h5tr2_chunkSize = 10000

contains

subroutine h5tr2_init

  type(h5_dataField), allocatable :: setFields(:)

  h5tr2_nElem = 0

  allocate(setFields(9))

  setFields(1) = h5_dataField(name="PID",  type=h5_typeInt)
  setFields(2) = h5_dataField(name="TURN", type=h5_typeInt)
  setFields(3) = h5_dataField(name="S",    type=h5_typeReal)
  setFields(4) = h5_dataField(name="X",    type=h5_typeReal)
  setFields(5) = h5_dataField(name="XP",   type=h5_typeReal)
  setFields(6) = h5_dataField(name="Y",    type=h5_typeReal)
  setFields(7) = h5_dataField(name="YP",   type=h5_typeReal)
  setFields(8) = h5_dataField(name="DE/E", type=h5_typeReal)
  setFields(9) = h5_dataField(name="TYPE", type=h5_typeInt)

  call alloc(h5tr2_iData, 2, h5tr2_chunkSize, 0,         "h5tr2_iData")
  call alloc(h5tr2_rData, 6, h5tr2_chunkSize, 0.0_fPrec, "h5tr2_rData")

  call h5_createFormat("tracksTwo", setFields, h5tr2_fmtID)
  call h5_createDataSet("tracks2", h5_rootID, h5tr2_fmtID, h5tr2_setID, h5tr2_chunkSize)

  deallocate(setFields)

end subroutine h5tr2_init

subroutine h5tr2_finalise

  call h5tr2_saveData

  call dealloc(h5tr2_iData, "h5tr2_iData")
  call dealloc(h5tr2_rData, "h5tr2_rData")

end subroutine h5tr2_finalise

subroutine h5tr2_writeLine(pid, turn, s, x, xp, y, yp, dee, type)

  integer,          intent(in) :: pid, turn, type
  real(kind=fPrec), intent(in) :: s, x, xp, y, yp, dee

  h5tr2_nElem = h5tr2_nElem + 1

  if(h5tr2_nElem >= h5tr2_chunkSize) then
    call h5tr2_saveData
  end if

  h5tr2_iData(:,h5tr2_nElem) = (/pid,turn,type/)
  h5tr2_rData(:,h5tr2_nElem) = (/s,x,xp,y,yp,dee/)

end subroutine h5tr2_writeLine

subroutine h5tr2_saveData

  call h5_prepareWrite(h5tr2_setID, h5tr2_nElem)
  call h5_writeData(h5tr2_setID, 1, h5tr2_nElem, h5tr2_iData(1,1:h5tr2_nElem))
  call h5_writeData(h5tr2_setID, 2, h5tr2_nElem, h5tr2_iData(2,1:h5tr2_nElem))
  call h5_writeData(h5tr2_setID, 3, h5tr2_nElem, h5tr2_rData(1,1:h5tr2_nElem))
  call h5_writeData(h5tr2_setID, 4, h5tr2_nElem, h5tr2_rData(2,1:h5tr2_nElem))
  call h5_writeData(h5tr2_setID, 5, h5tr2_nElem, h5tr2_rData(3,1:h5tr2_nElem))
  call h5_writeData(h5tr2_setID, 6, h5tr2_nElem, h5tr2_rData(4,1:h5tr2_nElem))
  call h5_writeData(h5tr2_setID, 7, h5tr2_nElem, h5tr2_rData(5,1:h5tr2_nElem))
  call h5_writeData(h5tr2_setID, 8, h5tr2_nElem, h5tr2_rData(6,1:h5tr2_nElem))
  call h5_writeData(h5tr2_setID, 9, h5tr2_nElem, h5tr2_iData(3,1:h5tr2_nElem))
  call h5_finaliseWrite(h5tr2_setID)

  if(h5_debugOn) then
    write(lout,"(a,i0,a)") "HDF5> DEBUG Wrote ",h5tr2_nElem," lines to tracks2."
  end if

  h5tr2_nElem      = 1
  h5tr2_iData(:,:) = 0
  h5tr2_rData(:,:) = 0.0_fPrec

end subroutine h5tr2_saveData

end module hdf5_tracks2

! END IFDEF HDF5
#endif
