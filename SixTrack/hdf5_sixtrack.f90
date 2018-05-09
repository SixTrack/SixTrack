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
  
contains

subroutine h5lin_init

  h5lin_nElem(1) = 0
  h5lin_nElem(2) = 1000

  call alloc(h5lin_iData,               h5lin_nElem(2), 0,                            "h5lin_iData")
  call alloc(h5lin_cData, max_name_len, h5lin_nElem(2), repeat(char(0),max_name_len), "h5lin_cData")
  call alloc(h5lin_rData, 17,           h5lin_nElem(2), 0.0_fPrec,                    "h5lin_rData")

end subroutine h5lin_init

subroutine h5lin_writeLine(lineNo, elemType, arrData)
  
  integer,          intent(in) :: lineNo
  character(len=*), intent(in) :: elemType
  real(kind=fPrec), intent(in) :: arrData(17)
  
  h5lin_nElem(1) = h5lin_nElem(1) + 1
  
  if(h5lin_nElem(1) > h5lin_nElem(2)) then
    h5lin_nElem(2) = h5lin_nElem(2) + 1000
    call resize(h5lin_iData,               h5lin_nElem(2), 0,                            "h5lin_iData")
    call resize(h5lin_cData, max_name_len, h5lin_nElem(2), repeat(char(0),max_name_len), "h5lin_cData")
    call resize(h5lin_rData, 17,           h5lin_nElem(2), 0.0_fPrec,                    "h5lin_rData")
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
  setFields(2)  = h5_dataField(name="TYP",    type=h5_typeChar, size=max_name_len)
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
  call h5_createDataSet("linopt", h5_rootID, optFmt, optSet, 2000)
  
  call h5_prepareWrite(optSet, h5lin_nElem(1))
  call h5_writeData(optSet, 1, h5lin_nElem(1), h5lin_iData)
  call h5_writeData(optSet, 2, h5lin_nElem(1), h5lin_cData)
  do rIdx=1,17
    call h5_writeData(optSet, rIdx+2, h5lin_nElem(1), h5lin_rData(rIdx,:))
  end do
  call h5_finaliseWrite(optSet)
  
  deallocate(setFields)
  
  call dealloc(h5lin_iData,               "h5lin_iData")
  call dealloc(h5lin_cData, max_name_len, "h5lin_cData")
  call dealloc(h5lin_rData,               "h5lin_rData")
  
end subroutine h5lin_saveData

end module hdf5_linopt

#ifdef COLLIMAT
! ================================================================================================ !
!  Special HDF5 Module for Writing Tracks2 to HDF5
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-09
! ================================================================================================ !
module hdf5_tracks2
  
  use hdf5_output
  use mod_alloc
  use mod_commons
  
  implicit none
  
  integer,                       private, save :: h5tr2_nElem(2)
  integer,          allocatable, private, save :: h5tr2_iData(:,:)
  real(kind=fPrec), allocatable, private, save :: h5tr2_rData(:,:)
  
contains

subroutine h5tr2_init
  
  h5tr2_nElem(1) = 0
  h5tr2_nElem(2) = 500
  
  call alloc(h5tr2_iData, 2, h5tr2_nElem(2), 0,         "h5tr2_iData")
  call alloc(h5tr2_rData, 6, h5tr2_nElem(2), 0.0_fPrec, "h5tr2_rData")
  
end subroutine h5tr2_init

subroutine h5tr2_writeLine(pid, turn, s, x, xp, y, yp, dee, type)
  
  integer,          intent(in) :: pid, turn, type
  real(kind=fPrec), intent(in) :: s, x, xp, y, yp, dee
  
  h5tr2_nElem(1) = h5tr2_nElem(1) + 1
  
  if(h5tr2_nElem(1) > h5tr2_nElem(2)) then
    h5tr2_nElem(2) = h5tr2_nElem(2) + 500
    call resize(h5tr2_iData, 2, h5tr2_nElem(2), 0,         "h5tr2_iData")
    call resize(h5tr2_rData, 6, h5tr2_nElem(2), 0.0_fPrec, "h5tr2_rData")
  end if
  
  h5tr2_iData(:,h5tr2_nElem(1)) = (/pid,turn,type/)
  h5tr2_rData(:,h5tr2_nElem(1)) = (/s,x,xp,y,yp,dee/)
  
end subroutine h5tr2_writeLine

subroutine h5tr2_saveData
  
  type(h5_dataField), allocatable :: setFields(:)
  integer tr2Fmt, tr2Set, rIdx
  
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
  
  call h5_createFormat("tracksTwo", setFields, tr2Fmt)
  call h5_createDataSet("tracks2", h5_rootID, tr2Fmt, tr2Set, 1000)
  
  call h5_prepareWrite(tr2Set, h5tr2_nElem(1))
  call h5_writeData(tr2Set, 1, h5tr2_nElem(1), h5tr2_iData(1,:))
  call h5_writeData(tr2Set, 2, h5tr2_nElem(1), h5tr2_iData(2,:))
  call h5_writeData(tr2Set, 3, h5tr2_nElem(1), h5tr2_rData(1,:))
  call h5_writeData(tr2Set, 4, h5tr2_nElem(1), h5tr2_rData(2,:))
  call h5_writeData(tr2Set, 5, h5tr2_nElem(1), h5tr2_rData(3,:))
  call h5_writeData(tr2Set, 6, h5tr2_nElem(1), h5tr2_rData(4,:))
  call h5_writeData(tr2Set, 7, h5tr2_nElem(1), h5tr2_rData(5,:))
  call h5_writeData(tr2Set, 8, h5tr2_nElem(1), h5tr2_rData(6,:))
  call h5_writeData(tr2Set, 9, h5tr2_nElem(1), h5tr2_iData(3,:))
  call h5_finaliseWrite(tr2Set)
  
  deallocate(setFields)
  
  call dealloc(h5tr2_iData, "h5tr2_iData")
  call dealloc(h5tr2_rData, "h5tr2_rData")
  
end subroutine h5tr2_saveData

end module hdf5_tracks2

! END IFDEF COLLIMAT
#endif

! END IFDEF HDF5
#endif
