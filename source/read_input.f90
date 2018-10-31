! ================================================================================================ !
!  Read Input
! ~~~~~~~~~~~~
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-18
!  Module for parsing special input files like fort.33, fort.13, fort.10, etc.
! ================================================================================================ !

module read_input

  use crcoall

  implicit none

contains

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-10-31
!  Rewritten parsing of initial coordinates from fort.13. Moved from maincr.
! ================================================================================================ !
subroutine readFort13

  use parpro
  use mod_hions
  use mod_common
  use mod_commonmn

  implicit none

  integer ia

#ifdef CRLIBM
  integer nchars
  parameter (nchars=160)
  character(len=nchars) ch
  character(len=nchars+nchars) ch1
  real(kind=fPrec) round_near
#endif

        do ia=1,napx,2
#ifndef CRLIBM
          read(13,*,iostat=ierro) xv(1,ia),yv(1,ia),xv(2,ia),yv(2,ia),  &
     &sigmv(ia),dpsv(ia),xv(1,ia+1),yv(1,ia+1),xv(2,ia+1),yv            &
     &(2,ia+1), sigmv(ia+1),dpsv(ia+1),e0,ejv(ia),ejv(ia+1)
#endif
#ifdef CRLIBM
          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ xv(1,ia)]"
             call prror(-1)
          endif
          xv(1,ia) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV xv(1,ia)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ yv(1,ia)]"
             call prror(-1)
          endif
          yv(1,ia) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV yv(1,ia)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ xv(2,ia)]"
             call prror(-1)
          endif
          xv(2,ia) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV xv(2,ia)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ yv(2,ia)]"
             call prror(-1)
          endif
          yv(2,ia) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV yv(2,ia)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ sigmv(ia)]"
             call prror(-1)
          endif
          sigmv(ia) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV sigmv(ia)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ dpsv(ia)]"
             call prror(-1)
          endif
          dpsv(ia) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV dpsv(ia)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ xv(1,ia+1)]"
             call prror(-1)
          endif
          xv(1,ia+1) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV xv(1,ia+1)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ yv(1,ia+1)]"
             call prror(-1)
          endif
          yv(1,ia+1) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV yv(1,ia+1)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ xv(2,ia+1)]"
             call prror(-1)
          endif
          xv(2,ia+1) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV xv(2,ia+1)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ yv(2,ia+1)]"
             call prror(-1)
          endif
          yv(2,ia+1) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV yv(2,ia+1)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)                                              &
     &            "Error when reading fort.13 [READ sigmv(ia+1)]"
             call prror(-1)
          endif
          sigmv(ia+1) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)                                              &
     &            "Error when reading fort.13 [CONV sigmv(ia+1)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)                                              &
     &            "Error when reading fort.13 [READ dpsv(ia+1)]"
             call prror(-1)
          endif
          dpsv(ia+1) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)                                              &
     &            "Error when reading fort.13 [CONV dpsv(ia+1)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ e0]"
             call prror(-1)
          endif
          e0 = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV e0]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ ejv(ia)]"
             call prror(-1)
          endif
          ejv(ia) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV ejv(ia)]"
             call prror(-1)
          endif

          read(13,'(a)', iostat=ierro) ch
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [READ ejv(ia+1)]"
             call prror(-1)
          endif
          ejv(ia+1) = round_near(ierro,nchars,ch)
          if(ierro.gt.0) then
             write(lout,*)"Error when reading fort.13 [CONV ejv(ia+1)]"
             call prror(-1)
          endif

#endif
          if(ierro.ne.0) call prror(56)
          mtc(ia)=one
          mtc(ia+1)=one
          nucm(ia)=nucm0
          nucm(ia+1)=nucm0
          e0f=sqrt(e0**2-nucm0**2)                                         !hr05
          ejfv(ia)=sqrt(ejv(ia)**2-nucm(ia)**2)                               !hr05
          ejfv(ia+1)=sqrt(ejv(ia+1)**2-nucm(ia+1)**2)                           !hr05
          oidpsv(ia)=one/(one+dpsv(ia))
          oidpsv(ia+1)=one/(one+dpsv(ia+1))
          moidpsv(ia)=mtc(ia)/(one+dpsv(ia))
          moidpsv(ia+1)=mtc(ia+1)/(one+dpsv(ia+1))
          omoidpsv(ia)=c1e3*((one-mtc(ia))*oidpsv(ia))
          omoidpsv(ia+1)=c1e3*((one-mtc(ia+1))*oidpsv(ia+1))
        end do
end subroutine readFort13

! ================================================================================================ !
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Last Modified: 2018-05-18
!  Reading fort.33. Moved from maincr/mainda.
! ================================================================================================ !
subroutine readFort33

  use string_tools
  use mod_common,     only : clo6, clop6
  use sixtrack_input, only : sixin_echoVal
  use mod_settings
  use mod_units
  use parpro

  implicit none

  character(len=:), allocatable :: lnSplit(:)
  character(len=mInputLn)       :: inLine
  integer nSplit, ioStat
  logical spErr, fErr

  call units_openUnit(33,"fort.33",.true.,"r",fErr)

  read(33,"(a)",iostat=ioStat) inLine
  if(ioStat /= 0) then
    write(lout,"(a,i0)") "READ33> ERROR Failed to read line from 'fort.33'. iostat = ",ioStat
    call prror(-1)
  end if

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "READ33> ERROR Failed to parse line from 'fort.33'"
    call prror(-1)
  end if

  if(nSplit > 0) call chr_cast(lnSplit(1),clo6(1), spErr)
  if(nSplit > 1) call chr_cast(lnSplit(2),clop6(1),spErr)
  if(nSplit > 2) call chr_cast(lnSplit(3),clo6(2), spErr)
  if(nSplit > 3) call chr_cast(lnSplit(4),clop6(2),spErr)
  if(nSplit > 4) call chr_cast(lnSplit(5),clo6(3), spErr)
  if(nSplit > 5) call chr_cast(lnSplit(6),clop6(3),spErr)

  if(st_debug) then
    call sixin_echoVal("clo6(1)", clo6(1), "READ33",1)
    call sixin_echoVal("clop6(1)",clop6(1),"READ33",1)
    call sixin_echoVal("clo6(2)", clo6(2), "READ33",1)
    call sixin_echoVal("clop6(2)",clop6(2),"READ33",1)
    call sixin_echoVal("clo6(3)", clo6(3), "READ33",1)
    call sixin_echoVal("clop6(3)",clop6(3),"READ33",1)
  end if

end subroutine readFort33

end module read_input
