module mod_splitSingletrack
  use iso_fortran_env, only : output_unit
  implicit none

  ! File units used for I/O;
  ! can be overridden USEing if clashing
  ! with units used in rest of the code.
  integer :: inUnit = 90
  integer :: outUnit = 91

  integer :: textOutUnit = output_unit

  ! Temp arrays for I/O
  character(len=80), private :: sixtit,commen
  character(len=8),  private :: cdate,ctime,progrm
  integer, private           :: i,ia,ia0,ifipa,ilapa,itopa,j,ntwin,numl,n
  double precision, private  :: dpscor,sigcor
  integer, private           ::  icode
  double precision, private  :: b,c,c1,d,d1,dizu0,dmmac,dnms,dnumlr,dummy,e,e1,f,f1,g,g1,h,h1,p,p1
  double precision, private  :: qwc(3)
  double precision, private  :: clo(3)
  double precision, private  :: clop(3)
  double precision, private  :: di0(3)
  double precision, private  :: dip0(3)
  double precision, private  :: ta(6,6)

contains
  function numSTFpairs(ifname)
    character(len=*), intent(in) :: ifname
    integer numSTFpairs
    integer iostat

    open(unit=inUnit, file=ifname, form='UNFORMATTED', status='OLD', action='READ', err=510, iostat=iostat)

    read(inUnit,end=511,err=520,iostat=iostat) &
         sixtit,commen,cdate,ctime,       &
         progrm,ifipa,ilapa,itopa,icode,numl,qwc(1),qwc(2),qwc(3), clo(1), &
         clop(1),clo(2),clop(2),clo(3),clop(3), di0(1),dip0(1),di0(2),dip0(2), &
         dummy,dummy, ta(1,1),ta(1,2),ta(1,3),ta(1,4),ta(1,5),ta(1,6), &
         ta(2,1),ta(2,2),ta(2,3),ta(2,4),ta(2,5),ta(2,6), ta(3,1),ta(3,2), &
         ta(3,3),ta(3,4),ta(3,5),ta(3,6), ta(4,1),ta(4,2),ta(4,3),ta(4,4), &
         ta(4,5),ta(4,6), ta(5,1),ta(5,2),ta(5,3),ta(5,4),ta(5,5),ta(5,6), &
         ta(6,1),ta(6,2),ta(6,3),ta(6,4),ta(6,5),ta(6,6), dmmac,dnms,dizu0,&
         dnumlr,sigcor,dpscor
    close(inUnit)

    if(itopa <= 0) then
      write(textOutUnit,'(a)')    "ERR in numSTFpairs: Corrupt itopa in first header of file '" // ifname // "'"
      write(textOutUnit,'(a,i0,a)') "itopa = ", itopa, " <= 0"
      stop 1
    end if

    if (modulo(itopa,2) /= 0) then
      write(textOutUnit,'(a)')    "ERR in numSTFpairs: Number of particles not divisible by 2."
      write(textOutUnit,'(a,i6)') "ERR in numSTFpairs: itopa = ", itopa
      stop 1
    end if

    numSTFpairs = itopa/2
    return

    !Error handling
510 continue
    write(textOutUnit,'(a,i6)') "ERR while opening file; iostat=",iostat
    stop 1
511 continue
    write(textOutUnit,'(a,i6)') "END while reading header; iostat=",iostat
    stop 1
520 continue
    write(textOutUnit,'(a,i6)') "ERR while reading header; iostat=",iostat
    stop 1
  end function numSTFpairs

  subroutine convertSTFpair(pairIdx, ifname, ofname)
    integer, intent(in)          :: pairIdx !Which pair to extract
    character(len=*), intent(in) :: ifname  !Name of the singletrackfile to split
    character(len=*), intent(in) :: ofname  !Name of the output file

    integer iostat

    integer iPair

    open(unit=inUnit,  file=ifname, form='UNFORMATTED', status='OLD',     action='READ',  err=410, iostat=iostat)
    open(unit=outUnit, file=ofname, form='UNFORMATTED', status='REPLACE', action='WRITE', err=411, iostat=iostat)

    headers: do
      read(inUnit,end=511,err=520,iostat=iostat) &
           sixtit,commen,cdate,ctime,       &
           progrm,ifipa,ilapa,itopa,icode,numl,qwc(1),qwc(2),qwc(3), clo(1), &
           clop(1),clo(2),clop(2),clo(3),clop(3), di0(1),dip0(1),di0(2),dip0(2), &
           dummy,dummy, ta(1,1),ta(1,2),ta(1,3),ta(1,4),ta(1,5),ta(1,6), &
           ta(2,1),ta(2,2),ta(2,3),ta(2,4),ta(2,5),ta(2,6), ta(3,1),ta(3,2), &
           ta(3,3),ta(3,4),ta(3,5),ta(3,6), ta(4,1),ta(4,2),ta(4,3),ta(4,4), &
           ta(4,5),ta(4,6), ta(5,1),ta(5,2),ta(5,3),ta(5,4),ta(5,5),ta(5,6), &
           ta(6,1),ta(6,2),ta(6,3),ta(6,4),ta(6,5),ta(6,6), dmmac,dnms,dizu0,&
           dnumlr,sigcor,dpscor

      iPair = (ifipa-1)/2+1

      if (iPair == pairIdx) then
        write(outUnit,err=620,iostat=iostat) &
             sixtit,commen,cdate,ctime,       &
             progrm,ifipa,ilapa,itopa,icode,numl,qwc(1),qwc(2),qwc(3), clo(1), &
             clop(1),clo(2),clop(2),clo(3),clop(3), di0(1),dip0(1),di0(2),dip0(2), &
             dummy,dummy, ta(1,1),ta(1,2),ta(1,3),ta(1,4),ta(1,5),ta(1,6), &
             ta(2,1),ta(2,2),ta(2,3),ta(2,4),ta(2,5),ta(2,6), ta(3,1),ta(3,2), &
             ta(3,3),ta(3,4),ta(3,5),ta(3,6), ta(4,1),ta(4,2),ta(4,3),ta(4,4), &
             ta(4,5),ta(4,6), ta(5,1),ta(5,2),ta(5,3),ta(5,4),ta(5,5),ta(5,6), &
             ta(6,1),ta(6,2),ta(6,3),ta(6,4),ta(6,5),ta(6,6), dmmac,dnms,dizu0,&
             dnumlr,sigcor,dpscor

        ntwin = ilapa-ifipa+1
        if (.not. (ntwin == 1 .or. ntwin == 2)) then
          write(textOutUnit,'(a)' ) "Error when reading header; ifipa=",ifipa, "ilapa=",ilapa,"ntwin=",ntwin
          stop 1
        end if

      end if

      if (ifipa >= itopa-1) then
        exit !Got all the headers
      end if
    end do headers

    if(ntwin == 1) then
      data1: do
        read(inUnit,end=530,err=531, iostat=iostat) &
             ia,ifipa,b,c,d,e,f,g,h,p
        iPair = (ia-1)/2+1

        if (iPair == pairIdx) then
          write(outUnit, err=621, iostat=iostat) &
               ia,ifipa,b,c,d,e,f,g,h,p
        end if
      end do data1
    else if (ntwin == 2) then
      data2: do
        read(inUnit,end=530,err=531, iostat=iostat) &
             ia,ifipa,b,c,d,e,f,g,h,p, ilapa,b,c1,d1,e1,f1,g1,h1,p1
        iPair = (ifipa-1)/2+1

        if (iPair == pairIdx) then
          write(outUnit, err=621, iostat=iostat) &
               ia,ifipa,b,c,d,e,f,g,h,p, ilapa,b,c1,d1,e1,f1,g1,h1,p1
        end if
      end do data2
    end if

    !Should never reach this point
    write(textOutUnit,'(a)') "Control reached an unreachable point?!?"
    stop 1

    !Error handling
410 continue
    write(textOutUnit,'(a,i6)') "ERR while opening input file; iostat=",iostat
    stop 1

411 continue
    write(textOutUnit,'(a,i6)') "ERR while opening output file; iostat=",iostat
    stop 1

511 continue
    write(textOutUnit,'(a,i6)') "END while reading header; iostat=",iostat
    stop 1

520 continue
    write(textOutUnit,'(a,i6)') "ERR while reading header; iostat=",iostat
    stop 1

530 continue
    !END while reading DATA
    close (inUnit)
    close (outUnit)
    return

531 continue
    write(textOutUnit,'(a,i6)') "ERR while reading data ; iostat=",iostat
    stop 1

620 continue
    write(textOutUnit,'(a,i6)') "ERR while writing header; iostat=",iostat
    stop 1

621 continue
    write(textOutUnit,'(a,i6)') "ERR while writing data; iostat=",iostat
    stop 1

  end subroutine convertSTFpair

  subroutine convertSTFfile(ifname, oldnames)
    character(len=*), intent(in) :: ifname   !Name of the singletrackfile to split
    logical, intent(in)          :: oldnames !Use the fort.91-pairIdx names?
                                             !Only allowed if <= 32 pairs

    integer nPairs
    integer iPair
    character(len=100) ofnameBuffer !TODO: More inteligent scaling (or at least a check) would be nice
                                    ! However the actual names are hard coded so whatever.

    nPairs = numSTFpairs(ifname)

    if (oldnames .and. nPairs > 32) then
      write(textOutUnit,'(a,i6,a)') "ERR in convertSTFfile: Oldpairs=true but nPairs = ", npairs, " > 32"
      stop 1
    end if

    pairs: do iPair = 1,nPairs
      if (oldnames) then
        write(ofnameBuffer,'(a,i0)')     "fort.",91-iPair
      else
        write(ofnameBuffer,'(a,a,i0.6)') ifname,'.',iPair
      end if

      call convertSTFpair(iPair, ifname, ofnameBuffer)

    end do pairs

  end subroutine convertSTFfile

end module mod_splitSingletrack
