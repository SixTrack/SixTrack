! =================================================================================================
!  END SIXTRACK
!  Last modified: 2018-06-25
! =================================================================================================

subroutine prror(ier)

  use crcoall
  use parpro,     only : mcor,mmul,mran,nbb,nran,nrco, npart, nelb, nele, nblo, nzfz, ntr, nper, nema, mNameLen
  use mod_common, only : ierro,errout_status, erbez

#ifdef FLUKA
  use mod_fluka
#endif

  implicit none

  integer ier

#ifdef FLUKA
  integer fluka_con
#endif

  save

  errout_status = ier

  write(lout,"(a)") "    +++++++++++++++++++++++++++++"
  write(lout,"(a)") "    +      ERROR DETECTED       +"
  write(lout,"(a)") "    + RUN TERMINATED ABNORMALLY +"
  if(ier > 0) then
    write(lout,"(a,i3,a)") "    +      ERROR CODE: ",ier,"      +"
  end if
  write(lout,"(a)") "    +++++++++++++++++++++++++++++"

  select case (ier)
   case (1)
     write(lout,10010)
   case (2)
     write(lout,10020) nele
   case (3)
     write(lout,10030)
  case (4)
    write(lout,10040)
  case (5)
    write(lout,10050)
   case (6)
     write(lout,10060)
   case (7)
     write(lout,10070)
  case (8)
    write(lout,10080)
  case (9)
    write(lout,10090)
   case (10)
     write(lout,10100)
  case (11)
    write(lout,10110)
  case (12)
    write(lout,10120)
  case (13)
    write(lout,10130)
  case (14)
    write(lout,10140)
   case (15)
     write(lout,10150)
   case (16)
     write(lout,10160) nele
   case (17)
     write(lout,10170) nper
   case (18)
     write(lout,10180) nblo
   case (19)
     write(lout,10190) erbez
   case (20)
     write(lout,10200) erbez
   case (21)
     write(lout,10210)
   case (22)
     write(lout,10220)
   case (23)
     write(lout,10230)
   case (24)
     write(lout,10240)
   case (25)
     write(lout,10250)
   case (26)
     write(lout,10260) nelb
   case (27)
     write(lout,10270)
   case (28)
     write(lout,10280)
   case (29)
     write(lout,10290)
  case (30)
    write(lout,10300) nran
  case (31)
    write(lout,10310)
  case (32)
    write(lout,10320)
  case (33)
    write(lout,10330)
  case (34)
    write(lout,10340) mran
  case (35)
    write(lout,10350)
  case (36)
    write(lout,10360)
   case (37)
     write(lout,10370)
  case (38)
    write(lout,10380)
   case (39)
     write(lout,10390)
   case (40)
     write(lout,10400)
   case (41)
     write(lout,10410)
   case (42)
     write(lout,10420)
   case (43)
     write(lout,10430) nzfz
   case (44)
     write(lout,10440)
   case (45)
     write(lout,10450)
   case (46)
     write(lout,10460) nrco
   case (47)
     write(lout,10470)
   case (48)
     write(lout,10480)
   case (49)
     write(lout,10490)
   case (50)
     write(lout,10500)
   case (51)
     write(lout,10510)
   case (52)
     write(lout,10520) nema
   case (53)
     write(lout,10530)
   case (54)
     write(lout,10540) npart
!   case (55)
!     write(lout,10550) nmac
  case (56)
    write(lout,10560) ierro
   case (57)
     write(lout,10570) ierro
  case (58)
    write(lout,10580) ierro
   case (59)
     write(lout,10590) ierro
   case (60)
     write(lout,10600) ierro
   case (61)
     write(lout,10610) ierro
   case (62)
     write(lout,10620)
   case (63)
     write(lout,10630)
  case (64)
    write(lout,10640)
   case (65)
     write(lout,10650) mcor
   case (66)
     write(lout,10660)
   case (67)
     write(lout,10670)
   case (68)
     write(lout,10680)
   case (69)
     write(lout,10690)
   case (70)
     write(lout,10700)
   case (71)
     write(lout,10710)
   case (72)
     write(lout,10720)
   case (73)
     write(lout,10730)
   case (74)
     write(lout,10740)
   case (75)
     write(lout,10750)
   case (76)
     write(lout,10760)
   case (77)
     write(lout,10770)
   case (78)
     write(lout,10780)
  case (79)
    write(lout,10790)
   case (80)
     write(lout,10800)
   case (81)
     write(lout,10810)
   case (82)
     write(lout,10820)
   case (83)
     write(lout,10830)
  case (84)
    write(lout,10840)
   case (85)
     write(lout,10850) mmul
   case (86)
     write(lout,10860)
   case (87)
     write(lout,10870)
  case (88)
    write(lout,10880)
  case (89)
    write(lout,10890)
  case (90)
    write(lout,10900)
  case (91)
    write(lout,10910)
  case (92)
    write(lout,10920)
  case (93)
    write(lout,10930)
  case (94)
    write(lout,10940)
  case (95)
    write(lout,10950)
  case (96)
    write(lout,10960)
  case (97)
    write(lout,10970)
   case (98)
     write(lout,10980)
   case (99)
     write(lout,10990)
   case (100)
     write(lout,11000) ntr
  case (101)
    write(lout,11010)
  case (102)
    write(lout,11020) nbb
  case (103)
    write(lout,11030)
   case (104)
     write(lout,11040) mNameLen
   case (105)
     write(lout,11050) mmul
  end select

#ifdef FLUKA
  call fluka_close
#endif

#ifdef CR
  call abend('                                                  ')
#else
  call closeUnits
  stop 1
#endif

10010 format(t4,'WRONG MODE DEFINITION')
10020 format(t4,'NOT MORE THAN ',i4,' POSITIONS FOR RESONANCE-COMPENSATION ALLOWED')
10030 format(t4,'ELEMENT FOR RESONANCE-COMPENSATION IS NOT IN THE ELEMENT LIST')
10040 format(t4,'UNSTABLE CLOSED ORBIT DURING INITIAL DISPERSION CALCULATION'/&
             t4,'INSTABILITY OCCURRED FOR SMALL RELATIVE ENERGY DEVIATION')
10050 format(t4,'UNSTABLE CLOSED ORBIT FOR ZERO ENERGY DEVIATION')
10060 format(t4,'UNSTABLE CLOSED ORBIT DURING DISPERSION CALCULATION AFTER ORBIT SCALING'/&
             t4,'INSTABILITY OCCURRED FOR SMALL RELATIVE ENERGY DEVIATION')
10070 format(t4,'UNSTABLE CLOSED ORBIT AFTER ORBIT SCALING')
10080 format(t4,'ELEMENTS SPECIFIED FOR TUNE VARIATION ARE NOT QUADRUPOLES')
10090 format(t4,'UNSTABLE CLOSED ORBIT DURING TUNE VARIATION')
10100 format(t4,'NO OPTICAL SOLUTION DURING TUNE VARIATION')
10110 format(t4,'ELEMENTS SPECIFIED FOR CHROMATICITY CORRECTION ARE NOT SEXTUPOLES')
10120 format(t4,'UNSTABLE CLOSED ORBIT DURING CHROMATICITY CORRECTION')
10130 format(t4,'NO OPTICAL SOLUTION DURING CHROMATICITY CORRECTION')
10140 format(t4,'ELEMENTS OF DIFFERENT TYPES ARE COMBINED IN DATA BLOCK COMBINATION OF ELEMENTS')
10150 format(t4,'UNKNOWN BLOCK SPECIFICATION')
10160 format(t4,'NO. OF SINGLE ELEMENTS EXCEEDS THE MAXIMUM ALLOWED VALUE: ',i4)
10170 format(t4,'NO. OF SUPER PERIODS LARGER THAN: ',i4)
10180 format(t4,'NO. OF DIFFERENT BLOCKS EXCEEDS THE MAXIMUM ALLOWED VALUE: ',i5)
10190 format(t4,'UNKNOWN SINGLE ELEMENT ',a16,' IN THE BLOCK DEFINITION')
10200 format(t4,'UNKNOWN BLOCK NAME OR INSERTION NAME ',a16,' IN THE STRUCTURE INPUT')
10210 format(t4,'MAXIMUM NUMBER OF STRUCTURE ELEMENTS SURPASSED')
10220 format(t4,'NO SOLUTION FOR ORBIT SCALING - POSSIBLE REASONS:'/&
             t4,'--> DIPOLE STRENGTHS OF NON-CORRECTOR ELEMENTS TO HIGH'/&
             t4,'--> NON-LINEARITIES TOO STRONG, TRY TO INCREASE INITIAL CORRECTOR STRENGTHS'/&
             t4,'--> USE ALL DIPOLE ELEMENTS FOR SCALING'/)
10230 format(t4,'NO OPTICAL SOLUTION')
10240 format(t4,'NO SOLUTION FOR DISPERSION')
10250 format(t4,'--> PLEASE INCLUDE LENGTH OF MACHINE IN THE <SYNCHROTRON> BLOCK')
10260 format(t4,'ONE BLOCK CAN NOT HAVE MORE THAN ',i4,' ELEMENTS')
10270 format(t4,'KINETIC ENERGY OF THE PARTICLE IS LESS OR EQUAL TO ZERO')
10280 format(t4,'EITHER YOUR RF FREQUENCY IS SHIFTED BY 180 DEGREES,'/,&
             t4,'THEN CHANGE THE SIGN OF <ITION> IN THE <SYNCHROTRON>-INPUTBLOCK,',/&
             t4,'OR YOUR ALFA-P IS WRONGLY INTRODUCED IN THE SAME INPUTBLOCK')
10290 format(t4,'MULTIPOLE COEFFICIENTS CANNOT BE SET EQUAL')
10300 format(t4,'THE RANDOM NUMBER: ',i15,' FOR THE INITIAL STRUCTURE IS TOO SMALL')
10310 format(t4,'ELEMENTS THAT NEED RANDOM NUMBERS HAVE A KZ > 0')
10320 format(t4,'THERE ARE NOT ENOUGH RANDOM NUMBERS FOR THE INSERTED ELEMENTS')
10330 format(t4,'TO USE THE SAME RANDOM NUMBERS FOR 2 ELEMENTS, THE INSERTED ELEMENT ',&
                 'MUST NOT NEED MORE OF SUCH NUMBERS THAN THE REFERENCE ELEMENT')
10340 format(t4,'NOT MORE THAN',i4,' OF EACH TYPE OF INSERTED ELEMENTS CAN BE USED')
10350 format(t4,'PROBLEMS DURING MATRIX-INVERSION IN QMOD')
10360 format(t4,'NO CONVERGENCE IN RMOD')
10370 format(t4,'CHOSEN ORDERS OF RESONANCES CAN NOT BE CALCULATED')
10380 format(t4,'PROBLEMS DURING MATRIX-INVERSION IN RMOD')
10390 format(t4,'WITH THE SPECIFIED ELEMENTS THE RESONANCE CANNOT BE COMPENSATED'/&
             t4,'RESONANCEORDER AND ELEMENTTYP # MUST BE THE SAME')
10400 format(t4,'NOT MORE THAN 2 PARTICLES CAN BE TRACKED')
10410 format(t4,'GEOMETRY AND STRENGTH FILE (UNIT 2) IS EMPTY OR DESTROYED')
10420 format(t4,'TRACKING PARAMETER FILE (UNIT 3) IS EMPTY OR NON-EXISTING')
10430 format(t4,'NOT MORE THAN ',i4,' RANDOM NUMBERS CAN BE USED')
10440 format(t4,'FOR THE INPUTBLOCK - ORBIT CORRECTION - ONLY CORRECTORS WITH THE KEYWORDS ( HCOR= ; VCOR= )'/&
             t4,'AND MONITORS WITH THE KEYWORDS ( HMON= ; VMON= ) ARE ALLOWED')
10450 format(t4,'FOR THE INPUTBLOCK - LINEAR OPTICS - ONLY THE KEYWORD ( ELEMENT ) AND ( BLOCK ) ARE ALLOWED')
10460 format(t4,'ORDER OF COMPENSATION CAN NOT BE LARGER THAN: ',i4)
10470 format(t4,'ONLY UP TO 3 RESONANCES CAN BE COMPENSATED')
10480 format(t4,'RESONANCE TYPE IS OUT OF THE RANGE OF THE RESONANCE ORDER')
10490 format(t4,'ONLY UP TO 3 SUB-RESONANCES CAN BE COMPENSATED')
10500 format(t4,'THE MULTIPOLE ORDER FOR THE SUB-RESONANCE COMPENSATION SHOULD NOT EXCEED THE VALUE 9')
10510 format(t4,'PROBLEMS WITH FILE 3 WITH DYNAMIC KICKS')
10520 format(t4,'MAXIMUM ORDER OF THE ONE TURN MAP IS ',i4,/,' NEMA HAS TO BE LARGER THAN NORD')
10530 format(t4,'# OF VARIABLES -NV- OF THE ONE TURN MAP IS NOT IN THE ALLOWED RANGE [2 <= NV <= 5]')
10540 format(t4,'MAXIMUM NUMBER OF PARTICLES FOR VECTORIZATION IS ',i7)
10550 format(t4,'MAXIMUM NUMBER OF DIFFERENT SEEDS FOR VECTORIZATION IS ',i4)
10560 format(t4,'PROBLEMS WITH FILE 13 WITH INITIAL COORDINATES - ERROR CODE: ',i10)
10570 format(t4,'PROBLEMS WITH FILE 2 WITH ACCELERATOR STRUCTURE - ERROR CODE: ',i10)
10580 format(t4,'PROBLEMS WITH FILE 3 WITH TRACKING PARAMETERS - ERROR CODE: ',i10)
10590 format(t4,'PROBLEMS WITH FILE 11 FOR CRAY INPUT - ERROR CODE: ',i10)
10600 format(t4,'PROBLEMS WITH FILE 99 FOR BINARY OUTPUT - ERROR CODE: ',i10)
10610 format(t4,'PROBLEMS WITH FILE 12 FOR END COORDINATES - ERROR CODE: ',i10)
10620 format(t4,'ELEMENTS SPECIFIED FOR DECOUPLING ARE NOT SKEW QUADRUPOLES')
10630 format(t4,'THERE ARE THE APPROPRIATE ELEMENTS FOR THE DECOUPLING OR SIMULTANEOUS TUNE ADJUSTMENT')
10640 format(t4,'PROBLEMS DURING MATRIX-INVERSION IN DECOUP')
10650 format(t4,'MAXIMUM NUMBER OF EXTRA PARAMETERS IS: ',i4)
10660 format(t4,'EXTRA PARAMETERS FOR THE MAP DOES NOT EXIST')
10670 format(t4,'ONLY SINGLE KICK ELEMENTS ALLOWED FOR MAP CALCULATION')
10680 format(t4,'THE ORDER OF THE NORMAL FORM IS TOO HIGH. CHECK THE DIFFERENTIAL ALGEBRA PARAMETERS')
10690 format(t4,'TOO MANY VARIABLES SPECIFIED. CHECK THE DIFFERENTIAL ALGEBRA PARAMETERS')
10700 format(t4,'NO CORRECTORS SPECIFIED')
10710 format(t4,'BOTH AMPLITUDE AND MOMENTUM ORDER ARE ZERO!')
10720 format(t4,'BOTH AMPLITUDE AND MOMENTUM ORDER ARE DIFFERENT FROM ZERO!')
10730 format(t4,'AMPLITUDE ORDER OUTSIDE RANGE [0,2]')
10740 format(t4,'MOMENTUM ORDER OUTSIDE RANGE [0,3] (ONE EXCLUDED!)')
10750 format(t4,'MINIMUM ORDER OUTSIDE RANGE [2,3]')
10760 format(t4,'MINIMUM ORDER GREATER THAN MAXIMUM!')
10770 format(t4,'MAXIMUM ORDER OUTSIDE RANGE [2,3]')
10780 format(t4,'NORMAL FORMS ANALYSIS IMPOSSIBLE',/,t4,'THE TRANSFER MAP DOES NOT EXIST!')
10790 format(t4,'ZERO OR NEGATIVE ENERGY DOES NOT MAKE MUCH SENSE')
10800 format(t4,'PROBLEM READING EXTERNAL MULTIPOLE ERRORS')
10810 format(t4,'TOO MANY ELEMENTS FOR LINEAR OPTICS WRITE-OUT')
10820 format(t4,'FOR CLOSED ORBIT CORRECTORS ONLY DIPOLES OF LEGTH ZERO OR MULTIPOLE LENSES ALLOWED')
10830 format(t4,'AN ELEMENT FOR CLOSED ORBIT CORRECTION CAN BE ONLY EITHER A HORIZONTAL MONITOR',/,&
             t4,'OR A VERTICAL MONITOR OR A HORIZONTAL CORRECTOR OR A VERTICAL CORRECTOR')
10840 format(t4,'NUMBER OF ORBIT CORRECTORS IS ZERO')
10850 format(t4,'THE ORDER OF MULTIPOLES MMUL: ',i4,' HAS TO BE LARGER THAN 10 BUT SMALLER THAN 20')
10860 format(t4,'PROBLEM READING EXTERNAL MISALIGNMENTS')
10870 format(t4,'PROBLEM READING FROM FILE 30 (SINGLE KICKS AND MISALIGNMENTS')
10880 format(t4,'BEAM_BEAM: EITHER NORMALIZED EMITTANCES OR THE RESULTING SIGMA VALUES EQUAL TO ZERO')
10890 format(t4,'BEAM_BEAM: AT EACH INTERACTION POINT THE BEAM MUST BE EITHER ROUND OR ELLIPTICAL FOR ALL PARTICLES')
10900 format(t4,'QUADRUPOLES ARE NOT SUITED TO ADJUST THE TUNES')
10910 format(t4,'ORDER AND NUMBER OF VARIABLES HAVE TO BE LARGER THAN ZERO TO CALCULATE A DIFFERENTIAL ALGEBRA MAP')
10920 format(t4,'YOU CANNOT COMBINE AN ELEMENT WITH ITSELF')
10930 format(t4,'INVERTED LINEAR BLOCKS NOT ALLOWED')
10940 format(t4,'  NUMBER OF NORMAL FORM VARIABLES HAVE TO BE: 2, 4, 5, 6 + PARAMETERS')
10950 format(t4,'  DA CORRECTIONS IMPLEMENTED FOR 4-D AND 6-D ONLY')
10960 format(t4,'SEXTUPOLES ARE NOT SUITED TO ADJUST THE CHROMATICITY')
10970 format(t4,'UNSTABLE CLOSED ORBIT IN DA CALCULATION')
10980 format(t4,'TROMBONE ELEMENT NOT IN LIST OF SINGLE ELEMENTS')
10990 format(t4,'INCOMPLETE PARAMETERS FOR TROMBONE ELEMENT')
11000 format(t4,'MAXIMUM NUMBER OF TROMBONES EXCEEDED : NTR = ',i4)
11010 format(t4,'AMPLITUDES EXCEED THE MAXIMUM VALUES IN UMLAUF')
11020 format(t4,'MAXIMUM ELEMENT NUMBER FOR BEAM_BEAM WITH COUPLING EXCEEDED:  NBB = ',i4)
11030 format(t4,'6D BEAM-BEAM WITH TILT NOT POSSIBLE')
11040 format(t4,'SINGLE ELEMENT NAME LONGER THAN ', i3, ' CHARACTERS')
11050 format(t4,'THE INPUT ORDER OF MULTIPOLES IS LARGER THAN THE MAXIMUM ALLOWED ORDER MMUL: ',i4)

end subroutine prror

#ifdef CR
! =================================================================================================
subroutine abend(cstring)

  use numerical_constants

  use, intrinsic :: iso_fortran_env, only : error_unit

  use crcoall
  use parpro
  use mod_common
  use checkpoint_restart
  use string_tools
  use mod_units
  use mod_version
  use mod_time

  implicit none

  integer i,lstring,j,itot,ttot
  character(len=*) cstring
  character(len=256) filename
  real(kind=fPrec) sumda(60)
  logical fopen, rErr, fErr
  character(len=8192) ch
  character(len=25) ch1
  integer errno,l1,l2
  integer ich

  save

  write(93,*) 'SIXTRACR STOP/ABEND called and closing files'
  endfile(93,iostat=ierro)
  backspace(93,iostat=ierro)

  ! Calling close to be very safe.......96 calls to abend
  ! Easier than adding the call on every abend
  call closeUnits

  ! If fort.10 is inexistent (physics error or some other problem)
  ! we try and write a 0d0 file with a turn number and CPU time
  write(93,"(a)") 'SIXTRACR STOP/ABEND checking fort.10'
  endfile(93,iostat=ierro)
  backspace(93,iostat=ierro)

  call units_openUnit(unit=10,fileName="fort.10",formatted=.true.,mode="rw",err=fErr,recl=8195)
  if(fErr) goto 11

  ! Now we try and read fort.10 i.e. is it empty?
  read(10,'(a1024)',end=11,err=11,iostat=ierro) arecord
  ! Seems to be OK
  goto 12

11 continue
  ! Now we try and write a fort.10
  ! We put some CPU for Igor, a version, and turn number 0
  write(93,*) 'SIXTRACR STOP/ABEND writing a fort.10'
  endfile(93,iostat=ierro)
  backspace(93,iostat=ierro)

  ! The version
  itot=0
  ttot=0
  do i=1,8
    if (version(i:i).ne.' ') then
      if (version(i:i).ne.'.') then
        itot=itot*10+ichar(version(i:i))-ichar('0')
      else
        ttot=ttot*10**2+itot
        itot=0
      end if
    end if
  end do
  ttot=ttot*10**2+itot
  do i=1,60
    sumda(i)=zero
  end do
  sumda(52)=real(ttot,fPrec)
  ! The CPU
  call time_timerCheck(time1)
  trtime=time1-time0
#ifdef CR
  trtime=trtime+crtime3
#endif
  sumda(60)=real(trtime,fPrec)
  ! Note it COULD happen that napxo is 0 for a very very early error and even napx!!!
  if (napxo.eq.0.and.napx.eq.0) napxo=1
  write(93,*) 'SIXTRACR STOP/ABEND writing fort.10, lines',napxo,'/',napx
  endfile(93,iostat=ierro)
  backspace(93,iostat=ierro)
  if (napxo.eq.0.and.napx.eq.0) napxo=1
  if (napxo.eq.0) napxo=napx
  do j=1,napxo
#ifndef CRLIBM
    write(ch,*,iostat=ierro) (sumda(i),i=1,60)
    do ich=8192,1,-1
      if(ch(ich:ich).ne.' ') goto 707
    end do
707 write(10,'(a)',iostat=ierro) ch(:ich)
#else
    l1=1
    do i=1,60
      ! We return the length of the string (always 24)
      call chr_fromReal(sumda(i),ch1,19,2,rErr)
      ch(l1:l1+25)=' '//ch1(1:25)
      l1=l1+26
    end do
    write(10,'(a)',iostat=ierro) ch(1:l1-1)
#endif
    if(ierro.ne.0) then
      write(lout,"(a,i0)") "ABEND> ERROR Problems writing to fort.10. ierro = ",ierro
    end if
  end do

12 continue
  close (10,iostat=ierro)
#ifdef CR
#ifdef DEBUG
  ! call system('../crend   >> crlog')
#endif
  close(91,err=4)
4 continue
  close(94,err=5)
5 continue
  close(95,err=6)
6 continue
  close(96,err=7)
7 continue
  if (lout.eq.92) then
    write(93,*) 'SIXTRACR STOP/ABEND copying fort.92 to fort.6'
    endfile(93,iostat=ierro)
    backspace(93,iostat=ierro)
    rewind 92
3   read(92,'(a1024)',end=1,err=8,iostat=ierro) arecord
    lstring=1024
    do i=1024,2,-1
      lstring=i
      if (arecord(i:i).ne.' ') goto 2
      lstring=lstring-1
    end do
2   write(6,'(a)',iostat=ierro) arecord(1:lstring)
    goto 3
  end if

1 write(6,"(a)",iostat=ierro) 'SIXTRACR> Stop '//cstring
  close(6,iostat=ierro)
  ! and get rid of fort.92
  rewind 92
  endfile(92,iostat=ierro)
  close(92)
  write(93,*) 'SIXTRACR stop '//cstring
  write(93,*)
#ifdef DEBUG
  ! call system('../crend   >> crlog')
#endif
#ifdef BOINC
  do i=2,120
    inquire(i,opened=fopen)
    write(93,*) 'UNIT ',i,' opened ',fopen
  end do
  ! call boinc_zipitall()
  ! call boinc_finish_graphics()
  if(errout_status.ne.0) then
    close(93)
    call boincrf('fort.93',filename)
    call print_lastlines_to_stderr(93,filename)
    call boincrf('fort.6',filename)
    call print_lastlines_to_stderr(6,filename)
  end if
  call boinc_finish(errout_status) !This call does not return
#else
  if(errout_status.ne.0) then
    close(93)
    call print_lastlines_to_stderr(93,"fort.93")
    call print_lastlines_to_stderr(6,"fort.6")

    write(error_unit,'(a,i5)') "Stopping, errout_status=",errout_status
    stop 1
  else
    ! No error
    stop
  end if
#endif

!!!!!! In case of errors when copying fort.92 (lout) -> fort.6 !!!!!!!!!

8 write(93,*) 'SIXTRACR CR ABEND *** ERROR *** reading fort.92, iostat=',ierro
  close(93)
  write(6,*) 'SIXTRACR CR ABEND *** ERROR *** reading fort.92, iostat=',ierro
#ifdef DEBUG
  ! call system('../crend   >> crlog')
#endif
#ifdef BOINC
  do i=2,120
    inquire(i,opened=fopen)
    write(6,*) 'UNIT ',i,' opened ',fopen
  end do
  close(6,err=31)
31 continue
  ! call boinc_zipitall()
  ! call boinc_finish_graphics()
  if(errout_status.ne.0) then
    close(93)
    call boincrf('fort.93',filename)
    call print_lastlines_to_stderr(93,filename)
    call boincrf('fort.6',filename)
    call print_lastlines_to_stderr(6,filename)
  end if
  call boinc_finish(errout_status) !This call does not return
#else
  if(errout_status.ne.0) then
    close(93)
    call print_lastlines_to_stderr(93,"fort.93")
    call print_lastlines_to_stderr(6,"fort.6")

    write(error_unit,'(a,i5)') "Stopping, errout_status=",errout_status
    stop 1
  else
    ! No error
    stop
  end if
#endif
#else
  ! This one should probably remain as write(*,*) or use output_unit
  write(*,*) 'SIXTRACK STOP/ABEND '//cstring
#ifdef DEBUG
  ! call system('../crend   >> crlog')
#endif
  ! No fort.6 and 93 if not CR -> don't do print_lastlines_to_stderr()
  if(errout_status.ne.0) then
    write(error_unit,'(a,i5)') "Stopping, errout_status=",errout_status
    stop 1
  else
    ! No error
    stop
  end if
#endif

end subroutine abend
#endif

! =================================================================================================
!  K. Sjobak, BE-ABP-HSS
!  Last modified: June 2017
!  Subroutine to copy the last nlines lines from a file to stderr
!  It is mainly used just before exiting SixTrack in case there was an error.
!  This is useful since STDERR is often returned from batch systems and BOINC.
! =================================================================================================
subroutine print_lastlines_to_stderr(file_unit, file_name)

  use, intrinsic :: iso_fortran_env, only : error_unit
  use string_tools

  implicit none

  integer,          intent(in) :: file_unit
  character(len=*), intent(in) :: file_name

  integer, parameter :: nlines = 40

  character(len=1024) fileBuff (nlines)
  integer fileBuff_idx
  integer i,j
  integer printLines

  logical lopen
  integer ierro

  ! Clear the buffer
  do i=1,nlines
    fileBuff(i)=''
  end do

  ! Open the file
  inquire(unit=file_unit,opened=lopen)
  if(lopen) then
    write(error_unit,'(a,1x,i5,1x,a)') "Error when opening unit #",file_unit,": The file is already open."
    return
  end if

  open(file_unit,file=file_name,form="formatted",status="old",iostat=ierro)
  if(ierro .ne. 0) then
    write(error_unit,'(a,a,a,1x,i5,1x,a,1x,i5)') &
      "Error when opening file '",trim(file_name),"' on unit #", file_unit,", iostat =",ierro
    return
  end if

  ! Read into fileBuff as a ring buffer.
  fileBuff_idx = 1
  j = 0

1 read(file_unit,'(a1024)',end=3,err=2,iostat=ierro) fileBuff(fileBuff_idx)
  ! write(error_unit,*) fileBuff_idx,":",trim(fileBuff(fileBuff_idx))
  fileBuff_idx = fileBuff_idx+1
  if (fileBuff_idx.ge.nlines) fileBuff_idx = 1
  j = j+1
  goto 1

2 continue ! An error occured
  write(error_unit,'(a,1x,i5)') "An error occured while reading the file, iostat=",ierro


3 continue ! EOF or error; anyway we're done.
  close(file_unit)

  ! Print stuff back out from the buffer
  if (j .lt. nlines) then
    printLines = j
    fileBuff_idx=1
  else
    printLines = nlines
  end if
  write(error_unit,'(a,1x,i5,1x,a,a,a)') "******* Last",printLines,"lines of file '",trim(file_name),"': *******"

  i = fileBuff_idx  ! Position in buffer (we have already incremented i)
  j = 0             ! How many have we printed

10 continue
  if(i.gt.nlines) i=1      ! i is wrapping
    write(error_unit,'(a)') trim(fileBuff(i))
    i = i+1
    j = j+1                   ! j just counts
    if(j.lt.printLines) goto 10
      write(error_unit,'(a,a,a)') "******* Done writing tail of file '",trim(file_name),"' to stderr *******"

end subroutine print_lastlines_to_stderr
