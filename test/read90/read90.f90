      program read90
      use, intrinsic :: iso_fortran_env, only : stdout=>output_unit
!-----------------------------------------------------------------------
!  Read a fort.90 and print with correct binary/decimal conversion
!-----------------------------------------------------------------------
      implicit none
      integer errno
      character(len=80) sixtit,commen
      integer icode
      double precision dpscor,sigcor
      integer i,ia,ia0,ifipa,ilapa,itopa,j,nfile,ntwin,numl,n
      double precision b,c,c1,clo,clop,d,d1,di0,dip0,dizu0,dmmac,dnms,dnumlr,  &
           dummy,e,e1,f,f1,g,g1,h,h1,p,p1,qwc,ta
      character(len=8) cdate,ctime,progrm
      character(len=25) ch1
      dimension qwc(3),clo(3),clop(3),di0(2),dip0(2)
      dimension ta(6,6)
      
      interface
         integer function dtostr(x,results)
           double precision, intent(in) :: x
           character(*), intent(out) :: results
         end function dtostr
      end interface
      
      integer :: stat
      
      INTEGER :: cmdarg_i, cmdarg_length, cmdarg_status
      CHARACTER(len=100) :: cmdarg_arg
      logical STF
      integer singleParticle, firstParticle, lastParticle

      character(len=100) :: fname
      character(len=100) :: ofname
      logical ofoutput
      logical hasInputFile
      integer ounit
!----------------------------------------------------------------------
      ofoutput = .FALSE.
      ounit = stdout
      fname = "fort.190"
      
      STF = .false. !Use the new STF format?
      singleParticle = -1 !If set, only print records for particle i (and i+1)
      firstParticle = -1
      lastParticle  = -1
      cmdarg_i = 0
      do 
         call get_command_argument(cmdarg_i, cmdarg_arg,cmdarg_length,cmdarg_status)
         if (len_trim(cmdarg_arg)==0) EXIT

         if (cmdarg_i.gt.0) then !Skip first argument (command name)
            if (cmdarg_arg .eq. "--STF") then
               STF=.true.
            else if (cmdarg_arg .eq. "--SP") then
               if(.not.STF) then
                  write(*,*) "--SP flag only valid following a --STF flag."
                  flush(stdout)
                  stop 2
               endif
               if(firstParticle.ne.-1 .or. lastParticle.ne.-1) then
                  write(*,*) "--SP and --PR flags are incompatible."
                  flush(stdout)
                  stop 8
               endif
               
               !Read the number
               cmdarg_i = cmdarg_i+1
               call get_command_argument(cmdarg_i, cmdarg_arg,cmdarg_length,cmdarg_status)
               if (len_trim(cmdarg_arg)==0) then
                  write(*,*) "No number found following --SP flag?"
                  flush(stdout)
                  stop 3
               endif

               read(cmdarg_arg,*,iostat=stat) singleParticle
               if (singleParticle.lt.1) then
                  write(*,*) "singleParticle=",singleParticle
                  write(*,*) "Did you specify an integer>0?"
                  flush(stdout)
                  stop 4
               endif
               if (mod(singleParticle,2).ne.1) then
                  write(*,*) "singleParticle=",singleParticle, "; expected odd number."
                  flush(stdout)
                  stop 5
               endif
            else if (cmdarg_arg .eq. "--PR") then
               if(.not.STF) then
                  write(*,*) "--PR flag only valid following a --STF flag."
                  flush(stdout)
                  stop 9
               endif
               if(singleParticle.ne.-1) then
                  write(*,*) "--SP and --PR flags are incompatible."
                  flush(stdout)
                  stop 10
               endif

               !Read the first number
               cmdarg_i = cmdarg_i+1
               call get_command_argument(cmdarg_i, cmdarg_arg,cmdarg_length,cmdarg_status)
               if (len_trim(cmdarg_arg)==0) then
                  write(*,*) "No number found following --PR flag?"
                  flush(stdout)
                  stop 11
               endif

               read(cmdarg_arg,*,iostat=stat) firstParticle
               if (firstParticle.lt.1) then
                  write(*,*) "firstParticle=",firstParticle
                  write(*,*) "Did you specify an integer>0?"
                  flush(stdout)
                  stop 12
               endif
               if (mod(firstParticle,2).ne.1) then
                  write(*,*) "firstParticle=",firstParticle, "; expected odd number."
                  flush(stdout)
                  stop 13
               endif

               !Read the second number
               cmdarg_i = cmdarg_i+1
               call get_command_argument(cmdarg_i, cmdarg_arg,cmdarg_length,cmdarg_status)
               if (len_trim(cmdarg_arg)==0) then
                  write(*,*) "No second number found following --PR flag?"
                  flush(stdout)
                  stop 14
               endif

               read(cmdarg_arg,*,iostat=stat) lastParticle
               if (lastParticle.lt.1) then
                  write(*,*) "lastParticle=",lastParticle
                  write(*,*) "Did you specify an integer>0?"
                  flush(stdout)
                  stop 15
               endif
               if (mod(lastParticle,2).ne.1) then
                  write(*,*) "lastParticle=",lastParticle, "; expected odd number."
                  flush(stdout)
                  stop 16
               endif
               if (.not. lastParticle .gt. firstParticle) then
                  write(*,*) "Expected lastParticle > firstParticle, got:"
                  write(*,*) "firstParticle=",firstParticle,"lastParticle=",lastParticle
                  flush(stdout)
                  stop 17
               endif
               
            else if (cmdarg_arg .eq. "--fname") then
               !Read the filename
               cmdarg_i = cmdarg_i+1
               call get_command_argument(cmdarg_i, cmdarg_arg,cmdarg_length,cmdarg_status)
               if (len_trim(cmdarg_arg)==0) then
                  write(*,*) "No filename found following --fname flag?"
                  flush(stdout)
                  stop 6
               endif

               if(cmdarg_status.eq.-1) then
                  write(*,*) "Filename was truncated to '"//cmdarg_arg//"'"
                  flush(stdout)
                  stop 7
               end if

               fname=cmdarg_arg
               
            else if (cmdarg_arg .eq. "--ofname") then
               !Read the filename
               cmdarg_i = cmdarg_i+1
               call get_command_argument(cmdarg_i, cmdarg_arg,cmdarg_length,cmdarg_status)
               if (len_trim(cmdarg_arg)==0) then
                  write(*,*) "No filename found following --ofname flag?"
                  flush(stdout)
                  stop 6
               endif

               if(cmdarg_status.eq.-1) then
                  write(*,*) "Filename was truncated to '"//cmdarg_arg//"'"
                  flush(stdout)
                  stop 7
               end if

               ofname=cmdarg_arg
               ofoutput = .TRUE.

            else
               write(*,*) "Did not recognize argument '"//trim(cmdarg_arg)//"'"
               write(*,*)
               write(*,*) "USAGE: read90 (--STF) (--SP <number> | --PR <number> <number>) (--fname <name>) (--ofname <name>)"
               write(*,*) "The --STF flag indicates that the file is in the new STF format"
               write(*,*) "The --SP and --PR flag can only be used together with the --STF flag,"
               write(*,*) " in order to extract only single particle (pair) data in a way that emulates the old format."
               write(*,*) " The number following should be a odd number, i.e. 1,3,5 etc. for particle pair 1+2/3+4/5+6 etc."
               write(*,*) " The --PR flag differs from the --SP flag in that it can be used to extract a range of particle pairs,"
               write(*,*) " first and last pair inclusive. Expects first number < last number"
               write(*,*) "The --fname flag can be used to specify the name of the file to read from (max 100 characters)."
               write(*,*) " If nothing is specified, it defaults to 'fort.190'."
               write(*,*) "The --ofname flag can be used to specify the name of the file to write to (max 100 characters)."
               write(*,*) " If nothing is specified, it defaults to stdout."
               flush(stdout)
               stop 1
            end if
         end if
         !write (*,*) cmdarg_i, cmdarg_arg
         cmdarg_i = cmdarg_i+1
      end do
      
!--open fort.190 (the input file)
      nfile=190
      n=0
      
      INQUIRE(FILE=fname,EXIST=hasInputFile)
      if (.not. hasInputFile) then
         write(*,'(a,a,a)') "Error in read90 - file '"//trim(fname)//"' was not found; try --help?"
         flush(stdout)
         stop 19
      endif
      
      open(nfile,file=fname,form='UNFORMATTED',status='OLD')

!--open the optional output file
      if (ofoutput .eqv. .TRUE.) then
         open(191,file=ofname,action='WRITE',status='REPLACE')
         ounit = 191
      end if

100   read(nfile,end=511,err=520) sixtit,commen,cdate,ctime,       &
     &progrm,ifipa,ilapa,itopa,icode,numl,qwc(1),qwc(2),qwc(3), clo(1), &
     &clop(1),clo(2),clop(2),clo(3),clop(3), di0(1),dip0(1),di0(2),dip0 &
     &(2),dummy,dummy, ta(1,1),ta(1,2),ta(1,3),ta(1,4),ta(1,5),ta(1,6), &
     &ta(2,1),ta(2,2),ta(2,3),ta(2,4),ta(2,5),ta(2,6), ta(3,1),ta(3,2), &
     &ta(3,3),ta(3,4),ta(3,5),ta(3,6), ta(4,1),ta(4,2),ta(4,3),ta(4,4), &
     &ta(4,5),ta(4,6), ta(5,1),ta(5,2),ta(5,3),ta(5,4),ta(5,5),ta(5,6), &
     &ta(6,1),ta(6,2),ta(6,3),ta(6,4),ta(6,5),ta(6,6), dmmac,dnms,dizu0,&
     &dnumlr,sigcor,dpscor

      if (STF .and. singleParticle.ne.-1) then
         if (singleParticle.ne.ifipa) then
            if (ifipa .eq. itopa-1) goto 210 !OK, we're done reading headers.
            goto 100
         endif
      endif
      if (STF .and. firstParticle.ne.-1) then
         if (.not. (ifipa.ge.firstParticle .and. ifipa.le.lastParticle)) then
            !We are outside of the interval of interest.
            if (ifipa .eq. itopa-1) goto 210 !OK, we're done  reading headers.
            goto 100
         endif
      endif
      
      write (ounit,'(1x,a,i10,a)') 'Read record = ', n, " (header)"
      ntwin=1
      if(ilapa.ne.ifipa) ntwin=2
      n=n+1 !Increase record number
      
      write (ounit,'(2x,a,i1,1x,i10,1x,i10)') 'ntwin, first_particle, last_particle : ',ntwin,ifipa,ilapa
      ! ifipa=0
      write(ounit,'(2x,a,a)')   "Title           : ", sixtit
      write(ounit,'(2x,a,a)')   "Comment         : ", commen
!      write(ounit,'(2x,a,a)')   "cdate           : ", cdate !Skip the date and time; won't match between different runs when diff
!      write(ounit,'(2x,a,a)')   "ctime           : ", ctime
      write(ounit,'(2x,a,a)')   "Program         : ", progrm
      write(ounit,'(2x,a,i10)') "Total particles : ", itopa
      write(ounit,'(2x,a,i10)') "icode           : ", icode
      write(ounit,'(2x,a,i10)') "numl            : ", numl
      errno=dtostr(qwc(1),ch1)
      write(ounit,'(2x,a,a)') "Qx:      ",ch1
      errno=dtostr(qwc(2),ch1)
      write(ounit,'(2x,a,a)') "Qy:      ",ch1
      errno=dtostr(qwc(3),ch1)
      write(ounit,'(2x,a,a)') "Qz:      ",ch1
      errno=dtostr(clo(1),ch1)
      write(ounit,'(2x,a,a)') "x_clo:   ", ch1
      errno=dtostr(clop(1),ch1)
      write(ounit,'(2x,a,a)') "x'_clo:  ",ch1
      errno=dtostr(clo(2),ch1)
      write(ounit,'(2x,a,a)') "y_clo:   ", ch1
      errno=dtostr(clop(2),ch1)
      write(ounit,'(2x,a,a)') "y'_clo:  ", ch1
      errno=dtostr(clo(3),ch1)
      write(ounit,'(2x,a,a)') "clo(3):  ", ch1
      errno=dtostr(clop(3),ch1)
      write(ounit,'(2x,a,a)') "clop(3): ", ch1
      errno=dtostr(di0(1),ch1)
      write(ounit,'(2x,a,a)') "di0(1):  ", ch1
      errno=dtostr(dip0(1),ch1)
      write(ounit,'(2x,a,a)') "dip0(1): ", ch1
      errno=dtostr(di0(2),ch1)
      write(ounit,'(2x,a,a)') "di0(2):  ", ch1
      errno=dtostr(dip0(2),ch1)
      write(ounit,'(2x,a,a)') "dip0(2): ",ch1
      do i=1,6
        do j=1,6
          errno=dtostr(ta(i,j),ch1)
          write (ounit,'(2x, A, I1, A, I1, A, A)') 'tas(',i,',',j,'): ',ch1
        enddo
      enddo
      errno=dtostr(dmmac,ch1)
      write(ounit,'(2x,a,a)') 'dmmac:   ', ch1
      errno=dtostr(dnms,ch1)
      write(ounit,'(2x,a,a)') 'dnms:    ', ch1
      errno=dtostr(dizu0,ch1)
      write(ounit,'(2x,a,a)') 'dizu0:   ', ch1
      errno=dtostr(dnumlr,ch1)
      write(ounit,'(2x,a,a)') 'dnumlr:  ', ch1
      errno=dtostr(sigcor,ch1)
      write(ounit,'(2x,a,a)') 'sigcor:  ', ch1
      errno=dtostr(dpscor,ch1)
      write(ounit,'(2x,a,a)') 'dpscor:  ', ch1

      if (STF) then
         if (ifipa .lt. itopa-1) goto 100 !Read more headers
      endif

      !Code for reading the first particle pair (turn)
210   if(ntwin.eq.1) read(nfile,end=530,err=510) ia,ifipa,b,c,d,e, &
     &f,g,h,p
      if(ntwin.eq.2) read(nfile,end=530,err=510) ia,ifipa,b,c,d,e, &
     &f,g,h,p, ilapa,b,c1,d1,e1,f1,g1,h1,p1
      if(ifipa.lt.1) goto 210
!     if(progrm.eq.'MAD') then
!     endif
      if (STF .and. singleParticle.ne.-1) then
         if (singleParticle.ne.ifipa) goto 210
      endif
      if (STF .and. firstParticle.ne.-1) then
         !Are we are outside of the interval of interest?
         if (.not. (ifipa.ge.firstParticle .and. ifipa.le.lastParticle)) goto 210
      endif

!--KEEP THE FIRST TURN NUMBER : IA0
      ia0=ia
      goto 212
      
  211 continue !Code for reading further pairs/turns
      if(ntwin.eq.1) read(nfile,end=540,err=510) ia,ifipa,b,c,d,e, &
     &f,g,h,p
      if(ntwin.eq.2) read(nfile,end=540,err=510) ia,ifipa,b,c,d,e, &
     &f,g,h,p, ilapa,b,c1,d1,e1,f1,g1,h1,p1
      if (STF .and. singleParticle.ne.-1) then
         if (singleParticle.ne.ifipa) goto 211
      endif
      if (STF .and. firstParticle.ne.-1) then
         !Are we are outside of the interval of interest?
         if (.not. (ifipa.ge.firstParticle .and. ifipa.le.lastParticle)) goto 211
      endif
      n=n+1
      
  212 continue
      write (ounit,'(1x,a,i10,a)') 'Read record = ', n, " (particle(s))"
      write (ounit,'(2x,a,i10,a,i10,a)') 'Turn ',ia,'   Particle ',ifipa, " (1st)"
      errno=dtostr(b,ch1)
      write(ounit,'(2x,a,a)') "dam:     ", ch1
      errno=dtostr(c,ch1)
      write(ounit,'(2x,a,a)') "x mm:    ", ch1
      errno=dtostr(d,ch1)
      write(ounit,'(2x,a,a)') "x' mrad: ", ch1
      errno=dtostr(e,ch1)
      write(ounit,'(2x,a,a)') "y mm:    ", ch1
      errno=dtostr(f,ch1)
      write(ounit,'(2x,a,a)') "y' mrad: ", ch1
      errno=dtostr(g,ch1)
      write(ounit,'(2x,a,a)') "ds mm    ", ch1
      errno=dtostr(h,ch1)
      write(ounit,'(2x,a,a)') "dp/p0:   ", ch1
      errno=dtostr(p,ch1)
      write(ounit,'(2x,a,a)') "Energy:  ", ch1
      if(ntwin.eq.2) then
         write (ounit,'(2x,a,i10,a,i10,a)') 'Turn ',ia,'   Particle ',ilapa, " (2nd)"
         errno=dtostr(b,ch1)
         write(ounit,'(2x,a,a)') "dam:     ", ch1
         errno=dtostr(c1,ch1)
         write(ounit,'(2x,a,a)') "x mm:    ", ch1
         errno=dtostr(d1,ch1)
         write(ounit,'(2x,a,a)') "x' mrad: ", ch1
         errno=dtostr(e1,ch1)
         write(ounit,'(2x,a,a)') "y mm:    ", ch1
         errno=dtostr(f1,ch1)
         write(ounit,'(2x,a,a)') "y' mrad: ", ch1
         errno=dtostr(g1,ch1)
         write(ounit,'(2x,a,a)') "ds mm    ", ch1
         errno=dtostr(h1,ch1)
         write(ounit,'(2x,a,a)') "dp/p0:   ", ch1
         errno=dtostr(p1,ch1)
         write(ounit,'(2x,a,a)') "Energy:  ", ch1
      endif
      goto 211
  510 continue
      write(ounit,*) nfile,'FILE CORRUPTED record',n
      goto 550
  511 continue
      write(ounit,*) nfile,'NO HEADER, EMPTY FILE!!!'
      goto 550
  520 continue
      write(ounit,*) nfile,'HEADER CORRUPTED'
      goto 550
  530 continue
      write(ounit,*) nfile,'NO DATA'
      goto 550
  540 continue
  550 continue
      if(ofoutput .eqv. .TRUE.) then
         close(ounit)
      endif
!----------------------------------------------------------------------
      end
      integer function dtostr(x,results)
! Uses the dtoa_c.c version of dtoa via the dtoaf.c interface copied from roundctl
      use, intrinsic :: iso_fortran_env, only : stdout=>output_unit
      implicit none
      double precision, intent(in) :: x
      character(*), intent(out) :: results
      integer dtoaf 
      integer ilen,mode,ndigits,decpoint,mysign
      integer i,l,d,e
      character(len=1) str(999)
      character(len=24) lstr
      character(len=3) e3
      mode=2
      ndigits=17
      ilen=dtoaf(x,mode,ndigits,decpoint,mysign,str(1))
      if (ilen.le.0.or.ilen.gt.17) then
! Always returns 17 or less characters as requested
      write (*,10000)
      write (*,*) 'Routine dtoa[f] returned string length ',ilen,'!!!'
      flush(stdout)
      stop 18
10000 format(5x///t10,'++++++++++++++++++++++++'/ t10,                  &
     &'+++++ERROR DETECTED+++++'/ t10,'++++++++++++++++++++++++'/ t10)
      endif
      lstr=' '
      do i=1,ilen
        lstr(i:i)=str(i)
      enddo
! Now try my formatting
      d=decpoint
      e=0
      l=1
      lstr=' '
      if (mysign.ne.0) then
        lstr(l:l)='-'
      endif
      if (decpoint.eq.9999) then
! Infinity or Nan
        do i=1,ilen
          lstr(l+i:l+i)=str(i)
        enddo
      else
! Pad with zeros
        do i=ilen+1,17
          str(i)='0'
        enddo
        if (decpoint.le.0) then
          e=decpoint-1
          d=1
        else
! I am using 17 as decision point to avoid dddd.e+eee
! but rather d.ddde+eee
          if (decpoint.ge.17) then
            e=decpoint-1
            d=1
          else
            d=decpoint
          endif
        endif
! and copy with the decimal point
        do i=1,17
          lstr(l+i:l+i)=str(i)
          if (i.eq.d) then
            l=l+1
            lstr(l+i:l+i)='.'
          endif
        enddo
! and add exponent e+/-nnn
        l=20
        lstr(l:l)='e'
        l=21
        lstr(l:l)='+'
        if (e.lt.0) then
          lstr(l:l)='-'
          e=-e
        endif
        l=22
        write (e3,'(I3.3)') e
        lstr(l:l+2)=e3(1:3)
      endif  
      results=lstr
      dtostr=24
      return
      end
