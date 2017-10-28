      program compf10
      use, intrinsic :: iso_fortran_env, only : output_unit
      implicit none
      double precision prob(60),prob1(60),eps(60)
      integer line,word,i
      logical diff,diffs
      logical hasInputFile
! Now compare the closed orbit in 53-58 as well
      do i=1,60
        prob(i) = 0
        prob1(i) = 0
        eps(i) = 0
      enddo
      line=0
      diff=.false.
      diffs=.false.
      
      hasInputFile = .false.
      INQUIRE(file="fort.20",EXIST=hasInputFile)
      if (.not. hasInputFile) then
         write(*,'(a,a)') "Error in compf10 - file 'fort.20'"//
     &        " was not found"
         flush(output_unit)
         stop 1
      endif
      hasInputFile = .false.
      INQUIRE(file="fort.21",EXIST=hasInputFile)
      if (.not. hasInputFile) then
         write(*,'(a,a)') "Error in compf10 - file 'fort.21'"//
     &        " was not found"
         flush(output_unit)
         stop 2
      endif
      
      open(20,status='OLD', file="fort.20")
      open(21,status='OLD', file="fort.21")
      
    1 read (20,*,end=100,err=98) prob
      line=line+1
      read (21,*,end=99,err=97) prob1
      if (diffs) diff=.true.
      diffs=.false.
      do word=1,51
        eps(word)=abs(prob(word))-abs(prob1(word))
        if (eps(word).ne.0d0) then
          if (.not.diffs) then
            write (*,*) "compf10_DIFF fort.10, line",line
            diffs=.true.
          endif
          write (*,*) 
     & "compf10_DIFF",word,prob(word),prob1(word),eps(word)
          if (abs(eps(word)).ge.1d-14) then
            write (*,*) 'HUGE!',abs(eps(word))
          endif
!       else
!         write (*,*) "compf10_SAME",word,prob(word)
        endif
      enddo 
      do word=53,58
        eps(word)=abs(prob(word))-abs(prob1(word))
        if (eps(word).ne.0d0) then
          if (.not.diffs) then
            write (*,*) "compf10_DIFF fort.10, line",line
            diffs=.true.
          endif
          write (*,*) 
     & "compf10_DIFF",word,prob(word),prob1(word),eps(word)
          if (abs(eps(word)).ge.1d-14) then
            write (*,*) 'HUGE!',abs(eps(word))
          endif
!       else
!         write (*,*) "compf10_SAME",word,prob(word)
        endif
      enddo 
      go to 1
 99   continue
      write (*,*) "Comparing VERSION ",prob(52)," to ",prob1(52)
      write (*,*)
     & "compf10_DIFF I/O error, wrong no of lines!! line no ",line
      flush(output_unit)
      stop 2
 98   continue
      write (*,*) "Comparing VERSION ",prob(52)," to ",prob1(52)
      write (*,*) "compf10_DIFF I/O error!! fort.20 line no ",line
      flush(output_unit)
      stop 3
 97   continue
      write (*,*) "Comparing VERSION ",prob(52)," to ",prob1(52)
      write (*,*) "compf10_DIFF I/O error!! fort.21 line no ",line
      flush(output_unit)
      stop 4
 100  continue
      if (line.eq.0) go to 99
      write (*,*) "Comparing VERSION ",prob(52)," to ",prob1(52)
      if (diff) then
        write (*,*) "compf10_DIFF after comparing ",line ,"lines"
        flush(output_unit)
        stop 1
      else
        write (*,*) "compf10_SAME after comparing ",line ,"lines"
        flush(output_unit)
        stop 0
      endif
      end
