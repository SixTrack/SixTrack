      program checkf1014
      implicit none
      double precision prob(60),prob1(60),lprob(60),lprob1(60)
      integer line,word,i
      logical diff,diffs
      character*80 buffer
! Now compare the closed orbit in 53-60 as well
      line=0
      diff=.false.
      diffs=.false.
    1 read (20,*,end=100,err=98) lprob
      do i=1,60
        write (buffer,'(E19.12)') lprob
        read (buffer,'(E19.12)') prob
      enddo
      line=line+1
      read (21,*,end=99,err=97) lprob1
      do i=1,60
        write (buffer,'(E19.12)') lprob1
        read (buffer,'(E19.12)') prob1
      enddo
      if (diffs) diff=.true.
      diffs=.false.
      do word=1,51
        if (prob(word).ne.prob1(word)) diffs=.true.
      enddo 
      do word=53,60
        if (prob(word).ne.prob1(word)) diffs=.true.
      enddo 
      if (diffs) then
        write (*,*)
        write (*,*) "checkf1014_DIFF fort.10, line",line
        do word=1,51
          if (prob(word).ne.prob1(word)) then
            write (*,*) "checkf1014_DIFF",word,prob(word),prob1(word)
          else
            write (*,*) "checkf1014_SAME",word,prob(word)
          endif
        enddo
        do word=53,60
          if (prob(word).ne.prob1(word)) then
            write (*,*) "checkf1014_DIFF",word,prob(word),prob1(word)
          else
            write (*,*) "checkf1014_SAME",word,prob(word)
          endif
        enddo
        write (*,*)
      else
        write (*,*) "checkf1014_SAME fort.10, line",line
      endif
      go to 1
 99   continue
      write (*,*) "Comparing VERSION ",prob(52)," to ",prob1(52)
      write (*,*)
     & "checkf1014_DIFF I/O error, wrong no of lines!! line no ",line
      stop
 98   continue
      write (*,*) "Comparing VERSION ",prob(52)," to ",prob1(52)
      write (*,*) "checkf1014_DIFF I/O error!! fort.20 line no ",line
      stop
 97   continue
      write (*,*) "Comparing VERSION ",prob(52)," to ",prob1(52)
      write (*,*) "checkf1014_DIFF I/O error!! fort.21 line no ",line
      stop
 100  continue
      if (line.eq.0) go to 99
      write (*,*) "Comparing VERSION ",prob(52)," to ",prob1(52)
      if (diff) then
        write (*,*) "checkf1014_DIFF after comparing ",line ,"lines"
      else
        write (*,*) "checkf1014_SAME after comparing ",line ,"lines"
      endif
      end
