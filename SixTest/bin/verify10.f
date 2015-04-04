      program verify10
      implicit none
      double precision prob(60),prob1(60),probv,probv1
      integer line,word
      logical diff,diffs
! Now checks the closed orbit as well as Total Turns
      line=0
      diff=.false.
    1 read (20,*,end=100,err=99) prob
      line=line+1
      read (21,*,end=99,err=99) prob1
      probv=prob(52)
      probv1=prob1(52)
      diffs=.false.
      do word=1,51
        if (prob(word).ne.prob1(word)) diffs=.true.
      enddo 
      do word=53,58
        if (prob(word).ne.prob1(word)) diffs=.true.
      enddo 
      if (diffs) then
        diff=.true.
        write (*,*)
        write (*,*) "verify10_DIFF fort.10, line",line
        do word=1,51
          if (prob(word).ne.prob1(word)) then
            write (*,*) "verify10_DIFF",word,prob(word),prob1(word)
          else
            write (*,*) "verify10_SAME",word,prob(word)
          endif
        enddo
        do word=53,58
          if (prob(word).ne.prob1(word)) then
            write (*,*) "verify10_DIFF",word,prob(word),prob1(word)
          else
            write (*,*) "verify10_SAME",word,prob(word)
          endif
        enddo
      else
        write (*,*) "verify10_SAME fort.10, line",line
      endif
      go to 1
 99   continue
      write (*,*) "Comparing VERSION ",probv," to ",probv1
      write (*,*)
     & "verify10_DIFF I/O error, wrong no of lines!!?? line no ",line
      call exit(1)
      stop
 100  continue
      if (line.eq.0) go to 99
      if (diff) then
        write (*,*) "Comparing VERSION ",probv," to ",probv1
        write (*,*) "Different after comparing ",line," lines"
        call exit(2)
      else
        write (*,*) "Comparing VERSION ",probv," to ",probv1
        write (*,*) "verify10_SAME after comparing ",line," lines"
        call exit(0) 
      endif
      end
