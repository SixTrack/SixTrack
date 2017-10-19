      program astuce
 
************************************************************************
*                                                                      *
*   Program to extract a FORTRAN version from a source form in which   *
*   common decks, source decks, and flags are defined similar to CDC   *
*   UPDATE or Historian. The program is intended for usage with MAD    *
*   only.                                                              *
*                                                                      *
*   Usage.                                                             *
*   The SOURCE input file is presented on logical unit 11, the FORTRAN *
*   output file will be written onto unit 12. Command input unit 5,    *
*   questions and statistics unit 6.                                   *
*                                                                      *
*   Commands.                                                          *
*   ASTUCE will use the first character of the first input line in the *
*   SOURCE file if it is either '*' or '+', otherwise it will use      *
*   the standard master character  '*'. To set it EXPLICITLY to        *
*   another character, e.g. '+' as in MAD, use                         *
*   *=+                                                                *
*   Flags are defined in the following way:                            *
*   DF FLAG1,FLAG2,FLAG3                                               *
*   DF FLAG4,...             as many lines as necessary.               *
*   Deck ranges:                                                       *
*   E DECK1,DECK2.DECK3,DECK4,... etc.                                 *
*   meaning                                                            *
*   "extract DECK1, all decks between DECK2 and DECK3 inclusive,       *
*   DECK4,... ".                                                       *
*                                                                      *
*                                                                      *
*   Last command: EX                                                   *
*   will start the extraction. If no deck range is specified, all      *
*   decks will be extracted.                                           *
*                                                                      *
*   Author: H. Grote, CERN/DD     Feb. 25, 1988                        *
*   Copyright CERN 2014                                                *
*   This software is distributed under the terms of the GNU            *
*   Lesser General Public License version 2.1, copied verbatim in the  *
*   file ``COPYING''.                                                  *
*
*  In applying this licence, CERN does not waive the privileges and    *
*  immunities granted to it by virtue of its status as an              *
*  Intergovernmental Organization or submit itself to any jurisdiction.*
*                                                                      *
*                                                                      *
************************************************************************

      call astdef
      call astini
      call inisum
      call pass1
      call pass2
      call summup
      call termin
      end
      subroutine astdef
*---------*---------*---------*---------*---------*---------*---------*-
      parameter (mline = 80, mtab = 8, irunit = 5, ipunit = 6,
     +inunit = 11, iounit = 12, maxflg = 100, mudeck = 1000,
     +mxdeck = 10000, mxdcom = 1000, mxlcom = 10000, marg = 4)
      common / aux / nlin, nlout, ndcomm, nlcomm, nudeck, nuflag,
     +ndflag, nddeck, iffull, idfull, iasepl(marg),
     +isc(mxdcom), ilc(mxdcom), idchec(mudeck,2), ifchec(maxflg)
 
      common / saux / suflag(maxflg+1), sdflag(maxflg+1),
     +sargin(2*marg), inpnam, outnam, smast
      character suflag * (mtab), sdflag * (mtab), sargin * (mline),
     +inpnam * 60, outnam * 60, smast * 1
 
      common / sbuff / sbcomm(mxlcom), sncomm(mxdcom), sddeck(mxdeck),
     +sudeck(mudeck,2)
      character sbcomm * (mline), sncomm * (mtab), sddeck * (mtab),
     +sudeck * (mtab)
*---------*---------*---------*---------*---------*---------*---------*-
 
      smast='*'
      nuflag = 0
      ndflag = 0
      nudeck = 0
      nddeck = 0
      ndcomm = 0
      nlcomm = 0
      nrdeck = 0
      nwdeck = 0
      iffull = 0
      idfull = 0
      do 10 i = 1, 2
         do 10 j = 1, mudeck
   10 idchec(j,i) = 0
      do 20 i = 1, maxflg
   20 ifchec(i) = 0
      do 30  i = 1, marg
   30 iasepl(i) = 0
      end
      subroutine astini
*---------*---------*---------*---------*---------*---------*---------*-
      parameter (mline = 80, mtab = 8, irunit = 5, ipunit = 6,
     +inunit = 11, iounit = 12, maxflg = 100, mudeck = 1000,
     +mxdeck = 10000, mxdcom = 1000, mxlcom = 10000, marg = 4)
      common / aux / nlin, nlout, ndcomm, nlcomm, nudeck, nuflag,
     +ndflag, nddeck, iffull, idfull, iasepl(marg),
     +isc(mxdcom), ilc(mxdcom), idchec(mudeck,2), ifchec(maxflg)
 
      common / saux / suflag(maxflg+1), sdflag(maxflg+1),
     +sargin(2*marg), inpnam, outnam, smast
      character suflag * (mtab), sdflag * (mtab), sargin * (mline),
     +inpnam * 60, outnam * 60, smast * 1
 
      common / sbuff / sbcomm(mxlcom), sncomm(mxdcom), sddeck(mxdeck),
     +sudeck(mudeck,2)
      character sbcomm * (mline), sncomm * (mtab), sddeck * (mtab),
     +sudeck * (mtab)
*---------*---------*---------*---------*---------*---------*---------*-
      character sline*(mline)
      call box('ASTUCE - start of execution')
*--- get arguments
      call gettem
      if (iasepl(1) .eq. 0)  then
        write(ipunit,*) ' please enter input (SOURCE) file name'
        read(irunit,'(A)',end=999) inpnam
      else
        inpnam = sargin(iasepl(1))
      endif
      write(ipunit,*)  inpnam
      open(inunit,file=inpnam,status='OLD')
      if (iasepl(2) .eq. 0)  then
        write(ipunit,*) ' please enter output (ftn) file name'
        read(irunit,'(A)',end=999) outnam
      else
        outnam = sargin(iasepl(2))
      endif
      write(ipunit,*)  outnam
      open(iounit,file=outnam,status='UNKNOWN')
*--- get master character from first line
      read(inunit,'(A)',end=1) sline
      if (sline(1:1) .eq. '+')  smast = '+'
      rewind inunit
      goto 2
    1 continue
      print *, ' '
      print *, ' +++ Warning: input file (SOURCE) is empty, STOP'
      print *, ' '
      stop 1
    2 continue
      write(ipunit,10000)  smast
   10 continue
      if (iasepl(3) .eq. 0)  then
        read(irunit,'(A)',end=999) sline
        call upper(sline)
        if(sline(:2).eq.'EX') goto 999
        l=lastnb(sline)
        if(l.ne.0)  then
          k=lfirnb(sline(:l))
          kast=index(sline,'*')
          if(kast.gt.0) then
            k1=lfirnb(sline(kast+1:l))
            if(k1.ne.0) then
              k1=k1+kast
              if(sline(k1:k1).eq.'=') then
                k2=lfirnb(sline(k1+1:l))
                if(k2.ne.0) then
                  k2=k2+k1
                  smast=sline(k2:k2)
                endif
              endif
            endif
          elseif(sline(k:k+1).eq.'DF')  then
            call flags1(sline,k,l)
          elseif(sline(k:k).eq.'E')  then
            call decks(sline,k,l)
          endif
        endif
        goto 10
      else
        sline = 'df ' // sargin(iasepl(3))
        call upper(sline)
        l = lastnb(sline)
        call flags1(sline, 1, l)
        if (iasepl(4) .ne. 0)  then
          sline = 'e ' // sargin(iasepl(4))
          call upper(sline)
          l = lastnb(sline)
          call decks(sline, 1, l)
        endif
      endif
10000 format(//' change the master character (current =',a1,
     +') with a line'//
     +' *=...'//' where "..." stands for the new one.'//
     +' Set flags and deck names in HISTE fashion, i.e.'//
     +' DF flag1,flag2,etc. and'/ ' E deck1,deck2.deck3,deck4,...'//
     +' Defaults: no flag given: no flag set; no deck given: ALL decks'/
     +' EX or E.O.F. or (under VM) empty <CR> will terminate'//)
  999 end
      subroutine box(string)
*---------*---------*---------*---------*---------*---------*---------*-
      parameter (mline = 80, mtab = 8, irunit = 5, ipunit = 6,
     +inunit = 11, iounit = 12, maxflg = 100, mudeck = 1000,
     +mxdeck = 10000, mxdcom = 1000, mxlcom = 10000, marg = 4)
      common / aux / nlin, nlout, ndcomm, nlcomm, nudeck, nuflag,
     +ndflag, nddeck, iffull, idfull, iasepl(marg),
     +isc(mxdcom), ilc(mxdcom), idchec(mudeck,2), ifchec(maxflg)
 
      common / saux / suflag(maxflg+1), sdflag(maxflg+1),
     +sargin(2*marg), inpnam, outnam, smast
      character suflag * (mtab), sdflag * (mtab), sargin * (mline),
     +inpnam * 60, outnam * 60, smast * 1
 
      common / sbuff / sbcomm(mxlcom), sncomm(mxdcom), sddeck(mxdeck),
     +sudeck(mudeck,2)
      character sbcomm * (mline), sncomm * (mtab), sddeck * (mtab),
     +sudeck * (mtab)
*---------*---------*---------*---------*---------*---------*---------*-
 
      character *(*)  string
      write(ipunit,10000) string
10000 format(/t5,15('++++')/t5,'+',t64,'+'/t5,'+',t15,a,t64,'+'/ t5,'+',
     +t64,'+'/t5,15('++++')/)
      end
      subroutine decks(sline,k1,l)
*---------*---------*---------*---------*---------*---------*---------*-
      parameter (mline = 80, mtab = 8, irunit = 5, ipunit = 6,
     +inunit = 11, iounit = 12, maxflg = 100, mudeck = 1000,
     +mxdeck = 10000, mxdcom = 1000, mxlcom = 10000, marg = 4)
      common / aux / nlin, nlout, ndcomm, nlcomm, nudeck, nuflag,
     +ndflag, nddeck, iffull, idfull, iasepl(marg),
     +isc(mxdcom), ilc(mxdcom), idchec(mudeck,2), ifchec(maxflg)
 
      common / saux / suflag(maxflg+1), sdflag(maxflg+1),
     +sargin(2*marg), inpnam, outnam, smast
      character suflag * (mtab), sdflag * (mtab), sargin * (mline),
     +inpnam * 60, outnam * 60, smast * 1
 
      common / sbuff / sbcomm(mxlcom), sncomm(mxdcom), sddeck(mxdeck),
     +sudeck(mudeck,2)
      character sbcomm * (mline), sncomm * (mtab), sddeck * (mtab),
     +sudeck * (mtab)
*---------*---------*---------*---------*---------*---------*---------*-
      character sline*(mline),sname*(mtab)
 
      if(idfull.eq.0)  then
         call getapo(sline,k1,l,k)
         if(k.gt.0) then
   10       continue
            call getnam(sline,k,l,kfch,klch)
            if(kfch.ne.0) then
               sname=sline(kfch:klch)
               call namtab(sname,sudeck,nudeck,ipos)
               if(ipos.gt.0) then
                  if(nudeck.eq.mudeck) then
                     call errms2 (
     +               '++++++ WARNING - user deck table full at:',
     +               nudeck)
                     idfull=1
                     stop 1
                     goto 999
                  else
                     do 20 i=nudeck,ipos,-1
   20                sudeck(i+1,2)=sudeck(i,2)
                     nudeck=nudeck+1
                     if(klch.lt.l) then
                        k=klch+1
                        if(sline(k:k).eq.'.') then
                           call getnam(sline,k+1,l,kfch,klch)
                           if(kfch.ne.0) then
                              sname=sline(kfch:klch)
                           endif
                        endif
                     endif
                     sudeck(ipos,2)=sname
                  endif
               endif
               call getapo(sline,klch,l,k)
               if(k.gt.0) goto 10
            endif
         endif
      endif
  999 end
      subroutine errms1(text,sline)
*---------*---------*---------*---------*---------*---------*---------*-
      parameter (mline = 80, mtab = 8, irunit = 5, ipunit = 6,
     +inunit = 11, iounit = 12, maxflg = 100, mudeck = 1000,
     +mxdeck = 10000, mxdcom = 1000, mxlcom = 10000, marg = 4)
      common / aux / nlin, nlout, ndcomm, nlcomm, nudeck, nuflag,
     +ndflag, nddeck, iffull, idfull, iasepl(marg),
     +isc(mxdcom), ilc(mxdcom), idchec(mudeck,2), ifchec(maxflg)
 
      common / saux / suflag(maxflg+1), sdflag(maxflg+1),
     +sargin(2*marg), inpnam, outnam, smast
      character suflag * (mtab), sdflag * (mtab), sargin * (mline),
     +inpnam * 60, outnam * 60, smast * 1
 
      common / sbuff / sbcomm(mxlcom), sncomm(mxdcom), sddeck(mxdeck),
     +sudeck(mudeck,2)
      character sbcomm * (mline), sncomm * (mtab), sddeck * (mtab),
     +sudeck * (mtab)
*---------*---------*---------*---------*---------*---------*---------*-
 
      character *(*)  text,sline
      write(ipunit,10000) text,nlin,sline
10000 format(/1x,a,' at line:',i6/1x,a)
      end
      subroutine errms2(string,n)
*---------*---------*---------*---------*---------*---------*---------*-
      parameter (mline = 80, mtab = 8, irunit = 5, ipunit = 6,
     +inunit = 11, iounit = 12, maxflg = 100, mudeck = 1000,
     +mxdeck = 10000, mxdcom = 1000, mxlcom = 10000, marg = 4)
      common / aux / nlin, nlout, ndcomm, nlcomm, nudeck, nuflag,
     +ndflag, nddeck, iffull, idfull, iasepl(marg),
     +isc(mxdcom), ilc(mxdcom), idchec(mudeck,2), ifchec(maxflg)
 
      common / saux / suflag(maxflg+1), sdflag(maxflg+1),
     +sargin(2*marg), inpnam, outnam, smast
      character suflag * (mtab), sdflag * (mtab), sargin * (mline),
     +inpnam * 60, outnam * 60, smast * 1
 
      common / sbuff / sbcomm(mxlcom), sncomm(mxdcom), sddeck(mxdeck),
     +sudeck(mudeck,2)
      character sbcomm * (mline), sncomm * (mtab), sddeck * (mtab),
     +sudeck * (mtab)
*---------*---------*---------*---------*---------*---------*---------*-
      character *(*)  string
 
      write(ipunit,'(1X,A,I6)')  string,n
      end
      subroutine expand(string)
*---------*---------*---------*---------*---------*---------*---------*-
      parameter (mline = 80, mtab = 8, irunit = 5, ipunit = 6,
     +inunit = 11, iounit = 12, maxflg = 100, mudeck = 1000,
     +mxdeck = 10000, mxdcom = 1000, mxlcom = 10000, marg = 4)
      common / aux / nlin, nlout, ndcomm, nlcomm, nudeck, nuflag,
     +ndflag, nddeck, iffull, idfull, iasepl(marg),
     +isc(mxdcom), ilc(mxdcom), idchec(mudeck,2), ifchec(maxflg)
 
      common / saux / suflag(maxflg+1), sdflag(maxflg+1),
     +sargin(2*marg), inpnam, outnam, smast
      character suflag * (mtab), sdflag * (mtab), sargin * (mline),
     +inpnam * 60, outnam * 60, smast * 1
 
      common / sbuff / sbcomm(mxlcom), sncomm(mxdcom), sddeck(mxdeck),
     +sudeck(mudeck,2)
      character sbcomm * (mline), sncomm * (mtab), sddeck * (mtab),
     +sudeck * (mtab)
*---------*---------*---------*---------*---------*---------*---------*-
      character string*(*),sline*(mline)
      character sexp(mxdcom)*(mtab)
      integer iki(mxdcom),ikj(mxdcom),iknl(mxdcom),iknn(mxdcom)
      logical ifdecl,eldecl,eidecl,cddecl,dkdecl,cadecl
      ifdecl(sline)=sline(2:3).eq.'IF'.and.(sline(4:4).eq.' '.or.
     +sline(4:4).eq.',')
      eldecl(sline)=sline(2:3).eq.'EL'.and.(sline(4:4).eq.' '.or.
     +sline(4:4).eq.',').or.sline(2:5).eq.'ELSE'.and.(sline(6:6)
     +.eq.' '.or.sline(6:6).eq.',')
      eidecl(sline)=sline(2:3).eq.'EI'.and.(sline(4:4).eq.' '.or.
     +sline(4:4).eq.',').or.sline(2:6).eq.'ENDIF'.and.(sline(7:7)
     +.eq.' '.or.sline(7:7).eq.',')
      dkdecl(sline)=sline(2:3).eq.'DK'.and.(sline(4:4).eq.' '.or.
     +sline(4:4).eq.',').or.sline(2:5).eq.'DECK'.and.(sline(6:6)
     +.eq.' '.or.sline(6:6).eq.',')
      cddecl(sline)=sline(2:3).eq.'CD'.and.(sline(4:4).eq.' '.or.
     +sline(4:4).eq.',').or.sline(2:8).eq.'COMDECK'.and.(sline(9:9)
     +.eq.' '.or.sline(9:9).eq.',')
      cadecl(sline)=sline(2:3).eq.'CA'.and.(sline(4:4).eq.' '.or.
     +sline(4:4).eq.',').or.sline(2:5).eq.'CALL'.and.(sline(6:6)
     +.eq.' '.or.sline(6:6).eq.',')
      sline=string
      lev=0
      ipt=0
      npt=0
   10 continue
      if(ipt+34.gt.mxdcom)  then
         call errms1(
     +   ' ++++++ WARNING - limit of COMDECK expansion (loop ?)',
     +   string)
         stop 1
         goto 999
      endif
      call getall(sline,sexp,ipt,npt)
      i=ipt
   20 i=i+1
      if(i.gt.npt)  then
         if(lev.eq.0) goto 999
         i=iki(lev)
         npt=iknn(lev)
         j=ikj(lev)
         nl=iknl(lev)
         lev=lev-1
         goto 30
      endif
      call namsrc(sexp(i),sncomm,ndcomm,ipos,last)
      if(ipos.eq.0)  then
         call errms1( ' ++++++ WARNING - comdeck '//sexp(i)//
     +   ' not found',string)
         stop 1
         goto 20
      endif
      j=isc(ipos)-1
      nl=ilc(ipos)
      if (nl .eq. 0)  goto 20
   30 j=j+1
      if(j.gt.nl) goto 20
      sline=sbcomm(j)
      if(cadecl(sline))  then
         if(lev.eq.mxdcom) then
            call errms1(
     +      ' ++++++ WARNING - limit of COMDECK expansion (loop ?)',
     +      sline)
            stop 1
            goto 999
         endif
         lev=lev+1
         iki(lev)=i
         ikj(lev)=j
         iknn(lev)=npt
         iknl(lev)=nl
         ipt=npt
         goto 10
      endif
      lp=max(1,lastnb(sline))
      write(iounit,'(A)')  sline(:lp)
      nlout=nlout+1
      goto 30
  999 end
      subroutine flags1(sline,k1,l)
*---------*---------*---------*---------*---------*---------*---------*-
      parameter (mline = 80, mtab = 8, irunit = 5, ipunit = 6,
     +inunit = 11, iounit = 12, maxflg = 100, mudeck = 1000,
     +mxdeck = 10000, mxdcom = 1000, mxlcom = 10000, marg = 4)
      common / aux / nlin, nlout, ndcomm, nlcomm, nudeck, nuflag,
     +ndflag, nddeck, iffull, idfull, iasepl(marg),
     +isc(mxdcom), ilc(mxdcom), idchec(mudeck,2), ifchec(maxflg)
 
      common / saux / suflag(maxflg+1), sdflag(maxflg+1),
     +sargin(2*marg), inpnam, outnam, smast
      character suflag * (mtab), sdflag * (mtab), sargin * (mline),
     +inpnam * 60, outnam * 60, smast * 1
 
      common / sbuff / sbcomm(mxlcom), sncomm(mxdcom), sddeck(mxdeck),
     +sudeck(mudeck,2)
      character sbcomm * (mline), sncomm * (mtab), sddeck * (mtab),
     +sudeck * (mtab)
*---------*---------*---------*---------*---------*---------*---------*-
      character sline*(mline),sflag*(mtab)
 
      if(iffull.eq.0)  then
         call getapo(sline,k1,l,k)
         if(k.gt.0) then
   10       continue
            call getnam(sline,k,l,kfch,klch)
            if(kfch.ne.0) then
               sflag=sline(kfch:klch)
               call namtab(sflag,suflag,nuflag,ipos)
               if(ipos.gt.0) then
                  if(nuflag.eq.maxflg) then
                     call errms2 (
     +               ' ++++++ WARNING - user flag table full at:',
     +               nuflag)
                     iffull=1
                     stop 1
                     goto 999
                  else
                     nuflag=nuflag+1
                  endif
               endif
               call getapo(sline,klch,l,k)
               if(k.ne.0) goto 10
            endif
         endif
      endif
  999 end
      subroutine flags2(sline,l)
*---------*---------*---------*---------*---------*---------*---------*-
      parameter (mline = 80, mtab = 8, irunit = 5, ipunit = 6,
     +inunit = 11, iounit = 12, maxflg = 100, mudeck = 1000,
     +mxdeck = 10000, mxdcom = 1000, mxlcom = 10000, marg = 4)
      common / aux / nlin, nlout, ndcomm, nlcomm, nudeck, nuflag,
     +ndflag, nddeck, iffull, idfull, iasepl(marg),
     +isc(mxdcom), ilc(mxdcom), idchec(mudeck,2), ifchec(maxflg)
 
      common / saux / suflag(maxflg+1), sdflag(maxflg+1),
     +sargin(2*marg), inpnam, outnam, smast
      character suflag * (mtab), sdflag * (mtab), sargin * (mline),
     +inpnam * 60, outnam * 60, smast * 1
 
      common / sbuff / sbcomm(mxlcom), sncomm(mxdcom), sddeck(mxdeck),
     +sudeck(mudeck,2)
      character sbcomm * (mline), sncomm * (mtab), sddeck * (mtab),
     +sudeck * (mtab)
*---------*---------*---------*---------*---------*---------*---------*-
      character sline*(mline),sflag*(mtab)
 
      data ifull/0/
 
      if(ifull.ne.0) goto 999
      call getapo(sline,2,l,k1)
      if(k1.eq.0) goto 999
   10 continue
      call mygetlog(sline,k1,l,kfch,klch,it)
      if(kfch.gt.0)  then
         sflag=sline(kfch:klch)
         call namtab(sflag,sdflag,ndflag,ipos)
         if(ipos.gt.0)  then
            if(ndflag.eq.maxflg) then
               call errms2(
     +         ' ++++++ WARNING - flag table full at limit =' ,ndflag)
               stop 1
               goto 999
            endif
            ndflag=ndflag+1
         endif
         k1=klch+1
         call getopt(sline,k1,l,kfch,klch)
         if(kfch.eq.0) goto 999
         k1=klch+1
         if(k1.lt.l) goto 10
      endif
  999 end
      subroutine getall(sline,sexp,ipt,npt)
*---------*---------*---------*---------*---------*---------*---------*-
      parameter (mline = 80, mtab = 8, irunit = 5, ipunit = 6,
     +inunit = 11, iounit = 12, maxflg = 100, mudeck = 1000,
     +mxdeck = 10000, mxdcom = 1000, mxlcom = 10000, marg = 4)
      common / aux / nlin, nlout, ndcomm, nlcomm, nudeck, nuflag,
     +ndflag, nddeck, iffull, idfull, iasepl(marg),
     +isc(mxdcom), ilc(mxdcom), idchec(mudeck,2), ifchec(maxflg)
 
      common / saux / suflag(maxflg+1), sdflag(maxflg+1),
     +sargin(2*marg), inpnam, outnam, smast
      character suflag * (mtab), sdflag * (mtab), sargin * (mline),
     +inpnam * 60, outnam * 60, smast * 1
 
      common / sbuff / sbcomm(mxlcom), sncomm(mxdcom), sddeck(mxdeck),
     +sudeck(mudeck,2)
      character sbcomm * (mline), sncomm * (mtab), sddeck * (mtab),
     +sudeck * (mtab)
*---------*---------*---------*---------*---------*---------*---------*-
      character sexp(*)*(mtab),sline*(mline)
 
      npt=ipt
      l=lastnb(sline)
      call getapo(sline,2,l,k1)
      if(k1.gt.0)  then
   10    continue
         call getnam(sline,k1,l,kfch,klch)
         if(kfch.gt.0) then
            npt=npt+1
            sexp(npt)=sline(kfch:klch)
            k1=klch+1
            if(k1.lt.l) then
               if(sline(k1:k1).eq.',') then
                  k1=k1+1
                  goto 10
               endif
            endif
         endif
      endif
      end
      subroutine getapo(string,k1,k2,kfch)
 
      character string*(*)
      character alphan*36,alphae*37
 
      data alphan/'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'/
      data alphae/'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.'/
 
      kfch=0
      do 10 i1=k1,k2
         if(index(alphan,string(i1:i1)).eq.0) goto 20
   10 continue
      goto 999
   20 continue
      do 30 i=i1,k2
         if(index(alphae,string(i:i)).ne.0)  then
            kfch=i
            goto 999
         endif
   30 continue
  999 end
      subroutine gettem
*---------*---------*---------*---------*---------*---------*---------*-
      parameter (mline = 80, mtab = 8, irunit = 5, ipunit = 6,
     +inunit = 11, iounit = 12, maxflg = 100, mudeck = 1000,
     +mxdeck = 10000, mxdcom = 1000, mxlcom = 10000, marg = 4)
      common / aux / nlin, nlout, ndcomm, nlcomm, nudeck, nuflag,
     +ndflag, nddeck, iffull, idfull, iasepl(marg),
     +isc(mxdcom), ilc(mxdcom), idchec(mudeck,2), ifchec(maxflg)
 
      common / saux / suflag(maxflg+1), sdflag(maxflg+1),
     +sargin(2*marg), inpnam, outnam, smast
      character suflag * (mtab), sdflag * (mtab), sargin * (mline),
     +inpnam * 60, outnam * 60, smast * 1
 
      common / sbuff / sbcomm(mxlcom), sncomm(mxdcom), sddeck(mxdeck),
     +sudeck(mudeck,2)
      character sbcomm * (mline), sncomm * (mtab), sddeck * (mtab),
     +sudeck * (mtab)
*---------*---------*---------*---------*---------*---------*---------*-
 
      m = min(2 * marg, command_argument_count())
      do 10  i = 1, m
        call get_command_argument(i, sargin(i))
   10 continue
      do 20  i = 1, m - 1
        if(sargin(i)(:2) .eq. '-s')  then
          iasepl(1) = i + 1
        elseif(sargin(i)(:2) .eq. '-f')  then
          iasepl(2) = i + 1
        elseif(sargin(i)(:2) .eq. '-d')  then
          iasepl(3) = i + 1
        elseif(sargin(i)(:2) .eq. '-e')  then
          iasepl(4) = i + 1
        endif
   20 continue
      end
      subroutine mygetlog(string,k1,k2,kfch,klch,it)
 
      character string*(*)
      character alphan*36,alphae*37
 
      data alphan/'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'/
      data alphae/'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789.'/
 
*--- find first + last alphanum.
*    and store in IT the values of a preceding '.NOT.'
      kfch=0
      do 10 i2=k1,k2
         if(index(alphae,string(i2:i2)).ne.0)  then
            if(i2.ne.k1) goto 999
            if(string(i2:i2).eq.'.')  then
               call getopt(string,i2,k2,kopt1,kopt2)
               if(kopt1.eq.0) goto 999
               if(string(kopt1:kopt2).ne.'.NOT.') goto 999
               it=0
               i3=i2+5
            else
               it=1
               i3=i2
            endif
            goto 20
         endif
   10 continue
      goto 999
   20 continue
      do 30 i2=i3,k2
         if(index(alphan,string(i2:i2)).ne.0)  then
            if(kfch.eq.0)  then
               if(i2.ne.i3) goto 999
               kfch=i2
            endif
            klch=i2
         elseif(kfch.ne.0)  then
            goto 999
         endif
   30 continue
  999 end
      subroutine getnam(string,k1,k2,kfch,klch)
 
      character string*(*)
      character alphan*37
 
      data alphan/'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_'/
 
*--- find first + last alphanum.
      kfch=0
      do 10 i2=k1,k2
         if(index(alphan,string(i2:i2)).ne.0)  then
            if(kfch.eq.0)  then
               if(i2.ne.k1) goto 999
               kfch=i2
            endif
            klch=i2
         elseif(kfch.ne.0)  then
            goto 999
         endif
   10 continue
  999 end
      subroutine getopt(string,k1,k2,kfch,klch)
 
      character string*(*)
 
      kfch=0
      if(k2.le.k1) goto 999
      if(string(k1:k1).ne.'.') goto 999
      do 10 i=k1+1,k2
         if(string(i:i).eq.'.')  then
            kfch=k1
            klch=i
            goto 999
         endif
   10 continue
  999 end
      logical function ifval(string)
*---------*---------*---------*---------*---------*---------*---------*-
      parameter (mline = 80, mtab = 8, irunit = 5, ipunit = 6,
     +inunit = 11, iounit = 12, maxflg = 100, mudeck = 1000,
     +mxdeck = 10000, mxdcom = 1000, mxlcom = 10000, marg = 4)
      common / aux / nlin, nlout, ndcomm, nlcomm, nudeck, nuflag,
     +ndflag, nddeck, iffull, idfull, iasepl(marg),
     +isc(mxdcom), ilc(mxdcom), idchec(mudeck,2), ifchec(maxflg)
 
      common / saux / suflag(maxflg+1), sdflag(maxflg+1),
     +sargin(2*marg), inpnam, outnam, smast
      character suflag * (mtab), sdflag * (mtab), sargin * (mline),
     +inpnam * 60, outnam * 60, smast * 1
 
      common / sbuff / sbcomm(mxlcom), sncomm(mxdcom), sddeck(mxdeck),
     +sudeck(mudeck,2)
      character sbcomm * (mline), sncomm * (mtab), sddeck * (mtab),
     +sudeck * (mtab)
*---------*---------*---------*---------*---------*---------*---------*-
      character string*(*),sname*(mtab)
      logical laux
 
      ifval=.false.
      ior=1
      l=len(string)
      call getapo(string,2,l,k)
      if(k.eq.0) goto 999
   10 continue
      call mygetlog(string,k,l,kfch,klch,it)
      if(kfch.eq.0) goto 999
      sname=string(kfch:klch)
      call namsrc(sname,suflag,nuflag,ipos,last)
      laux=ipos.gt.0.and.it.eq.1.or.ipos.eq.0.and.it.eq.0
      if(ior.eq.0)  then
         ifval=ifval.and.laux
      else
         ifval=ifval.or.laux
      endif
      call getopt(string,klch+1,l,kfch,klch)
      if(kfch.gt.0)  then
         if(string(kfch:klch).eq.'.OR.')  then
            ior=1
         else
            ior=0
         endif
         k=klch+1
         goto 10
      endif
  999 end
      subroutine inisum
*---------*---------*---------*---------*---------*---------*---------*-
      parameter (mline = 80, mtab = 8, irunit = 5, ipunit = 6,
     +inunit = 11, iounit = 12, maxflg = 100, mudeck = 1000,
     +mxdeck = 10000, mxdcom = 1000, mxlcom = 10000, marg = 4)
      common / aux / nlin, nlout, ndcomm, nlcomm, nudeck, nuflag,
     +ndflag, nddeck, iffull, idfull, iasepl(marg),
     +isc(mxdcom), ilc(mxdcom), idchec(mudeck,2), ifchec(maxflg)
 
      common / saux / suflag(maxflg+1), sdflag(maxflg+1),
     +sargin(2*marg), inpnam, outnam, smast
      character suflag * (mtab), sdflag * (mtab), sargin * (mline),
     +inpnam * 60, outnam * 60, smast * 1
 
      common / sbuff / sbcomm(mxlcom), sncomm(mxdcom), sddeck(mxdeck),
     +sudeck(mudeck,2)
      character sbcomm * (mline), sncomm * (mtab), sddeck * (mtab),
     +sudeck * (mtab)
*---------*---------*---------*---------*---------*---------*---------*-
 
      call box('User request summary')
      if(smast.ne.'*') write(ipunit,10000) smast
      write(ipunit,10010) nuflag
      write(ipunit,10020) (suflag(i),i=1,nuflag)
      write(ipunit,10030) nudeck
      do 10 i=1,nudeck
         if(sudeck(i,1).eq.sudeck(i,2)) then
            write(ipunit,'(1X,A)') sudeck(i,1)
         else
            write(ipunit,'(1X,A,'' to '',A)') sudeck(i,1),sudeck(i,2)
         endif
   10 continue
10000 format(/' master character changed to: ',a1/)
10010 format(' no. of flags:',i5/)
10020 format(1x,5a10)
10030 format(/' no. of decks:',i5/)
      end
      logical function inrang(snin)
*---------*---------*---------*---------*---------*---------*---------*-
      parameter (mline = 80, mtab = 8, irunit = 5, ipunit = 6,
     +inunit = 11, iounit = 12, maxflg = 100, mudeck = 1000,
     +mxdeck = 10000, mxdcom = 1000, mxlcom = 10000, marg = 4)
      common / aux / nlin, nlout, ndcomm, nlcomm, nudeck, nuflag,
     +ndflag, nddeck, iffull, idfull, iasepl(marg),
     +isc(mxdcom), ilc(mxdcom), idchec(mudeck,2), ifchec(maxflg)
 
      common / saux / suflag(maxflg+1), sdflag(maxflg+1),
     +sargin(2*marg), inpnam, outnam, smast
      character suflag * (mtab), sdflag * (mtab), sargin * (mline),
     +inpnam * 60, outnam * 60, smast * 1
 
      common / sbuff / sbcomm(mxlcom), sncomm(mxdcom), sddeck(mxdeck),
     +sudeck(mudeck,2)
      character sbcomm * (mline), sncomm * (mtab), sddeck * (mtab),
     +sudeck * (mtab)
*---------*---------*---------*---------*---------*---------*---------*-
      character snin*(*),sname*(mtab),swoff*(mtab)
      save swoff
      logical defaul,switch
 
      data defaul/.true./,switch/.false./
 
      if(nudeck.eq.0)  then
         inrang=defaul
      else
         sname=snin
         call namsrc(sname,sudeck,nudeck,ipos,last)
         if(ipos.gt.0)  then
            inrang=.true.
            swoff=sudeck(ipos,2)
            switch=swoff.ne.sname
         elseif(switch)  then
            inrang=.true.
            switch=swoff.ne.sname
         else
            inrang=.false.
         endif
      endif
      end
      subroutine keepcl(sline,ipos)
*---------*---------*---------*---------*---------*---------*---------*-
      parameter (mline = 80, mtab = 8, irunit = 5, ipunit = 6,
     +inunit = 11, iounit = 12, maxflg = 100, mudeck = 1000,
     +mxdeck = 10000, mxdcom = 1000, mxlcom = 10000, marg = 4)
      common / aux / nlin, nlout, ndcomm, nlcomm, nudeck, nuflag,
     +ndflag, nddeck, iffull, idfull, iasepl(marg),
     +isc(mxdcom), ilc(mxdcom), idchec(mudeck,2), ifchec(maxflg)
 
      common / saux / suflag(maxflg+1), sdflag(maxflg+1),
     +sargin(2*marg), inpnam, outnam, smast
      character suflag * (mtab), sdflag * (mtab), sargin * (mline),
     +inpnam * 60, outnam * 60, smast * 1
 
      common / sbuff / sbcomm(mxlcom), sncomm(mxdcom), sddeck(mxdeck),
     +sudeck(mudeck,2)
      character sbcomm * (mline), sncomm * (mtab), sddeck * (mtab),
     +sudeck * (mtab)
*---------*---------*---------*---------*---------*---------*---------*-
      character*(mline)  sline
 
      if(nlcomm.eq.mxlcom)  then
         call errms1(' ++++++ FATAL ERROR - execution terminated',sline)
         call errms2(' limit of COMDECK buffer reached =',nlcomm)
         call summup
      endif
      nlcomm=nlcomm+1
      if(isc(ipos).eq.0)  isc(ipos)=nlcomm
      ilc(ipos)=nlcomm
      sbcomm(nlcomm)=sline
      end
      function lastnb(string)
      character *(*)  string
      do 10 i=len(string),1,-1
         if(string(i:i).ne.' ')  then
            lastnb=i
            goto 999
         endif
   10 continue
      lastnb=0
  999 end
      function lfirnb(string)
      character string*(*)
 
      do 10 i=1,len(string)
         if(string(i:i).ne.' ') then
            lfirnb=i
            goto 999
         endif
   10 continue
  999 end
      subroutine namsrc(sname,slist,nlist,ipos,last)
*-----------------------------------------------------------------------
*
*   finds name in alphabetic table (binary search).
*
*   Input
*   SNAME           name to be looked up
*   SLIST           table
*   NLIST           length of table
*
*   Output
*   IPOS            = 0: name not in table
*                   > 0: position in table
*   LAST            for IPOS=0, position behind which name belongs
*
*-----------------------------------------------------------------------
      character *(*) sname,slist(*)
      ipos=0
      last=0
      n=nlist
      if(n.gt.0)  then
         kpos=0
   10    m=(n+1)/2
         last=kpos+m
         if (sname.lt.slist(last))  then
            n=m
            last=last-1
            if (n.gt.1) goto 10
         elseif (sname.gt.slist(last))  then
            kpos=last
            n=n-m
            if (n.gt.0) goto 10
         else
            ipos=last
         endif
      endif
      end
      subroutine namtab(sname,slist,nlist,ipos)
*-----------------------------------------------------------------------
*
*   enters a name in an alphabetic table, or gives position if already in.
*
*   input
*   SNAME                   name to be entered
*   SLIST                   name list
*   NUMTAB                  reference list to be updated (integers)
*   NLIST                   no. of names in SLIST
*   Output
*   IPOS                    <0: -pos of name already in table
*                           =0: NLIST <0
*                           >0: pos of newly entered name in table
*
*+++++++++++ IMPORTANT
*   In case the name has been entered, the user must increase the list
*   length himself.
*-----------------------------------------------------------------------
      character *(*) sname,slist(*)
      if(nlist.lt.0)  then
         ipos=0
      elseif(nlist.eq.0)  then
         ipos=1
         slist(1)=sname
      else
         call namsrc(sname,slist,nlist,kpos,last)
         if (kpos.eq.0)  then
*--- name not yet in table
            ipos=last+1
            do 10 i=nlist,ipos,-1
               slist(i+1)=slist(i)
   10       continue
            slist(ipos)=sname
         else
            ipos=-kpos
         endif
      endif
      end
      subroutine pass1
*---------*---------*---------*---------*---------*---------*---------*-
      parameter (mline = 80, mtab = 8, irunit = 5, ipunit = 6,
     +inunit = 11, iounit = 12, maxflg = 100, mudeck = 1000,
     +mxdeck = 10000, mxdcom = 1000, mxlcom = 10000, marg = 4)
      common / aux / nlin, nlout, ndcomm, nlcomm, nudeck, nuflag,
     +ndflag, nddeck, iffull, idfull, iasepl(marg),
     +isc(mxdcom), ilc(mxdcom), idchec(mudeck,2), ifchec(maxflg)
 
      common / saux / suflag(maxflg+1), sdflag(maxflg+1),
     +sargin(2*marg), inpnam, outnam, smast
      character suflag * (mtab), sdflag * (mtab), sargin * (mline),
     +inpnam * 60, outnam * 60, smast * 1
 
      common / sbuff / sbcomm(mxlcom), sncomm(mxdcom), sddeck(mxdeck),
     +sudeck(mudeck,2)
      character sbcomm * (mline), sncomm * (mtab), sddeck * (mtab),
     +sudeck * (mtab)
*---------*---------*---------*---------*---------*---------*---------*-
      character sline*(mline)
      logical ifval,cdflag
      logical ifdecl,eldecl,eidecl,cddecl,dkdecl,cadecl
      ifdecl(sline)=sline(2:3).eq.'IF'.and.(sline(4:4).eq.' '.or.
     +sline(4:4).eq.',')
      eldecl(sline)=sline(2:3).eq.'EL'.and.(sline(4:4).eq.' '.or.
     +sline(4:4).eq.',').or.sline(2:5).eq.'ELSE'.and.(sline(6:6)
     +.eq.' '.or.sline(6:6).eq.',')
      eidecl(sline)=sline(2:3).eq.'EI'.and.(sline(4:4).eq.' '.or.
     +sline(4:4).eq.',').or.sline(2:6).eq.'ENDIF'.and.(sline(7:7)
     +.eq.' '.or.sline(7:7).eq.',')
      dkdecl(sline)=sline(2:3).eq.'DK'.and.(sline(4:4).eq.' '.or.
     +sline(4:4).eq.',').or.sline(2:5).eq.'DECK'.and.(sline(6:6)
     +.eq.' '.or.sline(6:6).eq.',')
      cddecl(sline)=sline(2:3).eq.'CD'.and.(sline(4:4).eq.' '.or.
     +sline(4:4).eq.',').or.sline(2:8).eq.'COMDECK'.and.(sline(9:9)
     +.eq.' '.or.sline(9:9).eq.',')
      cadecl(sline)=sline(2:3).eq.'CA'.and.(sline(4:4).eq.' '.or.
     +sline(4:4).eq.',').or.sline(2:5).eq.'CALL'.and.(sline(6:6)
     +.eq.' '.or.sline(6:6).eq.',')
      call box('pass one starting (comdecks...)')
      nlin=0
      cdflag=.false.
   10 continue
      read(inunit,'(A)',end=999) sline
      nlin=nlin+1
      if(sline(1:1).eq.smast)  then
*--- skip comment lines
        if(sline(2:2).eq.'/')  goto 10
         call upper(sline)
         l=lastnb(sline)
         if(ifdecl(sline))  call flags2(sline,l)
         if(cddecl(sline)) then
            call getapo(sline,2,l,k)
            if(k.eq.0) then
               call errms1(' ++++++ WARNING - empty COMDECK declaration'
     +         ,sline)
               stop 1
               goto 10
            endif
            call getnam(sline,k,l,kfch,klch)
            call namtab(sline(kfch:klch),sncomm,ndcomm,iposd)
            if(iposd.lt.0) then
               call errms1(' ++++++ WARNING - duplicate COMDECK ignored'
     +         ,sline)
               cdflag=.false.
               stop 1
               goto 10
            endif
            if(ndcomm.eq.mxdcom) then
               call errms1(' ++++++ FATAL ERROR - execution terminated',
     +         sline)
               call errms2(' limit of COMDECKs reached =',ndcomm)
               call summup
               stop 1
            endif
            ndcomm=ndcomm+1
            do 20 ii=ndcomm,iposd+1,-1
               isc(ii)=isc(ii-1)
               ilc(ii)=ilc(ii-1)
   20       continue
            isc(iposd)=0
            cdflag=.true.
            ilev=0
         elseif(dkdecl(sline)) then
            cdflag=.false.
         elseif(cdflag) then
            if(ifdecl(sline)) then
               if(ilev.gt.0) then
                  ilev=ilev+1
               elseif(.not.ifval(sline)) then
                  ilev=ilev+1
               endif
            elseif(eldecl(sline)) then
               if(ilev.eq.0) then
                  ilev=1
               elseif(ilev.eq.1) then
                  ilev=0
               endif
            elseif(eidecl(sline)) then
               ilev=max(0,ilev-1)
            elseif(ilev.eq.0) then
               call keepcl(sline,iposd)
            endif
         endif
      elseif(cdflag.and.ilev.eq.0)  then
         call keepcl(sline,iposd)
      endif
      goto 10
  999 end
      subroutine pass2
*---------*---------*---------*---------*---------*---------*---------*-
      parameter (mline = 80, mtab = 8, irunit = 5, ipunit = 6,
     +inunit = 11, iounit = 12, maxflg = 100, mudeck = 1000,
     +mxdeck = 10000, mxdcom = 1000, mxlcom = 10000, marg = 4)
      common / aux / nlin, nlout, ndcomm, nlcomm, nudeck, nuflag,
     +ndflag, nddeck, iffull, idfull, iasepl(marg),
     +isc(mxdcom), ilc(mxdcom), idchec(mudeck,2), ifchec(maxflg)
 
      common / saux / suflag(maxflg+1), sdflag(maxflg+1),
     +sargin(2*marg), inpnam, outnam, smast
      character suflag * (mtab), sdflag * (mtab), sargin * (mline),
     +inpnam * 60, outnam * 60, smast * 1
 
      common / sbuff / sbcomm(mxlcom), sncomm(mxdcom), sddeck(mxdeck),
     +sudeck(mudeck,2)
      character sbcomm * (mline), sncomm * (mtab), sddeck * (mtab),
     +sudeck * (mtab)
*---------*---------*---------*---------*---------*---------*---------*-
      character sline*(mline),sname*(mtab)
      logical ifval,dkflag,inrang
      logical ifdecl,eldecl,eidecl,cddecl,dkdecl,cadecl
      ifdecl(sline)=sline(2:3).eq.'IF'.and.(sline(4:4).eq.' '.or.
     +sline(4:4).eq.',')
      eldecl(sline)=sline(2:3).eq.'EL'.and.(sline(4:4).eq.' '.or.
     +sline(4:4).eq.',').or.sline(2:5).eq.'ELSE'.and.(sline(6:6)
     +.eq.' '.or.sline(6:6).eq.',')
      eidecl(sline)=sline(2:3).eq.'EI'.and.(sline(4:4).eq.' '.or.
     +sline(4:4).eq.',').or.sline(2:6).eq.'ENDIF'.and.(sline(7:7)
     +.eq.' '.or.sline(7:7).eq.',')
      dkdecl(sline)=sline(2:3).eq.'DK'.and.(sline(4:4).eq.' '.or.
     +sline(4:4).eq.',').or.sline(2:5).eq.'DECK'.and.(sline(6:6)
     +.eq.' '.or.sline(6:6).eq.',')
      cddecl(sline)=sline(2:3).eq.'CD'.and.(sline(4:4).eq.' '.or.
     +sline(4:4).eq.',').or.sline(2:8).eq.'COMDECK'.and.(sline(9:9)
     +.eq.' '.or.sline(9:9).eq.',')
      cadecl(sline)=sline(2:3).eq.'CA'.and.(sline(4:4).eq.' '.or.
     +sline(4:4).eq.',').or.sline(2:5).eq.'CALL'.and.(sline(6:6)
     +.eq.' '.or.sline(6:6).eq.',')
      call box('pass two starting (extract)')
      rewind inunit
      nlin=0
      nlout=0
* Eric added this for nagfor -C=undefined
      ilev=0
      dkflag=.false.
   10 continue
      read(inunit,'(A)',end=999) sline
      nlin=nlin+1
      lp=max(1,lastnb(sline))
      if(sline(1:1).eq.smast)  then
*--- skip comment lines
        if(sline(2:2).eq.'/')  goto 10
         call upper(sline)
         l=lastnb(sline)
         if(dkdecl(sline)) then
            call getapo(sline,2,l,k)
            if(k.eq.0) then
               call errms1(' ++++++ WARNING - empty DECK declaration'
     +         ,sline)
               stop 1
               goto 10
            endif
            call getnam(sline,k,l,kfch,klch)
            sname=sline(kfch:klch)
            dkflag=inrang(sname)
            ilev=0
            if(nddeck.eq.mxdeck)  then
               call errms2(
     +         ' ++++++ WARNING - deck table reset at limit =',nddeck)
               nddeck=0
               stop 1
            endif
            nddeck=nddeck+1
            sddeck(nddeck)=sname
         elseif(cddecl(sline)) then
            dkflag=.false.
         elseif(dkflag) then
            if(ifdecl(sline)) then
               if(ilev.gt.0) then
                  ilev=ilev+1
               elseif(.not.ifval(sline)) then
                  ilev=ilev+1
               endif
            elseif(eldecl(sline)) then
               if(ilev.eq.0) then
                  ilev=1
               elseif(ilev.eq.1) then
                  ilev=0
               endif
            elseif(eidecl(sline)) then
               ilev=max(0,ilev-1)
            elseif(ilev.eq.0) then
               if(cadecl(sline))  then
                  call expand(sline)
               else
                  write(iounit,'(A)')  sline(:lp)
                  nlout=nlout+1
               endif
            endif
         endif
      elseif(dkflag.and.ilev.eq.0)  then
         write(iounit,'(A)')  sline(:lp)
         nlout=nlout+1
      endif
      goto 10
  999 end
      subroutine summup
*---------*---------*---------*---------*---------*---------*---------*-
      parameter (mline = 80, mtab = 8, irunit = 5, ipunit = 6,
     +inunit = 11, iounit = 12, maxflg = 100, mudeck = 1000,
     +mxdeck = 10000, mxdcom = 1000, mxlcom = 10000, marg = 4)
      common / aux / nlin, nlout, ndcomm, nlcomm, nudeck, nuflag,
     +ndflag, nddeck, iffull, idfull, iasepl(marg),
     +isc(mxdcom), ilc(mxdcom), idchec(mudeck,2), ifchec(maxflg)
 
      common / saux / suflag(maxflg+1), sdflag(maxflg+1),
     +sargin(2*marg), inpnam, outnam, smast
      character suflag * (mtab), sdflag * (mtab), sargin * (mline),
     +inpnam * 60, outnam * 60, smast * 1
 
      common / sbuff / sbcomm(mxlcom), sncomm(mxdcom), sddeck(mxdeck),
     +sudeck(mudeck,2)
      character sbcomm * (mline), sncomm * (mtab), sddeck * (mtab),
     +sudeck * (mtab)
*---------*---------*---------*---------*---------*---------*---------*-
 
      write(ipunit,10000)
      write(ipunit,10010) nlin,inunit,nlout,iounit
      write(ipunit,10020) 'flags',ndflag
      write(ipunit,10030) (sdflag(i),i=1,ndflag)
      write(ipunit,10020) 'COMDECKs',ndcomm
      write(ipunit,10030) (sncomm(i),i=1,ndcomm)
      write(ipunit,10020) 'DECKs',nddeck
      write(ipunit,10030) (sddeck(i),i=1,nddeck)
10000 format(//' ++++++ SUMMARY ++++++'/)
10010 format(i8,' lines read from input unit ',i4/ i8,
     +' lines written onto unit    ',i4)
10020 format(/' no. of ',a,t20,' :',i5)
10030 format(1x,8a9)
      end
      subroutine termin
      call box('ASTUCE - normal job termination')
      end
      subroutine upper(sl)
************************************************************************
*
*   Converts all characters in SL into upper case.
*
************************************************************************
 
      character * (*) sl
      character alfbet(2) * 26
 
      data alfbet(1) / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      data alfbet(2) / 'abcdefghijklmnopqrstuvwxyz'/
 
      do 10  i=1,len(sl)
        k = index(alfbet(2),sl(i:i))
        if (k .ne. 0)  sl(i:i) = alfbet(1)(k:k)
   10 continue
      end
 
 
 
