      program client
      implicit none

      integer PORT
      integer NPARTS, NTURNS
      parameter (PORT = 14999)
      parameter (NPARTS = 20000, NTURNS = 150)

      integer i, j, n
      integer cid

      integer mtype
      integer turn, id, gen
      double precision wgt, x, y, z, tx, ty, tz, m, pc, t
      integer aa, zz
      integer ntconnect, ntsendp, ntsendeob, ntsendipt
      integer ntsendeoc, ntrecv, ntwait, ntend
      integer ntsendnpart, ntnpart
      integer ntsendbrhono, ntbrho

      integer N_PART, N_EOB, N_END
      parameter (N_PART = 1)
      parameter (N_EOB = 2)
      parameter (N_END  = 3)

      integer npart
      double precision brho

      turn = 0
      id   = 0
      gen  = 1
      wgt  = 1.0D+0
      x    = -6.319439221D-01
      y    = 1.147679085D+00
      z    = 1.0D+0
      tx   = -3.856529999D-05
      ty   = -2.004572999D-06
      tz   =  5.004483739D-06
      aa   = 1
      zz   = 1
      m    = 450.0D+0
      pc   = 4520.0D+0
      t   = 5.0D-3

      call ntinit()
      cid = ntconnect("localhost", PORT)
      if (cid.eq.-1) then
        write(*,*) "Error connecting to server"
        goto 1000
      end if

      write(*,*) "Connected to server"

      npart = 56000
      n = ntsendnpart(cid, npart)
      brho = 1.2345d0
      n = ntsendbrhono(cid, brho)

      do j = 1,NTURNS

        turn = j
c       Send ipt
        write(*,*) "IPT ", j
        n = ntsendipt(cid, j, 0)
        if(n.lt.0) then
          write(*,*) "Error sending IPT"
          goto 900
        end if

        do i = 1,NPARTS
          id = i

c         Send particle
          n = ntsendp(cid, id, gen, wgt,
     &                x, y, z, tx, ty, tz, aa, zz, m, pc, t)
          if(n.lt.0) then
            write(*,*) "Error sending P"
            write(*,*) id, wgt, x, y, z, tx, ty, tz, aa, zz, m, pc, t
            goto 900
          end if

c         Read all incoming particles to release server out buffer
100       continue
          n = ntrecv(cid, mtype, id, gen, wgt,
     &               x, y, z, tx, ty, tz, aa, zz, m, pc, t)
          if(n.gt.0) then
            goto 100
          end if
c         End of reading incoming particles

        end do ! Done sending current turn particles

c       Send end of turn
        write(*,*) "End of turn ", j
        n = ntsendeob(cid)
        if(n.lt.0) then
          write(*,*) "Error sending T"
          goto 900
        end if

c     + Wait until end of turn (Synchronize)
c     |
200     continue
          n = ntwait(cid, mtype, id, gen, wgt,
     &               x, y, z, tx, ty, tz, aa, zz, m, pc, t)
          if(n.eq.-1) then
            write(*,*) "Server timed out while waiting end of turn"
            goto 900
          end if
          if(mtype.ne.N_EOB) then
            goto 200
          end if
c     |
c     + Finished waiting end of turn


      end do ! Turn loop

      write(*,*) "Done sending particles"
      write(*,*) "Closing connection..."

c     Send end of computation
      n = ntsendeoc(cid)
      if(n.lt.0) then
        write(*,*) "Error sending F"
        goto 900
      end if

c     Wait end of comp
      n = ntwait(cid, mtype, id, gen, wgt,
     &           x, y, z, tx, ty, tz, aa, zz, m, pc, t)
      if(n.eq.-1) then
        write(*,*) "Server timed out while waiting end of comp"
        goto 900
      end if
      if(mtype.ne.N_END) then
        write(*,*) "Unexpected message received"
      end if
c     Both ends agreed to disconnect

900   continue ! Error while connected

c     Finish connection
      n = ntend(cid)

1000  continue ! Error connecting

      stop
      end
