      program server
      implicit none

      integer sid ! Server id
      integer cid ! Connection id
      integer n
      integer count
      integer total
      integer PORT
      parameter (PORT = 14999)

      integer mtype
      integer turn, id, gen
      double precision wgt, x, y, z, tx, ty, tz, m, pc, t
      integer aa, zz
      integer ntaccept, ntstart, ntsendp, ntsendeob, ntserver
      integer ntsendeoc, ntshdwn, ntwait, ntend
      integer ntsendnpart, ntnpart
      integer ntsendbrhono, ntbrho

      integer N_PART, N_EOB, N_END, N_IPT
      parameter (N_PART = 1)
      parameter (N_EOB  = 2)
      parameter (N_END  = 3)
      parameter (N_IPT  = 5)

      integer npart
      double precision brho

      ! Initialize fluka connections
      call ntinit()

      sid = ntserver()
      if(sid.eq.-1) then
        write(*,*) "Error creating the server"
        stop
      end if
      ! Listen
      n = ntstart(sid, PORT)
      if(n.eq.-1) then
        write(*,*) "Error cannot start server"
        stop
      end if

      write(*,*) "Listening on port ", n


100   continue ! while(1)
        write(*,*) "** Waiting for new connection"
        cid = ntaccept(sid)

        npart = ntnpart(cid)
        write(*,*) ">>>>>>>>> npart=",npart
        n  = ntbrho(cid, brho)
        write(*,*) ">>>>>>>>> brho=",brho, n

        if(.not.cid.lt.0) then
          write(*,*) "** New connection accepted"
!          n = nttimeout(cid, 20)

            total = 0
            count = 0

200         continue ! while(1)

              n = ntwait(cid,
     &              mtype, id, gen, wgt,
     &              x, y, z, tx, ty, tz, aa, zz, m, pc, t)
              if(n.lt.0) then
                write(*,*) "Client timeout"
                goto 900
              end if

              if(mtype.eq.N_PART) then
                count = count + 1
!                write(*,*) "<"
                n = ntsendp(cid,
     &                id, gen, wgt,
     &                x, y, z, tx, ty, tz, aa, zz, m, pc, t)

                if(n.lt.0) then
                  write(*,*) "Error sending"
                  goto 900
                end if
!                write(*,*) ">"
              else if(mtype.eq.N_IPT) then
                !ipt = id
                turn = gen
              else if(mtype.eq.N_EOB) then
!                write(*,*) "T"
                write(*,*) "  End of turn ", turn, " received ", count,
     &                     " particles"
                n = ntsendeob(cid)
                total = total + count
                count = 0
              else if(mtype.eq.N_END) then
!                write(*,*) "F"
                n = ntsendeoc(cid)
!                write(*,*) "f"
                write(*,*) "** Connection closed"
                goto 900
              end if

            goto 200

900         continue

            n = ntend(cid)
            cid = -1

        end if

      goto 100

! Finished shutdown server
      n = ntshdwn(sid)

      end
