#ifdef FLUKA

subroutine check_coupling_integrity
!
!-----------------------------------------------------------------------
!     A.Mereghetti, for the FLUKA Team
!     last modified: 17-07-2013
!     check that an entrance MARKER is followed by an exit one in the
!        accelerator sequence, even though not immediately
!     inserted in main code by the 'fluka' compilation flag
!-----------------------------------------------------------------------

      use floatPrecision
      use mod_fluka
      use crcoall
      use parpro
      use mod_common
      use mod_commont
      use aperture
      implicit none

!     temporary variables
      integer i1 , i2
      integer ix1, ix2
      integer istart, istop
      logical lerr, lfound, lcurturn

      lerr = .false.

      write(lout,*) ''
      write(lout,10010)
      write(lout,*) ' CALL TO CHECK_COUPLING_INTEGRITY '
      write(lout,*) ' NB: only entrance/exit markers are listed;'
      write(lout,*) '     a single entry is by definition righteous'
      write(lout,10010)
      write(lout,*) ''
      write(lout,*) ''
      write(lout,*) '        keys to FLUKA types:'
      write(lout,*) FLUKA_ELEMENT,' --> simple element'
      write(lout,*) FLUKA_ENTRY,' --> entrance point'
      write(lout,*) FLUKA_EXIT,' --> exit point'
      write(lout,*) ''
      write(lout,*) ''

      i1=1
      do while ( i1.le.iu )
        if(ktrack(i1).ne.1.and.ic(i1).gt.nblo) then
!         SINGLE ELEMENT
          ix1=ic(i1)-nblo
          if ( fluka_type(ix1).eq.FLUKA_ENTRY ) then
            write(lout,*) ''
            write(lout,*) ''
            write(lout,10020) 'entry type', 'name', 'ID SING EL ID struct', 'ID geom'
            write(lout,10030) fluka_type(ix1), bez(ix1), ix1, i1, fluka_geo_index(ix1)
            istart = i1+1
            istop  = iu
            lcurturn = .true.
            lfound = .false.
 1982       continue
            do i2=istart,istop
              if(ktrack(i2).ne.1.and.ic(i2).gt.nblo) then
!               SINGLE ELEMENT
                ix2=ic(i2)-nblo
                if ( fluka_type(ix2).eq.FLUKA_EXIT ) then
                  if(fluka_geo_index(ix1).eq.fluka_geo_index(ix2))then
                    write(lout,10030) fluka_type(ix2), bez(ix2), ix2, i2, fluka_geo_index(ix2)
                    call contour_aperture_markers( i1, i2, .true.  )
                    i1 = i2 + 1
                    lfound = .true.
                    if ( lcurturn ) then
                      exit
                    else
                      goto 1983
                    endif
                  else
                    write(lout,10030) fluka_type(ix2), bez(ix2), ix2, i2, fluka_geo_index(ix2)
                    write(lout,*) 'ERROR! un-matched geo index'
                    write(lout,*) ''
                    lerr = .true.
                  endif
                elseif ( fluka_type(ix2).ne.FLUKA_NONE ) then
                  write(lout,*) 'ERROR! non-exit point when entrance is on'
                  write(lout,*) ''
                  write(lout,10030) fluka_type(ix2), bez(ix2), ix2, i2, fluka_geo_index(ix2)
                  lerr = .true.
                endif
              endif
            enddo
            if ( .not. lfound ) then
              if ( lcurturn ) then
!               the exit marker was not found: maybe the fluka insertion
!                  is across the end/beginning of the accelerator structure;
!               -> restart the research, upstream of the entrance marker:
                istart = 1
                istop  = i1
                lcurturn = .false.
                goto 1982
              else
!               failing research:
!               NB: in principle, this should never happen, but let's be picky
                write(lout,*)'ERROR! entrance point does not have the exit'
                write(lout,*)''
                lerr = .true.
              endif
            endif
          endif
        endif

!       go to next accelerator entry
        i1 = i1+1
      enddo

 1983 continue
      if ( lerr ) then
        write(lout,*) ' at least one inconsistency in flagging elements'
        write(lout,*) '    for coupling: please check carefully...'
        call prror(-1)
      endif

!     au revoir:
      return
10010 format(132('-'))
10020 format(1X,A10,1X,A4,12X,3(1X,A10))
10030 format(1X,I10,1X,A16,3(1X,I10))
end subroutine check_coupling_integrity

subroutine kernel_fluka_element( nturn, i, ix )
!
!-----------------------------------------------------------------------
!     A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!     last modified: 07-02-2014
!     'transport' subroutine, for a Fluka insertion corresponding to a
!       single SINGLE ELEMENT, of non-zero length
!     inserted in main code by the 'fluka' compilation flag
!-----------------------------------------------------------------------
!
      use floatPrecision
      use mod_fluka
      use numerical_constants, only : zero, one, c1e3, c1m3
      use crcoall
      use parpro
      use mod_common
      use mod_commont
      use mod_commonmn

      use mod_hions

#ifdef ROOT
      use root_output
#endif

      implicit none

!     interface variables:
      integer nturn, i, ix

!     temporary variables
      integer ret, j, k
      real(kind=fPrec) eltot
      integer pid_q               ! ph: hisix
      save

      eltot = fluka_synch_length( ix )

      if (fluka_debug) then
!        where are we?
         write(fluka_log_unit,*)'# In fluka element of type ', fluka_type(ix)
         write(fluka_log_unit,*)'#   i=', i
         write(fluka_log_unit,*)'#   ix=', ix
         write(fluka_log_unit,*)'#   fluka_geo_index=',fluka_geo_index(ix)
         write(fluka_log_unit,*)'#   eltot=',eltot
      end if

!     PH hisix compute the number of nucleons sent to FLUKA
!     PH hisix compute the total ion energy sent to FLUKA
      nnuc0 = 0
      ien0  = zero
      do j=1,napx
        nnuc0   = nnuc0 + naa(j)
        ien0    = ien0 + ejv(j)
        ! array of particle ids sent to FLUKA
        pids(j) = fluka_uid(j)
      end do


      ret = fluka_send_receive( nturn, fluka_geo_index(ix), eltot, napx, xv1, xv2, yv1, yv2, sigmv, ejv, naa, nzz, nucm )

      if (ret.eq.-1) then
         write(lout,*)'[Fluka] Error in Fluka communication in kernel_fluka_element...'
         write(fluka_log_unit,*)'# Error in Fluka communication in kernel_fluka_element...'
         call prror(-1)
      end if

      nnuc1 = 0                 ! hisix: number of nucleons leaving the collimato
      ien1  = zero              ! hisix: total energy leaving the collimator
!     particles to be tracked
      do j=1,napx
!        Update values related to losses
         partID(j) = j
         pstop (j) = .false.
!        Update variables depending on total energy
!         ejfv  (j) = sqrt((ejv(j)-pma)*(ejv(j)+pma))
!         rvv   (j) = (ejv(j)*e0f)/(e0*ejfv(j))
!         dpsv  (j) = (ejfv(j)-e0f)/e0f
!         oidpsv(j) = one/(one+dpsv(j))
         ejfv  (j) = sqrt((ejv(j)-nucm(j))*(ejv(j)+nucm(j)))   ! hisix: ion mass
         rvv   (j) = (ejv(j)*e0f)/(e0*ejfv(j))                 ! hisix: remains unchanged
         dpsv  (j) = (ejfv(j)*(nucm0/nucm(j))-e0f)/e0f         ! hisix: new delta
         oidpsv(j) = one/(one+dpsv(j))
         dpsv1 (j) = (dpsv(j)*c1e3)*oidpsv(j)
         mtc     (j) = (nzz(j)*nucm0)/(zz0*nucm(j))            ! hisix: mass to charge
         moidpsv (j) = mtc(j)*oidpsv(j)                        ! hisix
         omoidpsv(j) = c1e3*((one-mtc(j))*oidpsv(j))           ! hisix
         nnuc1       = nnuc1 + naa(j)                          ! outcoming nucleons
         ien1        = ien1  + ejv(j)                          ! outcoming energy
      end do

!     hisix: compute the nucleon and energy difference
!              reduce by factor 1e-3 to get the energy in GeV
      if((ien0-ien1).gt.one) then
        write(208,*) fluka_geo_index(ix), nnuc0-nnuc1, c1m3*(ien0-ien1)
#ifdef ROOT
        if(root_flag .and. root_FLUKA .eq. 1) then
          call root_FLUKA_EnergyDeposition(fluka_geo_index(ix), nnuc0-nnuc1, c1m3*(ien0-ien1))
        end if
#endif
      ! hisix debugging:
      ! write out the particle distribution after the primary
        if(fluka_geo_index(ix).eq.11) then
          do j=1,napx
            write(210,*) naa(j), nzz(j), nucm(j),ejfv(j),mtc(j),dpsv(j)
          end do
        end if
      end if

!     hisix: check which particle ids have not been sent back
!            write their ids to fort.209
      do j=1,npart                                             ! loop over all pids possible
        pid_q = zero

        do k=1,napx                                            ! loop over pids received
          if(pids(j).eq.fluka_uid(k)) then
            pid_q = one
          end if
        end do

        if(pid_q.eq.zero.and.pids(j).ne.zero) then
          write(209,*) fluka_geo_index(ix), pids(j)
        end if
      end do

!     empty places
      do j=napx+1,npart
!        Update values related to losses
         partID(j) = j
         pstop (j) = .true.
!        Update values related to momentum
         ejv   (j) = zero
         rvv   (j) = one
         ejfv  (j) = zero
         dpsv  (j) = zero
         oidpsv(j) = one
         dpsv1 (j) = zero
         mtc   (j) = one            ! hiSix
         naa   (j) = aa0            ! hiSix
         nzz   (j) = zz0            ! hiSix
         nucm  (j) = nucm0          ! hiSix
         moidpsv (j) = one          ! hiSix
         omoidpsv(j) = zero         ! hiSix
      end do

!     au revoir:
      return
end subroutine kernel_fluka_element

subroutine kernel_fluka_entrance( nturn, i, ix )
!
!-----------------------------------------------------------------------
!     A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!     last modified: 07-02-2014
!     'transport' subroutine, for the marker starting a Fluka insertion
!     inserted in main code by the 'fluka' compilation flag
!-----------------------------------------------------------------------
!
      use floatPrecision
      use mod_fluka
      use numerical_constants, only : zero
      use crcoall
      use parpro
      use mod_common
      use mod_commonmn

      use mod_hions

      implicit none


!     interface variables:
      integer nturn, i, ix

!     temporary variables
      integer ret, j
      real(kind=fPrec) eltot
      integer k                   ! ph: hisix
      integer pid_q               ! ph: hisix

      save

      eltot = zero

      if (fluka_debug) then
!        where are we?
         write(fluka_log_unit,*)'# In fluka element of type ', fluka_type(ix)
         write(fluka_log_unit,*)'#   i=', i
         write(fluka_log_unit,*)'#   ix=', ix
         write(fluka_log_unit,*)'#   fluka_geo_index=',fluka_geo_index(ix)
         write(fluka_log_unit,*)'#   eltot=',eltot
      end if

      ! P. HERMES for hisix
      ! send also A,Z,m to FLUKA

!     PH hisix compute the number of nucleons sent to FLUKA
!     PH hisix compute ion energy sent to FLUKA
!     PH initialize array of particle ids
      nnuc0 = 0
      ien0  = zero
      do j=1,npart
        pids(j) = 0
      end do

      do j=1,napx
        nnuc0   = nnuc0 + naa(j)
        ien0    = ien0  + ejv(j)
        pids(j) = fluka_uid(j)   ! array of particle ids sent to FLUKA
!    write(*,*),'PH:',pids(j)
      end do

      ret = fluka_send( nturn, fluka_geo_index(ix), eltot, napx, xv1, xv2, yv1, yv2, sigmv, ejv, naa, nzz, nucm )

      if (ret.eq.-1) then
         write(lout,*)'[Fluka] Error in Fluka communication in kernel_fluka_entrance...'
         write(fluka_log_unit,*)'# Error in Fluka communication in kernel_fluka_entrance...'
         call prror(-1)
      end if

!     au revoir:
      return
end subroutine kernel_fluka_entrance

subroutine kernel_fluka_exit( nturn, i, ix )
!
!-----------------------------------------------------------------------
!     A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!     last modified: 07-02-2014
!     'transport' subroutine, for the marker closing a Fluka insertion
!     inserted in main code by the 'fluka' compilation flag
!-----------------------------------------------------------------------
!
      use floatPrecision
      use mod_fluka
      use numerical_constants, only : zero, one, c1e3, c1m3
      use crcoall
      use parpro
      use mod_common
      use mod_commont
      use mod_commonmn

      use mod_hions

#ifdef ROOT
      use root_output
#endif

      implicit none

!     interface variables:
      integer nturn, i, ix

!     temporary variables
      integer ret, j, k
      real(kind=fPrec) eltot
      integer pid_q               ! ph: hisix

      save

      eltot = fluka_synch_length( ix )

      if (fluka_debug) then
!        where are we?
         write(fluka_log_unit,*)'# In fluka element of type ', fluka_type(ix)
         write(fluka_log_unit,*)'#   i=', i
         write(fluka_log_unit,*)'#   ix=', ix
         write(fluka_log_unit,*)'#   fluka_geo_index=',fluka_geo_index(ix)
         write(fluka_log_unit,*)'#   eltot=',eltot
      end if

      ret = fluka_receive( nturn, fluka_geo_index(ix), eltot, napx, xv1, xv2, yv1, yv2, sigmv, ejv, naa, nzz, nucm )

      if (ret.eq.-1) then
         write(lout,*)'[Fluka] Error in Fluka communication in kernel_fluka_exit...'
         write(fluka_log_unit,*)'# Error in Fluka communication in kernel_fluka_exit...'
         call prror(-1)
      end if

      nnuc1 = 0                 ! hisix: number of nucleons leaving the collimator
      ien1  = zero              ! hisix: total energy leaving the collimator
!     particles to be tracked
      do j=1,napx
!        Update values related to losses
         partID(j) = j
         pstop (j) = .false.
!        Update variables depending on total energy
!         ejfv  (j) = sqrt((ejv(j)-pma)*(ejv(j)+pma))
!         rvv   (j) = (ejv(j)*e0f)/(e0*ejfv(j))
!         dpsv  (j) = (ejfv(j)-e0f)/e0f
!         oidpsv(j) = one/(one+dpsv(j))
!         dpsv1 (j) = (dpsv(j)*c1e3)*oidpsv(j)
         ejfv  (j) = sqrt((ejv(j)-nucm(j))*(ejv(j)+nucm(j)))   ! hisix: ion mass
         rvv   (j) = (ejv(j)*e0f)/(e0*ejfv(j))                 ! hisix: remains unchanged
         dpsv  (j) = (ejfv(j)*(nucm0/nucm(j))-e0f)/e0f         ! hisix: new delta
         oidpsv(j) = one/(one+dpsv(j))
         dpsv1 (j) = (dpsv(j)*c1e3)*oidpsv(j)
         mtc     (j) = (nzz(j)*nucm0)/(zz0*nucm(j))            ! hisix: mass to charge
         moidpsv (j) = mtc(j)*oidpsv(j)                        ! hisix
         omoidpsv(j) = c1e3*((one-mtc(j))*oidpsv(j))           ! hisix
         nnuc1       = nnuc1 + naa(j)                          ! outcoming nucleons
         ien1        = ien1  + ejv(j)                          ! outcoming energy
      end do

!       hisix: compute the nucleon and energy difference
!              reduce by factor 1e-3 to get the energy in GeV
        if((ien0-ien1).gt.one) then
          write(208,*) fluka_geo_index(ix), nnuc0-nnuc1, c1m3*(ien0-ien1)
#ifdef ROOT
          if(root_flag .and. root_FLUKA .eq. 1) then
            call root_FLUKA_EnergyDeposition(fluka_geo_index(ix), nnuc0-nnuc1, c1m3*(ien0-ien1))
          end if
#endif

          ! hisix debugging:
          ! write out the particle distribution after the primary
          if (fluka_geo_index(ix).eq.11) then
            do j=1,napx
              write(210,*) naa(j), nzz(j), nucm(j),ejfv(j),mtc(j),dpsv(j)
            end do
          end if
        end if
!
!     hisix: check which particle ids have not been sent back
!            write their ids to fort.209
      do j=1,npart                                       ! loop over all pids possible
        pid_q = zero
        do k=1,napx                                    ! loop over pids received
          if(pids(j).eq.fluka_uid(k)) then
            pid_q = one
          end if
        end do
        if(pid_q.eq.zero.and.pids(j).ne.zero) then
          write(209,*) fluka_geo_index(ix), pids(j)
        end if
      end do

!     empty places
      do j=napx+1,npart
!        Update values related to losses
         partID(j) = j
         pstop (j) = .true.
!        Update values related to momentum
         ejv   (j) = zero
         rvv   (j) = one
         ejfv  (j) = zero
         dpsv  (j) = zero
         oidpsv(j) = one
         dpsv1 (j) = zero
         mtc   (j) = one            ! hiSix
         naa   (j) = aa0            ! hiSix
         nzz   (j) = zz0            ! hiSix
         nucm  (j) = nucm0          ! hiSix
         moidpsv (j) = one          ! hiSix
         omoidpsv(j) = zero         ! hiSix
      end do

!     au revoir:
      flush(208)
      flush(209)

      return
end subroutine kernel_fluka_exit

#endif
