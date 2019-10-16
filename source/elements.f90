! ================================================================================================ !
!  Routines related to the initialisation of elements
!  Need to be in a separate file due to dependency issues
! ================================================================================================ !

! ================================================================================================ !
!  K. Sjobak, A. Santamaria, V.K. Berglyd Olsen, BE-ABP-HSS
!  Created: 2016-12-23
!  Updated: 2019-08-29
!
!  Initialise a lattice element with index elIdx, such as done when reading fort.2 and in DYNK.
!
!  Never delete an element from the lattice, even if it is not making a kick.
!  If the element is not recognized, do nothing (for now).
!  If trying to initialize an element (not lfirst) which is disabled, print an error and exit.
! ================================================================================================ !
subroutine initialise_element(ix,lfirst)

  use crcoall
  use floatPrecision
  use mathlib_bouncer
  use numerical_constants

  use parpro
  use parbeam
  use mod_common
  use mod_common_main
  use mod_common_track
  use mod_utils

  use cheby, only : cheby_kz
  use dynk,  only : dynk_elemData

  integer, intent(in) :: ix
  logical, intent(in) :: lfirst

  integer i,m,k,im,nmz,izu,ibb,ii,j
  real(kind=fPrec) r0,r0a,bkitemp,sfac1,sfac2,sfac2s,sfac3,sfac4,sfac5,crkveb_d,cikveb_d,rho2b_d,   &
    tkb_d,rb_d,rkb_d,xrb_d,zrb_d,cbxb_d,cbzb_d,crxb_d,crzb_d,xbb_d,zbb_d,napx0

  ! Nonlinear Elements
  if(abs(kz(ix)) >= 1 .and. abs(kz(ix)) <= 10) then
    if(.not.lfirst) then
      do i=1,iu
        if(ic(i)-nblo == ix) then
          if(ktrack(i) == 31) goto 100 !ERROR
          sm(ix)  = ed(ix)          ! Also done in envar() which is called from clorb()
          smiv(i) = sm(ix)+smizf(i) ! Also done in program maincr
          smi(i)  = smiv(i)         ! Also done in program maincr
          call setStrack(abs(kz(ix)),i)
        end if
      end do
    end if

  ! Multipoles
  elseif(kz(ix) == 11) then
    if(lfirst) then
      if(abs(el(ix)+one) <= pieni) then
        dki(ix,1) = ed(ix)
        dki(ix,3) = ek(ix)
        ed(ix) = one
        ek(ix) = one
        el(ix) = zero
      else if(abs(el(ix)+two) <= pieni) then
        dki(ix,2) = ed(ix)
        dki(ix,3) = ek(ix)
        ed(ix) = one
        ek(ix) = one
        el(ix) = zero
      end if
    else
      do i=1,iu
        if(ic(i)-nblo == ix) then
          nmz = nmu(ix)
          im  = irm(ix)
          do k=1,nmz
            aaiv(k,i) = scalemu(im)*(ak0(im,k)+amultip(k,i)*aka(im,k))
            bbiv(k,i) = scalemu(im)*(bk0(im,k)+bmultip(k,i)*bka(im,k))
          end do
        end if
      end do
    end if

  ! Cavities (ktrack = 2 for thin)
  elseif(abs(kz(ix)) == 12) then
    dynk_elemData(ix,3) = el(ix)
    phasc(ix) = el(ix)*rad
    el(ix) = zero
    if(lfirst) then
      if(abs(ed(ix)) > pieni .and. abs(ek(ix)) > pieni) then
        ncy2   = ncy2 + 1
        kp(ix) = 6
      end if
    else
      hsyc(ix) = ((twopi)*ek(ix))/tlen                             ! SYNC block
      hsyc(ix) = (c1m3*hsyc(ix)) * real(sign(1,kz(ix)),kind=fPrec) ! trauthin/trauthck
    end if

  ! Wire
  else if(kz(ix) == 15) then
    ed(ix) = zero
    ek(ix) = zero
    el(ix) = zero

  ! BEAM-BEAM
  elseif(kz(ix) == 20) then

    if(nbeam == 0 .and. .not. lfirst) then
      write(lerr,"(a)") "BEAMBEAM> ERROR Beam-beam element encountered, but no BEAM block in '"//trim(fort3)//"'"
      call prror
    end if

    if(lfirst) then
      ptnfac(ix)  = el(ix)
      el(ix)      = zero
      parbe(ix,5) = ed(ix)
      ed(ix)      = zero
      parbe(ix,6) = ek(ix)
      ek(ix)      = zero
    end if
    ! This is to inialize all the beam-beam element before the tracking (or to update it for DYNK).
    if(.not.lfirst) then
      do i=1,iu
        if(ic(i)-nblo == ix) then
          ibb = imbb(i)
          if(parbe(ix,2) > zero) then
            if(beam_expflag == 1) then
              bbcu(ibb,1)  = parbe(ix,7)
              bbcu(ibb,4)  = parbe(ix,8)
              bbcu(ibb,6)  = parbe(ix,9)
              bbcu(ibb,2)  = parbe(ix,10)
              bbcu(ibb,9)  = parbe(ix,11)
              bbcu(ibb,10) = parbe(ix,12)
              bbcu(ibb,3)  = parbe(ix,13)
              bbcu(ibb,5)  = parbe(ix,14)
              bbcu(ibb,7)  = parbe(ix,15)
              bbcu(ibb,8)  = parbe(ix,16)
              do ii=1,10
                bbcu(ibb,ii) = bbcu(ibb,ii)*c1m6
              end do
            end if
            ktrack(i)   = 44
            parbe(ix,4) = (((-one*crad)*ptnfac(ix))*half)*c1m6
            if(ibeco == 1) then
              track6d(1,1) = parbe(ix,5)*c1m3
              track6d(2,1) = zero
              track6d(3,1) = parbe(ix,6)*c1m3
              track6d(4,1) = zero
              track6d(5,1) = zero
              track6d(6,1) = zero
              napx0 = napx
              napx  = 1
              call beamint(napx,track6d,parbe,sigz,bbcu,imbb(i),ix,ibtyp,ibbc, mtc)
              beamoff(1,imbb(i)) = track6d(1,1)*c1e3
              beamoff(2,imbb(i)) = track6d(3,1)*c1e3
              beamoff(3,imbb(i)) = track6d(5,1)*c1e3
              beamoff(4,imbb(i)) = track6d(2,1)*c1e3
              beamoff(5,imbb(i)) = track6d(4,1)*c1e3
              beamoff(6,imbb(i)) = track6d(6,1)
              napx = napx0
            end if

          else if(parbe(ix,2) == zero) then
            if(beam_expflag == 1) then
              bbcu(ibb,1) = parbe(ix,1)
              bbcu(ibb,2) = parbe(ix,3)
              bbcu(ibb,3) = parbe(ix,13)
            end if
            if(ibbc == 1) then
              sfac1  = bbcu(ibb,1)+bbcu(ibb,2)
              sfac2  = bbcu(ibb,1)-bbcu(ibb,2)
              sfac2s = one
              if(sfac2 < zero) sfac2s = -one
              sfac3 = sqrt(sfac2**2+(four*bbcu(ibb,3))*bbcu(ibb,3))
              if(sfac3 > sfac1) then
                write(lerr,"(a)") "BEAMBEAM> ERROR 6D beam-beam with tilt not possible."
                call prror
              end if
              sfac4 = (sfac2s*sfac2)/sfac3
              sfac5 = (((-one*sfac2s)*two)*bbcu(ibb,3))/sfac3
              sigman(1,ibb) = sqrt(((sfac1+sfac2*sfac4)+(two*bbcu(ibb,3))*sfac5)*half)
              sigman(2,ibb) = sqrt(((sfac1-sfac2*sfac4)-(two*bbcu(ibb,3))*sfac5)*half)
              bbcu(ibb,11)  = sqrt(half*(one+sfac4))
              bbcu(ibb,12)  = (-one*sfac2s)*sqrt(half*(one-sfac4))
              if(bbcu(ibb,3) < zero) bbcu(ibb,12) = -one*bbcu(ibb,12)
            else
              bbcu(ibb,11)  = one
              sigman(1,ibb) = sqrt(bbcu(ibb,1))
              sigman(2,ibb) = sqrt(bbcu(ibb,2))
            end if

            ! Round beam
            nbeaux(imbb(i)) = 0
            if(sigman(1,imbb(i)) == sigman(2,imbb(i))) then
              if(nbeaux(imbb(i)) == 2 .or. nbeaux(imbb(i)) == 3) then
                write(lerr,"(a)") "BEAMBEAM> ERROR At each interaction point the beam must be either "//&
                  "round or elliptical for all particles"
                call prror
              else
                nbeaux(imbb(i)) = 1
                sigman2(1,imbb(i)) = sigman(1,imbb(i))**2
              end if
            end if

            ! Elliptic beam x>z
            if(sigman(1,imbb(i)) > sigman(2,imbb(i))) then
              if(nbeaux(imbb(i)) == 1 .or. nbeaux(imbb(i)) == 3) then
                write(lerr,"(a)") "BEAMBEAM> ERROR At each interaction point the beam must be either "//&
                  "round or elliptical for all particles"
                call prror
              else
                nbeaux(imbb(i)) = 2
                ktrack(i)       = 42
                sigman2(1,imbb(i)) = sigman(1,imbb(i))**2
                sigman2(2,imbb(i)) = sigman(2,imbb(i))**2
                sigmanq(1,imbb(i)) = sigman(1,imbb(i))/sigman(2,imbb(i))
                sigmanq(2,imbb(i)) = sigman(2,imbb(i))/sigman(1,imbb(i))
              end if
            end if

            ! Elliptic beam z>x
            if(sigman(1,imbb(i)) < sigman(2,imbb(i))) then
              if(nbeaux(imbb(i)) == 1 .or. nbeaux(imbb(i)) == 2) then
                write(lerr,"(a)") "BEAMBEAM> ERROR At each interaction point the beam must be either "//&
                  "round or elliptical for all particles"
                call prror
              else
                nbeaux(imbb(i)) = 3
                ktrack(i)       = 43
                sigman2(1,imbb(i)) = sigman(1,imbb(i))**2
                sigman2(2,imbb(i)) = sigman(2,imbb(i))**2
                sigmanq(1,imbb(i)) = sigman(1,imbb(i))/sigman(2,imbb(i))
                sigmanq(2,imbb(i)) = sigman(2,imbb(i))/sigman(1,imbb(i))
              end if
            end if

            strack(i) = crad*ptnfac(ix)
            if(ibbc == 0) then
              crkveb_d = parbe(ix,5)
              cikveb_d = parbe(ix,6)
            else
              crkveb_d = parbe(ix,5)*bbcu(imbb(i),11)+parbe(ix,6)*bbcu(imbb(i),12)
              cikveb_d = parbe(ix,6)*bbcu(imbb(i),11)-parbe(ix,5)*bbcu(imbb(i),12)
            end if

            if(nbeaux(imbb(i)) == 1) then
              ktrack(i) = 41
              if(ibeco == 1) then
                rho2b_d = crkveb_d**2+cikveb_d**2
                tkb_d   = rho2b_d/(two*sigman2(1,imbb(i)))
                beamoff(4,imbb(i)) = ((strack(i)*crkveb_d)/rho2b_d)*(one-exp_mb(-one*tkb_d))
                beamoff(5,imbb(i)) = ((strack(i)*cikveb_d)/rho2b_d)*(one-exp_mb(-one*tkb_d))
              end if
            end if

            if(ktrack(i) == 42) then
              if(ibeco == 1) then
                rb_d  = sqrt(two*(sigman2(1,imbb(i))-sigman2(2,imbb(i))))
                rkb_d = (strack(i)*pisqrt)/rb_d
                xrb_d = abs(crkveb_d)/rb_d
                zrb_d = abs(cikveb_d)/rb_d
                tkb_d = (crkveb_d**2/sigman2(1,imbb(i))+cikveb_d**2/sigman2(2,imbb(i)))*half
                xbb_d = sigmanq(2,imbb(i))*xrb_d
                zbb_d = sigmanq(1,imbb(i))*zrb_d
                if(ibtyp == 0) then
                  call errf(xrb_d,zrb_d,crxb_d,crzb_d)
                  call errf(xbb_d,zbb_d,cbxb_d,cbzb_d)
                else
                  call wzsub(xrb_d,zrb_d,crxb_d,crzb_d)
                  call wzsub(xbb_d,zbb_d,cbxb_d,cbzb_d)
                end if
                beamoff(4,imbb(i)) = (rkb_d*(crzb_d-exp_mb(-one*tkb_d)*cbzb_d))*sign(one,crkveb_d)
                beamoff(5,imbb(i)) = (rkb_d*(crxb_d-exp_mb(-one*tkb_d)*cbxb_d))*sign(one,cikveb_d)
              end if
            end if

            if(ktrack(i) == 43) then
              if(ibeco == 1) then
                rb_d  = sqrt(two*(sigman2(2,imbb(i))-sigman2(1,imbb(i))))
                rkb_d = (strack(i)*pisqrt)/rb_d
                xrb_d = abs(crkveb_d)/rb_d
                zrb_d = abs(cikveb_d)/rb_d
                tkb_d = (crkveb_d**2/sigman2(1,imbb(i))+cikveb_d**2/sigman2(2,imbb(i)))*half
                xbb_d = sigmanq(2,imbb(i))*xrb_d
                zbb_d = sigmanq(1,imbb(i))*zrb_d
                if(ibtyp == 0) then
                  call errf(zrb_d,xrb_d,crzb_d,crxb_d)
                  call errf(zbb_d,xbb_d,cbzb_d,cbxb_d)
                else
                  call wzsub(zrb_d,xrb_d,crzb_d,crxb_d)
                  call wzsub(zbb_d,xbb_d,cbzb_d,cbxb_d)
                end if
                beamoff(4,imbb(i)) = (rkb_d*(crzb_d-exp_mb(-one*tkb_d)*cbzb_d))*sign(one,crkveb_d)
                beamoff(5,imbb(i)) = (rkb_d*(crxb_d-exp_mb(-one*tkb_d)*cbxb_d))*sign(one,cikveb_d)
              end if
            end if
          end if
        end if
      end do
    end if

  ! Crab Cavities
  ! Note: If setting something else than el(),
  ! DON'T call initialise_element on a crab, it will reset the phase to 0.
  elseif(abs(kz(ix)) == 23) then
    crabph(ix) = el(ix)
    el(ix)     = zero

  ! CC Mult kick order 2
  elseif(abs(kz(ix)) == 26) then
    crabph2(ix) = el(ix)
    el(ix)      = zero

  ! CC Mult kick order 3
  elseif(abs(kz(ix)) == 27) then
    crabph3(ix) = el(ix)
    el(ix)      = zero

  ! CC Mult kick order 4
  else if(abs(kz(ix)) == 28) then
    crabph4(ix) = el(ix)
    el(ix)      = zero

  ! e-lens
  else if(kz(ix) == 29) then
    ed(ix) = zero
    ek(ix) = zero
    el(ix) = zero

  ! Chebyshev lens
  else if(kz(ix) == cheby_kz) then
    ed(ix) = zero
    ek(ix) = zero
    el(ix) = zero
  end if

  return

  ! Error handlers
100 continue
  write(lerr,"(a)") "INITELEM> ERROR Tried to set the strength of an element which is disabled. bez = '",trim(bez(ix)),"'"
  call prror

end subroutine initialise_element

! ================================================================================================ !
!  Calculate strack for magnet types
!  Code merged from include files stra01.f90 to stra14.f90
!  V.K. Berglyd Olsen, BE-ABP-HSS
!  Updated: 2019-09-20
! ================================================================================================ !
subroutine setStrack(skz, i)

  use parpro
  use crcoall
  use mod_common,       only : tiltc, tilts, ic, dki
  use mod_common_main,  only : smiv
  use mod_common_track, only : strack, strackc, stracks
  use numerical_constants

  integer, intent(in) :: skz
  integer, intent(in) :: i

  integer ix

  ix = ic(i) - nblo

  select case(skz)
  case(1)
    strack(i) = smiv(i)*c1e3
  case(2)
    strack(i) = smiv(i)
  case(3)
    strack(i) = smiv(i)*c1m3
  case(4)
    strack(i) = smiv(i)*c1m6
  case(5)
    strack(i) = smiv(i)*c1m9
  case(6)
    strack(i) = smiv(i)*c1m12
  case(7)
    strack(i) = smiv(i)*c1m15
  case(8)
    strack(i) = smiv(i)*c1m18
  case(9)
    strack(i) = smiv(i)*c1m21
  case(10)
    strack(i) = smiv(i)*c1m24
  case(11)
    strack(i) = dki(ix,1)/dki(ix,3)
  case(12)
    strack(i) = dki(ix,1)
  case(13)
    strack(i) = dki(ix,2)/dki(ix,3)
  case(14)
    strack(i) = dki(ix,2)
  case default
    write(lerr,"(a,i0,a)") "TRACKING> ERROR Setting strack for type ",skz," not possible. This is a bug."
    call prror
  end select

  strackc(i) = strack(i)*tiltc(i)
  stracks(i) = strack(i)*tilts(i)

end subroutine setStrack
