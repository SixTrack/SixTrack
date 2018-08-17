! start include/multl12.f90
#ifndef TILT
  if(iv.eq.2) ekk=(((((bb(4)+four*(bb(5)*cr(2)+aa(5)*ci(2)))+10d0*(bb(6)*cr(3)+aa(6)*ci(3)))+20d0*&
                  (bb(7)*cr(4)+aa(7)*ci(4)))+35d0*(bb(8)*cr(5)+aa(8)*ci(5)))+56d0*(bb(9)*cr(6)+aa(9)*ci(6)))+84d0*&
                  (bb(10)*cr(7)+aa(10)*ci(7))
  if(iv.eq.3) ekk=(((bb(6)+6d0*(bb(7)*cr(2)+aa(7)*ci(2)))+21d0*(bb(8)*cr(3)+aa(8)*ci(3)))+56d0*&
                  (bb(9)*cr(4)+aa(9)*ci(4)))+126d0*(bb(10)*cr(5)+aa(10)*ci(5))
  if(iv.eq.4) ekk=(bb(8)+8d0*(bb(9)*cr(2)+aa(9)*ci(2)))+36d0*(bb(10)*cr(3)+aa(10)*ci(3))
  if(iv.eq.5) ekk=bb(10)
  call detune(iv,ekk,ep,beta,dtu,dtup,dfac)
#else
  tiltck=tiltc(k)**2-tilts(k)**2                               !hr03
  tiltsk=(two*tiltc(k))*tilts(k)                               !hr03
  tiltckuk=tiltck*tiltc(k)-tiltsk*tilts(k)
  tiltsk=tiltck*tilts(k)+tiltsk*tiltc(k)
  tiltck=tiltckuk
  tiltckuk=tiltck*tiltc(k)-tiltsk*tilts(k)
  tiltsk4=tiltck*tilts(k)+tiltsk*tiltc(k)
  tiltck4=tiltckuk
  tiltckuk=tiltck4*tiltc(k)-tiltsk4*tilts(k)
  tiltsk=tiltck4*tilts(k)+tiltsk4*tiltc(k)
  tiltck=tiltckuk
  tiltckuk=tiltck*tiltc(k)-tiltsk*tilts(k)
  tiltsk6=tiltck*tilts(k)+tiltsk*tiltc(k)
  tiltck6=tiltckuk
  tiltckuk=tiltck6*tiltc(k)-tiltsk6*tilts(k)
  tiltsk=tiltck6*tilts(k)+tiltsk6*tiltc(k)
  tiltck=tiltckuk
  tiltckuk=tiltck*tiltc(k)-tiltsk*tilts(k)
  tiltsk8=tiltck*tilts(k)+tiltsk*tiltc(k)
  tiltck8=tiltckuk
  tiltckuk=tiltck8*tiltc(k)-tiltsk8*tilts(k)
  tiltsk=tiltck8*tilts(k)+tiltsk8*tiltc(k)
  tiltck=tiltckuk
  tiltckuk=tiltck*tiltc(k)-tiltsk*tilts(k)
  tiltsk10=tiltck*tilts(k)+tiltsk*tiltc(k)
  tiltck10=tiltckuk
  if(iv.eq.2) then
    ekk= ((((((((((((tiltck4* bb(4)                     -     &!hr03
    tiltsk4*(             aa(4)       ))+                     &!hr03
    four *tiltck4*(bb(5) *cr(2)+aa(5) *ci(2)))+               &!hr03
    four *tiltsk4*(bb(5) *ci(2)-aa(5) *cr(2)))+               &!hr03
    10d0 *tiltck4*(bb(6) *cr(3)+aa(6) *ci(3)))+               &!hr03
    10d0 *tiltsk4*(bb(6) *ci(3)-aa(6) *cr(3)))+               &!hr03
    20d0 *tiltck4*(bb(7) *cr(4)+aa(7) *ci(4)))+               &!hr03
    20d0 *tiltsk4*(bb(7) *ci(4)-aa(7) *cr(4)))+               &!hr03
    35d0 *tiltck4*(bb(8) *cr(5)+aa(8) *ci(5)))+               &!hr03
    35d0 *tiltsk4*(bb(8) *ci(5)-aa(8) *cr(5)))+               &!hr03
    56d0 *tiltck4*(bb(9) *cr(6)+aa(9) *ci(6)))+               &!hr03
    84d0 *tiltck4*(bb(10)*cr(7)+aa(10)*ci(7)))+               &!hr03
    84d0 *tiltck4*(bb(10)*cr(7)+aa(10)*ci(7)))+               &!hr03
    84d0 *tiltsk4*(bb(10)*ci(7)-aa(10)*cr(7))
  endif
  if(iv.eq.3) then
    ekk= ((((((((tiltck6* bb(6)                     -         &!hr03
    tiltsk6*(             aa(6)       ))+                     &!hr03
    6d0  *tiltck6*(bb(7) *cr(2)+aa(7) *ci(2)))+               &!hr03
    6d0  *tiltsk6*(bb(7) *ci(2)-aa(7) *cr(2)))+               &!hr03
    21d0 *tiltck6*(bb(8) *cr(3)+aa(8) *ci(3)))+               &!hr03
    21d0 *tiltsk6*(bb(8) *ci(3)-aa(8) *cr(3)))+               &!hr03
    56d0 *tiltck6*(bb(9) *cr(4)+aa(9) *ci(4)))+               &!hr03
    56d0 *tiltsk6*(bb(9) *ci(4)-aa(9) *cr(4)))+               &!hr03
    126d0*tiltck6*(bb(10)*cr(5)+aa(10)*ci(5)))+               &!hr03
    126d0*tiltsk6*(bb(10)*ci(5)-aa(10)*cr(5))
  endif
  if(iv.eq.4) then
    ekk= ((((tiltck8* bb(8)                     -             &!hr03
    tiltsk8*(             aa(8)       ))+                     &!hr03
    8d0  *tiltck8*(bb(9) *cr(2)+aa(9) *ci(2)))+               &!hr03
    8d0  *tiltsk8*(bb(9) *ci(2)-aa(9) *cr(2)))+               &!hr03
    36d0 *tiltck8*(bb(10)*cr(3)+aa(10)*ci(3)))+               &!hr03
    36d0 *tiltsk8*(bb(10)*ci(3)-aa(10)*cr(3))
  endif
  if(iv.eq.5) then
    ekk= tiltck10*bb(10)-tiltsk10*aa(10)
  endif
  call detune(iv,ekk,ep,beta,dtu,dtup,dfac)
#endif
! end include/multl12.f90
