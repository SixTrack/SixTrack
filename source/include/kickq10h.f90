!start kickq10h.f90
#ifndef TILT
qu=(nine*ekk)*crkve                                               !hr02
qv=(nine*ekk)*cikve                                               !hr02
#else
tiltck=tiltc(k)**2-tilts(k)**2                                   !hr08
tiltsk=(two*tiltc(k))*tilts(k)                                   !hr02
qu=(nine*ekk)*(tiltck*crkve+tiltsk*cikve)                         !hr02
qv=(nine*ekk)*(tiltck*cikve-tiltsk*crkve)                         !hr02
#endif

!end kickq10h.f90
