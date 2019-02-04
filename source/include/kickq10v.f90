! start include/kickq10v.f90
#ifndef TILT
qu=(nine*ekk)*cikve                                               !hr02
qv=(-nine*ekk)*crkve                                              !hr02
#else
tiltck=tiltc(k)**2-tilts(k)**2                                   !hr08
tiltsk=(two*tiltc(k))*tilts(k)                                   !hr02
qu=(nine*ekk)*(tiltck*cikve-tiltsk*crkve)                         !hr02
qv=(-nine*ekk)*(tiltck*crkve+tiltsk*cikve)                        !hr02
#endif
! end include/kickq10v.f90
