#ifndef TILT
y(1,1)=y(1,1)+ekk
#else
y(1,1)=y(1,1)+ekk*tiltc(k)
y(1,2)=y(1,2)+ekk*tilts(k)
#endif
