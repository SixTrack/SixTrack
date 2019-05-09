! start include/kickvho.f90
! xlv and zlv is set in include/alignva.f90
crkveuk = crkve*xlv - cikve*zlv
cikve   = crkve*zlv + cikve*xlv
crkve   = crkveuk
! end include/kickvho.f90
