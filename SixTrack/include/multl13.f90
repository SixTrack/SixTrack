l1=l-1
izu=izu+1
aa(l)=ak0(im,l)+zfz(izu)*aka(im,l)
aa(l)=(benkr*aa(l))/r0a                                        !hr03
izu=izu+1
bb(l)=bk0(im,l)+zfz(izu)*bka(im,l)
bb(l)=(benkr*bb(l))/r0a                                        !hr03
r0a=r0a*r0
if(l.gt.2) then
  cxzyrr=cxzyr*cxzr-cxzyi*cxzi
  cxzyi=cxzyr*cxzi+cxzyi*cxzr
  cxzyr=cxzyrr
  cr(l)=cxzyr
  ci(l)=cxzyi
endif
dyy1=(dyy1+bb(l)*cr(l))+aa(l)*ci(l)                            !hr03
dyy2=(dyy2-bb(l)*ci(l))+aa(l)*cr(l)                            !hr03
if(l.gt.1.and.ium.ne.1) then
  qu=qu+real(l1,fPrec)*(bb(l)*cr(l1)+aa(l)*ci(l1))                   !hr03
  qv=qv+real(l1,fPrec)*(bb(l)*ci(l1)-aa(l)*cr(l1))                   !hr03
endif
