#ifdef DEBUG

subroutine warr(vname,value,i,j,k,l)

  use floatPrecision
  use parpro

  implicit none

  character(len=*) vname
  real(kind=fPrec) value

  integer i,j,k,l
  integer ierro
  character(len=mNameLen) myname

  myname=vname
  write(100) myname,value,i,j,k,l
  endfile (100,iostat=ierro)
  backspace (100,iostat=ierro)

  return

end subroutine warr

subroutine dumpbl1(dumpname,n,i)

  use floatPrecision
  use parpro
  use mod_common

  implicit none

  integer n,i
  character(len=*) dumpname
  save

  write(99,*) dumpname,'   Turn ',n,' Element ',i
  write(99,100) 'bl1 ',bl1
  write(99,100) 'bl2 ',bl2
  endfile (99,iostat=ierro)
  backspace (99,iostat=ierro)

100  format (a10,(Z20))

end subroutine dumpbl1

subroutine dumpzfz(dumpname,n,i)

  use floatPrecision
  use parpro
  use mod_common

  implicit none

  integer n,i
  integer j
  character(len=*) dumpname
  character(len=10) mydump,myzfz
  save

  mydump=dumpname
  myzfz='zfz'
  write(101) mydump,n,i
  write(101) myzfz
  do j=1,nzfz
    write(101) zfz(j)
  end do
  endfile (101,iostat=ierro)
  backspace (101,iostat=ierro)

end subroutine dumpzfz

subroutine dumpxy(dumpname,n,i,k)

  use floatPrecision
  use parpro
  use mod_common
  use mod_commont
  use mod_commonmn

  implicit none

  integer n,i,j,k
  character(len=*) dumpname
  save

  write(99,*) dumpname,'   Turn ',n,' Element ',i
  write(99,*) (xv1(j),j=1,k),(xv2(j),j=1,k),(yv1(j),j=1,k),(yv2(j),j=1,k),&
    (sigmv(j),j=1,k),(ejv(j),j=1,k),(ejfv(j),j=1,k),(rvv(j),j=1,k),           &
    (dpsv(j),j=1,k),(dpsv1(j),j=1,k),(oidpsv(j),j=1,k)
  endfile (99,iostat=ierro)
  backspace (99,iostat=ierro)

end subroutine dumpxy

subroutine dumpsynu(dumpname,n,i)

  use floatPrecision
  use parpro
  use mod_common
  use mod_commonmn
  use mod_commons

  implicit none

  integer n,i,j,l,m,k
  character(len=*) dumpname
  save

  write(99,*) dumpname,'   Turn ',n,' Element ',i
  write(99,*) (aek(j),j=1,napxo)
  write(99,*) (afok(j),j=1,napxo)
  write(99,*) (as3(j),j=1,napxo)
  write(99,*) (as4(j),j=1,napxo)
  write(99,*) (as6(j),j=1,napxo)
  write(99,*) (co(j),j=1,napxo)
  write(99,*) (dpd(j),j=1,napxo)
  write(99,*) (dpsq(j),j=1,napxo)
  write(99,*) (fi(j),j=1,napxo)
  write(99,*) (fok(j),j=1,napxo)
  write(99,*) (fok1(j),j=1,napxo)
  write(99,*) (fokqv(j),j=1,napxo)
  write(99,*) (g(j),j=1,napxo)
  write(99,*) (gl(j),j=1,napxo)
  write(99,*) (hc(j),j=1,napxo)
  write(99,*) (hi(j),j=1,napxo)
  write(99,*) (hi1(j),j=1,napxo)
  write(99,*) (hm(j),j=1,napxo)
  write(99,*) (hp(j),j=1,napxo)
  write(99,*) (hs(j),j=1,napxo)
  write(99,*) (rho(j),j=1,napxo)
  write(99,*) (rhoc(j),j=1,napxo)
  write(99,*) (rhoi(j),j=1,napxo)
  write(99,*) (si(j),j=1,napxo)
  write(99,*) (siq(j),j=1,napxo)
  write(99,*) (sm1(j),j=1,napxo)
  write(99,*) (sm12(j),j=1,napxo)
  write(99,*) (sm2(j),j=1,napxo)
  write(99,*) (sm23(j),j=1,napxo)
  write(99,*) (sm3(j),j=1,napxo)
  write(99,*) (wf(j),j=1,napxo)
  write(99,*) (wfa(j),j=1,napxo)
  write(99,*) (wfhi(j),j=1,napxo)
  write(99,*) ((((al(k,m,j,l),l=1,il),j=1,napxo),m=1,2),k=1,6)
  write(99,*) ((((as(k,m,j,l),l=1,il),j=1,napxo),m=1,2),k=1,6)

  endfile (99,iostat=ierro)
  backspace (99,iostat=ierro)

end subroutine dumpsynu

subroutine dump(dumpname,n,i)

  use floatPrecision
  use parpro
  use mod_common
  use mod_commonmn
  use mod_commons
  use mod_commont
  use mod_commond
  use checkpoint_restart

  implicit none

  integer n,i
  character(len=*) dumpname
  save

  write(99,*) dumpname,'   Turn ',n,' Element ',i

  ! my cr variables
  write(99,*) 'time0 ',time0
  write(99,*) 'time1 ',time1
  write(99,*) 'sixrecs ',sixrecs
  write(99,*) 'binrec ',binrec
  write(99,*) 'binrecs ',binrecs
  write(99,*) 'numlcr ',numlcr
  write(99,*) 'rerun ',rerun
  write(99,*) 'restart ',restart
  write(99,*) 'checkp ',checkp
  write(99,*) 'fort95 ',fort95
  write(99,*) 'fort96 ',fort96
  write(99,*) 'arecord ',arecord
  write(99,*) 'stxt ',stxt
  write(99,*) 'runtim ',runtim

  ! mycrio variables
  write(99,*) 'crnumlcr',crnumlcr
  write(99,*) 'crnuml',crnuml
  write(99,*) 'crsixrecs',crsixrecs
  write(99,*) 'crbinrec',crbinrec
  write(99,*) 'crbinrecs',crbinrecs
  write(99,*) 'crsythck',crsythck
  write(99,*) 'crtime3',crtime3
  write(99,*) 'crnapxo',crnapxo
  write(99,*) 'crnapx',crnapx
  write(99,*) 'cre0',cre0
  write(99,*) 'crnumxv(npart)',crnumxv
  write(99,*) 'crnnumxv(npart)',crnnumxv
  write(99,*) 'crpartID(npart)',crpartID
  write(99,*) 'crparentID(npart)',crparentID
  write(99,*) 'crpstop(npart)',crpstop
  write(99,*) 'crxv',crxv
  write(99,*) 'cryv',cryv
  write(99,*) 'crsigmv',crsigmv
  write(99,*) 'crdpsv',crdpsv
  write(99,*) 'crdpsv1',crdpsv1
  write(99,*) 'crejv',crejv
  write(99,*) 'crejfv',crejfv

  ! some tracking stuff
  write(99,*) 'nwri',nwri
  write(99,*) 'ktrack',ktrack
  write(99,*) 'strack',strack
  write(99,*) 'strackc',strackc
  write(99,*) 'stracks',stracks
  write(99,*) 'dpsv1',dpsv1

  write(99,*) 'ierro ',ierro
  write(99,*) 'erbez ',erbez
  write(99,*) 'pi ',pi
  write(99,*) 'pi2 ',pi2
  write(99,*) 'pisqrt ',pisqrt
  write(99,*) 'rad ',rad
  write(99,*) 'il ',il
  write(99,*) 'mper ',mper
  write(99,*) 'mblo ',mblo
  write(99,*) 'mbloz ',mbloz
  write(99,*) 'msym ',msym
  write(99,*) 'kanf ',kanf
  write(99,*) 'iu ',iu
  write(99,*) 'ic ',ic
  write(99,*) 'ed ',ed
  write(99,*) 'el ',el
  write(99,*) 'ek ',ek
  write(99,*) 'sm ',sm
  write(99,*) 'kz ',kz
  write(99,*) 'kp ',kp
  write(99,*) 'xpl ',xpl
  write(99,*) 'xrms ',xrms
  write(99,*) 'zpl ',zpl
  write(99,*) 'zrms ',zrms
  write(99,*) 'mel ',mel
  write(99,*) 'mtyp ',mtyp
  write(99,*) 'mstr ',mstr
  write(99,*) 'a ',a
  write(99,*) 'bl1 ',bl1
  write(99,*) 'bl2 ',bl2
  write(99,*) 'rvf ',rvf
  write(99,*) 'idfor ',idfor
  write(99,*) 'napx ',napx
  write(99,*) 'napxo ',napxo
  write(99,*) 'numlr ',numlr
  write(99,*) 'nde ',nde
  write(99,*) 'nwr ',nwr
  write(99,*) 'ird ',ird
  write(99,*) 'imc ',imc
  write(99,*) 'irew ',irew
  write(99,*) 'ntwin ',ntwin
  write(99,*) 'iclo6 ',iclo6
  write(99,*) 'iclo6r ',iclo6r
  write(99,*) 'iver ',iver
  write(99,*) 'qs ',qs
  write(99,*) 'e0 ',e0
  write(99,*) 'pma ',pma
  write(99,*) 'ej ',ej
  write(99,*) 'ejf ',ejf
  write(99,*) 'phas0 ',phas0
  write(99,*) 'phas ',phas
  write(99,*) 'hsy ',hsy
  write(99,*) 'crad ',crad
  write(99,*) 'hsyc ',hsyc
  write(99,*) 'phasc ',phasc
  write(99,*) 'dppoff ',dppoff
  write(99,*) 'sigmoff ',sigmoff
  write(99,*) 'tlen ',tlen
  write(99,*) 'iicav ',iicav
  write(99,*) 'itionc ',itionc
  write(99,*) 'ition ',ition
  write(99,*) 'idp ',idp
  write(99,*) 'ncy ',ncy
  write(99,*) 'ixcav ',ixcav
  write(99,*) 'dpscor ',dpscor
  write(99,*) 'sigcor ',sigcor
  write(99,*) 'icode ',icode
  write(99,*) 'idam ',idam
  write(99,*) 'its6d ',its6d
  write(99,*) 'bk0 ',bk0
  write(99,*) 'ak0 ',ak0
  write(99,*) 'bka ',bka
  write(99,*) 'aka ',aka
  write(99,*) 'benki ',benki
  write(99,*) 'benkc ',benkc
  write(99,*) 'r00 ',r00
  write(99,*) 'irm ',irm
  write(99,*) 'nmu ',nmu
  write(99,*) 'zfz ',zfz
  write(99,*) 'iorg ',iorg
  write(99,*) 'mzu ',mzu
  write(99,*) 'bezr ',bezr
  write(99,*) 'izu0 ',izu0
  write(99,*) 'mcut ',mcut
! write(99,*) 'exterr ',exterr
! write(99,*) 'extalign ',extalign
  write(99,*) 'tiltc ',tiltc
  write(99,*) 'tilts ',tilts
  write(99,*) 'mout2 ',mout2
  write(99,*) 'icext ',icext
  write(99,*) 'icextal ',icextal
  write(99,*) 'aper ',aper
  write(99,*) 'di0 ',di0
  write(99,*) 'dip0 ',dip0
  write(99,*) 'ta ',ta
  write(99,*) 'dma ',dma
  write(99,*) 'dmap ',dmap
  write(99,*) 'dkq ',dkq
  write(99,*) 'dqq ',dqq
  write(99,*) 'de0 ',de0
  write(99,*) 'ded ',ded
  write(99,*) 'dsi ',dsi
  write(99,*) 'dech ',dech
  write(99,*) 'dsm0 ',dsm0
  write(99,*) 'itco ',itco
  write(99,*) 'itcro ',itcro
  write(99,*) 'itqv ',itqv
  write(99,*) 'qw0 ',qw0
  write(99,*) 'iq ',iq
  write(99,*) 'iqmod ',iqmod
  write(99,*) 'kpa ',kpa
  write(99,*) 'iqmod6 ',iqmod6
  write(99,*) 'bez ',bez
  write(99,*) 'elbe ',elbe
  write(99,*) 'bezb ',bezb
  write(99,*) 'ilin ',ilin
  write(99,*) 'nt ',nt
  write(99,*) 'iprint ',iprint
  write(99,*) 'ntco ',ntco
  write(99,*) 'eui ',eui
  write(99,*) 'euii ',euii
  write(99,*) 'nlin ',nlin
  write(99,*) 'bezl ',bezl
  write(99,*) 'betam ',betam
  write(99,*) 'pam ',pam
  write(99,*) 'betac ',betac
  write(99,*) 'pac ',pac
  write(99,*) 'bclorb ',bclorb
  write(99,*) 'nhmoni ',nhmoni
  write(99,*) 'nhcorr ',nhcorr
  write(99,*) 'nvmoni ',nvmoni
  write(99,*) 'nvcorr ',nvcorr
  write(99,*) 'ncororb ',ncororb
  ! write(99,*) 'apx ',apx
  ! write(99,*) 'apz ',apz
  write(99,*) 'sigma0 ',sigma0
  write(99,*) 'iclo ',iclo
  write(99,*) 'ncorru ',ncorru
  write(99,*) 'ncorrep ',ncorrep
  write(99,*) 'icomb0 ',icomb0
  write(99,*) 'icomb ',icomb
  write(99,*) 'ratio ',ratio
  write(99,*) 'ratioe ',ratioe
  write(99,*) 'iratioe ',iratioe
  write(99,*) 'icoe ',icoe
  write(99,*) 'ise ',ise
  write(99,*) 'mesa ',mesa
  write(99,*) 'mp ',mp
  write(99,*) 'm21 ',m21
  write(99,*) 'm22 ',m22
  write(99,*) 'm23 ',m23
  write(99,*) 'ise1 ',ise1
  write(99,*) 'ise2 ',ise2
  write(99,*) 'ise3 ',ise3
  write(99,*) 'isea ',isea
  write(99,*) 'qxt ',qxt
  write(99,*) 'qzt ',qzt
  write(99,*) 'tam1 ',tam1
  write(99,*) 'tam2 ',tam2
  write(99,*) 'isub ',isub
  write(99,*) 'nta ',nta
  write(99,*) 'nte ',nte
  write(99,*) 'ipt ',ipt
  write(99,*) 'totl ',totl
  write(99,*) 'rtc ',rtc
  write(99,*) 'rts ',rts
  write(99,*) 'ire ',ire
  write(99,*) 'ipr ',ipr
  write(99,*) 'irmod2 ',irmod2
  write(99,*) 'dtr ',dtr
  write(99,*) 'nre ',nre
  write(99,*) 'nur ',nur
  write(99,*) 'nch ',nch
  write(99,*) 'nqc ',nqc
  write(99,*) 'npp ',npp
  write(99,*) 'nrr ',nrr
  write(99,*) 'nu ',nu
  write(99,*) 'dphix ',dphix
  write(99,*) 'dphiz ',dphiz
  write(99,*) 'qx0 ',qx0
  write(99,*) 'qz0 ',qz0
  write(99,*) 'dres ',dres
  write(99,*) 'dfft ',dfft
  write(99,*) 'cma1 ',cma1
  write(99,*) 'cma2 ',cma2
  write(99,*) 'nstart ',nstart
  write(99,*) 'nstop ',nstop
  write(99,*) 'iskip ',iskip
  write(99,*) 'iconv ',iconv
  write(99,*) 'imad ',imad
  write(99,*) 'ipos ',ipos
  write(99,*) 'iav ',iav
  write(99,*) 'iwg ',iwg
  write(99,*) 'ivox ',ivox
  write(99,*) 'ivoz ',ivoz
  write(99,*) 'ires ',ires
  write(99,*) 'ifh ',ifh
  write(99,*) 'toptit ',toptit
  write(99,*) 'kwtype ',kwtype
  write(99,*) 'itf ',itf
  write(99,*) 'icr ',icr
  write(99,*) 'idis ',idis
  write(99,*) 'icow ',icow
  write(99,*) 'istw ',istw
  write(99,*) 'iffw ',iffw
  write(99,*) 'nprint ',nprint
  write(99,*) 'ndafi ',ndafi
  write(99,*) 'qwsk ',qwsk
  write(99,*) 'betx ',betx
  write(99,*) 'betz ',betz
  write(99,*) 'alfx ',alfx
  write(99,*) 'alfz ',alfz
  write(99,*) 'iskew ',iskew
  write(99,*) 'nskew ',nskew
  write(99,*) 'hmal ',hmal
  write(99,*) 'sixtit ',sixtit
  write(99,*) 'commen ',commen
  write(99,*) 'ithick ',ithick
  write(99,*) 'clo6 ',clo6
  write(99,*) 'clop6 ',clop6
  write(99,*) 'dki ',dki
  write(99,*) 'sigman ',sigman
  write(99,*) 'sigman2 ',sigman2
  write(99,*) 'sigmanq ',sigmanq
  write(99,*) 'clobeam ',clobeam
  write(99,*) 'beamoff ',beamoff
  write(99,*) 'parbe ',parbe
  write(99,*) 'track6d ',track6d
  write(99,*) 'ptnfac ',ptnfac
  write(99,*) 'sigz ',sigz
  write(99,*) 'sige ',sige
  write(99,*) 'partnum ',partnum
  write(99,*) 'parbe14 ',parbe14
  write(99,*) 'emitx ',emitx
  write(99,*) 'emity ',emity
  write(99,*) 'emitz ',emitz
  write(99,*) 'gammar ',gammar
  write(99,*) 'nbeam ',nbeam
  write(99,*) 'ibbc ',ibbc
  write(99,*) 'ibeco ',ibeco
  write(99,*) 'ibtyp ',ibtyp
  write(99,*) 'lhc ',lhc
  write(99,*) 'cotr ',cotr
  write(99,*) 'rrtr ',rrtr
  write(99,*) 'imtr ',imtr
  write(99,*) 'bbcu ',bbcu
  write(99,*) 'ibb6d ',ibb6d
  write(99,*) 'imbb ',imbb
  write(99,*) 'as ',as
  write(99,*) 'al ',al
  write(99,*) 'sigm ',sigm
  write(99,*) 'dps ',dps
  write(99,*) 'idz ',idz
  write(99,*) 'dp1 ',dp1
  write(99,*) 'itra ',itra
  write(99,*) 'x ',x
  write(99,*) 'y ',y
  write(99,*) 'bet0 ',bet0
  write(99,*) 'alf0 ',alf0
  write(99,*) 'clo ',clo
  write(99,*) 'clop ',clop
  write(99,*) 'cro ',cro
  write(99,*) 'is ',is
  write(99,*) 'ichrom ',ichrom
  write(99,*) 'nnumxv ',nnumxv
  write(99,*) 'xsi ',xsi
  write(99,*) 'zsi ',zsi
  write(99,*) 'smi ',smi
  write(99,*) 'ampt ',ampt
  write(99,*) 'tlim ',tlim
  write(99,*) 'tasm ',tasm
  write(99,*) 'preda ',preda
  write(99,*) 'idial ',idial
  write(99,*) 'nord ',nord
  write(99,*) 'nvar ',nvar
  write(99,*) 'nvar2 ',nvar2
  write(99,*) 'nsix ',nsix
  write(99,*) 'ncor ',ncor
  write(99,*) 'ipar ',ipar
  write(99,*) 'nordf ',nordf
  write(99,*) 'nvarf ',nvarf
  write(99,*) 'nord1 ',nord1
  write(99,*) 'ndimf ',ndimf
  write(99,*) 'idptr ',idptr
  write(99,*) 'inorm ',inorm
  write(99,*) 'imod1 ',imod1
  write(99,*) 'imod2 ',imod2
  write(99,*) 'ekv ',ekv
  write(99,*) 'fokqv ',fokqv
  write(99,*) 'aaiv ',aaiv
  write(99,*) 'bbiv ',bbiv
  write(99,*) 'smiv ',smiv
  write(99,*) 'zsiv ',zsiv
  write(99,*) 'xsiv ',xsiv
  write(99,*) 'xsv ',xsv
  write(99,*) 'zsv ',zsv
  write(99,*) 'qw ',qw
  write(99,*) 'qwc ',qwc
  write(99,*) 'clo0 ',clo0
  write(99,*) 'clop0 ',clop0
  write(99,*) 'eps ',eps
  write(99,*) 'epsa ',epsa
  write(99,*) 'ekk ',ekk
  write(99,*) 'cr ',cr
  write(99,*) 'ci ',ci
  write(99,*) 'xv1 ',xv1
  write(99,*) 'yv1 ',yv1
  write(99,*) 'xv2 ',xv2
  write(99,*) 'yv2 ',yv2
  write(99,*) 'dam ',dam
  write(99,*) 'ekkv ',ekkv
  write(99,*) 'sigmv ',sigmv
  write(99,*) 'dpsv ',dpsv
  write(99,*) 'dp0v ',dp0v
  write(99,*) 'sigmv6 ',sigmv6
  write(99,*) 'dpsv6 ',dpsv6
  write(99,*) 'ejv ',ejv
  write(99,*) 'ejfv ',ejfv
  write(99,*) 'xlv ',xlv
  write(99,*) 'zlv ',zlv
  write(99,*) 'pstop ',pstop
  write(99,*) 'rvv ',rvv
  write(99,*) 'ejf0v ',ejf0v
  write(99,*) 'numxv ',numxv
  write(99,*) 'nms ',nms
  write(99,*) 'partID ',partID
  write(99,*) 'parentID ',parentID
  write(99,*) 'dpd ',dpd
  write(99,*) 'dpsq ',dpsq
  write(99,*) 'fok ',fok
  write(99,*) 'rho ',rho
  write(99,*) 'fok1 ',fok1
  write(99,*) 'si ',si
  write(99,*) 'co ',co
  write(99,*) 'g ',g
  write(99,*) 'gl ',gl
  write(99,*) 'sm1 ',sm1
  write(99,*) 'sm2 ',sm2
  write(99,*) 'sm3 ',sm3
  write(99,*) 'sm12 ',sm12
  write(99,*) 'as3 ',as3
  write(99,*) 'as4 ',as4
  write(99,*) 'as6 ',as6
  write(99,*) 'sm23 ',sm23
  write(99,*) 'rhoc ',rhoc
  write(99,*) 'siq ',siq
  write(99,*) 'aek ',aek
  write(99,*) 'afok ',afok
  write(99,*) 'hp ',hp
  write(99,*) 'hm ',hm
  write(99,*) 'hc ',hc
  write(99,*) 'hs ',hs
  write(99,*) 'wf ',wf
  write(99,*) 'wfa ',wfa
  write(99,*) 'wfhi ',wfhi
  write(99,*) 'rhoi ',rhoi
  write(99,*) 'hi ',hi
  write(99,*) 'fi ',fi
  write(99,*) 'hi1 ',hi1
  write(99,*) 'xvl ',xvl
  write(99,*) 'yvl ',yvl
  write(99,*) 'ejvl ',ejvl
  write(99,*) 'dpsvl ',dpsvl
  write(99,*) 'oidpsv ',oidpsv
  write(99,*) 'sigmvl ',sigmvl
  write(99,*) 'iv ',iv
  write(99,*) 'aperv ',aperv
  write(99,*) 'ixv ',ixv
  write(99,*) 'ampv ',ampv
  write(99,*) 'clo6v ',clo6v
  write(99,*) 'clop6v ',clop6v
  write(99,*) 'hv ',hv
  write(99,*) 'bl1v ',bl1v
  write(99,*) 'tas ',tas
  write(99,*) 'qwcs ',qwcs
  write(99,*) 'di0xs ',di0xs
  write(99,*) 'di0zs ',di0zs
  write(99,*) 'dip0xs ',dip0xs
  write(99,*) 'dip0zs ',dip0zs
  write(99,*) 'xau ',xau
  write(99,*) 'cloau ',cloau
  write(99,*) 'di0au ',di0au
  write(99,*) 'tau ',tau
  write(99,*) 'tasau ',tasau
  write(99,*) 'wx ',wx
  write(99,*) 'x1 ',x1
  write(99,*) 'x2 ',x2
  write(99,*) 'fake ',fake

  write(99,*) 'e0f ',e0f
  write(99,*) 'numx ',numx
  write(99,*) 'cotr ',cotr
  write(99,*) 'rrtr ',rrtr
  write(99,*) 'imtr ',imtr

  ! these other values???
  write(99,*) 'numl ',numl
  write(99,*) 'niu ',niu
  write(99,*) 'amp0 ',amp0
  write(99,*) 'amp ',amp
  write(99,*) 'damp ',damp
  write(99,*) 'chi0 ',chi0
  write(99,*) 'chid ',chid
  write(99,*) 'rat ',rat
  write(99,*) 'exz ',exz
  write(99,*) 'time0 ',time0
  write(99,*) 'time1 ',time1

  endfile (99,iostat=ierro)
  backspace (99,iostat=ierro)

end subroutine dump

subroutine dumpbin(dumpname,n,i)

  use floatPrecision
  use parpro
  use mod_common
  use mod_commonmn
  use mod_commons
  use mod_commont
  use mod_commond
  use checkpoint_restart

  implicit none

  integer n,i
  character(len=*) dumpname
  character(len=10) mydump
  save

  mydump=dumpname

  write(99) mydump
  write(99) n
  write(99) i

  ! my cr variables
  write(99) time0
  write(99) time1
  write(99) sixrecs
  write(99) binrec
  write(99) binrecs
  write(99) numlcr
  write(99) rerun
  write(99) restart
  write(99) checkp
  write(99) fort95
  write(99) fort96
  write(99) arecord
  write(99) stxt
  write(99) runtim

  ! mycrio variables
  write(99) crnumlcr
  write(99) crnuml
  write(99) crsixrecs
  write(99) crbinrec
  write(99) crbinrecs
  write(99) crsythck
  write(99) crtime3
  write(99) crnapxo
  write(99) crnapx
  write(99) cre0
  write(99) crnumxv
  write(99) crnnumxv
  write(99) crpartID
  write(99) crparentID
  write(99) crpstop
  write(99) crxv
  write(99) cryv
  write(99) crsigmv
  write(99) crdpsv
  write(99) crdpsv1
  write(99) crejv
  write(99) crejfv

  ! some tracking stuff
  write(99) nwri
  write(99) ktrack
  write(99) strack
  write(99) strackc
  write(99) stracks
  write(99) dpsv1

  write(99) ierro
  write(99) erbez
  write(99) pi
  write(99) pi2
  write(99) pisqrt
  write(99) rad
  write(99) il
  write(99) mper
  write(99) mblo
  write(99) mbloz
  write(99) msym
  write(99) kanf
  write(99) iu
  write(99) ic
  write(99) ed
  write(99) el
  write(99) ek
  write(99) sm
  write(99) kz
  write(99) kp
  write(99) xpl
  write(99) xrms
  write(99) zpl
  write(99) zrms
  write(99) mel
  write(99) mtyp
  write(99) mstr
  write(99) a
  write(99) bl1
  write(99) bl2
  write(99) rvf
  write(99) idfor
  write(99) napx
  write(99) napxo
  write(99) numlr
  write(99) nde
  write(99) nwr
  write(99) ird
  write(99) imc
  write(99) irew
  write(99) ntwin
  write(99) iclo6
  write(99) iclo6r
  write(99) iver
  write(99) qs
  write(99) e0
  write(99) pma
  write(99) ej
  write(99) ejf
  write(99) phas0
  write(99) phas
  write(99) hsy
  write(99) crad
  write(99) hsyc
  write(99) phasc
  write(99) dppoff
  write(99) sigmoff
  write(99) tlen
  write(99) iicav
  write(99) itionc
  write(99) ition
  write(99) idp
  write(99) ncy
  write(99) ixcav
  write(99) dpscor
  write(99) sigcor
  write(99) icode
  write(99) idam
  write(99) its6d
  write(99) bk0
  write(99) ak0
  write(99) bka
  write(99) aka
  write(99) benki
  write(99) benkc
  write(99) r00
  write(99) irm
  write(99) nmu
  write(99) zfz
  write(99) iorg
  write(99) mzu
  write(99) bezr
  write(99) izu0
  write(99) mcut
! write(99) exterr
! write(99) extalign
  write(99) tiltc
  write(99) tilts
  write(99) mout2
  write(99) icext
  write(99) icextal
  write(99) aper
  write(99) di0
  write(99) dip0
  write(99) ta
  write(99) dma
  write(99) dmap
  write(99) dkq
  write(99) dqq
  write(99) de0
  write(99) ded
  write(99) dsi
  write(99) dech
  write(99) dsm0
  write(99) itco
  write(99) itcro
  write(99) itqv
  write(99) qw0
  write(99) iq
  write(99) iqmod
  write(99) kpa
  write(99) iqmod6
  write(99) bez
  write(99) elbe
  write(99) bezb
  write(99) ilin
  write(99) nt
  write(99) iprint
  write(99) ntco
  write(99) eui
  write(99) euii
  write(99) nlin
  write(99) bezl
  write(99) betam
  write(99) pam
  write(99) betac
  write(99) pac
  write(99) bclorb
  write(99) nhmoni
  write(99) nhcorr
  write(99) nvmoni
  write(99) nvcorr
  write(99) ncororb
  ! write(99) apx
  ! write(99) apz
  write(99) sigma0
  write(99) iclo
  write(99) ncorru
  write(99) ncorrep
  write(99) icomb0
  write(99) icomb
  write(99) ratio
  write(99) ratioe
  write(99) iratioe
  write(99) icoe
  write(99) ise
  write(99) mesa
  write(99) mp
  write(99) m21
  write(99) m22
  write(99) m23
  write(99) ise1
  write(99) ise2
  write(99) ise3
  write(99) isea
  write(99) qxt
  write(99) qzt
  write(99) tam1
  write(99) tam2
  write(99) isub
  write(99) nta
  write(99) nte
  write(99) ipt
  write(99) totl
  write(99) rtc
  write(99) rts
  write(99) ire
  write(99) ipr
  write(99) irmod2
  write(99) dtr
  write(99) nre
  write(99) nur
  write(99) nch
  write(99) nqc
  write(99) npp
  write(99) nrr
  write(99) nu
  write(99) dphix
  write(99) dphiz
  write(99) qx0
  write(99) qz0
  write(99) dres
  write(99) dfft
  write(99) cma1
  write(99) cma2
  write(99) nstart
  write(99) nstop
  write(99) iskip
  write(99) iconv
  write(99) imad
  write(99) ipos
  write(99) iav
  write(99) iwg
  write(99) ivox
  write(99) ivoz
  write(99) ires
  write(99) ifh
  write(99) toptit
  write(99) kwtype
  write(99) itf
  write(99) icr
  write(99) idis
  write(99) icow
  write(99) istw
  write(99) iffw
  write(99) nprint
  write(99) ndafi
  write(99) qwsk
  write(99) betx
  write(99) betz
  write(99) alfx
  write(99) alfz
  write(99) iskew
  write(99) nskew
  write(99) hmal
  write(99) sixtit
  write(99) commen
  write(99) ithick
  write(99) clo6
  write(99) clop6
  write(99) dki
  write(99) sigman
  write(99) sigman2
  write(99) sigmanq
  write(99) clobeam
  write(99) beamoff
  write(99) parbe
  write(99) track6d
  write(99) ptnfac
  write(99) sigz
  write(99) sige
  write(99) partnum
  write(99) parbe14
  write(99) emitx
  write(99) emity
  write(99) emitz
  write(99) gammar
  write(99) nbeam
  write(99) ibbc
  write(99) ibeco
  write(99) ibtyp
  write(99) lhc
  write(99) cotr
  write(99) rrtr
  write(99) imtr
  write(99) bbcu
  write(99) ibb6d
  write(99) imbb
  write(99) as
  write(99) al
  write(99) sigm
  write(99) dps
  write(99) idz
  write(99) dp1
  write(99) itra
  write(99) x
  write(99) y
  write(99) bet0
  write(99) alf0
  write(99) clo
  write(99) clop
  write(99) cro
  write(99) is
  write(99) ichrom
  write(99) nnumxv
  write(99) xsi
  write(99) zsi
  write(99) smi
  write(99) ampt
  write(99) tlim
  write(99) tasm
  write(99) preda
  write(99) idial
  write(99) nord
  write(99) nvar
  write(99) nvar2
  write(99) nsix
  write(99) ncor
  write(99) ipar
  write(99) nordf
  write(99) nvarf
  write(99) nord1
  write(99) ndimf
  write(99) idptr
  write(99) inorm
  write(99) imod1
  write(99) imod2
  write(99) ekv
  write(99) fokqv
  write(99) aaiv
  write(99) bbiv
  write(99) smiv
  write(99) zsiv
  write(99) xsiv
  write(99) xsv
  write(99) zsv
  write(99) qw
  write(99) qwc
  write(99) clo0
  write(99) clop0
  write(99) eps
  write(99) epsa
  write(99) ekk
  write(99) cr
  write(99) ci
  write(99) xv1
  write(99) yv1
  write(99) xv2
  write(99) yv2
  write(99) dam
  write(99) ekkv
  write(99) sigmv
  write(99) dpsv
  write(99) dp0v
  write(99) sigmv6
  write(99) dpsv6
  write(99) ejv
  write(99) ejfv
  write(99) xlv
  write(99) zlv
  write(99) pstop
  write(99) rvv
  write(99) ejf0v
  write(99) numxv
  write(99) nms
  write(99) parentID
  write(99) dpd
  write(99) dpsq
  write(99) fok
  write(99) rho
  write(99) fok1
  write(99) si
  write(99) co
  write(99) g
  write(99) gl
  write(99) sm1
  write(99) sm2
  write(99) sm3
  write(99) sm12
  write(99) as3
  write(99) as4
  write(99) as6
  write(99) sm23
  write(99) rhoc
  write(99) siq
  write(99) aek
  write(99) afok
  write(99) hp
  write(99) hm
  write(99) hc
  write(99) hs
  write(99) wf
  write(99) wfa
  write(99) wfhi
  write(99) rhoi
  write(99) hi
  write(99) fi
  write(99) hi1
  write(99) xvl
  write(99) yvl
  write(99) ejvl
  write(99) dpsvl
  write(99) oidpsv
  write(99) sigmvl
  write(99) iv
  write(99) aperv
  write(99) ixv
  write(99) ampv
  write(99) clo6v
  write(99) clop6v
  write(99) hv
  write(99) bl1v
  write(99) tas
  write(99) qwcs
  write(99) di0xs
  write(99) di0zs
  write(99) dip0xs
  write(99) dip0zs
  write(99) xau
  write(99) cloau
  write(99) di0au
  write(99) tau
  write(99) tasau
  write(99) wx
  write(99) x1
  write(99) x2
  write(99) fake
  write(99) e0f
  write(99) numx
  write(99) cotr
  write(99) rrtr
  write(99) imtr

  ! these other values???
  write(99) numl
  write(99) niu
  write(99) amp0
  write(99) amp
  write(99) damp
  write(99) chi0
  write(99) chid
  write(99) rat
  write(99) exz
  write(99) time0
  write(99) time1

  endfile (99,iostat=ierro)
  backspace (99,iostat=ierro)

end subroutine dumpbin

subroutine dumphex(dumpname,n,i)

  use floatPrecision
  use parpro
  use mod_common
  use mod_commonmn
  use mod_commons
  use mod_commont
  use mod_commond
  use checkpoint_restart

  implicit none

  integer n,i
  character(len=*) dumpname
  save

  write(99,*) dumpname,'   Turn ',n,' Element ',i

  ! my cr variables
  write(99,100) 'time0 ',time0
  write(99,100) 'time1 ',time1
  write(99,100) 'sixrecs ',sixrecs
  write(99,100) 'binrec ',binrec
  write(99,100) 'binrecs ',binrecs
  write(99,100) 'numlcr ',numlcr
  write(99,100) 'rerun ',rerun
  write(99,100) 'restart ',restart
  write(99,100) 'checkp ',checkp
  write(99,100) 'fort95 ',fort95
  write(99,100) 'fort96 ',fort96
  write(99,100) 'arecord ',arecord
  write(99,100) 'stxt ',stxt
  write(99,100) 'runtim ',runtim

  ! mycrio variables
  write(99,100) 'crnumlcr',crnumlcr
  write(99,100) 'crnuml',crnuml
  write(99,100) 'crsixrecs',crsixrecs
  write(99,100) 'crbinrec',crbinrec
  write(99,100) 'crbinrecs',crbinrecs
  write(99,100) 'crsythck',crsythck
  write(99,100) 'crtime3',crtime3
  write(99,100) 'crnapxo',crnapxo
  write(99,100) 'crnapx',crnapx
  write(99,100) 'cre0',cre0
  write(99,100) 'crnumxv(npart)',crnumxv
  write(99,100) 'crnnumxv(npart)',crnnumxv
  write(99,100) 'crpartID(npart)',crpartID
  write(99,100) 'crparentID(npart)',crparentID
  write(99,100) 'crpstop(npart)',crpstop
  write(99,100) 'crxv',crxv
  write(99,100) 'cryv',cryv
  write(99,100) 'crsigmv',crsigmv
  write(99,100) 'crdpsv',crdpsv
  write(99,100) 'crdpsv1',crdpsv1
  write(99,100) 'crejv',crejv
  write(99,100) 'crejfv',crejfv

  ! some tracking stuff
  write(99,100) 'nwri',nwri
  write(99,100) 'ktrack',ktrack
  write(99,100) 'strack',strack
  write(99,100) 'strackc',strackc
  write(99,100) 'stracks',stracks
  write(99,100) 'dpsv1',dpsv1

  write(99,100) 'ierro ',ierro
  write(99,100) 'erbez ',erbez
  write(99,100) 'pi ',pi
  write(99,100) 'pi2 ',pi2
  write(99,100) 'pisqrt ',pisqrt
  write(99,100) 'rad ',rad
  write(99,100) 'il ',il
  write(99,100) 'mper ',mper
  write(99,100) 'mblo ',mblo
  write(99,100) 'mbloz ',mbloz
  write(99,100) 'msym ',msym
  write(99,100) 'kanf ',kanf
  write(99,100) 'iu ',iu
  write(99,100) 'ic ',ic
  write(99,100) 'ed ',ed
  write(99,100) 'el ',el
  write(99,100) 'ek ',ek
  write(99,100) 'sm ',sm
  write(99,100) 'kz ',kz
  write(99,100) 'kp ',kp
  write(99,100) 'xpl ',xpl
  write(99,100) 'xrms ',xrms
  write(99,100) 'zpl ',zpl
  write(99,100) 'zrms ',zrms
  write(99,100) 'mel ',mel
  write(99,100) 'mtyp ',mtyp
  write(99,100) 'mstr ',mstr
  write(99,100) 'a ',a
  write(99,100) 'bl1 ',bl1
  write(99,100) 'bl2 ',bl2
  write(99,100) 'rvf ',rvf
  write(99,100) 'idfor ',idfor
  write(99,100) 'napx ',napx
  write(99,100) 'napxo ',napxo
  write(99,100) 'numlr ',numlr
  write(99,100) 'nde ',nde
  write(99,100) 'nwr ',nwr
  write(99,100) 'ird ',ird
  write(99,100) 'imc ',imc
  write(99,100) 'irew ',irew
  write(99,100) 'ntwin ',ntwin
  write(99,100) 'iclo6 ',iclo6
  write(99,100) 'iclo6r ',iclo6r
  write(99,100) 'iver ',iver
  write(99,100) 'qs ',qs
  write(99,100) 'e0 ',e0
  write(99,100) 'pma ',pma
  write(99,100) 'ej ',ej
  write(99,100) 'ejf ',ejf
  write(99,100) 'phas0 ',phas0
  write(99,100) 'phas ',phas
  write(99,100) 'hsy ',hsy
  write(99,100) 'crad ',crad
  write(99,100) 'hsyc ',hsyc
  write(99,100) 'phasc ',phasc
  write(99,100) 'dppoff ',dppoff
  write(99,100) 'sigmoff ',sigmoff
  write(99,100) 'tlen ',tlen
  write(99,100) 'iicav ',iicav
  write(99,100) 'itionc ',itionc
  write(99,100) 'ition ',ition
  write(99,100) 'idp ',idp
  write(99,100) 'ncy ',ncy
  write(99,100) 'ixcav ',ixcav
  write(99,100) 'dpscor ',dpscor
  write(99,100) 'sigcor ',sigcor
  write(99,100) 'icode ',icode
  write(99,100) 'idam ',idam
  write(99,100) 'its6d ',its6d
  write(99,100) 'bk0 ',bk0
  write(99,100) 'ak0 ',ak0
  write(99,100) 'bka ',bka
  write(99,100) 'aka ',aka
  write(99,100) 'benki ',benki
  write(99,100) 'benkc ',benkc
  write(99,100) 'r00 ',r00
  write(99,100) 'irm ',irm
  write(99,100) 'nmu ',nmu
  write(99,100) 'zfz ',zfz
  write(99,100) 'iorg ',iorg
  write(99,100) 'mzu ',mzu
  write(99,100) 'bezr ',bezr
  write(99,100) 'izu0 ',izu0
  write(99,100) 'mcut ',mcut
! write(99,100) 'exterr ',exterr
! write(99,100) 'extalign ',extalign
  write(99,100) 'tiltc ',tiltc
  write(99,100) 'tilts ',tilts
  write(99,100) 'mout2 ',mout2
  write(99,100) 'icext ',icext
  write(99,100) 'icextal ',icextal
  write(99,100) 'aper ',aper
  write(99,100) 'di0 ',di0
  write(99,100) 'dip0 ',dip0
  write(99,100) 'ta ',ta
  write(99,100) 'dma ',dma
  write(99,100) 'dmap ',dmap
  write(99,100) 'dkq ',dkq
  write(99,100) 'dqq ',dqq
  write(99,100) 'de0 ',de0
  write(99,100) 'ded ',ded
  write(99,100) 'dsi ',dsi
  write(99,100) 'dech ',dech
  write(99,100) 'dsm0 ',dsm0
  write(99,100) 'itco ',itco
  write(99,100) 'itcro ',itcro
  write(99,100) 'itqv ',itqv
  write(99,100) 'qw0 ',qw0
  write(99,100) 'iq ',iq
  write(99,100) 'iqmod ',iqmod
  write(99,100) 'kpa ',kpa
  write(99,100) 'iqmod6 ',iqmod6
  write(99,100) 'bez ',bez
  write(99,100) 'elbe ',elbe
  write(99,100) 'bezb ',bezb
  write(99,100) 'ilin ',ilin
  write(99,100) 'nt ',nt
  write(99,100) 'iprint ',iprint
  write(99,100) 'ntco ',ntco
  write(99,100) 'eui ',eui
  write(99,100) 'euii ',euii
  write(99,100) 'nlin ',nlin
  write(99,100) 'bezl ',bezl
  write(99,100) 'betam ',betam
  write(99,100) 'pam ',pam
  write(99,100) 'betac ',betac
  write(99,100) 'pac ',pac
  write(99,100) 'bclorb ',bclorb
  write(99,100) 'nhmoni ',nhmoni
  write(99,100) 'nhcorr ',nhcorr
  write(99,100) 'nvmoni ',nvmoni
  write(99,100) 'nvcorr ',nvcorr
  write(99,100) 'ncororb ',ncororb
  ! write(99,100) 'apx ',apx
  ! write(99,100) 'apz ',apz
  write(99,100) 'sigma0 ',sigma0
  write(99,100) 'iclo ',iclo
  write(99,100) 'ncorru ',ncorru
  write(99,100) 'ncorrep ',ncorrep
  write(99,100) 'icomb0 ',icomb0
  write(99,100) 'icomb ',icomb
  write(99,100) 'ratio ',ratio
  write(99,100) 'ratioe ',ratioe
  write(99,100) 'iratioe ',iratioe
  write(99,100) 'icoe ',icoe
  write(99,100) 'ise ',ise
  write(99,100) 'mesa ',mesa
  write(99,100) 'mp ',mp
  write(99,100) 'm21 ',m21
  write(99,100) 'm22 ',m22
  write(99,100) 'm23 ',m23
  write(99,100) 'ise1 ',ise1
  write(99,100) 'ise2 ',ise2
  write(99,100) 'ise3 ',ise3
  write(99,100) 'isea ',isea
  write(99,100) 'qxt ',qxt
  write(99,100) 'qzt ',qzt
  write(99,100) 'tam1 ',tam1
  write(99,100) 'tam2 ',tam2
  write(99,100) 'isub ',isub
  write(99,100) 'nta ',nta
  write(99,100) 'nte ',nte
  write(99,100) 'ipt ',ipt
  write(99,100) 'totl ',totl
  write(99,100) 'rtc ',rtc
  write(99,100) 'rts ',rts
  write(99,100) 'ire ',ire
  write(99,100) 'ipr ',ipr
  write(99,100) 'irmod2 ',irmod2
  write(99,100) 'dtr ',dtr
  write(99,100) 'nre ',nre
  write(99,100) 'nur ',nur
  write(99,100) 'nch ',nch
  write(99,100) 'nqc ',nqc
  write(99,100) 'npp ',npp
  write(99,100) 'nrr ',nrr
  write(99,100) 'nu ',nu
  write(99,100) 'dphix ',dphix
  write(99,100) 'dphiz ',dphiz
  write(99,100) 'qx0 ',qx0
  write(99,100) 'qz0 ',qz0
  write(99,100) 'dres ',dres
  write(99,100) 'dfft ',dfft
  write(99,100) 'cma1 ',cma1
  write(99,100) 'cma2 ',cma2
  write(99,100) 'nstart ',nstart
  write(99,100) 'nstop ',nstop
  write(99,100) 'iskip ',iskip
  write(99,100) 'iconv ',iconv
  write(99,100) 'imad ',imad
  write(99,100) 'ipos ',ipos
  write(99,100) 'iav ',iav
  write(99,100) 'iwg ',iwg
  write(99,100) 'ivox ',ivox
  write(99,100) 'ivoz ',ivoz
  write(99,100) 'ires ',ires
  write(99,100) 'ifh ',ifh
  write(99,100) 'toptit ',toptit
  write(99,100) 'kwtype ',kwtype
  write(99,100) 'itf ',itf
  write(99,100) 'icr ',icr
  write(99,100) 'idis ',idis
  write(99,100) 'icow ',icow
  write(99,100) 'istw ',istw
  write(99,100) 'iffw ',iffw
  write(99,100) 'nprint ',nprint
  write(99,100) 'ndafi ',ndafi
  write(99,100) 'qwsk ',qwsk
  write(99,100) 'betx ',betx
  write(99,100) 'betz ',betz
  write(99,100) 'alfx ',alfx
  write(99,100) 'alfz ',alfz
  write(99,100) 'iskew ',iskew
  write(99,100) 'nskew ',nskew
  write(99,100) 'hmal ',hmal
  write(99,100) 'sixtit ',sixtit
  write(99,100) 'commen ',commen
  write(99,100) 'ithick ',ithick
  write(99,100) 'clo6 ',clo6
  write(99,100) 'clop6 ',clop6
  write(99,100) 'dki ',dki
  write(99,100) 'sigman ',sigman
  write(99,100) 'sigman2 ',sigman2
  write(99,100) 'sigmanq ',sigmanq
  write(99,100) 'clobeam ',clobeam
  write(99,100) 'beamoff ',beamoff
  write(99,100) 'parbe ',parbe
  write(99,100) 'track6d ',track6d
  write(99,100) 'ptnfac ',ptnfac
  write(99,100) 'sigz ',sigz
  write(99,100) 'sige ',sige
  write(99,100) 'partnum ',partnum
  write(99,100) 'parbe14 ',parbe14
  write(99,100) 'emitx ',emitx
  write(99,100) 'emity ',emity
  write(99,100) 'emitz ',emitz
  write(99,100) 'gammar ',gammar
  write(99,100) 'nbeam ',nbeam
  write(99,100) 'ibbc ',ibbc
  write(99,100) 'ibeco ',ibeco
  write(99,100) 'ibtyp ',ibtyp
  write(99,100) 'lhc ',lhc
  write(99,100) 'cotr ',cotr
  write(99,100) 'rrtr ',rrtr
  write(99,100) 'imtr ',imtr
  write(99,100) 'bbcu ',bbcu
  write(99,100) 'ibb6d ',ibb6d
  write(99,100) 'imbb ',imbb
  write(99,100) 'as ',as
  write(99,100) 'al ',al
  write(99,100) 'sigm ',sigm
  write(99,100) 'dps ',dps
  write(99,100) 'idz ',idz
  write(99,100) 'dp1 ',dp1
  write(99,100) 'itra ',itra
  write(99,100) 'x ',x
  write(99,100) 'y ',y
  write(99,100) 'bet0 ',bet0
  write(99,100) 'alf0 ',alf0
  write(99,100) 'clo ',clo
  write(99,100) 'clop ',clop
  write(99,100) 'cro ',cro
  write(99,100) 'is ',is
  write(99,100) 'ichrom ',ichrom
  write(99,100) 'nnumxv ',nnumxv
  write(99,100) 'xsi ',xsi
  write(99,100) 'zsi ',zsi
  write(99,100) 'smi ',smi
  write(99,100) 'ampt ',ampt
  write(99,100) 'tlim ',tlim
  write(99,100) 'tasm ',tasm
  write(99,100) 'preda ',preda
  write(99,100) 'idial ',idial
  write(99,100) 'nord ',nord
  write(99,100) 'nvar ',nvar
  write(99,100) 'nvar2 ',nvar2
  write(99,100) 'nsix ',nsix
  write(99,100) 'ncor ',ncor
  write(99,100) 'ipar ',ipar
  write(99,100) 'nordf ',nordf
  write(99,100) 'nvarf ',nvarf
  write(99,100) 'nord1 ',nord1
  write(99,100) 'ndimf ',ndimf
  write(99,100) 'idptr ',idptr
  write(99,100) 'inorm ',inorm
  write(99,100) 'imod1 ',imod1
  write(99,100) 'imod2 ',imod2
  write(99,100) 'ekv ',ekv
  write(99,100) 'fokqv ',fokqv
  write(99,100) 'aaiv ',aaiv
  write(99,100) 'bbiv ',bbiv
  write(99,100) 'smiv ',smiv
  write(99,100) 'zsiv ',zsiv
  write(99,100) 'xsiv ',xsiv
  write(99,100) 'xsv ',xsv
  write(99,100) 'zsv ',zsv
  write(99,100) 'qw ',qw
  write(99,100) 'qwc ',qwc
  write(99,100) 'clo0 ',clo0
  write(99,100) 'clop0 ',clop0
  write(99,100) 'eps ',eps
  write(99,100) 'epsa ',epsa
  write(99,100) 'ekk ',ekk
  write(99,100) 'cr ',cr
  write(99,100) 'ci ',ci
  write(99,100) 'xv ',xv
  write(99,100) 'yv ',yv
  write(99,100) 'dam ',dam
  write(99,100) 'ekkv ',ekkv
  write(99,100) 'sigmv ',sigmv
  write(99,100) 'dpsv ',dpsv
  write(99,100) 'dp0v ',dp0v
  write(99,100) 'sigmv6 ',sigmv6
  write(99,100) 'dpsv6 ',dpsv6
  write(99,100) 'ejv ',ejv
  write(99,100) 'ejfv ',ejfv
  write(99,100) 'xlv ',xlv
  write(99,100) 'zlv ',zlv
  write(99,100) 'pstop ',pstop
  write(99,100) 'rvv ',rvv
  write(99,100) 'ejf0v ',ejf0v
  write(99,100) 'numxv ',numxv
  write(99,100) 'nms ',nms
  write(99,100) 'partID ',partID
  write(99,100) 'parentID ',parentID
  write(99,100) 'dpd ',dpd
  write(99,100) 'dpsq ',dpsq
  write(99,100) 'fok ',fok
  write(99,100) 'rho ',rho
  write(99,100) 'fok1 ',fok1
  write(99,100) 'si ',si
  write(99,100) 'co ',co
  write(99,100) 'g ',g
  write(99,100) 'gl ',gl
  write(99,100) 'sm1 ',sm1
  write(99,100) 'sm2 ',sm2
  write(99,100) 'sm3 ',sm3
  write(99,100) 'sm12 ',sm12
  write(99,100) 'as3 ',as3
  write(99,100) 'as4 ',as4
  write(99,100) 'as6 ',as6
  write(99,100) 'sm23 ',sm23
  write(99,100) 'rhoc ',rhoc
  write(99,100) 'siq ',siq
  write(99,100) 'aek ',aek
  write(99,100) 'afok ',afok
  write(99,100) 'hp ',hp
  write(99,100) 'hm ',hm
  write(99,100) 'hc ',hc
  write(99,100) 'hs ',hs
  write(99,100) 'wf ',wf
  write(99,100) 'wfa ',wfa
  write(99,100) 'wfhi ',wfhi
  write(99,100) 'rhoi ',rhoi
  write(99,100) 'hi ',hi
  write(99,100) 'fi ',fi
  write(99,100) 'hi1 ',hi1
  write(99,100) 'xvl ',xvl
  write(99,100) 'yvl ',yvl
  write(99,100) 'ejvl ',ejvl
  write(99,100) 'dpsvl ',dpsvl
  write(99,100) 'oidpsv ',oidpsv
  write(99,100) 'sigmvl ',sigmvl
  write(99,100) 'iv ',iv
  write(99,100) 'aperv ',aperv
  write(99,100) 'ixv ',ixv
  write(99,100) 'ampv ',ampv
  write(99,100) 'clo6v ',clo6v
  write(99,100) 'clop6v ',clop6v
  write(99,100) 'hv ',hv
  write(99,100) 'bl1v ',bl1v
  write(99,100) 'tas ',tas
  write(99,100) 'qwcs ',qwcs
  write(99,100) 'di0xs ',di0xs
  write(99,100) 'di0zs ',di0zs
  write(99,100) 'dip0xs ',dip0xs
  write(99,100) 'dip0zs ',dip0zs
  write(99,100) 'xau ',xau
  write(99,100) 'cloau ',cloau
  write(99,100) 'di0au ',di0au
  write(99,100) 'tau ',tau
  write(99,100) 'tasau ',tasau
  write(99,100) 'wx ',wx
  write(99,100) 'x1 ',x1
  write(99,100) 'x2 ',x2
  write(99,100) 'fake ',fake
  write(99,100) 'e0f ',e0f
  write(99,100) 'numx ',numx
  write(99,100) 'cotr ',cotr
  write(99,100) 'rrtr ',rrtr
  write(99,100) 'imtr ',imtr

  ! these other values???
  write(99,100) 'numl ',numl
  write(99,100) 'niu ',niu
  write(99,100) 'amp0 ',amp0
  write(99,100) 'amp ',amp
  write(99,100) 'damp ',damp
  write(99,100) 'chi0 ',chi0
  write(99,100) 'chid ',chid
  write(99,100) 'rat ',rat
  write(99,100) 'exz ',exz
  write(99,100) 'time0 ',time0
  write(99,100) 'time1 ',time1

  endfile (99,iostat=ierro)
  backspace (99,iostat=ierro)

100 format (a10,(Z20))

end subroutine dumphex
#endif
