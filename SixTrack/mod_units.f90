module mod_units
  
  implicit none
  
  type, public :: unitSpec
    integer,            public :: unit     ! Unit number
    character(len=256), public :: filename ! Standard filename
    integer,            public :: form     ! 0 = unformatted, 1 = formatted
    logical,            public :: sixtrack ! Default open in standard SixTrack
    logical,            public :: sixda    ! Default open in differential algebra version
    logical,            public :: boinc    ! Open in BOINC if enabled
    integer,            public :: recl=0   ! RECL flag, if needed
  end type unitSpec
  
  type(unitSpec), private :: units_uList(200)
  integer,        private :: units_nList
  
contains

subroutine units_initUnits
  units_nList = 0
end subroutine units_initUnits

subroutine units_openDefaults
  
  implicit none
  
  integer i
  character(len=7) fileName
  
  call units_openUnits(unit=2, fileName="fort.2", formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=3, fileName="fort.3", formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=4, fileName="fort.4", formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=7, fileName="fort.7", formatted=.true., boinc=.true.,fio=.true.,recl=303)
  call units_openUnits(unit=8, fileName="fort.8", formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=9, fileName="fort.9", formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=11,fileName="fort.11",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=12,fileName="fort.12",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=13,fileName="fort.13",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=14,fileName="fort.14",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=15,fileName="fort.15",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=16,fileName="fort.16",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=17,fileName="fort.17",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=18,fileName="fort.18",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=19,fileName="fort.19",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=20,fileName="fort.20",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=21,fileName="fort.21",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=22,fileName="fort.22",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=23,fileName="fort.23",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=24,fileName="fort.24",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=25,fileName="fort.25",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=26,fileName="fort.26",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=27,fileName="fort.27",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=28,fileName="fort.28",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=29,fileName="fort.29",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=30,fileName="fort.30",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=31,fileName="fort.31",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=32,fileName="fort.32",formatted=.false.,boinc=.true.)
  call units_openUnits(unit=33,fileName="fort.33",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=34,fileName="fort.34",formatted=.true., boinc=.true.,fio=.true.)
  call units_openUnits(unit=35,fileName="fort.35",formatted=.true., boinc=.true.,fio=.true.)
  
#ifndef BNLELENS
#ifdef STF
  ! Open Single Track File
  call units_openUnits(unit=90,fileName="singletrackfile.dat",formatted=.false.,boinc=.true.)
#else
  ! Open binary files 59 to 90 for particle pair 1 to 32
  do i=59,90
    write(fileName,"(a5,i2)") "fort.",i
    call units_openUnits(unit=i,fileName=fileName,formatted=.false.,boinc=.true.)
  end do
#endif
#else
#ifdef BOINC
  call units_openUnits(unit=10,fileName="fort.10",formatted=.true.,boinc=.true.,fio=.true.,recl=8195)
  call units_openUnits(unit=54,fileName="fort.54",formatted=.true.,boinc=.true.,fio=.true.)
#else
#ifdef CR
  call units_openUnits(unit=51,fileName="fort.51",formatted=.true.,boinc=.false.,fio=.true.)
  call units_openUnits(unit=52,fileName="fort.52",formatted=.true.,boinc=.false.,fio=.true.)
  call units_openUnits(unit=53,fileName="fort.53",formatted=.true.,boinc=.false.,fio=.true.)
  call units_openUnits(unit=54,fileName="fort.54",formatted=.true.,boinc=.false.,fio=.true.)
  call units_openUnits(unit=97,fileName="fort.97",formatted=.true.,boinc=.false.,fio=.true.)
#else
  call units_openUnits(unit=51,fileName="SixTwiss.dat",formatted=.true.,boinc=.false.,fio=.true.)
  call units_openUnits(unit=52,fileName="beambeam-output.dat",formatted=.true.,boinc=.false.,fio=.true.)
  call units_openUnits(unit=52,fileName="beambeam-lostID.dat",formatted=.true.,boinc=.false.,fio=.true.)
  call units_openUnits(unit=54,fileName="beambeamdist.dat",formatted=.true.,boinc=.false.,fio=.true.)
  call units_openUnits(unit=97,fileName="checkdist.dat",formatted=.true.,boinc=.false.,fio=.true.)
#endif
#endif
#endif

  call units_openUnits(unit=98,fileName="fort.98",formatted=.true.,boinc=.true.,fio=.true.)
  
  ! Eric for the DA coefficients in BINARY
  call units_openUnits(unit=110,fileName="fort.110",formatted=.false.,boinc=.false.)
  call units_openUnits(unit=111,fileName="fort.111",formatted=.false.,boinc=.false.)
  
#ifdef DEBUG
  call units_openUnits(unit=99 ,fileName="dump",  formatted=.false.,boinc=.true.)
  call units_openUnits(unit=100,fileName="arrays",formatted=.false.,boinc=.true.)
#endif
  
  ! Heavy Ion Output
  call units_openUnits(unit=208,fileName="fort.208",formatted=.false.,boinc=.true.) ! collimator losses (energy)
  call units_openUnits(unit=209,fileName="fort.209",formatted=.false.,boinc=.true.) ! collimator losses in function of particle i
  call units_openUnits(unit=210,fileName="fort.210",formatted=.false.,boinc=.true.) ! mtc after each collimator interaction
  
end subroutine units_openDefaults

subroutine units_openUnits(unit,fileName,formatted,boinc,fio,recl)
  
  implicit none
  
  integer,                    intent(in) :: unit
  character(len=*),           intent(in) :: fileName
  logical,                    intent(in) :: formatted
  logical,          optional, intent(in) :: boinc
  logical,          optional, intent(in) :: fio
  integer,          optional, intent(in) :: recl
  
  integer i, fRecl
  logical fBoinc, fFio
  
  if(present(boinc)) then
    fBoinc = boinc
  else
    fBoinc = .false.
  end if
  
  if(present(fio)) then
    fFio = fio
  else
    fFio = .false.
  end if
  
  if(present(recl)) then
    fRecl = recl
  else
    fRecl = 0
  end if
  
#ifdef BOINC
  if(fBoinc) call boincrf(unitList(i)%filename, fileName)
#endif
#ifndef FIO
  fFio = .false.
#endif
#ifndef NAGFOR
  fRecl = 0
#endif
  
  if(.not. formatted) fFio = .false.
  
  if(formatted) then
    if(fRecl > 0) then
      if(fFio) then
        open(unit,file=fileName,form="formatted",status="unknown",round="nearest",recl=fRecl)
      else
        open(unit,file=fileName,form="formatted",status="unknown",recl=fRecl)
      end if
    else
      if(fFio) then
        open(unit,file=fileName,form="formatted",status="unknown",round="nearest")
      else
        open(unit,file=fileName,form="formatted",status="unknown")
      end if
    end if
  else
    open(unit,file=fileName,form="unformatted",status="unknown")
  endif
  
end subroutine units_openUnits

end module mod_units
