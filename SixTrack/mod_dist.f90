! ================================================================================================ !
!  A.Mereghetti and D.Sinuela Pastor, for the FLUKA Team
!  last modified: 2018-06-02
!  read a beam distribution
!  always in main code
!
!  Format of the input file:
!    id:         unique identifier of the particle (integer)
!    gen:        parent ID (integer)
!    weight:     statistical weight of the particle (double: >0.0)
!    x,  y,  s:  particle position  [m]
!    xp, yp, zp: particle direction (tangents) []
!    aa, zz:     mass and atomic number
!    m:          rest mass [GeV/c2]
!    pc:         particle momentum [GeV/c]
!    dt:         time delay with respect to the reference particle [s]
!  - aa,zz and m are now taken into account for hisix
!
!  NOTA BENE:
!  - id, gen and weight are assigned by the fluka_mod_init subroutine;
!  - z and zp are actually useless (but we never know);
!  - aa, zz and m are not stored at the moment (safer decision from the code
!    point of view, until a decision about ion tracking is taken);
!    the subroutine fluka_send is then responsible for using the corresponding
!    values for protons through the interface, whereas the subroutine fluka_receive
!    simply ignores the values passed through the FlukaIO interface;
!
!  variables in input to routine:
!  - napx: number of protons to be tracked (from fort.3 file);
!  - npart: max number of protons that can be tracked (array dimensioning);
!  - enom: nominal total energy of the beam (ie of synch particle) [MeV];
!  - pnom: nominal linear momentum of the beam (ie of synch particle) [MeV/c];
!  - clight: speed of light [m/s];
!  NB: in case the file contains less particle than napx, napx is
!      re-assigned
!
!  output variables:
!    all other variables in the interface (6D tracking variables);
!
! ================================================================================================ !
module mod_dist

  use crcoall
  use floatPrecision

  implicit none

  logical,           public, save :: dist_enable      ! DIST input block given
  logical,           public, save :: dist_echo        ! echo the read distribution?
  character(len=16), public, save :: dist_filename    !
  integer,           public, save :: dist_read_unit   ! unit for reading the distribution
  integer,           public, save :: dist_echo_unit   ! unit for echoing the distribution

contains

subroutine dist_parseInputLine(inLine, iLine, iErr)

  use string_tools
  use file_units

  implicit none

  character(len=*), intent(in)    :: inLine
  integer,          intent(inout) :: iLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable   :: lnSplit(:)
  integer nSplit
  logical spErr

  call chr_split(inLine, lnSplit, nSplit, spErr)
  if(spErr) then
    write(lout,"(a)") "DIST> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  if(nSplit == 0) return

  select case(lnSplit(1))

  case("ECHO")
    dist_echo = .true.

  case("RDUN")
    write(lout,"(a)") "DIST> INFO RDUN is deprecated. A unit will be assigned automatically."

  case("ECUN")
    if(nSplit <= 2) then
      write(lout,"(a,i0)") "DIST> ERROR ECUN must have 1 values, got ",(nSplit-1)
      iErr = .true.
      return
    end if
    call chr_cast(lnSplit(2),dist_echo_unit,iErr)

  case("READ")
    if(nSplit <= 2) then
      write(lout,"(a,i0)") "DIST> ERROR READ must have 1 values, got ",(nSplit-1)
      iErr = .true.
      return
    end if
    dist_filename = trim(lnSplit(2))
    call funit_requestUnit(dist_filename,dist_read_unit)
    if(.not.dist_enable) dist_enable = .true.

  end select

end subroutine dist_parseInputLine

subroutine dist_readdis(napx, npart, enom, pnom, clight, x, y, xp, yp, s, pc, aa, zz, m)

  use numerical_constants
  use, intrinsic :: iso_fortran_env, only : int16
  implicit none

! interface variables:
  integer napx, npart
  real(kind=fPrec) enom, pnom, clight
  real(kind=fPrec) x, y, xp, yp, s, pc
  real(kind=fPrec) m0

! temporary variables:
  integer id, gen
  integer(kind=int16) aa, zz
  real(kind=fPrec) weight, z, zp, m, dt
  integer jj
  character(240) tmp_line

  character comment_char
  parameter ( comment_char = '*' )

  dimension x ( npart ), y ( npart )
  dimension xp( npart ), yp( npart )
  dimension s ( npart ), pc( npart )
  dimension dt( npart )
! P. HERMES hisix
  dimension aa( npart ), zz( npart )
  dimension m ( npart )

  write(lout,*) ''
  write(lout,*) "Reading particles from ", dist_filename

! initialise tracking variables:
  do jj=1,npart
    x (jj) = zero
    y (jj) = zero
    xp(jj) = zero
    yp(jj) = zero
    pc(jj) = zero
    s (jj) = zero
    aa(jj) = 0            ! hisix
    zz(jj) = 0            ! hisix
    m (jj) = zero         ! hisix
  end do

! initialise particle counter
  jj = 0

  open( unit=dist_read_unit, file=dist_filename )

! cycle on lines in file:
1981 continue
  read(dist_read_unit,'(A)',end=1983,err=1982) tmp_line
  if( tmp_line(1:1).eq.comment_char ) goto 1981
  jj = jj+1

  if( jj.gt.npart ) then
    write(lout,*) 'Error while reading particles'
    write(lout,*) 'not enough memory for all particles in file'
    write(lout,*) 'please increase the npart parameter and recompile'
    write(lout,*) 'present value:', npart
    jj = npart
    goto 1984
  else if( jj.gt.napx ) then
    write(lout,*) ''
    write(lout,*) 'Stopping reading file, as already ', napx
    write(lout,*) ' particles have been read, as requested by the user'
    write(lout,*) ' in fort.3 file'
    write(lout,*) ''
    jj = napx
    goto 1983
  end if

  read( tmp_line, *, err=1982 ) id, gen, weight, x(jj), y(jj), z, xp(jj), yp(jj), zp, aa(jj), zz(jj), m(jj), pc(jj), dt(jj)
  goto 1981

! error while parsing file:
1982 continue
  write(lout,*) 'Error while reading particles at line:'
  write(lout,*) tmp_line
  goto 1984

1983 continue
  if( jj.eq.0 ) then
    write(lout,*) 'Error while reading particles'
    write(lout,*) 'no particles read from file'
    goto 1984
  end if

  close(dist_read_unit)
  write(lout,*) "Number of particles read = ", jj

  if( jj.lt.napx ) then
    write(lout,*) ''
    write(lout,*) 'Warning: read a number of particles'
    write(lout,*) '         LOWER than the one requested for tracking'
    write(lout,*) '         requested:',napx
    write(lout,*) ''
    napx = jj
  end if

! fix units:
  do jj=1,napx
    x (jj) = x(jj)  * c1e3 ! [m]     -> [mm]
    y (jj) = y(jj)  * c1e3 ! [m]     -> [mm]
    xp(jj) = xp(jj) * c1e3 ! []      -> [1.0E-03]
    yp(jj) = yp(jj) * c1e3 ! []      -> [1.0E-03]
    pc(jj) = pc(jj) * c1e3 ! [GeV/c] -> [MeV/c]
    m (jj) = m(jj)  * c1e3 ! [GeV/c^2] -> [MeV/c^2] ! P. HERMES
    s (jj) = -pnom/enom * dt(jj)*clight * c1e3
  end do

  return

! exit with error
1984 continue
  close(dist_read_unit)
  call prror(-1)
  return
end subroutine dist_readdis

end module mod_dist
