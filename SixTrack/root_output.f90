module root_output
  use, intrinsic :: iso_c_binding
  use crcoall

  implicit none

  logical root_flag                                      !ROOT input block exists
  integer root_eos_enabled
  integer root_Accelerator
  integer root_ApertureCheck
  integer root_Collimation
  integer root_CollimationDB
  integer root_Optics
  integer root_RunNumber
  character(len=512) :: root_eos_server
  character(len=512) :: root_folder
  character(len=512) :: root_prefix

  save

interface

!General stuff
subroutine DoSixTrackRootInit(eos, run_number, eos_server, root_path, root_prefix, Accelerator, Optics, ApertureCheck, Collimation,&
& CollimationDB) bind(C,name="DoSixTrackRootInit")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(kind=C_INT), intent(in), value :: eos
  integer(kind=C_INT), intent(in), value :: run_number
  character(kind=C_CHAR,len=1), intent(in) :: eos_server
  character(kind=C_CHAR,len=1), intent(in) :: root_path
  character(kind=C_CHAR,len=1), intent(in) :: root_prefix
  integer(kind=C_INT), intent(in), value :: Accelerator
  integer(kind=C_INT), intent(in), value :: Optics
  integer(kind=C_INT), intent(in), value :: ApertureCheck
  integer(kind=C_INT), intent(in), value :: Collimation
  integer(kind=C_INT), intent(in), value :: CollimationDB
end subroutine DoSixTrackRootInit

subroutine SixTrackRootExit() bind(C,name="SixTrackRootExit")
  use, intrinsic :: iso_c_binding
  implicit none
end subroutine SixTrackRootExit

!collimation stuff
subroutine CollimatorLossRootWrite(icoll_in,db_name,db_name_len,impact_in,absorbed_in,caverage_in,csigma_in,length_in) &
 & bind(C,name="CollimatorLossRootWrite")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(kind=C_INT), intent(in), value :: icoll_in
  character(kind=C_CHAR,len=1), intent(in) :: db_name
  integer(kind=C_INT), intent(in), value :: db_name_len
  integer(kind=C_INT), intent(in), value :: impact_in
  integer(kind=C_INT), intent(in), value :: absorbed_in
  real(kind=C_DOUBLE), intent(in), value :: caverage_in
  real(kind=C_DOUBLE), intent(in), value :: csigma_in
  real(kind=C_DOUBLE), intent(in), value :: length_in
end subroutine CollimatorLossRootWrite

subroutine ApertureCheckWriteLossParticle(turn_in, i_in, ix_in, bez_in, bez_in_len, slos_in, ipart_in, x_in, xp_in, y_in, yp_in, &
& ct_in, e_in, dp_in) bind(C,name="ApertureCheckWriteLossParticle")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(kind=C_INT), intent(in), value :: turn_in
  integer(kind=C_INT), intent(in), value :: i_in
  integer(kind=C_INT), intent(in), value :: ix_in
  character(kind=C_CHAR,len=1), intent(in) :: bez_in
  integer(kind=C_INT), intent(in), value :: bez_in_len
  real(kind=C_DOUBLE), intent(in), value :: slos_in
  integer(kind=C_INT), intent(in), value :: ipart_in
  real(kind=C_DOUBLE), intent(in), value :: x_in
  real(kind=C_DOUBLE), intent(in), value :: xp_in
  real(kind=C_DOUBLE), intent(in), value :: y_in
  real(kind=C_DOUBLE), intent(in), value :: yp_in
  real(kind=C_DOUBLE), intent(in), value :: ct_in
  real(kind=C_DOUBLE), intent(in), value :: e_in
  real(kind=C_DOUBLE), intent(in), value :: dp_in
end subroutine ApertureCheckWriteLossParticle

subroutine SurvivalRootWrite(nturn_in, npart_in) bind(C,name="SurvivalRootWrite")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(kind=C_INT), intent(in), value :: nturn_in
  integer(kind=C_INT), intent(in), value :: npart_in
end subroutine SurvivalRootWrite

subroutine OpticsRootWriteLin(i_in, name_in, c_len, s_in, x_in, xp_in, y_in, yp_in, beta_x_in, beta_y_in, alpha_x_in, alpha_y_in, &
& dispersion_x_in, dispersion_y_in, dispersion_xp_in, dispersion_yp_in) bind(C,name="OpticsRootWriteLin")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(kind=C_INT), intent(in), value :: i_in
  character(kind=C_CHAR,len=1), intent(in) :: name_in
  integer(kind=C_INT), intent(in), value :: c_len
  real(kind=C_DOUBLE), intent(in), value :: s_in
  real(kind=C_DOUBLE), intent(in), value :: x_in
  real(kind=C_DOUBLE), intent(in), value :: xp_in
  real(kind=C_DOUBLE), intent(in), value :: y_in
  real(kind=C_DOUBLE), intent(in), value :: yp_in
  real(kind=C_DOUBLE), intent(in), value :: beta_x_in
  real(kind=C_DOUBLE), intent(in), value :: beta_y_in
  real(kind=C_DOUBLE), intent(in), value :: alpha_x_in
  real(kind=C_DOUBLE), intent(in), value :: alpha_y_in
  real(kind=C_DOUBLE), intent(in), value :: dispersion_x_in
  real(kind=C_DOUBLE), intent(in), value :: dispersion_y_in
  real(kind=C_DOUBLE), intent(in), value :: dispersion_xp_in
  real(kind=C_DOUBLE), intent(in), value :: dispersion_yp_in
end subroutine

subroutine OpticsRootWriteCpl(phi1,phi2,bexi,bexii,bezi,bezii,   &
 &                                  alxi,alxii,alzi,alzii,       &
 &                                  gaxi,gaxii,gazi,gazii,       &
 &                                  phxi,phxii,phzi,phzii,       &
 &                                  phxpi,phxpii,phzpi,phzpii,   &
 &                                  couuang,                     &
 &                                  t61,t62,t63,t64, &
 &                                  t11,t12,t13,t14) bind(C,name="OpticsRootWriteCpl")
  use, intrinsic :: iso_c_binding

  implicit none

  real(kind=C_DOUBLE), intent(in), value :: phi1
  real(kind=C_DOUBLE), intent(in), value :: phi2

  real(kind=C_DOUBLE), intent(in), value :: bexi
  real(kind=C_DOUBLE), intent(in), value :: bexii
  real(kind=C_DOUBLE), intent(in), value :: bezi
  real(kind=C_DOUBLE), intent(in), value :: bezii

  real(kind=C_DOUBLE), intent(in), value :: alxi
  real(kind=C_DOUBLE), intent(in), value :: alxii
  real(kind=C_DOUBLE), intent(in), value :: alzi
  real(kind=C_DOUBLE), intent(in), value :: alzii

  real(kind=C_DOUBLE), intent(in), value :: gaxi
  real(kind=C_DOUBLE), intent(in), value :: gaxii
  real(kind=C_DOUBLE), intent(in), value :: gazi
  real(kind=C_DOUBLE), intent(in), value :: gazii

  real(kind=C_DOUBLE), intent(in), value :: phxi
  real(kind=C_DOUBLE), intent(in), value :: phxii
  real(kind=C_DOUBLE), intent(in), value :: phzi
  real(kind=C_DOUBLE), intent(in), value :: phzii

  real(kind=C_DOUBLE), intent(in), value :: phxpi
  real(kind=C_DOUBLE), intent(in), value :: phxpii
  real(kind=C_DOUBLE), intent(in), value :: phzpi
  real(kind=C_DOUBLE), intent(in), value :: phzpii

  real(kind=C_DOUBLE), intent(in), value :: couuang

  real(kind=C_DOUBLE), intent(in), value :: t61
  real(kind=C_DOUBLE), intent(in), value :: t62
  real(kind=C_DOUBLE), intent(in), value :: t63
  real(kind=C_DOUBLE), intent(in), value :: t64
  real(kind=C_DOUBLE), intent(in), value :: t11
  real(kind=C_DOUBLE), intent(in), value :: t12
  real(kind=C_DOUBLE), intent(in), value :: t13
  real(kind=C_DOUBLE), intent(in), value :: t14
end subroutine

subroutine  OpticsRootWrite() bind(C,name="OpticsRootWrite")
  use, intrinsic :: iso_c_binding
end subroutine

!extern "C" void CollimatorDatabaseRootWrite(int j, char* db_name_in, int db_name_len, char* db_material_in, int db_material_len, double db_nsig_in, double db_length_in, double db_rotation_in, double db_offset_in)
subroutine  CollimatorDatabaseRootWrite(j, db_name_in, db_name_len, db_material_in, db_material_len, db_nsig_in, db_length_in, &
db_rotation_in, db_offset_in) bind(C,name="CollimatorDatabaseRootWrite")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(kind=C_INT), intent(in), value :: j
  character(kind=C_CHAR,len=1), intent(in) :: db_name_in
  integer(kind=C_INT), intent(in), value :: db_name_len
  character(kind=C_CHAR,len=1), intent(in) :: db_material_in
  integer(kind=C_INT), intent(in), value :: db_material_len
  real(kind=C_DOUBLE), intent(in), value :: db_nsig_in
  real(kind=C_DOUBLE), intent(in), value :: db_length_in
  real(kind=C_DOUBLE), intent(in), value :: db_rotation_in
  real(kind=C_DOUBLE), intent(in), value :: db_offset_in
end subroutine

!extern "C" void RunTimeRootWrite(Float_t pretime_in, Float_t trtime_in, Float_t, posttime_in);
subroutine  RunTimeRootWrite(pretime_in, trtime_in, posttime_in) bind(C,name="RunTimeRootWrite")
  use, intrinsic :: iso_c_binding
  implicit none
  real(kind=C_FLOAT), intent(in), value :: pretime_in
  real(kind=C_FLOAT), intent(in), value :: trtime_in
  real(kind=C_FLOAT), intent(in), value :: posttime_in
end subroutine

subroutine SixTrackRootWrite() bind(C,name="SixTrackRootWrite")
  implicit none
end subroutine

!extern "C" void AcceleratorRootWrite(char* name_in, int name_len, int ktrack_in, double value_in, double extra_in, double length_in);
subroutine  AcceleratorRootWrite(name_in, name_len, ktrack_in, value_in, extra_in, length_in) bind(C,name="AcceleratorRootWrite")
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=C_CHAR,len=1), intent(in) :: name_in
  integer(kind=C_INT), intent(in), value :: name_len
  integer(kind=C_INT), intent(in), value :: ktrack_in
  real(kind=C_DOUBLE), intent(in), value :: value_in
  real(kind=C_DOUBLE), intent(in), value :: extra_in
  real(kind=C_DOUBLE), intent(in), value :: length_in
end subroutine


subroutine  ConfigurationOutputRootSet_npart(napx_in) bind(C,name="ConfigurationOutputRootSet_npart")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(kind=C_INT), intent(in), value :: napx_in
end subroutine

subroutine  ConfigurationOutputRootSet_nturns(nturns_in) bind(C,name="ConfigurationOutputRootSet_nturns")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(kind=C_INT), intent(in), value :: nturns_in
end subroutine

subroutine ConfigurationRootWrite() bind(C,name="ConfigurationRootWrite")
  implicit none
end subroutine


end interface

contains

subroutine SixTrackRootInit
  implicit none
  call DoSixTrackRootInit(root_eos_enabled, root_RunNumber, root_eos_server, root_folder, root_prefix, root_Accelerator, &
& root_Optics, root_ApertureCheck, root_Collimation, root_CollimationDB)
end subroutine SixTrackRootInit

subroutine SixTrackRootFortranInit
  implicit none
  root_flag          = .false.
  root_eos_enabled   = 0
  root_Accelerator   = 0
  root_ApertureCheck = 0
  root_Collimation   = 0
  root_CollimationDB = 0
  root_Optics        = 0
  root_RunNumber     = 0
  root_eos_server    = C_NULL_CHAR
  root_folder        = C_NULL_CHAR
  root_prefix        = C_NULL_CHAR
end subroutine SixTrackRootFortranInit

subroutine daten_root(ch)

  use string_tools

  implicit none

  character(len=*), intent(in) :: ch
  character getfields_fields(getfields_n_max_fields)*(getfields_l_max_string) ! Array of fields
  integer   getfields_nfields                                                 ! Number of identified fields
  integer   getfields_lfields(getfields_n_max_fields)                         ! Length of each what:
  logical   getfields_lerr                                                    ! An error flag

!ROOT is enabled
  root_flag = .true.

  !Read filenames
  call getfields_split( ch, getfields_fields, getfields_lfields, getfields_nfields, getfields_lerr )
  if( getfields_lerr ) call prror(-1)

  if(getfields_nfields .ne. 2) then
     write(lout,'(a)')         "ERROR in ROOT input:"
     write(lout,'(a,1x,i3,a)') "Expected 2 entries per line, got", getfields_nfields, ", line=",ch
     call prror(-1)
  end if

!  For input debugging if needed
!  write(lout,*) '1: ', getfields_fields(1)(1:getfields_lfields(1))
!  write(lout,*) '2: ', getfields_fields(2)(1:getfields_lfields(2))

!EOS: e.g. eosuser.cern.ch/
  if(getfields_fields(1)(1:getfields_lfields(1)).eq.'EOS') then
    root_eos_enabled = 1
    root_eos_server = getfields_fields(2)(1:getfields_lfields(2)) // C_NULL_CHAR
!PATH e.g. /eos/user/u/username/
  else if(getfields_fields(1)(1:getfields_lfields(1)).eq.'PATH') then
    root_folder = getfields_fields(2)(1:getfields_lfields(2)) // C_NULL_CHAR
!PREFIX sixtrack_
  else if(getfields_fields(1)(1:getfields_lfields(1)).eq.'PREFIX') then
    root_prefix = getfields_fields(2)(1:getfields_lfields(2)) // C_NULL_CHAR
!RUN number
  else if(getfields_fields(1)(1:getfields_lfields(1)).eq.'RUN') then
    read(getfields_fields(2)(1:getfields_lfields(2)),*) root_RunNumber
!blocks to enable
!ENABLE
!COLL, APER, ALL
  else if(getfields_fields(1)(1:getfields_lfields(1)).eq.'ENABLE') then
    write(lout,*) getfields_fields(2)(1:getfields_lfields(2))
    if(getfields_fields(2)(1:getfields_lfields(2)).eq.'ALL') then
      root_ApertureCheck = 1
      root_Accelerator = 1
      root_Collimation = 1
      root_CollimationDB = 1
      root_Optics = 1
    else if(getfields_fields(2)(1:getfields_lfields(2)).eq.'ACCEL') then
      root_Accelerator = 1
    else if(getfields_fields(2)(1:getfields_lfields(2)).eq.'COLL') then
      root_Collimation = 1
    else if(getfields_fields(2)(1:getfields_lfields(2)).eq.'COLDB') then
      root_CollimationDB = 1
    else if(getfields_fields(2)(1:getfields_lfields(2)).eq.'APER') then
      root_ApertureCheck = 1
    else if(getfields_fields(2)(1:getfields_lfields(2)).eq.'OPTICS') then
      root_Optics = 1
    end if
  end if

!Check configuration

!check + enable flags

end subroutine daten_root

subroutine root_parseInputDone
  implicit none
  write(lout,*) 'ROOT parsing done'
  if(root_eos_enabled.eq.1) then
    write(lout,*) 'Will write to eos directly: ', root_eos_enabled
    write(lout,*) 'server: ', root_eos_server
  else
    write(lout,*) 'Will write to -'
  end if

  write(lout,*) 'path:   ', root_folder
  write(lout,*) 'prefix: ', root_prefix
  write(lout,*) 'run:    ', root_RunNumber
  write(lout,*) 'Accelerator:   ', root_Accelerator
  write(lout,*) 'Optics:        ', root_Optics
  write(lout,*) 'Collimation:   ', root_Collimation
  write(lout,*) 'CollimationDB: ', root_CollimationDB
  write(lout,*) 'Aperture:      ', root_ApertureCheck

end subroutine root_parseInputDone

end module root_output
