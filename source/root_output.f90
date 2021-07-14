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
  integer root_FLUKA
  integer root_DumpPipe
  integer root_RunNumber
  character(len=512) :: root_eos_server
  character(len=512) :: root_folder
  character(len=512) :: root_prefix

  save

interface

!General stuff
subroutine DoSixTrackRootInit(eos, run_number, eos_server, root_path, root_prefix, Accelerator, Optics, ApertureCheck, Collimation,&
& CollimationDB, FLUKA_f, ApertureDump) bind(C,name="DoSixTrackRootInit")
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
  integer(kind=C_INT), intent(in), value :: FLUKA_f
  integer(kind=C_INT), intent(in), value :: ApertureDump
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
& p_in, dp_in, ct_in, naa_in, nzz_in, nqq_in, pdgid_in) bind(C,name="ApertureCheckWriteLossParticle")
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
  real(kind=C_DOUBLE), intent(in), value :: p_in
  real(kind=C_DOUBLE), intent(in), value :: dp_in
  real(kind=C_DOUBLE), intent(in), value :: ct_in
  integer(kind=C_INT), intent(in), value :: naa_in
  integer(kind=C_INT), intent(in), value :: nzz_in
  integer(kind=C_INT), intent(in), value :: nqq_in
  integer(kind=C_INT), intent(in), value :: pdgid_in
end subroutine ApertureCheckWriteLossParticle

subroutine ApertureCheckWriteLossParticleF(turn_in, i_in, ix_in, bez_in, bez_in_len, slos_in, partID_in, parentID_in, &
& partWeight_in, x_in, xp_in, y_in, yp_in, p_in, dp_in, ct_in, naa_in, nzz_in, nqq_in, pdgid_in) &
& bind(C,name="ApertureCheckWriteLossParticleF")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(kind=C_INT), intent(in), value :: turn_in
  integer(kind=C_INT), intent(in), value :: i_in
  integer(kind=C_INT), intent(in), value :: ix_in
  character(kind=C_CHAR,len=1), intent(in) :: bez_in
  integer(kind=C_INT), intent(in), value :: bez_in_len
  real(kind=C_DOUBLE), intent(in), value :: slos_in
  integer(kind=C_INT32_T), intent(in), value :: partID_in
  integer(kind=C_INT32_T), intent(in), value :: parentID_in
  real(kind=C_DOUBLE), intent(in), value :: partWeight_in
  real(kind=C_DOUBLE), intent(in), value :: x_in
  real(kind=C_DOUBLE), intent(in), value :: xp_in
  real(kind=C_DOUBLE), intent(in), value :: y_in
  real(kind=C_DOUBLE), intent(in), value :: yp_in
  real(kind=C_DOUBLE), intent(in), value :: p_in
  real(kind=C_DOUBLE), intent(in), value :: dp_in
  real(kind=C_DOUBLE), intent(in), value :: ct_in
  integer(kind=C_INT), intent(in), value :: naa_in
  integer(kind=C_INT), intent(in), value :: nzz_in
  integer(kind=C_INT), intent(in), value :: nqq_in
  integer(kind=C_INT), intent(in), value :: pdgid_in
end subroutine ApertureCheckWriteLossParticleF

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
  real(kind=C_DOUBLE), intent(in), value :: pretime_in
  real(kind=C_DOUBLE), intent(in), value :: trtime_in
  real(kind=C_DOUBLE), intent(in), value :: posttime_in
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

subroutine  ConfigurationOutputRootSet_aperture_binsize(bin_size_in) bind(C,name="ConfigurationOutputRootSet_aperture_binsize")
  use, intrinsic :: iso_c_binding
  implicit none
  real(kind=C_DOUBLE), intent(in), value :: bin_size_in
end subroutine

subroutine  ConfigurationOutputRootSet_reference_energy(e0_in) bind(C,name="ConfigurationOutputRootSet_reference_energy")
  use, intrinsic :: iso_c_binding
  implicit none
  real(kind=C_DOUBLE), intent(in), value :: e0_in
end subroutine

subroutine  ConfigurationOutputRootSet_reference_mass(nucm0_in) bind(C,name="ConfigurationOutputRootSet_reference_mass")
  use, intrinsic :: iso_c_binding
  implicit none
  real(kind=C_DOUBLE), intent(in), value :: nucm0_in
end subroutine

subroutine ConfigurationRootWrite() bind(C,name="ConfigurationRootWrite")
  implicit none
end subroutine

subroutine root_EnergyDeposition(id_in, nucleons_in, energy_in) bind(C,name="root_EnergyDeposition")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(kind=C_INT),     intent(in), value :: id_in
  integer(kind=C_INT16_T), intent(in), value :: nucleons_in
  real(kind=C_DOUBLE),     intent(in), value :: energy_in
end subroutine

subroutine root_FLUKA_Names(id_in, name_in, name_len, ins_type) bind(C,name="root_FLUKA_Names")
  use, intrinsic :: iso_c_binding
  implicit none
  integer(kind=C_INT),    intent(in), value :: id_in
  character(kind=C_CHAR,len=1), intent(in)  :: name_in
  integer(kind=C_INT),    intent(in), value :: name_len
  integer(kind=C_INT),    intent(in), value :: ins_type
end subroutine

subroutine root_DumpAperture(apname_in, apname_len, aptype_in, aptype_len, s_in, ap1_in, ap2_in, ap3_in, ap4_in, ap5_in, ap6_in, &
& ap7_in, ap8_in, ap9_in, ap10_in, ap11_in) bind(C,name="root_DumpAperture")
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=C_CHAR,len=1), intent(in)  :: apname_in
  integer(kind=C_INT),    intent(in), value :: apname_len
  character(kind=C_CHAR,len=1), intent(in)  :: aptype_in
  integer(kind=C_INT),    intent(in), value :: aptype_len
  real(kind=C_DOUBLE), intent(in), value :: s_in
  real(kind=C_DOUBLE), intent(in), value :: ap1_in
  real(kind=C_DOUBLE), intent(in), value :: ap2_in
  real(kind=C_DOUBLE), intent(in), value :: ap3_in
  real(kind=C_DOUBLE), intent(in), value :: ap4_in
  real(kind=C_DOUBLE), intent(in), value :: ap5_in
  real(kind=C_DOUBLE), intent(in), value :: ap6_in
  real(kind=C_DOUBLE), intent(in), value :: ap7_in
  real(kind=C_DOUBLE), intent(in), value :: ap8_in
  real(kind=C_DOUBLE), intent(in), value :: ap9_in
  real(kind=C_DOUBLE), intent(in), value :: ap10_in
  real(kind=C_DOUBLE), intent(in), value :: ap11_in
end subroutine

subroutine root_BunchDumpInit() bind(C,name="root_BunchDumpInit")
  use, intrinsic :: iso_c_binding
  implicit none
end subroutine

subroutine root_DumpBunch( name_in, name_len, i,ix, turn_in, particleID_in, parentID_in, pdgID_in, q_in, weight_in, &
& s_in, x_in, xp_in, y_in, yp_in, z_in, dp_in, sx_in, sy_in, sz_in, m_in) bind(C,name="root_DumpBunch")
  use, intrinsic :: iso_c_binding
  implicit none
  character(kind=C_CHAR,len=1), intent(in)  :: name_in
  integer(kind=C_INT),    intent(in), value :: name_len
  integer(kind=C_INT),    intent(in), value :: i
  integer(kind=C_INT),    intent(in), value :: ix
  integer(kind=C_INT),    intent(in), value :: turn_in
  integer(kind=C_INT),    intent(in), value :: particleID_in
  integer(kind=C_INT),    intent(in), value :: parentID_in
  integer(kind=C_INT),    intent(in), value :: pdgID_in
  integer(kind=C_INT16_T),intent(in), value :: q_in
  real(kind=C_DOUBLE),    intent(in), value :: weight_in
  real(kind=C_DOUBLE),    intent(in), value :: s_in
  real(kind=C_DOUBLE),    intent(in), value :: x_in
  real(kind=C_DOUBLE),    intent(in), value :: xp_in
  real(kind=C_DOUBLE),    intent(in), value :: y_in
  real(kind=C_DOUBLE),    intent(in), value :: yp_in
  real(kind=C_DOUBLE),    intent(in), value :: z_in
  real(kind=C_DOUBLE),    intent(in), value :: dp_in
  real(kind=C_DOUBLE),    intent(in), value :: sx_in
  real(kind=C_DOUBLE),    intent(in), value :: sy_in
  real(kind=C_DOUBLE),    intent(in), value :: sz_in
  real(kind=C_DOUBLE),    intent(in), value :: m_in
end subroutine
end interface

contains

subroutine SixTrackRootInit
  implicit none
  if(root_flag)  then
    call DoSixTrackRootInit(root_eos_enabled, root_RunNumber, root_eos_server, root_folder, root_prefix, root_Accelerator, &
&                           root_Optics, root_ApertureCheck, root_Collimation, root_CollimationDB, root_FLUKA, root_DumpPipe)
  end if
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
  root_FLUKA         = 0
  root_DumpPipe      = 0
  root_RunNumber     = 0
  root_eos_server    = C_NULL_CHAR
  root_folder        = C_NULL_CHAR
  root_prefix        = C_NULL_CHAR
end subroutine SixTrackRootFortranInit

subroutine root_daten(inLine,iErr)

  use string_tools

  implicit none

  character(len=*), intent(in)    :: inLine
  logical,          intent(inout) :: iErr

  character(len=:), allocatable :: lnSplit(:)
  integer nSplit
  logical spErr, cErr

  call chr_split(inLine,lnSplit,nSplit,spErr)
  if(spErr) then
    write(lerr,"(a)") "ROOT> ERROR Failed to parse input line."
    iErr = .true.
    return
  end if

  !ROOT is enabled
  root_flag = .true.

  if(nSplit /= 2) then
    write(lerr,"(a,i0)") "ROOT> ERROR Expected 2 entries per line, got ",nSplit
    iErr = .true.
    return
  end if

!  For input debugging if needed
!  write(lout,*) '1: ', getfields_fields(1)(1:getfields_lfields(1))
!  write(lout,*) '2: ', getfields_fields(2)(1:getfields_lfields(2))

!EOS: e.g. eosuser.cern.ch/
  if(lnSplit(1) == 'EOS') then
    root_eos_enabled = 1
    root_eos_server = trim(lnSplit(2)) // C_NULL_CHAR

!PATH e.g. /eos/user/u/username/
  else if(lnSplit(1) == 'PATH') then
    root_folder = trim(lnSplit(2)) // C_NULL_CHAR

!PREFIX sixtrack_
  else if(lnSplit(1) == 'PREFIX') then
    root_prefix = trim(lnSplit(2)) // C_NULL_CHAR

!RUN number
  else if(lnSplit(1) == 'RUN') then
    call chr_cast(lnSplit(2),root_RunNumber,cErr)
!blocks to enable
!ENABLE
!COLL, APER, ALL
  else if(lnSplit(1) == 'ENABLE') then
    write(lout,"(a)") "ROOT> "//trim(lnSplit(2))
    if(lnSplit(2) == 'ALL') then
      root_ApertureCheck = 1
      root_Accelerator = 1
      root_Collimation = 1
      root_CollimationDB = 1
      root_Optics = 1
      root_DumpPipe = 1
#ifdef FLUKA
      root_FLUKA = 1
#endif
    else if(lnSplit(2) == 'ACCEL') then
      root_Accelerator = 1
    else if(lnSplit(2) == 'COLL') then
      root_Collimation = 1
    else if(lnSplit(2) == 'COLDB') then
      root_CollimationDB = 1
    else if(lnSplit(2) == 'APER') then
      root_ApertureCheck = 1
    else if(lnSplit(2) == 'OPTICS') then
      root_Optics = 1
    else if(lnSplit(2) == 'FLUKA') then
      root_FLUKA = 1
    else if(lnSplit(2) == 'PIPE') then
      root_DumpPipe = 1
    end if
  end if

!Check configuration

!check + enable flags

end subroutine root_daten

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
  write(lout,*) 'FLUKA:         ', root_FLUKA
  write(lout,*) 'PIPE:          ', root_DumpPipe

end subroutine root_parseInputDone

end module root_output
