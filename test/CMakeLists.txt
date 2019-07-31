
##
##  SixTrack Test Suite
## =====================
##

##
## SelfTests
##

add_test(NAME SelfTest_verify10 COMMAND verify10)
set_tests_properties(SelfTest_verify10 PROPERTIES PASS_REGULAR_EXPRESSION "Error in verify10 - file 'fort.2[0|1]' was not found\n.*[1|2]")

add_test(NAME SelfTest_checkf10 COMMAND checkf10)
set_tests_properties(SelfTest_checkf10 PROPERTIES PASS_REGULAR_EXPRESSION "Error in checkf10 - file 'fort.2[0|1]' was not found\n.*[1|2]")

add_test(NAME SelfTest_compf10 COMMAND compf10)
set_tests_properties(SelfTest_compf10 PROPERTIES PASS_REGULAR_EXPRESSION "Error in compf10 - file 'fort.2[0|1]' was not found\n.*[1|2]")

add_test(NAME SelfTest_checkf110 COMMAND checkf110)
set_tests_properties(SelfTest_checkf110 PROPERTIES PASS_REGULAR_EXPRESSION "Error in checkf110 - file 'fort.2[0|1]' was not found\n.*[1|2]")

add_test(NAME SelfTest_checkf1014 COMMAND checkf1014)
set_tests_properties(SelfTest_checkf1014 PROPERTIES PASS_REGULAR_EXPRESSION "Error in checkf1014 - file 'fort.2[0|1]' was not found\n.*[1|2]")

add_test(NAME SelfTest_read90 COMMAND read90)
set_tests_properties(SelfTest_read90 PROPERTIES PASS_REGULAR_EXPRESSION "Error in read90 - file 'fort.190' was not found; try --help?\n.*19")

add_test(NAME SelfTest_readDump3 COMMAND readDump3)
set_tests_properties(SelfTest_readDump3 PROPERTIES PASS_REGULAR_EXPRESSION "Expected 2 arguments: inputfile outputfile\n.*3")

if(ZLIB)
  add_test(NAME SelfTest_zlib COMMAND ${CMAKE_BINARY_DIR}/zlib/example)
  if(64BIT)
    add_test(NAME SelfTest_zlib_64 COMMAND ${CMAKE_BINARY_DIR}/zlib/example64)
  endif()
endif()

##
## Checks
##

add_test(NAME CheckTestInput_STF WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test/ COMMAND ${CMAKE_SOURCE_DIR}/test/CheckTestInputForSTF.sh $<TARGET_FILE:read90>)
set_tests_properties(CheckTestInput_STF PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR")

add_test(NAME CheckBuildManual COMMAND ${CMAKE_SOURCE_DIR}/test/CheckBuildManual.sh "${CMAKE_SOURCE_DIR}")
set_tests_properties(CheckBuildManual PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR")

##
## ERROR Tests
## Only run for non-CR version as CR version doesn't write directly to stderr
##

if(NOT CR)
  list(APPEND SIXTRACK_ERROR
    error_init
    error_trac
  )
endif(NOT CR)

##
## Needs checking: thick6dsingles (broken input files -- this test could probably be removed)
## thick6dblocks elensidealthck6d dynk_globalvars thick6ddynk - empty fort.10 (no POST block in fort.3) -- this is OK, we still test the tracking with fort.90 / STF.
##

##
# List of Tests
##

list(APPEND SIXTRACK_TESTS
  bb
  bb_ntwin1
  bbe51
  bbe52
  bbe571ib0
  beam-HO_4D-newstyle
  beam-HO_4D-newstyle-loadParticles
  beam-HO_4D-oldstyle
  beam-HO_6D-newstyle
  beam-HO_6D-oldstyle
  beam-HO_6D-simple-newstyle
  beam-HO_6D-simple-newstyle-orbit
  beam-HO_6D-simple-oldstyle
  beam-HO_LR-ibbc-newstyle
  beam-HO_LR-ibbc-oldstyle
  beam-HO_LR-newstyle
  beam-HO_LR-oldstyle
  beambeam4dDYNK
  beambeamDYNK
  chebythin6d
  chebythin6d_DYNK
  chebythin6d_ions
  chebythin6d_ramp_DYNK
  crabs
  dipedge
  dist_file_6d
  dist_file_6d_direct
  distance
  dump_2_all_ele_highprec
  dump_all_highprec
  dump_binary
  dump4
  dump5
  dump6
  dump7_fma
  dump7_fma_thin4d
  dump8_fma
  dump9
  dynk_filelin
  dynk_globalvars
  elensidealthck4d
  elensidealthck6d
  elensidealthin4d
  elensidealthin6d
  elensidealthin6d_DYNK
  elensidealthin6d_ions
  elensidealthin6d_ramp_DYNK
  eric
  exact
  fcc
  fma
  fma_binary
  fringe_field
  frs
  frs60
  generalized-rf-multipoles
  generalized-rf-multipoles-missalignments
  geometry
  geometry_multicol
  last_amp
  lost
  lostevery
  lostnotilt
  lostnumxv
  multipoles_DYNK
  notilt
  orbit6d-element-combinedFunction
  orbit6d-element-quadrupole
  orbit6d-element-solenoid
  orbit6d-element-trombone
  orbit6d-element-trombone-kick
  orbit6d-element-rotations
  orbit6d-ions
  orbit6d-ions-long
  part_state_binary_ions
  part_state_text_ions
  perf1
  prob1
  prob3
  rfMultipoles_2nd3rd4th_skew
  rfMultipoles_2ndOrder
  rfMultipoles_3rdOrder
  rfMultipoles_4thOrder
  s316
  scatter_bbelastic
  thick4
  thick6d_below_transit
  thick6dblocks
  thick6ddynk
  thin4d_ions
  thin6d_below_transit
  thin6d_ions
  tilt
  wire_thck4d
  wire_thck6d
  wire_thin4d
  wire_thin6d
  wire_thin6d_ibeco0
)

##
# Tests depending on flags
##

if(STF)
  list(APPEND SIXTRACK_TESTS
    beam-HO_6D-simple-newstyle-coupling
    dynk_IIR
    dynk_operators
    dynk_polynomial
    dynk_random
    scatter_aperture
  )
endif()

if(NOT CR)
  list(APPEND SIXTRACK_TESTS
    collimation_dist_1h
    collimation_dist_1v
    collimation_dist_2h
    collimation_dist_2v
    collimation_dist_3h
    collimation_dist_3v
    collimation_dist_4
    collimation_dist_5
    collimation_dist_6
    collimation_dist_radial
    collimation_k2
    collimation_new-db_new-block
    collimation_old-db_new-block
    collimation_old-db_old-block
    scatter_collimation
  )
endif()

if(LIBARCHIVE)
  list(APPEND SIXTRACK_TESTS
    elensidealthin6d_DYNK_ZIPF
  )
endif()

if(NAFF)
  list(APPEND SIXTRACK_TESTS
    fma_naff
  )
endif()

if(PYTHIA)
  list(APPEND SIXTRACK_TESTS
    scatter_pythia
  )
endif()

if(DISTLIB)
  list(APPEND SIXTRACK_TESTS
    distlib_normalised
    distlib_tmatrix
    distlib_twiss
  )
endif()

if(SIXDA)
  list(APPEND SIXTRACK_DA_TESTS
    frs_da
    thick4_da
  )
endif()

## UNUSED TESTS

# orbit6d-element-wire
# thick6dsingles (only contains one canonical, fort.6, which shows that it crashes. Not sure what it is supposed to show?)
# thick4in_da (copy of thick4_da that contains an extra huge MULT block in fort.3?)

## END UNUSED TESTS

##
# Test Configurations
##

# List of tests that do not produce a valid fort.10
set(TESTS_NOFORT10
  beam-HO_6D-simple-newstyle
  beam-HO_6D-simple-newstyle-coupling
  beam-HO_6D-simple-newstyle-orbit
  beam-HO_6D-simple-oldstyle
  beambeam4dDYNK
  beambeamDYNK
  collimation_dist_1h
  collimation_dist_1v
  collimation_dist_2h
  collimation_dist_2v
  collimation_dist_3h
  collimation_dist_3v
  collimation_dist_4
  collimation_dist_5
  collimation_dist_6
  collimation_dist_radial
  collimation_k2
  collimation_new-db_new-block
  collimation_old-db_new-block
  collimation_old-db_old-block
  dist_file_6d
  dist_file_6d_direct
  distlib_normalised
  distlib_tmatrix
  distlib_twiss
  dump_2_all_ele_highprec
  dump_all_highprec
  dump_binary
  dump4
  dump5
  dump6
  dynk_filelin
  dynk_globalvars
  dynk_IIR
  dynk_operators
  dynk_polynomial
  dynk_random
  elensidealthck6d
  exact
  fringe_field
  generalized-rf-multipoles
  generalized-rf-multipoles-missalignments
  geometry
  geometry_multicol
  multipoles_DYNK
  orbit6d-element-combinedFunction
  orbit6d-element-quadrupole
  orbit6d-element-solenoid
  orbit6d-element-trombone
  orbit6d-element-trombone-kick
  orbit6d-element-rotations
  orbit6d-ions
  orbit6d-ions-long
  part_state_binary_ions
  part_state_text_ions
  rfMultipoles_2nd3rd4th_skew
  rfMultipoles_2ndOrder
  rfMultipoles_3rdOrder
  rfMultipoles_4thOrder
  scatter_aperture
  scatter_bbelastic
  scatter_collimation
  scatter_pythia
  thick6d_below_transit
  thick6dblocks
  thick6ddynk
  thin6d_below_transit
)

# List of tests that do not produce singletrackfile.dat or fort.90
set(TESTS_NO_STF_FORT90
  beam-HO_6D-simple-newstyle-coupling
  chebythin6d
  chebythin6d_DYNK
  chebythin6d_ions
  chebythin6d_ramp_DYNK
  collimation_dist_1h
  collimation_dist_1v
  collimation_dist_2h
  collimation_dist_2v
  collimation_dist_3h
  collimation_dist_3v
  collimation_dist_4
  collimation_dist_5
  collimation_dist_6
  collimation_dist_radial
  collimation_k2
  collimation_new-db_new-block
  collimation_old-db_new-block
  collimation_old-db_old-block
  crabs
  orbit6d-element-rotations
  dipedge
  dist_file_6d
  dist_file_6d_direct
  distlib_normalised
  distlib_tmatrix
  distlib_twiss
  dump_all_highprec
  dump4
  dump5
  dump6
  dynk_IIR
  dynk_operators
  dynk_polynomial
  dynk_random
  elensidealthin6d_ions
  fma
  fma_binary
  fma_naff
  fringe_field
  geometry
  geometry_multicol
  part_state_binary_ions
  part_state_text_ions
  rfMultipoles_2ndOrder
  rfMultipoles_4thOrder
  scatter_aperture
  scatter_bbelastic
  scatter_collimation
  scatter_pythia
  thick6d_below_transit
  thick6ddynk
  thin4d_ions
  thin6d_below_transit
  thin6d_ions
)

# For CR builds, test the CR functions, so kill/resume the test
if(CR)
  set(CR_ENABLE 1)
  # Should at least be able to reach the standard TIMEOUT.
  set(CR_KILL_TIME "2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,16384,16384,16384,16384,16384,16384,16384,16384,16384")
  list(APPEND CR_MUST_RESTART
    chebythin6d_ions
    chebythin6d_ramp_DYNK
    dipedge
    dynk_filelin
    dynk_globalvars
    dynk_IIR
    dynk_operators
    dynk_polynomial
    dynk_random
    elensidealthck4d
    elensidealthck6d
    elensidealthin6d_ramp_DYNK
    scatter_aperture
    scatter_bbelastic
  )
else()
  set(CR_ENABLE 0)
  set(CR_KILL_TIME 0)
  set(CR_MUST_RESTART "")
endif(CR)

# Set flags for specific tests
foreach(TST IN ITEMS ${SIXTRACK_TESTS})

  # Check fort.10 by default
  set(FORT10 1)
  if(TST IN_LIST TESTS_NOFORT10)
    set(FORT10 0)
  else()
    set(FORT10 1)
  endif()

  # if STF, we need to use the STF flag on read90, and disable checking fort.90
  if(STF)
    set(STF_ENABLE 1)
    set(FORT90 0)
  else()
    set(STF_ENABLE 0)
    set(FORT90 1)
  endif(STF)

  if(TST IN_LIST TESTS_NO_STF_FORT90)
    set(STF_ENABLE 0)
    set(FORT90 0)
  endif()

  if(LIBARCHIVE)
    if(TST STREQUAL elensidealthin6d_DYNK_ZIPF)
      set(SIXOUT_NUMFILES 2)
    else()
      set(SIXOUT_NUMFILES 0)
    endif()
  else()
    set(SIXOUT_NUMFILES 0)
  endif(LIBARCHIVE)

  set(DUMP_BIN_FILE NONE)
  if(TST STREQUAL dump_binary)
    set(DUMP_BIN_FILE dump.dat,dumpStart.dat)
  elseif(TST STREQUAL fma_binary)
    set(DUMP_BIN_FILE START_DUMP_1,IP2_DUMP_1,IP2_DUMP_2,IP2_DUMP_3,IP3_DUMP_1,IP3_DUMP_2,IP3_DUMP_3)
  elseif(TST STREQUAL dump7_fma)
    set(DUMP_BIN_FILE IP3_DUMP_1)
  elseif(TST STREQUAL dump8_fma)
    set(DUMP_BIN_FILE START_DUMP_1,IP3_DUMP_1)
  elseif(TST STREQUAL dump7_fma_thin4d)
    set(DUMP_BIN_FILE dump_ip3)
  elseif(TST STREQUAL dump_all_highprec)
    set(DUMP_BIN_FILE dump3.dat,dump8.dat)
  endif()

  set(CR_NUMRESTARTS 0)
  if(TST IN_LIST CR_MUST_RESTART)
    set(CR_NUMRESTARTS 1)
  endif()

  add_test(
    NAME ${TST} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test/${TST}
    COMMAND SixTestWrapper $<TARGET_FILE:SixTrackCR> $<TARGET_FILE:checkf10> $<TARGET_FILE:read90>
      ${FORT10} ${FORT90} ${STF_ENABLE} ${SIXOUT_NUMFILES} ${CR_ENABLE} ${CR_KILL_TIME}
      $<TARGET_FILE:readDump3> ${DUMP_BIN_FILE} ${CR_NUMRESTARTS}
  )

  #some tests take a very long time (set a 2 day timeout)
  set_tests_properties(${TST} PROPERTIES TIMEOUT "172800")

endforeach()

# DA Tests
if(SIXDA)
  foreach(TST IN ITEMS ${SIXTRACK_DA_TESTS})

    set(FORT10 0)
    set(FORT90 0)
    set(STF_ENABLE 0)
    set(SIXOUT_NUMFILES 0)
    set(CR_ENABLE 0)
    set(CR_KILL_TIME 0)
    set(DUMP_BIN_FILE NONE)

    add_test(
      NAME ${TST} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test/${TST}
      COMMAND SixTestWrapper $<TARGET_FILE:SixTrackDA> $<TARGET_FILE:checkf10> $<TARGET_FILE:read90>
        ${FORT10} ${FORT90} ${STF_ENABLE} ${SIXOUT_NUMFILES} ${CR_ENABLE} ${CR_KILL_TIME}
        $<TARGET_FILE:readDump3> ${DUMP_BIN_FILE} 0
    )

    #some tests take a very long time (set a 2 day timeout)
    set_tests_properties(${TST} PROPERTIES TIMEOUT "172800")

  endforeach()
endif()

# Set up tests
foreach(TST IN ITEMS ${SIXTRACK_TESTS} ${SIXTRACK_DA_TESTS} ${SIXTRACK_ERROR})
  add_custom_target(test_copy_${TST})
  add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/test/${TST}/)

  file(GLOB TEST_INPUTS LIST_DIRECTORIES FALSE ${CMAKE_SOURCE_DIR}/test/${TST}/*)
  #message(STATUS "${CMAKE_SOURCE_DIR}/../test/${TST}/")

  if(NOT BOINC)
    foreach(TSFILE IN ITEMS ${TEST_INPUTS})
      get_filename_component(TFILE ${TSFILE} NAME)
      if(WIN32 OR NOSYMLINK)
        add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${TSFILE} ${CMAKE_BINARY_DIR}/test/${TST}/${TFILE})
      else()
        add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E create_symlink ${TSFILE} ${CMAKE_BINARY_DIR}/test/${TST}/${TFILE})
      endif()
    endforeach()
  else()
    # Create Sixin.zip files
    set(TST_TO_ZIP fort.3) # Required
    foreach(ZIP_INPUT IN ITEMS fort.2 fort.8 fort.13 fort.16) # Optional
      if(EXISTS ${CMAKE_SOURCE_DIR}/test/${TST}/${ZIP_INPUT})
        set(TST_TO_ZIP ${TST_TO_ZIP} ${ZIP_INPUT})
      endif()
    endforeach()
    if(EXISTS ${CMAKE_SOURCE_DIR}/test/${TST}/extra_inputs.txt)
      add_custom_command(
        TARGET test_copy_${TST} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E tar "cf" "${CMAKE_BINARY_DIR}/test/${TST}/Sixin.zip" --format=zip --files-from=${CMAKE_SOURCE_DIR}/test/${TST}/extra_inputs.txt -- ${TST_TO_ZIP}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test/${TST}/
      )
    else()
      add_custom_command(
        TARGET test_copy_${TST} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E tar "cf" "${CMAKE_BINARY_DIR}/test/${TST}/Sixin.zip" --format=zip -- ${TST_TO_ZIP}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test/${TST}/
      )
    endif()

    # Copy/link the canonicals
    file(GLOB TEST_INPUT_CANONICAL LIST_DIRECTORIES FALSE ${CMAKE_SOURCE_DIR}/test/${TST}/*canonical)
    foreach(CANONICAL_FILES IN ITEMS ${TEST_INPUT_CANONICAL})
      get_filename_component(CANONICAL_FILE ${CANONICAL_FILES} NAME)
      if(WIN32 OR NOSYMLINK)
        add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${CANONICAL_FILES} ${CMAKE_BINARY_DIR}/test/${TST}/${CANONICAL_FILE})
      else()
        add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E create_symlink ${CANONICAL_FILES} ${CMAKE_BINARY_DIR}/test/${TST}/${CANONICAL_FILE})
      endif()
    endforeach()

    #And extra_checks.txt
    file(GLOB TEST_INPUT_EXTRAS LIST_DIRECTORIES FALSE ${CMAKE_SOURCE_DIR}/test/${TST}/extra_checks.txt)
    foreach(EXTRA_FILES IN ITEMS ${TEST_INPUT_EXTRAS})
      get_filename_component(EXTRA_FILE ${EXTRA_FILES} NAME)
      if(WIN32 OR NOSYMLINK)
        add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${EXTRA_FILES} ${CMAKE_BINARY_DIR}/test/${TST}/${EXTRA_FILE})
      else()
        add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E create_symlink ${EXTRA_FILES} ${CMAKE_BINARY_DIR}/test/${TST}/${EXTRA_FILE})
      endif()
    endforeach()

  endif(NOT BOINC)

  add_dependencies(SixTrackCore test_copy_${TST})
endforeach()

# Set up error tests
foreach(TST IN ITEMS ${SIXTRACK_ERROR})
  add_test(
    NAME ${TST} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test/${TST}
    COMMAND "errorTests.py" $<TARGET_FILE:SixTrackCR>
  )
  set_tests_properties(${TST} PROPERTIES LABELS error)
  if(WIN32 OR NOSYMLINK)
    add_custom_command(
      TARGET test_copy_${TST} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/test/tools/errorTestTools.py ${CMAKE_BINARY_DIR}/test/${TST}/errorTestTools.py
    )
  else()
    add_custom_command(
      TARGET test_copy_${TST} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/test/tools/errorTestTools.py ${CMAKE_BINARY_DIR}/test/${TST}/errorTestTools.py
    )
  endif()
endforeach()

#add labels to tests
#Lets say fast tests take <60 seconds to run
set(FAST_TESTS
  bb
  bb_ntwin1
  beam-HO_4D-newstyle
  beam-HO_4D-newstyle-loadParticles
  beam-HO_4D-oldstyle
  beam-HO_6D-newstyle
  beam-HO_6D-oldstyle
  beam-HO_6D-simple-newstyle
  beam-HO_6D-simple-newstyle-coupling
  beam-HO_6D-simple-newstyle-orbit
  beam-HO_6D-simple-oldstyle
  beam-HO_LR-ibbc-newstyle
  beam-HO_LR-ibbc-oldstyle
  beam-HO_LR-newstyle
  beam-HO_LR-oldstyle
  beambeam4dDYNK
  beambeamDYNK
  chebythin6d
  chebythin6d_DYNK
  chebythin6d_ions
  chebythin6d_ramp_DYNK
  collimation_dist_1h
  collimation_dist_1v
  collimation_dist_2h
  collimation_dist_2v
  collimation_dist_3h
  collimation_dist_3v
  collimation_dist_4
  collimation_dist_5
  collimation_dist_6
  collimation_dist_radial
  collimation_old-db_new-block
  collimation_old-db_old-block
  collimation_old-db_old-block
  crabs
  dist_file_6d
  dist_file_6d_direct
  distlib_normalised
  distlib_tmatrix
  distlib_twiss
  dump_2_all_ele_highprec
  dump_all_highprec
  dump_binary
  dump4
  dump5
  dump6
  dump7_fma
  dump7_fma_thin4d
  dump8_fma
  dump9
  dynk_filelin
  dynk_globalvars
  dynk_IIR
  dynk_operators
  dynk_polynomial
  dynk_random
  elensidealthck4d
  elensidealthck6d
  elensidealthin4d
  elensidealthin6d
  elensidealthin6d_DYNK
  elensidealthin6d_DYNK_ZIPF
  elensidealthin6d_ions
  exact
  fcc
  frs_da
  generalized-rf-multipoles
  generalized-rf-multipoles-missalignments
  geometry
  geometry_multicol
  last_amp
  lost
  lostnotilt
  lostnumxv
  multipoles_DYNK
  notilt
  orbit6d-element-combinedFunction
  orbit6d-element-quadrupole
  orbit6d-element-solenoid
  orbit6d-element-trombone
  orbit6d-element-trombone-kick
  orbit6d-element-rotations
  orbit6d-ions
  orbit6d-ions-long
  part_state_binary_ions
  part_state_text_ions
  rfMultipoles_2nd3rd4th_skew
  rfMultipoles_2ndOrder
  rfMultipoles_3rdOrder
  rfMultipoles_4thOrder
  s316
  scatter_aperture
  scatter_bbelastic
  scatter_collimation
  scatter_pythia
  thick4
  thick4_da
  thick6d_below_transit
  thick6dblocks
  thick6ddynk
  thin4d_ions
  thin6d_below_transit
  thin6d_ions
  wire_thck4d
  wire_thck6d
  wire_thin4d
  wire_thin6d
  wire_thin6d_ibeco0
)

# Medium tests take < 30 mins to run
set(MEDIUM_TESTS
  bbe51
  bbe52
  bbe571ib0
  chebythin6d_ramp_DYNK
  collimation_k2
  dipedge
  distance
  elensidealthin6d_ramp_DYNK
  eric
  fma
  fma_binary
  fma_naff
  fringe_field
  frs
  frs60
  lostevery
  tilt
)

# Slow is everything longer
set(SLOW_TESTS
  prob1
  prob3
)

# Performance test
set(PERF_TESTS
  perf1
)

foreach(TST IN ITEMS ${SIXTRACK_TESTS} ${SIXTRACK_DA_TESTS})

  list(FIND FAST_TESTS   ${TST} LABEL_FAST)
  list(FIND MEDIUM_TESTS ${TST} LABEL_MEDIUM)
  list(FIND SLOW_TESTS   ${TST} LABEL_SLOW)
  list(FIND PERF_TESTS   ${TST} LABEL_PERF)

  if(${LABEL_FAST} GREATER -1)
    set_property(TEST ${TST} PROPERTY LABELS fast)
  elseif(${LABEL_MEDIUM} GREATER -1)
    set_property(TEST ${TST} PROPERTY LABELS medium)
  elseif(${LABEL_SLOW} GREATER -1)
    set_property(TEST ${TST} PROPERTY LABELS slow)
  elseif(${LABEL_PERF} GREATER -1)
    set_property(TEST ${TST} PROPERTY LABELS perf)
  else()
    message(STATUS "TEST ${TST} needs a label setting!")
    set_property(TEST ${TST} PROPERTY LABELS NeedsLabel)
  endif()
endforeach()

# Tools
set(SELFTEST_TOOLS
  SelfTest_verify10
  SelfTest_compf10
  SelfTest_checkf1014
  SelfTest_checkf10
  SelfTest_checkf110
  SelfTest_read90
  SelfTest_readDump3
)
foreach(TOOL IN ITEMS ${SELFTEST_TOOLS})
  set_property(TEST ${TOOL} PROPERTY LABELS tool fast)
endforeach()

# Zlib tests
if(ZLIB)
  set_property(TEST SelfTest_zlib PROPERTY LABELS tool fast)
  if(64BIT)
    set_property(TEST SelfTest_zlib_64 PROPERTY LABELS tool fast)
  endif()
endif()

# Inputs
set(SELFTEST_INPUT
  CheckTestInput_STF
)
foreach(TOOL IN ITEMS ${SELFTEST_INPUT})
  set_property(TEST ${TOOL} PROPERTY LABELS input fast)
endforeach()

# Other
set(SELFTEST_OTHER
  CheckBuildManual
)
foreach(TOOL IN ITEMS ${SELFTEST_OTHER})
  set_property(TEST ${TOOL} PROPERTY LABELS input fast)
endforeach()
