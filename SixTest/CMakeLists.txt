#Test that the testing tools run!
ADD_TEST(NAME SelfTest_verify10 COMMAND verify10)
set_tests_properties(SelfTest_verify10 PROPERTIES PASS_REGULAR_EXPRESSION "Error in verify10 - file 'fort.2[0|1]' was not found\n.*[1|2]")

ADD_TEST(NAME SelfTest_checkf10 COMMAND checkf10)
set_tests_properties(SelfTest_checkf10 PROPERTIES PASS_REGULAR_EXPRESSION "Error in checkf10 - file 'fort.2[0|1]' was not found\n.*[1|2]")

ADD_TEST(NAME SelfTest_compf10 COMMAND compf10)
set_tests_properties(SelfTest_compf10 PROPERTIES PASS_REGULAR_EXPRESSION "Error in compf10 - file 'fort.2[0|1]' was not found\n.*[1|2]")

ADD_TEST(NAME SelfTest_checkf110 COMMAND checkf110)
set_tests_properties(SelfTest_checkf110 PROPERTIES PASS_REGULAR_EXPRESSION "Error in checkf110 - file 'fort.2[0|1]' was not found\n.*[1|2]")

ADD_TEST(NAME SelfTest_checkf1014 COMMAND checkf1014)
set_tests_properties(SelfTest_checkf1014 PROPERTIES PASS_REGULAR_EXPRESSION "Error in checkf1014 - file 'fort.2[0|1]' was not found\n.*[1|2]")

ADD_TEST(NAME SelfTest_read90 COMMAND read90)
set_tests_properties(SelfTest_read90 PROPERTIES PASS_REGULAR_EXPRESSION "Error in read90 - file 'fort.190' was not found\n.*19")

ADD_TEST(NAME SelfTest_readDump3 COMMAND readDump3)
set_tests_properties(SelfTest_readDump3 PROPERTIES PASS_REGULAR_EXPRESSION "Expected 2 arguments: inputfile outputfile\n.*3")

if(BOINC)
	ADD_TEST(NAME CheckTestInputForBOINC.sh WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/SixTest/ COMMAND CheckTestInputForBOINC.sh)
	set_tests_properties(CheckTestInputForBOINC.sh PROPERTIES FAIL_REGULAR_EXPRESSION "[^-.*$]|[^diff]")
endif(BOINC)

#SET(SIXTRACK_TESTS bbe51 bbe52 bbe571ib0 bb bb_ntwin1 crabamp dipedge distance dynk_globalvars elensidealthck4d elensidealthck6d elensidealthin4d elensidealthin6d_DYNK elensidealthin6d eric exact fma frs60 frs javier lostevery lost lostnotilt notilt prob1 prob3 s316 thick4 thick6dblocks thick6ddynk tilt ${SIXTRACK_TESTS})

##
## Needs checking: thick6dsingles (broken input files -- this test could probably be removed)
## thick6dblocks elensidealthck6d dynk_globalvars thick6ddynk - empty fort.10 (no POST block in fort.3) -- this is OK, we still test the tracking with fort.90.
## "Totally untested": bnl
##

#List of tests which are STF compatible
SET(STF_TESTS bbe51 bbe52 bbe571ib0 bb_ntwin1 bb crabamp dipedge distance elensidealthck4d elensidealthin4d elensidealthin6d elensidealthck6d elensidealthin6d_DYNK wire_thck4d wire_thck6d wire_thin4d wire_thin6d wire_thin6d_ibeco0 eric exact fma frs60 frs javier_bignpart javier lostevery lostnotilt lost notilt prob1 prob3 s316 thick4 tilt)

if(NOT COLLIMAT)
	if(LIBARCHIVE)
		SET(SIXTRACK_TESTS elensidealthin6d_DYNK_ZIPF ${SIXTRACK_TESTS})
	endif()

	#only turn on the bnl test if we have BNLELENS enabled due to slight differences in the test
	if(BNLELENS)
		SET(SIXTRACK_TESTS bnl ${SIXTRACK_TESTS})
	endif(BNLELENS)
	
	if(BIGNPART OR HUGENPART)
		SET(SIXTRACK_TESTS javier_bignpart ${SIXTRACK_TESTS})
	endif()

	if(BIGNBLZ OR HUGENBLZ)
		SET(SIXTRACK_TESTS fcc ${SIXTRACK_TESTS})
	endif()
	
	SET(SIXTRACK_TESTS bbe51 bbe52 bbe571ib0 bb bb_ntwin1 crabamp dipedge distance dynk_globalvars elensidealthin4d elensidealthin6d elensidealthin6d_DYNK eric exact frs60 frs javier lostevery lost lostnumxv lostnotilt notilt prob1 prob3 s316 thick4 thick6dblocks tilt elensidealthck4d elensidealthck6d fma thick6ddynk wire_thck4d wire_thck6d wire_thin4d wire_thin6d wire_thin6d_ibeco0 dump_binary dump4 dump5 dump6 elensidealthin6d_DYNK_RANDU elensidealthin6d_DYNK_ONOFF beam-HO_4D-newstyle beam-HO_4D-oldstyle ${SIXTRACK_TESTS})
	
endif(NOT COLLIMAT)

#For CR builds, test the CR functions, so kill/resume the test
if(CR)
	SET(CR_ENABLE 1)
	#Should at least be able to reach the standard TIMEOUT.
	SET(CR_KILL_TIME "2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,16384,16384,16384,16384,16384,16384,16384,16384,16384")
else()
	SET(CR_ENABLE 0)
	SET(CR_KILL_TIME 0)
endif(CR)

#Test commands
foreach(TST IN ITEMS ${SIXTRACK_TESTS})

	#determine what we need to do:
	#Check fort.10 by default
	SET(FORT10 1)
	#.. But some tests are missing this file:
	if(TST MATCHES exact OR TST MATCHES thick6dblocks OR TST MATCHES elensidealthck6d OR TST MATCHES dynk_globalvars OR TST MATCHES thick6ddynk OR TST MATCHES dump_binary OR TST MATCHES dump4 OR TST MATCHES dump5 OR TST MATCHES dump6)
		SET(FORT10 0)
	endif()

	#check fort.90 by default
	SET(FORT90 1)
	# if(TST MATCHES lostnumxv )
	# 	SET(FORT90 0)
	# endif()

	#if STF, we need to use the STF flag on read90
	if(STF)
		LIST(FIND STF_TESTS ${TST} IS_STF)
		if(NOT ${IS_STF} MATCHES -1)
			SET(STF_ENABLE 1)
		else()
			SET(STF_ENABLE 0)
		endif()
		##Disable checking fort.90
		SET(FORT90 0)
	else()
		SET(STF_ENABLE 0)
	endif(STF)

	if(LIBARCHIVE)
		if(TST MATCHES elensidealthin6d_DYNK_ZIPF)
			SET(SIXOUT_NUMFILES 6)
		else()
			SET(SIXOUT_NUMFILES 0)
		endif()
	else()
		SET(SIXOUT_NUMFILES 0)
	endif(LIBARCHIVE)

	SET( DUMP_BIN_FILE NONE )
	if(TST MATCHES dump_binary)
		SET( DUMP_BIN_FILE dump.dat )
	endif()
	
	ADD_TEST(NAME ${TST} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/SixTest/${TST} COMMAND SixTestWrapper $<TARGET_FILE:SixTrack> $<TARGET_FILE:checkf10> $<TARGET_FILE:read90> ${FORT10} ${FORT90} ${STF_ENABLE} ${SIXOUT_NUMFILES} ${CR_ENABLE} ${CR_KILL_TIME} $<TARGET_FILE:readDump3> ${DUMP_BIN_FILE})

	#some tests take a very long time (set a 2 day timeout)
	set_tests_properties(${TST} PROPERTIES TIMEOUT "172800")

endforeach()

#test files
foreach(TST IN ITEMS ${SIXTRACK_TESTS})
	add_custom_target(test_copy_${TST})
	add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/SixTest/${TST}/)

	file(GLOB TEST_INPUTS LIST_DIRECTORIES FALSE ${CMAKE_SOURCE_DIR}/SixTest/${TST}/*)
	#message(STATUS "${CMAKE_SOURCE_DIR}/../SixTest/${TST}/")

	if(NOT BOINC)
		foreach(TSFILE IN ITEMS ${TEST_INPUTS})
			get_filename_component(TFILE ${TSFILE} NAME)

			if(NOT WIN32)
				#create_symlink old new
				add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E create_symlink ${TSFILE} ${CMAKE_BINARY_DIR}/SixTest/${TST}/${TFILE})
			else()
				add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${TSFILE} ${CMAKE_BINARY_DIR}/SixTest/${TST}/${TFILE})
			endif()
		endforeach()
	else()
		if(NOT WIN32)
			add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/SixTest/${TST}/Sixin.zip ${CMAKE_BINARY_DIR}/SixTest/${TST}/Sixin.zip)
		else()
			add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/SixTest/${TST}/Sixin.zip ${CMAKE_BINARY_DIR}/SixTest/${TST}/Sixin.zip)
		endif()

		#Also need the canonicals
		file(GLOB TEST_INPUT_CANONICAL LIST_DIRECTORIES FALSE ${CMAKE_SOURCE_DIR}/SixTest/${TST}/*canonical)
		foreach(CANONICAL_FILES IN ITEMS ${TEST_INPUT_CANONICAL})
			get_filename_component(CANONICAL_FILE ${CANONICAL_FILES} NAME)
			if(NOT WIN32)
				add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E create_symlink ${CANONICAL_FILES} ${CMAKE_BINARY_DIR}/SixTest/${TST}/${CANONICAL_FILE})
			else()
				add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${CANONICAL_FILES} ${CMAKE_BINARY_DIR}/SixTest/${TST}/${CANONICAL_FILE})
			endif()
		endforeach()

		#And extra_checks.txt
		file(GLOB TEST_INPUT_EXTRAS LIST_DIRECTORIES FALSE ${CMAKE_SOURCE_DIR}/SixTest/${TST}/extra_checks.txt)
		foreach(EXTRA_FILES IN ITEMS ${TEST_INPUT_EXTRAS})
			get_filename_component(EXTRA_FILE ${EXTRA_FILES} NAME)
			if(NOT WIN32)
				add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E create_symlink ${EXTRA_FILES} ${CMAKE_BINARY_DIR}/SixTest/${TST}/${EXTRA_FILE})
			else()
				add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy ${EXTRA_FILES} ${CMAKE_BINARY_DIR}/SixTest/${TST}/${EXTRA_FILE})
			endif()
		endforeach()

	endif(NOT BOINC)

	add_dependencies( SixTrack test_copy_${TST})
endforeach()

#add labels to tests
#Lets say fast tests take <60 seconds to run
SET(FAST_TESTS crabamp thick4 thick6dblocks thick6ddynk exact javier bb bb_ntwin1 notilt lost lostnotilt lostnumxv s316 dipedge dynk_globalvars elensidealthin4d elensidealthck4d elensidealthck6d elensidealthin6d elensidealthin6d_DYNK elensidealthin6d_DYNK_ZIPF wire_thck4d wire_thck6d wire_thin4d wire_thin6d wire_thin6d_ibeco0 dump_binary dump4 dump5 dump6 fcc elensidealthin6d_DYNK_RANDU elensidealthin6d_DYNK_ONOFF)

#Medium tests take < 30 mins to run
SET(MEDIUM_TESTS bbe571ib0 tilt lostevery fma frs frs60 eric distance bbe51 bbe52 javier_bignpart)

#Slow is everything longer
SET(SLOW_TESTS prob1 prob3)

foreach(TST IN ITEMS ${SIXTRACK_TESTS})
	
	list(FIND FAST_TESTS ${TST} LABEL_FAST)
	list(FIND MEDIUM_TESTS ${TST} LABEL_MEDIUM)
	list(FIND SLOW_TESTS ${TST} LABEL_SLOW)
	
	if(${LABEL_FAST} GREATER -1)
		set_property(TEST ${TST} PROPERTY LABELS fast)
	elseif(${LABEL_MEDIUM} GREATER -1)
		set_property(TEST ${TST} PROPERTY LABELS medium)
	elseif(${LABEL_SLOW} GREATER -1)
		set_property(TEST ${TST} PROPERTY LABELS slow)
	else()
		message(STATUS "TEST ${TST} needs a label setting!")
		set_property(TEST ${TST} PROPERTY LABELS NeedsLabel)
	endif()
endforeach()

#tools
SET(SELFTEST_TOOLS SelfTest_verify10 SelfTest_compf10 SelfTest_checkf1014 SelfTest_checkf10 SelfTest_checkf110 SelfTest_read90 SelfTest_readDump3)
if(BOINC)
	SET(SELFTEST_TOOLS ${SELFTEST_TOOLS} CheckTestInputForBOINC.sh)
endif(BOINC)

foreach(TOOL IN ITEMS ${SELFTEST_TOOLS})
	set_property(TEST ${TOOL} PROPERTY LABELS tool fast)
endforeach()

