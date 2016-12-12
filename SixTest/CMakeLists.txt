#Test that the testing tools run!
ADD_TEST(NAME SelfTest_verify10 COMMAND verify10)
set_tests_properties(SelfTest_verify10 PROPERTIES PASS_REGULAR_EXPRESSION "Comparing VERSION    0.0000000000000000       to    0.0000000000000000.*\n.*verify10_DIFF I\\/O error, wrong no of lines\\!\\!\\?\\? line no            0")

ADD_TEST(NAME SelfTest_checkf10 COMMAND checkf10)
set_tests_properties(SelfTest_checkf10 PROPERTIES PASS_REGULAR_EXPRESSION "Comparing VERSION    0.0000000000000000       to    0.0000000000000000.*\n.*checkf10_DIFF I\\/O error, wrong no of lines\\!\\! line no            0")

ADD_TEST(NAME SelfTest_compf10 COMMAND compf10)
set_tests_properties(SelfTest_compf10 PROPERTIES PASS_REGULAR_EXPRESSION "Comparing VERSION    0.0000000000000000       to    0.0000000000000000.*\n.*compf10_DIFF I\\/O error, wrong no of lines\\!\\! line no            0")

ADD_TEST(NAME SelfTest_checkf110 COMMAND checkf110)
set_tests_properties(SelfTest_checkf110 PROPERTIES PASS_REGULAR_EXPRESSION "Comparing VERSION    0.0000000000000000       to    0.0000000000000000.*\n.*DIFF I\\/O error, wrong no of lines\\!\\! line no            0")

ADD_TEST(NAME SelfTest_checkf1014 COMMAND checkf1014)
set_tests_properties(SelfTest_checkf1014 PROPERTIES PASS_REGULAR_EXPRESSION "Comparing VERSION    0.0000000000000000       to    0.0000000000000000.*\n.*checkf1014_DIFF I\\/O error, wrong no of lines\\!\\! line no            0")

ADD_TEST(NAME SelfTest_read90 COMMAND read90)
set_tests_properties(SelfTest_read90 PROPERTIES PASS_REGULAR_EXPRESSION "At line .* of file .*read90.f90 \\(unit = 190\\)\nFortran runtime error: Cannot open file 'fort.190': No such file or directory")

SET(SIXTRACK_TESTS ${SIXTRACK_TESTS})
#SET(SIXTRACK_TESTS bbe51 bbe52 bbe571ib0 bb bb_ntwin1 crabamp dipedge distance dynk_globalvars elensidealthck4d elensidealthck6d elensidealthin4d elensidealthin6d_DYNK elensidealthin6d eric exact fma frs60 frs javier lostevery lost lostnotilt notilt prob1 prob3 s316 thick4 thick6dblocks thick6ddynk tilt ${SIXTRACK_TESTS})

##
## Needs checking: lostnumxv thick6dsingles (missing files)
## thick6dblocks elensidealthck6d dynk_globalvars thick6ddynk empty fort.10
## FIX: javier_bignpart (and STF)
## Disabled: elensidealthin6d_DYNK
##

SET(STF_TESTS bbe51 bbe52 bbe571ib0 bb_ntwin1 bb crabamp dipedge distance elensidealthck4d elensidealthin4d elensidealthin6d eric exact fma frs60 frs javier_bignpart javier lostevery lostnotilt lost notilt prob1 prob3 s316 thick4 tilt)

#only turn on the bnl test if we have BNLELENS enabled due to slight differences in the test
if(BNLELENS)
	SET(SIXTRACK_TESTS bnl ${SIXTRACK_TESTS})
endif(BNLELENS)

if(NOT COLLIMAT)
	SET(SIXTRACK_TESTS bbe51 bbe52 bbe571ib0 bb bb_ntwin1 crabamp dipedge distance dynk_globalvars elensidealthin4d elensidealthin6d eric exact frs60 frs javier lostevery lost lostnotilt notilt prob1 prob3 s316 thick4 thick6dblocks tilt elensidealthck4d elensidealthck6d fma thick6ddynk ${SIXTRACK_TESTS})
endif(NOT COLLIMAT)

#For CR builds, test the CR functions, so kill/resume the test
if(CR)
	SET(CR_ENABLE 1)
	SET(CR_KILL_TIME 10)
else()
	SET(CR_ENABLE 0)
	SET(CR_KILL_TIME 0)
endif(CR)

#Test commands
foreach(TST IN ITEMS ${SIXTRACK_TESTS})

	#determine what we need to do:
	#Check fort.10 by default
	SET(FORT10 1)
	#no fort.10 for the "exact" test
	if(TST MATCHES exact OR TST MATCHES thick6dblocks OR TST MATCHES elensidealthck6d OR TST MATCHES dynk_globalvars OR TST MATCHES thick6ddynk)
		SET(FORT10 0)
	endif()

	#check fort.90 by default
	SET(FORT90 1)
	if(TST MATCHES lostnumxv )
		SET(FORT90 0)
	endif()

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

	ADD_TEST(NAME ${TST} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/SixTest/${TST} COMMAND SixTestWrapper $<TARGET_FILE:SixTrack> ${FORT10} ${FORT90} ${STF_ENABLE} ${CR_ENABLE} ${CR_KILL_TIME})

	#some tests take a very long time
	if(TST MATCHES prob1 OR TST MATCHES prob3)
		set_tests_properties(${TST} PROPERTIES TIMEOUT "1000000")
	endif()

endforeach()

#test files
foreach(TST IN ITEMS ${SIXTRACK_TESTS})
	add_custom_target(test_copy_${TST})
	add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/SixTest/${TST}/)

	#link in the test tools
	add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE:read90> ${CMAKE_BINARY_DIR}/SixTest/${TST}/read90)
	add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE:checkf10> ${CMAKE_BINARY_DIR}/SixTest/${TST}/checkf10)

	file(GLOB TEST_INPUTS LIST_DIRECTORIES FALSE ${CMAKE_SOURCE_DIR}/SixTest/${TST}/*)
	#message(STATUS "${CMAKE_SOURCE_DIR}/../SixTest/${TST}/")

	foreach(TSFILE IN ITEMS ${TEST_INPUTS})
		get_filename_component(TFILE ${TSFILE} NAME)

		#create_symlink old new
		add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E create_symlink ${TSFILE} ${CMAKE_BINARY_DIR}/SixTest/${TST}/${TFILE})
	endforeach()
	add_dependencies( SixTrack test_copy_${TST})
endforeach()

