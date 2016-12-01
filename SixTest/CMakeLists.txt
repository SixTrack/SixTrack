#if(NOT COLLIMAT)
	#SET(SIXTRACK_TESTS bnl exact crabamp javier thick4 lost tilt dipedge bb bb_ntwin1 notilt lostnotilt lostevery s316 eric frs frs60 bbe51 bbe52 bbe571ib0 distance prob3 prob1 fma elensidealthck4d elensidealthck6d elensidealthin4d elensidealthin6d)

#only turn on the bnl test if we have BNLELENS enabled due to slight differences in the test
if(BNLELENS)
	SET(SIXTRACK_TESTS bnl ${SIXTRACK_TESTS})
endif(BNLELENS)

	SET(SIXTRACK_TESTS bb bbe51 bbe52 bbe571ib0 bb_ntwin1 crabamp dipedge distance elensidealthck4d elensidealthck6d elensidealthin4d elensidealthin6d elensidealthin6d_DYNK eric exact fma frs frs60 javier javier_bignpart lost lostevery lostnotilt notilt prob1 prob3 s316 stf_compare stf_output thick4 thick6dblocks thick6ddynk thick6dsingles tilt ${SIXTRACK_TESTS})
	foreach(TEST IN ITEMS ${SIXTRACK_TESTS})
		ADD_TEST(NAME ${TEST} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/SixTest/${TEST} COMMAND $<TARGET_FILE:SixTrack> )
		set_tests_properties(${TEST} PROPERTIES FAIL_REGULAR_EXPRESSION ".*RUN TERMINATED ABNORMALLY !!!$")
	endforeach()

	foreach(TST IN ITEMS ${SIXTRACK_TESTS})
		add_custom_target(test_copy_${TST})
		add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/SixTest/${TST}/)

		file(GLOB TEST_INPUTS LIST_DIRECTORIES FALSE ${CMAKE_SOURCE_DIR}/SixTest/${TST}/*)
		#message(STATUS "${CMAKE_SOURCE_DIR}/../SixTest/${TST}/")

		foreach(TSFILE IN ITEMS ${TEST_INPUTS})
			get_filename_component(TFILE ${TSFILE} NAME)
			#message(STATUS "ITEM: ${TFILE}")
			#create_symlink old new
			add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E create_symlink ${TSFILE} ${CMAKE_BINARY_DIR}/SixTest/${TST}/${TFILE})
		endforeach()
		add_dependencies( SixTrack test_copy_${TST})
	endforeach()

	add_dependencies( SixTrack test_copy)
	add_custom_target(test_copy)
	add_custom_command(TARGET test_copy POST_BUILD COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/SixTest/try_pro ${CMAKE_BINARY_DIR}/SixTest/try_pro)

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

#endif(NOT COLLIMAT)
