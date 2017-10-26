SET(SIXTRACK_TESTS thick4da frs ${SIXTRACK_TESTS})

#Test commands
foreach(TST IN ITEMS ${SIXTRACK_TESTS})

	SET(FORT10 0)
	SET(FORT90 0)
	SET(STF_ENABLE 0)
	SET(SIXOUT_NUMFILES 0)
	SET(CR_ENABLE 0)
	SET(CR_KILL_TIME 0)
	SET(DUMP_BIN_FILE NONE)

	ADD_TEST(NAME ${TST} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/SixTest_da/${TST} COMMAND SixTestWrapper $<TARGET_FILE:SixTrack> $<TARGET_FILE:checkf10> $<TARGET_FILE:read90> ${FORT10} ${FORT90} ${STF_ENABLE} ${SIXOUT_NUMFILES} ${CR_ENABLE} ${CR_KILL_TIME} $<TARGET_FILE:readDump3> ${DUMP_BIN_FILE})

	#some tests take a very long time (set a 2 day timeout)
	set_tests_properties(${TST} PROPERTIES TIMEOUT "172800")
endforeach()

#test files
foreach(TST IN ITEMS ${SIXTRACK_TESTS})
	add_custom_target(test_copy_${TST})
	add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/SixTest_da/${TST}/)

	#link in the test tools
	#add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE:read90> ${CMAKE_BINARY_DIR}/SixTest_da/${TST}/read90)
	#add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE:checkf10> ${CMAKE_BINARY_DIR}/SixTest_da/${TST}/checkf10)

	file(GLOB TEST_INPUTS LIST_DIRECTORIES FALSE ${CMAKE_SOURCE_DIR}/SixTest_da/${TST}/*)

	foreach(TSFILE IN ITEMS ${TEST_INPUTS})
		get_filename_component(TFILE ${TSFILE} NAME)

		#create_symlink old new
		add_custom_command(TARGET test_copy_${TST} POST_BUILD COMMAND ${CMAKE_COMMAND} -E create_symlink ${TSFILE} ${CMAKE_BINARY_DIR}/SixTest_da/${TST}/${TFILE})
	endforeach()
	add_dependencies( SixTrack test_copy_${TST})
endforeach()

