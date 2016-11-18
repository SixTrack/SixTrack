#get ready for testing
include(CTest)
#enable_testing()

set(CTEST_PROJECT_NAME "SixTrack")

#   add_test(NAME <name> [CONFIGURATIONS [Debug|Release|...]]
#            [WORKING_DIRECTORY dir]
#            COMMAND <command> [arg1 [arg2 ...]])
	#ADD_TEST(NAME bb WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/SixTest/bb/ COMMAND SixTrack)

SET(SIXTRACK_TESTS bnl exact crabamp javier thick4 lost tilt dipedge bb bb_ntwin1 notilt lostnotilt lostevery s316 eric frs frs60 bbe51 bbe52 bbe571ib0 distance prob3 prob1 fma elensidealthck4d elensidealthck6d elensidealthin4d elensidealthin6d)
foreach(TEST IN ITEMS ${SIXTRACK_TESTS})
	ADD_TEST(NAME ${TEST} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/SixTest/${TEST} COMMAND SixTrack)
endforeach()
