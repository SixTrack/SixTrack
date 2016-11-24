#get ready for testing
#include(CTest)
#enable_testing()

#set(CTEST_PROJECT_NAME "SixTrack")
#set(CTEST_NIGHTLY_START_TIME "01:00:00 UTC")
#
#set(CTEST_DROP_METHOD "https")
#set(CTEST_DROP_SITE "localhost")
#set(CTEST_DROP_LOCATION "CDash/submit.php?project=Sixtrack")
#set(CTEST_DROP_SITE_CDASH TRUE)

#   add_test(NAME <name> [CONFIGURATIONS [Debug|Release|...]]
#            [WORKING_DIRECTORY dir]
#            COMMAND <command> [arg1 [arg2 ...]])
	#ADD_TEST(NAME bb WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/SixTest/bb/ COMMAND SixTrack)




   ENABLE_TESTING()
   INCLUDE(CTest)

set(CTEST_PROJECT_NAME "Sixtrack")
set(CTEST_NIGHTLY_START_TIME "01:00:00 UTC")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "abp-cdash.web.cern.ch")
set(CTEST_DROP_LOCATION "/abp-cdash/submit.php?project=SixTrack")
set(CTEST_DROP_SITE_CDASH TRUE)

SET(SIXTRACK_TESTS bnl exact crabamp javier thick4 lost tilt dipedge bb bb_ntwin1 notilt lostnotilt lostevery s316 eric frs frs60 bbe51 bbe52 bbe571ib0 distance prob3 prob1 fma elensidealthck4d elensidealthck6d elensidealthin4d elensidealthin6d)
foreach(TEST IN ITEMS ${SIXTRACK_TESTS})
	ADD_TEST(NAME ${TEST} WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/SixTest/${TEST} COMMAND SixTrack)
endforeach()
