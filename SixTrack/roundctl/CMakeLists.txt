SET(LIBERIC_SOURCES round_near.c dtoa_c.c dtoaf.c enable_xp.c disable_xp.c) #Library for converting between Double <-> ASCII
add_library(roundctl STATIC ${LIBERIC_SOURCES})

if(32BIT)
  target_compile_options(roundctl PRIVATE -m32)
endif(32BIT)
if(64BIT)
  target_compile_options(roundctl PRIVATE -m64)
endif(64BIT)

if(${CMAKE_C_COMPILER_ID} MATCHES "GNU")
  target_compile_options(roundctl PRIVATE -mfpmath=sse -msse2 -std=c99 -W -Wall -pedantic)
  if(COVERAGE)
    target_compile_options(roundctl PRIVATE -fprofile-arcs -ftest-coverage)
  endif(COVERAGE)
elseif(${CMAKE_C_COMPILER_ID} MATCHES "Intel")
  target_compile_options(roundctl PRIVATE -msse2 -std=c99 -Wall -pedantic)
else()
  target_compile_options(roundctl PRIVATE -mfpmath=sse -msse2 -std=c99 -W -Wall -pedantic)
endif()

if(COVERAGE)
  SET_TARGET_PROPERTIES(roundctl PROPERTIES LINK_FLAGS "${CMAKE_C_LINKER_FLAGS}")
endif(COVERAGE)
