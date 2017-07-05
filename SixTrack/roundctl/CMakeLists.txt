SET(LIBERIC_SOURCES round_near.c dtoa_c.c dtoaf.c enable_xp.c disable_xp.c) #Library for converting between Double <-> ASCII
add_library(roundctl STATIC ${LIBERIC_SOURCES})

if(${CMAKE_SYSTEM_PROCESSOR} MATCHES AMD64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES x86_64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES i686)
	if(32BIT)
	  target_compile_options(roundctl PRIVATE -m32)
	endif(32BIT)
	if(64BIT)
	  target_compile_options(roundctl PRIVATE -m64)
	endif(64BIT)
endif()

if(${CMAKE_C_COMPILER_ID} MATCHES "GNU")
	if(${CMAKE_SYSTEM_PROCESSOR} MATCHES AMD64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES x86_64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES i686)
		target_compile_options(roundctl PRIVATE -mfpmath=sse -msse2 -std=c99 -W -Wall -pedantic)
	else()
		target_compile_options(roundctl PRIVATE -std=c99 -W -Wall -pedantic)
	endif()

elseif(${CMAKE_C_COMPILER_ID} MATCHES "Intel")
	target_compile_options(roundctl PRIVATE -msse2 -std=c99 -Wall -pedantic)
else()
	target_compile_options(roundctl PRIVATE -mfpmath=sse -msse2 -std=c99 -W -Wall -pedantic)
endif()
