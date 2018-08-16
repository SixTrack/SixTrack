
SET(LIBFLUKAIO_SOURCES src/Connection.c src/FlukaIO.c src/FlukaIOHandshake.c src/FlukaIOServer.c src/FortranFlukaIO.c src/NetIO.c)

include_directories(${CMAKE_SOURCE_DIR}/source/FlukaIO/include)

add_library(flukaIO STATIC ${LIBFLUKAIO_SOURCES})

if(${CMAKE_SYSTEM_PROCESSOR} MATCHES AMD64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES x86_64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES i686)
	if(32BIT)
		target_compile_options(flukaIO PRIVATE -m32)
	endif(32BIT)

	if(64BIT)
		target_compile_options(flukaIO PRIVATE -m64)
	endif(64BIT)
endif()
