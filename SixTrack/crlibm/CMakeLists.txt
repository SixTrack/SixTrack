#CMakeLists.txt for crlibm

SET(CRLIBM1_SOURCES csh_fast.c rem_pio2.c trigo_fast.c tan.c sine.c cosine.c exp.c exp_fast.c log.c log_fast.c log10.c atan.c atan_fast.c)
SET(CRLIBM2_SOURCES addition_scs.c double2scs.c multiplication_scs.c division_scs.c scs2double.c zero_scs.c)

add_library(crlibm1 OBJECT ${CRLIBM1_SOURCES})
add_library(crlibm2 OBJECT ${CRLIBM2_SOURCES})

SET(CRLIBM_DEFS "-DLINUX_INLINE -DHAVE_CONFIG_H -I.")
SET(CRLIBM_DEFS_DEBUG "-DHAVE_CONFIG_H -I.")

if(CMAKE_BUILD_TYPE STREQUAL Debug)
	target_compile_definitions(crlibm1 PRIVATE ${CRLIBM_DEFS_DEBUG})
	target_compile_definitions(crlibm2 PRIVATE ${CRLIBM_DEFS_DEBUG})
else()
	target_compile_definitions(crlibm1 PRIVATE ${CRLIBM_DEFS})
	target_compile_definitions(crlibm2 PRIVATE ${CRLIBM_DEFS})
endif()

if(${CMAKE_SYSTEM_PROCESSOR} MATCHES AMD64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES x86_64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES i686)
	if(32BIT)
	  target_compile_options(crlibm1 PRIVATE -m32)
	  target_compile_options(crlibm2 PRIVATE -m32)
	endif(32BIT)
	if(64BIT)
	  target_compile_options(crlibm1 PRIVATE -m64)
	  target_compile_options(crlibm2 PRIVATE -m64)
	endif(64BIT)
endif()

if(${CMAKE_C_COMPILER_ID} MATCHES "GNU")
	# In old crlibm, suppress some warnings. It's too many of them to be usefull...
	set(lessWarnings  -Wno-sign-conversion -Wno-conversion -Wno-long-long)

	if(${CMAKE_SYSTEM_PROCESSOR} MATCHES AMD64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES x86_64 OR ${CMAKE_SYSTEM_PROCESSOR} MATCHES i686)  
		target_compile_options(crlibm1 PRIVATE -mfpmath=sse -msse2 -fPIC -std=c99 -Wall -Wshadow -Wpointer-arith -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline -pedantic -fno-strict-aliasing ${lessWarnings})
		target_compile_options(crlibm2 PRIVATE -mfpmath=sse -msse2 -Wall -Wshadow -Wpointer-arith -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline ${lessWarnings})
	else()
		target_compile_options(crlibm1 PRIVATE -fPIC -std=c99 -Wall -Wshadow -Wpointer-arith -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline -pedantic -fno-strict-aliasing ${lessWarnings})
		target_compile_options(crlibm2 PRIVATE -Wall -Wshadow -Wpointer-arith -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline ${lessWarnings})
	endif()
  
elseif(${CMAKE_C_COMPILER_ID} MATCHES "Intel")
	target_compile_options(crlibm1 PRIVATE -msse2 -fPIC -std=c99 -Wall -Wshadow -Wpointer-arith -Wconversion -Wstrict-prototypes -Winline -pedantic -fno-strict-aliasing)
	target_compile_options(crlibm2 PRIVATE -msse2 -Wall -Wshadow -Wpointer-arith -Wconversion -Wstrict-prototypes -Winline)
else()
	target_compile_options(crlibm1 PRIVATE -mfpmath=sse -msse2 -fPIC -std=c99 -Wall -Wshadow -Wpointer-arith -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline -pedantic -fno-strict-aliasing)
	target_compile_options(crlibm2 PRIVATE -mfpmath=sse -msse2 -Wall -Wshadow -Wpointer-arith -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline)
endif()

#glue the parts together
add_library(crlibm STATIC $<TARGET_OBJECTS:crlibm1> $<TARGET_OBJECTS:crlibm2>)

