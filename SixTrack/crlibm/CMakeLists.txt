#CMakeLists.txt for crlibm

SET(CRLIBM1_SOURCES csh_fast.c rem_pio2.c trigo_fast.c tan.c sine.c cosine.c exp.c exp_fast.c log.c log_fast.c log10.c atan.c atan_fast.c ericc.c disable_xp.c enable_xp.c)
SET(CRLIBM2_SOURCES addition_scs.c double2scs.c multiplication_scs.c division_scs.c scs2double.c zero_scs.c)
SET(CRLIBM3_SOURCES round_near.c dtoa_c.c dtoaf.c)

add_library(crlibm1 OBJECT ${CRLIBM1_SOURCES})
add_library(crlibm2 OBJECT ${CRLIBM2_SOURCES})
add_library(crlibm3 OBJECT ${CRLIBM3_SOURCES})

SET(CRLIBM_DEFS "-DLINUX_INLINE -DHAVE_CONFIG_H -I.")

target_compile_definitions(crlibm1 PRIVATE ${CRLIBM_DEFS})
target_compile_definitions(crlibm2 PRIVATE ${CRLIBM_DEFS})

if(32BIT)
	target_compile_options(crlibm1 PRIVATE -m32)
	target_compile_options(crlibm2 PRIVATE -m32)
	target_compile_options(crlibm3 PRIVATE -m32)
endif(32BIT)
if(64BIT)
	target_compile_options(crlibm1 PRIVATE -m64)
	target_compile_options(crlibm2 PRIVATE -m64)
	target_compile_options(crlibm3 PRIVATE -m64)
endif(64BIT)

if(${CMAKE_C_COMPILER_ID} MATCHES "GNU")
	target_compile_options(crlibm1 PRIVATE -mfpmath=sse -msse2 -fPIC -std=c99 -Wall -Wshadow -Wpointer-arith -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline -pedantic -fno-strict-aliasing)
	target_compile_options(crlibm2 PRIVATE -mfpmath=sse -msse2 -Wall -Wshadow -Wpointer-arith -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline)
	target_compile_options(crlibm3 PRIVATE -mfpmath=sse -msse2 -std=c99 -W -Wall -pedantic)

	if(COVERAGE)
		target_compile_options(crlibm1 PRIVATE -fprofile-arcs -ftest-coverage)
		target_compile_options(crlibm2 PRIVATE -fprofile-arcs -ftest-coverage)
		target_compile_options(crlibm3 PRIVATE -fprofile-arcs -ftest-coverage)
	endif(COVERAGE)
elseif(${CMAKE_C_COMPILER_ID} MATCHES "Intel")
	target_compile_options(crlibm1 PRIVATE -msse2 -fPIC -std=c99 -Wall -Wshadow -Wpointer-arith -Wconversion -Wstrict-prototypes -Winline -pedantic -fno-strict-aliasing)
	target_compile_options(crlibm2 PRIVATE -msse2 -Wall -Wshadow -Wpointer-arith -Wconversion -Wstrict-prototypes -Winline)
	target_compile_options(crlibm3 PRIVATE -msse2 -std=c99 -Wall -pedantic)
else()
	target_compile_options(crlibm1 PRIVATE -mfpmath=sse -msse2 -fPIC -std=c99 -Wall -Wshadow -Wpointer-arith -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline -pedantic -fno-strict-aliasing)
	target_compile_options(crlibm2 PRIVATE -mfpmath=sse -msse2 -Wall -Wshadow -Wpointer-arith -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes -Wnested-externs -Wlong-long -Winline)
	target_compile_options(crlibm3 PRIVATE -mfpmath=sse -msse2 -std=c99 -W -Wall -pedantic)
endif()
	
#glue the parts together
add_library(crlibm STATIC $<TARGET_OBJECTS:crlibm1> $<TARGET_OBJECTS:crlibm2> $<TARGET_OBJECTS:crlibm3>)

if(COVERAGE)
	SET_TARGET_PROPERTIES(crlibm PROPERTIES LINK_FLAGS "${CMAKE_C_LINKER_FLAGS}")
endif(COVERAGE)
