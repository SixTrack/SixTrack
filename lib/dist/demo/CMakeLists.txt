cmake_minimum_required(VERSION 3.5)

project(hello)
enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f2008 -fimplicit-none")
    set(bounds "-fbounds-check")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

#
# Compile.
#
file(GLOB_RECURSE sources  ../source/*.c ../source/*.h ../source/*.f90 ../demo/*.f90)
add_library(${PROJECT_NAME} SHARED ${sources})

# Install library
install(TARGETS ${PROJECT_NAME} DESTINATION lib/${PROJECT_NAME})

# Install library headers
#file(GLOB HEADERS include/*.h)
#install(FILES ${HEADERS} DESTINATION include/${PROJECT_NAME})

add_executable(prog ${sources})
