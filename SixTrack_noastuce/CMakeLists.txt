cmake_minimum_required(VERSION 2.6)

project(sixtrack Fortran)

# project version
set( SIXTRACK_VERSION_MAJOR 0 )
set( SIXTRACK_VERSION_MINOR 95 )
set( SIXTRACK_VERSION_PATCH 0 )
set( SIXTRACK_VERSION ${SIXTRACK_VERSION_MAJOR}.${SIXTRACK_VERSION_MINOR}.${SIXTRACK_VERSION_PATCH})

# add source dir as a place for CMake modules (e.g. FindGSL.cmake)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# With this true we automatically include current src/binary directory
set(CMAKE_INCLUDE_CURRENT_DIR TRUE)

set(CMAKE_Fortran_COMPILE_FLAGS "${CMAKE_Fortran_COMPILE_FLAGS} -fno-second-underscore -funroll-loops")

include(common_defs)

add_subdirectory(dafor)
add_subdirectory(sixtrack)
add_subdirectory(lielib)
add_subdirectory(dabnew)
