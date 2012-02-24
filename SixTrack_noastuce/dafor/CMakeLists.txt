
add_executable(dafor dafor.f)

# this macro can be used to use dafor..
macro(daforize basename)
   set(_INCFLAGS "-I${CMAKE_CURRENT_SOURCE_DIR} -I${CMAKE_CURRENT_BINARY_DIR}")
add_custom_command(OUTPUT ${basename}.f
   COMMAND ${CMAKE_SOURCE_DIR}/cmake/precomp.sh "${CMAKE_Fortran_COMPILER}" "${_INCFLAGS}" "${CMAKE_CURRENT_SOURCE_DIR}/${basename}.F" "${CMAKE_CURRENT_BINARY_DIR}/${basename}.f"
   DEPENDS ${basename}.F
   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_custom_command(OUTPUT ${basename}.daf.f
   COMMAND cp ${basename}.f fort.1
   COMMAND dafor
   COMMAND cp fort.2 ${basename}.daf.f
   DEPENDS ${basename}.f
   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endmacro()
