# Set up header macro generations
include(sixtrack_defs)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/sixtrack.h.in ${CMAKE_CURRENT_BINARY_DIR}/sixtrack.h)

if(DYNAPER)
   daforize(sixda)
   daforize(sixscfox)
   set(srcfiles sixda.daf.f sixsc.F sixscfox.daf.f beamgas.F90)
else()
   daforize(sixvefox)
   set(srcfiles track.F sixve.F sixvefox.daf.f beamgas.F90)
endif()

add_executable(sixtrack ${srcfiles})

target_link_libraries(sixtrack lielib)
