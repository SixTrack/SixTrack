# Set up header macro generations
if(DYNAPER)
   option(_BIG "Defined by DYNAPER" ON)
   mark_as_advanced(_BIG)
else()
   option(_SMALL "Defined by DYNAPER OFF" ON)
   mark_as_advanced(_SMALL)
endif()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/dabnew.h.in ${CMAKE_CURRENT_BINARY_DIR}/dabnew.h)



add_library(dabnew dabnew.F)
