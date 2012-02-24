# Set up header macro generations
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lielib.h.in ${CMAKE_CURRENT_BINARY_DIR}/lielib.h)

add_library(lielib lielib.F)

target_link_libraries(lielib dabnew)
