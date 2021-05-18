include_directories(${ROOT_INCLUDE_DIRS})
add_definitions(${ROOT_DEFINITIONS})
set(CMAKE_CXX_FLAGS ${ROOT_CXX_FLAGS} ${CMAKE_CXX_FLAGS})

set(ROOT_OUTPUT_SOURCES SixTrack_root.cpp
ApertureCheck_root.cpp
Collimation_root.cpp
Optics_root.cpp
AcceleratorOutput_root.cpp
RunTime_root.cpp
ConfigurationOutput_root.cpp
ApertureDump_root.cpp
BunchDump_root.cpp
)
add_library(root_output STATIC ${ROOT_OUTPUT_SOURCES})
target_link_libraries(root_output ${ROOT_LIBRARIES})
