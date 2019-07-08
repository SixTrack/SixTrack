include_directories(${Geant4_INCLUDE_DIRS})         # Add -I type paths
add_definitions(${Geant4_DEFINITIONS})              # Add -D type defs
set(CMAKE_CXX_FLAGS ${Geant4_CXX_FLAGS} ${CMAKE_CXX_FLAGS})

set(G4COLLIMATION_SOURCES G4Interface.cpp CollimationGeometry.cpp CollimationParticleGun.cpp CollimationStackingAction.cpp CollimationMaterials.cpp CollimationEventAction.cpp CollimationTrackingAction.cpp)
add_library(g4collimation STATIC ${G4COLLIMATION_SOURCES})
