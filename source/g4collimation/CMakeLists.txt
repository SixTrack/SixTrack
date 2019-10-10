include_directories(${Geant4_INCLUDE_DIRS})         # Add -I type paths
add_definitions(${Geant4_DEFINITIONS})              # Add -D type defs
set(CMAKE_CXX_FLAGS "${Geant4_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")

if(ROOT)
	include_directories(${ROOT_INCLUDE_DIRS})
	add_definitions("${ROOT_DEFINITIONS} -DUSE_ROOT_FLAG")
	set(CMAKE_CXX_FLAGS "${ROOT_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
endif(ROOT)

set(G4COLLIMATION_SOURCES G4Interface.cpp CollimationGeometry.cpp CollimationParticleGun.cpp CollimationStackingAction.cpp CollimationMaterials.cpp CollimationEventAction.cpp CollimationTrackingAction.cpp  )

if(ROOT)
	set(G4COLLIMATION_SOURCES ${G4COLLIMATION_SOURCES} CollimationJawSD.cpp CollimationJawHit.cpp RootEDepHist.cpp RootEnergyDeposition.cpp)
endif(ROOT)

add_library(g4collimation STATIC ${G4COLLIMATION_SOURCES})
