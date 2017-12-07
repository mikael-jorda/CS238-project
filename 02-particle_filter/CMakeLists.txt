# create an executable
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BRACING2D_BINARY_DIR}/02-particle_filter)
ADD_EXECUTABLE (app02 app.cpp ${BRACING2D_COMMON_SOURCE}) # taco_library/include/taco/optimizer/QuadProgRT.cpp)
#include_directories(taco_library/include)

# and link the library against the executable
TARGET_LINK_LIBRARIES (app02 ${BRACING2D_COMMON_LIBRARIES})

# export resources such as model files.
# NOTE: this requires an install build
SET(RESOURCE_DIR ${BRACING2D_BINARY_DIR}/resources)
FILE(MAKE_DIRECTORY ${RESOURCE_DIR})
SET(APP_RESOURCE_DIR ${RESOURCE_DIR}/02-particle_filter)
FILE(MAKE_DIRECTORY ${APP_RESOURCE_DIR})
FILE(COPY world.urdf DESTINATION ${APP_RESOURCE_DIR})
FILE(COPY 4pbot_fixed.urdf DESTINATION ${APP_RESOURCE_DIR})
