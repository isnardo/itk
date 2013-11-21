# This is the root ITK CMakeLists file.
cmake_minimum_required(VERSION 2.4)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND CMAKE_POLICY)


# This project is designed to be built outside the Insight source tree.
project(OpticalFlow)

# Find ITK
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(OpticalFlow OpticalFlow.cxx)

target_link_libraries(OpticalFlow ${ITK_LIBRARIES})

