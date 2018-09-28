cmake_minimum_required(VERSION 2.8)
include($ENV{ROS_ROOT}/core/rosbuild/rosbuild.cmake)

rosbuild_init()

IF(NOT ROS_BUILD_TYPE)
  SET(ROS_BUILD_TYPE Release)
ENDIF()

MESSAGE("Build type: " ${ROS_BUILD_TYPE})

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -march=native -std=c++11")
set(CMAKE_CXX_FLAGS "${SSE_FLAGS} -O2 -g -std=c++0x -march=native -fno-omit-frame-pointer")

LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../cmake)

# required libraries
find_package(SuiteParse REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost COMPONENTS system thread)
find_package(OpenCV 3.0 QUIET)
if(NOT OpenCV_FOUND)
 find_package(OpenCV 2.4.3 REQUIRED)
endif()
find_package(Pangolin REQUIRED)

include_directories(
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/../../
  ${PROJECT_SOURCE_DIR}/../../src
  ${PROJECT_SOURCE_DIR}/../../thirdparty/Sophus
  ${EIGEN3_INCLUDE_DIR}
  ${Pangolin_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIR}
  ${CSPARSE_INCLUDE_DIR}
  ${CHOLMOD_INCLUDE_DIR}
)

set(LIBS 
  ${OpenCV_LIBS} 
  ${Pangolin_LIBRARIES}
  ${Boost_LIBRARIES}
  boost_system cxsparse boost_thread
  ${PROJECT_SOURCE_DIR}/../../build/lib/libdso.so
)

# Node for monocular camera
rosbuild_add_executable(DSO
  src/main_dso_ros.cc
)
target_link_libraries(DSO ${LIBS})



