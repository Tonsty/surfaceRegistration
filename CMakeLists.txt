CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(PROJECT_NAME SurfaceRegistration)

SET(HEADER_FILES 	
	algorithm.h
	preProcessing.h
	common/inputUtilities.h
	common/PCAStatistics.h
	common/visualizationUtilities.h)
SET(SOURCE_FILES 
	main.cpp
	algorithm.cpp
	preProcessing.cpp
	common/inputUtilities.cpp
	common/PCAStatistics.cpp
	common/visualizationUtilities.cpp)

PROJECT(${PROJECT_NAME})

SOURCE_GROUP(common REGULAR_EXPRESSION ".*/common/.*" )

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

find_package(Eigen3 REQUIRED)
include_directories(EIGEN3_INCLUDE_DIR)

SET(Boost_USE_MULTITHREAD ON)
FIND_PACKAGE(Boost REQUIRED)
INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${LINK_DIRECTORIES} ${Boost_LIBRARY_DIRS})

ADD_EXECUTABLE(${PROJECT_NAME} ${SOURCE_FILES} ${HEADER_FILES})

if(VTK_LIBRARIES)
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${VTK_LIBRARIES})
else()
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} vtkHybrid vtkWidgets)
endif()