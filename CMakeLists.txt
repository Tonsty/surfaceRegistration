CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(PROJECT_NAME SurfaceRegistration)

SET(LIB_DIRECTORY lib)

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

include_directories(${LIB_DIRECTORY})

SOURCE_GROUP(common REGULAR_EXPRESSION ".*/common/.*" )

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(${PROJECT_NAME} ${SOURCE_FILES} ${HEADER_FILES})

if(VTK_LIBRARIES)
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${VTK_LIBRARIES})
else()
  TARGET_LINK_LIBRARIES(${PROJECT_NAME} vtkHybrid vtkWidgets)
endif()
