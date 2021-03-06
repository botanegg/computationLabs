cmake_minimum_required(VERSION 2.8.4)
project(computationLabs)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(${CMAKE_SOURCE_DIR})

#lib
set(LIB_SOURCE_FILES
		share/matrix.cpp
		share/utils.cpp
		share/vector.cpp
		share/derivable.cpp
		)
add_library(computation ${LIB_SOURCE_FILES})
link_libraries(computation)

#lab 01
set(SOURCE_FILES lab01/lab01_rot.cpp)
add_executable(lab_01_rot ${SOURCE_FILES})

set(SOURCE_FILES lab01/lab01_lu.cpp)
add_executable(lab_01_lu ${SOURCE_FILES})

#lab 02
set(SOURCE_FILES lab02/lab02.cpp)
add_executable(lab_02 ${SOURCE_FILES})

#lab 03
set(SOURCE_FILES lab03/lab03.cpp)
add_executable(lab_03 ${SOURCE_FILES})

#lab 04
set(SOURCE_FILES lab04/lab04.cpp)
add_executable(lab_04 ${SOURCE_FILES})

#lab 05
set(SOURCE_FILES lab05/lab05.cpp)
add_executable(lab_05 ${SOURCE_FILES})


#########################################################
# FIND GLUT
#########################################################
find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIRS})
link_directories(${GLUT_LIBRARY_DIRS})
add_definitions(${GLUT_DEFINITIONS})
if (NOT GLUT_FOUND)
	message(ERROR " GLUT not found!")
endif (NOT GLUT_FOUND)
#########################################################
# FIND OPENGL
#########################################################
find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
link_directories(${OpenGL_LIBRARY_DIRS})
add_definitions(${OpenGL_DEFINITIONS})
if (NOT OPENGL_FOUND)
	message(ERROR " OPENGL not found!")
endif (NOT OPENGL_FOUND)


#lab 06
set(SOURCE_FILES lab06/lab06.cpp)
add_executable(lab_06 ${SOURCE_FILES})

target_link_libraries(lab_06
		${OPENGL_LIBRARIES}
		${GLUT_LIBRARY}
		)

#lab 07
set(SOURCE_FILES lab07/lab07.cpp)
add_executable(lab_07 ${SOURCE_FILES})

#lab 08
set(SOURCE_FILES lab08/lab08.cpp)
add_executable(lab_08 ${SOURCE_FILES})

#lab 09
set(SOURCE_FILES lab09/lab09.cpp)
add_executable(lab_09 ${SOURCE_FILES})

#lab 10
set(SOURCE_FILES lab10/lab10.cpp)
add_executable(lab_10 ${SOURCE_FILES})

#lab 11
set(SOURCE_FILES lab11/lab11.cpp)
add_executable(lab_11 ${SOURCE_FILES})

#lab 12
set(SOURCE_FILES lab12/lab12.cpp)
add_executable(lab_12 ${SOURCE_FILES})

#lab 13
set(SOURCE_FILES lab13/lab13.cpp)
add_executable(lab_13 ${SOURCE_FILES})
