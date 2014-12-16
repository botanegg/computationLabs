cmake_minimum_required(VERSION 2.8.4)
project(computationLabs)

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

#lab 04
set(SOURCE_FILES lab05/lab05.cpp)
add_executable(lab05 ${SOURCE_FILES})

#gen_gilbert_input
set(SOURCE_FILES gen_gilbert_input.cpp)
add_executable(gen_gilbert_input ${SOURCE_FILES})
