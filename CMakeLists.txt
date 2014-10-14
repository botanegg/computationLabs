cmake_minimum_required(VERSION 2.8.4)
project(computationLabs)

include_directories(${CMAKE_SOURCE_DIR})

#lib
set(LIB_SOURCE_FILES
    share/matrix.cpp
    share/utils.cpp
    share/vector.cpp
)
add_library(computation ${LIB_SOURCE_FILES})
link_libraries(computation)

#lab 01
set(SOURCE_FILES lab01/lab01_rot.cpp)
add_executable(lab01_rot ${SOURCE_FILES})

set(SOURCE_FILES lab01/lab01_lu.cpp)
add_executable(lab01_lu ${SOURCE_FILES})

#lab 02
set(SOURCE_FILES lab02/lab02.cpp)
add_executable(lab02 ${SOURCE_FILES})

#gen_gilbert_input
set(SOURCE_FILES gen_gilbert_input.cpp)
add_executable(gen_gilbert_input ${SOURCE_FILES})