cmake_minimum_required(VERSION 2.8.4)
project(computationLabs)

#lab 01
set(SOURCE_FILES lab01/lab01_rot.cpp)
add_executable(lab01_rot ${SOURCE_FILES})

set(SOURCE_FILES lab01/lab01_lu.cpp)
add_executable(lab01_lu ${SOURCE_FILES})