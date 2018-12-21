cmake_minimum_required(VERSION 3.6)
project(Assignment)

set(CMAKE_CXX_COMPILER icpc)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11")

set(SOURCE_FILES assignment-2018.c)

SET_SOURCE_FILES_PROPERTIES(${SOURCE_FILES} PROPERTIES LANGUAGE CXX)

add_executable(Assignment ${SOURCE_FILES})