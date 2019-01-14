cmake_minimum_required(VERSION 3.6)
project(nbody)

set(CMAKE_CXX_COMPILER icpc)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11")

set(STEP_1 solution-step1.c)
set(STEP_2 solution-step2.c)
set(STEP_6 solution-step6.c)

SET_SOURCE_FILES_PROPERTIES(${STEP_1} PROPERTIES LANGUAGE CXX)
SET_SOURCE_FILES_PROPERTIES(${STEP_2} PROPERTIES LANGUAGE CXX)
SET_SOURCE_FILES_PROPERTIES(${STEP_6} PROPERTIES LANGUAGE CXX)

add_executable(step1 ${STEP_1})
add_executable(step2 ${STEP_2})
add_executable(step6 ${STEP_6})