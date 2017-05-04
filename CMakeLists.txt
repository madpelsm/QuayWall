cmake_minimum_required(VERSION 2.6)
project(Kaaimuur)
#if(CMAKE_VERSION VERSION_GREATER 3)
    #cmake_policy(SET CMP0054 NEW)
#endif() 
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")



include_directories(${PROJECT_SOURCE_DIR}/ext/)
include_directories(${PROJECT_SOURCE_DIR}/src/)

message("source dir = ${PROJECT_SOURCE_DIR}")

file(GLOB KAAIMUUR_SRC_FILES "./src/*.cpp")
file(GLOB KAAIMUUR_HDR_FILES "./src/*.h")

ADD_EXECUTABLE(Kaaimuur ${KAAIMUUR_SRC_FILES} ${KAAIMUUR_HDR_FILES} ) 
target_link_libraries(Kaaimuur)
