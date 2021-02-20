project(geo1004)

cmake_minimum_required(VERSION 3.2)

add_definitions(-std=c++11)

include_directories( ${CMAKE_SOURCE_DIR}/voxelizer )

FILE(GLOB SRC_FILES voxeliser/*.cpp)
add_executable(geo1004 ${SRC_FILES})

if (WIN32)
    set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT "geo1004")
endif (WIN32)
