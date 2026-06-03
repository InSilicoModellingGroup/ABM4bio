cmake_minimum_required(VERSION 3.19.3)

project(ABM4bio)

find_package(BioDynaMo REQUIRED)
include(${BDM_USE_FILE})
include_directories("src")

file(GLOB_RECURSE HEADERS src/*.h)
file(GLOB_RECURSE SOURCES src/*.cc)

bdm_add_executable(ABM4bio
                   HEADERS ${HEADERS}
                   SOURCES ${SOURCES}
                   LIBRARIES ${BDM_REQUIRED_LIBRARIES})
