cmake_minimum_required(VERSION 2.8.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(flute)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Choose the type of build, options are: Debug Release."
    FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# Embed lut.ICCAD2015 so readLUT_embedded() needs no runtime data files.
add_executable(flute_embed_gen ${CMAKE_CURRENT_SOURCE_DIR}/flute_embed_gen.cpp)
set(FLUTE_LUT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lut.ICCAD2015)
set(FLUTE_LUT_DATA ${CMAKE_CURRENT_BINARY_DIR}/flute_lut_data.cpp)
add_custom_command(
  OUTPUT ${FLUTE_LUT_DATA}
  COMMAND flute_embed_gen ${FLUTE_LUT_DIR}/POWV9.dat ${FLUTE_LUT_DIR}/POST9.dat ${FLUTE_LUT_DATA}
  DEPENDS flute_embed_gen ${FLUTE_LUT_DIR}/POWV9.dat ${FLUTE_LUT_DIR}/POST9.dat
  VERBATIM)

set(SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/flute.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/flute_embedded.cpp
  ${FLUTE_LUT_DATA})

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_library(${PROJECT_NAME} STATIC ${SOURCES})

install(TARGETS ${PROJECT_NAME} DESTINATION thirdparty/${PROJECT_NAME})
install(DIRECTORY lut.ICCAD2015 DESTINATION thirdparty/${PROJECT_NAME})
