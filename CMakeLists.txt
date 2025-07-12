cmake_minimum_required(VERSION 3.12)
project(PicoHFAnalysis)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ROOT
find_package(ROOT REQUIRED COMPONENTS RIO Tree Hist MathCore)
include(${ROOT_USE_FILE})

# Look for STAR/StRoot installation (use STAR environment variable if set)
if(DEFINED ENV{STAR})
    message(STATUS "Using STAR software from $ENV{STAR}")
    list(APPEND CMAKE_PREFIX_PATH "$ENV{STAR}")
endif()

find_path(STAR_INCLUDE_DIR StMaker.h PATH_SUFFIXES StRoot)
find_library(STAR_LIB St_base PATH_SUFFIXES lib)

if(NOT STAR_INCLUDE_DIR OR NOT STAR_LIB)
    message(FATAL_ERROR "Could not find STAR/StRoot headers or libraries. Set $STAR or specify paths manually.")
endif()

include_directories(${STAR_INCLUDE_DIR})
link_directories(${STAR_LIB})

file(GLOB_RECURSE SRC_FILES src/*.cxx)
add_library(PicoHFAnalysis SHARED ${SRC_FILES})

target_link_libraries(PicoHFAnalysis ${ROOT_LIBRARIES} St_base StPicoDst)

# Install step
install(TARGETS PicoHFAnalysis DESTINATION lib)
install(DIRECTORY src/ DESTINATION include FILES_MATCHING PATTERN "*.h")
