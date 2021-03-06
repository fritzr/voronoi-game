cmake_minimum_required(VERSION 2.8.12)

#set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")
#cmake_policy(PUSH)
#include(cotire)
#cmake_policy(POP)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake/Modules)

project(voronoi-game)
if(NOT WIN32 AND NOT WIN64 AND NOT APPLE)
  set(CMAKE_CXX_COMPILER "g++")
endif()

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  execute_process(
    COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
  string(STRIP ${GCC_VERSION} GCC_VERSION)
  if (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
    message(FATAL_ERROR
      "${PROJECT_NAME} requires g++ 4.7 or greater for C++11 support"
      " (got ${CMAKE_CXX_COMPILER} == ${GCC_VERSION})")
  endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
else()
  message (FATAL_ERROR "${CMAKE_CXX_COMPILER} does not appear to support C++11")
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
if (${DEBUG})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -DDEBUG")
else()
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()

# Unfortunately FIST is 32-bit so we must build everything in 32 bits
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m32")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -m32")
set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -m32")
set (CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} -m32")

set(CXX_STATIC "-static-libgcc -static-libstdc++ ")
if (GCC_VERSION VERSION_GREATER 5.0 AND GCC_VERSION VERSION_LESS 5.4)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_STATIC}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CXX_STATIC}")
endif()

include_directories(include)
add_subdirectory(src)
