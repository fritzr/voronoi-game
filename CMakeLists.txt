cmake_minimum_required(VERSION 2.8.12)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")
cmake_policy(PUSH)
include(cotire)
cmake_policy(POP)

project(VDNative)
if(NOT WIN32 AND NOT WIN64)
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

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -std=c++11")

find_package(OpenCV REQUIRED)
find_package(Boost REQUIRED)

include_directories(
  ${OpenCV_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  )

if (OpenCV_VERSION VERSION_LESS 2.13)
  add_definitions(-DOPENCV2_PRE13)
endif()

set(CXX_STATIC "-static-libgcc -static-libstdc++ ")
if (GCC_VERSION VERSION_GREATER 5.0 AND GCC_VERSION VERSION_LESS 5.4)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_STATIC}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CXX_STATIC}")
endif()

# for debugging to dump cmake environment variables
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#  message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

set(PROG vgtest)

add_executable(${PROG} vgtest.cc vgame.h )

target_link_libraries(${PROG}
  ${OpenCV_LIBS}
  ${Boost_LIBRARIES}
  )

set_target_properties(${PROG} PROPERTIES
  COTIRE_ENABLE_PRECOMPILED_HEADER TRUE
  )
cotire(${PROG})


install (TARGETS ${PROG} DESTINATION bin)
