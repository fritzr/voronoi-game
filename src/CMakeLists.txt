find_package(OpenCV REQUIRED)
find_package(Boost REQUIRED)

# require boost 1.47.0 for polygon and geometry packages
if (Boost_VERSION LESS 105200)
  message (FATAL_ERROR "${PROJECT_NAME} requires boost 1.52.0")
endif()

set(CXX_STATIC "-static-libgcc -static-libstdc++ ")
if (GCC_VERSION VERSION_GREATER 5.0 AND GCC_VERSION VERSION_LESS 5.4)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_STATIC}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CXX_STATIC}")
endif()

# annoying warnings embedded in boost library code for recent gcc
if (GCC_VERSION VERSION_GREATER 7.0)
  add_definitions(-Wno-class-memaccess)
endif()

# for debugging to dump cmake environment variables
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#  message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

set(LIB_INCLUDE "${CMAKE_SOURCE_DIR}/lib")

include_directories(
  ${OpenCV_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${LIB_INCLUDE}
  )

set(TEST_PROG vgtest)

add_executable(${TEST_PROG}
  voronoi.h voronoi.cc
  sweep.h
  components.h components.cc
  main.h main.cc
)
#add_executable(l23test l23test.cc)

# these are just annoying
add_definitions(-Wno-parentheses)

target_link_libraries(${TEST_PROG}
  ${OpenCV_LIBS}
  ${Boost_LIBRARIES}
  )

#target_link_libraries(l23test ${Boost_LIBRARIES} )

#cotire(${TEST_PROG})

install (TARGETS ${TEST_PROG} DESTINATION bin)
#install (TARGETS l23test DESTINATION bin)
