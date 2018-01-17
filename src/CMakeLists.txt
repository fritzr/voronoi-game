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

# for debugging to dump cmake environment variables
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#  message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

include_directories(
  ${OpenCV_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  )

set(TEST_PROG vgtest)

add_executable(${TEST_PROG} voronoi.h voronoi.cc main.h main.cc )

# use slow method for user/site mapping until we fix the fast method
#add_definitions(-DCELLS_SLOW)

# this is just annoying
add_definitions(-Wno-parentheses)

target_link_libraries(${TEST_PROG}
  ${OpenCV_LIBS}
  ${Boost_LIBRARIES}
  )

set_target_properties(${TEST_PROG} PROPERTIES
  COTIRE_ENABLE_PRECOMPILED_HEADER TRUE
  )
cotire(${TEST_PROG})

install (TARGETS ${TEST_PROG} DESTINATION bin)