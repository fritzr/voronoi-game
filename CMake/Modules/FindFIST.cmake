set(FIST_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/include)

find_path(FIST_INCLUDE_DIR
  NAMES header.h basic.h defs.h graphics_header.h local.h numerics.h vertex.h
        bv_tree.h ext_appl.h header.h martin.h sgo.h data.h graphics.h ipe_io.h
        matrix.h triangulate.h
  PATHS ${CMAKE_SOURCE_DIR}/include
  PATH_SUFFIXES FIST
)

if (APPLE)
  set(FIST_LIBNAME libFIST-osx.a)
elseif (WIN32)
  set(FIST_LIBNAME libFIST.lib)
else() # UNIX
  set(FIST_LIBNAME libFIST-linux.a)
endif()

find_library(FIST_LIBRARY
  NAMES ${FIST_LIBNAME}
  PATHS ${CMAKE_SOURCE_DIR}/lib
  PATH_SUFFIXES FIST
)
get_filename_component(FIST_LIBRARY_DIR ${FIST_LIBRARY} DIRECTORY)

set(FIST_VERSION 1.0)

find_package_handle_standard_args(FIST
  FOUND_VAR FIST_FOUND
  REQUIRED_VARS
    FIST_LIBRARY
    FIST_INCLUDE_DIR
    FIST_LIBRARY_DIR
  VERSION_VAR FIST_VERSION
)

if(FIST_FOUND)
  set(FIST_INCLUDE_DIRS ${FIST_INCLUDE_DIR})
endif()

if(FIST_FOUND AND NOT TARGET FIST::FIST)
  add_library(${FIST_LIBNAME} STATIC IMPORTED)
  set_target_properties(${FIST_LIBNAME} PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${FIST_INCLUDE_DIR}"
  )
endif()
