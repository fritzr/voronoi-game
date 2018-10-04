find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${mpfr_INCLUDE_DIRS}
  ${MPFR_INCLUDES}
  ${GMP_INCLUDES}
)


add_library(mathtool SHARED
  Basic.cpp
  Basic.h
  Gauss.cpp
  Gauss.h
  Matrix.h
  Point.h
  Quaternion.cpp
  Quaternion.h
  Vector.h
  mpreal.h
)

target_link_libraries(mathtool
  ${MPFR_LIBRARIES}
  ${GMP_LIBRARIES}
)

install (TARGETS mathtool DESTINATION lib)
