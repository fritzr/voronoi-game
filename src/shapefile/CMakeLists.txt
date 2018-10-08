find_package(OpenCV REQUIRED)
find_package(FIST REQUIRED)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src/shapelib
  ${FIST_INCLUDE_DIR}
  ${OpenCV_INCLUDE_DIRS}
)

add_library(shapefile SHARED
  gettime.h
  intersection.h
  polygon.cpp
  polygon.h
  shpReader.cpp
  shpReader.h
)

target_link_libraries(shapefile
  ${FIST_LIBRARY}
  ${OpenCV_LIBS}
)

install (TARGETS shapefile DESTINATION lib)
