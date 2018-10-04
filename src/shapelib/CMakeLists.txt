
add_library(shapelib SHARED
  dbfopen.cpp
  safileio.cpp
  shapefil.h
  shpopen.cpp
)

install (TARGETS shapelib DESTINATION lib)
