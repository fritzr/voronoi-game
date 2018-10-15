/* Adapt OpenCV types to be usable with boost::polygon/boost:geometry.  */
#pragma once

#include <opencv2/core/core.hpp>

#include <boost/polygon/polygon.hpp>
#include <boost/polygon/point_concept.hpp>
#include <boost/polygon/rectangle_concept.hpp>
#include <boost/polygon/voronoi.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/register/box.hpp>
#include <boost/geometry/geometries/adapted/boost_polygon.hpp>

#ifndef BOOST_POLY_REGISTER_POINT
#define BOOST_POLY_REGISTER_POINT(ptype, ctype) \
namespace boost { namespace polygon { \
template <> struct \
geometry_concept<ptype> \
{ \
  typedef point_concept type; \
}; \
template <> struct \
point_traits<ptype> \
{ \
  typedef ctype coordinate_type; \
  static inline coordinate_type get(const ptype &pt, orientation_2d orient) \
  { \
    return (orient == HORIZONTAL) ? pt.x : pt.y; \
  } \
}; \
}} // end namespace boost::polygon
#endif

#ifndef BOOST_POLY_REGISTER_POINT
#define BOOST_POLY_REGISTER_POINT(ptype, ctype) \
namespace boost { namespace polygon { \
template <> struct \
geometry_concept<ptype> \
{ \
  typedef point_concept type; \
}; \
template <> struct \
point_traits<ptype> \
{ \
  typedef ctype coordinate_type; \
  static inline coordinate_type get(const ptype &pt, orientation_2d orient) \
  { \
    return (orient == HORIZONTAL) ? pt.x : pt.y; \
  } \
}; \
}} // end namespace boost::polygon
#endif

#define _EXPAND(...) __VA_ARGS__

#ifndef BOOST_POLY_REGISTER_MPOINT
#define BOOST_POLY_REGISTER_MPOINT(ptype, ctype) \
_EXPAND(BOOST_POLY_REGISTER_POINT(ptype, ctype)) \
namespace boost { namespace polygon { \
template <> struct \
point_mutable_traits<ptype> \
{ \
  typedef ctype coordinate_type; \
  static inline void set(ptype& point, orientation_2d orient, \
      coordinate_type value) \
  { \
    switch(orient.to_int()) { \
    case VERTICAL: point.y = value; break; \
    case HORIZONTAL: default: point.x = value; break; \
    } \
  } \
  static inline ptype construct(coordinate_type xval, coordinate_type yval) \
    { return ptype(xval, yval); } \
}; \
}} // end namespace boost::polygon
#endif

// Tell boost about cv points
BOOST_POLY_REGISTER_MPOINT(cv::Point, int);
BOOST_POLY_REGISTER_MPOINT(cv::Point2f, float);
BOOST_POLY_REGISTER_MPOINT(cv::Point2d, double);

BOOST_GEOMETRY_REGISTER_POINT_2D(cv::Point, int, cs::cartesian, x, y);
BOOST_GEOMETRY_REGISTER_POINT_2D(cv::Point2f, float, cs::cartesian, x, y);
BOOST_GEOMETRY_REGISTER_POINT_2D(cv::Point2d, double, cs::cartesian, x, y);
