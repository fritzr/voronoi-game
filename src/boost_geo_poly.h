/* Adapt boost::geometry types to be usable with boost::polygon.  */
#pragma once

#include <boost/polygon/polygon.hpp>
#include <boost/polygon/point_concept.hpp>
#include <boost/polygon/rectangle_concept.hpp>
#include <boost/polygon/voronoi.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
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

#ifndef BOOST_POLY_REGISTER_POINT_CONST
#define BOOST_POLY_REGISTER_POINT_CONST(ptype, ctype, GetX, GetY) \
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
    return (orient == HORIZONTAL) ? pt. GetX : pt. GetY; \
  } \
}; \
}} // end namespace boost::polygon
#endif

#define _EXPAND(...) __VA_ARGS__

#ifndef BOOST_POLY_REGISTER_MPOINT_GET_SET
#define BOOST_POLY_REGISTER_MPOINT_GET_SET(ptype, ctype, GetX, GetY, SetX, SetY) \
_EXPAND(BOOST_POLY_REGISTER_POINT_CONST(ptype, ctype, GetX, GetY)) \
namespace boost { namespace polygon { \
template <> struct \
point_mutable_traits<ptype> \
{ \
  typedef ctype coordinate_type; \
  static inline void set(ptype& point, orientation_2d orient, \
      coordinate_type value) \
  { \
    switch(orient.to_int()) { \
    case VERTICAL: point. SetY (value); break; \
    case HORIZONTAL: default: point. SetX (value); break; \
    } \
  } \
  static inline ptype construct(coordinate_type xval, coordinate_type yval) \
    { return ptype(xval, yval); } \
}; \
}} // end namespace boost::polygon
#endif

// Tell boost about geometry points
BOOST_POLY_REGISTER_MPOINT_GET_SET(
    boost::geometry::model::d2::point_xy<float>, float, x(), y(), x, y);
BOOST_POLY_REGISTER_MPOINT_GET_SET(
    boost::geometry::model::d2::point_xy<double>, double, x(), y(), x, y);

namespace bp = boost::polygon;
namespace bg = boost::geometry;
namespace bgm = boost::geometry::model;
