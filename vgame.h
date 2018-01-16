#include <opencv2/core.hpp>

#ifndef OPENCV2_4_13
#include <opencv2/contrib/contrib.hpp>
typedef int ColormapTypes;
#endif

#include <boost/unordered_map.hpp>
#include <boost/polygon/polygon.hpp>
#include <boost/polygon/point_concept.hpp>
#include <boost/polygon/rectangle_concept.hpp>
#include <boost/polygon/voronoi.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/register/box.hpp>
#include <boost/geometry/geometries/adapted/boost_polygon.hpp>

namespace bp = boost::polygon;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

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

// Tell boost about cv points
BOOST_POLY_REGISTER_POINT(cv::Point, int);
BOOST_POLY_REGISTER_POINT(cv::Point2f, float);
BOOST_POLY_REGISTER_POINT(cv::Point2d, double);

BOOST_GEOMETRY_REGISTER_POINT_2D(cv::Point, int, cs::cartesian, x, y);
BOOST_GEOMETRY_REGISTER_POINT_2D(cv::Point2f, float, cs::cartesian, x, y);
BOOST_GEOMETRY_REGISTER_POINT_2D(cv::Point2d, double, cs::cartesian, x, y);

typedef bp::rectangle_data<float> Rectf;
typedef bp::rectangle_data<double> Rectd;

typedef bp::voronoi_diagram<double> voronoi_diagram;
typedef voronoi_diagram::cell_type   cell_t;
typedef voronoi_diagram::vertex_type vertex_t;
typedef voronoi_diagram::edge_type   edge_t;
typedef voronoi_diagram::coordinate_type coordinate_t;

#if !defined(_WIN32) && !defined(_WIN64)
#include <unistd.h>
#include <getopt.h>
extern "C" {
	extern int opterr;    /* if error message should be printed */
	extern	int optind;    /* index into parent argv vector */
	extern int optopt;    /* character checked for validity */
	extern int optreset;  /* reset getopt */
	extern char    *optarg;
	extern int getopt(int nargc, char * const nargv[], const char *ostr);
};
#else
#error no getopt support
#endif
