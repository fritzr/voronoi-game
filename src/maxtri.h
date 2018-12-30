#pragma once

#include <set>
#include <queue>
#include <vector>
#include <unordered_set>
#include <functional>
#include <boost/array.hpp>

#ifdef MAXTRI_DEBUG
#include <iostream>
#include <iomanip>
#endif

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp> // point_xy, ring
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/algorithms/is_convex.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/register/ring.hpp>

#include "util.h"
#include "user.h"
#include "polygon.h"
#include "intersection.h"
#include "sweep.h"
#include "nn1.h"

#ifdef MAXTRI_DEBUG
struct ioptr
{
  const uintptr_t ptrval;

  ioptr(const void *p) : ptrval(reinterpret_cast<uintptr_t>(p)) {}

  inline friend
  std::ostream &operator<<(std::ostream &os, ioptr const& iop)
  {
    os << std::hex << std::setw(2*__SIZEOF_POINTER__) << std::setfill('0')
      << iop.ptrval;
    return os;
  }
};
#endif

namespace cfla { namespace tri
{

// namespaces

namespace b  = boost;
namespace bp = boost::polygon;
namespace bg = boost::geometry;
namespace bgm = boost::geometry::model;

template<typename Tp_>
  using point_xy = boost::geometry::model::d2::point_xy<Tp_>;


} } // end namespace cfla::tri

// Adapt point_xy to boost::polygon point_concept
namespace boost { namespace polygon
{

  template<typename Tp_> struct
    geometry_concept<cfla::tri::point_xy<Tp_> >
    { typedef point_concept type; };

  template<typename Tp_> struct
    point_traits<cfla::tri::point_xy<Tp_> >
    {
      typedef typename bg::traits::coordinate_type<cfla::tri::point_xy<Tp_> >
        ::type coordinate_type;

      static inline coordinate_type get(const cfla::tri::point_xy<Tp_> &p,
          orientation_2d orient)
        { return (orient == HORIZONTAL) ? p.x() : p.y(); }
    };

  template<typename Tp_> struct
    point_mutable_traits<cfla::tri::point_xy<Tp_> >
    {
      typedef typename bg::traits::coordinate_type<cfla::tri::point_xy<Tp_> >
        ::type coordinate_type;

      static inline void set(cfla::tri::point_xy<Tp_> &p,
          orientation_2d orient, coordinate_type value) {
        switch (orient.to_int()) {
          case VERTICAL: p.y(value); break;
          case HORIZONTAL: p.x(value); break;
        }
      }

      static inline cfla::tri::point_xy<Tp_>
        construct(coordinate_type x, coordinate_type y) {
          return cfla::tri::point_xy<Tp_> (x, y);
        }
    };

} } // end namespace boost::polygon


namespace cfla { namespace tri
{

// forward declarations

template<typename Tp_>
class Triangle;

template<typename Tp_>
class Edge;

template<typename Tp_>
class EdgePoint;

template<typename Tp_>
class SolutionCell;

template<typename Tp_>
struct compare_status;

template<typename Tp_>
struct event_point_compare_y;

template<typename Tp_>
class TriEventPoint;

template<typename Tp_>
class IXPoint;

template<typename Tp_>
class EventPoint;

template<class Tp_>
class MaxTri;

template<class Tp_>
struct edge_data;

template<class Tp_>
struct triangle_data;

template<typename Tp_>
struct event_data;

template<typename Tp_>
class StatusSegment;

template<typename Tp_>
struct traits
{
  typedef Tp_                            coordinate_type;
  typedef point_xy<Tp_>                  point_type;

  typedef edge_data<Tp_>                 edge_data_type;
  typedef triangle_data<Tp_>             triangle_data_type;
  typedef event_data<Tp_>                event_data_type;

  typedef Edge<Tp_>                      edge_type;
  typedef EdgePoint<Tp_>                 edge_point_type;
  typedef Triangle<Tp_>                  triangle_type;

  typedef TriEventPoint<Tp_>             tri_point_type;
  typedef IXPoint<Tp_>                   isect_point_type;
  typedef EventPoint<Tp_>                event_point_type;
  typedef event_point_compare_y<Tp_>     event_point_ycompare;

  typedef StatusSegment<Tp_>             status_seg_type;
  typedef compare_status<Tp_>            status_compare;

  typedef c_polygon<point_type>          polygon_type;
  typedef SolutionCell<Tp_>              solution_ref_type;
  typedef bgm::ring<point_type>          solution_cell_type;
  typedef MaxTri<coordinate_type>        solver_type;
};

#define INHERIT_TRAITS(Targ) \
  typedef cfla::tri::traits<Targ> traits; \
  typedef typename traits::coordinate_type coordinate_type; \
  typedef typename traits::point_type point_type; \
  \
  typedef typename traits::edge_data_type edge_data_type; \
  typedef typename traits::triangle_data_type triangle_data_type; \
  typedef typename traits::event_data_type event_data_type; \
  \
  typedef typename traits::edge_type edge_type; \
  typedef typename traits::edge_point_type edge_point_type; \
  typedef typename traits::triangle_type triangle_type; \
  \
  typedef typename traits::tri_point_type tri_point_type; \
  typedef typename traits::isect_point_type isect_point_type; \
  typedef typename traits::event_point_type event_point_type; \
  typedef typename traits::event_point_ycompare event_point_ycompare; \
  \
  typedef typename traits::status_seg_type status_seg_type; \
  typedef typename traits::status_compare status_compare; \
  \
  typedef typename traits::polygon_type polygon_type; \
  typedef typename traits::solution_ref_type solution_ref_type; \
  typedef typename traits::solution_cell_type solution_cell_type; \
  typedef typename traits::solver_type solver_type; \

// functions

template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
    almost_equal(T x, T y, int dlp=4)
{
  // the machine epsilon has to be scaled to the magnitude of the values used
  // and multiplied by the desired precision in DLPs (decimals in the last place)
  return std::abs(x-y)
    <= std::numeric_limits<T>::epsilon() * std::abs(x+y) * (dlp << 2)
    // unless the result is subnormal
    || std::abs(x-y) < std::numeric_limits<T>::min();
}

/* Note that solution cells, which are formed from the intersection of many
 * triangles, will always be a convex shape with a single boundary.  */
namespace boost { namespace geometry {
  template<typename Tp_>
    inline bool is_convex(typename traits<Tp_>::solution_cell_type const& c)
    {
      return true;
    }
} } // end namespace boost::geometry

/* Intersect many geometries into out.  */
template<typename Iter, typename Polygon>
bool
intersection(Iter begin, Iter end, Polygon& out)
{
  if (begin == end)
    return false;

  // Case A. one geometry, return it
  auto const& first = *begin++;
  if (begin == end)
  {
    bg::assign(out, *first);
    return true;
  }

  // Case B. two geometries, return their intersection
  auto const& second = *begin++;
  if (!bg::intersection(*first, *second, out))
    return false;

  // Case C. many geometries, return the cumulative intersection with A&B
  while (begin != end)
    // XXX is providing out twice okay??
    if (!bg::intersection(*(*begin++), out, out))
      return false;

  return true;
}

template<class Point>
inline typename bg::coordinate_type<Point>::type
getx(Point const& p) {
  return bg::get<0>(p);
}

template<class Point>
inline typename bg::coordinate_type<Point>::type
gety(Point const& p) {
  return bg::get<1>(p);
}

// classes

template<class Tp_>
struct triangle_data
{
  INHERIT_TRAITS(Tp_);

  typedef edge_type edge_array[3];
  typedef std::array<point_type, 3> point_array;

  // members
  polygon_type const& poly; // parent polygon
  point_array points;
  edge_array edges;

  template <class Tri>
  triangle_data(polygon_type const& parent, Tri const& input_tri)
    : poly(parent), points()
      , edges({ // to be filled in by Triangle constructor
          { this, 0, input_tri[0], input_tri[0], bp::RIGHT },
          { this, 1, input_tri[0], input_tri[0], bp::LEFT },
          { this, 2, input_tri[0], input_tri[0], bp::LEFT },
        })
  {
  }

  // parent edge
  inline edge_type const& edge(int index) const {
    if (index < 0)
      throw std::runtime_error("triangle_data: negative edge index!");
    return edges[index % 3];
  }

  inline point_type const& point(int index) const {
    if (index < 0)
      throw std::runtime_error("triangle_data: negative point index!");
    return points[index % 3];
  }
};

template<class Tp_>
struct edge_data
{
  INHERIT_TRAITS(Tp_);

  // members
  edge_point_type first, second;
  triangle_data_type *tri; // parent triangle -- weak reference
  int index; // index of edge in parent triangle (0, 1, or 2)
  bp::direction_1d dir; // LEFT or RIGHT
  int depth; // intersection depth

  edge_data(triangle_data_type *tridata, int edge_index,
      point_type start, point_type end, bp::direction_1d d)
    : first(this, bp::LOW, start), second(this, bp::HIGH, end), tri(tridata)
      , index(edge_index), dir(d), depth(0)
  {}

  // parent edge
  inline edge_type const& edge(void) const {
    return tri->edges[index];
  }
};

template<class Tp_>
class EdgePoint : public point_xy<Tp_>
{
public:
  // typedefs
  typedef point_xy<Tp_> super;
  INHERIT_TRAITS(Tp_);

  // members
  edge_data<Tp_> *edata; // parent edge data -- weak reference
  bp::direction_1d dir; // LOW or HIGH, meaning first or second

  // constructors
  EdgePoint(edge_data<Tp_> *parent, bp::direction_1d d, super pt)
    : super(pt), edata(parent), dir(d)
  {}

  EdgePoint(EdgePoint const& o) : super(o), edata(o.edata), dir(o.dir) {}

  // methods

  void set(point_type const& p)
  {
    super::x(p.x());
    super::y(p.y());
  }

  // parent edge
  inline edge_type const& edge(void) const {
    return edata->tri->edges[edata->index];
  }

  // return the point at the other end of this edge
  inline edge_point_type const& other(void) const {
    return (dir == bp::LOW) ? edata->second : edata->first;
  }

  inline const triangle_data_type *tridata(void) const {
    return edata->tri;
  }

#ifdef MAXTRI_DEBUG
  inline friend std::ostream& operator<<(std::ostream& os, EdgePoint const& p) {
    os << "(" << p.x() << " , " << p.y() << ")";
    return os;
  }
#endif
};

// Point of a Triangle (index 0, 1, or 2)
template<class Tp_>
class TriEventPoint
{
public:
  INHERIT_TRAITS(Tp_);

  const triangle_data_type *tridata;
  int index;

  TriEventPoint(const triangle_data_type *tdata, int point_index)
    : tridata(tdata), index(point_index)
  {}

  point_type const& value(void) const { return tridata->point(index); }

 /* When looking into the center of the triangle from a TriEventPoint,
  * the two edges that share it are considered the "left" and "right" edges
  * with respect to the point according to the turn direction from the
  * bisector of the two edges to the other point of the edge.
  * The indexes are consistent regardless of the triangle's actual shape
  * according to their definitions in the comments for the Triangle class.
  * Here is the same example diagram:
  *
  *       [0]               [0]                               >
  *        *                    left edge: E1               / left edge
  * V E0  /  \  E1 V           right edge: E0   exterior  /
  *      /     \            [1]                         /
  * [1] *- _    \               left edge: E0         * - - - - - > bisector
  *         - _   \            right edge: E2           \  interior
  *             - _\        [2]                           \
  *        E2 >    * [2]        left edge: E2               \ right edge
  *                            right edge: E1                 >
  */

  inline edge_type const& left_edge(void) const
  {
    static const int left_edge_indexes[] = {
      /* point index -> left edge index */
      /* 0 -> */ 1,
      /* 1 -> */ 0,
      /* 2 -> */ 2,
    };
    return tridata->edges[left_edge_indexes[index]];
  }

  inline edge_type const& right_edge(void) const
  {
    static const int right_edge_indexes[] = {
      /* point index -> right edge index */
      /* 0 -> */ 0,
      /* 1 -> */ 2,
      /* 2 -> */ 1,
    };
    return tridata->edges[right_edge_indexes[index]];
  }

#ifdef MAXTRI_DEBUG
  inline friend
  std::ostream& operator<<(std::ostream& os, TriEventPoint const& p)
  {
    os << "[" << p.index << "] joining "
      << p.left_edge() << " + " << p.right_edge();
    return os;
  }
#endif
};

// Intersection point
template<typename Tp_>
class IXPoint : public point_xy<Tp_>
{
public:
  INHERIT_TRAITS(Tp_);

  typedef point_xy<Tp_> super;

  const edge_data<Tp_> *e1, *e2;

  using super::x;
  using super::y;

  IXPoint()
    : super(), e1(nullptr), e2(nullptr)
  {}

  template<typename Pt_>
  IXPoint(Pt_ const& point,
      const edge_data<Tp_> *edge1,
      const edge_data<Tp_> *edge2)
    : super(point), e1(edge1), e2(edge2)
  {}

#ifdef MAXTRI_DEBUG
  inline friend
  std::ostream& operator<<(std::ostream& os, IXPoint const& p)
  {
    os << static_cast<super const&>(p)
      << " from " << p.edge1() << " and " << p.edge2();
    return os;
  }
#endif

  inline edge_type const& edge1(void) const { return e1->edge(); }
  inline edge_type const& edge2(void) const { return e2->edge(); }

  // Assign to the segment a-b such that a is the intersection point and
  // b is the lowest endpoint of the preferred edge. The direction refers to
  // the edge direction past the point of intersection.
  // bp::LEFT means the left-most angled segment (converse to bp::RIGHT).
  status_seg_type make_segment(bp::direction_1d preferred_direction) const
  {
    point_type a(*this);

    const edge_type *left_edge = &edge2();
    const edge_type *right_edge = &edge1();
    if (!leftTurn(a, left_edge->second(), right_edge->second()))
    {
      const edge_type *left_swp = left_edge;
      left_edge = right_edge;
      right_edge = left_swp;
    }

    if (preferred_direction == bp::LEFT)
      return status_seg_type(*left_edge, a);

    else // (preferred_direction == bp::RIGHT)
      return status_seg_type(*right_edge, a);
  }
};

enum EventType
{
    TRIPOINT = 0
  , INTERSECTION
};

template<typename Tp_>
struct event_data
{
  INHERIT_TRAITS(Tp_);

  union {
    tri_point_type point;
    isect_point_type isect;
  } u;

  EventType type;

  event_data(tri_point_type const& e)
    : u{.point=e}, type(TRIPOINT) {}

  event_data(isect_point_type const& i)
    : u{.isect=i}, type(INTERSECTION) {}

  ~event_data()
  {
    switch (type)
    {
    case TRIPOINT:
      (u.point).~TriEventPoint();
    case INTERSECTION:
      (u.isect).~IXPoint();
    default:
      // should be unreachable!
      break;
    }
  }

};

template<typename Tp_>
class EventPoint
{
private:
  typedef event_data<Tp_> m_type;
  std::unique_ptr<m_type> m_union;

public:
  INHERIT_TRAITS(Tp_);

  EventPoint(tri_point_type const& point)
    : m_union(new m_type(point)) {}

  EventPoint(IXPoint<Tp_> const& ipoint)
    : m_union(new m_type(ipoint)) {}

  inline EventType type(void) const {
    return m_union->type;
  }

  inline tri_point_type const& point(void) const {
    if (type() != TRIPOINT)
      throw std::runtime_error("EventPoint: type mismatch (exp. TRIPOINT)");
    return m_union->u.point;
  }

  inline isect_point_type const& isect(void) const {
    if (type() != INTERSECTION)
      throw std::runtime_error("EventPoint: type mismatch (exp. INTERSECTION)");
    return m_union->u.isect;
  }

  inline coordinate_type x(void) const
  {
    switch (type())
    {
    case TRIPOINT:
      return bg::get<0>(point().value());
    case INTERSECTION:
      return bg::get<0>(isect());
    default:
      // should be unreachable!
      return -1;
    }
  }

  inline coordinate_type y(void) const
  {
    switch (type())
    {
    case TRIPOINT:
      return bg::get<1>(point().value());
    case INTERSECTION:
      return bg::get<1>(isect());
    default:
      // should be unreachable!
      return -1;
    }
  }

};

template<typename Tp_>
inline bool
_event_point_ycompare(
    typename traits<Tp_>::event_point_type const& p1,
    typename traits<Tp_>::event_point_type const& p2)
{
  // We've just found the points are equivalent by value, so compare by...

  // event type
  if (p1.type() < p2.type())
    return true;
  if (p1.type() != p2.type())
    return false;

  // then by point index (for triangle points)
  if (p1.type() == TRIPOINT)
    return p1.point().index < p2.point().index;

  // else // p1.type() == INTERSECTION
  // otherwise we really don't care, but if it comes to this
  // then std::set will explode
  return false;
}

// Unfortunately this can't be a member function before C++17, so it's best to
// delegate to a namespace-bound function for now.
template<typename Pt1_, typename Pt2_>
inline bool
_specialized_ycompare(Pt1_ const& p1, Pt2_ const& p2)
{
  return false;
}

template<>
inline bool
_specialized_ycompare<
    traits<float>::event_point_type
  , traits<float>::event_point_type>
  (typename traits<float>::event_point_type const& p1,
   typename traits<float>::event_point_type const& p2)
{
  return _event_point_ycompare<float>(p1, p2);
}

template<>
inline bool
_specialized_ycompare<
    traits<double>::event_point_type
  , traits<double>::event_point_type>
  (typename traits<double>::event_point_type const& p1,
   typename traits<double>::event_point_type const& p2)
{
  return _event_point_ycompare<double>(p1, p2);
}

template<typename Tp_>
class event_point_compare_y
{
public:
  INHERIT_TRAITS(Tp_);

private:


public:

  // When comparing two edge points A and C in edges E1 = <A,B> and E2 = <C,D>
  // the points are lexicographically sorted by following criteria:
  //   #.  A.y < C.y
  //   #.  A is right edge < C is left edge
  //   #.  A.x < C.x
  //   #.  B.y < D.y
  //   #.  B is right edge < D is left edge
  //   #.  B.x < D.x
  template<typename Pt1_, typename Pt2_>
  inline bool
  operator()(Pt1_ const& p1, Pt2_ const& p2) const
  {
    // We can do all this in one big long condition a && (b || (x && d || ...
    // but this is much more readable, and should be equivalent

    // Order by y
    if (!AlmostEqualV(gety(p1), gety(p2)))
      return gety(p1) < gety(p2);

    // Then by x
    if (!AlmostEqualV(getx(p1), getx(p2)))
      return getx(p1) < getx(p2);

    // Type-specific fall-back comparisons.
    return _specialized_ycompare(p1, p2);
  }

};

// Element of the status tree -- like an Edge, but it's not always a full
// triangle edge, e.g. when one point is formed by the intersection of two
// Edges.
template<typename Tp_>
class StatusSegment
{
public:
  INHERIT_TRAITS(Tp_);

  StatusSegment(const edge_data_type *edge_data, point_type a, point_type b)
    : edata_(edge_data), first_(a), second_(b)
  {}

  StatusSegment(const edge_type &edge)
    : edata_(edge.data()), first_(edge.first()), second_(edge.second())
  {}

  StatusSegment(const edge_type &edge, point_type const& new_top)
    : edata_(edge.data()), first_(new_top), second_(edge.second())
  {}

  edge_type const& edge(void) const { return edata_->edge(); }
  point_type const& first(void) const { return first_; }
  point_type const& second(void) const { return second_; }

#ifdef MAXTRI_DEBUG
  inline friend
  std::ostream& operator<<(std::ostream& os, StatusSegment const& s)
  {
    os << s.first() << " -> " << s.second();
    if (!AlmostEqual(s.first(), s.edge().first())
        || !AlmostEqual(s.second(), s.edge().second()))
      os << ", fragment of " << s.edge();
    else
      os << " (@tri=0x" << ioptr(s.edge().tridata()) << ")";
    return os;
  }
#endif

private:
  const edge_data_type *edata_;
  point_type first_, second_;
};


template<typename Tp_>
class compare_status
{
public:
  INHERIT_TRAITS(Tp_);

private:
  compare_status();

public:
  // We maintain a reference to the sweep line algorithm object so we can
  // dynamically sort based on the current sweep position. This is a bit sneaky
  // since std::set expects an invariable comparator, but we manually fix the
  // events at intersection points, which is the only time this would change
  // change the sorting.
  solver_type const& sweeper;

  compare_status(solver_type const& parent)
    : sweeper(parent) {}

  // Return the point on the segment which intersects the sweep line at y_s.
  static point_type
    isect_sweep(status_seg_type const& seg, coordinate_type y_s)
  {
    const coordinate_type
        x_0 = getx(seg.first()), y_0 = gety(seg.first())
      , x_1 = getx(seg.second()), y_1 = gety(seg.second());

    // derived from: y=m(x-xa)+ya where m=(yb-ya)/(xb-xa)
    coordinate_type x_s;

    // Degenerate case; segment is parallel to and intersects sweep line.
    // In this case the y value better be the same as the sweep coordinate...
    if (y_1 == y_0)
      x_s = x_0;

    else
      x_s = x_0 + ((y_s - y_0) * (x_1 - x_0) / (y_1 - y_0));

    return point_type(x_s, y_s);
  }

  // Compare two segments which share a top point by their orientation,
  // where the left-most-facing segment has lower priority matching the x-sort
  // order.
  static coordinate_type
    compare_orientation(status_seg_type const& a, status_seg_type const& b)
  {
    /* As a precondition, we assume A and B share a top point.
     * The compare result ("a less than b") is computed as follows,
     * and returned according to the rule that a lefter edge is lesser:
     *
     *                        *                 *
     *                       / \               / \
     *                   A  / r \  B       B  / r \  A
     *                     /     \           /     \
     *                    <       \         /       >
     *                   q         p       p         q
     *
     *  result?             A < B            B < A
     *  leftTurn(p,r,q)?    true             false
     *
     */
    point_type p = b.second(), r = a.first(), q = a.second();
    return orientation(p, r, q);
  }

  // Whether two segments are identical except for their triangle.
  inline bool same_segment(status_seg_type const& a, status_seg_type const& b)
    const
  {
    // I don't see a good way to factor these checks in common with the same
    // ones in operator(), so for now just duplicate the code here.
    const coordinate_type sweepy = sweeper.current_y();
    const point_type sort_a = isect_sweep(a, sweepy);
    const point_type sort_b = isect_sweep(b, sweepy);

    return AlmostEqualV(getx(sort_a), getx(sort_b))
      && AlmostEqualV(gety(sort_a), gety(sort_b))
      && AlmostEqualV(compare_orientation(a, b), coordinate_type(0));
  }

  // Compare two segments by their top point.
  // If the top point is equal it's a point of intersection, so compare by
  // orientation.
  // Note that this must dynamically use the sweep line's current y value.
  // Don't worry! We manually swap elements in the set before the comparison
  // would be invalidated, since this can only happen at intersection points.
  inline bool operator()(status_seg_type const& a, status_seg_type const& b)
    const
  {
    const coordinate_type sweepy = sweeper.current_y();

    const point_type sort_a = isect_sweep(a, sweepy);
    const point_type sort_b = isect_sweep(b, sweepy);

    // First compare by x/y value
    if (!AlmostEqualV(getx(sort_a), getx(sort_b)))
      return getx(sort_a) < getx(sort_b);

    if (!AlmostEqualV(gety(sort_a), gety(sort_b)))
      return gety(sort_a) < gety(sort_b);

    // If the points are equal, sort by orientation;
    // the left-heading edge should be less to match our x-coordinate check
    coordinate_type turn = compare_orientation(a, b);

    // If the points are parallel (same slope), compare by owning triangle
    if (abs(turn) < SMALLNUMBER)
      return a.edge().tridata() < b.edge().tridata();
    else
      return turn > 0; // left turn
  }

};

} } // end namespace cfla::tri

BOOST_GEOMETRY_REGISTER_POINT_2D_GET_SET(
    cfla::tri::EdgePoint<float>, float, cs::cartesian, x, y, x, y);
BOOST_GEOMETRY_REGISTER_POINT_2D_GET_SET(
    cfla::tri::EdgePoint<double>, double, cs::cartesian, x, y, x, y);

BOOST_GEOMETRY_REGISTER_POINT_2D_GET_SET(
    cfla::tri::IXPoint<float>, float, cs::cartesian, x, y, x, y);
BOOST_GEOMETRY_REGISTER_POINT_2D_GET_SET(
    cfla::tri::IXPoint<double>, double, cs::cartesian, x, y, x, y);

BOOST_GEOMETRY_REGISTER_POINT_2D_CONST(
    cfla::tri::TriEventPoint<float>, float, cs::cartesian, value().x(), value().y());
BOOST_GEOMETRY_REGISTER_POINT_2D_CONST(
    cfla::tri::TriEventPoint<double>, double, cs::cartesian, value().x(), value().y());

BOOST_GEOMETRY_REGISTER_POINT_2D_CONST(
    cfla::tri::EventPoint<float>, float, cs::cartesian, x(), y());
BOOST_GEOMETRY_REGISTER_POINT_2D_CONST(
    cfla::tri::EventPoint<double>, double, cs::cartesian, x(), y());

namespace cfla { namespace tri {

template<class Tp_>
class Edge
{
private:
  typedef edge_data<Tp_> m_type;
  std::shared_ptr<m_type> m_data;

public:

  // typedefs
  INHERIT_TRAITS(Tp_);

  // constructors
  Edge(triangle_data<Tp_> * parent, int edge_index,
      point_type p1, point_type p2, bp::direction_1d d)
    : m_data(new m_type(parent, edge_index, p1, p2, d))
  {}

  // methods

  inline bool operator==(Edge const& e) const {
    return (&tridata() == &tridata()) && (dir() == e.dir())
      //&& almost_equal(this->first, e.first)
      //&& almost_equal(this->second, e.second)
      ;
  }

  inline const edge_data<Tp_> *data(void) const { return m_data.get(); }
  inline edge_data<Tp_> *data(void) { return m_data.get(); }
  inline const triangle_data<Tp_> *tridata(void) const { return m_data->tri; }

  polygon_type const& poly(void) const { return tridata()->poly; }

  inline edge_point_type& first(void) { return m_data->first; }
  inline edge_point_type const& first(void) const { return m_data->first; }

  inline edge_point_type& second(void) { return m_data->second; }
  inline edge_point_type const& second(void) const { return m_data->second; }

  inline int index(void) const { return m_data->index; }

  inline bp::direction_1d dir(void) const { return m_data->dir; }
  inline void dir(bp::direction_1d d) { m_data->dir = d; }

  inline int depth(void) const { return m_data->depth; }
  inline void depth(int d) { m_data->depth = d; }

  // return the next/previous linked edge
  inline Edge const& nextEdge(void) const {
    return tridata()->edge((index() + 1) % 3); // +1 mod 3
  }
  inline Edge const& prevEdge(void) const {
    return tridata()->edge((index() + 2) % 3); // -1 mod 3
  }

#ifdef MAXTRI_DEBUG
  inline friend std::ostream& operator<<(std::ostream& os, Edge const& e) {
    os << "[ " << e.first() << " => " << e.second() << " ]"
      " (@tri=0x" << ioptr(e.tridata()) << ")";
    return os;
  }
#endif
};

} } // end namespace cfla::tri

// Trait specializations for Edge as a segment concept
namespace boost { namespace geometry { namespace traits {

  template<class Tp_>
  struct tag<cfla::tri::Edge<Tp_> > {
    typedef segment_tag type;
  };

  template<class Tp_>
  struct point_type<cfla::tri::Edge<Tp_> > {
    typedef typename cfla::tri::traits<Tp_>::edge_point_type type;
  };

  template <typename Tp_, std::size_t Dimension>
  struct indexed_access<cfla::tri::Edge<Tp_>, 0, Dimension>
  {
    typedef cfla::tri::Edge<Tp_> segment_type;
    typedef typename cfla::tri::traits<Tp_>::coordinate_type coordinate_type;
    static inline coordinate_type get(segment_type const& s) {
      return geometry::get<Dimension>(s.first);
    }
    static inline void set(segment_type& s, coordinate_type const& value) {
      geometry::set<Dimension>(s.first, value);
    }
  };

} } } // end namespace boost::geometry::traits

namespace cfla { namespace tri {

/* A triangle's vertices are always ordered counterclockwise starting at the
 * highest point. The edges are ordered lexicographically from top to bottom
 * then left to right. All edges are either left edges or right edges. The
 * bottom edge E2 is either <1,2> (left) or <2,1> (right); the second endpoint
 * of E2 is always the point at the lowest height.
 *
 *       [0]
 * (0,1)  *
 * V E0  /  \  E1 V
 *      /     \  (0,2)
 * [1] *- _    \
 *         - _   \
 *             - _\
 *        E2 >     * [2]
 *        (1,2) or (2,1)
 *
 * In this case:
 *     E0 === (0, 1) ; dir = LEFT
 *     E1 === (0, 2) ; dir = RIGHT
 *     E2 === (1, 2) ; dir = LEFT
 *
 */

template<class Tp_>
class Triangle
{
private:
  typedef triangle_data<Tp_> m_type;
  std::shared_ptr<m_type> m_data;

public:
  // typedefs
  INHERIT_TRAITS(Tp_);

  typedef typename m_type::edge_array edge_array;
  typedef typename m_type::point_array point_array;

  typedef typename point_array::iterator iterator;
  typedef typename point_array::const_iterator const_iterator;
  typedef typename point_array::size_type size_type;
  typedef typename point_array::difference_type difference_type;
  typedef typename point_array::value_type value_type;

  // constructors
  template <class Tri>
  Triangle(polygon_type const& parent, Tri const& input_tri)
    : m_data(new m_type(parent, input_tri))
  {
    // Find the highest point: this is index 0.
    int maxy_index = 0;
    event_point_ycompare ycmp;

    if (ycmp(input_tri[maxy_index], input_tri[1]))
      maxy_index = 1;
    if (ycmp(input_tri[maxy_index], input_tri[2]))
      maxy_index = 2;
    point_type maxy_point = input_tri[maxy_index];

    // Find the order of points that form a left turn.
    int left_index = (maxy_index+2) % 3; // equivalent to -1, but not negative
    int right_index = (maxy_index+1) % 3;
    if (leftTurn(input_tri[left_index], maxy_point, input_tri[right_index]))
    {
      int left_swp = left_index;
      left_index = right_index;
      right_index = left_swp;
    }

    // Now we have the point order.
    points()[0] = input_tri[maxy_index];
    points()[1] = input_tri[left_index];
    points()[2] = input_tri[right_index];

    // Find the lowest of (left, right): this is the second index of E2.
    int e2_first = 1, e2_second = 2;
    bp::direction_1d e2_dir = bp::LEFT;
    if (ycmp(points()[1], points()[2]))
    {
      e2_first = 2;
      e2_second = 1;
      e2_dir = bp::RIGHT;
    }

    // Finally we can finish the edges.
    bg::assign(edges()[0].first(), points()[0]);
    bg::assign(edges()[0].second(), points()[1]);
    bg::assign(edges()[1].first(), points()[0]);
    bg::assign(edges()[1].second(), points()[2]);
    bg::assign(edges()[2].first(), points()[e2_first]);
    bg::assign(edges()[2].second(), points()[e2_second]);
    edges()[2].dir(e2_dir);
  }

private:
  inline point_array& points(void) { return m_data->points; }
  inline edge_array& edges(void) { return m_data->edges; }

public:
  inline point_array const& points(void) const { return m_data->points; }
  inline edge_array const& edges(void) const { return m_data->edges; }
  inline polygon_type const& poly(void) const { return m_data->poly; }
  inline const triangle_data_type *data(void) const { return m_data.get(); }

  // methods
  inline edge_type const& edge(int i) const { return edges()[i % size()]; }
  inline point_type const& operator[](int i) const { return points()[i]; }
  inline size_t size(void) const { return 3u; }

  inline iterator begin() { return points().begin(); }
  inline iterator end() { return points().end(); }
  inline const_iterator begin() const { return points().cbegin(); }
  inline const_iterator end() const { return points().cend(); }
};

} } // end namespace cfla::tri

BOOST_GEOMETRY_REGISTER_RING_TEMPLATED(cfla::tri::Triangle);

namespace cfla { namespace tri {

template<class Tp_>
class SolutionCell
{
public:
  INHERIT_TRAITS(Tp_);

  // triangles which form the maximal intersection
  typedef typename std::unordered_set<const triangle_type*> triangle_ref_set;

  typedef typename triangle_ref_set::const_iterator iterator;
  typedef typename triangle_ref_set::size_type size_type;

private:
  triangle_ref_set source_tris_;
  int depth_;

public:
  // A solution cell must be formed from at least two triangles.
  SolutionCell(triangle_type const& t1, triangle_type const& t2)
    : source_tris_({ &t1, &t2 }), depth_(-1)
  {}

  inline int depth(void) const { return depth_; }

  inline iterator begin(void) const { return source_tris_.cbegin(); }
  inline iterator end(void) const { return source_tris_.cend(); }
  inline size_type size(void) const { return source_tris_.size(); }

  // Intersect the input triangles to form the final solution cell.
  solution_cell_type cell(void) const
  {
    solution_cell_type out;
    if (size() == 0)
      return out;
    intersection(source_tris_.begin(), source_tris_.end(), out);
    return out;
  }

};

template<typename Tp_>
struct make_sla_traits
{
  INHERIT_TRAITS(Tp_);

  typedef typename traits::polygon_type *value_type; // isochromes
  typedef typename traits::solution_ref_type solution_type;

  typedef typename traits::event_point_type event_type;
  typedef typename traits::event_point_ycompare event_compare;
  typedef typename std::vector<event_type> event_container;
  typedef EventType       event_id_type;
  typedef make_sla_traits etraits;

  typedef sla_traits<value_type, event_type, solution_type, etraits> type;

  inline static EventType get_type(event_type const& e) {
    return e.type();
  }
};

template<class Tp_>
class MaxTri
  : public SweepLineAlgorithm<typename make_sla_traits<Tp_>::type>
{
public:
  typedef typename make_sla_traits<Tp_>::type sla;
  typedef SweepLineAlgorithm<sla> super_type;

  INHERIT_TRAITS(Tp_);

  typedef typename sla::value_type value_type;
  typedef std::vector<std::unique_ptr<polygon_type> > poly_container;
  typedef std::vector<triangle_type> triangle_container;
  typedef std::set<status_seg_type, status_compare> event_status;

  typedef typename event_status::iterator status_iterator;
  typedef typename event_status::reverse_iterator status_riterator;

  typedef typename super_type::solution_type solution_type;
  typedef typename super_type::solution_container solution_container;
  typedef typename solution_container::const_iterator solution_iterator;

private:
  triangle_container tris;
  poly_container polys;

  // sweep line status: points of current edges, ordered by X coordinate
  event_status status;
  int max_depth;
  coordinate_type _lasty;

  void insert_edge(edge_type const& e);
  void remove_edge(edge_type const& e);

  void handle_intersection(isect_point_type const& event);
  void handle_tripoint(tri_point_type const& p);

  using super_type::solutions;
  using super_type::queue;

  // Quick check for whether the ray given by the directed
  // segment (first->second) intersects the test segment.
  // This means the line segment e *might* intersect it.
  // This assumes already that the highest point of test is in the proper
  // direction. This is used for the loop condition in insert_edge(), so we
  // already know that test is sorted by x in the correct direction.
  inline
  bool can_intersect(status_seg_type const& seg, status_seg_type const& test)
  {
    bool test_dir1 = getx(seg.second()) - getx(test.first()) < 0;
    bool test_dir2 = getx(seg.second()) - getx(test.second()) < 0;

    // One point must be to the left, the other to the right.
    return test_dir1 != test_dir2;
  }

  // Look through the status in both directions to find which edges intersect
  // the given edge. Queue the intersection points.
  void check_intersections(const status_iterator center);

  // Like check_intersections, but only search in the direction of ++edge.
  template<typename Iter>
  void check_intersection(const Iter edge, Iter neighbor, const Iter end);

  // Leftwards check_intersection: convert arguments to reverse iterators
  // and pass to check_intersection above. Disambiguated by a reverse end()
  // iterator.
  inline void check_intersection(status_iterator edge, status_iterator left,
      status_riterator rend) {
    check_intersection(status_riterator(edge), status_riterator(left), rend);
  }

  // Intersect the two line segments defined by Edges.
  // If there is an intersection (return code other than '0'), sets the fields
  // of ix appropriately.
  char intersect(edge_type const& e1, edge_type const& e2, isect_point_type &ix)
    const;

#ifdef MAXTRI_DEBUG
  void dump_status_to_octave(std::ostream& os);
#endif

protected:
  // override
  void handle_event(EventType type, event_point_type const& event);
  void initialize(void);
  void finalize(void);

public:
  ~MaxTri();

  // No inputs yet.
  MaxTri()
    : tris(), status(status_compare(*this)), max_depth(-1), _lasty()
  {}

  // Construct from a list of input polygons (may be in no particular order).
  template<class Iter>
  MaxTri(Iter begin, Iter end)
    : tris(), status(status_compare(*this)), max_depth(-1), _lasty()
  {
    super_type::insert(begin, end);
  }

  coordinate_type current_y(void) const {
    return _lasty;
  }

  // Add triangles from polygon.
  // Same as the constructor form, if you're lazy and want to do it
  // after construction. (Overridden from SLA.)
  void add_event(value_type const& ply)
  {
    polys.emplace_back(ply);
    if (ply->triangles_size() == 0u)
      throw std::runtime_error("must triangulate the polygon first!");

    for (auto ptri = ply->triangles_begin(); ptri != ply->triangles_end();
        ++ptri)
    {
      point_type triangle[3] = {
        (*ply)[ptri->v[0]]->getPos(),
        (*ply)[ptri->v[1]]->getPos(),
        (*ply)[ptri->v[2]]->getPos(),
      };

      // Avoid filthy degenerates, as they will gain us nothing
      if (Collinear(triangle[0], triangle[1], triangle[2]))
        continue;

      tris.emplace_back(*ply, triangle);

#ifdef MAXTRI_DEBUG
      std::cout << "QUEUE triangle" << std::endl
        << "    [0] " << tris.back().edge(0) << std::endl
        << "    [1] " << tris.back().edge(1) << std::endl
        << "    [2] " << tris.back().edge(2) << std::endl;
#endif

      const triangle_data_type *const tri_data = tris.back().data();
      queue().emplace(tri_point_type(tri_data, 0));
      queue().emplace(tri_point_type(tri_data, 1));
      queue().emplace(tri_point_type(tri_data, 2));
    }
  }

  inline int depth(void) const { return max_depth; }

  // Return the solution cells.
  inline size_t size(void) const { return solutions().size(); }
  inline solution_iterator begin(void) const { return solutions().cbegin(); }
  inline solution_iterator end(void) const { return solutions().cend(); }

  inline solution_type const& solution(int idx) const { return *(begin()+idx); }
};

/* Adapter for MaxTri which conforms to the cfla_traits::solver_type
 * interface. That is, defines Point operator()(facilities...) and
 * iterator methods that iterate over input user points.
 */
template<class Tp_>
class MaxTriSolver
{
public:
  INHERIT_TRAITS(Tp_);

  typedef User<point_type> user_type;
  typedef std::vector<user_type> user_list;
  typedef typename user_list::iterator       iterator;
  typedef typename user_list::const_iterator const_iterator;

  typedef TTNN1<user_type, MaxTriSolver> nn1_type;
  typedef typename nn1_type::size_type size_type;

private:
  solver_type solver_;
  user_list users_;

  // Return the ideal location for player 2 among some player 1 facilities.
  template<typename point_filter>
  point_type compute(const nn1_type &nn1, point_filter filter)
  {
    // To solve for player 2, build isochromes from the opposing player's
    // facilities to each user point, then triangulate the polygons (for
    // simplicity) and solve for max-depth triangular intersections.
    for (auto userp = begin(); userp != end(); ++userp)
    {
      // Only build rects for certain points.
      if (!filter(*userp))
        continue;

      // Compute an isochrome representing the fixed travel time distance to
      // the nearest facility.
      point_type nearest_facility(nn1(*userp));

      polygon_type *p = new polygon_type;
      userp->isochrome(*p, nearest_facility);
      solver_.add_event(p);
    }
    solver_.compute();

    // We have potentially multiple solution cells.
    if (solver_.size() == 0)
    {
      std::cerr << "warning: empty solution" << std::endl;
      return point_type();
    }

    // Choose a random point in a random solution cell.
    size_type chosen_one = randrange(size_type(0), solver_.size()-1);
    auto solution = solver_.solution(chosen_one).cell();

    bgm::box<point_type> bbox;
    bg::envelope(solution, bbox);
    point_type out;
    coordinate_type randx, randy;
    do {
      randx = randrange(bg::get<bg::min_corner, 0>(bbox),
                        bg::get<bg::max_corner, 0>(bbox));
      randy = randrange(bg::get<bg::min_corner, 1>(bbox),
                        bg::get<bg::max_corner, 1>(bbox));
      out = point_type(randx, randy);
    } while (!bg::within(out, solution));

    return out;
  }

public:

  template<typename UserIter>
    MaxTriSolver(UserIter begin, UserIter end)
      : solver_(), users_(begin, end)
    {}

  const_iterator begin(void) const { return users_.begin(); }
  const_iterator end(void) const { return users_.end(); }
  //const_iterator cbegin(void) const { return users_.cbegin(); }
  //const_iterator cend(void) const { return users_.cend(); }

  // Return the "optimal solution" for the CFL game.
  template<typename point_filter>
  point_type operator()(const nn1_type &facilities, point_filter filter) {
    return compute(facilities, filter);
  }

  static point_type user_point(user_type const& u) {
    return u.center();
  }

  inline static coordinate_type
    distance(user_type const& user, point_type const& p) {
      return user.travelTime(p);
    }

};

extern template class MaxTri<double>;
extern template class MaxTri<float>;

#ifdef MAXTRI_DEBUG
template<typename U>
std::ostream& operator<<(std::ostream& os, Edge<U> const& e) {
  os << "<[" << std::setw(2) << std::setfill(' ') << e.rect_index
    << "] " << (e.dir == bp::LEFT ? "LEFT " : "RIGHT")
    << " " << e.coord << " d=" << e.depth << ">";
  return os;
}

template<typename U>
std::ostream& operator<<(std::ostream& os, SolutionCell<U> const& e)
{
  os << "solution (depth " << e.depth << ", " << e.size() << ")" << std::endl;
  unsigned int idx = 0u;
  for (const auto& sol : e)
  {
    os << "    [" << std::setw(2) << std::setfill(' ') << idx << "] "
      << sol << std::endl;
    idx++;
  }
  return os;
}
#endif

#undef INHERIT_TRAITS

} } // end namespace cfla::tri
