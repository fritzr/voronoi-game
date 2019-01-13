#pragma once

#include <set>
#include <unordered_map>
#include <queue>
#include <vector>
#include <functional>
#include <boost/array.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp> // point_xy, ring
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/ring.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/algorithms/is_convex.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/register/ring.hpp>

#if BOOST_VERSION == 106000
// https://svn.boost.org/trac10/ticket/11880
// adjacency_matrix.hpp uses ice_not after it has been deprecated
// from type_traits.hpp
#include <boost/type_traits/ice.hpp>
#endif

#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/adjacency_matrix.hpp>

#if BOOST_VERSION < 106700
#include <boost/functional/hash.hpp>
#else
#include <boost/container_hash/hash.hpp>
#endif

#include "util.h"
#include "user.h"
#include "polygon.h"
#include "intersection.h"
#include "sweep.h"
#include "nn1.h"

#ifdef MAXTRI_DEBUG

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

// With debug enabled, don't inline stuff so it's easier to debug.
#define MTINLINE

#define MTFAIL(...) do { \
  std::stringstream ss; \
  ss << __PRETTY_FUNCTION__ << ": " << __VA_ARGS__; \
  throw std::runtime_error(ss.str()); \
} while(0)

#define MTDEBUG(...) __VA_ARGS__

#else

#define MTINLINE inline
#define MTFAIL(...)
#define MTDEBUG(...)

#endif // MAXTRI_DEBUG

#ifdef MAXTRI_DEBUG
struct ioptr
{
  const uintptr_t ptrval;

  ioptr(const void *p) : ptrval(reinterpret_cast<uintptr_t>(p)) {}

  MTINLINE friend
  std::ostream &operator<<(std::ostream &os, ioptr const& iop)
  {
    os << "0x" << std::hex << std::setw(2*sizeof(iop.ptrval))
      << std::setfill('0') << iop.ptrval;
    return os;
  }
};
#endif // MAXTRI_DEBUG

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

      static MTINLINE coordinate_type get(const cfla::tri::point_xy<Tp_> &p,
          orientation_2d orient)
        { return (orient == HORIZONTAL) ? p.x() : p.y(); }
    };

  template<typename Tp_> struct
    point_mutable_traits<cfla::tri::point_xy<Tp_> >
    {
      typedef typename bg::traits::coordinate_type<cfla::tri::point_xy<Tp_> >
        ::type coordinate_type;

      static MTINLINE void set(cfla::tri::point_xy<Tp_> &p,
          orientation_2d orient, coordinate_type value) {
        switch (orient.to_int()) {
          case VERTICAL: p.y(value); break;
          case HORIZONTAL: p.x(value); break;
        }
      }

      static MTINLINE cfla::tri::point_xy<Tp_>
        construct(coordinate_type x, coordinate_type y) {
          return cfla::tri::point_xy<Tp_> (x, y);
        }
    };

} } // end namespace boost::polygon

namespace boost {
  template <typename Coord>
  struct hash<cfla::tri::point_xy<Coord> >
  {
    inline size_t operator()(cfla::tri::point_xy<Coord>  const& p) const
    {
      std::size_t hval = 0u;
      hash_combine(hval, geometry::get<0>(p));
      hash_combine(hval, geometry::get<1>(p));
      return hval;
    }
  };

} // end namespace boost

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
struct compare_status;

template<typename Tp_>
struct compare_segment_unique;

class event_point_compare_y;

template<typename Tp_>
class TriEventPoint;

template<typename Tp_>
class IXEvent;

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
  typedef IXEvent<Tp_>                   isect_event_type;
  typedef EventPoint<Tp_>                event_point_type;
  typedef event_point_compare_y          event_point_ycompare;

  typedef StatusSegment<Tp_>             status_seg_type;
  typedef compare_status<Tp_>            status_compare;
  typedef compare_segment_unique<Tp_>    status_equal_compare;

  typedef c_polygon<point_type>          polygon_type;
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
  typedef typename traits::isect_event_type isect_event_type; \
  typedef typename traits::event_point_type event_point_type; \
  typedef typename traits::event_point_ycompare event_point_ycompare; \
  \
  typedef typename traits::status_seg_type status_seg_type; \
  typedef typename traits::status_compare status_compare; \
  typedef typename traits::status_equal_compare status_equal_compare; \
  \
  typedef typename traits::polygon_type polygon_type; \
  typedef typename traits::solution_cell_type solution_cell_type; \
  typedef typename traits::solver_type solver_type; \

// functions

using ::operator<<;

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
    bg::assign(out, first);
    return true;
  }

  // Case B. two geometries, return their intersection
  auto const& second = *begin++;
  if (!bg::intersection(first, second, out))
    return false;

  // Case C. many geometries, return the cumulative intersection with A&B
  while (begin != end)
  {
    // XXX is providing out twice okay??
    if (!bg::intersection(*begin++, out, out))
      return false;
  }

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

// These are used to describe triangle points, and follow the definitions in
// Triangle.
enum PointHeight {
    TOP=0
  , MIDDLE
  , BOTTOM
};

template<class Tp_>
struct triangle_data
{
  INHERIT_TRAITS(Tp_);

  typedef edge_type edge_array[3];
  typedef std::array<point_type, 3> point_array;

  typedef typename point_array::iterator iterator;
  typedef typename point_array::const_iterator const_iterator;
  typedef typename point_array::size_type size_type;
  typedef typename point_array::difference_type difference_type;
  typedef typename point_array::value_type value_type;

  inline iterator begin(void) { return points.begin(); }
  inline iterator end(void) { return points.end(); }
  inline const_iterator begin(void) const { return points.begin(); }
  inline const_iterator end(void) const { return points.end(); }
  inline const_iterator cbegin(void) const { return points.cbegin(); }
  inline const_iterator cend(void) const { return points.cend(); }

  // members
  int id; // index of the triangle in the triangle list; serves as unique ID
  polygon_type const& poly; // parent polygon
  point_array points;
  edge_array edges;
  int middle_point_index; // whether the middle point is point 1 or 2

  template <class Tri>
  triangle_data(polygon_type const& parent, Tri const& input_tri, int idx)
    : id(idx), poly(parent), points()
      , edges({ // to be filled in by Triangle constructor
          { this, 0, input_tri[0], input_tri[0], bp::RIGHT },
          { this, 1, input_tri[0], input_tri[0], bp::LEFT },
          { this, 2, input_tri[0], input_tri[0], bp::LEFT },
        })
      , middle_point_index(1)
  {
  }

  // parent edge
  MTINLINE edge_type const& edge(int index) const {
#ifdef MAXTRI_DEBUG_INTERSECT
    if (index < 0)
      MTFAIL("negative edge index!");
#endif
    return edges[index % 3];
  }

  MTINLINE point_type const& point(int index) const {
#ifdef MAXTRI_DEBUG_INTERSECT
    if (index < 0)
      MTFAIL("negative point index!");
    if (index > 2)
      MTFAIL("out of bounds index!");
#endif
    return points[index];
  }

  MTINLINE point_type const& top(void) const { return point(TOP); }
  MTINLINE point_type const& middle(void) const { return point(MIDDLE); }
  MTINLINE point_type const& bottom(void) const { return point(BOTTOM); }

  MTINLINE int top_index(void) const { return 0; }
  MTINLINE int middle_index(void) const { return middle_point_index; }
  MTINLINE int bottom_index(void) const { return (middle_point_index * 2) % 3; }
  /* middle_point_index is either 1 or 2.
   * multiplying by 2 mod 3 maps 1 <-> 2, so we get the bottom index
   * without using a conditional:
   *  0 * 2 = 0, 0 % 3 = 0
   *  1 * 2 = 2, 2 % 3 = 2
   *  2 * 2 = 4, 4 % 3 = 1
   */
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
  MTINLINE edge_type const& edge(void) const {
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
  MTINLINE edge_type const& edge(void) const {
    return edata->tri->edges[edata->index];
  }

  // return the point at the other end of this edge
  MTINLINE edge_point_type const& other(void) const {
    return (dir == bp::LOW) ? edata->second : edata->first;
  }

  MTINLINE const triangle_data_type *tridata(void) const {
    return edata->tri;
  }

#ifdef MAXTRI_DEBUG
  MTINLINE friend std::ostream& operator<<(std::ostream& os, EdgePoint const& p) {
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
  PointHeight height;

  TriEventPoint()
    : tridata(nullptr), index(-1), height(TOP)
  {}

  TriEventPoint(const triangle_data_type *tdata, int point_index, PointHeight h)
    : tridata(tdata), index(point_index), height(h)
  {}

  MTINLINE point_type const& value(void) const {
    return tridata->point(index);
  }

 /* From the perspective of the algorithm, we need to visit the edges in the
  * correct order. For the top point we insert both edges, and for the bottom
  * point we remove both edges. Since the operation is the same on both edges
  * it doesn't matter which is which. But for the middle point (which can be
  * point [1] or [2]) we need to remove the edge we inserted from [0], and
  * insert the other edge which will be removed by the bottom point.
  * Therefore we define a TOP edge and BOTTOM edge for each point with respect
  * to the point indexes as defined by Triangle:
  *
  *       [0]               [0]
  *        *                     top edge: E0
  * V E0  /  \  E1 V          bottom edge: E1
  *      /     \            [1]
  * [1] *- _    \                top edge: E0
  *         - _   \           bottom edge: E2
  *             - _\        [2]
  *        E2 >    * [2]         top edge: E1
  *                           bottom edge: E2
  */

  MTINLINE edge_type const& top_edge(void) const
  {
    static const int top_edge_indexes[] = {
      /* point index -> top edge index */
      /* 0 -> */ 0,
      /* 1 -> */ 0,
      /* 2 -> */ 1,
    };
    return tridata->edges[top_edge_indexes[index]];
  }

  MTINLINE edge_type const& bottom_edge(void) const
  {
    static const int bottom_edge_indexes[] = {
      /* point index -> bottom edge index */
      /* 0 -> */ 1,
      /* 1 -> */ 2,
      /* 2 -> */ 2,
    };
    return tridata->edges[bottom_edge_indexes[index]];
  }

#ifdef MAXTRI_DEBUG
  MTINLINE friend
  std::ostream& operator<<(std::ostream& os, TriEventPoint const& p)
  {
    os << "[" << p.index << "] joining "
      << p.top_edge() << " and " << p.bottom_edge();
    return os;
  }
#endif
};

/* Intersection point.
 *
 * An intersection point is classically the intersection of TWO segments...
 * But we need to handle many segments intersecting at the same point.
 * Generally this sounds rare, but actually in the case of several collinear
 * (parallel overlapping) segments, a single independent segment that
 * intersects any of them would otherwise form a pairwise intersection with
 * each of them, resulting in many separate intersection events at the same
 * point.
 *
 * Here is an example where B and C are collinear, and both intersect A at P.
 * Let's say B comes before C by nature of the order in which it was inserted.
 * Classically, the status tree progresses as shown on the right:
 *
 *    A      B C
 *     \     /
 *       \  /      L1. A B C    ISECT A,B: SWAP A,B
 *      P  X       L2. B A C    ISECT A,C: SWAP A,C
 *        /  \     L3. B C A    ...
 *       /     \
 *
 * However, there is no way to represent L2. without directly accessing
 * the status tree nodes; in C++ we would have to modify the comparator so that
 * we could remove A and insert it again, and the set would order A after
 * B but before C. This is a problem because B and C are the same segment as
 * far as the comparator is concerned. Thus, any such collection of collinear
 * segments like B, C must be handled atomically with respect to some other
 * intersecting segment like A. That is, the tree must transition directly
 * from state L1. to L3. in our example above.
 *
 * This means each intersection point needs to track a set of segments!
 * To do this we have to keep a[n unordered] map of intersection sets keyed
 * on position for fast lookup. When we identify that two segments intersect,
 * we have to lookup the intersection set by position and add to it, then
 * handle the intersection exactly once.
 */

// Specialization of hash for StatusSegment: hash first based on endpoints
// then on triangle data.
template<typename Tp_>
struct hash_segment
{
  INHERIT_TRAITS(Tp_);

  typedef status_seg_type argument_type;
  typedef size_t result_type;

  hash_segment() {}

  inline bool operator()(argument_type const& seg) const
  {
    size_t hval = 0u;
    ::boost::hash_combine(hval, seg.first());
    ::boost::hash_combine(hval, seg.second());
    ::boost::hash_combine(hval, seg.edge().tridata());
    return hval;
  }
};

template<typename Tp_>
class IXEvent : public point_xy<Tp_>
{
public:
  INHERIT_TRAITS(Tp_);

  typedef point_xy<Tp_> super;

  // Tracks the segments which form this point. They cannot be in any
  // particular order until the instant the intersection point is handled.  The
  // only reason the main status tree can do this is because the order is
  // corrected whenever segments are swapped due to an intersection. From this
  // class' perspective, we don't know if one of our segments is swapped with
  // another one due to a different intersection point, so we can't keep our
  // set up-to-date. That's okay -- we just don't want to have any segment
  // added multiple times. Thus, we use an unordered_set.
  typedef std::unordered_set<
    /* Key */ status_seg_type
    /* Hash */, hash_segment<Tp_>
    /* KeyEqual */, status_equal_compare
    > segment_set;

  typedef typename segment_set::iterator iterator;
  typedef typename segment_set::const_iterator const_iterator;
  typedef typename segment_set::size_type size_type;

  using super::x;
  using super::y;

  IXEvent()
    : super(), segments_(status_equal_compare())
  {}

  ~IXEvent()
  {}

  IXEvent(IXEvent const& other)
    : super(getx(other), gety(other)), segments_(other.segments_)
  {}

  // Construct from point.
  template<typename Pt_>
  IXEvent(Pt_ const& point)
    : super(point), segments_()
  {}

  // Record that two segments intersect at this point.
  inline void insert(status_seg_type const& a, status_seg_type const& b)
  {
    segments_.insert(a);
    segments_.insert(b);
  }

  // Methods.

  inline segment_set const& segments(void) const { return segments_; }
  inline const_iterator begin(void) const { return segments().cbegin(); }
  inline const_iterator end(void) const { return segments().cend(); }
  inline const_iterator cbegin(void) const { return segments().cbegin(); }
  inline const_iterator cend(void) const { return segments().cend(); }
  inline size_type size(void) const { return segments().size(); }

#ifdef MAXTRI_DEBUG
  MTINLINE friend
  std::ostream& operator<<(std::ostream& os, IXEvent const& p)
  {
    os << "intersection at " << static_cast<super const&>(p) << " from:"
      << std::endl;
    int index = 0;
    for (auto const& segment : p.segments())
    {
      os << "    [" << std::setw(2) << std::setfill(' ') << index++
        << "] " << segment << std::endl;
    }
    return os;
  }
#endif

private:
  segment_set segments_;
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
    point_type isect;
  } u;

  EventType type;

  event_data(tri_point_type const& e)
    : u{.point=e}, type(TRIPOINT) {}

  event_data(point_type const& ix)
    : u{.isect=ix}, type(INTERSECTION) {}

  ~event_data()
  {
    switch (type)
    {
    case TRIPOINT:
      (u.point).~TriEventPoint();
    case INTERSECTION:
      (u.isect).~point_type();
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

  EventPoint(IXEvent<Tp_> const& ipoint)
    : m_union(new m_type(ipoint)) {}

  MTINLINE EventType type(void) const {
    return m_union->type;
  }

  MTINLINE tri_point_type const& point(void) const {
#ifdef MAXTRI_DEBUG_INTERSECT
    if (type() != TRIPOINT)
      MTFAIL("type mismatch, got " << type() << ", expected TRIPOINT");
#endif
    return m_union->u.point;
  }

  MTINLINE point_type const& isect(void) const {
#ifdef MAXTRI_DEBUG_INTERSECT
    if (type() != INTERSECTION)
      MTFAIL("type mismatch, got " << type() << ", expected INTERSECTION");
#endif
    return m_union->u.isect;
  }

  MTINLINE coordinate_type x(void) const
  {
    switch (type())
    {
    case TRIPOINT:
      return getx(point().value());
    case INTERSECTION:
      return getx(isect());
    default:
      // should be unreachable!
      return -1;
    }
  }

  MTINLINE coordinate_type y(void) const
  {
    switch (type())
    {
    case TRIPOINT:
      return gety(point().value());
    case INTERSECTION:
      return gety(isect());
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

class event_point_compare_y
{
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
  MTINLINE bool
  operator()(Pt1_ const& p1, Pt2_ const& p2) const
  {
    // We can do all this in one big long condition a && (b || (x && d || ...
    // but this is much more readable, and should be equivalent

    // Order by y
    if (!AlmostEqualV(gety(p1), gety(p2)))
      return gety(p1) < gety(p2);

    // Then by x (in reverse)
    if (!AlmostEqualV(getx(p1), getx(p2)))
      return getx(p2) < getx(p1);

    // Type-specific fall-back comparisons.
    return _specialized_ycompare(p1, p2);
  }

};

// Element of the status tree -- reference to an Edge.
template<typename Tp_>
class StatusSegment
{
public:
  INHERIT_TRAITS(Tp_);

  StatusSegment(const StatusSegment& other)
    : edata_(other.edge().data())
  {}

  StatusSegment(const edge_type &edge)
    : edata_(edge.data())
  {}

  edge_type const& edge(void) const { return edata_->edge(); }
  point_type const& first(void) const { return edge().first(); }
  point_type const& second(void) const { return edge().second(); }

#ifdef MAXTRI_DEBUG
  MTINLINE friend
  std::ostream& operator<<(std::ostream& os, StatusSegment const& s)
  {
    os << "&" << s.edge();
    return os;
  }
#endif

private:
  const edge_data_type *edata_;
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

    coordinate_type x_s;

    // Degenerate case; segment is parallel to and intersects sweep line.
    // In this case the y value better be the same as the sweep coordinate...
    if (AlmostEqualV(y_1, y_0))
      x_s = x_0;

    // derived from: y=m(x-xa)+ya where m=(yb-ya)/(xb-xa)
    else
      x_s = x_0 + ((y_s - y_0) * (x_1 - x_0) / (y_1 - y_0));

    return point_type(x_s, y_s);
  }

  /* Returns the orientation from the end of B through a pivot point where A
   * and B intersect to the end of A.
   *
   * The orientation value is such that right turn < parallel =~ 0 < left turn:
   *
   *          *   <- pivot ->   *
   *         / \               / \
   *     A  / q \  B       B  / q \  A
   *       /     \           /     \
   *      <       \         /       >
   *     r         p       p         r
   *
   *      left turn        right turn
   *        A < B             B < A
   */
  static coordinate_type
    compare_orientation(status_seg_type const& a, status_seg_type const& b,
        point_type pivot)
  {
    // As a precondition, we assume A and B share the pivot point.
    point_type p = b.second(), r = a.second();
    return orientation(p, pivot, r);
  }

  // Compare two segments by their intersection with the sweep line.
  // At intersection points we compare by orientation.
  //
  // Note that this must dynamically use the sweep line's current y value.
  // Don't worry! The order only changes at intersection points, at which time
  // we remove the segments, change the ordering (sweep line position), and
  // re-insert the only segments which could have been invalidated (the
  // segments forming the intersection).
  MTINLINE bool operator()(status_seg_type const& a, status_seg_type const& b)
    const
  {
    const coordinate_type sweepy = sweeper.current_y();

    const point_type sort_a = isect_sweep(a, sweepy);
    const point_type sort_b = isect_sweep(b, sweepy);

#ifdef MAXTRI_DEBUG_INTERSECT
    std::cerr
      << "        comparing " << sort_a << " | " << a << std::endl
      << "                < " << sort_b << " | " << b << "..." << std::endl
      << "            -> ";
#endif

    // First compare by x/y value
    if (!AlmostEqualV(getx(sort_a), getx(sort_b)))
    {
      bool ret = getx(sort_a) < getx(sort_b);
#ifdef MAXTRI_DEBUG_INTERSECT
      std::cerr << ret << ":"
        << " x_a (" << getx(sort_a) << ") " << (ret ? " <" : ">=")
        << " x_b (" << getx(sort_b) << ")" << std::endl;
#endif
      return ret;
    }

    if (!AlmostEqualV(gety(sort_a), gety(sort_b)))
    {
      bool ret = gety(sort_a) < gety(sort_b);
#ifdef MAXTRI_DEBUG_INTERSECT
      std::cerr << ret << ":"
        << " y_a (" << gety(sort_a) << ") " << (ret ? " <" : ">=")
        << " y_b (" << gety(sort_b) << ")" << std::endl;
#endif
      return ret;
    }

    // If the points are equal, sort by orientation;
    coordinate_type turn = compare_orientation(a, b, sort_a);

    // If the segments are collinear (parallel + overlap) they are unsortable.
    // They will be put in the same equal_range in the status multiset and
    // must be dealt with specially.
    // Equals means we should return false associatively.
    if (AlmostEqualV(turn, 0))
    {
#ifdef MAXTRI_DEBUG_INTERSECT
      std::cerr << false << ": collinear (equal)" << std::endl;
#endif
      return false;
    }

    // The left-heading edge should be less so it is first in a left-to-right
    // (low-to-high) traversal of the sweep status.
    bool ret = turn > 0; // true for left turn
#ifdef MAXTRI_DEBUG_INTERSECT
    std::cerr << ret << ": " << (ret ? " left" : "right") << " turn "
      << b.second() << " -> " << sort_a << " -> " << a.second() << std::endl;
#endif
    return ret;
  }

};

// This is used as KeyEqual for unordered sorting of segments in the IXEvent
// segments set.  All we care about here is to uniquely identify segments.
template<typename Tp_>
struct compare_segment_unique
{
  INHERIT_TRAITS(Tp_);

  inline bool operator()(status_seg_type const& a, status_seg_type const& b)
    const
  {
    // Both points have to be equal of course...
    // Then disambiguate by the owning triangle.
    return AlmostEqual(a.first(), b.first())
      && AlmostEqual(a.second(), b.second())
      && a.edge().tridata() == b.edge().tridata();
  }
};

} } // end namespace cfla::tri

BOOST_GEOMETRY_REGISTER_POINT_2D_GET_SET(
    cfla::tri::EdgePoint<float>, float, cs::cartesian, x, y, x, y);
BOOST_GEOMETRY_REGISTER_POINT_2D_GET_SET(
    cfla::tri::EdgePoint<double>, double, cs::cartesian, x, y, x, y);

BOOST_GEOMETRY_REGISTER_POINT_2D_GET_SET(
    cfla::tri::IXEvent<float>, float, cs::cartesian, x, y, x, y);
BOOST_GEOMETRY_REGISTER_POINT_2D_GET_SET(
    cfla::tri::IXEvent<double>, double, cs::cartesian, x, y, x, y);

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

  MTINLINE bool operator==(Edge const& e) const {
    return (&tridata() == &tridata()) && (dir() == e.dir())
      //&& almost_equal(this->first, e.first)
      //&& almost_equal(this->second, e.second)
      ;
  }

  MTINLINE const edge_data<Tp_> *data(void) const { return m_data.get(); }
  MTINLINE edge_data<Tp_> *data(void) { return m_data.get(); }
  MTINLINE const triangle_data<Tp_> *tridata(void) const { return m_data->tri; }

  polygon_type const& poly(void) const { return tridata()->poly; }

  MTINLINE edge_point_type& first(void) { return m_data->first; }
  MTINLINE edge_point_type const& first(void) const { return m_data->first; }

  MTINLINE edge_point_type& second(void) { return m_data->second; }
  MTINLINE edge_point_type const& second(void) const { return m_data->second; }

  MTINLINE int index(void) const { return m_data->index; }

  MTINLINE bp::direction_1d dir(void) const { return m_data->dir; }
  MTINLINE void dir(bp::direction_1d d) { m_data->dir = d; }

  MTINLINE int depth(void) const { return m_data->depth; }
  MTINLINE void depth(int d) { m_data->depth = d; }

  // return the next/previous linked edge
  MTINLINE Edge const& nextEdge(void) const {
    return tridata()->edge((index() + 1) % 3); // +1 mod 3
  }
  MTINLINE Edge const& prevEdge(void) const {
    return tridata()->edge((index() + 2) % 3); // -1 mod 3
  }

#ifdef MAXTRI_DEBUG_INTERSECT
  MTINLINE friend std::ostream& operator<<(std::ostream& os, Edge const& e) {
    os << "[ " << e.first() << " => " << e.second() << " ]"
      " (@tri=" << ioptr(e.tridata()) << ")";
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
    static MTINLINE coordinate_type get(segment_type const& s) {
      return geometry::get<Dimension>(s.first);
    }
    static MTINLINE void set(segment_type& s, coordinate_type const& value) {
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
 *       [0]                TOP
 * (0,1)  *
 * V E0  /  \  E1 V
 *      /     \  (0,2)
 * [1] *- _    \            MIDDLE
 *         - _   \
 *             - _\
 *        E2 >     * [2]    BOTTOM
 *        (1,2) or (2,1)
 *
 * In this case:
 *     E0 === (0, 1) ; dir = LEFT
 *     E1 === (0, 2) ; dir = RIGHT
 *     E2 === (1, 2) ; dir = LEFT
 *
 * Additionally, we describe the vertices as "top", "middle", and "bottom".
 * Top is always point[0]; then middle and bottom are either [1] or [2]
 * depending on which is higher. These follow the event queue sort order.
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

  bg::order_selector point_order = bg::counterclockwise;

  typedef typename point_array::iterator iterator;
  typedef typename point_array::const_iterator const_iterator;
  typedef typename point_array::size_type size_type;
  typedef typename point_array::difference_type difference_type;
  typedef typename point_array::value_type value_type;

  // constructors
  template <class Tri>
  Triangle(polygon_type const& parent, Tri const& input_tri, int index)
    : m_data(new m_type(parent, input_tri, index))
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
    bp::direction_1d e2_side = bp::LEFT;
    if (ycmp(points()[1], points()[2]))
    {
      e2_first = 2;
      e2_second = 1;
      e2_side = bp::RIGHT;
    }

    m_data->middle_point_index = e2_first;

    // Finally we can finish the edges.
    bg::assign(edges()[0].first(), points()[0]);
    bg::assign(edges()[0].second(), points()[1]);
    bg::assign(edges()[1].first(), points()[0]);
    bg::assign(edges()[1].second(), points()[2]);
    bg::assign(edges()[2].first(), points()[e2_first]);
    bg::assign(edges()[2].second(), points()[e2_second]);
    edges()[2].dir(e2_side);
  }

private:
  MTINLINE point_array& points(void) { return m_data->points; }
  MTINLINE edge_array& edges(void) { return m_data->edges; }

public:
  MTINLINE point_array const& points(void) const { return m_data->points; }
  MTINLINE edge_array const& edges(void) const { return m_data->edges; }
  MTINLINE polygon_type const& poly(void) const { return m_data->poly; }
  MTINLINE const triangle_data_type *data(void) const { return m_data.get(); }

  MTINLINE point_type const& point(int i) const { return m_data->point(i); }
  // Top, middle, and bottom points as described above.
  MTINLINE point_type const& point(PointHeight h) const {
    switch(h) {
      case BOTTOM: return bottom();
      case MIDDLE: return middle();
      default:
      case TOP: return top();
    }
  }

  MTINLINE point_type const& top(void) const { return m_data->top(); }
  MTINLINE point_type const& middle(void) const { return m_data->middle(); }
  MTINLINE point_type const& bottom(void) const { return m_data->bottom(); }

  MTINLINE edge_type const& edge(int i) const { return edges()[i % size()]; }
  MTINLINE point_type const& operator[](int i) const { return points()[i]; }
  MTINLINE size_t size(void) const { return 3u; }

  // Push the three triangle points (as tri_point_type objects)
  template<typename InsertIter>
  MTINLINE void insert(InsertIter it) const
  {
    *it++ = tri_point_type(data(), m_data->top_index(), TOP);
    *it++ = tri_point_type(data(), m_data->middle_index(), MIDDLE);
    *it++ = tri_point_type(data(), m_data->bottom_index(), BOTTOM);
  }

  inline iterator begin() { return points().begin(); }
  inline iterator end() { return points().end(); }
  inline const_iterator begin() const { return points().cbegin(); }
  inline const_iterator end() const { return points().cend(); }
  inline const_iterator cbegin() const { return points().cbegin(); }
  inline const_iterator cend() const { return points().cend(); }
};

} } // end namespace cfla::tri

BOOST_GEOMETRY_REGISTER_RING_TEMPLATED(cfla::tri::Triangle);
BOOST_GEOMETRY_REGISTER_RING_TEMPLATED(cfla::tri::triangle_data);

namespace cfla { namespace tri {

template<typename Tp_>
struct make_sla_traits
{
  INHERIT_TRAITS(Tp_);

  typedef typename traits::polygon_type *value_type; // isochrones
  typedef typename traits::solution_cell_type solution_type;

  typedef typename traits::event_point_type event_type;
  typedef typename traits::event_point_ycompare event_compare;
  typedef typename std::vector<event_type> event_container;
  typedef EventType       event_id_type;
  typedef make_sla_traits etraits;

  typedef sla_traits<value_type, event_type, solution_type, etraits> type;

  MTINLINE static EventType get_type(event_type const& e) {
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
  typedef std::multiset<status_seg_type, status_compare> event_status;

  typedef typename event_status::iterator status_iterator;
  typedef typename event_status::reverse_iterator status_riterator;

  // See IXEvent for a detailed description of what this is and why we need it.
  typedef std::unordered_map<
    /* Key */ point_type
    /* T */, isect_event_type
    /* Hash */, ::boost::hash<point_type>
    /* KeyEqual */, almost_equal_to<point_type>
    > isect_event_map;

  typedef typename super_type::solution_type solution_type;
  typedef typename super_type::solution_container solution_container;
  typedef typename solution_container::const_iterator solution_iterator;

  /* Keep an adjacency list of triangles which overlap so we can find the
   * solution from the max clique in the end. Each vertex of the graph maps 1-1
   * with a triangle.
   *
   * This is not as efficient as it is for rectangles (squares) but it works
   * for now. */
  typedef b::adjacency_matrix<b::undirectedS> components_graph;

private:
  triangle_container tris;
  poly_container polys;

  // sweep line status: points of current edges, ordered by X coordinate
  event_status status;
  isect_event_map intersections;
  int max_depth;
  coordinate_type _lasty;

  // connected components graph
  components_graph *graph;

  void insert_edge(edge_type const& e);
  void remove_edge(edge_type const& e);

  void handle_intersection(point_type const& event);
  void handle_tripoint(tri_point_type const& p);

  using super_type::solutions;
  using super_type::queue;

  // Look through the status in both directions to find which edges intersect
  // the given edge. Queue the intersection points.
  void check_intersections(status_iterator center);

  // Check for intersections between edge and every segment in [begin, end).
  template<typename Iter>
  void intersect_range(status_seg_type const& seg, Iter begin, Iter end);

  // Intersect the two line segments defined by Edges.
  // If there is an intersection (return code other than '0'), sets the fields
  // of ix appropriately.
  char intersect(edge_type const& a, edge_type const& b, point_type &pt) const {
    return SegSegInt(a.first(), a.second(), b.first(), b.second(), pt);
  }

  // Find a segment unambiguously in the status tree.
  // We need a special check for this since the status is a multiset.
  inline status_iterator find_unique(status_seg_type const& query) const
  {
    auto range = status.equal_range(query);
    status_iterator it = range.first;
    status_equal_compare equal;
    while (it != range.second && it != status.end() && !equal(query, *it))
      ++it;
    if (it == range.second)
      it = status.end();

    // We may or may not have found what the user wanted. We'll never know.
    return it;
  }

  MTINLINE void update_sweep(coordinate_type new_y) {
#ifdef MAXTRI_DEBUG_INTERSECT
    std::cerr << "SWEEP " << _lasty << " -> " << new_y << std::endl;
#endif
    _lasty = new_y;
  }

  template<typename Point>
  MTINLINE void update_sweep(Point new_y_point) {
    update_sweep(gety(new_y_point));
  }

#ifdef MAXTRI_DEBUG
  void dump_tris_to_octave(std::ostream& os) const;
  void dump_solution_to_octave(std::ostream& os,
      std::set<int> const& indexes, solution_cell_type const& sol) const;
#endif

#ifdef MAXTRI_DEBUG_INTERSECT
  void dump_status_to_octave(std::ostream& os) const;
  void debug_status(void) const;
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
      , graph(nullptr)
  {}

  // Construct from a list of input polygons (may be in no particular order).
  template<class Iter>
  MaxTri(Iter begin, Iter end)
    : tris(), status(status_compare(*this)), max_depth(-1), _lasty()
      , graph(nullptr)
  {
    super_type::insert(begin, end);
  }

  MTINLINE coordinate_type current_y(void) const {
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

      // new triangle_type
      tris.emplace_back(*ply, triangle, tris.size());

#ifdef MAXTRI_DEBUG_INTERSECT
      std::cerr << "QUEUE triangle" << std::endl
        << "  POINTS" << std::endl
        << "    [0] " << tris.back()[0] << std::endl
        << "    [1] " << tris.back()[1] << std::endl
        << "    [2] " << tris.back()[2] << std::endl
        << "  EDGES" << std::endl
        << "    [0] " << tris.back().edge(0) << std::endl
        << "    [1] " << tris.back().edge(1) << std::endl
        << "    [2] " << tris.back().edge(2) << std::endl;
#endif

      // queue triangle points
      tris.back().insert(emplace_inserter(queue()));
    }
  }

  MTINLINE int depth(void) const { return max_depth; }

  // Return the solution cells.
  MTINLINE size_t size(void) const { return solutions().size(); }
  MTINLINE solution_iterator begin(void) const { return solutions().cbegin(); }
  MTINLINE solution_iterator end(void) const { return solutions().cend(); }

  MTINLINE solution_type const& solution(int idx) const { return *(begin()+idx); }
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
    // To solve for player 2, build isochrones from the opposing player's
    // facilities to each user point, then triangulate the polygons (for
    // simplicity) and solve for max-depth triangular intersections.
    for (auto userp = begin(); userp != end(); ++userp)
    {
      // Only build rects for certain points.
      if (!filter(*userp))
        continue;

      // Compute an isochrone representing the fixed travel time distance to
      // the nearest facility.
      point_type nearest_facility(nn1(*userp));

      polygon_type *p = new polygon_type;
      userp->isochrone(*p, nearest_facility);
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
    auto solution = solver_.solution(chosen_one);

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

#ifdef MAXTRI_DEBUG
    // Append the final solution point to the solution script
    std::ofstream sol("./scripts/solution.m", std::ios_base::app);
    sol << "spoint = [" << randx << " , " << randy << "];" << std::endl;
    sol.close();
#endif

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

  MTINLINE static coordinate_type
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

#endif

#undef INHERIT_TRAITS

} } // end namespace cfla::tri
