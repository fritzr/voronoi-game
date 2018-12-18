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

#include "user.h"
#include "polygon.h"
#include "intersection.h"
#include "sweep.h"
#include "nn1.h"

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
struct event_point_compare_x;

template<typename Tp_>
struct event_point_compare_y;

template<typename Tp_>
class IXPoint;

template<typename Tp_>
class EventPoint;

template<typename Tp_>
struct traits
{
  typedef Tp_                            coordinate_type;
  typedef point_xy<Tp_>                  point_type;
  typedef Edge<Tp_>                      edge_type;
  typedef EdgePoint<Tp_>                 edge_point_type;
  typedef event_point_compare_x<Tp_>     event_point_xcompare;
  typedef event_point_compare_y<Tp_>     event_point_ycompare;
  typedef IXPoint<Tp_>                   isect_point_type;
  typedef EventPoint<Tp_>                event_point_type;
  typedef Triangle<Tp_>                  triangle_type;
  // XXX could be point_type
  typedef c_polygon<point_type>          polygon_type;
  typedef SolutionCell<Tp_>              solution_ref_type;
  typedef bgm::ring<point_type>          solution_cell_type;

};

#define INHERIT_TRAITS(Targ) \
  typedef cfla::tri::traits<Targ> traits; \
  typedef typename traits::coordinate_type coordinate_type; \
  typedef typename traits::point_type point_type; \
  typedef typename traits::edge_type edge_type; \
  typedef typename traits::edge_point_type edge_point_type; \
  typedef typename traits::event_point_xcompare event_point_xcompare; \
  typedef typename traits::event_point_ycompare event_point_ycompare; \
  typedef typename traits::isect_point_type isect_point_type; \
  typedef typename traits::event_point_type event_point_type; \
  typedef typename traits::triangle_type triangle_type; \
  typedef typename traits::polygon_type polygon_type; \
  typedef typename traits::solution_ref_type solution_ref_type; \
  typedef typename traits::solution_cell_type solution_cell_type; \

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
          { this, 0, input_tri[0], input_tri[0], bp::HIGH },
          { this, 1, input_tri[0], input_tri[0], bp::LOW },
          { this, 2, input_tri[0], input_tri[0], bp::LOW },
        })
  {
  }

  // parent edge
  inline edge_type const& edge(int index) const {
    return edges[index % 3];
  }
};

template<class Tp_>
struct edge_data
{
  INHERIT_TRAITS(Tp_);

  // members
  edge_point_type first, second;
  triangle_data<Tp_> *tri; // parent triangle -- weak reference
  int index; // index of edge in parent triangle (0, 1, or 2)
  bp::direction_1d dir; // LOW (left) or HIGH (right)
  int depth; // intersection depth

  edge_data(triangle_data<Tp_> *tridata, int edge_index,
      point_type start, point_type end, bp::direction_1d d)
    : first(this, bp::LOW, start), second(this, bp::HIGH, end), tri(tridata)
      , index(edge_index), dir(d), depth(-1)
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

  inline const triangle_data<Tp_> *tridata(void) const {
    return edata->tri;
  }

#ifdef MAXTRI_DEBUG
  inline friend std::ostream& operator<<(std::ostream& os, EdgePoint const& p) {
    os << "(" << p.x() << " , " << p.y() << ")";
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

  template<typename Pt_>
  IXPoint(Pt_ const& point,
      const edge_data<Tp_> *edge1,
      const edge_data<Tp_> *edge2)
    : super(point), e1(edge1), e2(edge2)
  {}
};

template<typename Tp_>
struct event_data
{
  INHERIT_TRAITS(Tp_);

  union {
    edge_point_type edge;
    isect_point_type isect;
  } u;

  bool intersection;

  event_data(edge_point_type const& e)
    : u{.edge=e}, intersection(false) {}
  //{ new (&u.edge) edge_point_type(e); }

  event_data(isect_point_type const& i)
    : u{.isect=i}, intersection(true) {}
  // { new (&u.isect) isect_point_type(i); }

  ~event_data()
  {
    if (intersection)
      (u.edge).~EdgePoint();
    else
      (u.isect).~IXPoint();
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

  EventPoint(edge_point_type const& point)
    : m_union(new m_type(point)) {}

  EventPoint(IXPoint<Tp_> const& ipoint)
    : m_union(new m_type(ipoint)) {}

  inline bool intersection(void) const {
    return m_union->intersection;
  }

  inline edge_point_type const& edge(void) const {
    if (intersection())
      throw std::runtime_error("EventPoint: requested intersection; is edge");
    return m_union->u.edge;
  }

  inline isect_point_type const& isect(void) const {
    if (!intersection())
      throw std::runtime_error("EventPoint: requested edge; is intersection");
    return m_union->u.isect;
  }

  inline coordinate_type x(void) const {
    return intersection() ? isect().x() : edge().x();
  }

  inline coordinate_type y(void) const {
    return intersection() ? isect().y() : edge().y();
  }

};

template<typename Tp_>
struct event_point_compare_y
{
  INHERIT_TRAITS(Tp_);

  // When comparing two edge points A and C in edges E1 = <A,B> and E2 = <C,D>
  // the points are lexicographically sorted by following criteria:
  //   #.  A.y < C.y
  //   #.  A is right edge < C is left edge
  //   #.  A.x < C.x
  //   #.  B.y < D.y
  //   #.  B is right edge < D is left edge
  //   #.  B.x < D.x
  inline bool operator()(event_point_type const& p1, event_point_type const& p2,
      bool recursing=false) const
  {
    bool y1 = p1.y() < p2.y();
    bool y2 = p1.y() == p2.y();
    bool d1 = !p1.intersection() && !p2.intersection()
      && p1.edge().dir == bp::HIGH && p2.edge().dir == bp::LOW;
    bool d2 = (p1.intersection() || p2.intersection())
      || (!p1.intersection() && !p2.intersection()
          && (p1.edge().dir == p2.edge().dir || p1.edge().dir == bp::LOW));
    bool x = p1.x() < p2.x();

    return y1 || (y2 && (d1 || (d2 && x)));
      // || (!recursing && (*this)(p1.edge().other(), p2.edge().other(), true))
  }

};

template<typename Tp_>
struct event_point_compare_x
{
  INHERIT_TRAITS(Tp_);

  inline bool operator()(event_point_type const& p1, event_point_type const& p2)
    const
  {
    bool x1 = p1.x() < p2.x();
    bool x2 = p1.x() == p2.x();
    bool d1 = !p1.intersection() && !p2.intersection()
      && p1.edge().dir == bp::HIGH && p2.edge().dir == bp::LOW;
    bool d2 = (p1.intersection() || p2.intersection())
      || (!p1.intersection() && !p2.intersection()
          && (p1.edge().dir == p2.edge().dir || p1.edge().dir == bp::LOW));
    bool y = p1.y() < p2.y();
    return x1 || (x2 && (d1 || (d2 && y)));
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
  inline const triangle_data<Tp_> *tridata(void) const { return m_data->tri; }

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
    os << "[ " << e.first() << " => " << e.second() << " ]";
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
 *        *
 * V E0  /  \  E1 V
 *      /     \
 * [1] *- _     \
 *         - _   \
 *             - _* [2]
 *        E2 >
 *
 * E0 === 0 -> 1 ; dir = LOW (left)
 * E1 === 0 -> 2 ; dir = HIGH (right)
 * E2 === 1 -> 2 ; dir = LOW (left)
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
    int maxy_index = 0;
    point_type maxy_point = input_tri[0];

    // Find the highest point: this is index 0.
    if (input_tri[1].y() > maxy_point.y()) maxy_index = 1;
    if (input_tri[2].y() > maxy_point.y()) maxy_index = 2;
    maxy_point = input_tri[maxy_index];

    // Find the order of points that form a left turn.
    int left_index = (maxy_index+2) % 3; // equivalent to -1, but not negative
    int right_index = (maxy_index+1) % 3;
    if (leftTurn(input_tri[left_index], maxy_point, input_tri[right_index]))
    {
      int swp = left_index;
      left_index = right_index;
      right_index = swp;
    }

    // Now we have the point order.
    points()[0] = input_tri[maxy_index];
    points()[1] = input_tri[left_index];
    points()[2] = input_tri[right_index];

    // Find the lowest of (left, right): this is the second index of E2.
    int e2_first = 1, e2_second = 2;
    bp::direction_1d e2_dir = bp::LOW;
    if (points()[1].y() < points()[2].y())
    {
      e2_first = 2;
      e2_second = 1;
      e2_dir = bp::HIGH;
    }

    // Finally we can finish the edges.
    edges()[0].first().set(points()[0]);
    edges()[0].second().set(points()[1]);
    edges()[1].first().set(points()[0]);
    edges()[1].second().set(points()[2]);
    edges()[2].first().set(points()[e2_first]);
    edges()[2].second().set(points()[e2_second]);
    edges()[2].dir(e2_dir);
  }

private:
  inline point_array& points(void) { return m_data->points; }
  inline edge_array& edges(void) { return m_data->edges; }

public:
  inline point_array const& points(void) const { return m_data->points; }
  inline edge_array const& edges(void) const { return m_data->edges; }
  inline polygon_type const& poly(void) const { return m_data->poly; }

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

  typedef typename traits::polygon_type value_type; // polygons are the input
  typedef typename traits::solution_ref_type solution_type;

  typedef typename traits::event_point_type event_type;
  typedef typename traits::event_point_ycompare event_compare;
  typedef typename std::vector<event_type> event_container;
  typedef bool       event_id_type;
  typedef make_sla_traits etraits;

  typedef sla_traits<value_type, event_type, solution_type, etraits> type;

  inline static bool get_type(event_type const& e) {
    return e.intersection();
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

  typedef std::vector<triangle_type> triangle_container;
  typedef std::set<edge_point_type, event_point_xcompare> edge_point_status;

  typedef typename edge_point_status::iterator edge_point_iterator;

  typedef typename super_type::solution_type solution_type;
  typedef typename super_type::solution_container solution_container;
  typedef typename solution_container::const_iterator solution_iterator;

private:
  triangle_container tris;

  // sweep line status: points of current edges, ordered by X coordinate
  edge_point_status edge_points;
  int max_depth;

  void insert_edge(edge_type const& e);
  void handle_intersection(isect_point_type const& i);
  void remove_edge(edge_type const& e);
  using super_type::solutions;
  using super_type::queue;

protected:
  // override
  void handle_event(bool intersection_event, event_point_type const& event);
  void initialize(void);
  void finalize(void);

public:
  ~MaxTri();

  // No inputs yet.
  MaxTri();

  // Construct from a list of input polygons (may be in no particular order).
  template<class Iter>
    MaxTri(Iter begin, Iter end)
    : tris(), edge_points(), max_depth(-1)
    {
      super_type::insert(begin, end);
    }

  // Add triangles from polygon.
  // Same as the constructor form, if you're lazy and want to do it
  // after construction. (Overridden from SLA.)
  void add_event(polygon_type const& ply)
  {
    if (ply.triangles_size() == 0u)
      throw std::runtime_error("must triangulate the polygon first!");

    for (auto ptri = ply.triangles_begin(); ptri != ply.triangles_end(); ++ptri)
    {
      point_type triangle[3] = {
        ply[ptri->v[0]]->getPos(),
        ply[ptri->v[1]]->getPos(),
        ply[ptri->v[2]]->getPos(),
      };
      tris.emplace_back(ply, triangle);

#ifdef MAXTRI_DEBUG
      std::cout << "QUEUE edges" << std::endl
        << "    [0] " << tris.back().edge(0) << std::endl
        << "    [1] " << tris.back().edge(1) << std::endl
        << "    [2] " << tris.back().edge(2) << std::endl;
#endif

      queue().emplace(tris.back().edge(2).second());
      queue().emplace(tris.back().edge(2).first());
      queue().emplace(tris.back().edge(1).second());
      queue().emplace(tris.back().edge(1).first());
      queue().emplace(tris.back().edge(0).second());
      queue().emplace(tris.back().edge(0).first());
    }
  }

  // Safe, non-const version.
  void add_event(polygon_type& ply)
  {
    ply.triangulate();
    add_event(const_cast<polygon_type const&>(ply));
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
  typedef MaxTri<coordinate_type> solver_type;

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

      // Compute an iso representing the travel time
      point_type nearest_facility(nn1(*userp));
      solver_.add_event(userp->isochrome(nearest_facility));
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
    << "] " << (e.dir == bp::LOW ? "LOW " : "HIGH")
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
