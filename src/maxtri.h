#pragma once

#include <set>
#include <vector>
#include <unordered_set>
#include <functional>

#ifdef DEBUG
#include <iostream>
#include <iomanip>
#endif

#include <opencv2/core/core.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp> // point_xy, ring
#include <boost/geometry/algorithms/is_convex.hpp>

#include <boost/graph/adjacency_list.hpp>

#include "user.h"
#include "polygon.h"
#include "sweep.h"
#include "nn1.h"

namespace cfla
{

namespace tri
{

// namespaces

namespace b  = boost;
namespace bp = boost::polygon;
namespace bg = boost::geometry;
namespace bgm = boost::geometry::model;

template<typename Tp_>
  using point_xy = boost::geometry::model::d2::point_xy<Tp_>;

// forward declarations

template<typename Tp_>
struct Triangle;

template<typename Tp_>
struct Edge;

template<typename Tp_>
struct EdgePoint;

template<typename Tp_>
struct SolutionCell;

template<typename Tp_>
struct traits
{
  typedef Tp_                            coordinate_type;
  typedef point_xy<Tp_> cv::Point_<Tp_>  point_type;
  typedef Edge<Tp_>                      edge_type;
  typedef EdgePoint<Tp_>                 edge_point_type;
  typedef Triangle<Tp_>                  triangle_type;
  typedef c_polygon<Tp_>                 polygon_type;
  typedef SolutionCell<Tp_>              solution_ref_type;
  typedef bgm::ring<point_xy<Tp_> >      solution_cell_type;
};

#define INHERIT_TRAITS(Targ) \
  typedef cfla::tri::traits<Targ> traits; \
  typedef typename traits::coordinate_type coordinate_type; \
  typedef typename traits::point_type point_type; \
  typedef typename traits::edge_type edge_type; \
  typedef typename traits::edge_point_type edge_point_type; \
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
    inline bool is_convex(typename traits<Tp_>::solution_cell_type const& c) {
      return true;
    }
} } // end namespace boost::geometry

/* Intersect many geometries into out.  */
template<typename Iter, typename Box, typename Polygon>
bool
intersection(Iter begin, Iter end, Box const& bbox, Polygon& out)
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
    // XXX is providing out twice okay??
    if (!bg::intersection(*begin++, out, out))
      return false;

  return true;
}

// classes

template<class Tp_>
  struct EdgePoint : public point_xy<Tp_>
{
  // typedefs
  typedef point_xy<Tp_> super;
  INHERIT_TRAITS(Tp_);

  // members
  edge_type &const edge; // parent edge
  const bp::direction_1d dir; // LOW or HIGH, meaning first or second

  // constructors
  EdgePoint() : super() {}

  EdgePoint(edge_type& parent, bp::direction_1d dir,
      coordinate_type x, coordinate_type y)
    : super(x, y), edge(parent), dir(d)
  {}

  EdgePoint(EdgePoint const& o, edge_type const& new_parent, bp::direction_1d d)
    : super(o), edge(new_parent), dir(d)
  {}

  // methods

  // Edge points are lex. sorted by x then y
  inline bool operator<(EdgePoint const& o) {
    return (x() < o.x()) || (x() == o.x() && y() < o.y());
  }

  // return the point at the other end of this edge
  inline other(void) const {
    return (dir == bp::LOW) ? edge.second : edge.first;
  }
};

template<class Tp_>
struct Edge : public bgm::segment<EdgePoint<Tp_> >
{
  // typedefs
  typedef bgm::segment<EdgePoint<Tp_> > super;
  INHERIT_TRAITS(Tp_);

  // members
  triangle_type &const tri; // parent triangle
  int index = -1; // index of edge in parent triangle (0, 1, or 2)
  typename bp::direction_1d dir; // LOW or HIGH, meaning left or right edge
  int depth; // intersection depth

  // constructors
  Edge(triangle_type& parent, int tidx,
      point_type p1, point_type p2, bp::direction_1d d)
    : super(p1, p2), tri(parent), index(tidx), dir(d), depth(-1)
  {}

  template<class Segment>
  Edge(triangle_type& parent, int tidx, Segment s, bp::direction_1d d)
    : super(s), tri(parent), index(tidx), dir(d), depth(-1)
  {}

  Edge(triangle_type& parent, int tidx, Edge const& other)
    : super(other), tri(parent), index(tidx), dir(other.dir), depth(other.depth)
  {}

  // methods

  inline bool operator==(Edge const& e) const {
    return (&tri == &e.tri) && (dir == e.dir)
      && almost_equal(first, e.first) && almost_equal(second, e.second)
  }

  // return the next/previous linked edge
  Edge const& nextEdge(void) const { return tri.edge(index + 1); }
  Edge const& prevEdge(void) const { return tri.edge(index - 1); }
};

template<class Tp_>
struct Triangle
{
  // typedefs
  typedef Tp_ coordinate_type;
  INHERIT_TRAITS(Tp_);

  typedef std::array<edge_type, 3> edge_array;

  // members
  polygon_type &const poly; // parent polygon
  edge_array edges;

  // constructors
  template <class Tri>
  Triangle(polygon_type const& parent, Tri const& t)
    : poly(p), edges({ {*this, t[0]}, {*this, t[1]}, {*this, t[2]} })
  {}

  // methods
  inline edge_type const& operator[](int i) const { return edges[i % size()]; }
  inline size_t size(void) const { return edges.size(); }

  inline edge_type const& edge(void) const { return edges[i % size()]; }
};


template<class Tp_>
struct SolutionCell
{
public:
  INHERIT_TRAITS(Tp_);

  // triangles which form the maximal intersection
  typedef typename std::unordered_set<
    std::reference_wrapper<const triangle_type> > triangle_ref_set;

  typedef typename triangle_ref_set::const_iterator iterator;
  typedef typename triangle_ref_set::size_type size_type;

private:
  triangle_ref_set source_tris_;
  int depth_;

public:
  // A solution cell must be formed from at least two triangles.
  SolutionCell(triangle_type const& t1, triangle_type const& t2)
    : source_tris_({ {t1}, {t2} }), depth_(-1)
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
  typedef cfla::tri::traits<Tp_> traits;
  typedef typename traits::polygon_type value_type; // polygons are the input
  typedef typename traits::solution_ref_type solution_type;

  typedef typename traits::edge_point_type    event_type;
  typedef typename traits::edge_point_compare event_compare;
  typedef typename traits:edge_type     event_id_type;
  typedef make_sla_traits etraits;

  typedef sla_traits<value_type, event_type, solution_type, etraits> type;

  inline static event_id_type get_type(event_type const& p) { return p.edge; }
};

template<class Tp_>
class MaxTri
  : public SweepLineAlgorithm<typename make_sla_traits<Tp_>::type>
{
public:
  typedef SweepLineAlgorithm<typename make_sla_traits<Tp_>::type> super_type;

  INHERIT_TRAITS(Tp_);

  typedef typename std::vector<triangle_type> triangle_container;
  typedef typename std::set<edge_type, edge_compare> edge_point_set;

  // Vertexes of the adjacency list are triangles.
  typedef typename b::property<b::vertex_index_t, int> VertexProps;
  typedef typename b::adjacency_list<
      b::hash_setS, b::listS, b::undirectedS, VertexProps
    > components_graph;
  typedef typename b::graph_traits<components_graph>::vertex_descriptor
    vertex_descriptor;
  typedef typename b::property_map<components_graph, b::vertex_index_t>::type
    id_map;
  typedef typename std::vector<vertex_descriptor> descriptor_list;

  typedef typename edge_point_set::iterator edge_point_iterator;

  typedef typename super_type::solution_container solution_container;
  typedef typename super_type::solution_iterator solution_iterator;

private:
  triangle_container tris;

  // sweep line status: points of current edges, ordered by X coordinate
  edge_point_set edge_points;

  // data needed for the adjacency list
  components_graph graph;
  id_map component_ids;
  descriptor_list vertexes;
  int max_depth;

  inline vertex_descriptor vd(int idx) const { return vertexes[idx]; }

  void insert_edge(edge_type const& e);
  void handle_intersection(edge_type const& e1, edge_type const& e2);
  void remove_edge(edge_type const& e);
  using super_type::solutions;
  using super_type::queue;

protected:
  // override
  void handle_event(edge_type const& edge, edge_point_type const& event);
  void initialize(void);
  void finalize(void);

public:
  ~MaxTri();

  // No inputs yet.
  MaxTri();

  // Construct from a list of input polygons (may be in no particular order).
  template<class Iter>
    MaxTri(Iter begin, Iter end)
    : tris(), edge_points(), graph(),
      component_ids(b::get(b::vertex_index_t(), graph)),
      vertexes(), max_depth(-1)
    {
      super_type::insert(begin, end);
    }

  inline int index(vertex_descriptor v) const {
    return b::get(component_ids, v);
  }

  components_graph const& adj_graph(void) const { return graph; }

  // Add triangles from polygon.
  // Same as the constructor form, if you're lazy and want to do it
  // after construction. (Overridden from SLA.)
  void add_event(polygon_type const& ply)
  {
    if (ply.triangles_size() == 0u)
      throw std::runtime_error("must triangulate the polygon first!");

    auto it = ply.triangles_begin();
    while (it != ply.triangles_end())
    {
      tris.emplace_back(ply, *it);
      ++it;
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

  typedef User<Tp_> user_type;
  typedef std::vector<user_type> user_list;
  typedef typename user_list::iterator       iterator;
  typedef typename user_list::const_iterator const_iterator;

  typedef cfla::L1NN1<coordinate_type, point_type> nn1_type;
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
      point_type nearest_facility(nn1(userp->center()));
      solver_.add_event(userp->isochrome(nearest_facility));
    }
    solver_.compute();

    // We have potentially multiple solution cells.
    if (solver_.size() == 0)
    {
      std::cerr << "warning: empty solution" << std::endl;
      return point_type(0, 0);
    }

    // Choose a random point in a random solution cell.
    size_type chosen_one = randrange(size_type(0), solver_.size()-1);
    auto solution = solver_.cell(chosen_one);

    coordinate_type randx = randrange(
        bp::get(solution, bp::HORIZONTAL, bp::LOW),
        bp::get(solution, bp::HORIZONTAL, bp::HIGH));
    coordinate_type randy = randrange(
        bp::get(solution, bp::VERTICAL, bp::LOW),
        bp::get(solution, bp::VERTICAL, bp::HIGH));
    point_type out(randx, randy);

    return out;
  }


public:

  template<typename UserIter>
    MaxTriSolver(UserIter begin, UserIter end)
      : solver_(), users_(begin, end)
    {}

  iterator begin(void) const { return users_.begin(); }
  iterator end(void) const { return users_.end(); }
  //const_iterator cbegin(void) const { return users_.cbegin(); }
  //const_iterator cend(void) const { return users_.cend(); }

  // Return the "optimal solution" for the CFL game.
  template<typename point_filter>
  point_type operator()(const nn1_type &facilities, point_filter filter) {
    return compute(facilities, filter);
  }

  inline static coordinate_type
    distance(user_type const& user, point_type const& p)
      return user.travelTime(p);
    }

};

extern template class MaxTri<double>;
extern template class MaxTri<float>;

#ifdef DEBUG
template<typename U>
std::ostream& operator<<(std::ostream& os, cfla::Edge<U> const& e) {
  os << "<[" << std::setw(2) << std::setfill(' ') << e.rect_index
    << "] " << (e.dir == boost::polygon::LOW ? "LOW " : "HIGH")
    << " " << e.coord << " d=" << e.depth << ">";
  return os;
}

template<typename U>
std::ostream& operator<<(std::ostream& os, cfla::SolutionEdge<U> const& e)
{
  os << "[" << std::setw(2) << std::setfill(' ') << e.solution << "] from "
    << e.edge;
  return os;
}
#endif

#undef INHERIT_TRAITS

} // end namespace cfla

#endif
