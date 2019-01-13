#pragma once

// Tree structures with custom node data.
#include <set>
#include <vector>
#include <unordered_set>

#ifdef DEBUG
#include <iostream>
#include <iomanip>
#endif

// XXX not all of the classes here have been instantiated recently.
// some things may have been broken since refactoring for the maxtri code.

#include <boost/polygon/polygon.hpp>
#include <boost/polygon/point_concept.hpp>
#include <boost/polygon/rectangle_concept.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>

#include <boost/graph/adjacency_list.hpp>

#include "boost_geo_poly.h"
#include "sweep.h"
#include "nn1.h"


namespace cfla { namespace rect
{

namespace b  = boost;
namespace bp = boost::polygon;
namespace bg = boost::geometry;

template<typename Tp_>
  using point_xy = boost::geometry::model::d2::point_xy<Tp_>;

template<class Tp_>
struct Edge
{
  typedef typename bp::rectangle_traits<Tp_> traits;
  typedef Tp_ coordinate_type;

  coordinate_type coord; // actual coordinate
  typename bp::direction_1d dir; // LOW or HIGH
  int rect_index; // parent rectangle
  mutable int depth; // does not affect sorting

  Edge() : coord(), dir(), rect_index(-1), depth(-1) {}
  Edge(coordinate_type const& c, bp::direction_1d const& d, int idx)
    : coord(c), dir(d), rect_index(idx), depth(-1) {}
  Edge(Edge const& other, int new_idx)
    : coord(other.coord), dir(other.dir), rect_index(new_idx),
      depth(other.depth) {}

  inline bool operator<(Edge const& e) const
  {
    return coord < e.coord
      || ((coord == e.coord) && (rect_index < e.rect_index));
  }

  inline bool operator==(Edge const& e) const {
    return rect_index == e.rect_index
      && dir == e.dir
      && coord == e.coord; // epsilon compare?
    // depth compare?
  }

};

template<class Tp_>
struct RectWrapper : bp::rectangle_data<Tp_>
{
  typedef typename bp::rectangle_data<Tp_> super_type;
  typedef typename bp::rectangle_traits<Tp_> traits;
  typedef Edge<Tp_> edge_type;

  size_t index;
  //int component; // ID of the component to which this rect belongs.

  template <class HInterval, class VInterval>
    RectWrapper(const HInterval& hrange, const VInterval& vrange, size_t idx)
      : super_type(hrange, vrange), index(idx)/*, component(-1)*/ {}

  inline edge_type edge(bp::orientation_2d orient, bp::direction_1d dir) const
  {
    return edge_type(bp::get(*this, orient, dir), dir, index);
  };

  template<class OutIter> void
    add_edges(OutIter out, bp::orientation_2d orient) {
      *out++ = edge(orient, bp::HIGH);
      *out++ = edge(orient, bp::LOW);
    }
};

} } // we interrupt your regularly scheduled cfla namespace to bring you...

namespace boost { namespace polygon {
  template <typename Tp_>
    struct geometry_concept<cfla::rect::RectWrapper<Tp_> >
      { typedef rectangle_concept type; };
}}


namespace cfla { namespace rect
{

// Node to track 'cell depths'. These form the tree T(V_k) from the paper,
// the leaves of which represent g_k[j]: the depth of the cell at position
// j in the current sweep partition.
/*
template<class Tp_>
class DepthNode : public tree::avl_node<DepthNode<Tp_> >
{
public:
  typedef Tp_ coordinate_type;
  typedef DepthNode<Tp_> my_type;
  typedef tree::avl_node<my_type> base_type;
  DECLARE_TRAITS(tree::avl_node_traits<my_type>);
  typedef Edge<coordinate_type> edge_type;
  typedef depth_compare<my_type> compare;

  // Internal depth value of the node.
  int depth;
  // If we are a leaf node, which X edge we correspond to.
  edge_type leaf_edge;

  DepthNode() : base_type(), depth(0), leaf_edge() {}

  // Whether we are a leaf node.
  inline bool leaf(void) const {
    return !node_traits::get_left(this) && !node_traits::get_right(this);
  }
};
*/

template<class Tp_>
struct SolutionCompare;

template<class Tp_, class EdgeType=Edge<Tp_> >
struct SolutionEdge
{
  typedef EdgeType edge_type;
  typedef typename edge_type::coordinate_type coordinate_type;
  typedef SolutionEdge<Tp_> sedge_type;

  edge_type edge;
  int solution;

  SolutionEdge(edge_type const& e, int sidx=-1)
    : edge(e), solution(sidx)
  {
  }

  inline bool operator<(sedge_type const& other) const {
    return edge < other.edge
      || (!(other.edge < edge) && solution < other.solution);
  }
  inline bool operator<(edge_type const& other) const {
    return edge < other;
  }
  inline bool operator==(sedge_type const& other) const {
    return !(*this < other) && !(other < *this);
  }
  inline bool operator==(edge_type const& other) const {
    return !(*this < other) && !(other < this->edge);
  }

};

template<class Tp_>
struct SolutionCell;

template<typename Tp_>
struct MaxDepthRectTraits
{
  typedef Tp_                            coordinate_type;
  typedef point_xy<Tp_>       point_type;
  typedef RectWrapper<Tp_>             rect_type;
  typedef typename rect_type::super_type pure_rect_type;
  typedef Edge<Tp_>                      edge_type;
  typedef std::less<edge_type>           edge_compare;

  typedef SolutionCell<Tp_>                  solution_type;
  typedef typename solution_type::sedge_type sedge_type;
  typedef std::less<sedge_type>              sedge_compare;
};

template<class Tp_>
struct SolutionCell
{
public:
  typedef MaxDepthRectTraits<Tp_> traits;
  typedef typename traits::coordinate_type coordinate_type;
  typedef typename traits::point_type      point_type;
  typedef typename traits::rect_type       rect_type;
  typedef typename traits::edge_type       edge_type;
  typedef typename traits::edge_compare    edge_compare;
  typedef typename traits::pure_rect_type  pure_rect_type;

  typedef SolutionEdge<Tp_> sedge_type;
  typedef std::less<sedge_type> sedge_compare;

  // indexes of rects which form the maximal intersection
  typedef typename std::unordered_set<size_t> index_set;
  typedef typename index_set::const_iterator index_iterator;
  typedef typename index_set::size_type size_type;

private:
  index_set source_rects_;
  int top_, bot_, left_, right_;
  bool hit_left_, hit_right_;
  int depth_;

public:
  SolutionCell(int top, int left, int right, int depth)
    : source_rects_(),
      top_(top), bot_(-1), left_(left), right_(right),
      hit_left_(false), hit_right_(false), depth_(depth)
  {
    source_rects_.insert(top_);
    source_rects_.insert(left_);
    source_rects_.insert(right_);
#ifdef MAXRECT_DEBUG
    std::cerr << "new max from rect " << top_
      << " and edges: " << left_ << " and " << right_ << std::endl;
#endif
  }

  inline int depth(void) const { return depth_; }

  inline void found(bp::direction_1d dir) {
    if (dir == bp::LEFT)  hit_left_  = true;
    //else
    if (dir == bp::RIGHT) hit_right_ = true;
  }
  inline bool marked(void) const { return hit_left_ || hit_right_; }
  inline bool marked(int bot) {
    if (marked())
    {
      set_bot(bot);
      return true;
    }
    return false;
  }

  inline int top(void) const { return top_; }
  inline int bot(void) const { return bot_; }
  inline int left(void) const { return left_; }
  inline int right(void) const { return right_; }
  inline void set_bot(int new_bot) {
    // only update bottom once
    if (bot_ < 0 && new_bot >= 0)
      source_rects_.insert(bot_ = new_bot);
  }

  inline index_iterator begin(void) const { return source_rects_.cbegin(); }
  inline index_iterator end(void) const { return source_rects_.cend(); }
  inline size_type size(void) const { return source_rects_.size(); }

  // solution cell
  template<class RectArray>
  pure_rect_type cell(RectArray const& rects, size_t rects_size) const {
    auto rit = source_rects_.begin();
    if (rit == source_rects_.end())
      return pure_rect_type();
    int idx = *rit++;
    if (idx < 0 || static_cast<size_t>(idx) >= rects_size)
      return pure_rect_type();

    pure_rect_type retcell = rects[idx];
    while (rit != source_rects_.end())
    {
      idx = *rit++;
      if (idx >= 0 && static_cast<size_t>(idx) < rects_size)
        /*assert(*/bp::intersect(retcell, rects[idx])/*)*/;
    }
    return retcell;
  }
};

template<typename Tp_>
struct make_sla_traits
{
  typedef MaxDepthRectTraits<Tp_> traits;
  typedef typename traits::pure_rect_type value_type;
  typedef typename traits::solution_type solution_type;

  typedef typename traits::edge_type    event_type;
  typedef typename traits::edge_compare event_compare;
  typedef typename bp::direction_1d     event_id_type;
  typedef make_sla_traits etraits;

  typedef sla_traits<value_type, event_type, solution_type, etraits> type;

  inline static event_id_type get_type(event_type const& e) { return e.dir; }
};

template<class Tp_>
class MaxRect
  : public SweepLineAlgorithm<typename make_sla_traits<Tp_>::type>
{
public:
  typedef SweepLineAlgorithm<typename make_sla_traits<Tp_>::type> super_type;

  typedef MaxDepthRectTraits<Tp_> traits;
  typedef typename traits::coordinate_type coordinate_type;
  typedef typename traits::point_type      point_type;
  typedef typename traits::rect_type       rect_type;
  typedef typename traits::edge_type       edge_type;
  typedef typename traits::edge_compare    edge_compare;
  typedef typename rect_type::super_type   pure_rect_type;

  typedef typename traits::solution_type   solution_type;
  typedef typename traits::sedge_type      sedge_type;
  typedef typename traits::sedge_compare   sedge_compare;

  typedef typename std::vector<edge_type> edge_container;
  typedef typename std::vector<rect_type> rect_container;
  typedef typename std::set<edge_type, edge_compare> edge_set;
  // T(V_k) from the paper - TODO
  //typedef typename tree::make_avltree<DepthNode<Tp_> >::type depth_tree_type;
  // for now we take a short-cut and keep the entire array gk[j]

  // Vertexes are rect_type objects.
  typedef typename b::property<b::vertex_index_t, int> VertexProps;
  typedef typename b::adjacency_list<
      b::hash_setS, b::listS, b::undirectedS, VertexProps
    > components_graph;
  typedef typename b::graph_traits<components_graph>::vertex_descriptor
    vertex_descriptor;
  typedef typename b::property_map<components_graph, b::vertex_index_t>::type
    id_map;
  typedef typename std::vector<vertex_descriptor> descriptor_list;

  typedef typename edge_set::iterator edge_iterator;

  typedef typename std::set<sedge_type, sedge_compare> solution_edge_set;
  typedef typename super_type::solution_container solution_container;
  typedef typename super_type::solution_container::const_iterator
    solution_iterator;

  typedef typename std::vector<pure_rect_type> pure_rect_list;
  typedef typename pure_rect_list::const_iterator pure_rect_iterator;

private:
  rect_container rects;
  edge_set edges_x;   // vertical rect edges within each line event (status)
  components_graph graph;
  id_map component_ids;
  descriptor_list vertexes;
  int max_depth;

  solution_edge_set solution_edges;
  pure_rect_list solution_cells;

  inline vertex_descriptor vd(int idx) const { return vertexes[idx]; }

  inline rect_type &rect(edge_type const& edge) {
    return rects[edge.rect_index];
  }

  template<typename EdgeSetIter>
    int check_max_depth(rect_type const& r, EdgeSetIter edge_lb, int depth);

  void insert_rect(rect_type const& r);
  void remove_rect(rect_type const& r);
  using super_type::solutions;
  using super_type::queue;

protected:
  void handle_event(bp::direction_1d dir, edge_type const& event);
  void initialize(void);
  void finalize(void);

public:
  ~MaxRect();

  // No inputs yet.
  MaxRect();

  // Construct from a list of input rectangles (may be in no particular order).
  template<class RectIter>
    MaxRect(RectIter begin, RectIter end)
    : rects(), edges_x(), graph(),
      component_ids(b::get(b::vertex_index_t(), graph)),
      vertexes(), max_depth(-1)
    {
      super_type::insert(begin, end);
    }

  inline int index(vertex_descriptor v) const {
    return b::get(component_ids, v);
  }

  components_graph const& adj_graph(void) const { return graph; }

  void add_rect(rect_type const& rect);

  // Add rectangles. Same as the constructor form, if you're lazy and want to
  // do it sometime after construction.
  template<class RectType>
    void add_event(RectType const& rect)
    {
      // We must make sure the horizontal/vertical intervals are in the proper
      // order, otherwise our algorithm will explode since we assume LOW < HIGH.
      auto hivl = bp::get(rect, bp::HORIZONTAL);
      auto vivl = bp::get(rect, bp::VERTICAL);
      auto hl = bp::get(hivl, bp::LOW), hh = bp::get(hivl, bp::HIGH);
      auto vl = bp::get(vivl, bp::LOW), vh = bp::get(vivl, bp::HIGH);
      hivl = bp::construct<decltype(hivl)>(std::min(hl,hh), std::max(hl,hh));
      vivl = bp::construct<decltype(vivl)>(std::min(vl,vh), std::max(vl,vh));
      add_rect(rect_type(hivl, vivl, rects.size()));
#ifdef MAXRECT_DEBUG
      std::cerr << "ADD [" << (rects.size()-1) << "]"
        << " h(" << bp::get(hivl, bp::LOW) << "," << bp::get(hivl, bp::HIGH)
        << ")"
        << " v(" << bp::get(vivl, bp::LOW) << "," << bp::get(vivl, bp::HIGH)
        << ")" << std::endl;
#endif
    }

  // override
  void add_event(pure_rect_type const& rect) {
    add_event<pure_rect_type>(rect);
  }

  inline int depth(void) const { return max_depth; }

  // Return the solution cells.
  inline size_t size(void) const { return solution_cells.size(); }

  inline pure_rect_iterator begin(void) const { return solution_cells.cbegin(); }
  inline pure_rect_iterator end(void) const { return solution_cells.cend(); }
  inline pure_rect_type const& cell(size_t i) { return solution_cells[i]; }

  inline solution_iterator sol_begin(void) const { return solutions().cbegin(); }
  inline solution_iterator sol_end(void) const { return solutions().cend(); }
  inline solution_type const& solution(size_t i) const { return solutions()[i]; }
};

/* Adapter for MaxRect which conforms to the cfla_traits::solver_type
 * interface. That is, defines Point operator()(facilities...) and
 * iterator methods that iterate over input user points.
 */
template<class Tp_,
  typename PointContainer=std::list<typename MaxDepthRectTraits<Tp_>::Point>
>
class MaxRectSolver
{
public:
  typedef MaxDepthRectTraits<Tp_> traits;
  typedef typename traits::coordinate_type coordinate_type;
  typedef typename traits::point_type      point_type;
  typedef typename traits::rect_type       rect_type;
  typedef typename traits::edge_type       edge_type;
  typedef typename traits::edge_compare    edge_compare;
  typedef typename rect_type::super_type   pure_rect_type;

  typedef PointContainer point_list;
  typedef typename point_list::iterator       iterator;
  typedef typename point_list::const_iterator const_iterator;

  typedef cfla::L1NN1<coordinate_type, point_type> nn1_type;
  typedef typename nn1_type::size_type size_type;

private:
  typedef MaxRect<coordinate_type> solver_type;

  solver_type solver_;
  point_list users_;

  // Build an "L1 ball" (square) for each user point to pass to the maxrect
  // algorithm. The rectangles (squares) are formed using the L1 distance from
  // the user to its nearest facility.
  template<typename point_filter>
  void build_rects(const nn1_type& nn1, point_filter filter)
  {
    for (auto userp = begin(); userp != end(); ++userp)
    {
      // Only build rects for certain points.
      if (!filter(*userp))
        continue;

      point_type user = *userp;
      point_type site = nn1(user); /* nearest neighbor */

      /* This is actually the distance to the corner points from   o      o o
       * the center.  The left and right points become top-left  o---o =>  \
       * and bottom-right when rotated.                            o      o o */
      auto l1dist = nn1.distance(site, user);
      point_type tl = point_type(user.x - l1dist, user.y);
      point_type br = point_type(user.x + l1dist, user.y);
#if 0
//#ifdef DEBUG
      point_type tr = point_type(user.x, user.y + l1dist);
      point_type bl = point_type(user.x, user.y - l1dist);
      if (img_) {
        cv::line(*img_, tl, tr, cv::Scalar(0xcc, 0xcc, 0xcc), 1);
        cv::line(*img_, tr, br, cv::Scalar(0xcc, 0xcc, 0xcc), 1);
        cv::line(*img_, br, bl, cv::Scalar(0xcc, 0xcc, 0xcc), 1);
        cv::line(*img_, bl, tl, cv::Scalar(0xcc, 0xcc, 0xcc), 1);
        cv::imshow("rectangles", *img_);
        cv::waitKey(0);
      }
#endif
      tl = rotateZ2f_pos(tl);
      br = rotateZ2f_pos(br);

      // Build the upright rectangle.
      pure_rect_type l1rect(bp::construct<pure_rect_type>(tl.x, tl.y, br.x, br.y));
      solver_.add_event(l1rect);
    }
  }

  // Return the ideal location for player 2 among some player 1 facilities.
  template<typename point_filter>
  point_type compute(const nn1_type &nn1, point_filter filter)
  {
    // To solve for player 2, build L1 (rotated) rectangles from the
    // opposing player's facilities and solve for max-depth region(s).
    // We first rotate the space about the origin so we solve with axis-up
    // rectangles.
    build_rects(nn1, filter);
    solver_.compute();

    // We have potentially multiple solution cells.
    // Pick a randome one and return its center.
    if (solver_.size() == 0)
    {
      std::cerr << "warning: empty solution" << std::endl;
      return point_type(0, 0);
    }

    size_type chosen_one = randrange(size_type(0), solver_.size()-1);
#ifdef DEBUG
    // Dump all possible solutions.
    std::cerr << "solutions:" << std::endl;
    size_type cellidx = size_type(0);
    for (auto cellp = solver_.begin(); cellp != solver_.end(); ++cellp)
    {
      pure_rect_type solution = *cellp;
      point_type ctmp(0, 0);
      bp::center(ctmp, solution);
      std::cerr << "[" << std::setw(2) << std::setfill(' ') << std::dec
        << cellidx << "] (depth=" << solver_.solution(cellidx).depth()
        << ") " << ctmp << std::endl;
    }
#endif

    // Instead of the center, choose a random point in the cell.
    pure_rect_type solution = solver_.cell(chosen_one);

    coordinate_type randx = randrange(
        bp::get(solution, bp::HORIZONTAL, bp::LOW),
        bp::get(solution, bp::HORIZONTAL, bp::HIGH));
    coordinate_type randy = randrange(
        bp::get(solution, bp::VERTICAL, bp::LOW),
        bp::get(solution, bp::VERTICAL, bp::HIGH));
    point_type next(randx, randy);;
#ifdef DEBUG
    // Dump our chosen solution.
    std::cerr << "chose [" << std::setw(2) << std::setfill(' ') << std::dec
      << chosen_one << "] " << next << std::endl;
#endif

    // Make sure to rotate the solution point back to the original space, since
    // the solution space itself is rotated.
    return rotateZ2f_neg(next);
  }


public:

  template<typename UserIter>
    MaxRectSolver(UserIter begin, UserIter end)
      : solver_(), users_(begin, end)
    {}

  iterator begin(void) const { return users_.begin(); }
  iterator end(void) const { return users_.end(); }
  //const_iterator cbegin(void) const { return users_.cbegin(); }
  //const_iterator cend(void) const { return users_.cend(); }

  template<typename point_filter>
  point_type operator()(const nn1_type &facilities, point_filter filter) {
    return compute(facilities, filter);
  }

  inline static coordinate_type
    distance(const point_type &p1, const point_type &p2) {
      return nn1_type::distance(p1, p2);
    }

};

extern template class MaxRect<double>;
extern template class MaxRect<float>;

#ifdef DEBUG
template<typename U>
std::ostream& operator<<(std::ostream& os, cfla::rect::Edge<U> const& e) {
  os << "<[" << std::setw(2) << std::setfill(' ') << e.rect_index
    << "] " << (e.dir == boost::polygon::LOW ? "LOW " : "HIGH")
    << " " << e.coord << " d=" << e.depth << ">";
  return os;
}

template<typename U>
std::ostream& operator<<(std::ostream& os, cfla::rect::SolutionEdge<U> const& e)
{
  os << "[" << std::setw(2) << std::setfill(' ') << e.solution << "] from "
    << e.edge;
  return os;
}
#endif

} } // end namespace cfla::rect

template<typename Tp_>
std::ostream& operator<<(std::ostream& os,
    boost::polygon::rectangle_data<Tp_> const& r)
{
  namespace bp = boost::polygon;
  os << "|" << bp::get(r, bp::HORIZONTAL, bp::LOW) << ","
            << bp::get(r, bp::VERTICAL, bp::HIGH)
     << "|";
  return os;
}

