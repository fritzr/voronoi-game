
// Tree structures with custom node data.
#include <set>
#include <queue>
#include <vector>
#include <unordered_set>

#ifdef DEBUG
#include <iostream>
#include <iomanip>
#endif

#include <opencv2/core/core.hpp>

#include <boost/polygon/polygon.hpp>
#include <boost/polygon/point_concept.hpp>
#include <boost/polygon/rectangle_concept.hpp>
#include <boost/geometry.hpp>

#include <boost/graph/adjacency_list.hpp>

//#include "tree.h"

template <class Container>
  class push_insert_iterator:
    public std::iterator<std::output_iterator_tag,void,void,void,void>
{
protected:
  Container* container;

public:
  typedef Container container_type;
  explicit push_insert_iterator(Container& x) : container(&x) {}
  push_insert_iterator<Container>& operator=(
      typename Container::const_reference value)
    { container->push(value); return *this; }
  push_insert_iterator<Container>& operator* (){ return *this; }
  push_insert_iterator<Container>& operator++ (){ return *this; }
  push_insert_iterator<Container> operator++ (int){ return *this; }
};

template<typename Container>
push_insert_iterator<Container> push_inserter(Container& container){
  return push_insert_iterator<Container>(container);
}

namespace components
{

namespace b  = boost;
namespace bp = boost::polygon;
namespace bg = boost::geometry;

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

#ifdef DEBUG
  template<typename U>
    friend std::ostream& operator<<(std::ostream& os, Edge<U> const& e) {
      os << "<d=" << e.depth << " " << (e.dir == bp::LOW ? "LOW " : "HIGH")
         << " [" << std::setw(2) << std::setfill(' ') << e.rect_index
         << "] " << e.coord << ">";
      return os;
    }
#endif
};

template<class Tp_>
struct RectComponent : bp::rectangle_data<Tp_>
{
  typedef typename bp::rectangle_data<Tp_> super_type;
  typedef typename bp::rectangle_traits<Tp_> traits;
  typedef Edge<Tp_> edge_type;

  size_t index;
  //int component; // ID of the component to which this rect belongs.

  template <class HInterval, class VInterval>
    RectComponent(const HInterval& hrange, const VInterval& vrange, size_t idx)
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

} // we interrupt your regularly scheduled components namespace to bring you...

namespace boost { namespace polygon {
  template <typename Tp_>
    struct geometry_concept<components::RectComponent<Tp_> >
      { typedef rectangle_concept type; };
}}

#ifdef DEBUG
template<typename T, typename C>
  std::ostream& operator<<(std::ostream& os, std::set<T, C> const& s) {
    os << "set< ";
    for (auto it = s.begin(); it != s.end(); ++it)
      os << *it << ", ";
    os << " >";
    return os;
  }

template<typename T>
  std::ostream& operator<<(std::ostream& os, std::unordered_set<T> const& s) {
    os << "set{ ";
    for (auto it = s.begin(); it != s.end(); ++it)
      os << *it << ", ";
    os << " }";
    return os;
  }

template<typename T>
  std::ostream& operator<<(std::ostream& os, std::vector<T> const& v) {
    os << "[ ";
    for (auto it = v.begin(); it != v.end(); ++it)
      os << *it << ", ";
    os << " ]";
    return os;
  }

template<typename T>
  std::ostream& operator<<(std::ostream& os, std::list<T> const& l) {
    os << "[ ";
    for (auto it = l.begin(); it != l.end(); ++it)
      os << *it << ", ";
    os << " ]";
    return os;
  }
#endif

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


namespace components
{

template<class Tp_>
struct EdgeCompare
{
  typedef Edge<Tp_> edge_type;
  bool operator()(edge_type const& e1, edge_type const& e2)
  {
    return e1.coord < e2.coord
      || ((e1.coord == e2.coord) && (e1.rect_index < e2.rect_index));
  }
};

template<class Node_>
struct depth_compare
{
  typedef Node_ node;
  inline bool operator()(const node &d1, const node &d2) const
  {
    return d1.depth < d2.depth;
  }
};

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

template<class Iter1, class Iter2>
Iter1 sync_iters(Iter1 begin1, Iter2 begin2, Iter2 end2) {
  while (begin2++ != end2)
    ++begin1;
  return begin1;
}

template<class Iter1, class Iter2>
Iter1 rsync_iters(Iter1 end1, Iter2 end2, Iter2 begin2) {
  while (end2-- != begin2)
    --end1;
  return end1;
}

template<class Tp_>
class ConnectedComponents
{
public:
  typedef Tp_ coordinate_type;
  typedef typename cv::Point_<Tp_> point_type;
  typedef RectComponent<Tp_> rect_type;
  typedef typename rect_type::super_type pure_rect_type;
  typedef Edge<Tp_> edge_type;
  typedef EdgeCompare<Tp_> edge_comparator;

  typedef typename std::vector<edge_type> edge_container;
  typedef typename std::vector<rect_type> rect_container;
  typedef typename std::set<edge_type, edge_comparator>
    edge_set;
  typedef std::priority_queue<edge_type, edge_container, edge_comparator>
    edge_queue;
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

  typedef typename std::unordered_set<size_t> index_set;
  typedef typename index_set::const_iterator index_iterator;

  typedef typename edge_set::iterator edge_iterator;

private:
  rect_container rects;
  edge_queue edges_y; // horizontal rect edges which are the sweep events
  edge_set edges_x;   // vertical rect edges within each sweep line event
  //depth_tree_type depth_tree;
  components_graph graph;
  id_map component_ids;
  descriptor_list vertexes;
  int max_depth = -1;
  bool max_flag = false;
  // indexes of rects which form the maximal intersection
  index_set max_rects;
  pure_rect_type max_rect;

  inline vertex_descriptor vd(int idx) const { return vertexes[idx]; }

  inline rect_type &rect(edge_type const& edge) {
    return rects[edge.rect_index];
  }

  template<typename EdgeSetIter>
    int check_max_depth(rect_type const& r, EdgeSetIter edge_lb, int depth);

  void insert_rect(rect_type const& r);
  void remove_rect(rect_type const& r);

public:
  // No inputs yet.
  ConnectedComponents();

  inline int index(vertex_descriptor v) const {
    return b::get(component_ids, v);
  }

  // Construct from a list of input rectangles (may be in no particular order).
  template<class RectIter>
    ConnectedComponents(RectIter begin, RectIter end)
    : rects(), edges_y(), edges_x(), graph(),
      component_ids(b::get(b::vertex_index_t(), graph))
    {
      add_rects(begin, end);
    }

  components_graph const& adj_graph(void) const { return graph; }

  // Add rectangles. Same as the constructor form, if you're lazy and want to
  // do it sometime after construction.
  template<class RectIter>
    void add_rects(RectIter begin, RectIter end)
    {
      size_t idx = rects.size();
      while (begin != end)
      {
        // Construct our custom rectangle wrappers.
        //vertex_descriptor vd = b::add_vertex(graph);
        vertexes.push_back(b::add_vertex(graph));
        b::put(component_ids, vd(idx), b::vertex_index_t(idx));
        // Make sure our intervals are in the proper order, otherwise
        // the algorithm WILL fail.
        auto hivl = bp::get(*begin, bp::HORIZONTAL);
        auto hl = bp::get(hivl, bp::LOW), hh = bp::get(hivl, bp::HIGH);
        auto vivl = bp::get(*begin, bp::VERTICAL);
        auto vl = bp::get(vivl, bp::LOW), vh = bp::get(vivl, bp::HIGH);
        hivl = bp::construct<decltype(hivl)>(std::min(hl,hh), std::max(hl,hh));
        vivl = bp::construct<decltype(vivl)>(std::min(vl,vh), std::max(vl,vh));
        rects.emplace_back(hivl, vivl, idx++);
        ++begin;
        // Queue up the horizontal edges. Vertical edges go in at each event.
        rects.back().add_edges(push_inserter(edges_y), bp::HORIZONTAL);
      }
    }

  // Run the algorithm and compute the connected components.
  void compute(void);

  inline int depth(void) const { return max_depth; }

  // Return the rectangle with maximal depth.
  inline pure_rect_type max(void) const { return max_rect; }

  inline index_iterator begin(void) const { return max_rects.cbegin(); }
  inline index_iterator end(void) const { return max_rects.cend(); }
};

extern template class ConnectedComponents<double>;

} // end namespace components
