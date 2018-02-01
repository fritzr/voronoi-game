
// Tree structures with custom node data.
#include <set>
#include <queue>
#include <vector>

#include <opencv2/core/core.hpp>

#include <boost/polygon/polygon.hpp>
#include <boost/polygon/point_concept.hpp>
#include <boost/polygon/rectangle_concept.hpp>
#include <boost/geometry.hpp>

#include "tree.h"

namespace bp = boost::polygon;
namespace bg = boost::geometry;

template<class Tp_>
struct Edge
{
  typedef typename bp::rectangle_traits<Tp_> traits;
  typedef Tp_ coordinate_type;

  coordinate_type coord; // actual coordinate
  typename bp::direction_1d dir; // LOW or HIGH
  size_t rect_index; // parent rectangle

  Edge(coordinate_type const& c, bp::direction_1d const& d, size_t idx)
    : coord(c), dir(d), rect_index(idx) {}
};

template<class Tp_>
struct RectComponent : bp::rectangle_data<Tp_>
{
  typedef typename bp::rectangle_data<Tp_> super_type;
  typedef typename bp::rectangle_traits<Tp_> traits;
  typedef Edge<Tp_> edge_type;

  size_t index;
  int component; // ID of the component to which this rect belongs.

  template <class HInterval, class VInterval>
    RectComponent(const HInterval& hrange, const VInterval& vrange, size_t idx)
      : super_type(hrange, vrange), index(idx), component(-1) {}

  inline edge_type edge(bp::orientation_2d orient, bp::direction_1d dir)
  {
    return edge_type(bp::get(*this, orient, dir), dir, index);
  };

  template<class OutIter> void
    add_edges(OutIter out, bp::orientation_2d orient) {
      *out++ = edge(orient, bp::LOW);
      *out++ = edge(orient, bp::HIGH);
    }
};

namespace boost { namespace polygon {
  template <typename Tp_> struct geometry_concept<RectComponent<Tp_> > {
    typedef rectangle_concept type;
  };
}}


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

// Node to track 'cell depths'. These form the tree T(V_k) from the paper,
// the leaves of which represent g_k[j]: the depth of the cell at position
// j in the current sweep partition.
template<class Tp_>
class DepthNode : public tree::avl_node<DepthNode<Tp_> >
{
public:
  DECLARE_TRAITS(tree::avl_node_traits<DepthNode>);
  typedef Edge<Tp_> edge_type;

  // Internal depth value of the node.
  size_t depth;
  // Whether we are a leaf node.
  bool leaf(void) const {
    return (node_traits::get_left(this) == nullptr)
      && (node_traits::get_right(this) == nullptr);
  }
  // If we are a leaf node, which X edge we correspond to.
  edge_type leaf_edge;
};

template<class Tp_>
class ConnectedComponents
{
public:
  typedef Tp_ coordinate_type;
  typedef typename cv::Point_<Tp_> point_type;
  typedef RectComponent<Tp_> rect_type;
  typedef Edge<Tp_> edge_type;
  typedef EdgeCompare<Tp_> edge_comparator;

  typedef typename std::vector<edge_type> edge_container;
  typedef typename std::vector<rect_type> rect_container;
  typedef typename std::set<edge_type, edge_comparator>
    edge_set;
  typedef std::priority_queue<edge_type, edge_container, edge_comparator>
    edge_queue;
  /*
  typedef typename tree::avl_tree<DepthNode, typename DepthNode::compare>
    depth_tree; // T(V_k) from the paper
  */

private:
  rect_container rects;
  edge_queue edges_y; // horizontal rect edges which are the sweep events
  edge_set edges_x;   // vertical rect edges within each sweep line event

  inline rect_type &rect(edge_type const& edge) {
    return rects[edge.rect_index];
  }

public:
  // No inputs yet.
  ConnectedComponents();

  // Construct from a list of input rectangles (may be in no particular order).
  template<class RectIter>
    ConnectedComponents(RectIter begin, RectIter end);

  // Add rectangles. Same as the constructor form, if you're lazy and want to
  // do it sometime after construction.
  template<class RectIter>
    void add_rects(RectIter begin, RectIter end);

  // Run the algorithm and compute the connected components.
  void compute(void);
};

extern template class ConnectedComponents<double>;
