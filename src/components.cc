#include <stdexcept>
#include <algorithm>
#include <functional>

#include "components.h"

using namespace std;

namespace components
{

template<class Tp_>
ConnectedComponents<Tp_>::ConnectedComponents(void)
  : rects(), edges_y(), edges_x(),
    component_ids(b::get(b::vertex_index_t(), graph))
{
}

template<class Tp_>
template<class EdgeSetIter>
void ConnectedComponents<Tp_>::
check_max_depth(rect_type const& r, EdgeSetIter edge_lb, int new_depth)
{
  if (new_depth <= max_depth)
    return;
  max_depth = new_depth;
  // TODO - once we find a max-depth... then what
}

template<class Tp_>
void ConnectedComponents<Tp_>::
insert_rect(rect_type const& r)
{
  // Insert edges in-place.
  auto lb = edges_x.insert(r.edge(bp::VERTICAL, bp::LOW)).first;
  auto ub = edges_x.insert(r.edge(bp::VERTICAL, bp::HIGH)).first;
  // Find the equivalent location in the depth list to insert
  // and insert the depth of the left edge, which is one more than before.
  int prev_depth = 0;
  auto depthit = sync_iters(depths.begin(), edges_x.begin(), lb);
  if (depthit != depths.end())
    prev_depth = *depthit++;
  depthit = ++depths.insert(depthit, prev_depth+1);
  // Now increment all depths covered by the rect.
  // TODO - the efficient way to do this is to use the tree T(V_k).
  // For now we are skipping the tree and using a list to avoid the work of
  // writing the custom 2-3 tree data-structure.
  while (depthit != depths.end() && lb != ub) {
    check_max_depth(r, lb, ++(*depthit++));
    b::add_edge(vd(r.index), vd((*lb).rect_index), graph);
    ++lb;
  }
  // Now insert the depth of the right edge.
  // The depth of all unvisited edges to the right are unchanged.
  depths.insert(depthit, prev_depth+1);
}

template<class Tp_>
void ConnectedComponents<Tp_>::
remove_rect(rect_type const& r)
{
  auto lb = edges_x.find(r.edge(bp::VERTICAL, bp::LOW));
  auto ub = edges_x.find(r.edge(bp::VERTICAL, bp::HIGH));
  assert(lb != edges_x.end() && ub != edges_x.end());
  // Find the first location in our depth list to remove
  // and remove the corresponding depth and lower bound.
  auto depthit = sync_iters(depths.begin(), edges_x.begin(), lb);
  depths.erase(depthit++);
  edges_x.erase(lb++);
  // Now subtract one from all depths until the end.
  // Again, this should really be done by manipulating internal tree nodes.
  while (lb != ub)
  {
    --(*depthit++);
    // ??? b::add_edge(vd(r.index), vd((*eit).rect_index), graph);
    ++lb;
  }
  // Now erase the right-edges.
  depths.erase(depthit);
  edges_x.erase(ub);
}

template<class Tp_>
void ConnectedComponents<Tp_>::
compute(void)
{
  // Run through the queue of sorted horizontal edges.
  while (!edges_y.empty())
  {
    typename edge_queue::const_reference edge = edges_y.top();

    switch (edge.dir.to_int())
    {
    case bp::LOW:
      // When we first encounter a rectangle, insert its vertical edges in the
      // sweep status.
      insert_rect(rect(edge));
      break;
    case bp::HIGH:
      // When we encounter the top of a rectangle, we can remove its vertical
      // edges from our sweep status.
      remove_rect(rect(edge));
      break;
    default:
      throw runtime_error("unreachable");
    }

    edges_y.pop();
  }
}

template class ConnectedComponents<double>;

} // end namespace components
