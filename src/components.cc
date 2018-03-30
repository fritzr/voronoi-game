#include <stdexcept>
#include <algorithm>
#include <functional>

#ifdef DEBUG
#include <iostream>
#include <iomanip>
#endif

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
int ConnectedComponents<Tp_>::
check_max_depth(rect_type const& r, EdgeSetIter edge_lb, int new_depth)
{
  if (new_depth <= max_depth)
    return new_depth;
  max_depth = new_depth;
  max_rects.clear();
  max_bottom_rect = -1;
  max_rects.insert(r.index);
  max_left = *edge_lb;
  max_rects.insert(max_left.rect_index);
#ifdef DEBUG
  cerr << "new max from rect " << r.index
    << " and edges: " << max_left.rect_index
    << (max_left.dir == bp::LOW ? "-" : "+");
#endif
  ++edge_lb;
  max_right = *edge_lb;
  max_rects.insert(max_right.rect_index);
#ifdef DEBUG
  cerr << " and " << max_right.rect_index
    << (max_right.dir == bp::LOW ? "-" : "+");
  cerr << endl;
#endif
  return new_depth;
}

template<class Tp_>
void ConnectedComponents<Tp_>::
insert_rect(rect_type const& r)
{
#ifdef DEBUG
  cerr << "inserting" << endl;
#endif
  // Insert edges in-place.
  edge_iterator elb = edges_x.insert(r.edge(bp::HORIZONTAL, bp::LOW)).first;
  edge_iterator ebefore = --edge_iterator(elb);
  edge_iterator eub = edges_x.insert(r.edge(bp::HORIZONTAL, bp::HIGH)).first;

  int depth = 0;
  if (elb != edges_x.begin())
  {
    depth = ebefore->depth;
    // If we were not inserted after a left edge, decrement depth
    // since we will erroneously increment it in the loop below.
    if (ebefore->dir != bp::LOW)
      --depth;
  }
  elb->depth = ++depth;
  ebefore = elb++;

  // Now increment depths of all edges inside us.
  while (elb != eub)
  {
    depth = ++(elb->depth);
    elb->depth = check_max_depth(r, elb, depth);
    // Mark that the rectangles intersect.
    b::add_edge(vd(r.index), vd((*elb).rect_index), graph);
    ++elb;
    ++ebefore;
  }

  // End edge gets the last depth. Decrement if we are after a right edge.
  if (ebefore->dir != bp::LOW)
    --depth;
  eub->depth = depth;

#ifdef DEBUG
  cerr << "inserted edges, status: " << edges_x << endl;
#endif
}

template<class Tp_>
void ConnectedComponents<Tp_>::
remove_rect(rect_type const& r)
{
#ifdef DEBUG
  cerr << "removing" << endl;
#endif

  auto lb = edges_x.find(r.edge(bp::HORIZONTAL, bp::LOW));
  auto ub = edges_x.find(r.edge(bp::HORIZONTAL, bp::HIGH));
  assert(lb != edges_x.end() && ub != edges_x.end());

  // Find the first location in our depth list to remove
  // and remove the corresponding depth and lower bound.
  edges_x.erase(lb++);
  // Now subtract one from all depths until the end.
  // Again, this should really be done by manipulating internal tree nodes.
  bool found_left=false, found_right=false;
  while (lb != ub)
  {
    // Now if we just identified a max-depth rectangle, we may not have seen
    // its bottom yet, since the bottom may be the top of this rectangle. The
    // first top-edge to intersect both the left/right edges of the current
    // max-rect must be the bottom of the actual max cell.
    if (max_bottom_rect < 0)
    {
      if (!found_left && (*lb == max_left))
        found_left = true;
      if (!found_right && (*lb == max_right))
        found_right = true;
    }

    --(lb->depth);
    ++lb;
    // ??? b::add_edge(vd(r.index), vd((*eit).rect_index), graph);
  }
  // Now erase the right edge.
  edges_x.erase(ub);

  // If we intersected the max-cell's vertical edges, we must be its bottom.
  if (found_left && found_right) {
    max_bottom_rect = r.index;
    max_rects.insert(max_bottom_rect);
  }

#ifdef DEBUG
  cerr << "removed edges, status: " << edges_x << endl;
#endif
}

template<class Tp_>
void ConnectedComponents<Tp_>::
compute(void)
{
  // Run through the queue of sorted horizontal edges.
  while (!edges_y.empty())
  {
    typename edge_queue::const_reference edge = edges_y.top();
#ifdef DEBUG
    cerr << "reading edge " << edge << endl;
#endif

    switch (edge.dir.to_int())
    {
    case bp::HIGH:
      // When we first encounter a rectangle, insert its vertical edges in the
      // sweep status.
      insert_rect(rect(edge));
      break;
    case bp::LOW:
      // When we encounter the top of a rectangle, we can remove its vertical
      // edges from our sweep status.
      remove_rect(rect(edge));
      break;
    default:
      throw runtime_error("unreachable");
    }

    edges_y.pop();
  }

  // Now if we found any max intersections, close them.
  if (!max_rects.empty())
  {
    auto rit = max_rects.begin();
    int idx = *rit++;
    max_rect = rects[idx];
    while (rit != max_rects.end())
    {
      idx = *rit++;
      /*assert(*/bp::intersect(max_rect, rects[idx])/*)*/;
    }
  }

#ifdef DEBUG
  cerr << "max rects are: " << max_rects << endl;
  cerr << "max cell is: " << max_rect << endl;
#endif
}

template class ConnectedComponents<double>;

} // end namespace components
