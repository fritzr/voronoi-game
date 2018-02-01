#include <stdexcept>
#include <algorithm>
#include <functional>

#include "components.h"

using namespace std;

template<class Tp_>
ConnectedComponents<Tp_>::ConnectedComponents(void)
  : rects(), edges_y(), edges_x()
{
}

template<class Tp_>
template<class RectIter>
ConnectedComponents<Tp_>::ConnectedComponents(RectIter begin, RectIter end)
  : rects(), edges_y(), edges_x()
{
  add_rects(begin, end);
}

template<class Tp_>
template<class RectIter>
void ConnectedComponents<Tp_>::
add_rects(RectIter it, RectIter end)
{
  size_t idx = rects.size();
  while (it != end)
  {
    // Construct our custom rectangle wrappers.
    rects.emplace_back(
        bp::get(*it, bp::HORIZONTAL),
        bp::get(*it, bp::VERTICAL),
        idx++);
    ++it;
    // Queue up the horizontal edges. Vertical edges go in at each event.
    rects.back().add_edges(back_inserter(edges_y), bp::HORIZONTAL);
  }
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
      rect(edge).add_edges(inserter(edges_x, edges_x.end()), bp::VERTICAL);
      break;
    case bp::HIGH:
      // When we encounter the top of a rectangle, we can remove its vertical
      // edges from our sweep status.
      edges_x.erase(rect(edge).edge(bp::VERTICAL, bp::LOW));
      edges_x.erase(rect(edge).edge(bp::VERTICAL, bp::HIGH));
      break;
    default:
      throw runtime_error("unreachable");
    }

    edges_y.pop();
  }
}

template class ConnectedComponents<double>;
