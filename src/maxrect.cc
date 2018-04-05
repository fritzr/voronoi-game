#include <stdexcept>
#include <algorithm>
#include <functional>

#ifdef DEBUG
#include <iostream>
#include <iomanip>
#endif

#include "maxrect.h"

using namespace std;

// competitive facility location algorithms
namespace cfla
{

template<class Tp_>
MaxRect<Tp_>::MaxRect(void)
  : rects(), edges_x(),
    component_ids(b::get(b::vertex_index_t(), graph)),
    vertexes(), max_depth(-1)
{
}

template<class Tp_>
MaxRect<Tp_>::~MaxRect(void)
{
}

template<class Tp_>
template<class EdgeSetIter>
int MaxRect<Tp_>::
check_max_depth(rect_type const& r, EdgeSetIter leftp, int new_depth)
{
  // Only bother with max depth cells, AND only if the edge of interest is
  // a left edge.
  if (new_depth <= 1 || new_depth < max_depth || leftp->dir != bp::LOW)
    return new_depth;

  // If we encounter a deeper cell, everything we've ever known until now
  // has been a lie.
  if (new_depth > max_depth)
  {
    max_depth = new_depth;
#ifdef DEBUG
    cerr << "solutions cleared, found new max depth " << max_depth << endl;
#endif
    solutions().clear();
    solution_edges.clear();
  }

  EdgeSetIter rightp = ++EdgeSetIter(leftp);
  // Track the new edges as they belong to the solution cell.
  // This is necessary so we can "efficiently" (see remarks on T(V_k))
  // see whether we have found the bottom of a solution cell in remove_rect().
  const int next_solution_index = solutions().size();
  solution_edges.emplace(*leftp, next_solution_index);
  solution_edges.emplace(*rightp, next_solution_index);
  solutions().emplace_back(r.index, leftp->rect_index, rightp->rect_index);

  return new_depth;
}

template<class Tp_>
void MaxRect<Tp_>::
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
  elb->depth = depth; // will be incremented in the loop below.

  // Now increment depths of all edges inside us.
  while (elb != eub)
  {
    depth = ++(elb->depth);
    elb->depth = check_max_depth(r, elb, depth);
    // Mark that the rectangles intersect.
    // The nature of the graph removes parallel (duplicate) edges.
    b::add_edge(vd(r.index), vd((*elb).rect_index), graph);
    ++elb;
  }

  // End edge gets the last depth. Decrement if we are after a right edge.
  ebefore = --edge_iterator(elb);
  if (ebefore->dir != bp::LOW)
    --depth;
  eub->depth = depth;

#ifdef DEBUG
  cerr << "inserted edges, status: " << edges_x << endl;
  cerr << "solution status: " << solution_edges << endl;
#endif
}

template<class Tp_>
void MaxRect<Tp_>::
remove_rect(rect_type const& r)
{
#ifdef DEBUG
  cerr << "removing" << endl;
#endif

  const edge_type left_edge  = r.edge(bp::HORIZONTAL, bp::LOW);
  const edge_type right_edge = r.edge(bp::HORIZONTAL, bp::HIGH);

  // Check any edges that lie inside the rect to remove.
  auto lb = edges_x.find(left_edge);
  auto ub = edges_x.find(right_edge);
  assert(lb != edges_x.end() && ub != edges_x.end());

  // Also check the edges of solution cells.
  auto sol_lb = solution_edges.lower_bound(left_edge);
  auto sol_ub = solution_edges.upper_bound(right_edge);

  // First find the closest solution cell to the real edge.
  while (sol_lb != sol_ub && *sol_lb < *lb)
    ++sol_lb;
  auto sol_it = sol_lb;

  // Subtract one from the depths of all edges we intersected.
  // Again, this should really be done by manipulating internal tree nodes.
  auto it = lb;
  while (it != ub)
  {
    // Check each solution cell to see whether the bottom edge of this
    // rectangle forms the bottom edge of the solution cell. We know this is
    // the case if the edges we are removing intersects
    if (sol_it != solution_edges.end())
    {
      // Mark any solution edge(s) that we intersect. This bottom edge is now
      // known as the bottom edge of the corresponding solution cell(s).
      if (*sol_it == *it) {
        solutions()[sol_it->solution].found(sol_it->edge.dir);
        ++sol_it;
      }
    }

    --(it->depth);
    ++it;
    // ??? b::add_edge(vd(r.index), vd((*eit).rect_index), graph);
  }
  // Now erase the left and right edges.
  edges_x.erase(lb);
  edges_x.erase(ub);

  // Now go through all the edges of solution cells we checked and see if we
  // are their bottom edge (we intersected one of their edges).
  sol_ub = sol_it;
  sol_it = sol_lb;
  while (sol_it != sol_ub)
  {
    solution_type& cell = solutions()[sol_it->solution];
    // We do not check our right edge in the loop above, so check now.
    if (*sol_it == right_edge)
      cell.found(sol_it->edge.dir);
    // If we did indeed intersect the solution cell, remove its edges
    // from the solution edge set since our sweep line has now passed the cell.
    if (cell.marked(r.index)) {
#ifdef DEBUG
      cerr << "found bottom for solution ["
        << setw(2) << setfill(' ') << sol_it->solution << "] from rect ["
        << setw(2) << setfill(' ') << r.index << "] at "
        << r.edge(bp::VERTICAL, bp::LOW) << endl;
#endif
      sol_it = solution_edges.erase(sol_it);
    }
    else
      ++sol_it;
  }

#ifdef DEBUG
  cerr << "removed edges, status: " << edges_x << endl;
  cerr << "solution status: " << solution_edges << endl;
#endif
}

template<class Tp_>
void MaxRect<Tp_>::
handle_event(bp::direction_1d const& dir, edge_type const& edge)
{
  switch (dir.to_int())
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
}

template<class Tp_>
void MaxRect<Tp_>::
add_rect(rect_type const& rect)
{
  size_t next_idx = rects.size();
  // Construct our custom rectangle wrappers.
  //vertex_descriptor vd = b::add_vertex(graph);
  vertexes.push_back(b::add_vertex(graph));
  b::put(component_ids, vd(next_idx), b::vertex_index_t(next_idx));
  // Add the rect and queue up the horizontal edges.
  // Nb. a horizontal edge consists of (VERTICAL,LOW) and (VERTICAL,HIGH)
  // coordinates in boost nonmenclature.
  rects.push_back(rect);
  rects.back().add_edges(push_inserter(queue()), bp::VERTICAL);
}

template<class Tp_>
void MaxRect<Tp_>::
initialize(void)
{
  // We can never have more solutions than rects so just reserve this upper
  // bound to prevent implicitly resizing the solutions vector.
  solutions().reserve(rects.size());
}

template<class Tp_>
void MaxRect<Tp_>::
finalize(void)
{
  // Cache the solution cells.
  solution_cells.clear();
  solution_cells.reserve(solutions().size());
  for (auto sol_it = solutions().begin(); sol_it != solutions().end(); ++sol_it)
    solution_cells.emplace_back(sol_it->cell(rects, rects.size()));
}

template class MaxRect<double>;

} // end namespace cfla
