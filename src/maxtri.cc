#include <stdexcept>
#include <algorithm>
#include <functional>
#include <iterator>
#include <cmath>

#ifdef DEBUG
#include <iostream>
#include <iomanip>
#endif

#include "maxtri.h"
#include "intersection.h"

using namespace std;

// competitive facility location algorithms
namespace cfla { namespace tri
{

using ::operator<<;


template<class Tp_>
MaxTri<Tp_>::~MaxTri(void)
{
  status.clear();
  tris.clear();
  polys.clear();
}

template<class Tp_>
char MaxTri<Tp_>::
intersect(edge_type const& a, edge_type const& b, isect_point_type &result)
  const
{
  char c = '0';
  // Only consider intersections between edges from triangles belonging to
  // different polygons.
  if (&a.poly() != &b.poly())
  {
    c = SegSegInt(a.first(), a.second(), b.first(), b.second(), result);
    if (c != '0')
    {
      // Note the position is already set by SegSegInt.
      result.e1 = a.data();
      result.e2 = b.data();
    }
  }
  return c;
}

template<class Tp_>
template<typename Iter>
void
MaxTri<Tp_>::
check_intersection(const Iter segment, const Iter neighbor, const Iter end)
{
  if (segment != end && neighbor != end)
  {
    // Check for intersections and queue any we find.
    isect_point_type intersection;
    char code = intersect(segment->edge(), neighbor->edge(), intersection);
    if ('0' != code && 'v' != code && 'e' != code)
    {
#ifdef MAXTRI_DEBUG
      cerr << "    found intersection '" << code << "' "
        << intersection << " with " << *neighbor << endl;
#endif
      queue().push(intersection);
    }
  }
}

template<class Tp_>
void MaxTri<Tp_>::
insert_edge(edge_type const& e)
{
  // Insert the edge by its top point.
  status_iterator elb = status.emplace(e).first;

  // Look for intersections with edges in the status to the...

  // ... left
  check_intersection(status_riterator(elb), ++status_riterator(elb),
        status.rend());

  // ... right
  check_intersection(elb, ++status_iterator(elb), status.end());
}

template<class Tp_>
void MaxTri<Tp_>::
remove_edge(edge_type const& e)
{
  // Remove the edge from the status.
  // If this is not the bottom edge ([2]) in the triangle, we will immediately
  // be adding another edge, but it is handled by a separate event
  status_seg_type edge_segment(e);
  status_iterator elb = status.find(edge_segment);

#ifdef MAXTRI_DEBUG
  // this should always be true: if it's not, handle_tripoint will catch it
  if (elb != status.end())
#else
  assert (elb != status.end())
#endif
    status.erase(elb);
}

template<class Tp_>
void MaxTri<Tp_>::
handle_intersection(isect_point_type const& isect)
{
  // Swap the edges that intersect in the status.
  // The old left/right edges are now the new right/left edges.
  status_seg_type lseg_new = isect.make_segment(bp::LEFT);
  status_seg_type rseg_old(lseg_new.edge());
  status_seg_type rseg_new = isect.make_segment(bp::RIGHT);
  status_seg_type lseg_old(rseg_new.edge());

  status_iterator lb = status.find(lseg_old);
  status_iterator rb = status.find(rseg_old);
#ifdef MAXTRI_DEBUG
  if (lb == status.end())
    throw runtime_error("failed to find old left edge!");
  if (rb == status.end())
    throw runtime_error("failed to find old right edge!");

  size_t status_size = status.size();
#else
  assert(lb != status.end());
  assert(rb != status.end());
#endif

  // Perform the swap.
  lb = status.erase(lb);
#ifdef MAXTRI_DEBUG
  if (status.size() != status_size - 1)
    throw runtime_error("failed to remove old right edge from status!");
#endif
  --lb;

  rb = status.erase(rb);
#ifdef MAXTRI_DEBUG
  if (status.size() != status_size - 2)
    throw runtime_error("failed to remove old left edge from status!");
#endif

  rb = lb = status.insert(lb, lseg_new);
#ifdef MAXTRI_DEBUG
  if (status.size() != status_size - 1)
    throw runtime_error("failed to insert new left edge into status!");
#endif

  rb = status.insert(rb, rseg_new);
#ifdef MAXTRI_DEBUG
  if (status.size() != status_size)
    throw runtime_error("failed to insert new right edge into status!");
#endif

  // See if the new left edge intersects with its new left neighbor,
  // ditto for the new right edge and its right neighbor
  status_riterator lend = status.rend();
  check_intersection(status_riterator(lb), ++status_riterator(lb), lend);

  status_iterator rend = status.end();
  check_intersection(status_iterator(rb), ++status_iterator(rb), rend);
}

template<class Tp_>
void MaxTri<Tp_>::
handle_tripoint(tri_point_type const& tpoint)
{
  // We can do this by just queueing all points twice, but this way is clearer
  // and simplifies event ordering. See the diagrams in the header for
  // reference on which points are which.
#ifdef MAXTRI_DEBUG
  size_t status_size = status.size();
#endif
  switch (tpoint.index)
  {
  case 0:
    insert_edge(tpoint.left_edge());
#ifdef MAXTRI_DEBUG
    if (status.size() != status_size + 1)
      throw runtime_error("failed to uniquely insert left edge!");
#endif

    insert_edge(tpoint.right_edge());
#ifdef MAXTRI_DEBUG
    if (status.size() != status_size + 2)
      throw runtime_error("failed to uniquely insert right edge!");
#endif
    break;

  case 1:
    remove_edge(tpoint.left_edge());
    insert_edge(tpoint.right_edge());
#ifdef MAXTRI_DEBUG
    if (status.size() < status_size)
      throw runtime_error("failed to insert right edge!");
    if (status.size() > status_size)
      throw runtime_error("failed to remove left edge!");
#endif
    break;

  case 2:
    remove_edge(tpoint.left_edge());
#ifdef MAXTRI_DEBUG
    if (status.size() != status_size - 1)
      throw runtime_error("failed to remove left edge!");
#endif
    remove_edge(tpoint.right_edge());
#ifdef MAXTRI_DEBUG
    if (status.size() != status_size - 2)
      throw runtime_error("failed to remove right edge!");
#endif
    break;

  default:
    throw runtime_error("handle_tripoint: bad tpoint index!");
  }
}

template<class Tp_>
void MaxTri<Tp_>::
handle_event(EventType type, event_point_type const& event)
{
  _lasty = gety(event);

  switch (type)
  {
  case INTERSECTION:
    {
      isect_point_type const& isect(event.isect());

#ifdef MAXTRI_DEBUG
      cerr << "handling intersection " << isect << endl;
#endif

      handle_intersection(isect);
    }
    break;

  case TRIPOINT:
    {
      tri_point_type const& tri_point(event.point());

#ifdef MAXTRI_DEBUG
      cerr << "handling point " << tri_point << endl;
#endif

      handle_tripoint(tri_point);
    }
    break;

  default:
    throw std::runtime_error("unhandled event type!");
  }

#ifdef MAXTRI_DEBUG
  cerr << "status:" << endl;
  int i = 0;
  for (auto const& segment : status)
  {
    cerr << "    [" << setw(2) << dec << setfill(' ') << i << "] "
      << segment << endl;
    ++i;
  }
#endif
}

template<class Tp_>
void MaxTri<Tp_>::
initialize(void)
{
  // TODO
}

template<class Tp_>
void MaxTri<Tp_>::
finalize(void)
{
  // TODO
}

template class MaxTri<double>;
template class MaxTri<float>;

} } // end namespace cfla::tri
