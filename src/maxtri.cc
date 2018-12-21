#include <stdexcept>
#include <algorithm>
#include <functional>
#include <iterator>

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
check_intersections(const Iter segment, const Iter begin, const Iter end)
{
  for (Iter sit = begin; sit != end && can_intersect(*segment, *sit); ++sit)
  {
    // Check for intersections and queue any we find.
    isect_point_type intersection;
    char code;
    if ('0' != (code = intersect(segment->edge(), sit->edge(), intersection)))
    {
#ifdef MAXTRI_DEBUG
      cerr << "    found intersection '" << code << "' "
        << intersection << " with " << *sit << endl;
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
  status_iterator eit = elb;

  // Look for intersections with edges in the status to the...

  // ... left
  if (getx(e.second()) < getx(e.first()))
    check_intersections(status_riterator(elb), ++status_riterator(elb),
        status.rend());

  // ... right
  else
    check_intersections(elb, ++status_iterator(elb), status.end());
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
  assert(elb != status.end());
  status.erase(elb);
}

template<class Tp_>
void MaxTri<Tp_>::
handle_intersection(isect_point_type const& isect)
{
  // Swap the edges that intersect in the status. TODO
  status_seg_type lseg = isect.make_segment(bp::LEFT);
  status_seg_type rseg = isect.make_segment(bp::RIGHT);

  status_iterator lb = status.lower_bound(lseg);
}

template<class Tp_>
void MaxTri<Tp_>::
handle_tripoint(tri_point_type const& tpoint)
{
  // We can do this by just queueing all points twice, but this way is clearer
  // and simplifies event ordering. See the diagrams in the header for
  // reference on which points are which.
  switch (tpoint.index)
  {
  case 0:
    insert_edge(tpoint.left_edge());
    insert_edge(tpoint.right_edge());
    break;
  case 1:
    remove_edge(tpoint.left_edge());
    insert_edge(tpoint.right_edge());
    break;
  case 2:
    remove_edge(tpoint.left_edge());
    remove_edge(tpoint.right_edge());
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
    cerr << "    [" << setw(2) << i << "] " << segment << endl;
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
