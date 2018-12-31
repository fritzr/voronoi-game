#include <stdexcept>
#include <algorithm>
#include <functional>
#include <iterator>
#include <cmath>

#ifdef DEBUG
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
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
check_intersection(const Iter segment, Iter neighbor, const Iter end)
{
  if (segment == end || neighbor == end)
    return;

  const Iter first_neighbor = neighbor;
  point_type current_point = status_compare::isect_sweep(*segment, current_y());
  do
  {
    // Check for intersections and queue any we find.
    // Do not queue filthy degenerates or points we've already visited.
    isect_point_type intersection;
    char code = intersect(segment->edge(), neighbor->edge(), intersection);
    if ('0' != code && 'v' != code && 'e' != code
        && event_point_ycompare()(intersection, current_point))
    {
#ifdef MAXTRI_DEBUG
      cerr << "    found intersection '" << code << "' " << intersection
        << endl;
#endif
      // We can only queue intersection points we haven't visited yet.
      queue().push(intersection);
    }

    // Typically this isn't a loop, BUT we can have the exact same
    // segment across several different triangles. Such identical segments will
    // all be sorted adjacent, so we can loop through them linearly.
    // We could use a multimap, in which case we'd remove the triangle pointer
    // sorting criterion and pass iterators from equal_range() to this function.
    // But then we'd still have to write additional code to implement find()
    // (in other functions) based on the owning triangle.
  } while (status.key_comp().same_segment(*first_neighbor, *++neighbor));
}

template<class Tp_>
void MaxTri<Tp_>::
check_intersections(const status_iterator center)
{
  // Look left
  status_riterator rcenter(center);
  check_intersection(rcenter, ++status_riterator(rcenter), status.rend());

  // Look right
  check_intersection(center, ++status_iterator(center), status.end());
}

template<class Tp_>
void MaxTri<Tp_>::
insert_edge(edge_type const& e)
{
  // Insert the edge by its top point.
  status_iterator elb = status.emplace(e).first;

  // Look for intersections with adjacent edges in the status.
  check_intersections(elb);
}

template<class Tp_>
void MaxTri<Tp_>::
remove_edge(edge_type const& e)
{
  // Remove the edge from the status.
  status_seg_type edge_segment(e);
  status_iterator eub, elb = status.find(edge_segment);

  // Nb. If this is not the bottom edge ([2]) in the triangle, we will
  // immediately be adding another edge, but it is already done in
  // handle_tripoint.

#ifndef MAXTRI_DEBUG
  assert (elb != status.end());
#else
  // this should always be true: if it's not, handle_tripoint should catch it
  if (elb != status.end())
#endif
    eub = elb = status.erase(elb);
  --elb;

  // Now elb, eub refer to the edges immediately left and right of the removed
  // edge. Theoretically we only need to check for intersections between these
  // two edges; but because multiple consecutive edges can be from the same
  // triangle, we need to run the check (which involves a loop) in each
  // direction. However, we don't want to detect the same point twice.

  // Check to the right of LB, including LB & UB.
  check_intersection(elb, eub, status.end());

  // Check to the left of UB, skipping LB so we don't check LB & UB twice.
  if (elb != status.begin())
  {
    --elb;
    if (elb != status.begin())
    {
      --elb;
      check_intersection(eub, elb, status.rend());
    }
  }
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

#ifdef MAXTRI_DEBUG
  cerr << "    old  left: " << lseg_old << endl;
  cerr << "    old right: " << rseg_old << endl;
  cerr << "    new  left: " << lseg_new << endl;
  cerr << "    new right: " << rseg_new << endl;
#endif

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

  // Remove the old edges.
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

  // Update the sweep line to the intersection point now that we've removed the
  // intersecting edges. Now when we insert them they should be swapped
  // relative to their old positions.
  update_sweep(isect);

  // Perform the swap by inserting the edges again.
  // This time they will follow the new world order.
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
  // Update the sweep line to the point of this event.
  update_sweep(tpoint.value());

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

#ifdef MAXTRI_DEBUG

template<class Tp_>
void MaxTri<Tp_>::
dump_status_to_octave(ostream& os) const
{
  vector<coordinate_type> xs, ys;
  coordinate_type minx = numeric_limits<coordinate_type>::max()
                , maxx = numeric_limits<coordinate_type>::lowest();
  xs.reserve(status.size());
  ys.reserve(status.size());
  for (auto const& segment : status)
  {
    coordinate_type x1 = getx(segment.first()), x2 = getx(segment.second());

    if (x1 < minx)
      minx = x1;
    if (x1 > maxx)
      maxx = x1;
    if (x2 < minx)
      minx = x2;
    if (x2 > maxx)
      maxx = x2;

    xs.push_back(x1);
    xs.push_back(x2);
    ys.push_back(gety(segment.first()));
    ys.push_back(gety(segment.second()));
  }

  os << "figure();" << endl
    << "hold on;" << endl;
  // draw the sweep line itself if we have anything
  if (status.size() > 0)
  {
    os << "plot([ " << minx << " ; " << maxx << " ], "
               "[ " << current_y() << " ; " << current_y() << " ]"
               ", 'Color', 'red', 'LineWidth', 2);" << endl
       << "xlim([ " << minx << " ; " << maxx << " ]);" << endl;
  }

  // scatter plot for all endpoints
  os << "scatter([";
  for (auto const& x : xs)
    os << x << " ; ";
  os << "] , [";
  for (auto const& y : ys)
    os << y << " ; ";
  os << "]);" << endl;

  // also, plot each segment separately
  for (auto const& segment : status)
  {
    os << "plot(["
      << getx(segment.first()) << " ; "
      << getx(segment.second())
      << "], ["
      << gety(segment.first()) << " ; "
      << gety(segment.second())
      << "]);" << endl;
  }
}

template<class Tp_>
void MaxTri<Tp_>::
debug_status(void) const
{
  static int status_cnt = 0;

  stringstream ofname_s;
  ofname_s << "status_" << dec << setw(2) << setfill('0') << status_cnt << ".m";
  string ofname(ofname_s.str());

  cerr << ofname << ":" << endl;
  int i = 0;
  for (auto const& segment : status)
  {
    cerr << "    [" << setw(2) << dec << setfill(' ') << i << "] "
      << segment << endl;
    ++i;
  }

  ofstream ofstatus(ofname);
  dump_status_to_octave(ofstatus);
  ofstatus.close();

  ++status_cnt;
}

#endif

template<class Tp_>
void MaxTri<Tp_>::
handle_event(EventType type, event_point_type const& event)
{
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
  debug_status();
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
