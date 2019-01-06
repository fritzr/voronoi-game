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

template<class Tp_>
MaxTri<Tp_>::~MaxTri(void)
{
  status.clear();
  tris.clear();
  polys.clear();
}

template<class Tp_>
template<typename Iter>
void
MaxTri<Tp_>::
intersect_range(status_seg_type const& segment, Iter begin, Iter end)
{
  // Check intersection with each segment from [neighbor, end).
  point_type current_point = status_compare::isect_sweep(segment, current_y());

  // See the comments on IXEvent for why this is a loop.
  while (begin != end)
  {
#ifdef MAXTRI_DEBUG
    cerr << "    checking intersection with " << *begin << endl;
#endif

    // Check for intersections and queue any we find.
    // Do not queue filthy degenerates or points we've already visited.
    point_type pt;
    char code = intersect(segment.edge(), begin->edge(), pt);
    if ('0' != code && 'v' != code && 'e' != code
        && event_point_ycompare()(pt, current_point))
    {
#ifdef MAXTRI_DEBUG
      cerr << "    found intersection '" << code << "' " << pt << endl;
#endif

      // We need to look up the intersection in the intersection map to see if
      // we've already found an intersection here.
      auto ixit = intersections.find(pt);

      // If it doesn't exist yet, this is the first intersection at this point.
      if (ixit == intersections.end())
      {
        auto insert = intersections.insert(make_pair(pt, isect_event_type(pt)));
        assert(insert.second);
        ixit = insert.first; // ref to the newly inserted pair.

        // Queue the intersection the first time we find it.
        // Subsequent intersections at the same point will just add segments.
        queue().emplace(ixit->second);
      }

      // Record that these segments also form this intersection.
      ixit->second.insert(segment, *begin);
    }

    ++begin;
  }
}

// Check new intersections for a segment based on the current status.
template<class Tp_>
void MaxTri<Tp_>::
check_intersections(status_iterator center)
{
  // don't call us with a bad segment!
  assert(center != status.end());

#ifdef MAXTRI_DEBUG
  cerr << "    looking for intersections around " << *center << endl;
#endif

  auto range = status.equal_range(*center);

  // Innner bounds.
  status_iterator right = range.second;
  status_iterator left = range.first;

  /* The following diagram illustrates the iterators below referring to
   * segments in the sweep status, where the same letter indicates segments in
   * the same equal_range:
   *
   *      0  1  2      3  4  5      6  7  8      9 10 11     12 13 14
   *    ---------------------------------------------------------------
   *   |  X  X  X      L  L  L      C  C  C      R  R  R      Z  Z  Z  |
   *    ---------------------------------------------------------------
   *            ^      ^     ^      ^  ^  ^      ^            ^
   *     left_end      |     |   left  |  right  |            right_end
   *         left_last |     |       center      |
   *              left_begin |                   | right_begin
   *
   * Ultimately, we want to check *center for intersections with all the
   * left (L) and right (R) segments, but not the segments in its own range (C).
   * We use reverse iterators in the left direction, and forward iterators in
   * the right direction so the actual checking code in intersect_range can be
   * the same.
   */

  // Look for intersections with the left-adjacent edge(s).
  if (left != status.end() && left != status.begin())
  {
    // Move out of the equal_range on center into the range to the left.
    // Nb. converting to reverse_iterator subtracts 1 from the pointer.
    status_riterator left_begin = status_riterator(left);

    status_iterator left_last = status.lower_bound(*left_begin);
    status_riterator left_end = status.rend();
    if (left_last != status.end())
      left_end = status_riterator(left_last);

#ifdef MAXTRI_DEBUG
      cerr << "      left: from " << *left_begin << " up to ";
      if (left_end != status.rend())
        cerr << *left_end;
      else
        cerr << "REND" << endl;
#endif

    intersect_range(*center, left_begin, left_end);
  }

  // Look for intersections with the right-adjacent edge(s).
  if (right != status.end())
  {
    // Move out of the equal_range on center into the range to the right.
    status_iterator right_begin = right;
    ++right_begin;
    if (right_begin != status.end())
    {
      status_iterator right_end = status.upper_bound(*right_begin);
      if (right_end != status.end())
        ++right_end;

#ifdef MAXTRI_DEBUG
      cerr << "      right: from " << *right_begin << " up to ";
      if (right_end != status.end())
        cerr << *right_end;
      else
        cerr << "END" << endl;
#endif

      intersect_range(*center, right_begin, right_end);
    }
  }
}

template<class Tp_>
void MaxTri<Tp_>::
insert_edge(edge_type const& e)
{
#ifdef MAXTRI_DEBUG
  cerr << "  inserting " << e << endl;
#endif

  // Insert the edge by its top point.
  // It is always inserted as the upper bound of any equal_range.
  status_iterator eit = status.insert(e);
  assert(eit != status.end());

  check_intersections(eit);
}

template<class Tp_>
void MaxTri<Tp_>::
remove_edge(edge_type const& e)
{
#ifdef MAXTRI_DEBUG
  size_t status_size = status.size();
  cerr << "  removing " << e << endl;
#endif

  // Remove the edge from the status.
  status_seg_type edge_segment(e);
  status_iterator eub, elb = find_unique(edge_segment);

  // Nb. If this is not the bottom edge ([2]) in the triangle, we will
  // immediately be adding another edge, but it is already done in
  // handle_tripoint.

#ifndef MAXTRI_DEBUG
  assert (elb != status.end());
#else
  // this should always be true: if it's not, handle_tripoint should catch it
  if (elb == status.end())
    MTFAIL("edge not found!");
#endif

  eub = elb = status.erase(elb);
  --elb;

#ifdef MAXTRI_DEBUG
  if (status.size() != status_size - 1)
    MTFAIL("failed to remove edge!");
#endif

  // Now elb, eub refer to the edges immediately left and right of the removed
  // edge. Theoretically we only need to check for intersections between these
  // two edges, but pass them both to the helper anyway.
  // It's okay if we find the same intersection twice, because it should end up
  // collapsing in the unordered_set used by IXEvent.
  check_intersections(elb);
  check_intersections(eub);
}

template<class Tp_>
void MaxTri<Tp_>::
handle_intersection(point_type const& isect_point)
{
  // The first step is to reorder all the segments which form the intersection.
  //
  // Because we don't have access to the internal tree nodes of a C++ set,
  // we have to do this by removing the segments, modifying the comparator,
  // then re-inserting the segments so they'll be sorted correctly by the
  // updated comparator.
  //
  // This sounds like bad news, but it's actually fine, since the algorithm is
  // designed so that the change in comparison criteria happens only after each
  // intersection point, and affects only the intersecting segments! Therefore
  // we are not modifying the sort order of any nodes existing in the tree.

#ifdef MAXTRI_DEBUG
  size_t status_size;
#endif

  // Lookup the actual intersection event from its point.
  auto ixit = intersections.find(isect_point);
  assert(ixit != intersections.end());
  isect_event_type const& isect = ixit->second;

#ifdef MAXTRI_DEBUG
  cerr << "handling " << isect;
#endif

  // Remove the segments forming the intersection (so we can add them later).
  for (auto it = isect.begin(); it != isect.end(); ++it)
  {
#ifdef MAXTRI_DEBUG
    status_size = status.size();
    cerr << "    removing old segment " << *it << endl;
#endif

    status_iterator segit = find_unique(*it);

#ifdef MAXTRI_DEBUG
    if (segit == status.end())
      MTFAIL("failed to find old segment " << *it);
#endif

    status.erase(segit);

#ifdef MAXTRI_DEBUG
    if (status.size() != status_size - 1)
      MTFAIL("failed to remove old segment " << *it);
#endif
  }

  // Update the sweep line to the intersection point now that we've removed the
  // intersecting segments. This changes the sorting criteria, since the
  // comparator refers to us for the y-position of the sweep line.
  // Now when we insert the segments again they should be in the correct order.
  // Classically this equates to a single swap.
  update_sweep(isect_point);

  // Perform the reordering by inserting the intersecting segments again.
  // This time they will follow the new world order.
  // We need to remember each segment that is inserted, so we can check for
  // intersections after they're all added again.
  std::vector<status_iterator> new_segments;
  for (auto it = isect.begin(); it != isect.end(); ++it)
  {
#ifdef MAXTRI_DEBUG
    status_size = status.size();
    cerr << "    re-inserting segment " << *it << endl;
#endif

    status_iterator newit = status.insert(*it);
    new_segments.push_back(newit);

#ifdef MAXTRI_DEBUG
    if (status.size() != status_size + 1)
      MTFAIL("failed to insert new segment!");
#endif
  }

#ifdef MAXTRI_DEBUG
  cerr << "reordered ";
  debug_status();
#endif

  // Now that the sweep line has been updated to the intersection point and the
  // intersecting segments have been reordered, we have to check each of these
  // segments for their next point of intersection below the sweep line.
  for (auto const& segiter : new_segments)
    check_intersections(segiter);
}

template<class Tp_>
void MaxTri<Tp_>::
handle_tripoint(tri_point_type const& tpoint)
{
  // Update the sweep line to the point of this event.
  // Though this implicitly modifies the status comparator, no segments can be
  // reordered in the status as a result because of the invariants of the
  // algorithm.
  update_sweep(tpoint.value());

  // We could queue all points twice (once for the top of each segment, once
  // for the bottom of each segment) but handling each point of the triangle
  // is clearer and simplifies event ordering. The diagram from the header
  // is reproduced below for reference.
  switch (tpoint.height)               /*                                 */
  {                                    /*       [0]                TOP    */
  case TOP:                            /* (0,1)  *                        */
    insert_edge(tpoint.top_edge());    /* V E0  /  \  E1 V                */
    insert_edge(tpoint.bottom_edge()); /*      /     \  (0,2)             */
    break;                             /* [1] *- _    \            MIDDLE */
                                       /*         - _   \                 */
  case MIDDLE:                         /*             - _\                */
    remove_edge(tpoint.top_edge());    /*        E2 >     * [2]    BOTTOM */
    insert_edge(tpoint.bottom_edge()); /*        (1,2) or (2,1)           */
    break;

  case BOTTOM:
    remove_edge(tpoint.top_edge());
    remove_edge(tpoint.bottom_edge());
    break;

  default:
    MTFAIL("bad tpoint index!");
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
      point_type const& isect_point(event.isect());

      handle_intersection(isect_point);
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
    MTFAIL("unhandled event type!");
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
