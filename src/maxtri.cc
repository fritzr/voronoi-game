#include <stdexcept>
#include <algorithm>
#include <functional>

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
  polys.clear();
}

template<class Tp_>
MaxTri<Tp_>::MaxTri(void)
  : tris(), polys(), edge_points(), max_depth(-1)
{
}

template<class Tp_>
void MaxTri<Tp_>::
insert_edge(edge_type const& e)
{
  // Insert the edge by its top point.
  edge_point_iterator elb = edge_points.insert(e.first()).first;
  edge_point_iterator eit = elb;

  // See if we overlap the left/right edges.
  --eit;
  bool intersects = true;
  int incr = -1;
  while (intersects && eit != edge_points.end())
  {
    // Check for intersection, but only if the edges are from different
    // polygon triangulations.
    if (&e.tridata()->poly != &eit->tridata()->poly)
    {
      point_type p;
      char isect_code = SegSegInt(
          static_cast<point_type const&>(e.first()),
          static_cast<point_type const&>(e.second()),
          static_cast<point_type const&>(*eit),
          static_cast<point_type const&>(eit->other()),
          p);
      if (isect_code == '0')
        intersects = false;

      // Queue the intersection point.
      else
      {
#ifdef MAXTRI_DEBUG
        cerr << "    found intersection " << p << " with "
          << *eit << endl;
#endif
        queue().push(isect_point_type(p, e.data(), eit->edata));
      }
    }
    if (incr < 0)
      --eit;
    else
      ++eit;

    // Start over increasing from insertion point.
    if ((!intersects || eit == edge_points.end()) && incr < 0)
    {
      eit = elb;
      ++eit;
      incr = 1;
    }
  }

}

template<class Tp_>
void MaxTri<Tp_>::
remove_edge(edge_type const& e)
{
  // TODO
}

template<class Tp_>
void MaxTri<Tp_>::
handle_intersection(isect_point_type const& ix)
{
  // TODO
}

template<class Tp_>
void MaxTri<Tp_>::
handle_event(bool intersection_event, event_point_type const& event)
{

  if (intersection_event)
  {
    isect_point_type const& isect(event.isect());

#ifdef MAXTRI_DEBUG
    cerr << "handling intersection of " << isect.e1->edge()
      << " and " << isect.e2->edge() << endl;
#endif

    handle_intersection(isect);
  }

  else
  {
    edge_point_type const& edge_point(event.edge());
    edge_type const& edge(edge_point.edge());

#ifdef MAXTRI_DEBUG
  cerr << "handling "
    << setw(6) << setfill(' ') << (edge_point.dir == bp::LOW ? "FIRST":"SECOND")
    << " point from "
    << setw(6) << setfill(' ') << (edge.dir() == bp::LOW ? "LEFT":"RIGHT" )
    << " edge" << edge << endl;
#endif

    switch (edge_point.dir.to_int())
    {
      case bp::LOW: // first (top) point of edge
        insert_edge(edge); break;

      case bp::HIGH: // second (bottom) point of edge
        remove_edge(edge); break;

      default:
        // unreachable
        throw std::runtime_error("bad int value for EdgePoint::dir!");
    }
  }
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
