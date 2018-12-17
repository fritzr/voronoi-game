#include <stdexcept>
#include <algorithm>
#include <functional>

#ifdef DEBUG
#include <iostream>
#include <iomanip>
#endif

#include "maxtri.h"

using namespace std;

// competitive facility location algorithms
namespace cfla { namespace tri
{

template<class Tp_>
MaxTri<Tp_>::~MaxTri(void)
{
}

template<class Tp_>
MaxTri<Tp_>::MaxTri(void)
  : tris(), edge_points(), max_depth(-1)
{
}

template<class Tp_>
void MaxTri<Tp_>::
insert_edge(edge_type const& e)
{
  // Insert the edge by its top point.
  edge_point_iterator elb = edge_points.insert(e.first()).first;

  // Then we need to queue an intersection point.
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
