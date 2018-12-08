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
  : tris(), edge_points(), graph(),
    component_ids(b::get(b::vertex_index_t(), graph)),
    vertexes(), max_depth(-1)
{
}

template<class Tp_>
void MaxTri<Tp_>::
insert_edge(edge_type const& e)
{
  // TODO
}

template<class Tp_>
void MaxTri<Tp_>::
remove_edge(edge_type const& e)
{
  // TODO
}

template<class Tp_>
void MaxTri<Tp_>::
handle_intersection(edge_type const& e1, edge_type const& e2)
{
  // TODO
}

template<class Tp_>
void MaxTri<Tp_>::
handle_event(edge_type const& edge, edge_point_type const& event)
{
  // TODO
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
