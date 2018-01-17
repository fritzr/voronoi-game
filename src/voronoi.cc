#include "voronoi.h"

#include <boost/foreach.hpp>

using namespace std;
using namespace cv;

template<class Tp_>
void
VoronoiDiagram<Tp_>::clip_infinite_edge(const edge_type& edge,
    point_type& p0, point_type& p1)
  const
{
  const cell_type& cell1 = *edge.cell();
  const cell_type& cell2 = *edge.twin()->cell();
  point_type origin, direction;

  point_type s1 = sites[cell1.source_index()];
  point_type s2 = sites[cell2.source_index()];
  origin.x = (s1.x + s2.x) * 0.5;
  origin.y = (s1.y + s2.y) * 0.5;
  direction.x = (s1.y - s2.y);
  direction.y = (s2.x - s1.x);

  coordinate_type koef = boundsx
    / (std::max)(fabs(direction.x), fabs(direction.y));
  if (edge.vertex0() == NULL) {
    p0 = point_type(
        origin.x - direction.x * koef,
        origin.y - direction.y * koef);
  } else {
    p0 = point_type(edge.vertex0()->x(), edge.vertex0()->y());
  }
  if (edge.vertex1() == NULL) {
    p1 = point_type(
        origin.x + direction.x * koef,
        origin.y + direction.y * koef);
  } else {
    p1 = point_type(edge.vertex1()->x(), edge.vertex1()->y());
  }
}

// Function object to filter-out users whose sites are already known.
template <class Value>
struct is_unknown {
  const vector<int> &index_set;
  is_unknown(const vector<int> &s) : index_set(s) {}
  inline bool operator()(Value const& p) const {
    return index_set[p.second] == -1;
  }
};

template <class Value, class Tp_>
struct inside_cell {
  typedef typename VoronoiDiagram<Tp_>::cell_type cell_type;
  const VoronoiDiagram<Tp_> &vd;
  const cell_type &cell;
  inside_cell(const VoronoiDiagram<Tp_> &v, const cell_type& c)
    : vd(v), cell(c) {}
  inline bool operator()(Value const& p) const {
    return vd.is_inside(cell, p.first); // point inside cell
  }
};

// This is an efficient method for mapping user points to their cells using
// voronoi cells and a range-tree to quickly find the site to which each user
// point belongs.
// This runs in O(m log m + n log n) for m sites and n users.
// Currently this method doesn't work - the cell containment check using clip
// lines is incomplete since some cells are infinite.
// TODO use winding checks since Voronoi edge set is O(m) and cells are convex.
template<class Tp_>
void
VoronoiDiagram<Tp_>::build_fast(void)
{
#ifdef CELLS_SLOW
  throw runtime_error("Voronoi::build_fast(): method unimplmented");
#else
  // u2s maps user index to nearest facility index
  boost::polygon::construct_voronoi(sites.begin(), sites.end(), &vd);
  u2s = vector<int>(users.size(), -1);

  // Construct an R-tree for range queries on the user points.
  // We use this to efficiently find the user points that lie within each
  // voronoi cell using its bounding box.
  typedef pair<point_type, int> rvalue_t;
  vector<rvalue_t> uvals;
  int user_idx = 0;
  for (auto user = users.begin(); user != users.end(); ++user, ++user_idx)
    uvals.push_back(make_pair(point_type(user->x, user->y), user_idx));
  bgi::rtree<rvalue_t, bgi::quadratic<16> > utree(uvals.begin(), uvals.end());

  for (auto it = vd.cells().begin(); it != vd.cells().end(); ++it)
  {
    const cell_type &cell = *it;
    // Map the users that lie inside our cell to this site.
    rect_type bbox = boundingRect(cell);
    vector<rvalue_t> query_users;
    utree.query(
           bgi::satisfies(is_unknown<rvalue_t>(u2s))
        && bgi::within(bbox)
        && bgi::satisfies(inside_cell<rvalue_t, Tp_>(*this, cell)),
        back_inserter(query_users));
    BOOST_FOREACH(rvalue_t const& v, query_users)
    {
      int user_idx = v.second;
      u2s[user_idx] = cell.source_index();
    }
  }
#endif // !CELLS_SLOW
}

// This is a brute-force (slow) method for mapping user points to their cells
// simply by associating each user to the site to which it is nearest.
// This runs in O(m*n) for m sites and n users.
template<class Tp_>
void
VoronoiDiagram<Tp_>::build_slow(void)
{
#ifndef CELLS_SLOW
  throw runtime_error("Voronoi::build_fast(): method unimplmented");
#else // CELLS_SLOW
  u2s = vector<int>(users.size(), -1);
  for (size_t user_idx = 0u; user_idx < users.size(); ++user_idx)
  {
    const Point &user = users[user_idx];
    double mindist = std::numeric_limits<double>::infinity();
    int closest_site_idx = -1;
    for (size_t site_idx = 0u; site_idx < sites.size(); ++site_idx)
    {
      const Point &site = sites[site_idx];
      double distance = cv::norm(site - user);
      if (distance < mindist) {
        mindist = distance;
        closest_site_idx = static_cast<int>(site_idx);
      }
    }
    u2s[user_idx] = closest_site_idx;
  }
#endif // CELLS_SLOW
}

template<class Tp_>
void
VoronoiDiagram<Tp_>::build(void)
{
#ifdef CELLS_SLOW
  this->build_slow();
#else
  this->build_fast();
#endif
}

// Common instantiations.
template class VoronoiDiagram<double>;
