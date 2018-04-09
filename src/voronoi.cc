#include "voronoi.h"

#include <boost/foreach.hpp>

using namespace std;
using namespace cv;

namespace voronoi
{

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

template<class Tp_>
void
VoronoiDiagram<Tp_>::map_knn(void)
{
  // Construct an R-tree for NN(1) queries on the sites.
  // We use this to efficiently find the site that is nearest to each user.
  typedef pair<point_type, int> rvalue_t;
  vector<rvalue_t> svals;
  int site_idx = 0;
  for (auto site = sites.begin(); site != sites.end(); ++site, ++site_idx)
    svals.push_back(make_pair(point_type(site->x, site->y), site_idx));
  typedef typename bgi::rtree<rvalue_t, bgi::quadratic<16> > voronoi_tree;
  voronoi_tree vtree(svals.begin(), svals.end());

  // The nearest-neighbor is of course the voronoi cell.
  for (size_t user_idx = 0u; user_idx < users.size(); user_idx++)
  {
    const point_type &user = users[user_idx];
    u2s[user_idx] = vtree.qbegin(bgi::nearest(user, 1))->second;
  }
}

// This is an efficient method for mapping user points to their nearest site,
// using a voronoi diagram and range-tree.
// This runs in O(m log m + (m + n) log n) on average for m sites and n users:
//   O(m log m) worst-case to build the voronoi diagram;
//   O(n log n) worst-case to build the user range tree;
//   O(m log n) average-case to find user points in each cell (see below).
// Since Voronoi cells are convex and the Voronoi edge set is O(m), we can
// use the winding check for containment of a query point in amortized O(1)
// time per cell.
// Each containment check in the worst case must consider O(n) users.
// We improve this to O(log n) on average by using a range tree to filter each
// user query down to points which lie in the cell's bounding box. Of course
// this only works when the cells and users are distributed randomly.
// In the worst case, all Voronoi cells may share the same bounding box, and
// all O(n) user points will be checked for all O(m) cells anyway. However, the
// containment check is short-circuited for those users who have already been
// assigned a cell yet, so the overhead is small for the extra comparisons.
template<class Tp_>
void
VoronoiDiagram<Tp_>::map_quick(void)
{
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
}

// This is a brute-force (slow) method for mapping user points to their cells
// simply by associating each user to the site to which it is nearest.
// This runs in O(m*n) for m sites and n users.
template<class Tp_>
void
VoronoiDiagram<Tp_>::map_slow(void)
{
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
}

template<class Tp_>
void
VoronoiDiagram<Tp_>::build(SearchMethod sm)
{
  vd.clear();

  // Always construct the voronoi diagram from the sites so we have edges.
  boost::polygon::construct_voronoi(sites.begin(), sites.end(), &vd);

  // u2s maps user index to nearest facility index
  u2s = vector<int>(users.size(), -1);

  // Now map users to their cells through whichever method we want.
  switch (sm) {
    case SearchMethod::Slow:
      this->map_slow(); break;
    case SearchMethod::Quick:
      this->map_quick(); break;
    case SearchMethod::KNN:
    case SearchMethod::Default:
      this->map_knn(); break;
    default:
      throw runtime_error("VoronoiDiagram::build(): bad SearchMethod!");
  };
}

// Common instantiations.
template class VoronoiDiagram<double>;

}
