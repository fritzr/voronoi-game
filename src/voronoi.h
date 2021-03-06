#pragma once

#include <iterator>
#include <vector>
#include <algorithm>
#include <cmath>

#include "intersection.h" // leftTurn
#include "boost_geo_poly.h"

namespace voronoi
{

namespace bp = boost::polygon;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

template<class Tp_> class edge_iterator;

// A simple edge wrapper, consisting of two finite points.
// This simplifies edge_t, because
template <class Tp_>
struct Edge {
  typedef boost::geometry::model::d2::point_xy<Tp_> point_type;
  point_type p0, p1;
};

// Voronoi diagram wrapper for convenience with Voronoi game usage.
template<class Tp_>
class VoronoiDiagram
{
  friend class edge_iterator<Tp_>;
public:
  typedef typename Edge<Tp_>::point_type point_type;

  typedef typename bp::voronoi_diagram<Tp_> voronoi_diagram;
  typedef typename voronoi_diagram::edge_type edge_type;
  typedef typename voronoi_diagram::cell_type cell_type;
  typedef typename voronoi_diagram::vertex_type vertex_type;
  typedef typename voronoi_diagram::coordinate_type coordinate_type;
  typedef typename bp::rectangle_data<Tp_> rect_type;
  typedef typename rect_type::interval_type interval_type;
  typedef Edge<Tp_> finite_edge_type;

  typedef cell_type value_type;
  typedef ptrdiff_t difference_type;
  typedef size_t size_type;
  typedef value_type* pointer;
  typedef value_type& reference;
  typedef value_type const& const_reference;

  typedef typename voronoi_diagram::const_cell_iterator iterator;
  typedef edge_iterator<Tp_> finite_edge_iterator;

  typedef typename std::vector<point_type> container_type;
  typedef typename container_type::const_iterator point_iterator;

private:
  container_type sites;
  container_type users;
  // Used to clip infinite edges of the VD: y bounds are ignored
  double boundsx;

  // Map users to their closest site, or -1 if unknown.
  std::vector<int> u2s;

  voronoi_diagram vd;

  void clip_infinite_edge(const edge_type& edge, point_type& p0, point_type& p1)
    const;

  // Bounding rectangle for a cell.
  rect_type boundingRect(const cell_type &cell) const {
    rect_type bbox;
    bool bbox_init = false;
    for (auto eit = cell_begin(cell); eit != cell_end(cell); ++eit) {
      finite_edge_type e = *eit;
      if (!bbox_init) {
        bp::set_points(bbox, e.p0, e.p1);
        bbox_init = true;
      } else {
        bp::encompass(bbox, e.p0);
        bp::encompass(bbox, e.p1);
      }
    }
    return bbox;
  }

  // Workhorse(s) for build() (only one of these is implemented at a time)
  void map_knn(void);
  void map_quick(void);
  void map_slow(void);

public:
  // Initialize a Voronoi diagram from the given input facilities
  // (from which point_type must be constructable).
  // The max width and height are used to clip infinite edges of the diagram.
  // To actually use the diagram, you must call build().
  template<class SiteInputIter, class UserInputIter>
  VoronoiDiagram(SiteInputIter sites_begin, SiteInputIter sites_end,
      UserInputIter users_begin, UserInputIter users_end,
      double max_width=1e9, double max_height=1e9)
    : sites(sites_begin, sites_end), users(users_begin, users_end),
      boundsx(max_width), u2s(users.size(), -1), vd()
  {
  }

  template<class UserInputIter>
  VoronoiDiagram(UserInputIter users_begin, UserInputIter users_end,
      double max_width=1e9, double max_height=1e9)
    : sites(), users(users_begin, users_end),
      boundsx(max_width), u2s(users.size(), -1), vd()
  {
  }

  // Add sites or users one-by-one. Note that if you do this you wil need
  // to re-run build() before the answers will be correct.
  template<class Point>
    inline void add_site(Point const& site) { sites.push_back(site); }
  template<class InputIter> 
    inline void add_sites(InputIter begin, InputIter end)
    { sites.push_back(begin, end); }

  template<class Point>
    inline void add_user(Point const& user) {
      users.push_back(user);
      u2s.push_back(-1);
    }
  template<class InputIter>
    inline void add_users(InputIter begin, InputIter end) {
      users.push_back(begin, end);
      while (u2s.size() < users.size())
        u2s.push_back(-1);
    }

  // Return the coordinate of the site to which the given user belongs in O(1).
  inline point_type user_site(size_t user_index) const {
    return sites.at(u2s.at(user_index)); // bounds-checked
  }
  inline point_type user_site(point_iterator user_iter) const {
    return site(user_iter - users_begin());
  }

  // Return the index of the site to which the given user belongs in O(1).
  inline size_t site_index(size_t user_index) const {
    return u2s.at(user_index); // bounds-checked
  }
  inline size_t site_index(point_iterator user_iter) const {
    return site_index(user_iter - users_begin()); // bounds-checked
  }

  // Iterate over all cells.
  inline iterator begin() const { return vd.cells().begin(); }
  inline iterator end() const { return vd.cells().end(); }
  inline size_t size() const { return vd.num_cells(); }

  inline point_iterator sites_begin() const { return sites.begin(); }
  inline point_iterator sites_end() const { return sites.end(); }
  inline size_t sites_size() const { return sites.size(); }
  inline point_type const& site(size_t site_index) const {
    return sites.at(site_index); // bounds-checked
  }

  inline point_iterator users_begin() const { return users.begin(); }
  inline point_iterator users_end() const { return users.end(); }
  inline size_t users_size() const { return users.size(); }
  inline point_type const& user(size_t user_index) const {
    return users.at(user_index); // bounds-checked
  }
  inline point_type nearest_site(size_t user_index) const {
    return sites.at(site_index(user_index));
  }
  inline point_type nearest_site(point_iterator user_pt) const {
    return sites.at(site_index(user_pt));
  }

  // Iterate through finite points of the edges of a cell.
  inline finite_edge_iterator cell_begin(const cell_type &cell) const {
    return finite_edge_iterator(*this, cell.incident_edge());
  }
  inline finite_edge_iterator cell_end(const cell_type &cell) const {
    return finite_edge_iterator(*this, NULL);
  }

  // Whether the point is inside the cell.
  // Just make sure all the edges turn in the same direction towards pt.
  // This works because the cell must be convex.
  // The number of edges per cell are amortized O(1).
  inline bool is_inside(const cell_type& cell, const point_type& point) const
  {
    bool winding_known = false;
    bool winding = true;
    for (auto eit = cell_begin(cell); eit != cell_end(cell); ++eit) {
      finite_edge_type e = *eit;
      bool this_wind = leftTurn(e.p0, e.p1, point);
      if (winding_known && winding != this_wind)
        return false;
      winding = this_wind;
      winding_known = true;
    }
    return true;
  }

  enum SearchMethod {
    Default=0,
    Slow,
    Quick,
    KNN,
    SM_BAD,
  };

  // Build the Voronoi diagram and figure out which users belong to which
  // sites. This method performs two steps. The first step is to build
  // the Voronoi diagram out of the sites(), which takes O(m log m) time.
  // The second step is to map users to their closest sites. How the nearest
  // site is located depends on the SearchMethod parameter. The SearchMethod
  // algorithms have different run-times:
  // SearchMethod   Average case		Worst case
  //   Slow	    O(m * n)			O(m * n)
  //   Quick	    O((m + n) log n)		O(m * n)
  //   KNN*	    O(n (log n + log m))	O(n (log n + log m))
  // * The method marked with a '*' is the default.
  // For implementation details, see source code of the map_* methods.
  void build(SearchMethod method=Default);
};

template <class Tp_>
class edge_iterator
  : public std::iterator<std::bidirectional_iterator_tag, Edge<Tp_> >
{
private:
  typedef VoronoiDiagram<Tp_> parent_type;
  typedef typename parent_type::edge_type edge_type;
  typedef typename parent_type::finite_edge_type finite_edge_type;
  typedef finite_edge_type value_type;
  typedef ptrdiff_t difference_type;
  typedef size_t size_type;
  typedef value_type* pointer;
  typedef value_type const& reference;
  typedef value_type const& const_reference;

  const parent_type &parent;
  const edge_type *const first;
  const edge_type *current;

public:
  edge_iterator(const parent_type &p, const edge_type *e)
    : parent(p), first(e), current(e) {}

  inline bool operator==(edge_iterator const& o) const {
    return (current == NULL && o.current == NULL)
      || (first == o.first && current == o.current);
  }
  inline bool operator!=(edge_iterator const& o) const {
    return !(*this == o);
  }

  // postfix
  inline edge_iterator operator++(int) {
    edge_iterator it(parent, current);
    this->operator++();
    return it;
  }

  // prefix
  inline edge_iterator& operator++(void) {
    // If we reach the beginning again, set to NULL to identify the end
    if (current != NULL && (current = current->next()) == first)
      current = NULL;
    return *this;
  }

  inline value_type operator*(void) {
    value_type ret;
    if (current != NULL)
      parent.clip_infinite_edge(*current, ret.p0, ret.p1);
    return ret;
  }
};

// Common instantiations.
extern template class VoronoiDiagram<double>;

} // end namespace voronoi
