#pragma once

#ifdef DEBUG
#include <iostream>
#endif

#include <cmath> // std::abs

#include <boost/geometry/geometry.hpp>
#include <boost/geometry/algorithms/distance.hpp>
#include <boost/geometry/algorithms/comparable_distance.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "boost_geo_poly.h"

namespace cfla
{
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

template<class coordinate_type, class point_type>
inline coordinate_type
l1dist(const point_type &p1, const point_type &p2)
{
  return std::abs(bg::get<0>(p1) - bg::get<0>(p2))
       + std::abs(bg::get<1>(p1) - bg::get<1>(p2));
}

template<class coordinate_type, class point_type>
inline coordinate_type
l2dist(const point_type &p1, const point_type &p2)
{
  return bg::distance(p1, p2);
}

/* L2 (Euclidean) nearest-neighbor 1 adapter.  */
template<class coordinate_type, class point_type>
class L2NN1
{
public:
  typedef typename bgi::rtree<point_type, bgi::quadratic<16> > nn_tree;
  typedef typename nn_tree::value_type     value_type;
  typedef typename nn_tree::const_iterator const_iterator;
  typedef typename nn_tree::const_iterator iterator;
  typedef typename nn_tree::size_type      size_type;

  // Members
private:
  // Data structure for locating the nearest facility to each user point.
  nn_tree nn_;

public:
  template<class InputIter>
  L2NN1(InputIter ptbegin, InputIter ptend)
    : nn_(ptbegin, ptend)
  {}

  // Iterate over input points.
  inline const_iterator begin(void) const { return nn_.begin(); }
  inline const_iterator end(void) const { return nn_.end(); }
  inline size_type size(void) const { return nn_.size(); }
  inline void insert(point_type const& f) { nn_.insert(f); }

  /* TODO paramaterize on L(n) to make L1/L2/etc easier.  */
  inline point_type const& operator()(point_type const& pt) const {
    return *nn_.qbegin(bgi::nearest(pt, 1));
  }

  /* comparable_distance() uses pythagoras/euclidean (L2) by default. */
  inline static coordinate_type comparable_distance(point_type p1, point_type p2)
    { return bg::comparable_distance(p1, p2); }

  inline static coordinate_type distance(point_type p1, point_type p2)
    { return bg::distance(p1, p2); }

#ifdef DEBUG
  template <class Tp_, class Pt_>
  friend inline std::ostream& operator<<(std::ostream& os,
      L2NN1<Tp_, Pt_> const& p)
  {
    for (auto ptit = p.begin(); ptit != p.end(); ++ptit)
      std::cout << *ptit << " ";
    std::cout << std::endl;
    return os;
  }
#endif
}; // end class L2NN1

/* L1 (Manhattan) nearest-neighbor 1 adapter.  */
template<class coordinate_type, class point_type>
class L1NN1
{
public:
  typedef L2NN1<coordinate_type, point_type> nn1_type;
  typedef typename nn1_type::value_type     value_type;
  typedef typename nn1_type::const_iterator const_iterator;
  typedef typename nn1_type::const_iterator iterator;
  typedef typename nn1_type::size_type      size_type;

private:
  // FIXME figure out how to pass custom distance predicate to boost rtree.
  // For now use L2 for this as it's [hopefully] close enough.
  nn1_type nn_;

public:

  template<class InputIter>
  L1NN1(InputIter ptbegin, InputIter ptend)
    : nn_(ptbegin, ptend)
  {}

  // Iterate over input points.
  inline const_iterator begin(void) const { return nn_.begin(); }
  inline const_iterator end(void) const { return nn_.end(); }
  inline size_type size(void) const { return nn_.size(); }
  inline void insert(const point_type &p) { nn_.insert(p); }

  /* XXX see above.  */
  inline point_type const& operator()(point_type const& pt) const {
    return nn_(pt);
  }

  inline static coordinate_type distance(point_type p1, point_type p2)
    { return l1dist<coordinate_type>(p1, p2); }

  inline static coordinate_type comparable_distance(point_type p1, point_type p2)
    { return distance(p1, p2); }

#ifdef DEBUG
  template <class Tp_, class Pt_>
  friend inline std::ostream& operator<<(std::ostream& os,
      L1NN1<Tp_, Pt_> const& p)
  {
    for (auto ptit = p.begin(); ptit != p.end(); ++ptit)
      std::cout << *ptit << " ";
    std::cout << std::endl;
    return os;
  }
#endif
}; // end class L1NN1

} // end namespace cfla
