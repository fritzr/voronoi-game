#pragma once

template<typename T>
inline T deg2rad(T deg) {
  return (deg * static_cast<T>(M_PI)) / static_cast<T>(180.0);
}
template<typename T>
inline T rad2deg(T rad) {
  return (rad * static_cast<T>(180.0)) / static_cast<T>(M_PI);
}

// Pre-computed 45-deg Euler rotation matrices for Z-rotation
const float ANGLE_DEG = 45.0f;
const float pangle = deg2rad(ANGLE_DEG);
const float nangle = -pangle;

#include <opencv2/core/core.hpp>

inline cv::Point2f rotateZ2f_neg(cv::Point2f pt) {
  static const cv::Mat rotZn = cv::Mat(cv::Matx22f({
    cosf(nangle), -sinf(nangle), /* 0, */
    sinf(nangle),  cosf(nangle), /* 0, */
    /*         0,             0,    1, */
  }));
  cv::Mat ans = cv::Mat(cv::Matx12f({pt.x, pt.y})) * rotZn;
  return cv::Point2f(ans.at<float>(0u, 0u), ans.at<float>(0u, 1u));
}

inline cv::Point2f rotateZ2f_pos(cv::Point2f pt) {
  static const cv::Mat rotZp = cv::Mat(cv::Matx22f({
    cosf(pangle), -sinf(pangle), /* 0, */
    sinf(pangle),  cosf(pangle), /* 0, */
    /*         0,             0,    1, */
  }));
  cv::Mat ans = cv::Mat(cv::Matx12f({pt.x, pt.y})) * rotZp;
  return cv::Point2f(ans.at<float>(0u, 0u), ans.at<float>(0u, 1u));
}

#include <random>

template<class Int>
inline Int randint(Int min, Int max)
{
  static std::default_random_engine rand;
  std::uniform_int_distribution<Int> distribution(min, max);
  return distribution(rand);
}

#include <map>
#include <list>
#include <boost/geometry/geometry.hpp>
namespace bgi = boost::geometry::index;

template<class UserList, class FacilityList>
class NN1
{
public:
  typedef typename UserList::const_iterator user_citerator;
  typedef typename UserList::iterator user_iterator;

  typedef typename FacilityList::value_type point_type;
  typedef std::pair<point_type, int> fpair_type;
  typedef typename bgi::rtree<fpair_type, bgi::quadratic<16> > voronoi_tree;

private:
  UserList const& users_;
  FacilityList const& facilities_;
  std::list<fpair_type> fpairs_;
  std::map<int, int> u2s_;

public:
  NN1(UserList const& u, FacilityList const& f)
    : users_(u), facilities_(f), u2s_(users_.size(), -1)
  {
    // Construct an R-tree for NN(1) queries on the facilitys.
    // We use this to efficiently find the facility that is nearest to each user
    int facility_idx = 0;
    auto fbegin = f.begin();
    while (fbegin != f.end())
      fpairs_.emplace_back(*fbegin++, facility_idx++);
  }

  void add_site(point_type p)
  {
    facilities_.push_back(p);
    fpairs_.emplace_back(p, fpairs_.size());
  }

  template<class InputIter>
  void build(void)
  {
    // The nearest-neighbor is of course the voronoi cell.
    voronoi_tree vtree(fpairs_.begin(), fpairs_.end());
    int user_idx = 0;
    auto userp = users_.begin();
    while (userp != users_.end())
      u2s_[user_idx++] = vtree.qbegin(bgi::nearest(*userp++, 1))->second;
  }

  inline point_type nearest_facility(int up) const
    { return facilities_.at(u2s_[up]); }
  inline int nearest_index(int up) const
    { return u2s_[up]; }
  inline point_type nearest_facility(user_iterator up) const
    { return facilities_.at(u2s_[up - users_.begin()]); }
  inline int nearest_index(user_iterator up) const
    { return u2s_[up - users_.begin()]; }

};
