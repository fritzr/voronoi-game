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

template<class point_type, class user_container>
class NN1
{
public:

private:
  value_list site_indexes;
  user_container const& users;

  NN1(user_container const& u) : users(u) { }

  class NN1Bind {
    typedef std::pair<point_type, int> rvalue_t;
    typedef std::vector<rvalue_t> value_list;
    typedef typename bgi::rtree<rvalue_t, bgi::quadratic<16> > nn1tree;

    template<class facility_container>
    NN1Bind(facility_container const& facs) {
    }
  };

  template<class facility_container>
  NN1Bind bind(facility_container const& facs) {
    return NN1Bind(facs);
  }
};
