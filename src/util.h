#pragma once

#ifdef MSVC
#define __attribute__()
#endif

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
#include <ctime>
#include <type_traits>
#ifdef DEBUG
#include <iostream>
#endif

typedef std::default_random_engine rng_type;
extern rng_type* rng;

template<class Type>
inline Type randrange(Type min, Type max)
{
  if (rng == nullptr)
  {
    auto seed = time(NULL);
    rng = new rng_type(seed);
#ifdef DEBUG
    std::cerr << "using RNG seed: " << seed << std::endl;
#endif
  }

  typedef typename std::conditional<std::is_integral<Type>::value,
      std::uniform_int_distribution<Type>,
      std::uniform_real_distribution<Type>
    >::type distribution_type;

  distribution_type distribution(min, max);
  return distribution(*rng);
}
