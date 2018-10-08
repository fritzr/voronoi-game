/* Compatibility for different versions of OpenCV.  */

#include <opencv2/core/core.hpp>

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#if CV_VERSION_MAJOR < 3 && !defined(OPENCV2_4_13)
#include <opencv2/contrib/contrib.hpp>
typedef int ColormapTypes;
#endif

#if defined(CV_VERSION_EPOCH)
#include <opencv2/contrib/contrib.hpp>
#include <opencv2/legacy/legacy.hpp>
#if CV_VERSION_MINOR < 13
typedef int ColormapTypes;

template <typename It, typename Pt>
struct pointiter : public std::iterator<
    typename std::iterator_traits<It>::iterator_category,
    Pt>
{
  typedef typename std::iterator_traits<It> traits;
  typedef std::iterator<
    typename traits::iterator_category,
    Pt> iterator;

  typedef typename iterator::value_type value_type;
  typedef typename iterator::reference reference;
  typedef typename iterator::pointer pointer;
  typedef typename iterator::difference_type difference_type;
  typedef size_t size_type;

  It it_;
  pointiter(It it) : it_(it) {}

  inline pointiter& operator++() { ++it_; return *this; }
  inline pointiter operator++(int) { return pointiter(it_++); }
  inline bool operator==(pointiter o) const { return it_ == o.it_; }
  inline bool operator!=(pointiter o) const { return it_ != o.it_; }
  inline bool operator<(pointiter o) const { return it_ < o.it_; }
  inline pointiter& operator--() { --it_; return *this; }
  inline pointiter operator--(int) { return pointiter(it_--); }
  inline pointiter operator+(int i) { return pointiter(it_ + i); }
  inline difference_type operator-(pointiter o) { return it_ - o.it_; }
  inline pointiter& operator+=(int i) { it_ += i; return *this; }
  inline pointiter& operator-=(int i) { it_ -= i; return *this; }
  inline value_type operator*() const {
    auto pt = *it_;
    return Pt(pt.x, pt.y);
  }
};
#endif
#endif

#if defined(CV_VERSION_EPOCH) && CV_VERSION_MINOR < 13
// With older OpenCV, you can't construct Point2d from Point
typedef pointiter<typename std::vector<cv::Point>::iterator, cv::Point2d> p2di;
#else
typedef typename std::vector<cv::Point>::iterator p2di;
#endif

#if defined(CV_VERSION_EPOCH) && CV_VERSION_MINOR < 13

// Point order from RotatedRect::points()
#define RRBL 1
#define RRTL 0
#define RRTR 3
#define RRBR 2

#else

// Point order from RotatedRect::points()
#define RRBL 0
#define RRTL 1
#define RRTR 2
#define RRBR 3
#endif

using namespace cv;
