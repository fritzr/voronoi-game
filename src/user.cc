#include "user.h"

#include "opencv_compat.h"
#include "adapt_boost_poly.h"

using namespace cv;
using namespace std;

template<typename Pt_>
typename User<Pt_>::polygon_type
User<Pt_>::isochrome(coordinate_type time) const
{
  // TODO
  return polygon_type();
}

/* Find the travel time given fixed-travel-time (FTT) rings.  */
template<typename Pt_>
typename User<Pt_>::coordinate_type
User<Pt_>::travelTime(typename User<Pt_>::point_type query) const
{
  if (isolines_.empty())
    return std::numeric_limits<coordinate_type>::infinity();

  /* Find the upper/lower rings between which this point is contained.  */
  auto lower_ring = isolines_.end();
  auto upper_ring = isolines_.begin();

  /* Keep going until we find the first ring that encloses the query point. */
  while (upper_ring != isolines_.end() && !upper_ring->enclosed(query))
  {
    lower_ring = upper_ring;
    ++upper_ring;
  }

  /* This should already be captured by rings.empty() above. */
  if (upper_ring == isolines_.end() && lower_ring == isolines_.end())
    return HUGE_VAL;

  /* Interpolate the travel time using the FTT endpoints and distance.
   * We interpolate with:
   *
   *  L <- dl ->   <- du -> U
   *  L -------- F -------- U
   *  L                     U
   * fttl       ttt        fttu
   *
   * Where ttt = | dl / (du + dl) | * | fttu - fttl | .
   */
  coordinate_type dl = 0.0, du = 0.0, fttl = 0.0,
                  fttu = std::numeric_limits<coordinate_type>::infinity();

  /* 3 cases for the surrounding rings L, U:
   *    1. L=nul, F < U
   *    2.            L < F < U
   *    3.                    L < F, U=nul
   */

  /* Case 1. No lower ring -- lower distance is from center point (origin).  */
  if (lower_ring == isolines_.end() /* && upper_ring != rings.end() */)
  {
    dl = bg::distance(center(), query);
    fttl = 0.0;
    du = bg::distance(query, upper_ring->front());
    fttu = upper_ring->front().extra().ftt;
  }

  /* Case 2. Interpolate between the bounding rings.  */
  if (lower_ring != isolines_.end() && upper_ring != isolines_.end())
  {
    dl = bg::distance(lower_ring->front(), query);
    fttl = lower_ring->front().extra().ftt;
    du = bg::distance(query, upper_ring->front());
    fttu = upper_ring->front().extra().ftt;
  }

  /* Case 3. No upper bound -- extrapolate past lower bound.  */
  if (lower_ring != isolines_.end() && upper_ring == isolines_.end())
  {
    coordinate_type du = -bg::distance(lower_ring->front(), query);
    dl = bg::distance(center(), lower_ring->front()) - du;
    fttl = 0.0;
    fttu = lower_ring->front().extra().ftt;
  }

  /* Time-to-travel.  */
  return fttl + (abs(dl / (du + dl)) * abs(fttu - fttl));
}

template class User<cv::Point2d>;
