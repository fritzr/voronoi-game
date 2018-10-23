#include "user.h"

#include "opencv_compat.h"

using namespace cv;
using namespace std;

/* Find the travel time given fixed-travel-time (FTT) rings.  */
double
User::travelTime(Point2d query) const
{
  if (rings.empty())
    return HUGE_VAL;

  /* Find the upper/lower rings between which this point is contained.  */
  auto lower_ring = rings.end();
  auto upper_ring = rings.begin();

  /* Keep going until we find the first ring that encloses the query point. */
  while (upper_ring != rings.end() && !upper_ring->enclosed(query))
  {
    lower_ring = upper_ring;
    ++upper_ring;
  }

  /* This should already be captured by rings.empty() above. */
  if (upper_ring == rings.end() && lower_ring == rings.end())
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
  double dl = 0.0, du = 0.0;
  double fttu = HUGE_VAL, fttl = 0.0;

  /* 3 cases for the surrounding rings L, U:
   *    1. L=nul, F < U
   *    2.            L < F < U
   *    3.                    L < F, U=nul
   */

  /* Case 1. No lower ring -- lower distance is from center point (origin).  */
  if (lower_ring == rings.end() /* && upper_ring != rings.end() */)
  {
    dl = distance(center, query);
    fttl = 0.0;
    du = distance(query, upper_ring->front());
    fttu = upper_ring->front().extra().ftt;
  }

  /* Case 2. Interpolate between the bounding rings.  */
  if (lower_ring != rings.end() && upper_ring != rings.end())
  {
    dl = distance(lower_ring->front(), query);
    fttl = lower_ring->front().extra().ftt;
    du = distance(query, upper_ring->front());
    fttu = upper_ring->front().extra().ftt;
  }

  /* Case 3. No upper bound -- extrapolate past lower bound.  */
  if (lower_ring != rings.end() && upper_ring == rings.end())
  {
    dl = distance(center, query);
    fttl = 0.0;
    du = distance(query, lower_ring->front());
    fttu = lower_ring->front().extra().ftt;
  }

  /* Time-to-travel.  */
  return abs(dl / (du + dl)) * abs(fttu - fttl);
}
