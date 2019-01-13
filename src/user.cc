#include "user.h"

#include "adapt_boost_poly.h"

using namespace std;

template<typename Pt_>
bool
User<Pt_>::isochrone_bounds(point_type const& query,
    const_iterator &lb, const_iterator &ub) const
{
  if (isolines_.empty())
    return false;

  /* Find the upper/lower rings between which this point is contained.  */
  lb = isolines_.end();
  ub = isolines_.begin();

  /* Keep going until we find the first ring that encloses the query point. */
  while (ub != isolines_.end() && !ub->enclosed(query))
  {
    lb = ub;
    ++ub;
  }

  /* This should already be captured by rings.empty() above. */
  if (ub == isolines_.end() && lb == isolines_.end())
    return false;

  return true;
}

template<typename Pt_>
void
User<Pt_>::isochrone(polygon_type &result, point_type const& query) const
{
  const_iterator lower_ring, upper_ring;
  if (isochrone_bounds(query, lower_ring, upper_ring))
  {
    // TODO for now we just copy one of the bounds. We should really generate
    // an interpolated isochrone from travelTime(query).
    result.copy(*((upper_ring != isolines_.end()) ? upper_ring : lower_ring));
    result.triangulate(); // we will need this
  }
}

/* Find the travel time given fixed-travel-time (FTT) rings.  */
template<typename Pt_>
typename User<Pt_>::coordinate_type
User<Pt_>::travelTime(typename User<Pt_>::point_type query) const
{
  /* If the helper fails we have no way to know the answer.  */
  const_iterator lower_ring, upper_ring;
  if (!isochrone_bounds(query, lower_ring, upper_ring))
    return std::numeric_limits<coordinate_type>::max();

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
    du = -bg::distance(lower_ring->front(), query);
    dl = bg::distance(center(), lower_ring->front()) - du;
    fttl = 0.0;
    fttu = lower_ring->front().extra().ftt;
  }

  /* Time-to-travel.  */
  return fttl + (abs(dl / (du + dl)) * abs(fttu - fttl));
}

template class User<boost::geometry::model::d2::point_xy<double> >;
