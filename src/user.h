#pragma once

#include <vector>
#include <string>
// #include <set>

#ifdef DEBUG
#include <iostream>
#endif

#include "opencv_compat.h"

#include "polygon.h"
#include "shpReader.h"

/* A user has a center point and several layers of rings (like a topological
 * map). Each ring approximates a known fixed travel time (FTT) radius.
 * To identify the "distance" from a user, we interpolate between these
 * rings. We use time rather than distance, because it is a more realistic
 * measurement of a user's willingness to travel, which is a better heuristic
 * for determining the ideal facility location.
 *   x------------x
 *   |  x-------x |  For example; in the simplified diagram to the left, three
 *   |  |  x--x | |  rings are shown around a center point. Each ring is
 *   |  |1 x--x | |  associated with a FTT in seconds, shown as 1, 2, or 3min.
 *   |2 x-------x |  Along all points on the outer ring it would take the user
 * 3 x------------x
 */

#if 0
/* */
template<typename point>
struct angle_sort
{
  point m_dreference;

public:
  angle_sort(const point reference)
    : m_dreference(reference) {}

  bool operator()(const point *a, const point *b) const
  {
    const point da = *a, db = *b;
    const double detb = m_dreference.cross(db);

    // nothing is less than zero degrees
    if (detb == 0 && db.x * m_dreference.x + db.y * m_dreference.y >= 0) return false;

    const double deta = m_dreference.cross(da);

    // zero degrees is less than anything else
    if (deta == 0 && da.x * m_dreference.x + da.y * m_dreference.y >= 0) return true;

    if (deta * detb >= 0) {
      // both on same side of reference, compare to each other
      return da.cross(db) > 0;
    }

    // vectors "less than" zero degrees are actually large, near 2 pi
    return deta > 0;
  }
};

class user_ring : public c_polygon
{
private:
  typedef std::set<ply_vertex *, angle_sort<Point2d> > angular_vertex_set;

public:
  /* Take the first polygon from cpl (move semantics).  */
  template<typename It>
    user_ring(const Point2d &center, c_polygon &cpl)
    {
      c_ply &cp = cpl.front();
      cpl.pop_front();
      push_back(cp);
      indexing();
      sorted.insert(cp.begin(), cp.end());
    }

private:
  angular_vertex_set sorted;
};
#endif

template<typename Pt_>
class User
{
public:
  typedef Pt_ point_type;
  typedef c_ply<point_type> ring_type;
  typedef c_polygon<point_type> polygon_type;
  typedef typename polygon_type::coordinate_type coordinate_type;

private:
  /* Return the distance from one point to another.  */
  static inline coordinate_type distance(point_type p1, point_type p2) {
    return cv::norm(p2 - p1);
  }

  /* Return the distance from a center point to the line segment a-b.  */
  static inline coordinate_type
    distance(point_type p, point_type l1, point_type l2)
  {
    cv::Point3d aa(l1.x, l1.y, 1.0);
    cv::Point3d bb(l2.x, l2.y, 1.0);
    cv::Point3d l(aa.cross(bb));
    return coordinate_type(
        std::abs((p.x * l.x + p.y * l.y + l.z)) * 1.0
        / std::sqrt(double(l.x * l.x + l.y * l.y)));
  }

  /* Return the distance from a point to the closest point on a ring.  */
  static inline coordinate_type distance(point_type p, const ring_type &ring) {
    return cv::pointPolygonTest(ring.mat(), p, true); // measureDistance=true
  }

  /* Commutative with above.  */
  static inline coordinate_type distance(const ring_type &ring, point_type p) {
    return distance(p, ring);
  }

public:
  /* Each c_polygon must contain only one ring for now.  */
  User(const point_type& pt, std::vector<polygon_type> plys)
    : center(pt), isolines(plys)
  {
    // Go ahead and triangulate now -- this will be needed for computing
    // travel time later.
    for (auto rit = isolines.begin(); rit != isolines.end(); ++rit)
    {
      rit->triangulate();
      rit->getSize(); // index
    }
  }

  /* Interpolate the travel time to a point based on the nearest isoline(s)
   * with a known fixed travel time.  */
  coordinate_type travelTime(point_type pt) const;

#ifdef DEBUG
  template<typename Pt__>
  friend inline std::ostream& operator<<(std::ostream& os, const User<Pt__> &u)
  {
    os << u.center << " [ ftts: ";
    auto cpi = u.isolines.cbegin();
    while (cpi != u.isolines.cend())
    {
      if (cpi != u.isolines.cbegin())
        os << ", ";
      os << cpi->front().extra().ftt;
      ++cpi;
    }
    os << " ]";
    return os;
  }
#endif // DEBUG

private:
  point_type center;
  /* Fixed-travel-time (FTT) isolines.  */
  std::vector<polygon_type> isolines;
};

template<typename Pt_>
extern std::vector<User<Pt_> >
readUsers(const std::string &points_path, const std::string &rings_path)
{
  auto /* point_map */ upoints = shp::readIndexedPoints<Pt_>(points_path);
  auto /* polygon_map */ rings = shp::readPolygons<Pt_>(rings_path);
  std::vector<User<Pt_> > users;
  for (auto uit = upoints.begin(); uit != upoints.end(); ++uit)
  {
    unsigned int pointIndex = uit->first;
    typename User<Pt_>::point_type center = uit->second;
    /* Match a ring list with the user.  */
    auto ring_lookup = rings.find(pointIndex);
    if (ring_lookup != rings.end())
    {
      auto ring_list = ring_lookup->second;
      users.emplace_back(center, ring_list);
    }
    else
      users.emplace_back(center,
          std::vector<typename User<Pt_>::polygon_type>());
  }
  return users;
}

extern template class User<cv::Point2d>;
