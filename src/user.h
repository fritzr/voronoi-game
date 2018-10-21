#include <vector>
#include <set>
#include <array>

#include "opencv_compat.h"

#include "polygon.h"
#include "shpReader.h"

using namespace std;
using namespace cv;

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

class User
{
private:
  /* Return the distance from one point to another.  */
  static inline double distance(Point2d p1, Point2d p2) {
    return cv::norm(p2 - p1);
  }

  /* Return the distance from a center point to the line segment a-b.  */
  static inline double distance(Point2d p, Point2d l1, Point2d l2)
  {
    Point3d aa(l1.x, l1.y, 1.0);
    Point3d bb(l2.x, l2.y, 1.0);
    Point3d l(aa.cross(bb));
    return std::abs((p.x * l.x + p.y * l.y + l.z)) * 1.0
      / std::sqrt(double(l.x * l.x + l.y * l.y));
  }

  /* Return the distance from a point to the closest point on a ring.  */
  static inline double distance(Point2d p, const c_ply &ring) {
    return cv::pointPolygonTest(ring, p, true); // measureDistance=true
  }

  /* Commutative with above.  */
  static inline double distance(const c_ply &ring, Point2d p) {
    return distance(p, ring);
  }

public:
  /* Each c_polygon must contain only one ring.  */
  User(const Point2d& pt, vector<c_polygon> plys)
    : center(pt), rings(plys)
  {
#if 0
    for (auto pply = plys.begin(); pply != plys.end(); ++pply)
    {
      /* Take only the outermost ring of the polygon.  */
      c_polygon &polygon = *pply;
      c_ply &ply = polygon.front(); /* outer ring */
      rings.emplace_back(pt, ply); /* convert to user_ring */
    }
#endif
  }

  /* To identify the travel time to a point, we first have to draw a ray from
   * our center point to the query point. Then we must intersect our rings with
   * the ray to find the upper and lower bound points along two known rings.
   * Then we can interpolate between the FTT of the two rings based on the
   * linear distance between the two rings. This approximation should be "close
   * enough" - of course it is a tradeoff. The more rings, the more accurate
   * the travel time interpolation, but the higher the computation overhead
   * is in traversing the rings.
   */
  double travelTime(Point2d pt) const;

private:
  Point2d center;
  /* Fixed-travel-time (FTT) rings.  */
  vector<c_polygon> rings;
};

static vector<User>
readUsers(const string &points_path, const string &rings_path)
{
  shp::point_map upoints = shp::readIndexedPoints(points_path);
  shp::polygon_map rings = shp::readPolygons(rings_path);
  vector<User> users;
  for (auto uit = upoints.begin(); uit != upoints.end(); ++uit)
  {
    unsigned int pointIndex = uit->first;
    Point2d center = uit->second;
    /* Match a ring list with the user.  */
    auto ring_lookup = rings.find(pointIndex);
    if (ring_lookup != rings.end())
    {
      auto ring_list = ring_lookup->second;
      users.emplace_back(center, ring_list);
    }
    else
      users.emplace_back(center, vector<c_polygon>());
  }
  return users;
}
