#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cerrno>
#include <bitset>
#include <set>
#include <iomanip>
#include <array>
#include <algorithm>
#include <iterator>
#include <cmath>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#if !defined(CV_VERSION_MAJOR) || CV_VERSION_MAJOR < 3
#include <opencv2/contrib/contrib.hpp>
#include <opencv2/legacy/legacy.hpp>
#endif

#include <boost/range/join.hpp>
#include <boost/foreach.hpp>

#include "vgame.h"

#define i2f(...) static_cast<float>(__VA_ARGS__)
#define i2u(...) static_cast<unsigned int>(__VA_ARGS__)

using namespace std;
using namespace cv;

struct options_t {
  // Global configurable options
  bool drawEdges=false;
  bool dogrid=false;
  bool gridLabels=true;
  unsigned int ngridx=10u, ngridy=10u; // number of grid lines (if any)
  int gridThickness=1; // thickness of grid lines
  bool debug = false;
  ColormapTypes colormap = COLORMAP_HOT;
  bool fill_inputs = false;
  int screenWidth = 1920;
  int screenHeight = 1080;

  string sites_path;
  string users_path;
} options;

static options_t opts;

// Colors for solution polygons
static const size_t maxcolor = (1 << 8) - 1;
static const Scalar fillColor(maxcolor-1);
static const Scalar holeColor(0);
static const Scalar gridColor(fillColor);

// Default fonts
static const int FONT = FONT_HERSHEY_SCRIPT_SIMPLEX;
static const double FONT_SCALE = 0.5;
static const int FONT_THICKNESS = 1;
static const Scalar FONT_COLOR(Scalar::all(255));

static void usage(const char *prog, const char *errmsg=NULL)
#ifndef _MSC_VER
  __attribute__((noreturn))
#endif
	;

void
usage(const char *prog, const char *errmsg)
{
  cerr << "usage: " << prog << " [OPTIONS] <sites-file> <users-file>" << endl
    << "OPTIONS are:" << endl
    << "  -F			Fill input polygons (default: trace)" << endl
    << "  -C COLORMAP		Give a colormap index" << endl
    << "    0:autumn, 1:bone, 2:jet, 3:winter, 4:rainbow, 5:ocean," << endl
    << "    6:summer, 7:spring, 8:cool, 9:hsv, 10:pink, [11:hot], 12:parula"
    << endl
    << "  -X <N>		Number of X grid lines to draw (default: 10)"
    << endl
    << "  -Y <N>		Number of Y grid lines to draw (default: 10)"
    << endl
    << "  -T <N>		Thickness (pixels) of grid lines (default: 1)"
    << endl
    << "  -e			Draw edges of VD (default: no)" << endl
    << "  -E			Do not draw edges of VD (default)" << endl
    << "  -g			Draw a grid (default: no grid)" << endl
    << "  -G			Do not draw a grid (this is default)" << endl
    << "  -l			Label the grid lines (default: yes)" << endl
    << "  -L			Do not label grid lines" << endl
    << "  -d			Debug mode (extra output)" << endl
    << "  -W SIZE		Max width of the display (default: 1920)"
    << endl
    << "  -H SIZE		Max height of the display (default: 1080)"
    << endl
    ;
  if (errno) {
    cerr << endl;
    perror(NULL);
  }
  if (errmsg && errmsg[0]) {
    cerr << endl << "error: " << errmsg << endl;
  }
  exit(1);
}


template <typename T>
T getenum(int maxval, const char *instr, ostream &err, const char *errtype)
{
  istringstream in(instr);
  int ival = -1;
  in >> ival;
  if (!in.eof() || in.fail() || ival < 0 || ival >= maxval) {
    err << "bad value for " << errtype << ": '" << instr << "'";
    opterr = 1;
    return static_cast<T>(-1);
  }
  return static_cast<T>(ival);
}

static const char *sopts = "FC:X:Y:T:eEgGlLdW:H:hs:u:";

// Parses options and sets the global options structure.
static void
get_options(int argc, char *argv[], options_t &o)
{
  int opt;
  opterr = 0;
  ostringstream errstr;
  while ((opt = getopt(argc, argv, sopts)) >= 0 && !opterr) {
    switch (opt) {
    case 'd':
      o.debug = true; break;
    case 'F':
      o.fill_inputs = true; break;
    case 'C':
      o.colormap = getenum<ColormapTypes>(13, optarg, errstr, "colormap");
      break;
    case 'X':
      o.ngridx = getenum<int>(INT_MAX, optarg, errstr, "X grid lines"); break;
    case 'Y':
      o.ngridy = getenum<int>(INT_MAX, optarg, errstr, "Y grid lines"); break;
    case 'e':
      o.drawEdges = true; break;
    case 'E':
      o.drawEdges = false; break;
    case 'g':
      o.dogrid = true; break;
    case 'G':
      o.dogrid = false; break;
    case 'l':
      o.gridLabels = true; break;
    case 'L':
      o.gridLabels = false; break;
    case 'T':
      o.gridThickness = getenum<int>(INT_MAX, optarg, errstr, "grid thickness");
      break;
    case 'W':
      o.screenWidth = getenum<int>(INT_MAX, optarg, errstr, "screen width");
      break;
    case 'H':
      o.screenHeight = getenum<int>(INT_MAX, optarg, errstr, "screen height");
      break;
    case 'h':
      usage(argv[0]);
    case ':':
      errstr << "option -" << (char)opt << " missing argument";
      if (!opterr) opterr = 1;
      break;
    default:
      errstr << "unknown option -" << (char)opt;
      if (!opterr) opterr = 1;
      break;
    }
  }
  if (opterr) {
    usage(argv[0], errstr.str().c_str());
  }

  if (argc - optind < 2) {
    usage(argv[0], "not enough arguments");
  }

  // Unlimited screen size when set to zero.
  if (o.screenHeight == 0) o.screenHeight = INT_MAX;
  if (o.screenWidth == 0) o.screenWidth = INT_MAX;

  if (o.debug) {
    cout << "reading polygons from file '" << argv[optind] << "'" << endl;
  }

  o.sites_path = string(argv[optind]);
  o.users_path = string(argv[optind+1]);
}

static inline Scalar
colorKey(int i, int imax) {
  return ((i+1) * maxcolor)/ (imax+1);
}

// Draw grid lines.
static void
drawGrid(Mat img, int nx, int ny, const Scalar &color, int thicc, bool labels)
{
  const Size sz(img.size());
  const unsigned int dx = sz.width / nx;
  const unsigned int dy = sz.height / ny;
  const int xspace = 15;
  const int yspace = 15;

  // Vertical lines.
  for (unsigned int xpos=0; xpos < i2u(sz.width); xpos += dx) {
    cv::line(img, Point(xpos, 0), Point(xpos, sz.height), color, thicc);
    // Label the grid lines.
    if (labels) {
      ostringstream ss; ss << xpos;
      putText(img, ss.str(), Point(xpos + xspace, sz.height - yspace),
          FONT, FONT_SCALE, FONT_COLOR, FONT_THICKNESS);
    }
  }

  // Horizontal lines. Need to flip the labels.
  for (unsigned int ypos=0; ypos < i2u(sz.height); ypos += dy) {
    cv::line(img, Point(0, ypos), Point(sz.width, ypos), color, thicc);
    // Label the grid lines.
    if (labels) {
      ostringstream ss; ss << (sz.height - ypos - dy);
      putText(img, ss.str(), Point(xspace, ypos + dy - yspace),
          FONT, FONT_SCALE, FONT_COLOR, FONT_THICKNESS);
    }
  }
}

static Mat
translate(Mat img, float transx, float transy)
{
  Mat displace = (Mat_<float>(2,3)<< 1, 0, transx, 0, 1, transy);
  cv::warpAffine(img, img, displace, img.size());
  return img;
}

static Mat
resizeToFit(Mat img, Mat out, const Rect &bbox, const Size &screen)
{
  // First subtract origin (top-left) of the rectangle and crop.
  img = translate(img, -bbox.x, -bbox.y);
  img = img(Rect(0, 0, bbox.width, bbox.height));

  // Then see if we need to squash the image to fit our screen,
  // maintaining aspect ratio.
  Size nsz(bbox.width, bbox.height);
  if (nsz.height > screen.height) {
    nsz.width = static_cast<int>((screen.height/(float)nsz.height) * nsz.width);
    nsz.height = screen.height;
  }

  if (nsz.width > screen.width) {
    nsz.width = screen.height;
    nsz.width = static_cast<int>((screen.width/(float)nsz.width) * nsz.height);
  }

  cv::resize(img, out, nsz);
  return out;
}

// Read points from a file (one point per line) and write to out.
template <class It>
static size_t
readPoints(const string &path, It out)
{
  ifstream inf;
  inf.exceptions(ifstream::failbit | ifstream::badbit);
  inf.open(path);
  inf >> std::skipws;
  size_t npts = 0u;
  while (!inf.fail() && !inf.eof())
  {
    Point pt;
    inf >> pt.x;
    inf >> pt.y;
    if (!inf.fail())
    {
      *out++ = pt;
      ++npts;
    }
  }
  return npts;
}

template<class It>
static void
readOrDie(const string &path, It out)
{
  try {
    errno = 0;
    readPoints(path, out);
  } catch (exception) { // should be ios_base::failure; PR libstdc++/66145
    if (errno) {
      cerr << "error opening file '" << path << "': ";
      perror(NULL);
      exit(2);
    }
  }
}

struct voronoi_plane
{
private:
  const vector<Point> &sites_;
  const Size &bounds_;

  void clip_infinite_edge(const edge_t& edge, Point2d &p0, Point2d &p1) const
  {
    const cell_t& cell1 = *edge.cell();
    const cell_t& cell2 = *edge.twin()->cell();
    Point2d origin, direction;

    Point2d s1 = sites_[cell1.source_index()];
    Point2d s2 = sites_[cell2.source_index()];
    origin.x = (s1.x + s2.x) * 0.5;
    origin.y = (s1.y + s2.y) * 0.5;
    direction.x = (s1.y - s2.y);
    direction.y = (s2.x - s1.x);

    coordinate_t koef =
        bounds_.width / (std::max)(fabs(direction.x), fabs(direction.y));
    if (edge.vertex0() == NULL) {
      p0 = Point2d(
          origin.x - direction.x * koef,
          origin.y - direction.y * koef);
    } else {
      p0 = Point2d(edge.vertex0()->x(), edge.vertex0()->y());
    }
    if (edge.vertex1() == NULL) {
      p1 = Point2d(
          origin.x + direction.x * koef,
          origin.y + direction.y * koef);
    } else {
      p1 = Point2d(edge.vertex1()->x(), edge.vertex1()->y());
    }
  }

  void clip_line(const Point2d &inpt, const edge_t *clip, Point3d &line) const
  {
    // Get the endpoints of the vertex, which might be infinite.
    Point2d v0, v1;
    clip_infinite_edge(*clip, v0, v1);
    using std::abs;
    using std::numeric_limits;
    if (abs(v0.x - v1.x) < numeric_limits<double>::epsilon()
        && abs(v0.y - v1.y) < numeric_limits<double>::epsilon())
    {
      line = Point3d(0.0);
      return;
    }
    line.x = -(v0.y - v1.y);
    line.y = v0.x - v1.x;
    line.z = (v0.y - v1.y) * v0.x - (v0.x - v1.x) * v0.y;
    if (inpt.x * line.x + inpt.y * line.y + line.z > 0.0) {
      line *= -1.0;
    }
  }

public:
  voronoi_plane(const vector<Point> &sites, const Size &bounds)
    : sites_(sites), bounds_(bounds)
  {}

  Rectd boundingRect (const edge_t *const edges) const
  {
    Rectd bbox;
    bool bbox_init = false;
    const edge_t *e = edges;
    do {
      // Get the endpoints of the vertex, which might be infinite.
      Point2d v0, v1;
      clip_infinite_edge(*e, v0, v1);
      if (!bbox_init)
      {
        bp::set_points(bbox, v0, v1);
        bbox_init = true;
      }
      else {
        bp::encompass(bbox, v0);
        bp::encompass(bbox, v1);
      }
    } while ((e = e->next()) != edges);
    return bbox;
  }

  void get_clip_lines (const Point &inside, const edge_t *const edges,
      vector<Point3d> &clip_lines) const
  {
    const Point2d inpt(inside.x, inside.y);
    const edge_t *e = edges;
    do {
      Point3d cline;
      clip_line(inpt, e, cline);
      if (cline.x >= std::numeric_limits<double>::epsilon ()
          || abs(cline.y) >= std::numeric_limits<double>::epsilon ())
        clip_lines.push_back(cline);
    } while ((e = e->next()) != edges);
  }

  void draw_cell (Mat img, const edge_t *const edges) const
  {
    size_t site_idx = edges->cell()->source_index();
    const Scalar color = colorKey(site_idx, sites_.size());
    const edge_t *e = edges;
    do {
      Point2d v0, v1;
      clip_infinite_edge(*e, v0, v1);
      cv::line(img, v0, v1, color, 2);
    } while ((e = e->next()) != edges);
  }
};

static inline bool
is_inside(const Point2d &pt, const vector<Point3d> &clips)
{
  // Use clip lines to determine whether the point lies inside the polygon
  for (auto clip = clips.begin(); clip != clips.end(); ++clip)
    if (pt.x * clip->x + pt.y * clip->y + clip->z >= 0.0)
      return true;

  return false;
}

static void
voronoi_map(Mat img,
    const vector<Point> &sites, const vector<Point> &users, vector<int> &u2s,
    const Size &bounds)
{
  // u2s maps user index to nearest facility index
  voronoi_diagram vd;
  boost::polygon::construct_voronoi(sites.begin(), sites.end(), &vd);
  u2s = vector<int>(users.size(), -1);

  // Helper for clipping infinite edges etc...
  const voronoi_plane vp(sites, bounds);

  // Construct an R-tree for range queries on the user points.
  // We use this to efficiently find the user points that lie within each
  // voronoi cell using its bounding box.
  typedef pair<Point2d, int> rvalue_t;
  vector<rvalue_t> uvals;
  size_t user_idx = 0u;
  for (auto user = users.begin(); user != users.end(); ++user, ++user_idx)
    uvals.push_back(make_pair(Point2d(user->x, user->y), user_idx));
  bgi::rtree<rvalue_t, bgi::quadratic<16> > utree(uvals.begin(), uvals.end());

  for (auto it = vd.cells().begin(); it != vd.cells().end(); ++it)
  {
    const cell_t &cell = *it;
    Point2d cmin, cmax, v0, v1;
    const edge_t *edge = cell.incident_edge();
    const size_t cell_idx = cell.source_index();
    const Point &center = sites[cell_idx];

    if (opts.drawEdges) {
      vp.draw_cell(img, edge);
      if (opts.debug) {
        cv::imshow("output", img);
        cv::waitKey(0);
      }
    }

    // Get the cell bounding box for range query.
    Rectd bbox = vp.boundingRect (edge);

    // Generate clip lines.
    vector<Point3d> clip_lines;
    vp.get_clip_lines (center, edge, clip_lines);

    // Map the users that lie inside our cell to this site.
    vector<rvalue_t> query_users;
    utree.query(bgi::within(bbox), back_inserter(query_users));
    BOOST_FOREACH(rvalue_t const& v, query_users)
    {
      const Point2d &user = v.first;
      int user_idx = v.second;
      if (u2s[user_idx] == -1 && is_inside(user, clip_lines))
        u2s[user_idx] = cell_idx;
    }
  }
}

static void
draw_sites(Mat img, const vector<Point> &sites, const vector<Point> &users,
    vector<int> &u2s)
{
  // Fill sites
  for (size_t site_idx = 0u; site_idx < sites.size(); ++site_idx) {
    cv::circle(img, sites[site_idx], 10, colorKey(site_idx, sites.size()), -1);
  }

  // Don't fill users
  for (size_t user_idx = 0u; user_idx < users.size(); ++user_idx) {
    cv::circle(img, users[user_idx], 5, colorKey(u2s[user_idx], sites.size()));
  }
}

int main(int argc, char *argv[])
{
  get_options(argc, argv, opts);

  vector<Point> sites, users;

  readOrDie(opts.sites_path, back_inserter(sites));
  readOrDie(opts.users_path, back_inserter(users));

  if (sites.empty()) {
    cerr << "error: empty sites file (" << opts.sites_path << ")" << endl;
    exit(2);
  }

  // We got resolution (x, y) as (width, height); flip to (rows, cols)
  Size resolution(opts.screenWidth, opts.screenHeight);
  Mat img(resolution.height, resolution.width, CV_8U);

  if (opts.debug)
    cout << "canvas size " << img.cols << " x " << img.rows << endl;

  // Draw a grid first if we want.
  if (opts.dogrid)
    drawGrid(img, opts.ngridx, opts.ngridy, gridColor, opts.gridThickness,
        opts.gridLabels);
  // Flip vertical initially so text comes out up-right (since we flip later)
  cv::flip(img, img, 0);

  vector<int> u2s;
  voronoi_map(img, sites, users, u2s, resolution);
  draw_sites(img, sites, users, u2s);

  // Convert grayscale gradient to color gradient for plotting
  Mat colorimg;
  auto pts = boost::join(sites, users);
  Rect bbox = boundingRect(vector<Point>(pts.begin(), pts.end()));
  // squash to screen size, maintaining aspect ratio
  if ((bbox.width > 0 || bbox.height > 0)
      && (bbox.height > opts.screenHeight || bbox.width > opts.screenWidth)) {
    img = resizeToFit(img, img, bbox, resolution);
  }
  cv::flip(img, colorimg, 0); // flip vertical
  applyColorMap(colorimg, colorimg, opts.colormap);

  imshow("output", colorimg);
  waitKey(0);

  return 0;
}
