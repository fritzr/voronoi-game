#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cerrno>
#include <array>
#include <algorithm>
#include <iterator>
#include <cmath>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

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

// With older OpenCV, you can't construct Point2d from Point
#if defined(CV_VERSION_EPOCH) && CV_VERSION_MINOR < 13
typedef pointiter<typename std::vector<cv::Point>::iterator, cv::Point2d> p2di;

// Point order from RotatedRect::points()
#define RRBL 1
#define RRTL 0
#define RRTR 3
#define RRBR 2

#else
typedef typename std::vector<cv::Point>::iterator p2di;

// Point order from RotatedRect::points()
#define RRBL 0
#define RRTL 1
#define RRTR 2
#define RRBR 3
#endif

#include <boost/range/join.hpp>
#include <boost/foreach.hpp>

#include "main.h"
#include "voronoi.h"
#include "components.h"

#define i2f(...) static_cast<float>(__VA_ARGS__)
#define i2u(...) static_cast<unsigned int>(__VA_ARGS__)

using namespace std;
using namespace cv;
using namespace voronoi;

enum DrawRects {
  RECTS_NONE = 0,
  RECTS_ROTATED = 1,
  RECTS_UPRIGHT = 2,
  RECTS_BOTH = 3,
  RECTS_BAD = 4,
};

enum DrawCell {
  CELL_NONE = 0,
  CELL_ROTATED = 1,
  CELL_UPRIGHT = 2,
  CELL_BOTH = 3,
  CELL_BAD = 4,
};

struct options_t
{
  // Global configurable options
  bool drawEdges=true;
  bool dogrid=false;
  bool gridLabels=true;
  unsigned int ngridx=10u, ngridy=10u; // number of grid lines (if any)
  int gridThickness=1; // thickness of grid lines
  bool debug = false;
  ColormapTypes colormap = COLORMAP_HOT;
  bool fill_inputs = false;
  int screenWidth = 1920;
  int screenHeight = 1080;
  DrawCell computeCell = CELL_NONE;
  int drawRects = true;
  bool dumpGraph = false; // adjacency graph
  VoronoiDiagram<double>::SearchMethod queryType
    =VoronoiDiagram<double>::Default;

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
const int FONT_XSPACE = 10;
const int FONT_YSPACE = 10;

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
    << "  -e			Draw edges of VD (default: yes)" << endl
    << "  -E			Do not draw edges of VD (default: yes)" << endl
    << "  -g			Draw a grid (default: no grid)" << endl
    << "  -G			Do not draw a grid (this is default)" << endl
    << "  -l			Label the grid lines (default: yes)" << endl
    << "  -L			Do not label grid lines" << endl
    << "  -d			Debug mode (extra output)" << endl
    << "  -W SIZE		Max width of the display (default: 1920)"
    << endl
    << "  -H SIZE		Max height of the display (default: 1080)"
    << endl
    << "  -q TYPE		Query type for mapping users to sites:" << endl
    << "    [0:default], 1:brute force, 2:range search, 3:NN(1) (fastest)"
    << endl
    << "  -c			Compute connected components." << endl
    << "  -r			Draw rects (default: none):" << endl
    << "    [0:none], 1:rotated, 2:upright, 3:both" << endl
    << "  -j			Dump adJacency list (default: no)" << endl
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
static T
getenum(int maxval, const char *instr, ostream &err, const char *errtype)
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

static const char *sopts = "Fc:C:X:Y:T:eEgGlLdW:H:hs:u:q:r:";

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
    case 'q':
      o.queryType = getenum<VoronoiDiagram<double>::SearchMethod>(
        VoronoiDiagram<double>::SM_BAD, optarg, errstr, "query type");
      break;
    case 'r':
      o.drawRects = getenum<DrawRects>(RECTS_BAD, optarg, errstr, "rect type");
      break;
    case 'c':
      o.computeCell = getenum<DrawCell>(CELL_BAD, optarg, errstr, "cell type");
      break;
    case 'j':
      o.dumpGraph = true; break;
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

// Draw grid lines.
static void
drawGrid(Mat img, int nx, int ny, const Scalar &color, int thicc, bool labels)
{
  const Size sz(img.size());
  const unsigned int dx = sz.width / nx;
  const unsigned int dy = sz.height / ny;

  // Vertical lines.
  for (unsigned int xpos=0; xpos < i2u(sz.width); xpos += dx) {
    cv::line(img, Point(xpos, 0), Point(xpos, sz.height), color, thicc);
    // Label the grid lines.
    if (labels) {
      ostringstream ss; ss << xpos;
      putText(img, ss.str(), Point(xpos + FONT_XSPACE, sz.height - FONT_YSPACE),
          FONT, FONT_SCALE, FONT_COLOR, FONT_THICKNESS);
    }
  }

  // Horizontal lines. Need to flip the labels.
  for (unsigned int ypos=0; ypos < i2u(sz.height); ypos += dy) {
    cv::line(img, Point(0, ypos), Point(sz.width, ypos), color, thicc);
    // Label the grid lines.
    if (labels) {
      ostringstream ss; ss << (sz.height - ypos - dy);
      putText(img, ss.str(), Point(FONT_XSPACE, ypos + dy - FONT_YSPACE),
          FONT, FONT_SCALE, FONT_COLOR, FONT_THICKNESS);
    }
  }
}

#if 0
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
#endif

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
  } catch (exception&) { // should be ios_base::failure; PR libstdc++/66145
    if (errno) {
      cerr << "error opening file '" << path << "': ";
      perror(NULL);
      exit(2);
    }
  }
}

template<typename T>
static inline T deg2rad(T deg) {
  return (deg * static_cast<T>(M_PI)) / static_cast<T>(180.0);
}
template<typename T>
static inline T rad2deg(T rad) {
  return (rad * static_cast<T>(180.0)) / static_cast<T>(M_PI);
}

// Pre-computed 45-deg Euler rotation matrices
static const float ANGLE_DEG = 45.0f;
static const float pangle = deg2rad(ANGLE_DEG);
static const float nangle = -pangle;
static const Mat rotZn = Mat(Matx22f({
  cosf(nangle), -sinf(nangle), /* 0, */
  sinf(nangle),  cosf(nangle), /* 0, */
  /*         0,             0,    1, */
}));
static const Mat rotZp = Mat(Matx22f({
  cosf(pangle), -sinf(pangle), /* 0, */
  sinf(pangle),  cosf(pangle), /* 0, */
  /*         0,             0,    1, */
}));

static inline Point2f rotateZ2f_neg(Point2f pt) {
  Mat ans = Mat(Matx12f({pt.x, pt.y})) * rotZn;
  return Point2f(ans.at<float>(0u, 0u), ans.at<float>(0u, 1u));
}

static inline Point2f rotateZ2f_pos(Point2f pt) {
  Mat ans = Mat(Matx12f({pt.x, pt.y})) * rotZp;
  return Point2f(ans.at<float>(0u, 0u), ans.at<float>(0u, 1u));
}

static Point2f pbias, nbias;

template<class pt>
static inline pt rotp(pt p)
{
  return rotateZ2f_pos(p) - pbias;
}

template<class pt>
static inline pt rotn(pt p)
{
  return rotateZ2f_neg(p) - nbias;
}

template<class pt>
static inline pt check_rot(pt p)
{
  return (opts.drawRects == RECTS_UPRIGHT) ? rotp(p) : p;
}

static inline Scalar
colorKey(int i, int imax) {
  return Scalar(((i+1) * double(maxcolor))/ (imax+1));
}

template <class Tp_>
static void
draw_cells (Mat img, const VoronoiDiagram<Tp_> &vd)
{
  size_t site_idx = 0u;
  for (auto cell_iter=vd.begin(); cell_iter!=vd.end(); ++cell_iter, ++site_idx)
  {
    Scalar color = colorKey(site_idx, vd.sites_size());
    auto const& cell = *cell_iter;
    for (auto edge_iter = vd.cell_begin(cell); edge_iter != vd.cell_end(cell);
        ++edge_iter) {
      auto edge = *edge_iter;
      cv::line(img, check_rot(edge.p0), check_rot(edge.p1), color, 2);
    }
    // With debug, show each cell's borders as it is drawn
    if (opts.debug) {
      cv::imshow("output", img);
      cv::waitKey(0);
    }
  }
}

template <class Tp_>
static void
draw_sites(Mat img, const VoronoiDiagram<Tp_> &vd)
{
  // Fill sites
  for (size_t site_idx = 0u; site_idx < vd.sites_size(); ++site_idx) {
    Scalar color = colorKey(site_idx, vd.sites_size());
    cv::circle(img, check_rot(vd.site(site_idx)), 10, color, -1);
  }

  // Don't fill users
  for (size_t user_idx = 0u; user_idx < vd.users_size(); ++user_idx) {
    Scalar color = colorKey(vd.site_index(user_idx), vd.sites_size());
    cv::circle(img, check_rot(vd.user(user_idx)), 5, color);
  }
}

// rotated rect bounding box as rectangle data
static bp::rectangle_data<float>
getRRBBox(RotatedRect const& rr)
{
  // Older OpenCV doesn't have floating-point bounding rects
#if defined(CV_VERSION_EPOCH) && CV_VERSION_MINOR < 13
  Rect brect = rr.boundingRect();
#else
  Rect_<float> brect = rr.boundingRect2f();
#endif
  bp::rectangle_data<float> ret;
  bp::set_points(ret, brect.tl(), brect.br());
  return ret;
}

// rotated rect bounding box as rectangle data
template<class RRectIter>
static bp::rectangle_data<float>
getRRBBox(RRectIter begin, RRectIter end)
{
  auto bbox = getRRBBox(*begin++); // assumes non-empty iterator
  while (begin != end)
  {
    RotatedRect const& rr = *begin++;
    std::array<Point2f, 4> pts;
    rr.points(pts.begin());
    for (unsigned int i = 0u; i < 4; ++i)
      bp::encompass(bbox, pts[i]);
  }
  return bbox;
}

static void
draw_rotrect(Mat img, RotatedRect const& rrect, Scalar const& color,
    int thicc=2)
{
  std::array<Point2f, 4> pts;
  rrect.points(pts.begin());
  if (thicc < 0) {
    std::array<Point, 4> ipts;
    std::copy(pts.begin(), pts.end(), ipts.begin());
    cv::fillConvexPoly(img, ipts.begin(), ipts.size(), color);
  }
  else
    for (unsigned int i = 0u; i < pts.size(); ++i)
      cv::line(img, pts[i], pts[(i+1)%pts.size()], color, thicc);
}

template<class Tp_, class RectIter>
static void
draw_rotrects(Mat img, VoronoiDiagram<Tp_> const& vd, RectIter rrectp,
    int thicc=2)
{
  unsigned int rectidx = 0u;
  for (auto userp = vd.users_begin(); userp != vd.users_end(); ++userp)
  {
    Scalar color = colorKey(vd.site_index(userp), vd.sites_size());
    const RotatedRect& rrect = *rrectp++;
    draw_rotrect(img, rrect, color, thicc);
    putText(img, to_string(rectidx++),
        Point(rrect.center.x + FONT_XSPACE, rrect.center.y - FONT_YSPACE),
        FONT, FONT_SCALE, color, FONT_THICKNESS);
  }
}

template<class rect_type, class RectIter>
static void
build_rects(Mat img, vector<Point>& sites, vector<Point>& users,
    Size const& resolution, RectIter rects_out)
{
  VoronoiDiagram<double> vd(
      p2di(sites.begin()), p2di(sites.end()),
      p2di(users.begin()), p2di(users.end()),
      resolution.width, resolution.height);
  vd.build(opts.queryType);

  // First build the L1 rotated rects.
  vector<RotatedRect> rotrects;
  vd.build_rects(back_inserter(rotrects));

  // Draw rotated rects if we want.
  if (opts.drawRects & RECTS_ROTATED)
    draw_rotrects(img, vd, rotrects.begin());

  // Get the center of the rotated rect-space.
  auto bbox = getRRBBox(rotrects.begin(), rotrects.end());
  Point2f center;
  bp::center(center, bbox);
  //cv::circle(img, center, 10, fillColor);
  // set the global bias position so points appear back near the center
  pbias = rotateZ2f_pos(center) - center;
  nbias = rotateZ2f_neg(center) - center;

  // Now rotate the entire space to output axis-up rectangles,
  // biased by the center point.
  size_t user_idx = 0u;
  for (auto rrectp = rotrects.begin(); rrectp != rotrects.end(); ++rrectp)
  {
    const RotatedRect& rrect = *rrectp;
    array<Point2f, 4> pts;
    rrect.points(pts.begin());
    for (auto pp = pts.begin(); pp != pts.end(); ++pp)
      *pp = rotp(*pp); // rotate and bias

    // Write the axis-up rectangles to the output array.
     // BL, TL, TR, BR -> xl, yl, xh, yh
    *rects_out++ = bp::construct<rect_type>(
        pts[RRBL].x, pts[RRBL].y, pts[RRTR].x, pts[RRTR].y);

    // Draw the axis-up rects if we want.
    if (opts.drawRects & RECTS_UPRIGHT)
    {
      Scalar color = colorKey(vd.site_index(user_idx), vd.sites_size());
      cv::rectangle(img, Rect(pts[RRTL], pts[RRBR]), color, 2);
    }

    ++user_idx;
  }

  // Draw cells/sites now that we know the proper bias.
  if (opts.drawEdges) {
    draw_cells(img, vd);
  }
  draw_sites(img, vd);
}

static Mat
finish_img(Mat img, Rect const& bbox, Size const& resolution)
{
  // Convert grayscale gradient to color gradient for plotting
  Mat colorimg;
  // squash to screen size, maintaining aspect ratio
  if ((bbox.width > 0 || bbox.height > 0)
      && (bbox.height > opts.screenHeight || bbox.width > opts.screenWidth)) {
    //img = resizeToFit(img, img, bbox, resolution);
  }
  cv::flip(img, colorimg, 0); // flip vertical
  applyColorMap(colorimg, colorimg, opts.colormap);
  return colorimg;
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

  typedef typename components::ConnectedComponents<double> CComp;
  typedef typename CComp::coordinate_type coordinate_type;
  typedef typename bp::rectangle_data<coordinate_type> rect_type;
  vector<rect_type> rects;
  build_rects<rect_type>(img, sites, users, resolution, back_inserter(rects));

  if (opts.computeCell)
  {
    CComp comp(rects.begin(), rects.end());
    comp.compute();
    cout << "maximal depth: " << comp.depth() << endl;

    cout << "maximal intersections from: ";
    for (auto it = comp.begin(); it != comp.end(); ++it)
      cout << *it << ", ";
    cout << endl;

    auto maxrect = comp.max();
    cout << "maximal rect:  " << maxrect << endl;
    Point tl, br;
    bp::assign(tl, bp::ul(maxrect));
    bp::assign(br, bp::lr(maxrect));
    if (opts.computeCell & CELL_UPRIGHT)
    {
      Rect r(tl, br);
      cv::rectangle(img, r, fillColor, CV_FILLED);
    }
    if (opts.computeCell & CELL_ROTATED)
    {
      Point2f center((tl.x + br.x)/2.0f, (tl.y + br.y)/2.0f);
      center = rotn(center);
      Size2f size(br.x - tl.x, tl.y - br.y);
      RotatedRect rrect(center, size, ANGLE_DEG);
      draw_rotrect(img, rrect, fillColor, CV_FILLED);
    }

    if (opts.dumpGraph)
    {
      /* Dump adjacency list */
      auto graph = comp.adj_graph();
      //size_t i = 0u, nedges = boost::num_edges(graph);
      auto itpair = boost::edges(graph);
      auto eit = itpair.first;
      auto end = itpair.second;
      cout << "connected rects: " << endl;
      while (eit != end)
      {
        cout << "  " << comp.index(boost::source(*eit, graph))
          << " <-> " << comp.index(boost::target(*eit, graph)) << endl;
        ++eit;
      }
    }
  }

  auto pts = boost::join(sites, users);
  Rect bbox = boundingRect(vector<Point>(pts.begin(), pts.end()));
  img = finish_img(img, bbox, resolution);

  // Now that we've flipped the image draw the user point labels
  // so they are upright.
  unsigned int user_idx = 0u;
  for (auto uit = users.begin(); uit != users.end(); ++uit)
  {
    Point up = check_rot(*uit);
    up.y = resolution.height - up.y;
    putText(img, to_string(user_idx++),
        Point(up.x + FONT_XSPACE, up.y - FONT_YSPACE),
        FONT, FONT_SCALE, FONT_COLOR, FONT_THICKNESS);
  }

  imshow("output", img);
  waitKey(0);

  return 0;
}
