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

#if !defined(CV_VERSION_MAJOR) || CV_VERSION_MAJOR < 3
#include <opencv2/contrib/contrib.hpp>
#include <opencv2/legacy/legacy.hpp>
#endif

#include <boost/range/join.hpp>
#include <boost/foreach.hpp>

#include "main.h"
#include "voronoi.h"

#define i2f(...) static_cast<float>(__VA_ARGS__)
#define i2u(...) static_cast<unsigned int>(__VA_ARGS__)

using namespace std;
using namespace cv;

struct options_t {
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
  } catch (exception&) { // should be ios_base::failure; PR libstdc++/66145
    if (errno) {
      cerr << "error opening file '" << path << "': ";
      perror(NULL);
      exit(2);
    }
  }
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
      cv::line(img, edge.p0, edge.p1, color, 2);
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
    cv::circle(img, vd.site(site_idx), 10, color, -1);
  }

  // Don't fill users
  for (size_t user_idx = 0u; user_idx < vd.users_size(); ++user_idx) {
    Scalar color = colorKey(vd.site_index(user_idx), vd.sites_size());
    cv::circle(img, vd.user(user_idx), 5, color);
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

  VoronoiDiagram<double> vd(
      sites.begin(), sites.end(), users.begin(), users.end(),
      resolution.width, resolution.height);
  vd.build();

  if (opts.drawEdges) {
    draw_cells(img, vd);
  }
  draw_sites(img, vd);

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
