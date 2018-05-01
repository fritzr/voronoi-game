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
#include "maxrect.h"
#include "vgame.h"

#define i2f(...) static_cast<float>(__VA_ARGS__)
#define i2u(...) static_cast<unsigned int>(__VA_ARGS__)

using namespace std;
using namespace cv;
using namespace voronoi;

typedef typename cfla::MaxRect<double> MaxRect;
typedef typename MaxRect::coordinate_type coordinate_type;
typedef typename bp::rectangle_data<coordinate_type> rect_type;
typedef typename MaxRect::solution_type solution_type;

typedef typename cfla::cfla_traits<cv::Point2f> vgtraits;
typedef typename cfla::VGame<vgtraits> VGame;

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
  DrawRects drawRects = RECTS_NONE;
  bool dumpGraph = false; // adjacency graph
  unsigned int rounds = 1u;
  bool userLabels = false;
  bool interactive = false;

  string users_path;
  string p1sites_path;
  string p2sites_path;
};

static options_t opts;

// Default colors and fonts
static const Scalar gridColor(Scalar::all(0xff)); // white
static const int FONT = FONT_HERSHEY_SCRIPT_SIMPLEX;
static const double FONT_SCALE = 0.5;
static const int FONT_THICKNESS = 1;
static const Scalar FONT_COLOR(Scalar::all(0xff)); // white
const int FONT_XSPACE = 10;
const int FONT_YSPACE = 10;

static const Scalar P1COLOR(Scalar(255, 0, 0));
static const Scalar P2COLOR(Scalar(0, 0, 255));

static void usage(const char *prog, const char *errmsg=NULL)
#ifndef _MSC_VER
  __attribute__((noreturn))
#endif
	;

void
usage(const char *prog, const char *errmsg)
{
  cerr << "usage: " << prog << " [OPTIONS] <users> <sites>" << endl
    << "       " << prog << " [OPTIONS] -p<N> <users> <p1sites> <p2sites>"
    << endl
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
    << "  -u			Draw user labels (default: yes)" << endl
    << "  -U                    Do not draw user labels" << endl
    << "  -W SIZE		Max width of the display (default: 1920)"
    << endl
    << "  -H SIZE		Max height of the display (default: 1080)"
    << endl
    << "  -q TYPE		Query type for mapping users to sites:" << endl
    << "    [0:default], 1:brute force, 2:range search, 3:NN(1) (fastest)"
    << endl
    << "  -r			Draw rects (default: none):" << endl
    << "  -p N			Play N rounds of the game (default: 1)." << endl
    << "  -i			Play interactively (default: no)" << endl
    << "  -I			Do not play interactively (default)" << endl
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

static const char *sopts = "FC:X:Y:T:eEgGlLdW:H:hs:u:q:p:iI";

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
    case 'g':
      o.dogrid = true; break;
    case 'G':
      o.dogrid = false; break;
    case 'l':
      o.gridLabels = true; break;
    case 'L':
      o.gridLabels = false; break;
    case 'u':
      o.userLabels = true; break;
    case 'U':
      o.userLabels = false; break;
    case 'T':
      o.gridThickness = getenum<int>(INT_MAX, optarg, errstr, "grid thickness");
      break;
    case 'W':
      o.screenWidth = getenum<int>(INT_MAX, optarg, errstr, "screen width");
      break;
    case 'H':
      o.screenHeight = getenum<int>(INT_MAX, optarg, errstr, "screen height");
      break;
    case 'r':
      o.drawRects = getenum<DrawRects>(RECTS_BAD, optarg, errstr, "rect type");
      break;
    case 'j':
      o.dumpGraph = true; break;
    case 'p':
      o.rounds = strtoul(optarg, NULL, 0); break;
    case 'i':
      o.interactive = true; break;
    case 'I':
      o.interactive = false; break;
    case 'h':
      usage(argv[0]);
    case ':':
      errstr << "option -" << (char)opt << " missing argument";
      if (!opterr) opterr = 1;
      break;
    case '?':
      errstr << "unknown option -" << (char)optopt;
      if (!opterr) opterr = 1;
      break;
    default:
      errstr << "unknown return from getopt: " << opt;
      if (!opterr) opterr = 1;
      break;
    }
  }
  if (opterr) {
    usage(argv[0], errstr.str().c_str());
  }

  if (o.computeCell)
    usage(argv[0], "-c is no longer implemented");
  if (o.drawRects)
    usage(argv[0], "-r is no longer implemented");

  int nargs = 3;
  if (argc - optind < nargs) {
    usage(argv[0], "not enough arguments");
  }

  // Unlimited screen size when set to zero.
  if (o.screenHeight == 0) o.screenHeight = INT_MAX;
  if (o.screenWidth == 0) o.screenWidth = INT_MAX;

  o.users_path = string(argv[optind]);
  o.p1sites_path = string(argv[optind+1]);
  if (nargs > 2)
    o.p2sites_path = string(argv[optind+2]);
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

template<class InputIter>
static void
draw_facilities(Mat& img, InputIter begin, InputIter end, Scalar const& color)
{
  while (begin != end)
    cv::circle(img, *begin++, 10, color, -1);
}

static void
draw_users(Mat& img, VGame const& vd)
{
  // Don't fill users
  for (auto userp = vd.users_begin(); userp != vd.users_end(); ++userp) {
    typename VGame::player_type const& p = vd.owner(*userp);
    Scalar color = p.id() == 0 ? P1COLOR : P2COLOR;
    cv::circle(img, check_rot(*userp), 5, color);
  }
}

static void
show_score(Mat& img, VGame& vg)
{
  // Draw updated users and scores based on new owners.
  vg.score();
  draw_users(img, vg);
  cout << "scores:" << endl;
  for (unsigned int pid = 0u; pid < VGame::nplayers; ++pid)
  {
    const typename VGame::player_type& player = vg.player(pid);
    cout << "  player " << player.id() << " = "
      << setw(2) << setfill(' ') << player.score() << endl;
  }
}

int main(int argc, char *argv[])
{
  get_options(argc, argv, opts);

  vector<Point> users, p1sites, p2sites;

  readOrDie(opts.users_path, back_inserter(users));
  readOrDie(opts.p1sites_path, back_inserter(p1sites));
  if (p1sites.empty()) {
    cerr << "error: empty sites file (" << opts.p1sites_path << ")" << endl;
    exit(2);
  }
  readOrDie(opts.p2sites_path, back_inserter(p2sites));
  if (p2sites.empty()) {
    cerr << "error: empty sites file (" << opts.p2sites_path << ")" << endl;
    exit(2);
  }

  // We got resolution (x, y) as (width, height); flip to (rows, cols)
  Size resolution(opts.screenWidth, opts.screenHeight);
  Mat img(resolution.height, resolution.width, CV_8UC3);

  if (opts.debug)
    cout << "canvas size " << img.cols << " x " << img.rows << endl;

  // Draw a grid first if we want.
  if (opts.dogrid)
    drawGrid(img, opts.ngridx, opts.ngridy, gridColor, opts.gridThickness,
        opts.gridLabels);

  // Play the game one round at a time.
  VGame vg(users.begin(), users.end());
#ifdef DEBUG
  vg.set_img(img);
#endif
  vg.init_player(0, p1sites.begin(), p1sites.end());
  vg.init_player(1, p2sites.begin(), p2sites.end());
  const typename VGame::player_type& p1 = vg.player(0);
  const typename VGame::player_type& p2 = vg.player(1);
  draw_facilities(img, p1.begin(), p1.end(), P1COLOR);
  draw_facilities(img, p2.begin(), p2.end(), P2COLOR);
  putText(img, "P1", Point(FONT_XSPACE, 2*FONT_YSPACE),
      FONT, FONT_SCALE, P1COLOR, FONT_THICKNESS);
  putText(img, "P2", Point(FONT_XSPACE+40, 2*FONT_YSPACE),
      FONT, FONT_SCALE, P2COLOR, FONT_THICKNESS);

  unsigned int turns_remaining = (unsigned int)(opts.rounds != 0);
  unsigned int rounds_per_turn = opts.rounds;
  if (opts.interactive)
  {
    turns_remaining = opts.rounds;
    rounds_per_turn = 1u;
  }
  Scalar pcolor = P1COLOR;
  while (turns_remaining--)
  {
    show_score(img, vg);
    pcolor = (vg.next_round() % 2 == 0) ? P2COLOR : P1COLOR;
    unsigned int nextp = vg.next_player().id();
    Point2f next_facility = vg.play_round(rounds_per_turn);
    cout << "player " << nextp << " solution at " << next_facility << endl;
    if (opts.interactive)
    {
      // First show a white circle...
      circle(img, next_facility, 10, FONT_COLOR, -1);
      imshow("output", img);
      waitKey(0);
    }
    // Then show it as the right color once we've moved on.
    circle(img, next_facility, 10, pcolor, -1);
  }
  show_score(img, vg);

  if (opts.userLabels)
  {
    unsigned int user_idx = 0u;
    for (auto uit = users.begin(); uit != users.end(); ++uit)
    {
      Point up = check_rot(*uit);
      putText(img, to_string(user_idx++),
          Point(up.x + FONT_XSPACE, up.y - FONT_YSPACE),
          FONT, FONT_SCALE, FONT_COLOR, FONT_THICKNESS);
    }
  }

  imshow("output", img);
  waitKey(0);

  return 0;
}
