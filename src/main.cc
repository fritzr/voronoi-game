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
#include <type_traits>

#include "main.h"
#include "voronoi.h"
#include "maxrect.h"
#include "maxtri.h"
#include "vgame.h"
#include "shpReader.h"
#include "adapt_boost_poly.h"

#define i2f(...) static_cast<float>(__VA_ARGS__)
#define i2u(...) static_cast<unsigned int>(__VA_ARGS__)

using namespace std;
using namespace shp;

typedef double coordinate_type; // required precision of shapefiles
typedef typename boost::geometry::model::d2::point_xy<coordinate_type>
  point_type;

// Nearest-neighbor solvers for L1 (Manhattan) and L2 (Euclidean) distance
typedef typename cfla::L2NN1<coordinate_type, point_type> L2NN1;
typedef typename cfla::L1NN1<coordinate_type, point_type> L1NN1;

// L1 max-rect solver (Imai & Asano, 1983)
// XXX This hasn't been instantiated in a while.
// MaxTri is the main contribution.
typedef typename cfla::rect::MaxRectSolver<coordinate_type, L1NN1>
  rect_solver_type;
typedef typename cfla::cfla_traits<point_type, rect_solver_type> rect_traits;
typedef typename cfla::VGame<rect_traits> RectGame;

// Travel time max-tri solver
typedef User<point_type> User2d; // Customer points with isochrones
typedef typename cfla::tri::MaxTriSolver<coordinate_type> tri_solver_type;
typedef typename cfla::cfla_traits<point_type, tri_solver_type, User2d>
  tri_traits;
typedef typename cfla::VGame<tri_traits> TriGame;


// Enum to int: initialize positive types as zero, negative types as -1
template <typename T>
struct to_int {
  static const bool signum = std::is_signed<T>::value;

  typedef typename std::make_signed<T>::type signed_type;
  typedef typename std::make_unsigned<T>::type unsigned_type;

  typedef typename std::conditional<signum,
      typename std::make_signed<T>::type,
      typename std::make_unsigned<T>::type
        >::type
    type;

  static inline type initval(void) {
    if (signum)
      return -1;
    else
      return 0u;
  }
};

template <typename T>
static T
getenum(T maxval, const char *instr, ostream &err, const char *errtype)
{
  typedef typename to_int<T>::type int_type;
  istringstream in(instr);
  const int_type initval = to_int<T>::initval();
  int_type ival = initval;
  int_type imax = int_type(maxval);
  in >> ival;
  if (!in.eof() || in.fail() || ival < 0 || ival >= imax) {
    err << "bad value for " << errtype << ": '" << instr << "'";
    opterr = 1;
    return T(initval);
  }
  return T(ival);
}

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

enum SolverAlgorithm {
  ALG_RECT = 0,
  ALG_TRI = 1,
  ALG_BAD = 2,
};

struct options_t
{
  // Global configurable options
  bool debug = false;
  DrawCell computeCell = CELL_NONE;
  DrawRects drawRects = RECTS_NONE;
  bool dumpGraph = false; // adjacency graph
  unsigned int rounds = 1u;
  bool userLabels = false;
  bool interactive = false;
  string output_path;
  long unsigned int seed = 0u;
  SolverAlgorithm algorithm = ALG_TRI;

  string user_points_path;
  string user_rings_path;
  string p1sites_path;
  string p2sites_path;
};

static options_t opts;

static const char *sopts = "Dp:s:";

static void usage(const char *prog, const char *errmsg=NULL)
#ifndef _MSC_VER
  __attribute__((noreturn))
#endif
	;

void
usage(const char *prog, const char *errmsg)
{
  cerr << "usage: " << prog
    << " [OPTIONS] [-p<N>] <user_points> <user_rings> <p1_sites> <p2_sites>"
    << " <output>"
    << endl << endl
    << "Read users and their fixed-travel-time (FTT) rings, along with some "
       "known facilities/sites for players 1 and 2. Then successively find the"
       " optimal location to place sites for the second player. With -p, play"
       " N rounds, considering each player alternatively as 'player 2'. The"
       " output is a Point Shapefile containing the solution points."
    << endl << endl
    << "OPTIONS are:" << endl
    << "  -D			Debug mode (extra output)" << endl
    << "  -p N			Play N rounds of the game (default: 1)." << endl
    << "  -s SEED		Seed the RNG with the given number." << endl
    << "			Default is to use current Epoch time." << endl
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

// Parses options and sets the global options structure.
static void
get_options(int argc, char *argv[], options_t &o)
{
  int opt;
  opterr = 0;
  ostringstream errstr;
  while ((opt = getopt(argc, argv, sopts)) >= 0 && !opterr)
  {
    switch (opt)
    {
    case 'D':
      o.debug = true; break;
    case 'd':
      o.algorithm = getenum(ALG_BAD, optarg, errstr, "distance algorithm");
      break;
    case 'p':
      o.rounds = strtoul(optarg, NULL, 0); break;
    case 's':
      if (rng != nullptr)
        usage(argv[0], "-s given multiple times");
      o.seed = getenum(ULONG_MAX, optarg, errstr, "seed");
      rng = new rng_type(o.seed);
      break;
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

  if (o.algorithm == ALG_RECT)
    usage(argv[0], "-d: interface to the"
        " manhattan distance algorithm is currently unimplemented");

  int nargs = 5;
  if (argc - optind < nargs) {
    usage(argv[0], "not enough arguments");
  }

  o.user_points_path = string(argv[optind]);
  o.user_rings_path = string(argv[optind+1]);
  o.p1sites_path = string(argv[optind+2]);
  o.p2sites_path = string(argv[optind+3]);
  o.output_path = string(argv[optind+4]);
}


#ifdef DEBUG
static ostream&
dumpDistances(ostream& os, const vector<User2d> &users,
    const vector<point_type> &p1sites, const vector<point_type> &p2sites)
{
  os << "===== USERS =====" << endl;
  for (auto uit = users.begin(); uit != users.end(); ++uit)
  {
    os << "    " << *uit << endl;

    // output shortest distance to each ring as a fuzzy metric
    unsigned int nlayer = 0u;
    os << "      ";
    for (auto rit = uit->begin(); rit != uit->end(); ++rit)
    {
      if (nlayer != 0u)
        os << ", ";
      os << User2d::distance(uit->center(), rit->front());
      ++nlayer;
    }
    os << endl;
  }

  typedef vector<point_type>::const_iterator piterator;
  typedef pair<const piterator, const piterator> iter_range;
  typedef pair<const string, const iter_range> zip;
  array<zip, 2> players = {
      zip(string("===== PLAYER ONE ====="),
          iter_range(p1sites.cbegin(), p1sites.cend())),
      zip(string("===== PLAYER TWO ====="),
          iter_range(p2sites.cbegin(), p2sites.cend())),
  };

  /* Output both players..  */
  for (auto pn = players.cbegin(); pn != players.cend(); ++pn)
  {
    const string &title(pn->first);
    os << endl << title << endl;

    unsigned int idx = 0u;
    const iter_range &range(pn->second);
    for (auto pp = range.first; pp != range.second; ++pp)
    {
      os << "[" << setw(2) << setfill(' ') << idx++ << "]    " << *pp << endl;
      unsigned int uidx = 0u;
      for (auto uit = users.begin(); uit != users.end(); ++uit)
      {
        os << "    (" << setw(2) << uidx++ << ") "
          << *uit << " | FTT " << uit->travelTime(*pp) << endl;
      }
    }
  }

  return os;
}
#endif

template<class Game> /* VGame */
static void
play_game(vector<User2d> const& users,
    vector<point_type> const& p1sites, vector<point_type> const& p2sites,
    vector<point_type> &solutions)
{
  Game vg(users.begin(), users.end());

  vg.init_player(0, p1sites.begin(), p1sites.end());
  vg.init_player(1, p2sites.begin(), p2sites.end());

  //const typename TriGame::player_type& p1 = vg.player(0);
  //const typename TriGame::player_type& p2 = vg.player(1);

  unsigned int turns_remaining = opts.rounds;
  unsigned int rounds_per_turn = 1u;
  while (turns_remaining--)
  {
    vg.score();
#ifdef DEBUG
    unsigned int nextp = vg.next_player().id();
#endif
    point_type next_facility = vg.play_round(rounds_per_turn);
    solutions.push_back(next_facility);
#ifdef DEBUG
    cout << "player " << nextp << " solution at " << next_facility << endl;
#endif
  }
}

std::default_random_engine* rng = nullptr;

int main(int argc, char *argv[])
{
  get_options(argc, argv, opts);

  vector<User2d> users(readUsers<point_type>(
        opts.user_points_path, opts.user_rings_path));
  vector<point_type> p1sites(readPoints<point_type>(opts.p1sites_path));
  vector<point_type> p2sites(readPoints<point_type>(opts.p2sites_path));
  vector<point_type> solutions;

#ifdef DEBUG
  dumpDistances(cout, users, p1sites, p2sites);
#endif

  // Play the game one round at a time.
  if (opts.algorithm == ALG_TRI)
    play_game<TriGame>(users, p1sites, p2sites, solutions);
  else
    /*play_game<RectGame>(users, p1sites, p2sites, solutions)*/
    ;

  writePoints(opts.output_path, solutions.begin(), solutions.end());
  cerr << "Wrote " << opts.output_path << endl;

  if (rng != nullptr)
    delete rng;

  return 0;
}
