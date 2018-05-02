#pragma once

#include <utility>   // std::copy
#include <algorithm> // std::swap
#include <iterator>  // std::advance
#include <list>
#include <array>
#include <cmath>
#include <stdexcept>

#include <iostream>

#include <opencv2/core/core.hpp> // cv::norm

#include <boost/geometry/geometry.hpp>

#include "maxrect.h"
#include "util.h" // randrange

namespace cfla
{
namespace bgi = boost::geometry::index;

// Generic typedefs common for cfla stuff to ease the templating
// of the many classes which require these typedefs.
template<
  class Point, class PointContainer=std::list<Point>,
  unsigned int NumPlayers=2u,
  unsigned int DistanceNorm=1u
>
struct cfla_traits
{
  static const unsigned int nplayers = NumPlayers;
  static const unsigned int norm     = DistanceNorm;
  typedef Point                               point_type;
  typedef decltype(point_type::x)             coordinate_type;
  typedef typename bp::rectangle_data<coordinate_type> rect_type;
  typedef PointContainer                      point_list;
  typedef typename point_list::iterator       point_iterator;
  typedef typename point_list::const_iterator point_citerator;
  typedef typename point_list::size_type      size_type;
};

#define inherit_traits(t) \
public: \
  typedef t traits; \
  static const auto nplayers = t::nplayers; \
  static const auto norm     = t::norm; \
  typedef typename traits::point_type      point_type; \
  typedef typename traits::coordinate_type coordinate_type; \
  typedef typename traits::rect_type       rect_type; \
  typedef typename traits::point_list      point_list; \
  typedef typename traits::point_iterator  point_iterator; \
  typedef typename traits::point_citerator point_citerator; \
  typedef typename traits::size_type       size_type; \

template<class Traits>
class VPlayer
{
  inherit_traits(Traits);
  typedef typename bgi::rtree<point_type, bgi::quadratic<16> > nn_tree;
  typedef typename nn_tree::value_type     value_type;
  typedef typename nn_tree::const_iterator const_iterator;
  typedef typename nn_tree::const_iterator iterator;

  // Members
private:
  // Player ID.
  const int id_;
  // Data structure for locating the nearest facility to each user point.
  nn_tree nn_;
  // Number of customers we 'own' (closest facility is ours).
  int score_;

public:
  template<class UserList, class InputIter>
  VPlayer(int id, UserList const& users, InputIter fbegin, InputIter fend)
    : id_(id), nn_(fbegin, fend), score_(-1)
  {}

  inline int id(void) const { return id_; }

  // Iterate over facilities.
  inline const_iterator begin(void) const { return nn_.begin(); }
  inline const_iterator end(void) const { return nn_.end(); }
  inline size_type size(void) const { return nn_.size(); }

  inline void set_score(int new_score) { score_ = new_score; }
  inline int score(void) const { return score_; }
  inline void add_score(int to_add=1) { set_score(score()+to_add); }

  inline point_type const& nearest_facility(point_type const& user) const {
    return *nn_.qbegin(bgi::nearest(user, 1));
  }
  inline void insert(point_type const& f) { nn_.insert(f); }

#ifdef DEBUG
  template <class Traits_>
  friend inline std::ostream& operator<<(std::ostream& os,
      VPlayer<Traits_> const& p)
  {
    for (auto ptit = p.begin(); ptit != p.end(); ++ptit)
      std::cout << *ptit << " ";
    std::cout << std::endl;
    return os;
  }
#endif
}; // end class VPlayer

template<class Traits>
class VGame
{
  inherit_traits(Traits);
  typedef VPlayer<traits> player_type;
  typedef std::array<player_type*, nplayers> player_list;

  // Members
private:
  // List of customer (user) points.
  player_list players_;
  point_list users_;
  unsigned int current_round;
#ifdef DEBUG
  cv::Mat* img_ = nullptr;
#endif

public:
  VGame(point_list const& users)
    : players_({}), users_(users), current_round(0)
  {}

  template<class InputIter>
  VGame(InputIter users_begin, InputIter users_end)
    : players_({}), users_(users_begin, users_end), current_round(0)
  {}

  ~VGame()
  {
    for (size_type playerid = 0u; playerid < nplayers; ++playerid)
    {
      if (players_[playerid])
      {
        delete players_[playerid];
        players_[playerid] = nullptr;
      }
    }
  }

#ifdef DEBUG
  inline void set_img(cv::Mat& img) { img_ = &img; }
#endif

  inline point_type user(int user_idx) const {
    return std::advance(users_begin(), user_idx);
  }
  inline point_citerator users_begin(void) const { return users_.cbegin(); }
  inline point_citerator users_end(void) const { return users_.cend(); }
  inline size_type users_size(void) const { return users_.size(); }

  // Initialize a player with his facilities.
  template<class InputIter>
  inline void init_player(unsigned int id, InputIter fbegin, InputIter fend) {
    if (players_[id])
      delete players_[id];
    players_[id] = new player_type(id, users_, fbegin, fend);
  }

  // Return the player congruent to id modulo the number of players.
  // If the player has not been initialized this will segfault (hah).
  inline const player_type& player(int id) const {
    return *players_[id % players_.size()];
  }

  // Modifiable reference.
  inline player_type& player(int id) {
    return const_cast<player_type&>(const_cast<const VGame*>(this)->player(id));
  }

public:

  template<int Norm>
  static inline
  typename std::enable_if<Norm == 1u, coordinate_type>::type
    distance(const point_type &p1, const point_type &p2) {
      // L1 distance
      return std::abs(p1.x - p2.x) + std::abs(p1.y - p2.y);
    }

  template<int Norm>
  static inline
  typename std::enable_if<Norm == 2u, coordinate_type>::type
    distance(const point_type &p1, const point_type &p2) {
      // L2 distance
      return cv::norm(p2 - p1);
    }

  static inline coordinate_type
    distance(const point_type& p1, const point_type& p2)
      { return distance<VGame::norm>(p1, p2); }

  // To find who owns a point, get each player's closest facility and narrow
  // down to which one is really closer. This prevents us from having to
  // keep an extra range tree (since we already have one per player).
  const player_type& owner(point_type const& user) const
  {
    coordinate_type dist = std::numeric_limits<coordinate_type>::infinity();
    const player_type* min_player = &player(0);
    for (auto playerp = players_.begin(); playerp != players_.end(); ++playerp)
    {
      point_type facility = (*playerp)->nearest_facility(user);
      coordinate_type d = VGame::distance(facility, user);
      if (std::abs(d) < std::abs(dist)) {
        dist = d;
        min_player = *playerp;
      }
    }
    return *min_player;
  }

  inline player_type& owner(point_type const& user) {
    return const_cast<player_type&>(const_cast<const VGame*>(this)->owner(user));
  }

public:

  // Score all players.
  void score(void)
  {
    for (auto playerp = players_.begin(); playerp != players_.end(); ++playerp)
      (*playerp)->set_score(0);
    for (auto userp = users_begin(); userp != users_end(); ++userp)
      owner(*userp).add_score(1);
  }

  // Return the winning player. You must run score() first.
  const player_type& winner(void) const
  {
    int max_score = -1;
    int winning_player = -1;
    for (unsigned int playerid = 0u; playerid < nplayers; ++playerid)
    {
      int score = players_[playerid]->score();
      if (score > max_score)
      {
        max_score = score;
        winning_player = playerid;
      }
    }
    return winning_player;
  }

  inline player_type& winner(void) {
    return const_cast<player_type&>(winner());
  }

public:

  // The number of the next round to be played (the first round is round 0).
  inline unsigned int next_round(void) const { return current_round; }
  inline player_type const& next_player(void) const {
    return player(current_round + nplayers - 1);
  }

  // Play the game for the given number of rounds (default: play 1 round).
  // Note that one round means one player is solved (two rounds is symmetric).
  // Note also that this uses the player 2 nth-round solution in each round.
  // Return the last facility added.
  // If starting_player is negative, continue with the last player played.
  // Otherwise, reset such that "player two" is the given player index.
  point_type play_round(size_type rounds=1, int starting_player=-1)
  {
    point_type last_solution;
    while (rounds--) {
      unsigned int current_player = next_player().id();
      unsigned int prev_player = (current_player - 1) % nplayers;
#ifdef DEBUG
      std::cout << "round " << current_round << ":" << std::endl
        << prev_player << ": " << player(prev_player)
        << current_player << ": " << player(current_player);
#endif
      // Solve for player 2 and add the point to p2's list and the VD sites.
      last_solution = nth_round(player(prev_player), player(current_player));
      player(current_player).insert(last_solution);
      // Then switch players for next round.
      ++current_round;
    }
    return last_solution;
  }

  // Shortcut for init_player() and then play_round(rounds) with 2 players.
  template<class InputIter1, class InputIter2>
  inline void play2(
      InputIter1 p1f_begin, InputIter1 p1f_end,
      InputIter2 p2f_begin, InputIter2 p2f_end,
      size_type rounds=1) {
    init_player(0, p1f_begin, p1f_end);
    init_player(1, p2f_begin, p2f_end);
    play_round(rounds);
  }

private:

  // Build an "L1 ball" (square) for each user point to pass to the maxrect
  // algorithm. The rectangles (squares) are formed using the L1 distance from
  // the user to its nearest facility.
  void build_rects(const player_type& player, std::list<rect_type> &rects_out)
  {
    for (auto userp = users_begin(); userp != users_end(); ++userp)
    {
      // Only build rects for points owned by the selected player.
      if (owner(*userp).id() != player.id())
        continue;

      point_type user = *userp;
      point_type site = player.nearest_facility(user);

      /* This is actually the distance to the corner points from   o      o o
       * the center.  The left and right points become top-left  o---o =>  \
       * and bottom-right when rotated.                            o      o o */
      auto l1dist = VGame::distance(site, user);
      point_type tl = point_type(user.x - l1dist, user.y);
      point_type br = point_type(user.x + l1dist, user.y);
#if 0
//#ifdef DEBUG
      point_type tr = point_type(user.x, user.y + l1dist);
      point_type bl = point_type(user.x, user.y - l1dist);
      if (img_) {
        cv::line(*img_, tl, tr, cv::Scalar(0xcc, 0xcc, 0xcc), 1);
        cv::line(*img_, tr, br, cv::Scalar(0xcc, 0xcc, 0xcc), 1);
        cv::line(*img_, br, bl, cv::Scalar(0xcc, 0xcc, 0xcc), 1);
        cv::line(*img_, bl, tl, cv::Scalar(0xcc, 0xcc, 0xcc), 1);
        cv::imshow("rectangles", *img_);
        cv::waitKey(0);
      }
#endif
      tl = rotateZ2f_pos(tl);
      br = rotateZ2f_pos(br);

      // Build the upright rectangle.
      rect_type l1rect;
      l1rect = bp::construct<rect_type>(tl.x, tl.y, br.x, br.y);
      rects_out.push_back(l1rect);
    }
  }

  // Play a single round and return thei deal location for player 2.
  point_type nth_round(player_type const& p1, player_type const& p2)
  {
    // To solve for player 2, build L1 (rotated) rectangles from the
    // opposing player's facilities and solve for max-depth region(s).
    // We first rotate the space about the origin so we solve with axis-up
    // rectangles.
    std::list<rect_type> l1rects;
    build_rects(p1, l1rects);

    MaxRect<coordinate_type> maxrect(l1rects.begin(), l1rects.end());
    maxrect.compute();

    // We have potentially multiple solution cells.
    // Pick a randome one and return its center.
    if (maxrect.size() == 0)
    {
      std::cerr << "warning: empty solution" << std::endl;
      return point_type(0, 0);
    }

    size_type chosen_one = randrange(size_type(0), maxrect.size()-1);
#ifdef DEBUG
    // Dump all possible solutions.
    std::cerr << "solutions:" << std::endl;
    size_type cellidx = size_type(0);
    for (auto cellp = maxrect.begin(); cellp != maxrect.end(); ++cellp)
    {
      rect_type solution = *cellp;
      point_type ctmp(0, 0);
      bp::center(ctmp, solution);
      std::cerr << "[" << std::setw(2) << std::setfill(' ') << std::dec
        << cellidx << "] (depth=" << maxrect.solution(cellidx).depth()
        << ") " << ctmp << std::endl;
    }
#endif

    // Instead of the center, choose a random point in the cell.
    rect_type solution = maxrect.cell(chosen_one);

    coordinate_type randx = randrange(
        bp::get(solution, bp::HORIZONTAL, bp::LOW),
        bp::get(solution, bp::HORIZONTAL, bp::HIGH));
    coordinate_type randy = randrange(
        bp::get(solution, bp::VERTICAL, bp::LOW),
        bp::get(solution, bp::VERTICAL, bp::HIGH));
    point_type next(randx, randy);;
#ifdef DEBUG
    // Dump our chosen solution.
    std::cerr << "chose [" << std::setw(2) << std::setfill(' ') << std::dec
      << chosen_one << "] " << next << std::endl;
#endif

    // Make sure to rotate the solution point back to the original space, since
    // the solution space itself is rotated.
    return rotateZ2f_neg(next);
  }

}; // end class VGame

} // end namespace cfla
