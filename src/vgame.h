#pragma once

#include <utility>   // std::copy
#include <algorithm> // std::swap
#include <iterator>  // std::advance
#include <list>
#include <array>
#include <cmath>
#include <stdexcept>

#ifdef DEBUG
#include <iostream>
#endif

#include <opencv2/core/core.hpp> // cv::norm

#include <boost/geometry/geometry.hpp>

#include "maxrect.h"
#include "util.h" // randint

namespace cfla
{
namespace bgi = boost::geometry::index;

// Generic typedefs common for cfla stuff to ease the templating
// of the many classes which require these typedefs.
template<
  class Point, class PointContainer=std::list<Point>,
  unsigned int NumPlayers=2u
>
struct cfla_traits
{
  static const unsigned int nplayers = NumPlayers;
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
  unsigned int current_player = (nplayers-1); // starting with player 2

  // Members
private:
  // List of customer (user) points.
  player_list players_;
  point_list users_;

public:
  VGame(point_list const& users)
    : players_({}), users_(users)
  {}

  template<class InputIter>
  VGame(InputIter users_begin, InputIter users_end)
    : players_({}), users_(users_begin, users_end)
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

  static inline coordinate_type
    distance(const point_type &p1, const point_type &p2) {
      return cv::norm(p2 - p1);
    }

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
      coordinate_type d = distance(facility, user);
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
    for (int user_idx = 0; user_idx < users_size(); ++user_idx)
      owner(user_idx).add_score(1);
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

  // Play the game for the given number of rounds (default: play 1 round).
  // Note that one round means one player is solved (two rounds is symmetric).
  // Note also that this uses the player 2 nth-round solution in each round.
  // Return the last facility added.
  // If starting_player is negative, continue with the last player played.
  // Otherwise, reset such that "player two" is the given player index.
  point_type play_round(size_type rounds=1, int starting_player=-1)
  {
    point_type last_solution;
    unsigned int roundnum = 0u;
    while (rounds--) {
      unsigned int prev_player = current_player-1;
#ifdef DEBUG
      std::cout << "round " << roundnum++ << ":" << std::endl
        << prev_player << ": " << player(prev_player) << std::endl
        << current_player << ": " << player(current_player) << std::endl;
#endif
      // Solve for player 2 and add the point to p2's list and the VD sites.
      last_solution = nth_round(player(prev_player), player(current_player));
      player(current_player).insert(last_solution);
      // Then switch players for next round.
      ++current_player;
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
      point_type user = *userp;
      point_type site = player.nearest_facility(user);

      /* This is actually the distance to the corner points from   o      o o
       * the center.  The left and right points become top-left  o---o =>  \
       * and bottom-right when rotated.                            o      o o */
      auto l1dist = 2.0f * (std::abs(site.x-user.x) + std::abs(site.y-user.y));
      point_type tl = rotateZ2f_pos(point_type(user.x - l1dist, user.y));
      point_type br = rotateZ2f_pos(point_type(user.x + l1dist, user.y));

      // Build the upright rectangle.
      rects_out.push_back(bp::construct<rect_type>(tl.x, tl.y, br.x, br.y));
    }
  }

  // Play a single round and return thei deal location for player 2.
  point_type nth_round(player_type const& p1, player_type const& p2)
  {
    // Build L1 (rotated) rectangles, and rotate the space about the
    // origin to solve with axis-up rectangles.
    std::list<rect_type> l1rects;
    build_rects(p2, l1rects);

    MaxRect<coordinate_type> maxrect(l1rects.begin(), l1rects.end());
    maxrect.compute();

    // We have potentially multiple solution cells.
    // Pick a randome one and return its center.
    rect_type solution = maxrect.cell(randint(size_type(0), maxrect.size()));
    point_type center;
    bp::center(center, solution);
    // Make sure to rotate the solution point back to the original space, since
    // the solution space itself is rotated.
    return rotateZ2f_neg(center);
  }

}; // end class VGame

} // end namespace cfla
