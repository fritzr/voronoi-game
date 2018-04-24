#pragma once

#include <utility> // std::copy
#include <algorithm> // std::swap
#include <list>
#include <array>
#include <map>
#include <cmath>
#include <stdexcept>

#include <boost/range/join.hpp> // boost::join
#include <opencv2/core/core.hpp>

#include "maxrect.h"
#include "util.h"

namespace cfla
{

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
class Facility : public Traits::point_type
{
  inherit_traits(Traits);

  // Membrs
public:
  int player_id;

  // Methods
public:
  Facility(coordinate_type x, coordinate_type y)
    : point_type(x, y), player_id(-1) {}
  Facility(coordinate_type x, coordinate_type y, unsigned int pid)
    : point_type(x, y), player_id(pid) {}
};


template<class Traits>
class VPlayer
{
  inherit_traits(Traits);
  typedef point_citerator iterator;
  typedef point_citerator const_iterator;
  typedef point_type      value_type;
  typedef NN1<point_list, point_list> VNN;

  // Members
private:
  // Player ID.
  int id_;
  // List of facility locations.
  point_list flist_;
  // Voronoi diagram to locate our nearest facility to each user point.
  VNN nn_;
  // Number of customers we 'own' (closest facility is ours).
  int score_;

public:
  template<class UserList, class InputIter>
  VPlayer(int id, UserList const& users, InputIter fbegin, InputIter fend)
    : id_(id), flist_(fbegin, fend), nn_(users, flist_), score_(-1)
  {}

  inline int id(void) const { return id_; }

  // Iterate over facilities.
  inline iterator begin(void) const { return flist_.begin(); }
  inline iterator end(void) const { return flist_.end(); }
  inline size_type size(void) const { return flist_.size(); }

  inline void set_score(int new_score) { score_ = new_score; }
  inline int score(void) const { return score_; }
  inline void add_score(int to_add=1) { set_score(score()+to_add); }

  inline VNN const& nn(void) const { return nn_; }
  inline void build(void) { nn_.build(); }

  inline void add(point_type const& f) { nn_.add_site(f); }
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

public:
  VGame(point_list const& users)
    : users_(users)
  {}

  template<class InputIter>
  VGame(InputIter users_begin, InputIter users_end)
    : users_(users_begin, users_end)
  {}

  inline point_citerator users_begin(void) const { return users_.cbegin(); }
  inline point_citerator users_end(void) const { return users_.cend(); }
  inline size_type users_size(void) const { return users_.size(); }

  // Initialize a player with his facilities.
  template<class InputIter>
  inline void init_player(unsigned int id, InputIter fbegin, InputIter fend) {
    if (players_[id])
      delete players_[id];
    players_[id] = new player_type(users_, fbegin, fend);
  }

  // Return the player congruent to id modulo the number of players.
  // If the player has not been initialized this will segfault (hah).
  inline const player_type& player(int id) const {
    return *players_[id % players_.size()];
  }

protected:
  // Modifiable reference.
  inline player_type& player(int id) {
    return const_cast<player_type&>(player(id));
  }

public:

  static inline coordinate_type
    distance(const point_type &p1, const point_type &p2) {
      return cv::norm(p2 - p1);
    }

  // To find who owns a point, get each player's closest facility and narrow
  // down to which one is really closer. This prevents us from having to
  // keep an extra range tree (since we already have one per player).
  const player_type& owner(int user_index) const
  {
    coordinate_type dist = std::numeric_limits<coordinate_type>::infinity();
    player_type* min_player = &player(0);
    point_type user = users_.at(user_index);
    for (auto playerp = players_.begin(); playerp != players_.end(); ++playerp)
    {
      point_type facility = (*playerp)->nn().nearest_facility(user_index);
      coordinate_type d = distance(facility, user);
      if (std::abs(d) < std::abs(dist)) {
        dist = d;
        min_player = *playerp;
      }
    }
    return min_player;
  }

protected:
  inline player_type& owner(int user_index) {
    return const_cast<player_type&>(owner(user_index));
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

protected:
  inline player_type& winner(void) {
    return const_cast<player_type&>(winner());
  }

public:

  // Play the game for the given number of rounds (default: play 1 round).
  // Note that one round means one player is solved (two rounds is symmetric).
  // Note also that this uses the player 2 nth-round solution in each round.
  // Return the last facility added.
  point_type play_round(size_type rounds=1)
  {
    player_type* p1ptr = &players_[0];
    int player_id = 1; // starting with player 2
    point_type last_solution;
    while (rounds--) {
      // Each round we need a new NN search since we may have added a facility.
      player(player_id).build();
      // Solve for player 2 and add the point to p2's list and the VD sites.
      last_solution = nth_round(player(player_id-1), player(player_id));
      player(player_id).add(last_solution);
      // Then switch players for next round.
      ++player_id;
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
  void build_rects(std::list<rect_type> &rects_out)
  {
    for (auto userp = users_begin(); userp != users_end(); ++userp)
    {
      point_type user = *userp;
      point_type site = player(1).nn().nearest_facility(userp);

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
    build_rects(l1rects);

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
