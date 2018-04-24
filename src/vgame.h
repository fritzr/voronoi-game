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

#include "voronoi.h"
#include "maxrect.h"
#include "util.h"

namespace cfla
{

// Generic typedefs common for cfla stuff to ease the templating
// of the many classes which require these typedefs.
template<
  class Point, class PointContainer=std::list<Point>,
  unsigned int NumPlayers=2u,
  class VoronoiDiagram=voronoi::VoronoiDiagram<double>
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
  typedef VoronoiDiagram                      voronoi_diagram;
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
  typedef typename traits::voronoi_diagram voronoi_diagram; \

template<class Traits>
class Facility : public typename Traits::point_type
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

  // Members
private:
  // List of facility locations.
  point_list flist_;
  // Voronoi diagram to locate our nearest facility to each user point.
  voronoi_diagram vd_;

public:
  VPlayer() : {}

  template<class InputIter>
  VPlayer(InputIter facilities_begin, InputIter facilities_end)
    : flist_(facilities_begin, facilities_end)
  {}

  // Iterate over facilities.
  inline iterator begin(void) const { return flist_.begin(); }
  inline iterator end(void) const { return flist_.end(); }
  inline size_type size(void) const { return flist_.size(); }

  inline void add(point_type const& f) { return flist_.push_back(f); }
}; // end class VPlayer

template<class Traits>
class VGame
{
  inherit_traits(Traits);
  typedef VPlayer<traits> player_type;
  typedef std::array<player_type, nplayers> player_list;

  // Members
private:
  // List of customer (user) points.
  player_list players_;

public:
  VGame(point_list const& users)
    : vd_(users.begin(), users.end())
  {}

  template<class InputIter>
  VGame(InputIter users_begin, InputIter users_end)
    : vd_(users_begin, users_end)
  {}

  ~VGame()
  {
    if (vd_)
      delete vd_;
    vd_ = nullptr;
  }

  inline point_citerator users_begin(void) const { return users_.cbegin(); }
  inline point_citerator users_end(void) const { return users_.cend(); }
  inline size_type users_size(void) const { return users_.size(); }

  inline voronoi_diagram& voronoi(void) { return vd_; }

  // Initialize a player with his facilities.
  template<class InputIter>
  inline void init_player(unsigned int id, InputIter fbegin, InputIter fend) {
    players_[playerid] = player_type(fbegin, fend);
    auto player = &players_[playerid];
    vd.add_sites(player->begin(), players_->end());
  }

  // Return the player congruent to id modulo the number of players.
  inline player_type& player(int id) {
    return players_[id % players_.size()];
  }

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
      vd_.build();
      // Solve for player 2 and add the point to p2's list and the VD sites.
      last_solution = nth_round(player(player_id-1), player(player_id));
      player(player_id).add(last_solution);
      vd_.add_site(last_solution);
      // Then switch players for next round.
      ++player_id;
    }
    return p2_solution;
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
    for (auto userp = vd_.users_begin(); userp != vd_.users_end(); ++userp)
    {
      point_type user = *userp;
      point_type site = vd_.nearest_site(userp);

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
