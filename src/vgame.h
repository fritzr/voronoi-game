#include <utility> // std::copy
#include <algorithm> // std::swap
#include <boost/range/join.hpp> // boost::join
#include <list>

#include "voronoi.h"

namespace cfla
{

// Generic typedefs common for cfla stuff to ease the templating
// of the many classes which require these typedefs.
template<
  class Point, class PointContainer=std::list<Point>,
  class VoronoiDiagram=voronoi::VoronoiDiagram<double>
>
struct cfla_traits
{
  typedef Point                               point_type;
  typedef PointContainer                      point_list;
  typedef typename point_list::iterator       point_iterator;
  typedef typename point_list::const_iterator point_citerator;
  typedef typename point_list::size_type      size_type
  typedef VoronoiDiagram                      voronoi_diagram;
};

#define inherit_traits(t) \
public: \
  typedef t traits; \
  typedef typename traits::point_type      point_type; \
  typedef typename traits::point_list      point_list; \
  typedef typename traits::iterator        iterator; \
  typedef typename traits::const_iterator  const_iterator; \
  typedef typename traits::size_type       size_type \
  typedef typename traits::voronoi_diagram voronoi_diagram; \


template<class Traits>
class VPlayer
{
  inherit_traits(Traits);
  typedef point_citerator iterator;
  typedef point_type      value_type;

  // Members
private:
  // List of facility locations.
  point_list flist_;

public:
  // Iterate over facilities.
  inline iterator begin(void) const { return flist_.begin(); }
  inline iterator end(void) const { return flist_.begin(); }
  inline size_type size(void) const { return flist_.size(); }
}; // end class VPlayer

template<class Traits>
class VGame
{
  inherit_traits(Traits);
  typedef VPlayer<traits> player_type;

  // Members
private:
  // List of customer (user) points.
  point_list users_;
  player_type p1_, p2_;
  voronoi_diagram vd_;

  // Play a single round and return the ideal location for player 2.
  point_type nth_round(player_type const& p1, player_type const& p2);

public:
  VGame(point_list const& users)
    : users_(users),
  {}

  template<class InputIter>
  VGame(InputIter user_begin, InputIter user_end)
    : user_(begin, end),
  {}

  // Initialize the game with the given points for each player.
  void start(point_list& p1f, point_list& p2f)
  {
    p1_ = player_type(p1f);
    p2_ = player_type(p2f);
    // Initialize the Voronoi Diagram so we can quickly find the nearest
    // neighbor.
    auto pfacilities = facilities();
    vd_ = voronoi_diagram(
        pfacilities.begin(), pfacilities.end(),
        users_.begin(), users_.end());
  }

  // Play the game for the given number of rounds (default: play 1 round).
  // Note that one round means one player is solved (two rounds is symmetric).
  // Note also that this uses the player 2 nth-round solution in each round.
  void play_round(size_type rounds=1)
  {
    player_type* p1ptr = &p1_;
    player_type* p2ptr = &p2_;
    while (rounds--) {
      // Each round we need a new NN search since we may have added a facility.
      vd.build();
      // Solve for player 2 and add the point to p2's list and the VD sites.
      point_type p2_solution = nth_round(*p1ptr, *p2ptr);
      p2ptr->add(p2_solution);
      vd.add_site(p2_solution);
      std::swap(p1ptr, p2ptr);
    }
  }

  // Shortcut for start() and then play_round(rounds).
  template<class InputIter1, class InputIter2>
  inline void play(
      InputIter1 p1f_begin, InputIter1 p1f_end,
      InputIter2 p2f_begin, InputIter2 p2f_end,
      size_type rounds=1) {
    start(p1f_begin, p1f_end, p2f_begin, p2f_end);
    play_round(rounds);
  }

  inline player_type const& player1(void) const { return p1_; }
  inline player_type const& player2(void) const { return p2_; }
  inline auto facilities(void) const {
    return boost::join(player1(), player2());
  } // TODO auto?

}; // end class VGame

} // end namespace cfla
