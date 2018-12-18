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
#include <boost/geometry/algorithms/comparable_distance.hpp>

namespace cfla
{

// Generic typedefs common for cfla stuff to ease the templating
// of the many classes which require these typedefs.
//
//
// MaxDepthSolver solves the max-depth problem (P2 solution for the Voronoi
// Game). It should be constructible from an iterator over user Points. Then it
// should define `Point operator()(nn1_type, point_filter)`.
// The arguments are an `nn1_type` reference (see below) and a unary filter
// function object which accepts user points and returns whether or not the
// solver should consider them in the solution.
// The returned Point is the P2 solution.
// That is to say the point, when placed among the facilities, has a maximal
// number of users closer to it than any other facility.
//
// The solver type must provide a user_point() method to return the
// location of a user_type object, unless user_type and point_type are
// equivalent.
//
// Additionally, MaxDepthSolver must expose a typedef `nn1_type` used for
// finding the nearest point to a set of points. Some common acceptable types
// are the L(n)NN1 types in `nn1.h`.
//
// `nn1_type` should have the following interface:
//
// `nn1_type` should be constructible from an iterator over facility Points,
// and should define `point_type operator()(user_type)` to return the nearest
// facility point to a given user point. The user point need not be, and indeed
// usually is not, among the set of facility points.
//
// `nn1_type` should also be a container-like type which has begin(), end(),
// size(), and insert() methods that involve the facility points.

template<class Point
  , class MaxDepthSolver
  , class UserType=Point
  , unsigned int NumPlayers=2u
>
struct cfla_traits
{
  static const unsigned int nplayers = NumPlayers;
  typedef MaxDepthSolver                       solver_type;
  typedef typename solver_type::const_iterator user_iterator;

  typedef typename solver_type::nn1_type       nn1_type;
  typedef typename nn1_type::const_iterator    facility_iterator;

  typedef UserType                             user_type;
  typedef Point                                point_type;
  typedef typename bg::coordinate_type<point_type>::type coordinate_type;
  typedef size_t                               size_type;
  /*
  static const unsigned int norm     = DistanceNorm;
  typedef PointContainer                      point_list;
  typedef typename point_list::iterator       point_iterator;
  typedef typename point_list::const_iterator point_citerator;
  */
};

#define inherit_traits(t) \
public: \
  typedef t traits; \
  static const auto nplayers = t::nplayers; \
  typedef typename traits::user_iterator     user_iterator; \
  typedef typename traits::facility_iterator facility_iterator; \
  typedef typename traits::solver_type       solver_type; \
  typedef typename traits::nn1_type          nn1_type; \
  typedef typename traits::user_type         user_type; \
  typedef typename traits::point_type        point_type; \
  typedef typename traits::coordinate_type   coordinate_type; \
  typedef typename traits::size_type         size_type; \

  /*
  static const auto norm     = t::norm; \
  typedef typename traits::point_list      point_list; \
  typedef typename traits::point_iterator  point_iterator; \
  typedef typename traits::point_citerator point_citerator; \
  */

template<class Traits>
class VPlayer
{
  inherit_traits(Traits);

  // Members
private:
  // Player ID.
  const int id_;
  // Number of users we 'own' (closest facility is ours).
  int score_;
  // Used for locating the nearest facility to each user point.
  nn1_type nn_;

public:
  template<class FacilityIter>
  VPlayer(int id, FacilityIter fbegin, FacilityIter fend)
    : id_(id), score_(-1), nn_(fbegin, fend)
  {}

  inline int id(void) const { return id_; }

  // Iterate over facilities.
  inline facility_iterator begin(void) const { return nn_.begin(); }
  inline facility_iterator end(void) const { return nn_.end(); }
  inline size_type size(void) const { return nn_.size(); }
  inline void insert(const point_type &facility) { nn_.insert(facility); }

  inline void set_score(int new_score) { score_ = new_score; }
  inline int score(void) const { return score_; }
  inline void add_score(int to_add=1) { set_score(score()+to_add); }

  inline const nn1_type &nn1(void) const { return nn_; }

  inline point_type nearest_facility(user_type const& user) const {
    return nn_(user);
  }

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
  solver_type solver_;
  unsigned int current_round;
#ifdef DEBUG
  cv::Mat* img_ = nullptr;
#endif

public:
  template<class UserIter>
  VGame(UserIter users_begin, UserIter users_end)
    : players_({}), solver_(users_begin, users_end), current_round(0)
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

  inline user_iterator users_begin(void) const { return solver_.begin(); }
  inline user_iterator users_end(void) const { return solver_.end(); }
  inline size_type users_size(void) const { return solver_.size(); }
  inline user_type& user(int user_idx) const {
    return *std::advance(users_begin(), user_idx);
  }

  inline typename std::enable_if<
      !std::is_same<point_type, user_type>::value
    , point_type>::type
  user_point(user_type const&u) const {
    return solver_.user_point(u);
  }

  // Initialize a player with his facilities.
  template<class InputIter>
  inline void init_player(unsigned int id, InputIter fbegin, InputIter fend) {
    if (players_[id])
      delete players_[id];
    players_[id] = new player_type(id, fbegin, fend);
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

  // To find who owns a point, get each point's closest facility and narrow
  // down to which one is really closer. This prevents us from having to
  // keep an extra range tree (since we already have one per player).
  const player_type& owner(user_type const& user) const
  {
    coordinate_type dist = std::numeric_limits<coordinate_type>::infinity();
    const player_type* min_player = &player(0);
    for (auto playerp = players_.begin(); playerp != players_.end(); ++playerp)
    {
      point_type facility = (*playerp)->nearest_facility(user);
      coordinate_type d = solver_type::distance(user, facility);
      if (std::abs(d) < std::abs(dist)) {
        dist = d;
        min_player = *playerp;
      }
    }
    return *min_player;
  }

  // Modifiable reference.
  inline player_type& owner(user_type const& user) {
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
  const int winner_id(void) const
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

  // Return the winning player. You must run score() first.
  inline const player_type& winner(void) const { return player(winner_id()); }
  inline player_type& winner(void) { return player(winner_id()); }

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

  typedef std::unary_function<user_type, bool> point_filter;
  struct filter_users : public point_filter
  {
    const VGame &game;
    const int player_id;
    filter_users(const VGame &g, int id)
      : game(g), player_id(id)
    {}

    inline bool operator()(user_type const& p) const {
      return game.owner(p).id() == player_id;
    }
  };

  // Play a single round and return the ideal location for player 2.
  // Only visit user points which are owned by p1.
  inline point_type nth_round(player_type const& p1, player_type const& p2)
  {
    return solver_(p1.nn1(), filter_users(*this, p1.id()));
  }

}; // end class VGame

#undef inherit_traits

} // end namespace cfla
