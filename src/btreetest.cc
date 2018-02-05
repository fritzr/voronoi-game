#include <iostream>
#include <algorithm>
#include <utility>
#include <vector>
#include <iomanip>

#include <boost/geometry/geometry.hpp>

#include "btree/btree.h"

namespace bg = boost::geometry;
using namespace std;

inline static string d2s(bg::direction1d which) {
  return which == LOW ? string("LOW") : string("HIGH");
}

struct Node
{
  float x;      // x coordinate of the cell boundary - actual key value
  int depth;    // cell depth just before x - g_k[i] from the paper
  //int rect_idx; // parent rectangle
  string name; // parent's name
  bg::direction1d side; // left or right edge (LOW or HIGH)?

  Node(float pos, const string& n, bg::direction1d which)
    : x(pos), depth(0), name(n), side(which) {}

  friend ostream& operator<<(ostream& os, const Node& n) {
    return os
      << "{ " << setprecision(2) << setw(4) << setfill('0') << os.x
      << " " << name
      << " " << setw(4) << setfill(' ') << d2s(side)
      << " d=" << depth << " }";
  }
};

template <typename Key, typename Compare, int TargetNodeValues>
struct make_params
{
  typedef std::allocator<Key> Allocator;
  typedef typename btree::btree_common_params<
    Key, Compare, Allocator, 256, sizeof(Key)
    > base_params;

  // Reverse the size calculation of kNodeTargetValues from kTargetNodeSize
  // to get the proper size for the given number of nodes
  typedef typename btree::btree_node<base_params> base_node;
  typename base_params::size_type kTargetNodeSize
    = (sizeof(Key)*TargetNodeValues) + sizeof(base_node::base_fields),

  typedef typename btree::btree_common_params<
    Key, Compare, Allocator, kTargetNodeSize, sizeof(Key)
    > type;
};

static const int leaf_nodes = 3;
typedef make_params<Node, Node::compare, leaf_nodes>::type params;

typedef btree::btree<params> Tree;

int main(int argc, char* argv[])
{
  vector<Node> nodes = {
    { 1.0f, "x3", bg::LOW },
    { 1.5f, "x3", bg::HIGH },
    { 1.2f, "x2", bg::LOW },
    { 0.5f, "x1", bg::LOW },
    { 8f, "x1", bg::HIGH },
    { 3f, "x2", bg::HIGH },
    { 2f, "x4", bg::LOW },
    { 3f, "x4", bg::HIGH },
  };

  Tree t(nodes.begin(), nodes.end());

  cout << "nodes: " << endl;
  for (auto it = t.begin(); it != t.end(); ++it)
    cout << *it << endl;
  return 0;
}
