#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>

#include <boost/intrusive/avltree.hpp>
#include "tree.h"

namespace bi = boost::intrusive;
using namespace std;

// Node comparison function.
template <class node>
struct node_compare {
  inline bool operator()(const node& a, const node& b) const
    { return a.value_ < b.value_; }
  inline bool operator()(const node& b, float val) const
    { return b.value_ < val; }
  inline bool operator()(float val, const node& b) const
    { return val < b.value_; }
};

// internal nodes
class Node : public tree::avl_node<Node>
{
public:
  typedef node_compare<Node> compare;

  // Required members are inherited - here are extra members
  float value_;
  string name;

  Node() : value_(-1.0f), name("header") {}
  Node(float val, const string &n) : value_(val), name(n) {}

  static string nstr(const Node *n) {
    ostringstream ss;
    if (n) {
      ss << setw(6) << setprecision(4) << n->value_ << " "
        << setw(6) << n->name;
    } else {
      ss << "<null>";
    }
    return ss.str();
  }
  string str(void) const {
    ostringstream ss;
    ss << nstr(this) << " { " << nstr(parent_) << " | <"
      << nstr(left_) << " , " << nstr(right_) << "> }";
    return ss.str();
  }
};

typedef typename tree::make_avltree<Node>::type AvlTree;

int main()
{
  vector<Node> node_storage = {
    { 1.0f, "x3-" },
    { 3.0f, "x1-" },
    { 9.0f, "x3+" },
    { 1.5f, "x2-" },
    { 4.2f, "x1+" },
    { 2.3f, "x2+" },
  };

  // Initialize our AVL tree from the nodes above
  AvlTree tree;
  const size_t num_nodes = node_storage.size();
  for (size_t idx = 0; idx < num_nodes; ++idx) {
    tree.insert_equal(node_storage[idx]);
  }

  for (auto it = tree.begin(); it != tree.end(); ++it)
  {
    cout << it->str() << endl;
  }

  /*
  // Tree diagram - traverse breadth-first
  cout << "tree:" << endl;
  int last_depth = 0;
  for (auto it = tree.bfs_begin(); it != tree.bfs_end(); ++it) {
    auto node = *it;
    int depth = tree.depth(node);
    // When we move to a new row, end the line
    if (depth != last_depth) {
      cout << endl;
      last_depth = depth;
    }
    cout << "  " << setw(10) << setfill('*') << (node ? node->str() : "<null>");
  }
  cout << endl;

  // Sorted nodes
  cout << "inorder:" << endl;
  for (auto it = tree.begin(); it != tree.end(); ++it) {
    cout << it->str() << "  ";
  }
  cout << endl;

  // Traverse depth-first (postorder)
  cout << "dfs:" << endl;
  for (auto it = tree.dfs_begin(); it != tree.dfs_end(); ++it) {
    cout << (*it ? (*it)->str() : "<null>") << "  ";
  }
  cout << endl;

  // Bounded range
  const int bmin=1, bmax=2;
  cout << "bounded [" << bmin << ", " << bmax << "]:" << endl;
  auto range = tree.bounded_range(bmin, bmax, true, true);
  for (auto it = range.first; it != range.second; it = tree.next(it)) {
    cout << it->str() << "  ";
  }
  cout << endl;

  */
  return 0;
}
