#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>

#include "tree.h"

using namespace std;

// Node comparison function.
template <class Traits>
struct node_compare {
  DECLARE_TRAITS(Traits);
  bool operator()(const_node_ptr a, const_node_ptr b) {
    return a->value_ < b->value_;
  }
  bool operator()(const_node_ptr b, int val) {
    return b->value_ < val;
  }
  bool operator()(int val, const_node_ptr b) {
    return val < b->value_;
  }
};

// internal nodes
class Node : public tree::avl_node<Node>
{
public:
  typedef node_compare<typename tree::avl_node<Node>::node_traits> compare;
  // Required members are inherited - here are extra members
  int value_;
  string name;

  Node() : value_(-1), name("header") {}
  Node(int val, const string &n) : value_(val), name(n) {}

  string str(void) const {
    ostringstream ss;
    ss << dec << setw(2) << setfill(' ') << value_ << " " << name;
    return ss.str();
  }
};

typedef typename tree::avl_tree<Node, typename Node::compare>::type AvlTree;
typedef typename AvlTree::node_algorithms avl;

int main()
{
  vector<Node> node_storage = {
    { 0, "hi" },
    { 1, "x3-" },
    { 1, "x1-" },
    { 7, "x3+" },
    { 1, "x2-" },
    { 3, "x1+" },
    { 2, "x2+" },
  };

  // Initialize our AVL tree from the nodes above
  AvlTree tree;
  for (auto it = node_storage.begin(); it != node_storage.end(); ++it) {
    tree.insert_equal_upper_bound(&(*it));
  }

  // Tree diagram - traverse breadth-first
  cout << "tree:" << endl;
  int max_height = tree.height();
  int last_height = max_height;
  for (auto it = tree.bfs_begin(); it != tree.bfs_end(); ++it) {
    auto node = *it;
    int height = tree.height(node);
    // When we move to a new row, end the line
    if (height != last_height) {
      cout << endl;
      last_height = height;
    }
    cout << "  " << (node ? node->str() : "<null>");
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
  cout << "bounded (" << bmin << ", " << bmax << "):" << endl;
  auto range = tree.bounded_range(bmin, bmax, true, true);
  for (auto it = range.first; it != range.second; it = tree.next(it)) {
    cout << it->str() << "  ";
  }
  cout << endl;

  return 0;
}
