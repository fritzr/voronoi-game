#include <iostream>
#include <algorithm>
#include <utility>
#include <vector>
#include <iomanip>
#include <list>
#include <sstream>
#include <fstream>
#include <cstdlib>

#include <boost/geometry.hpp>
#include <boost/polygon/polygon.hpp>

#include "btree/btree.h"

namespace bp = boost::polygon;
using namespace std;

static inline string d2s(bp::direction_1d which) {
  return which == bp::LOW ? "-" : "+";
}

template <typename Tp>
static inline Tp randf(Tp min=0.0, Tp max=1.0) {
  return min + (random()*((max-min) / RAND_MAX));
}

struct flt
{
  inline friend ostream& operator<<(ostream& os, const flt& f)
  {
    os.setf(ios_base::fixed, ios_base::floatfield);
    return os << setprecision(3) << setw(4) << setfill('0');
  }
};

// These are the internal nodes of the tree
struct KeyData
{
  float x;      // x coordinate of the cell boundary - actual key value
  int depth;    // cell depth for this subtree

  KeyData() : x(0.0), depth(0) {}
  KeyData(float pos) : x(pos), depth(0) {}
  KeyData(float pos, int d) : x(pos), depth(d) {}

  inline string str(void) const {
    return static_cast<ostringstream&>(ostringstream()
        << flt() << x << "(" << depth << ")").str();
  }

  inline friend ostream& operator<<(ostream& os, const KeyData& n) {
    return os << n.str();
  }
};

// These are the data nodes of the tree (leaf nodes, if you will)
struct LeafData
{
  //int rect_idx; // parent rectangle
  string name; // parent's name
  bp::direction_1d dir; // left or right edge (LOW or HIGH)?

  LeafData() : name("<null>"), dir(bp::LOW) {}
  LeafData(const string& n, bp::direction_1d which)
    : name(n), dir(which) {}

  inline string str(void) const {
    return name + d2s(dir);
  }

  inline friend ostream& operator<<(ostream& os, const LeafData& n) {
    return os << n.str();
  }
};

struct Key_compare
{
  inline bool operator()(const KeyData& n1, const KeyData& n2) const {
    return (n1.x < n2.x);
  }
};

template<typename node_type>
static inline string
data_string(const node_type* node)
{
  int i = 0;
  ostringstream nodess;
  for (; i < node->count()-1; ++i)
    nodess << node->value(i).second << " | ";
  nodess << node->value(i).second;
  return nodess.str();
}

template<typename node_type>
static inline string
key_string(const node_type* node)
{
  int i = 0;
  ostringstream nodess;
  for (; i < node->count()-1; ++i)
    nodess << node->key(i) << " | ";
  nodess << node->key(i);
  return nodess.str();
}

template<typename node_type>
static inline string
node_string(const node_type* node)
{
  ostringstream ss;
  ss << '"' << key_string(node) << "\\n" << data_string(node) << '"';
  return ss.str();
}

template <typename Tree>
ostream&
operator<<(ostream &os, const typename Tree::node_type& node)
{
  return os << node_string(node);
}

template <typename Key, typename Data, typename Compare, int TargetNodeValues>
struct make_params
{
  typedef std::allocator<Key> Allocator;
  typedef typename btree::btree_map_params<Key, Data, Compare, Allocator, 256>
    base_params;

  // Reverse the size calculation of kNodeTargetValues from kTargetNodeSize
  // to get the proper size for the given number of nodes
  typedef typename btree::btree_node<base_params> base_node;
  static const typename base_params::size_type kTargetNodeSize
    = (TargetNodeValues * base_params::kValueSize)
      + sizeof(typename base_node::base_fields);

  typedef typename btree::btree_map_params<
    Key, Data, Compare, Allocator, kTargetNodeSize
    > type;
};

static const int leaf_nodes = 3;
typedef make_params<KeyData, LeafData, Key_compare, leaf_nodes>::type
  params;

typedef btree::btree<params> Tree;

template<typename Tree_type>
struct node_visitor
{
  typedef typename Tree_type::iterator iterator;
  inline void operator()(iterator it) const {
    cout << " => " << node_string(it.node) << endl;
  }
};

template<typename Tree>
int
write_tree(const string& filename, const Tree& tree)
{
  cout << "writing " << tree.nodes() << " nodes to " << filename << endl;
  const typename Tree::node_type* root = tree.root();

  // BFS traversal
  ofstream of(filename);
  of << "digraph G {" << endl;
  list<const typename Tree::node_type*> queue;
  queue.push_front(root);
  while (!of.fail() && !queue.empty())
  {
    const typename Tree::node_type* node = queue.front(); queue.pop_front();
    if (!node->leaf())
      for (int i = 0; i < node->count()+1; ++i)
      {
        const typename Tree::node_type* child = node->child(i);
        queue.push_back(child);
        of << node_string(node) << " -> " << node_string(child) << ";" << endl;
      }
  }
  of << "}" << endl;
  of.close();

  // now run dot to generate PS file
  string dot_cmd = string("dot -Tps -O ") + filename;
  return system(dot_cmd.c_str());
}

int main(int argc, char* argv[])
{
  vector<Tree::value_type> nodes;
  srand(time(NULL));

  int num_nodes = 10;
  if (argc > 1)
    num_nodes = atoi(argv[1]);
  // bad
  if (num_nodes <= 0)
    return 1;

  // generate requested number of nodes
  const float fmax = 10.0f;
  for (int i = 0; i < num_nodes; ++i)
  {
    string name = string("x") + to_string(i);
    float lowval = randf(0.0f, fmax);
    float highval = randf(lowval, fmax);
    nodes.push_back(make_pair(
          KeyData(lowval, 0), LeafData(name, bp::LOW)));
    nodes.push_back(make_pair(
          KeyData(highval, 0), LeafData(name, bp::HIGH)));
  }

  // create the tree by inserting the nodes
  Tree t(nodes.begin(), nodes.end());
  cout << "nodes: " << endl;
  for (auto it = t.begin(); it != t.end(); ++it)
    cout << node_string(it.node) << endl;

  // test path
  float keyval = 3.33;
  cout << "path to " << flt() << keyval << ":" << endl;
  node_visitor<Tree> v;
  v(t.lower_traverse(keyval, v));

  // write nodes to a PS file
  const string filename("btree.gv");
  return write_tree(filename, t);
}
