#include <iostream>
#include <algorithm>
#include <utility>
#include <vector>
#include <iomanip>
#include <list>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cstdio>

#include <boost/geometry.hpp>
#include <boost/polygon/polygon.hpp>

#include "l23tree.h"

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
template<typename Tp_>
struct KeyData
{
  typedef std::less<KeyData<Tp_> > compare;

  typedef Tp_ coordinate_type;
  coordinate_type x; // x coordinate of the cell boundary - actual key value
  int depth;         // cell depth for this subtree

  KeyData() : x(0), depth(0) {}
  KeyData(coordinate_type pos) : x(pos), depth(0) {}
  KeyData(coordinate_type pos, int d) : x(pos), depth(d) {}

  inline string str(void) const {
    return static_cast<ostringstream&>(ostringstream()
        << flt() << x << "(" << depth << ")").str();
  }

  inline bool operator<(const KeyData& o) const { return x < o.x; }

  inline friend ostream& operator<<(ostream& os, const KeyData& n) {
    return os << n.str();
  }
};

// These are the data nodes of the tree (leaf nodes, if you will)
struct LeafData
{
  //int rect_idx; // parent rectangle
  int id; // parent's name
  bp::direction_1d dir; // left or right edge (LOW or HIGH)?

  LeafData() : id(-1), dir(bp::LOW) {}
  LeafData(int index, bp::direction_1d which)
    : id(index), dir(which) {}

  inline string name(void) const {
    if (id >= 0)
      return string("x") + to_string(id);
    return string("<bad>");
  }

  inline string str(void) const {
    return name() + d2s(dir);
  }

  inline friend ostream& operator<<(ostream& os, const LeafData& n) {
    return os << n.str();
  }
};

template<typename Tp_>
struct make_tree_params
{
  typedef LeafData data_type;
  typedef KeyData<Tp_> key_type;
  typedef typename key_type::compare compare_type;

  typedef typename tree::l23_default_params<key_type, data_type, compare_type>
    type;
};

typedef make_tree_params<float>::type TreeParams;
typedef tree::l23_tree<TreeParams> Tree;


/*
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

template<typename Tree>
int
write_tree(const string& base_filename, const Tree& tree)
{
  string gv_filename = base_filename + ".gv";
  string ps_filename = base_filename + ".ps";
  cout << "writing " << tree.nodes() << " nodes to " << ps_filename << endl;
  const typename Tree::node_type* root = tree.root();

  // BFS traversal
  ofstream of(gv_filename);
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

  // now run dot to generate PS file; cleanup the gv file on success
  string dot_cmd = string("dot -Tps ") + gv_filename
    + string(" > ") + ps_filename;
  int ret = system(dot_cmd.c_str());
  if (ret == 0)
    remove(gv_filename.c_str());
  return ret;
}

template<typename Tree>
static void
init_depths(Tree& t)
{
  typedef typename Tree::data_type data_type;
  typedef typename Tree::key_type key_type;
  typedef typename Tree::iterator iterator;
  int depth = 0;
  set<int> current;
  for (iterator it = t.begin(); it != t.end(); ++it)
  {
    key_type& key = const_cast<key_type&>(it->first);
    const data_type& data = it->second;
    auto curid = current.find(data.id);
    if (curid == current.end())
    {
      current.insert(data.id);
      key.depth = ++depth;
    }
    else
    {
      current.erase(curid);
      key.depth = depth--;
    }
  }
}
*/

int main(int argc, char* argv[])
{
  vector<Tree::value_type> nodes;

  int num_nodes = 10;
  if (argc > 1)
    num_nodes = atoi(argv[1]);
  // bad
  if (num_nodes <= 0)
    return 1;

  if (argc > 2)
    srand(atoi(argv[2]));

  // generate requested number of nodes
  const float fmax = 10.0f;
  for (int i = 0; i < num_nodes; ++i)
  {
    float lowval = randf(0.0f, fmax);
    float highval = randf(lowval, fmax);
    nodes.push_back(make_pair(
          Tree::key_type(lowval, 0), Tree::data_type(i, bp::LOW)));
    nodes.push_back(make_pair(
          Tree::key_type(highval, 0), Tree::data_type(i, bp::HIGH)));
  }

  // create the tree by inserting the nodes
  Tree t(nodes.begin(), nodes.end());
  cout << "nodes: " << endl;
  //for (auto it = t.begin(); it != t.end(); ++it)
    //cout << node_string(*it) << endl;

  // test path
  /*
  float keyval = 3.33;
  cout << "path to " << flt() << keyval << ":" << endl;
  node_visitor<Tree> v;
  v(t.lower_traverse(keyval, v));

  if (argc > 3)
  {
    size_t rcount = atoi(argv[3]);
    while (rcount--)
      t.erase_unique(t.root()->key(0));
  }

  // write nodes to a PS file
  return write_tree("btree", t);
  */
}
