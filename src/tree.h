#include <list>
#include <iterator>

#include <boost/intrusive/avltree.hpp>
#include <boost/intrusive/rbtree.hpp>
#include <boost/intrusive/splaytree.hpp>
#include <boost/intrusive/treap.hpp>
#include <boost/intrusive/avltree_algorithms.hpp>
#include <boost/intrusive/rbtree_algorithms.hpp>
#include <boost/intrusive/splaytree_algorithms.hpp>
#include <boost/intrusive/treap_algorithms.hpp>

namespace tree {

#define _DECLARE_TRAITS(...) \
  DECLARE_TRAITS(__VA_ARGS__)
#define DECLARE_TRAITS(...) \
  typedef __VA_ARGS__ node_traits; \
  typedef typename node_traits::node node; \
  typedef typename node_traits::node_ptr node_ptr; \
  typedef typename node_traits::const_node_ptr const_node_ptr; \
  typedef tree::node_value_traits<node_traits> value_traits

#define DEFINE_NODE(...) \
  _DECLARE_TRAITS(__VA_ARGS__); \
  node_ptr parent_, left_, right_;

#define DEFINE_TRAITS(...) \
  typedef __VA_ARGS__ node; \
  typedef node* node_ptr; \
  typedef const node* const_node_ptr; \
  static inline node_ptr get_parent(const_node_ptr n) { return n->parent_; } \
  static inline void set_parent(node_ptr n, node_ptr p) { n->parent_ = p; } \
  static inline node_ptr get_left(const_node_ptr n) { return n->left_; } \
  static inline void set_left(node_ptr n, node_ptr l) { n->left_ = l; } \
  static inline node_ptr get_right(const_node_ptr n) { return n->right_; } \
  static inline void set_right(node_ptr n, node_ptr r) { n->right_ = r; }

namespace bi = boost::intrusive;

/* Generic tree structures */

template<class Nt_, bi::link_mode_type Link_t=bi::normal_link>
struct node_value_traits
{
   typedef Nt_                                                 node_traits;
   typedef typename node_traits::node_ptr                      node_ptr;
   typedef typename node_traits::const_node_ptr                const_node_ptr;
   typedef typename node_traits::node                          value_type;
   typedef node_ptr                                            pointer;
   typedef const_node_ptr                                      const_pointer;
   static const bi::link_mode_type link_mode = Link_t;
   static node_ptr to_node_ptr (value_type &value)
     {  return node_ptr(&value); }
   static const_node_ptr to_node_ptr (const value_type &value)
     {  return const_node_ptr(&value); }
   static pointer to_value_ptr(node_ptr n)
     {  return pointer(n); }
   static const_pointer to_value_ptr(const_node_ptr n)
     {  return const_pointer(n); }
};


template<class Nt_, class Node_algorithms>
class _node_iterator
  : public std::iterator<std::forward_iterator_tag, typename Nt_::node_ptr>
{
public:
  DECLARE_TRAITS(Nt_);
  typedef Node_algorithms node_algorithms;

protected:
  node_ptr header_;
  node_ptr current_;

public:
  _node_iterator(node_ptr header, node_ptr current)
    : header_(header), current_(current) {}
  inline bool operator==(const _node_iterator &o) const {
    return current_ == o.current_
      || ((current_ == header_) && (o.current_ == o.header_));
  }
  inline bool operator!=(const _node_iterator &o) const {return !(*this==o); }
  inline node_ptr operator*(void) {
    return current_;
  }
  inline node_ptr operator->(void) {
    return current_;
  }
};

template<class Nt_, class Node_algorithms>
  class _bfs_iterator : public _node_iterator<Nt_, Node_algorithms>
{
public:
  DECLARE_TRAITS(Nt_);
  typedef Node_algorithms node_algorithms;

private:
  std::list<node_ptr> q;

public:
  _bfs_iterator(node_ptr header, node_ptr current)
    : _node_iterator<node_traits, node_algorithms>(header, current)
  {
    if (current && current != header) {
      q.push_back(current);
      ++(*this);
    }
  }
  inline _bfs_iterator &operator++(void)
  {
    if (q.empty()) {
      this->current_ = this->header_;
      return *this;
    }
    this->current_ = q.front(); q.pop_front();
    node_ptr l = node_traits::get_left(this->current_);
    if (l) q.push_back(l);
    node_ptr r = node_traits::get_right(this->current_);
    if (r) q.push_back(r);
    return *this;
  }
  inline _bfs_iterator operator++(int) { // post-fix
    _bfs_iterator ret = *this;
    ++(*this);
    return ret;
  }
};

// postorder
template<class Nt_, class Node_algorithms>
  class _dfs_iterator : public _node_iterator<Nt_, Node_algorithms>
{
public:
  DECLARE_TRAITS(Nt_);
  typedef Node_algorithms node_algorithms;

private:
  enum Direction {
    Left,
    Right,
  };
  std::list<Direction> stack;

public:
  _dfs_iterator(node_ptr header, node_ptr current)
    : _node_iterator<node_traits, node_algorithms>(header, current)
  {
    if (this->current_ != this->header_)
    {
      // Push down the left side first to get to the left-deepest node.
      node_ptr left = node_traits::get_left(this->current_);
      while (left)
      {
        stack.push_back(Left);
        this->current_ = left;
        left = node_traits::get_left(this->current_);
      }
    }
  }

  node_ptr push_branch (node_ptr n)
  {
    node_ptr last = n;
    if (n)
    {
      node_ptr left = node_traits::get_left(last);
      while (left)
      {
        last = left;
        stack.push_back(Left);
        left = node_traits::get_left(last);
      }
    }
    return last;
  }

  inline _dfs_iterator &operator++(void)
  {
    if (stack.empty())
    {
      this->current_ = this->header_;
      return *this;
    }

    // We are always ready to back up to the parent.
    this->current_ = node_traits::get_parent(this->current_);
    Direction lastdir = stack.back(); stack.pop_back();

    // If we went left to get to this leaf, now traverse the right sub-tree
    // (starting from its left-most leaf).
    if (lastdir == Left)
    {
      node_ptr right = node_traits::get_right(this->current_);
      if (right)
      {
        stack.push_back(Right);
        this->current_ = push_branch(right);
      }
    }
    return *this;
  }

  inline _dfs_iterator operator++(int) { // post-fix
    _dfs_iterator ret = *this;
    ++(*this);
    return ret;
  }
};

template<class Nt_, class Node_algorithms>
  class _inorder_iterator : public _node_iterator<Nt_, Node_algorithms>
{
public:
  DECLARE_TRAITS(Nt_);
  typedef Node_algorithms node_algorithms;

  _inorder_iterator(node_ptr header, node_ptr current)
    : _node_iterator<node_traits, node_algorithms>(header, current) { }

  inline _inorder_iterator &operator++(void) {
    if (this->current_ != this->header_)
      this->current_ = node_algorithms::next_node(this->current_);
    return *this;
  }

  inline _inorder_iterator &operator--(void) {
    if (this->current_ != this->header_)
      this->current_ = node_algorithms::prev_node(this->current_);
    return *this;
  }

  inline _inorder_iterator operator++(int) { // post-fix
    _inorder_iterator ret = *this;
    ++(*this);
    return ret;
  }
};

// Like the boost::intrusive::*tree, but templated on node type,
// providing wrapper methods to invoke the boost::intrusive::*tree_algorithms
// capabilities which operate on nodes.
// This is actually not needed with proper manipulation of value_traits.
template<class Node_traits, class Node_compare, class Node_algorithms>
class generic_tree
{
public:
  DECLARE_TRAITS(Node_traits);
  typedef Node_compare node_compare;
  typedef Node_algorithms node_algorithms;

  // For iteration.
  typedef node_ptr value_type;
  typedef size_t size_type;
  typedef _inorder_iterator<node_traits, node_algorithms> inorder_iterator;
  typedef _dfs_iterator<node_traits, node_algorithms> dfs_iterator;
  typedef _bfs_iterator<node_traits, node_algorithms> bfs_iterator;
  typedef inorder_iterator iterator;

private:
  // Header, needed for the algorithms methods.
  node_ptr header_;
  node_compare cmp_;

  inline void init(void) {
    if (header_ == nullptr) {
      header_ = new node;
    }
    node_algorithms::init_header(header_);
  }

public:

  ~generic_tree() {
    if (header_ != nullptr) {
      delete header_;
      header_ = nullptr;
    }
  }

  // Initialize a default tree. Uses the default constructor for node_compare.
  generic_tree(void) : header_(nullptr), cmp_() { init(); }
  // Initialize a tree with a different key comparator.
  generic_tree(const node_compare &cmp)
    : header_(nullptr), cmp_(cmp) { init(); }

  // Node traversal (all iterators share a common 'end' iterator).
  // In-order traversal.
  inline iterator begin (void) {
    return inorder_iterator(header_, node_algorithms::begin_node(header_));
  }
  inline iterator end (void) {
    return inorder_iterator(header_, header_);
  }
  // Breadth-first traversal.
  inline bfs_iterator bfs_begin(void) {
    return bfs_iterator(header_, root());
  }
  inline bfs_iterator bfs_end(void) {
    return bfs_iterator(header_, header_);
  }
  // Post-order traversal.
  dfs_iterator dfs_begin(void) {
    return dfs_iterator(header_, root());
  }
  dfs_iterator dfs_end(void) {
    return dfs_iterator(header_, header_);
  }

  // Node methods.
  inline node_ptr header(void) const { return header_; }
  inline node_ptr root(void) const { return node_traits::get_parent(header_); }
  inline node_ptr left(void) const { return node_traits::get_left(header_); }
  inline node_ptr right(void) const { return node_traits::get_right(header_); }

  static inline node_ptr parent(node_ptr n) {
    return node_traits::get_parent(n);
  }
  static inline node_ptr left(node_ptr n) { return node_traits::get_left(n); }
  static inline node_ptr right(node_ptr n) { return node_traits::get_right(n); }
  static inline node_ptr next(node_ptr n) {
    return node_algorithms::next_node(n);
  }
  static inline node_ptr prev(node_ptr n) {
    return node_algorithms::prev_node(n);
  }

  // Search methods.
  template<typename KeyType> inline
    node_ptr lower_bound(const KeyType &k) {
      return node_algorithms::lower_bound(header_, k, cmp_);
    }
  template<typename KeyType, typename KeyNodePtrCompare> inline
    node_ptr lower_bound(const KeyType &k, KeyNodePtrCompare c) {
      return node_algorithms::lower_bound(header_, k, c);
    }
  template<typename KeyType> inline
    node_ptr upper_bound(const KeyType &k) {
      return node_algorithms::upper_bound(header_, k, cmp_);
    }
  template<typename KeyType, typename KeyNodePtrCompare> inline
    node_ptr upper_bound(const KeyType &k, KeyNodePtrCompare c) {
      return node_algorithms::upper_bound(header_, k, c);
    }
  template<typename KeyType> inline
    node_ptr find(const KeyType &k) {
      return node_algorithms::find(header_, k, cmp_);
    }
  template<typename KeyType, typename KeyNodePtrCompare> inline
    node_ptr find(const KeyType &k, KeyNodePtrCompare c) {
      return node_algorithms::find(header_, k, c);
    }
  template<typename KeyType> inline
    node_ptr equal_range(const KeyType &k) {
      return node_algorithms::equal_range(header_, k, cmp_);
    }
  template<typename KeyType, typename KeyNodePtrCompare> inline
    node_ptr equal_range(const KeyType &k, KeyNodePtrCompare c) {
      return node_algorithms::equal_range(header_, k, c);
    }
  template<typename KeyType> inline
    node_ptr count(const KeyType &k) {
      return node_algorithms::count(header_, k, cmp_);
    }
  template<typename KeyType, typename KeyNodePtrCompare> inline
    node_ptr count(const KeyType &k, KeyNodePtrCompare c) {
      return node_algorithms::count(header_, k, c);
    }
  template<typename KeyType>
    inline std::pair< node_ptr, node_ptr >
    bounded_range(const KeyType &kl, const KeyType &ku, bool cl, bool cu) {
      return node_algorithms::bounded_range(header_, kl, ku, cmp_, cl, cu);
    }
  template<typename KeyType, typename KeyNodePtrCompare>
    inline std::pair< node_ptr, node_ptr >
    bounded_range(const KeyType &kl, const KeyType &ku, KeyNodePtrCompare cmp,
        bool cl, bool cu) {
      return node_algorithms::bounded_range(header_, kl, ku, cmp, cl, cu);
    }

  // Insertion methods.
  inline node_ptr insert_equal_upper_bound(node_ptr const & n) {
    return node_algorithms::insert_equal_upper_bound(header_, n, cmp_);
  }
  template <typename NodeCompare>
  inline node_ptr insert_equal_upper_bound(node_ptr const & n, NodeCompare cmp) {
    return node_algorithms::insert_equal_upper_bound(header_, n, cmp);
  }
  inline node_ptr insert_equal_lower_bound(node_ptr const & n) {
    return node_algorithms::insert_equal_lower_bound(header_, n, cmp_);
  }
  template <typename NodeCompare>
  inline node_ptr insert_equal_lower_bound(node_ptr const & n, NodeCompare cmp) {
    return node_algorithms::insert_equal_lower_bound(header_, n, cmp);
  }
  inline node_ptr insert_equal(node_ptr const & n) {
    return node_algorithms::insert_equal(header_, parent(header_), n, cmp_);
  }
  inline node_ptr insert_equal(node_ptr const & hint, node_ptr const & n) {
    return node_algorithms::insert_equal(header_, hint, n, cmp_);
  }
  template <typename C>
  inline node_ptr insert_equal(node_ptr const & n, C cmp) {
    return node_algorithms::insert_equal(header_, parent(header_), n, cmp);
  }
  template <typename C>
  inline node_ptr insert_equal(node_ptr const & hint, node_ptr const & n, C cmp) {
    return node_algorithms::insert_equal(header_, hint, n, cmp);
  }

  // Manual insertions (without checking preconditions).
  inline node_ptr insert_before(node_ptr const & pos, node_ptr const & n) {
    return node_algorithms::insert_before(header_, pos, n);
  }
  inline void push_back(node_ptr const & n) {
    return node_algorithms::push_back(header_, n);
  }
  inline void push_front(node_ptr const & n) {
    return node_algorithms::push_front(header_, n);
  }

  // Removal methods.
  inline void unlink(node_ptr const & n) {
    return node_algorithms::erase(header_, n);
  }
  inline void erase(node_ptr const & n) {
    return node_algorithms::erase(header_, n);
  }

  // Property methods.
  inline size_type size(void) const {
    return node_algorithms::size(root()); // size of subtree at root (all)
  }
  inline static size_type size(node_ptr const & n) {
    return node_algorithms::size(n); // size of subtree at n
  }
  inline size_type depth(node_ptr const & n) const {
    // Depth of a subtree (distance from root).
    int depth = 0;
    node_ptr cur = n;
    while (cur != root()) {
      ++depth;
      cur = node_traits::get_parent(cur);
    }
    return depth;
  }
  inline static bool is_leaf(const_node_ptr n) {
    return (node_traits::get_left(n) == nullptr)
      && (node_traits::get_right(n) == nullptr);
  }
};

/* AVL_TREE */

enum default_balance {
  negative_t, zero_t, positive_t
};

template<class Node, class Balance=default_balance>
struct avl_node_traits
{
  DEFINE_TRAITS(Node);
  typedef Balance balance;
  static inline balance get_balance(const_node_ptr n) { return n->balance_; }
  static inline void set_balance(node_ptr n, balance b) { n->balance_ = b; }
  static inline balance negative(void) { return balance(-1); }
  static inline balance zero(void) { return balance(0); }
  static inline balance positive(void) { return balance(1); }
};

// Partial specialization
template<class Node>
struct avl_node_traits<Node, default_balance>
{
  DEFINE_TRAITS(Node);
  typedef default_balance balance;
  static inline balance get_balance(const_node_ptr n) { return n->balance_; }
  static inline void set_balance(node_ptr n, balance b) { n->balance_ = b; }
  static inline balance negative(void) { return balance::negative_t; }
  static inline balance zero(void) { return balance::zero_t; }
  static inline balance positive(void) { return balance::positive_t; }
};

template<class Derived, class Balance=default_balance>
class avl_node
{
public:
  DEFINE_NODE(avl_node_traits<Derived, Balance>);
  typename node_traits::balance balance_;
};

template<class Node, class Node_compare, class Balance=default_balance>
struct avl_tree
{
  typedef avl_node_traits<Node, Balance> node_traits;
  typedef typename bi::avltree_algorithms<node_traits> node_algorithms;
  typedef generic_tree<node_traits, Node_compare, node_algorithms> type;
};

template<class Node>
struct make_avltree
{
  typedef typename bi::value_traits<typename Node::value_traits> value_option;
  typedef typename bi::compare<typename Node::compare> compare_option;
  typedef typename bi::constant_time_size<false> size_option;
  typedef typename bi::avltree<Node, value_option, compare_option, size_option>
    type;
};


/* RB_TREE */

enum default_color {
  red_t, black_t
};

template<class Node, class Color=default_color>
struct rb_node_traits
{
  DEFINE_TRAITS(Node);
  typedef Color color;
  static inline color get_color(const_node_ptr n) { return n->color_; }
  static inline void set_color(node_ptr n, color c) { n->color_ = c; }
  static color black (void);
  static color red (void);
};

// Partial specialization
template<class Node>
struct rb_node_traits<Node, default_color>
{
  DEFINE_TRAITS(Node);
  typedef default_color color;
  static inline color get_color(const_node_ptr n) { return n->color_; }
  static inline void set_color(node_ptr n, color c) { n->color_ = c; }
  color black(void) { return color::black_t; }
  color red(void) { return color::red_t; }
};

template<class Derived, class Color=default_color>
class rb_node
{
public:
  DEFINE_NODE(rb_node_traits<Derived, Color>);
  typename node_traits::color color_;
};

template<class Node, class Node_compare, class Color=default_color>
struct rb_tree
{
  typedef rb_node_traits<Node, Color> node_traits;
  typedef typename bi::rbtree_algorithms<node_traits> node_algorithms;
  typedef generic_tree<node_traits, Node_compare, node_algorithms> type;
};

template<class Node>
struct make_rbtree
{
  typedef typename bi::value_traits<typename Node::value_traits> value_option;
  typedef typename bi::compare<typename Node::compare> compare_option;
  typedef typename bi::constant_time_size<false> size_option;
  typedef typename bi::rbtree<Node, value_option, compare_option, size_option>
    type;
};


/* SPLAY_TREE */

template<class Node>
struct splay_node_traits
{
  DEFINE_TRAITS(Node);
};

template<class Derived>
class splay_node
{
public:
  DEFINE_NODE(splay_node_traits<Derived>);
};

template<class Node, class Node_compare>
struct splay_tree
{
  typedef splay_node_traits<Node> node_traits;
  typedef typename bi::splaytree_algorithms<node_traits> node_algorithms;
  typedef generic_tree<node_traits, Node_compare, node_algorithms> type;
};

template<class Node>
struct make_splaytree
{
  typedef typename bi::value_traits<typename Node::value_traits> value_option;
  typedef typename bi::compare<typename Node::compare> compare_option;
  typedef typename bi::constant_time_size<false> size_option;
  typedef typename bi::splaytree<Node, value_option, compare_option, size_option>
    type;
};


/* TREAP_TREE */

template<class Node>
struct treap_node_traits
{
  DEFINE_TRAITS(Node);
};

template<class Derived>
class treap_node
{
public:
  DEFINE_NODE(treap_node_traits<Derived>);
};

template<class Node, class Node_compare>
struct treap
{
  typedef treap_node_traits<Node> node_traits;
  typedef typename bi::treap_algorithms<node_traits> node_algorithms;
  typedef generic_tree<node_traits, Node_compare, node_algorithms> type;
};

template<class Node>
struct make_treap
{
  typedef typename bi::value_traits<typename Node::value_traits> value_option;
  typedef typename bi::compare<typename Node::compare> compare_option;
  typedef typename bi::constant_time_size<false> size_option;
  typedef typename bi::treap<Node, value_option, compare_option, size_option>
    type;
};

}; // namespace tree
