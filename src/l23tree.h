//
// l23tree.h
//
// Implementation of a "leafy" 2-3 tree. This tree follows the "classical"
// definition of a 2-3 tree, where internal nodes contain only copies of the
// keys, and data nodes are only present in the leaves.
// We call it a "leafy 2-3 tree", or "l23tree", since data are only in leaves
// (and since C++ identifiers cannot start with numbers).
//

#include <utility> // pair
#include <array>
#include <list>
#include <algorithm>
#include <functional>
#include <stdexcept>

#include <cstdint>

namespace tree
{

  // Pointer helpers.
template<class Params>
struct mptr
{
  typedef typename Params::key_type* key_ptr;
  typedef typename Params::key_type& key_reference;
  typedef typename Params::data_type* data_ptr;
  typedef typename Params::data_type& data_reference;
  typedef typename Params::value_type* value_ptr;
  typedef typename Params::value_type& value_reference;
  typedef typename l23_inode<Params>::inode_type* inode_ptr;
  typedef typename l23_lnode<Params>::lnode_type* lnode_ptr;
};

template <class Params>
struct cptr
{
  typedef const typename Params::key_type* key_ptr;
  typedef const typename Params::key_type& key_reference;
  typedef const typename Params::data_type* data_ptr;
  typedef const typename Params::data_type& data_reference;
  typedef const typename Params::value_type* value_ptr;
  typedef const typename Params::value_type& value_reference;
  typedef const typename l23_inode<Params>::inode_type* inode_ptr;
  typedef const typename l23_lnode<Params>::lnode_type* lnode_ptr;
};

//! Tree parameters. You can make your own, but they should contain the members
//  that the base one here contains at least.
template<class Key, class Data, class KeyCompare,
  class Value=std::pair<Key, Data> >
struct l23_default_params
{
  typedef Key key_type;
  typedef Data data_type;
  typedef Value value_type;
  typedef KeyCompare key_compare;
  typedef uint8_t size_type;

  // Must implement for your Value class.
  static const key_type &key(const value_type& v) { return v.first; }
  static key_type &key(value_type& v) { return v.first; }
  static const data_type &data(const value_type& v) { return v.second; }
  static data_type &data(value_type& v) { return v.second; }
};

// Some forward class declarations.
template<class Params> class l23_inode;
template<class Params> class l23_lnode;
template<class Params, class Ptrs> class l23_iterator;
template<class Params> class l23_tree;

// Common node operations which operate equivalently for internal nodes with
// leaf children and internal nodes with internal children.
template<class Child_node_type, class Ptrs=mptr>
class l23_node_operations
{
public:
  // Typedefs.
  typedef typename node_type::parameters parameters;
  typedef typename parameters::key_type key_type;
  typedef typename parameters::key_type data_type;
  typedef typename parameters::value_type value_type;
  typedef typename parameters::key_type key_compare;
  typedef typename parameters::size_type size_type;
  typedef l23_inode<parameters> inode_type;
  typedef l23_lnode<parameters> lnode_type;
  typedef Child_node_type node_type;
  typedef node_type* node_ptr;
  typedef const node_type* cnode_ptr;

  typedef std::is_same<node_type, lnode_type> is_leaf;

  // Maybe-const, maybe not
  typedef typename Ptrs ptr;
  typedef typename ptr::key_ptr key_ptr;
  typedef typename ptr::key_reference key_reference;
  typedef typename ptr::data_ptr data_ptr;
  typedef typename ptr::data_reference data_reference;
  typedef typename ptr::value_ptr value_ptr;
  typedef typename ptr::value_reference value_reference;
  typedef typename ptr::inode_ptr inode_ptr;
  typedef typename ptr::lnode_ptr lnode_ptr;

  // Static helper methods.
  static inline inode_ptr toi(node_ptr p)
    { return reinterpret_cast<inode_ptr>(p); }
  static inline const inode_ptr toi(cnode_ptr p)
    { return reinterpret_cast<inode_ptr>(p); }
  static inline lnode_ptr tol(node_ptr p)
    { return reinterpret_cast<lnode_ptr>(p); }
  static inline const lnode_ptr tol(cnode_ptr p)
    { return reinterpret_cast<lnode_ptr>(p); }
  template<typename nptr>
    static inline node_ptr ton(nptr p)
      { return reinterpret_cast<node_ptr>(p); }
  template<typename nptr>
    static inline cnode_ptr ton(const nptr p)
      { return reinterpret_cast<cnode_ptr>(p); }

  static inline void set(inode_ptr parent, size_type index, lnode_ptr leaf)
    { parent->set_leaf(index, child); parent->leaves_ = true; }
  static inline void set(inode_ptr parent, size_type index, inode_ptr child)
    { parent->set_child(index, child); }

  static inline node_type child(inode_ptr parent, size_type index)
    { return ton(parent->child(index)); }

  static inline key_reference key(const lnode_ptr n) { return n->key(); }
  static inline key_reference key(const inode_ptr n) { return n->right_key(); }

  static inline void fix_keys(inode_ptr parent) {
    for (size_type kidx = 0u; kidx < parent->key_count(); ++kidx)
      keys_[kidx] = key(child(parent, kidx));
  }

  // Insert node at the given index. MUST be in sorted order,
  // or the tree will be FUBAR.
  static node_ptr insert(inode_ptr parent, node_ptr node, size_type index)
  {
    if (parent->full() || index >= inode_type::max_count())
      return nullptr;
    // Scoot over the nodes starting at index and insert the new leaf.
    // We MUST have a null pointer in the end slot which is overwritten
    // (checked by verifying we are not full() above).
    for (size_type idx = inode_type::max_count()-1; idx > index; --idx)
      set(parent, idx, child(parent, idx-1));
    set(parent, index, node);
    ++(parent->child_count_);
    // Now fix the keys. Note if the insert index was 2,
    // the caller will have to fix parent keys recursively.
    fix_keys(parent);
    return node;
  }

  static inline
  std::enable_if<is_leaf, node_ptr>::type
    insert(inode_ptr parent, size_type index, value_reference value)
  {
    if (parent->full() || index >= inode_type::max_count())
      return nullptr;
    node_type node = new node_type(value);
    return insert(parent, index, node);
  }

  // Append a leaf as the new middle or max leaf.
  // Note that if it is the new max leaf, parent keys must be fixed by the
  // caller.
  static inline node_type append(inode_ptr parent, node_ptr node)
    { return insert(parent, parent->child_count(), node); }
  static inline lnode_type* append(inode_ptr parent, value_reference value)
    { return insert(parent, parent->child_count(), value); }

  // Remove and return the current max leaf. Updates count.
  // Returns nullptr iff there are no leaves.
  // Note that if there were only two leaves, the node is now invalidated
  // (nodes should always have 2 or 3 children).
  static node_ptr pop(inode_ptr parent)
  {
    if (parent->child_count() == 0)
      return nullptr;
    node_ptr n = child(parent, --(parent->child_count_));
    if (n) n->link(nullptr);
    set(parent, parent->child_count(), nullptr);
    return n;
  }

  // Swap out a node with a new one and fix the keys.
  // You must ensure the sort invariant is maintained. Returns NULL for
  // out-of-bounds index, or if the leaf at index was already null.
  static node_ptr swap(inode_ptr parent, size_type index, node_ptr node)
  {
    if (index >= inode_type::max_count())
      return nullptr;
    node_ptr old = child(parent, index);
    set(parent, index, node);
    if (node->parent())
      set(node->parent(), index, old);
    // Unless we swapped the max key, we must fixup our internal keys.
    // If we are a max key, parent keys must be fixed by the caller.
    if (index != inode_type::max_count()-1)
      fix_keys();
    return old;
  }

};

//! 2-3 tree (i)nternal node.
template<class Params>
class l23_inode
{
public:
  // Typedefs.
  typedef Params parameters;
  typedef typename Params::key_type key_type;
  typedef typename Params::key_type data_type;
  typedef typename Params::value_type value_type;
  typedef typename Params::key_type key_compare;
  typedef typename Params::size_type size_type;
  typedef l23_inode<Params> inode_type;
  typedef l23_lnode<Params> lnode_type;
  typedef l23_node_operations<inode_type> iops;
  typedef l23_node_operations<lnode_type> lops;
  friend class l23_lnode<Params>;
  friend class l23_tree<Params>;
  friend class l23_node_operations<inode_type>;
  friend class l23_node_operations<lnode_type>;

  // Various array types. We never hold more than 3 of anything.
  typedef std::array<key_type, 2> key_array; // key copies for indexing
  typedef std::array<inode_type*, 3> node_array; // internal node children
  typedef std::array<lnode_type*, 3> leaf_array; // leaf children

public:
  // Methods.
  inline inode_type* parent(void) const { return parent_; }
  inline size_type pos(void) const { return ppos_; }
  static inline size_type max_count(void) { return 3u; }
  static inline size_type max_keys(void) { return 2u; }
  // Number of internal child nodes. This will still return the number of leaf
  // children if our children are leaves. Always 0, 2, or 3 for non-leaves.
  inline size_type child_count(void) const { return child_count_; }
  // Get internal child i. Unchecked: valid only for 0 <= i < child_count().
  inline const inode_type* child(size_type i) const { return u_.children[i]; }
  inline inode_type* child(size_type i) { return u_.children[i]; }

  inline inode_type* left_sibling(void) {
    if (parent() && pos() > 0)
      return parent()->child(pos()-1);
    return nullptr;
  }
  inline inode_type* right_sibling(void) {
    if (parent() && pos() < max_count()-1)
      return parent()->child(pos()+1);
    return nullptr;
  }

  inline inode_type* left_child(void) { return child(0); }
  inline inode_type* middle_child(void) { return child(1); }
  inline inode_type* right_child(void) { return child(2); }
  inline const inode_type* left_child(void) const { return child(0); }
  inline const inode_type* middle_child(void) const { return child(1); }
  inline const inode_type* right_child(void) const { return child(2); }

  // Whether our children are leaf nodes (true) or internal nodes (false).
  inline bool has_leaves(void) const { return leaves_; }
  // Number of leaves. 0 unless has_leaves() returns true. At most 3.
  inline size_type leaf_count(void) const
    { return has_leaves() ? child_count() : 0; }
  // Get leaf node i. Unchecked: valid only for 0 <= i < leaf_count().
  inline const lnode_type* leaf(size_type i) const { return u_.leaves[i]; }
  inline lnode_type* leaf(size_type i) { return u_.leaves[i]; }

  inline inode_type* left_leaf(void) { return leaf(0); }
  inline inode_type* middle_leaf(void) { return leaf(1); }
  inline inode_type* right_leaf(void) { return leaf(2); }
  inline const inode_type* left_leaf(void) const { return leaf(0); }
  inline const inode_type* middle_leaf(void) const { return leaf(1); }
  inline const inode_type* right_leaf(void) const { return leaf(2); }

private:
  // Set child or leaf. No bounds checking, and size is not updated.
  // However, the new node's parent pointer/pos are updated.
  inline void set_child(size_type index, inode_type* node) {
    if (node) node->link(this, index);
    u_.children[index] = node;
  }
  inline void set_leaf(size_type index, lnode_type* node) {
    if (node) node->link(this, index);
    u_.leaves[index] = node;
  }

public:
  // If we are full, we have no more room for another child.
  inline bool full(void) const { return child_count() == max_count(); }

  // Link to a new parent node.
  inline void link(inode_type* parent, size_type pos=0) {
    parent_ = parent;
    if (parent) ppos_ = pos;
    else ppos_ = -1;
  }

  // Boundary keys. All nodes in child(i) are less than key(i).
  // Generally all nodes always have exactly two keys, and either two or three
  // children. However, nodes which
  //   0 children    -> 0 keys;
  //   1, 2 children -> 1, 2 keys;
  //   3 children    -> 2 keys.
  inline size_type key_count(void) const
    { return std::min(child_count(), static_cast<size_type>(2)); }
  // Get key i. Unchecked: valid only for 0 <= i < key_count().
  inline key_type& key(size_type i) { return keys_[i]; }
  inline const key_type& key(size_type i) const { return keys_[i]; }
  inline key_type& left_key(void) { return key(0); }
  inline key_type& right_key(void) { return key(1); }
  inline const key_type& left_key(void) const { return key(0); }
  inline const key_type& right_key(void) const { return key(1); }

  // Unchecked.
  inline void set_key(size_type i, const key_type& k) { keys_[i] = k; }

public:
  // Constructors.
  l23_inode()
    : parent_(nullptr), keys_(),
      u_{.children={0}}, ppos_(-1), leaves_(false), child_count_(0)
  {}

  // Construct a leaf-container with two leaves in the given sorted order.
  l23_inode(lnode_type* left, lnode_type* middle)
    : parent_(nullptr), keys_{left->key(), middle->key()},
      u_{.leaves={left, middle}}, ppos_(-1), leaves_(true), child_count_(2)
  {
    left->link(this, 0);
    middle->link(this, 1);
  }

private:
  ////// Members.
  // Parent back-pointer. Only NULL for the root node.
  inode_type* parent_;

  // Copies of keys used to sort the children. We never have more than 2 keys,
  // since we never have more than three children.
  key_array keys_;

  // Pointers to our child nodes. Always have 0, 2, or 3 children.
  // If our child nodes are leaves, we instead store pointers to leaf nodes.
  union {
    node_array children;
    leaf_array leaves;
  } u_;

  // Position of this node in the parent.
  int ppos_;

  // If leaves_ is true, the children of this node are leaves.
  // Otherwise the children of this node are internal nodes.
  // This saves space by storing data only in the leaves, rather than
  // in internal nodes.
  bool leaves_;

  // Number of valid children (generally 2 or 3).
  size_type child_count_;
};

//! 2-3 tree (l)eaf node.
template<class Params>
class l23_lnode
{
public:
  // Typedefs.
  typedef typename Params::key_type key_type;
  typedef typename Params::key_type data_type;
  typedef typename Params::value_type value_type;
  typedef typename Params::key_type key_compare;
  typedef typename Params::size_type size_type;
  typedef l23_inode<Params> inode_type;
  typedef l23_lnode<Params> lnode_type;
  typedef l23_tree<Params> tree_type;
  friend class l23_inode<Params>;
  friend class l23_tree<Params>;

public:
  inline inode_type* parent(void) { return parent_; }
  inline const inode_type* parent(void) const { return parent_; }
  inline size_type pos(void) const { return ppos_; }

  inline value_type& value(void) { return value_; }
  inline const value_type& value(void) const { return value_; }

  inline key_type& key(void) { return Params::key(value()); }
  inline const key_type& key(void) const { return Params::key(value()); }

  inline data_type& data(void) { return Params::data(value()); }
  inline const data_type& data(void) const { return Params::data(value()); }

  inline lnode_type* left_sibling(void) {
    if (parent() && pos() > 0)
      return parent()->leaf(pos()-1);
    return nullptr;
  }
  inline lnode_type* right_sibling(void) {
    if (parent() && pos() < inode_type::max_count()-1)
      return parent()->leaf(pos()+1);
    return nullptr;
  }

  inline void link(inode_type* parent, size_type pos=0) {
    parent_ = parent;
    if (parent) ppos_ = pos;
    else ppos_ = -1;
  }

public:
  // Constructors.
  l23_lnode()
    : parent_(nullptr), value_()
  {}

  template<typename ValueType>
  l23_lnode(ValueType v)
    : parent_(nullptr), value_(v)
  {}

private:
  // Members. A leaf node just has a value and a back-pointer to its parent.
  inode_type* parent_;
  // Position in parent.
  int ppos_;
  value_type value_;
};

/*
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
*/

template<class Params, class Ptrs>
class l23_iterator
{
public:
  // Typedefs.
  typedef typename Params::key_type key_type;
  typedef typename Params::key_type data_type;
  typedef typename Params::key_type key_compare;
  typedef typename Params::size_type size_type;

  // Maybe-const, maybe not
  typedef typename Ptrs::key_ptr key_ptr;
  typedef typename Ptrs::key_reference key_reference;
  typedef typename Ptrs::data_ptr data_ptr;
  typedef typename Ptrs::data_reference data_reference;
  typedef typename Ptrs::value_ptr value_ptr;
  typedef typename Ptrs::value_reference value_reference;
  typedef typename Ptrs::inode_ptr inode_ptr;
  typedef typename Ptrs::lnode_ptr lnode_ptr;

public:
  // Constructors.

  // Iterator pointing to an internal node.
  l23_iterator(inode_ptr nptr)
    : node_(nptr), lpos_(-1) {}

  // Iterator pointing to a leaf node.
  l23_iterator(lnode_ptr lptr)
    : node_(lptr->parent()), lpos_(lptr->pos()) {}

  // Invalid iterator.
  l23_iterator()
    : node_(nullptr), lpos_(-1) {}

  // Copy constructor.
  l23_iterator(const l23_iterator& other)
    : node_(other.node_), lpos_(other.lpos_) {}

  // Whether this points to a leaf node.
  // Note that the actual leaf node pointed to may be NULL.
  inline bool is_leaf(void) const {
    return node_ && node_->has_leaves();
  }

  // Only valid when is_leaf() is true.
  inline lnode_ptr leaf(void) const { return node_->leaf(lpos_); }
  inline value_reference value(void) const { return leaf()->value(); }
  inline key_reference key(void) const { return leaf()->key(); }
  // Negative if is_leaf() is false, otherwise index of current leaf.
  inline int pos(void) const { return lpos_; }

  // Return the internal node pointed to by this iterator (if any).
  // Only valid when leaf() is false.
  inline inode_ptr node(void) { return node_; }
  // Return the parent of node().
  inline inode_ptr parent(void) { return node_->parent(); }

  // Whether this points to a valid node (false for tree.end()).
  inline operator bool(void) const { return node_ != nullptr; }

  // pre-fix
  l23_iterator& operator++(void)
  {
    // Short-circuit the null case.
    if (!*this)
      return *this;

    // Leaves move onto the next leaf until we run out.
    if (is_leaf() && ++lpos_ < node()->leaf_count())
      return *this;

    // Go to the next sibling.
    // Either we are an internal node, or we've just finished a leaf subtree.
    if ((node_ = node()->right_sibling()))
      return *this;

    // Once we run out of siblings, check the queue for the next node.
    if (q.empty())
    {
      node_ = nullptr;
      lpos_ = -1;
      return *this;
    }

    // If the next node from the queue has leaves, reset our leaf position.
    // Otherwise queue our children (only need to queue the left child, since
    // children will traverse to their right siblings as per above).
    node_ = q.front(); q.pop_front();
    if (node()->has_leaves())
      lpos_ = 0;
    else
      q.push_back(node()->left_child());
    return *this;
  }

  // post-fix
  inline l23_iterator operator++(int)
  {
    l23_iterator it(*this);
    ++(*this);
    return it;
  }

  inline inode_ptr operator*(void) const {
    if (!node()) throw std::runtime_error("l23_iterator: null dereference");
    return node();
  }

  inline bool operator==(const l23_iterator& other) const
    { return (lpos_ == other.lpos_ && node_ == other.node_); }
  inline bool operator!=(const l23_iterator& other) const
    { return !(*this == other); }

private:
  // Members.
  // Current node, or parent of current leaf node.
  inode_ptr node_;
  // Current leaf position, or -1 for internal nodes.
  int lpos_;
  // Node queue for BFS traversal.
  std::list<inode_ptr> q;
};

//! 2-3 tree with data in the leaves.
template<class Params>
class l23_tree
{
public:
  // Typedefs.
  typedef typename Params::key_type key_type;
  typedef typename Params::data_type data_type;
  typedef typename Params::value_type value_type;
  typedef typename Params::key_compare key_compare;
  typedef typename Params::size_type size_type;
  typedef l23_inode<Params> inode_type;
  typedef l23_lnode<Params> lnode_type;
  typedef l23_iterator<Params, cptr<Params> > const_iterator;
  typedef l23_iterator<Params, mptr<Params> > iterator;
  typedef l23_node_operations<inode_type> iops;
  typedef l23_node_operations<lnode_type> lops;

public:
  // Constructors.
  template<class KeyCompare>
  l23_tree(KeyCompare kcmp=key_compare())
    : root_(nullptr), compare_(kcmp)
  {}

  template<typename InputIter, typename KeyCompare=key_compare>
  l23_tree(InputIter begin, InputIter end, KeyCompare kcmp=key_compare())
    : root_(nullptr), compare_(kcmp)
  {
    while (begin != end)
      insert_unique(*begin++);
  }

  inode_type* root(void) { return root_; }
  const inode_type* root(void) const { return root_; }

  // Iterators.
  iterator begin(void) { return iterator(root()); }
  iterator end(void) { return iterator(); }
  const_iterator begin(void) const { return const_iterator(root()); }
  const_iterator end(void) const { return const_iterator(); }

  // Lookup.
public:
  template<class KeyType=const key_type&>
    const_iterator lower_bound(KeyType k) const;
  inline const_iterator lower_bound(const value_type& v) const {
    return lower_bound(Params::key(v));
  }

  template<class KeyType=const key_type&>
    const_iterator upper_bound(KeyType k) const;
  inline const_iterator upper_bound(const value_type& v) const {
    return upper_bound(Params::key(v));
  }

  template<class KeyType=const key_type&>
    const_iterator find_unique(KeyType k) const;
  inline const_iterator find_unique(const value_type& v) const {
    return find_unique(Params::key(v));
  }

  // Insert.
  iterator insert_unique(const value_type& v);

  // Remove.
  template<class KeyType=const key_type&>
    void erase_unique(KeyType k);
  inline void erase_unique(const value_type& v)
    { erase_unique(Params::key(v)); }
  void erase(iterator it);

private:
  // Simple helpers for comparing without having to remember the STL
  // comparators.
  template<class KeyType1, class KeyType2, class KeyCompare>
  inline bool compare(KeyType1 k1, KeyType2 k2, KeyCompare kcmp) const
    { return kcmp(k1, k2); }
  template<class KeyType1, class KeyType2>
  inline bool compare(KeyType1 k1, KeyType2 k2) const
    { return compare_(k1, k2); }
  template<class KeyType1, class KeyType2, class KeyCompare>
  inline bool equal(KeyType1 k1, KeyType2 k2, KeyCompare kcmp) const
    { return compare(k1, k2, kcmp) && compare(k2, k1, kcmp); }
  template<class KeyType1, class KeyType2>
  inline bool equal(KeyType1 k1, KeyType2 k2) const
    { return compare(k1, k2) && compare(k2, k1); }
  template<class KeyType1, class KeyType2, class KeyCompare>
  inline bool less_equal(KeyType1 k1, KeyType2 k2, KeyCompare kcmp) const
    { return !compare(k2, k1, kcmp); }
  template<class KeyType1, class KeyType2>
  inline bool less_equal(KeyType1 k1, KeyType2 k2) const
    { return !compare(k2, k1); }
  template<class KeyType1, class KeyType2, class KeyCompare>
  inline bool less(KeyType1 k1, KeyType2 k2, KeyCompare kcmp) const
    { return compare(k1, k2, kcmp); }
  template<class KeyType1, class KeyType2>
  inline bool less(KeyType1 k1, KeyType2 k2) const
    { return compare(k1, k2); }
  template<class KeyType1, class KeyType2, class KeyCompare>
  inline bool greater(KeyType1 k1, KeyType2 k2, KeyCompare kcmp) const
    { return compare(k2, k1); }
  template<class KeyType1, class KeyType2>
  inline bool greater(KeyType1 k1, KeyType2 k2) const
    { return compare(k2, k1); }
  template<class KeyType1, class KeyType2, class KeyCompare>
  inline bool greater_equal(KeyType1 k1, KeyType2 k2, KeyCompare kcmp) const
    { return !compare(k1, k2); }
  template<class KeyType1, class KeyType2>
  inline bool greater_equal(KeyType1 k1, KeyType2 k2) const
    { return !compare(k1, k2); }

  // Lower bound helpers.
  template<class KeyType>
    const_iterator lower_bound_checked(KeyType k) const;
  template<class KeyType>
    const_iterator lower_bound_generic(KeyType k) const;
  template<class KeyType>
    iterator lower_bound_insert(KeyType k);

  template<class KeyType>
    int child_index_generic(const inode_type* node, KeyType k) const;
  template<class KeyType>
    int child_index_insert(const inode_type* node, KeyType k) const;

  // Insert a leaf in the general case (all nodes have 2 or 3 nodes).
  iterator insert_leaf(const value_type& v);
  // Insert an internal node recursively.
  template<class Ops>
  void insert_internal(inode_type* parent, typename Ops::node_ptr new_node);

  // Fix keys along a branch from a new max key insertion.
  void fix_branch(const key_type& max_key, inode_type* pnode);

private:
  // Members.
  inode_type* root_;
  key_compare compare_;
};

// Tree implementations.
template<class Params> template<class KeyType>
typename l23_tree<Params>::const_iterator
l23_tree<Params>::lower_bound_checked(KeyType k) const
{
  // Special cases: no root or root has 1 child.
  if (!root() || !root()->child_count())
    return end();
  if (root()->child_count() == 1)
  {
    if (equal(k, root()->key(0)))
      return const_iterator(root(), 0);
    return end();
  }

  // Otherwise we can assume all nodes have 2 or 3 children.
  return lower_bound_generic(k);
}

template<class Params> template<class KeyType>
int
l23_tree<Params>::child_index_generic(const inode_type* node, KeyType k) const
{
  typename inode_type::size_type i = 0;
  for (; i < inode_type::max_count(); ++i)
    if (less_equal(k, node->key(i)))
      break;
  return i;
}

template<class Params> template<class KeyType>
int
l23_tree<Params>::child_index_insert(const inode_type* node, KeyType k) const
{
  typename inode_type::size_type i = 0;
  for (; i < node->child_count(); ++i)
    if (less(k, node->key(i)))
      break;
  return i;
}

template<class Params> template<class KeyType>
typename l23_tree<Params>::const_iterator
l23_tree<Params>::lower_bound_generic(KeyType k) const
{
  const inode_type* node = root();
  while (node && !node->has_leaves())
    node = node->child(child_index_generic(node, k));
  if (!node || !node->has_leaves())
    return end();
  // Now we must have landed on a non-null node containing 2 or 3 leaves.
  // Return an iterator pointing to the appropriate lower-bound leaf node.
  // Note that the leaf itself may not actually exist.
  return const_iterator(node->leaf(child_index_generic(node, k)));
}

template<class Params> template<class KeyType>
typename l23_tree<Params>::const_iterator
l23_tree<Params>::lower_bound(KeyType k) const
{
  const_iterator leafp = lower_bound_checked(k);
  if (leafp.leaf() == nullptr)
    return end();
  return leafp;
}

template<class Params> template<class KeyType>
typename l23_tree<Params>::const_iterator
l23_tree<Params>::upper_bound(KeyType k) const
{
  const_iterator leafp = lower_bound_checked(k);
  // Increment the iterator until we find a key value greater than k.
  while (leafp != end() && !compare(k, leafp.key()))
    ++leafp;
  return leafp; // may be end()
}

template<class Params> template<class KeyType>
typename l23_tree<Params>::const_iterator
l23_tree<Params>::find_unique(KeyType k) const
{
  const_iterator leafp = lower_bound_checked(k);
  if (leafp == end() || leafp.leaf() == nullptr)
    return end();
  // Only return the given node if its key is equal to the query key.
  // We already know compare(k, leafp.key()) is false by definition of
  // lower_bound (k NOT < leaf). If the converse is false we have equality
  // by STL convention.
  if (!compare(leafp.key(), k))
    return leafp;
  return end();
}

template<class Params> template<typename KeyType>
typename l23_tree<Params>::iterator
l23_tree<Params>::lower_bound_insert(KeyType k)
{
  inode_type* node = root();
  while (node && !node->has_leaves())
    node = node->child(child_index_insert(node, k));
  if (!node || !node->has_leaves())
    return end();
  // Now we must have landed on a non-null node containing 2 or 3 leaves.
  // Return an iterator pointing to the appropriate lower-bound leaf node.
  // Note that the leaf itself may not actually exist.
  return iterator(node->leaf(child_index_insert(node, k)));
}

template<class Params>
typename l23_tree<Params>::iterator
l23_tree<Params>::insert_unique(const value_type& value)
{
  // Special case: root doesn't exist yet, or root has only one node.
  // These cases only occur for the first two nodes.
  if (!root())
  {
    root_ = new inode_type();
    return iterator(lops::append(root_, new lnode_type(value)));
  }

  if (root()->child_count() == 1)
  {
    // Insert the node in the right spot.
    const key_type& key = Params::key(value);
    size_t pos = less(key, root()->left_leaf()->key()) ? 0 : 1;
    return iterator(lops::insert(root(), pos, value));
  }

  // Handle the general case.
  return insert_leaf(value);
}

template<class Params>
void
l23_tree<Params>::fix_branch(const key_type& max_key, inode_type* pnode)
{
  while (pnode)
  {
    if (greater(max_key, pnode->right_key()))
      pnode->set_key(1, max_key);
    pnode = pnode->parent();
  }
}

template<class Params>
typename l23_tree<Params>::iterator
l23_tree<Params>::insert_leaf(const value_type& value)
{
  // Without loss of generality, every node has either 2 or 3 children.
  // Find the proper parent node to contain the new leaf.
  iterator leafp = lower_bound_insert(Params::key(value));
  if (leafp == end())
    return end();

  return insert_internal<lops>(leafp.node(), new lnode_type(value));
}

template<class Params> template<class Ops>
void
l23_tree<Params>::insert_internal(inode_type* parent,
    typename Ops::node_ptr new_node)
{
  // If the node is not full, we can just add the node directly.
  int pos = child_index_insert(parent, Ops::key(new_node));
  if (Ops::insert(parent, pos, new_node))
  {
    // But, if we added a new max node, we must fix our ancestors' keys
    // because the new upper bound for this subtree may be higher than the
    // upper bound key from our parent.
    if (pos == inode_type::max_count() - 1)
      fix_branch(Ops::key(new_node), parent->parent());
    return iterator(new_node);
  }

  // Here the immediate parent of the new node is full, so we must create a new
  // node to hold the largest two children and push it upwards.

  // One of the two max nodes must be the current max node.
  node_ptr newnode_middle = Ops::pop(parent);
  // The other max node is either: (1) the new node,
  node_ptr newnode_left = nullptr;
  if (pos == inode_type::max_count()-1)
    newnode_left = new_leaf;
  // or (2) the current middle node, in which case the new node becomes the
  // middle node of this parent and we detach the original right two nodes.
  else
    newnode_left = Ops::swap(parent, (pos = 1), new_node);

  // Construct the new node with the two new leaves in the right order.
  if (greater(newnode_left->key(), newnode_middle->key()))
  {
    node_ptr tmp = newnode_left;
    newnode_left = newnode_middle;
    newnode_middle = tmp;
  }
  inode_ptr new_parent = new inode_type(newnode_left, newnode_middle);

  // Insert the new node into the parent recursively.
  insert_internal<iopt>(parent->parent(), new_parent);
}

} // end namespace tree
