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

namespace tree
{

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
  typedef size_t size_type;

  // Must implement for your Value class.
  static const key_type &key(const value_type& v) { return v.first; }
  static key_type &key(value_type& v) { return v.first; }
  static const data_type &data(const value_type& v) { return v.second; }
  static data_type &data(value_type& v) { return v.second; }
};

// Some forward class declarations.
template<class Params> class l23_inode;
template<class Params> class l23_lnode;
template<class Params> class l23_rootnode;
template<class Params> class l23_iterator;
template<class Params> class l23_tree;

//! 2-3 tree (i)nternal node.
template<class Params>
class l23_inode
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
  friend class lnode_type;
  friend class tree_type;

  // Various array types. We never hold more than 3 of anything.
  typedef std::array<key_type, 2> key_array; // key copies for indexing
  typedef std::array<inode_type*, 3> node_array; // internal node children
  typedef std::array<lnode_type*, 3> leaf_array; // leaf children

public:
  // Methods.
  inline inode_type* parent(void) const { return parent_; }
  static inline size_type max_count(void) { return 3; }
  // Number of internal child nodes. This will still return the number of leaf
  // children if our children are leaves. Always 0, 2, or 3 for non-leaves.
  inline size_type child_count(void) const { return child_count_; }
  // Get internal child i. Unchecked: valid only for 0 <= i < child_count().
  inline const inode_type* child(size_type i) const { return u_.children[i]; }
  inline inode_type* child(size_type i) { return u_.children[i]; }
  inline bool add_child(inode_type *new_child)
  {
    if (child_count() < max_count())
    {
      u_.children[(child_count_++ - 1)] = new_child;
      return true;
    }
    return false;
  }
  inline inode_type* set_child(size_type i, inode_type* new_child)
  {
    inode_type* old_child = u_.children[i];
    u_.children[i] = new_child;
    child_count_ = std::max(child_count_, i);
    return old_child;
  }

  inline inode_type* left_child(void) { return child(0); }
  inline inode_type* middle_child(void) { return child(1); }
  inline inode_type* right_child(void) { return child(2); }
  inline const inode_type* left_child(void) const { return child(0); }
  inline const inode_type* middle_child(void) const { return child(1); }
  inline const inode_type* right_child(void) const { return child(2); }

  // Whether our children are leaf nodes (true) or internal nodes (false).
  inline bool has_leaves(void) const { return leaf_; }
  // Number of leaves. 0 unless has_leaves() returns true. At most 3.
  inline size_type leaf_count(void) const
    { return has_leaves() ? child_count():0; }
  // Get leaf node i. Unchecked: valid only for 0 <= i < leaf_count().
  inline const lnode_type* leaf(size_type i) const
    { return static_cast<lnode_type*>(u_.children[i]); }
  inline lnode_type* leaf(size_type i)
    { return static_cast<lnode_type*>(u_.children[i]); }

  // Insert at leaf at the given index. MUST be in sorted order,
  // or the tree will be FUBAR.
  inline bool insert_leaf(lnode_type* new_leaf, size_type index)
  {
    if (full() || index >= max_count())
      return false;
    // Scoot over the nodes starting at index and insert the new leaf.
    // We MUST have a null pointer in the end slot which is overwritten
    // (checked by verifying we are not full() above).
    for (size_type lidx = max_count()-1; lidx > index; --lidx)
      u_.leaves[lidx] = u_.leaves[lidx-1];
    u_.leaves[index] = new_leaf;
    ++child_count_;
    // Now fix the keys: key(i) is a copy of key(leaf(i)).
    // If the insert index was 2, the caller will have to fix parent keys.
    for (size_type kidx = 0u; kidx < key_count(); ++kidx)
      keys_[kidx] = leaf(kidx)->key();
    return (leaves_ = true);
  }
  inline bool insert_leaf(const value_type& value, size_type index) {
    if (full() || index >= max_count())
      return false;
    lnode_type *new_leaf = new lnode_type(value);
    return insert_leaf(new_leaf, index);
  }

  inline bool add_leaf(const value_type& value)
    { return insert_leaf(value, leaf_count()); }
  inline bool add_leaf(lnode_type *new_leaf)
    { return insert_leaf(new_leaf, leaf_count()); }
  inline lnode_type* set_leaf(size_type i, lnode_type* new_leaf)
  {
    inode_type* old_child = u_.children[i];
    u_.children[i] = new_child;
    if (new_leaf != nullptr)
      leaves_ = true;
    return old_child;
  }

  inline inode_type* left_leaf(void) { return leaf(0); }
  inline inode_type* middle_leaf(void) { return leaf(1); }
  inline inode_type* right_leaf(void) { return leaf(2); }
  inline const inode_type* left_leaf(void) const { return leaf(0); }
  inline const inode_type* middle_leaf(void) const { return leaf(1); }
  inline const inode_type* right_leaf(void) const { return leaf(2); }

  // If we are full, we have no more room for another child.
  inline bool full(void) const { return child_count() == max_count(); }

  // Boundary keys. All nodes in child(i) are less than key(i).
  // The special case child(2) is just not less than key(1). For key_count():
  //   0 children    -> 0 keys;
  //   1, 2 children -> 1, 2 keys;
  //   3 children    -> 2 keys.
  inline size_type key_count(void) const { return std::max(child_count(), 2u); }
  // Get key i. Unchecked: valid only for 0 <= i < key_count().
  inline key_type& key(size_type i) { return keys_[i]; }
  inline const key_type& key(size_type i) const { return keys_[i]; }
  inline key_type& left_key(void) { return key(0); }
  inline key_type& right_key(void) { return key(1); }
  inline const key_type& left_key(void) const { return key(0); }
  inline const key_type& right_key(void) const { return key(1); }

public:
  // Constructors.
  l23_inode()
    : parent_(nullptr),
      key_count(0), keys_({0}),
      child_count(0), u_.children({0}),
      leaves_(false)
  {}

private:
  ////// Members.
  // Parent back-pointer. Only NULL for the root node.
  inode_type* parent_;

  // Copies of keys used to sort the children. We never have more than 2 keys,
  // since we never have more than three children.
  key_array keys_;

  // Pointers to our child nodes. Always have 0, 2, or 3 children.
  // If our child nodes are leaves, we instead store pointers to leaf nodes.
  size_type child_count_;
  union {
    node_array children;
    leaf_array leaves;
  } u_;

  // If leaves_ is true, the children of this node are leaves.
  // Otherwise the children of this node are internal nodes.
  // This saves space by storing data only in the leaves, rather than
  // in internal nodes.
  bool leaves_;
};

template<class Params>
class l23_rootnode : public l23_inode<Params>
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
  friend class lnode_type;
  friend class tree_type;

public:
  // Constructors.
  l23_rootnode()
    : inode_type(), leftmost(nullptr), rightmost(nullptr)
  {}

private:
  // Members.
  inode_type *leftmost, *rightmost;
  size_type tree_size;
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
  typedef l23_inode<Params> inode_type;
  typedef l23_lnode<Params> lnode_type;
  typedef l23_tree<Params> tree_type;
  friend class inode_type;
  friend class tree_type;

public:
  l23_lnode()
    : parent_(nullptr), value_()
  {}

  template<ValueType>
  l23_lnode(ValueType v)
    : parent_(nullptr), value_(v)
  {}

  inline inode_type* parent(void) { return parent_; }
  inline const inode_type* parent(void) const { return parent_; }

  inline value_type& value(void) { return value_; }
  inline const value_type& value(void) const { return value_; }

  inline key_type& key(void) { return Params::key(value()); }
  inline const key_type& key(void) const { return Params::key(value()); }

  inline data_type& data(void) { return Params::data(value()); }
  inline const data_type& data(void) const { return Params::data(value()); }

private:
  // Members. A leaf node just has a value and a back-pointer to its parent.
  inode_type* parent_;
  value_type value_;
};

template<class Params>
class l23_iterator
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
  typedef l23_rootnode<Params> root_type;

public:
  // Constructors.

  // Iterator pointing to nptr itself.
  l23_iterator(inode_type* nptr)
    : node_(nptr), pos_(-1)

  // Iterator pointing to the p'th child of nptr.
  l23_iterator(inode_type* nptr, size_type p)
    : node_(nptr), pos_(nptr ? p : -1)
  {}

  // Invalid iterator.
  l23_iterator()
    : node_(nullptr), pos_(-1)
  {}

  // Copy constructor.
  l23_iterator(const l23_iterator& other)
    : node_(other.node_), pos_(other.pos_)
  {}

  // Whether this points to a leaf node.
  // Note that the actual leaf node pointed to may be NULL.
  inline bool is_leaf(void) const {
    return !is_parent() && node_ && node_->has_leaves();
  }

  // Only valid when leaf() is true.
  inline lnode_type* leaf(void) { return node_->leaf(pos_); }
  inline value_type& value(void) { return leaf()->value(); }
  inline key_type& key(void) { return leaf()->key(); }

  // Return the node pointed to by this iterator (if any).
  // Only valid when leaf() is false.
  inline inode_type* node(void)
    { return is_parent() ? node_ : node_->child(pos_); }
  // Return the parent of node(). Useful in the case that node() is NULL, but
  // node's parent is not.
  inline inode_type* parent(void)
    { return is_parent() ? node_->parent() : node_; }

  // Whether this points to node (true) or a child of node (false).
  inline bool is_parent(void) const { return pos_ < 0; }

  // Whether this points to a valid node (false for tree.end()).
  inline operator bool(void) const { return node_ != nullptr; }

  // Position of the pointed-to node in the parent node.
  inline int pos(void) const { return pos_; }

  // pre-fix
  inline iterator& operator++(void)
  {
    // If we can visit the next node in the same parent, go ahead.
    if (pos() < node_->child_count() - 1)
      ++pos_;
    else {
      node_ = nullptr;
      pos_ = -1;
    }
    // TODO DFS by increment?
    return *this;
  }

  // post-fix
  inline iterator operator++(int)
  {
    iterator it(*this);
    ++(*this);
    return it;
  }

  inline bool operator==(const iterator& other) const
    { return (pos_ == other.pos_ && node_ == other.node_); }
  inline bool operator!=(const iterator& other) const
    { return !(*this == other); }

private:
  // Members.
  inode_type* node_;
  // If pos is negative, points to node (which is the root node).
  int pos_;
};

//! 2-3 tree with data in the leaves.
template<class Params>
class l23_tree
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
  typedef l23_rootnode<Params> root_type;
  typedef l23_iterator<Params> iterator;

public:
  // Constructors.
  template<class KeyCompare>
  l23_tree(KeyCompare kcmp=key_compare())
    : root_(), compare_(kcmp)
  {}

  root_type* root(void) { return &root_; }

  // Iterators.
  iterator begin(void) { return iterator(root, 0); }
  iterator end(void) { return iterator(nullptr, 0); }

  // Lookup.
public:
  template<class KeyType=const key_type&>
    iterator lower_bound(KeyType k);
  inline iterator lower_bound(const value_type& v) {
    return lower_bound(Params::key(v));
  }

  template<class KeyType=const key_type&>
    iterator upper_bound(KeyType k);
  inline iterator upper_bound(const value_type& v) {
    return upper_bound(Params::key(v));
  }

  template<class KeyType=const key_type&>
    iterator find_unique(KeyType k);
  inline iterator find_unique(const value_type& v) {
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
  template<class KeyType1, class KeyType2, class KeyCompare>
  inline compare(KeyType1 k1, KeyType2 k2, KeyCompare kcmp) const {
    return kcmp(k1, k2);
  }
  template<class KeyType1, class KeyType2>
  inline compare(KeyType1 k1, KeyType2 k2) const {
    return compare_(k1, k2);
  }

  // Return where a new key should go in sorted order as a child of the given
  // node by comparing against its (at most) two key values.
  // This version uses a "less than or equal to" comparison.
  template<class KeyType>
  struct cle
  {
    inline size_type operator()(KeyType key, inode_type* node) const
    {
      // Return the only node if we only have one node.
      // Special behavior: return -1 if we have zero nodes.
      if (node->count() == 0 || !compare(node->left_key(), k))
        return 0u; // no left or k <= left: return left child
      if (node->key_count() > 1 && compare(node->right_key(), k))
        return 2u; // k > right: return right child
      // left < k <= right or no right: return middle child
      return 1u;
    }
  };

  // Return where a new key should go in sorted order as a child of the given
  // node by comparing against its (at most) two key values.
  // This version uses a "strictly less than" comparison.
  template<class KeyType>
  struct clt
  {
    inline size_type operator()(KeyType key, inode_type* node) const
    {
      // Return the only node if we only have one node.
      // Special behavior: return -1 if we have zero nodes.
      if (node->count() == 0 || compare(k, node->left_key()))
        return 0u; // no left or k < left: return left child
      if (node->key_count() > 1 && compare(node->right_key(), k))
        return 2u; // k > right: return right child
      // left < k < right or no right: return middle child
      return 1u;
    }
  };

  // IndexFunctor must be one of 'cle' or 'clt'.
  template<class Params> template<class IndexFunctor, class KeyType>
    iterator internal_lower_bound(IndexFunctor f, KeyType k) const;

private:
  // Members.
  root_type root_;
  key_compare compare_;
};

// Tree implementations.
template<class Params> template<class IndexFunctor, class KeyType>
typename l23_tree<Params>::iterator
l23_tree<Params>::internal_lower_bound(IndexFunctor get_index, KeyType k) const
{
  if (!root()->count())
    return end();
  const inode_type* node = root();
  while (node && !node->has_leaves())
    node = node->child(get_index(k, node));
  // We must have landed on a non-null node pointing to leaves.
  if (!node || !node->has_leaves())
    return end();
  // Return an iterator pointing to the appropriate lower-bound leaf node.
  // Note that the child node may not actually exist.
  return iterator(node, get_index(k, node));
}

template<class Params> template<class KeyType>
typename l23_tree<Params>::iterator
l23_tree<Params>::lower_bound(KeyType k) const
{
  iterator leafp = internal_lower_bound(cle(), k);
  if (leafp.leaf() == nullptr)
    return end();
  return leafp;
}

template<class Params> template<class KeyType>
typename l23_tree<Params>::iterator
l23_tree<Params>::upper_bound(KeyType k) const
{
  iterator leafp = internal_lower_bound(cle(), k);
  // Increment the iterator until we find a key value greater than k.
  while (leafp != end() && !compare(k, leafp.key()))
    ++leafp;
  return leafp; // may be end()
}

template<class Params> template<class KeyType>
typename l23_tree<Params>::iterator
l23_tree<Params>::find_unique(KeyType k)
{
  iterator leafp = internal_lower_bound(cle(), k);
  if (leafp == end() || leafp.leaf() == nullptr)
    return end();
  // Only return the given node if its key is equal to the query key.
  // We already know compare(k, leafp.key()) is false by definition of
  // lower_bound (k NOT < leaf). If the converse is false we have equality
  // by STL convention.
  if (!compare(leafp.key(), k))
    return iterator(node, leaf_idx);
  return end();
}

template<class Params>
typename l23_tree<Params>::iterator
l23_tree<Params>::insert_unique(const value_type& value)
{
  // Special case: root has fewer than two child leaves.
  if (root()->count() == 0)
  {
    root()->add_leaf(value);
    return iterator(root(), 0);
  }
  if (root()->count() == 1)
  {
    // Insert the node in the right spot.
    const key_type& key = Params::key(value);
    size_t pos = compare(key, root()->leaf(0)->key()) ? 0 : 1;
    root()->insert_leaf(value, pos);
    return iterator(root(), pos);
  }

  // Now without loss of generality, every node has either 2 or 3 children.
  // Find the proper parent node to contain the new leaf.
  iterator leafp = internal_lower_bound(clt(), Params::key(value));
  if (leafp == end())
    return end();

  // If the node is not full, we can just add the leaf directly.
  inode_type* parent = leafp.parent();
  if (parent->insert_leaf(value, leafp.pos()))
  {
    // But, if we added a new max leaf, we must fix our ancestors' keys
    // because the new upper bound for this subtree may be higher than the
    // upper bound key from our parent.
    if (leafp.pos() == inode_type::max_count() - 1)
    {
      const key_type& new_key = Params::key(value);
      inode_type* pnode = parent->parent();
      while (pnode)
      {
        if (pnode->key_count() > 1 && compare(pnode->key(1), new_key))
          pnode->set_key(1, new_key);
        pnode = pnode->parent();
      }
    }
    return iterator(parent, leafp.pos());
  }

  // Here the parent of the new leaf is full, so we must rebalance.
  // TODO
}



} // end namespace tree
