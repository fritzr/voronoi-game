//------------------------------------------------------------------------------
//  Copyright 2007-2011 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _POLYGON_DATA_H_
#define _POLYGON_DATA_H_

#ifdef WIN32
#pragma warning(disable : 4786)
#endif

#include "opencv_compat.h"
#include "boost_cv_compat.h"

#include <boost/iterator/indirect_iterator.hpp>

#include <functional> // reference_wrapper
#include <iostream>
#include <list>
#include <vector>
#include <cassert>
#include <float.h>
#include <iterator>
using namespace cv;

#if CV_MAJOR_VERSION < 3
using namespace cv::traits;
#endif

namespace bp = boost::polygon;

typedef unsigned int uint;

template<typename T>
struct PI { };

template<>
struct PI<float> {
  float value = 6.2831853071794f;
};

template<>
struct PI<double> {
  double value = 6.2831853071794;
};

//
//a triangle, used in triangulating a polygon
//
struct triangle
{
    uint v[3]; // id to the vertices
};

// Forward declarations
template<typename Pt_>
class ply_vertex;

template <typename Pt, bool is_const>
class vertex_iterator;

template<typename Pt_>
class c_ply;

template<typename Pt_>
class c_plylist;

template<typename Pt_>
class c_polygon;

template<typename Pt_>
class PlyMat;

template <typename Pt_> std::istream &
  operator>>(std::istream&, c_ply<Pt_>&);

template <typename Pt_> std::ostream &
  operator<<(std::ostream &, const c_ply<Pt_> &);

template <typename Pt_> std::istream &
  operator>>(std::istream&, c_plylist<Pt_>&);

template <typename Pt_> std::ostream &
  operator<<(std::ostream &, const c_plylist<Pt_> &);

template <typename Pt_> std::ostream & \
  operator<<(std::ostream &, const c_polygon<Pt_> &);

template <typename Pt_> std::ostream & \
  operator<<(std::ostream &, const std::vector<c_polygon<Pt_> >&);

//
// extra information about the polygon
//

template<typename Pt_>
struct ply_extra_info
{
  typedef typename ply_vertex<Pt_>::coordinate_type coordinate_type;

  ply_extra_info()
  {
      ftt=0.0;
      has_box=false;
      box[0]=box[1]=box[2]=box[3]=0;
      pointIndex=UINT_MAX;
  }

  bool has_box;
  coordinate_type ftt; // fixed-travel-time for this isoline
  coordinate_type box[4];
  uint pointIndex; // index of the corresponding center point
};


//
// Vertex of polygon
//
template<typename Pt_>
class ply_vertex
{
public:
    typedef Pt_ point_type;
    typedef typename bp::point_traits<point_type>::coordinate_type
      coordinate_type;

    friend class vertex_iterator<Pt_, false>;
    friend class vertex_iterator<Pt_, true>;

    ///////////////////////////////////////////////////////////////////////////
    ply_vertex(){ init(); }
    ply_vertex( const point_type& p ){ pos=p; init(); }
    virtual ~ply_vertex();
    void setNext(ply_vertex * n){next=n; if(n!=NULL) n->pre=this; }
    void setPre(ply_vertex * n){pre=n; if(n!=NULL) n->next=this; }
    void computeExtraInfo();

    //negate the vertex
    void negate();

    //reverse the order
    void reverse();

    //copy
    void copy(ply_vertex * other);

    operator point_type(void) const { return getPos(); }

    ///////////////////////////////////////////////////////////////////////////
    void setPos(const point_type& p) { pos=p; }
    void setX(coordinate_type x) { pos.x=x; }
    void setY(coordinate_type y) { pos.y=y; }
    virtual const point_type& getPos() const { return pos; }

    void translate(const point_type& v){ pos=pos+v; }

    virtual ply_vertex * getNext() const { return next; }
    virtual ply_vertex * getPre() const { return pre; }

    const Vec<coordinate_type, 2>& getNormal() const { return normal; }
    bool isReflex() const { return reflex; }

    //get extra information
    uint getVID() const { return vid; }
    void setVID(uint id) {vid=id;}

private:

    void init(){
        next=pre=NULL; 
        reflex=false;
        vid=UINT_MAX;
    }

    //basic info
    point_type pos;       //position
    ply_vertex * next; //next vertex in the polygon
    ply_vertex * pre;  //previous vertex in the polygon
    //normal, the segment normal from this v to the next.
    Vec<coordinate_type, 2> normal;
    bool reflex;
    uint vid;
};

template <typename Pt, bool is_const>
  class vertex_iterator
    : public boost::iterator_facade<
        /* Derived */ vertex_iterator<Pt, is_const>
        /* Value */, typename std::conditional<is_const, const Pt&, Pt&>::type
        /* Category */, boost::random_access_traversal_tag
      >
{
public:
  typedef boost::iterator_facade<vertex_iterator,
          Pt, boost::random_access_traversal_tag,
          typename std::conditional<is_const, const Pt&, Pt&>::type>
    super;

  typedef typename super::value_type value_type;
  typedef typename super::reference reference;
  typedef typename super::difference_type difference_type;

  typedef c_ply<Pt> ring_type;
  typedef typename std::conditional<
    is_const, const ply_vertex<Pt>, ply_vertex<Pt>
    >::type vertex_type;

private:
  vertex_type *head_;
  vertex_type *cur_;
  size_t sz_;

  struct enabler{};

public:
  vertex_iterator() : head_(nullptr), cur_(nullptr), sz_(0u) {}
  vertex_iterator(const ring_type &ply)
    : head_(ply.head), cur_(NULL), sz_(ply.getSize()) {}
  vertex_iterator(vertex_type *head, size_t size)
    : head_(head), cur_(NULL), sz_(size) {}
  vertex_iterator(vertex_type *head, vertex_type *next, size_t size)
    : head_(head), cur_(next), sz_(size) {}

  // copy from const vertex iterator
  template<bool other_const>
  vertex_iterator(vertex_iterator<Pt, other_const> const& other,
      typename std::enable_if<other_const, enabler>::type = enabler())
    : head_(other.head_), cur_(other.head_), sz_(other.sz_) {}

private:
  friend class boost::iterator_core_access;
  friend class vertex_iterator<Pt, !is_const>;

  inline reference dereference(void) const {
    return cur_ ? cur_->pos : head_->pos;
  }
  template <bool other_const>
  inline bool equal(vertex_iterator<Pt, other_const> const& o) const {
    return head_ == o.head_ && cur_ == o.cur_;
  }
  inline vertex_iterator &increment() {
    cur_ = (cur_ ? cur_->getNext() : head_->getNext());
    return *this;
  }
  inline vertex_iterator &decrement() {
    cur_ = (cur_ ? cur_->getPre() : head_->getPre());
    return *this;
  }
  inline vertex_iterator &advance(difference_type n) {
    std::advance(*this, n); return *this;
  }

  template <bool other_const>
  inline difference_type distance_to(vertex_iterator<Pt, other_const> const& o)
    const
  {
    difference_type dist;
    const vertex_iterator  end = vertex_iterator(head_, head_, sz_);
    const vertex_iterator rend = vertex_iterator(head_, head_->getPre(), sz_);

    dist = 0;
    for (vertex_iterator copy = *this; copy != end; ++copy, ++dist)
    {
      if (copy == o)
        return dist;
    }

    dist = 0;
    for (vertex_iterator copy = *this; copy != rend; --copy, --dist)
      if (copy == o)
        return dist;

    return dist;
  }

};

//
// Polygon chain
//
template<typename Pt_>
class c_ply
{
public:
    typedef Pt_ point_type;
    typedef ply_vertex<Pt_> vertex_type;
    typedef typename vertex_type::coordinate_type coordinate_type;

    friend class vertex_iterator<Pt_, false>;
    friend class vertex_iterator<Pt_, true>;
    friend class c_polygon<Pt_>;
    friend class PlyMat<Pt_>;

    typedef vertex_iterator<point_type, true> const_iterator;
    typedef vertex_iterator<point_type, false> iterator;

    ///////
    // Iterate through points directly.
    iterator begin(void) { return iterator(*this); }
    iterator end(void) { return iterator(head, head, getSize()); }

    const_iterator begin(void) const { return const_iterator(*this); }
    const_iterator cbegin(void) const { return const_iterator(*this); }

    const_iterator end(void) const
      { return const_iterator(head, head, getSize()); }
    const_iterator cend(void) const
      { return const_iterator(head, head, getSize()); }

    enum POLYTYPE { UNKNOWN, PIN, POUT };

    ///////////////////////////////////////////////////////////////////////////
    c_ply(POLYTYPE t)
      : head(NULL), tail(NULL), all(), center(), radius(-1),
      area(std::numeric_limits<coordinate_type>::lowest()),
      type(t), extra_info(), triangulation()
    {}

    ///////////////////////////////////////////////////////////////////////////
    void copy(const c_ply<Pt_>& ply); //copy from the other ply
    void destroy();

    ///////////////////////////////////////////////////////////////////////////
    // create c_ply
    void beginPoly();
    void addVertex(coordinate_type x, coordinate_type y, bool remove_duplicate);
    inline void addVertex(coordinate_type x, coordinate_type y) {
      return addVertex(x, y, false);
    }
    void addVertex( vertex_type * v );

    void endPoly(bool remove_duplicate);
    inline void endPoly(void) { endPoly(false); }

    ///////////////////////////////////////////////////////////////////////////
    void negate();
    void reverse(); //reverse vertex order

    ///////////////////////////////////////////////////////////////////////////
    void translate(const point_type& p);

    ///////////////////////////////////////////////////////////////////////////
    //
    point_type findEnclosedPt(); //find a point that is enclosed by this polychain

    //triangulate the polygon
    void triangulate(std::vector<triangle>& tris);

    ///////////////////////////////////////////////////////////////////////////
    // Access functions
    vertex_type * getHead() const { return head; }

    //hole or external boundary
    POLYTYPE getType() const { return type; }

    void set(POLYTYPE t,vertex_type * h){ 
        type=t; head=h; 
        if(h!=NULL){ tail=h->getPre(); }
        else{ tail=NULL; }
    }

    int getSize() {
        if(all.empty()) indexing();
        return all.size();
    }

    int getSize() const {
        return all.size();
    }

    vertex_type * operator[](unsigned int id){
        if(all.empty()) indexing();
        return all[id];
    }

    //get com //center of mass
    const point_type& getCenter();

	//compute the Radius of the poly chain
	coordinate_type getRadius();
	
	//area
        coordinate_type getArea();
        inline coordinate_type getArea() const { return area; }

    // Convert to a cv::Matrix
    operator PlyMat<Pt_>() { return mat(); }
    inline PlyMat<Pt_> mat(void) const { return PlyMat<Pt_>(*this); }

    //
	// additional functions
	//

	//check if a point is enclosed
	//the behavior is unknown if pt is on the boundary of the polygon
	bool enclosed(const point_type& pt);

	//check if convex
	bool is_convex() const;

	//delete a vertex
	void delete_vertex(vertex_type * p);

	//get extra info
	ply_extra_info<Pt_>& extra(){ return extra_info; }
	const ply_extra_info<Pt_>& extra() const { return extra_info; }

    ///////////////////////////////////////////////////////////////////////////
    // Operator
    //check if give poly line is the same as this
    bool operator==( const c_ply<Pt_>& other ){ return other.head==head; }

    //
    //this is a more general check if other and this ply are identical
    //
    bool identical(const c_ply<Pt_>& other);

    friend std::istream &
      operator>> <Pt_>(std::istream&, c_ply<Pt_>&);

    friend std::ostream &
      operator<< <Pt_>(std::ostream &os, const c_ply<Pt_> &p);

protected:

    ///////////////////////////////////////////////////////////////////////////
    void doInit(); /*return # of vertice in this poly*/

    //indexing all elements and store in vector<ply_vertex*>
    void indexing();

private:

    vertex_type * head; //the head of vertex list
    vertex_type * tail; //end of the vertex list

    std::vector<vertex_type*> all; //all vertices
	
	//additional info
	point_type center;
	coordinate_type radius;
	coordinate_type area;

    //In, out or unknown.
    POLYTYPE type;

    //extrac info
    ply_extra_info<Pt_> extra_info;

    //triangulation
    std::vector<triangle> triangulation; //catched triangulation, calculated by triangulate
};


/* This is just a viewport into a c_ply to provide a CV adapter.  */
template<typename Pt_>
class PlyMat : public cv::Mat_<typename c_ply<Pt_>::point_type >
{
private:
  typedef cv::Mat_<typename c_ply<Pt_>::point_type > super;

public:
  typedef c_ply<Pt_> ring_type;
  typedef typename ring_type::vertex_type vertex_type;
  typedef typename vertex_type::point_type point_type;
  typedef typename vertex_type::coordinate_type coordinate_type;

private:
  const ring_type &ply_;

public:
  PlyMat(const ring_type &ply)
    : ply_(ply) {}

  // Copy assignment
  PlyMat &operator=(const PlyMat &m);

  const point_type &at(int pos) {
    return ply_.all[pos]->getPos();
  }
};

//a c_plylist is a list of c_ply
template <typename Pt_>
class c_plylist : public std::list<c_ply<Pt_> >
{
public:
    typedef std::list<c_ply<Pt_> > super;
    typedef c_ply<Pt_> ring_type;
    typedef typename ring_type::vertex_type vertex_type;
    typedef typename ring_type::point_type point_type;
    typedef typename vertex_type::coordinate_type coordinate_type;

    using super::iterator;
    using super::const_iterator;
    using super::empty;
    using super::size;
    using super::begin;
    using super::cbegin;
    using super::end;
    using super::cend;
    using super::front;
    using super::back;
    using super::clear;
    using super::push_back;

private:
    friend std::ostream&
      operator<< <Pt_> (std::ostream&, const c_plylist<Pt_>&);

    friend std::istream&
      operator>> <Pt_>(std::istream&, c_plylist<Pt_>& );

public:

    c_plylist()
    {
        box[0]=box[1]=box[2]=box[3]=0;
        is_buildboxandcenter_called=false;
    }

    void translate(const point_type& v);

    //access
    void buildBoxAndCenter();
    coordinate_type * getBBox() { assert(is_buildboxandcenter_called); return box; }
    const point_type& getCenter() { assert(is_buildboxandcenter_called); return center; }

protected:

    point_type center;
    coordinate_type box[4];

private:

    bool is_buildboxandcenter_called;
};


// Forward decl
namespace boost { namespace polygon {
  template<typename Pt> struct polygon_traits_general<c_polygon<Pt> >;
} }

template<typename RefIter>
class ring_view_iterator;

namespace std {
  template<typename RefIter>
  class iterator_traits<ring_view_iterator<RefIter> >
  {
    typedef ptrdiff_t difference_type;
    typedef typename RefIter::value_type::type value_type;
    typedef value_type& reference;
    typedef value_type* pointer;
  };
}

template<typename RefIter>
class ring_view_iterator
  : public boost::iterators::indirect_iterator<
      ring_view_iterator<RefIter>
    , typename RefIter::value_type::type
    , std::bidirectional_iterator_tag
    >
{
public:
  typedef boost::iterators::indirect_iterator<
      ring_view_iterator
    , typename RefIter::value_type::type
    , std::bidirectional_iterator_tag
    > super;

  typedef typename std::iterator_traits<ring_view_iterator> traits;
  typedef typename traits::value_type value_type;
  typedef typename traits::reference reference;
  typedef typename traits::pointer pointer;

  ring_view_iterator() : inner() {}
  ring_view_iterator(RefIter i) : inner(i) {}

  // override
  reference operator*() { return inner->get(); }
  pointer operator->() { return &inner->get(); }

private:
  RefIter inner;
};

/* XXX Copy and return a list of only our inner rings.  */
template<typename Pt>
class ring_view
  : public std::list<
      boost::reference_wrapper<const typename c_polygon<Pt>::ring_type>
    >
{
public:
  typedef c_polygon<Pt> polygon_type;
  typedef typename polygon_type::ring_type ring_type;

  typedef std::list<boost::reference_wrapper<const ring_type> > super;
  typedef typename super::const_iterator list_iterator;

  typedef typename super::value_type wrapper_type;
  typedef typename super::value_type::type value_type;

  typedef ring_view_iterator<typename super::const_iterator>
    const_iterator;

  template<typename Iter>
  ring_view(Iter begin, Iter end)
    : super(begin, end) {}

  const_iterator begin(void) const { return const_iterator(super::begin()); }
  const_iterator end(void) const { return const_iterator(); }
};

//
// a c_polygon is a restricted kind of c_plylist
// this defines a simple polygon so that
// the first element much be a POUT c_ply and
// the rest ply lines are holes
//
template<typename Pt>
class c_polygon : public c_plylist<Pt>
{
public:
    friend struct boost::polygon::polygon_traits_general<c_polygon>;

    typedef c_plylist<Pt> super;
    typedef Pt point_type;

    typedef c_ply<point_type> ring_type;
    typedef c_plylist<point_type> list_type;
    typedef typename ring_type::vertex_type vertex_type;
    typedef typename vertex_type::coordinate_type coordinate_type;

    typedef std::vector<triangle> tri_container;
    typedef tri_container::const_iterator tri_iterator;

    // For viewing the inner rings as required by boost
    typedef ring_view<Pt> ring_view_type;

    using super::iterator;
    using super::const_iterator;
    using super::empty;
    using super::size;
    using super::begin;
    using super::cbegin;
    using super::end;
    using super::cend;
    using super::front;
    using super::back;
    using super::clear;
    using super::push_back;

    c_polygon() : super() { area=0; }
    c_polygon(const ring_type &ply) : super() { super::push_back(ply); }

    bool valid() const; //check if this is a valid polygon

    //copy from the given polygon
    void copy(const c_polygon& other);

    //triangulate the polygon
    void triangulate(std::vector<triangle>& tris);

    // cache the triangulation for later
    void triangulate(void) {
      triangulate(triangulation);
    }

    tri_iterator triangles_begin(void) const { return triangulation.begin(); }
    tri_iterator triangles_end(void) const { return triangulation.begin(); }
    size_t triangles_size(void) const { return triangulation.size(); }

    void reverse(); //reverse the vertex order (not the list order)

    //access the vertices of the polygon as an array
    int getSize()
    {
        if(all.empty()) indexing();
        return all.size();
    }

    vertex_type * operator[](unsigned int id){
        if(all.empty()) indexing();
        return all[id];
    }

    vertex_type * operator[](unsigned int id) const {
      return all.at(id); // checked
    }

    coordinate_type getArea();
    inline coordinate_type getArea() const { return area; }

    // Remove the first c_ply.
    void pop_front();

    //destroy
    void destroy();

    //check if a point is enclosed
    //the behavior is unknown if pt is on the boundary of the polygon
    bool enclosed(const point_type& pt) const;

    // We need to triangulate first to get the right answer.
    bool enclosed(const point_type& pt) {
      if (triangulation.empty())
          triangulate(triangulation);
      return const_cast<const c_polygon*>(this)->enclosed(pt);
    }

    //find a point inside the polygon
    point_type findEnclosedPt();

    //get number of vertices
    uint getSize() const;

    bool is_convex() const;

    //check if polygons are identical (or not)
    bool identical(c_polygon& p);

    friend std::ostream &
      operator<< <Pt>(std::ostream &os, const c_polygon<Pt> &p);

    friend std::ostream &
      operator<< <Pt>(std::ostream &os, const std::vector<c_polygon<Pt> > &v);

    ring_view_type inner(void) const {
      return ring_view_type(std::next(this->begin()), this->end());
    }

private:

    //indexing the vertices and store them in vector<ply_vertex*> all
    void indexing();

    std::vector<vertex_type*> all; //all vertices

    //triangulation
    tri_container triangulation; //cached triangulation, calculated by triangulate

    coordinate_type area;
};


//
// Given radius and resolution (# of vertices), create a circle in polygon format
//
template <typename Pt_>
inline void create_circle(c_polygon<Pt_>& p,
    typename c_polygon<Pt_>::coordinate_type radius, uint res)
{
    typedef typename c_polygon<Pt_>::coordinate_type T;
    typedef typename c_polygon<Pt_>::ring_type ring_type;

    T delta = PI<T>::value / res;
    ring_type ply(ring_type::POUT);
    ply.beginPoly();
    for(uint i=0;i<res;i++){
        T r=delta*i;
        T x=cos(r)*radius;
        T y=sin(r)*radius;
        ply.addVertex(x,y);
    }//end for i
    ply.endPoly();
    p.push_back(ply);
}


// PLY_INSTANTIATE(point_type, [extern])
//
// Declare (with extern) or instantiate (with empty second arg) polygon classes
// using the given point type.
#define PLY_INSTANTIATE(Pt, E) \
  E template std::istream& \
    operator>> <Pt>(std::istream&, c_ply<Pt>& ); \
  E template std::ostream & \
    operator<< <Pt>(std::ostream &, const c_ply<Pt> &);\
  E template std::istream& \
    operator>> <Pt>(std::istream&, c_plylist<Pt>& ); \
  E template std::ostream & \
    operator<< <Pt>(std::ostream &, const c_plylist<Pt> &);\
  E template std::ostream & \
    operator<< <Pt>(std::ostream &, const c_polygon<Pt> &);\
  E template std::ostream & \
    operator<< <Pt>(std::ostream &, const std::vector<c_polygon<Pt> > &);\
  E template class ply_vertex<Pt>; \
  E template class c_ply<Pt>; \
  E template class c_plylist<Pt>; \
  E template class c_polygon<Pt>; \

PLY_INSTANTIATE(cv::Point2d, extern);
PLY_INSTANTIATE(cv::Point2f, extern);

#endif //_POLYGON_H_
