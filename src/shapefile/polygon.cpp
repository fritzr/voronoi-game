//------------------------------------------------------------------------------
//  Copyright 2007-2011 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "polygon.h"
#include "intersection.h"
#include <vector>
#include <cmath>

#include <boost/polygon/polygon.hpp>

#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/arithmetic/arithmetic.hpp>
#include <boost/geometry/arithmetic/dot_product.hpp>
#include <boost/geometry/algorithms/equals.hpp>
#include <boost/geometry/util/math.hpp>

#include "adapt_boost_poly.h"

using namespace std;
namespace bp = boost::polygon;
namespace bg = boost::geometry;
namespace bgmth = boost::geometry::math;

#ifdef WIN32
extern "C"{
#include "triangulate.h"
}
#else
#include "triangulate.h"
#endif

template<typename Pt_>
ply_vertex<Pt_>::~ply_vertex()
{
    //doing nothing for now
}

template<typename Pt_, typename T>
void
normalize(Pt_ &p, T &magnitude)
{
  typedef typename bg::coordinate_type<Pt_>::type coord_t;

  // normalize
  magnitude = bgmth::sqrt(bg::dot_product(p, p));
  if (!bgmth::equals(magnitude, coord_t(0)))
    bg::divide_value(p, magnitude);
}

template<typename Pt_>
void normalize(Pt_ &p)
{
  typedef typename bg::coordinate_type<Pt_>::type coord_t;

  coord_t magnitude_unused;
  normalize(p, magnitude_unused);
}

// - compute normal
// - check if the vertex is reflex or not
template<typename Pt_>
void
ply_vertex<Pt_>::computeExtraInfo()
{
    //compute normal direction
    point_type v= next->pos;
    bg::subtract_point(v, pos);
    if( bg::get<0>(v) == 0 ){
        if (bg::get<1>(v)>0) { bg::set<0>(normal,1); bg::set<1>(normal,0); }
        else { bg::set<0>(normal,-1); bg::set<1>(normal,0); }
    }
    else if( bg::get<0>(v)>0 ){
      bg::set<1>(normal, -1);
      bg::set<0>(normal, bg::get<1>(v) / bg::get<0>(v));
    }
    else{//get<0>(v)<0
      bg::set<1>(normal, 1);
      bg::set<0>(normal, -(bg::get<1>(v) / bg::get<0>(v)));
    }

    // normalize
    normalize(normal);

    // compute left or right turn
    reflex = !leftTurn(pre->pos, pos, next->pos);
}

template<typename Pt_>
void
ply_vertex<Pt_>::negate()
{
    bg::multiply_value(normal, -1);
    bg::multiply_value(pos, -1);
}

template<typename Pt_>
void
ply_vertex<Pt_>::reverse()
{
    swap(next,pre);
    //normal=-normal;
    computeExtraInfo();
}

template<typename Pt_>
void
ply_vertex<Pt_>::copy(ply_vertex * other)
{
    pos=other->pos;
    normal=other->normal;
    reflex=other->reflex;
    vid=other->vid;
}

///////////////////////////////////////////////////////////////////////////////

//copy from the given ply
template<typename Pt_>
void
c_ply<Pt_>::copy(const c_ply<Pt_>& other)
{
    destroy();//detroy myself first

    vertex_type* ptr=other.head;
    beginPoly();
    do{
        vertex_type * v=new vertex_type();
        assert(v); //check for memory
        v->copy(ptr);
        addVertex(v);
        ptr=ptr->getNext();
    }while( ptr!=other.head );

    //endPoly();
    //finish up
    tail->setNext(head);

    //copy extra info
    area=other.area;
    bg::assign_point(center, other.center);
    radius=other.radius;
    type=other.type;
    triangulation=other.triangulation;
    extra_info=other.extra_info;
}

// clean up the space allocated
template<typename Pt_>
void
c_ply<Pt_>::destroy()
{
    if( head==NULL ) return;
    vertex_type* ptr=head;
    do{
        vertex_type * n=ptr->getNext();
        delete ptr;
        ptr=n;
    }while( ptr!=head );
    head=tail=NULL;

    all.clear();
    triangulation.clear();
}

// Create a empty polygon
template<typename Pt_>
void
c_ply<Pt_>::beginPoly()
{
    head=tail=NULL;
    all.clear();
    triangulation.clear();
}

// Add a vertex to the polygonal chian
template<typename Pt_>
void
c_ply<Pt_>::addVertex(
    typename c_ply<Pt_>::coordinate_type x,
    typename c_ply<Pt_>::coordinate_type y,
    bool remove_duplicate )
{
    point_type pt(x,y);

    if(tail!=NULL){
        if(bg::equals(tail->getPos(), pt) && remove_duplicate) return; //don't add
    }

    vertex_type * v=new vertex_type(pt);
    if( tail!=NULL ){
        tail->setNext(v);
    }
    tail=v;
    if( head==NULL ) head=tail;
    v->setVID(all.size()); //id of the vertex in this ply
	all.push_back(v);

}

// Add a vertex to the polygonal chian
template<typename Pt_>
void
c_ply<Pt_>::addVertex( vertex_type * v )
{
    if( tail!=NULL ){
        tail->setNext(v);
    }
    tail=v;
    if( head==NULL ) head=tail;
    v->setVID(all.size()); //id of the vertex in this ply
    all.push_back(v);
}

// finish building the polygon
template<typename Pt_>
void
c_ply<Pt_>::endPoly(bool remove_duplicate)
{
    if(head!=NULL && tail!=NULL){
        if(remove_duplicate){
            if(bg::equals(head->getPos(), tail->getPos())){ //remove tail..
                delete tail;
                all.pop_back();
                tail=all.back();
            }
        }//
    }

    tail->setNext(head);
    doInit();
}

// initialize property of the this polychain
// Compute normals and find reflective vertices
template<typename Pt_>
void
c_ply<Pt_>::doInit()
{
    //compute area
    getArea();
    if(this->area<0 && type==POUT){
       //cerr<<"! Warning: polygon type is POUT but has negative area. Reverse the vertex ordering."<<endl;
       reverse();
    }
    else if(this->area>0 && type==PIN){
       //cerr<<"! Warning: polygon type is PIN but has positive area. Reverse the vertex ordering."<<endl;
       reverse();
    }

    //compute normals
    vertex_type* ptr=head;
    do{
        ptr->computeExtraInfo();
        ptr=ptr->getNext();
    }while( ptr!=head );
}

template<typename Pt_>
const typename c_ply<Pt_>::point_type&
c_ply<Pt_>::getCenter()
{
    if(radius<0){
        center = point_type(0,0);
        vertex_type * ptr=head;
        const point_type& first=ptr->getPos();
        uint size=0;
        do{
            size++;
            point_type v = ptr->getPos();
            bg::subtract_point(v, first);
            bg::add_point(center, v);
            ptr=ptr->getNext();
        }while(ptr!=head); //end while
        bg::divide_value(center, size);
        bg::add_point(center, first);

        radius=0;
    }

    return center;
}


///////////////////////////////////////////////////////////////////////////
template<typename Pt_>
void
c_ply<Pt_>::negate()
{
    vertex_type * ptr=head;
    do{
        ptr->negate();
        ptr=ptr->getNext();
    }while(ptr!=head); //end while
}

//reverse the order of the vertices
template<typename Pt_>
void
c_ply<Pt_>::reverse()
{
    vertex_type * ptr=head;
    do{
        ptr->reverse();
        ptr=ptr->getNext();
    }
    while(ptr!=head); //end while

    this->area=-this->area;
    all.clear();
    triangulation.clear();
}

///////////////////////////////////////////////////////////////////////////
template<typename Pt_>
void
c_ply<Pt_>::translate(const point_type& p)
{
    vertex_type * ptr=head;
    do{
        ptr->translate(p);
        ptr=ptr->getNext();
    }while(ptr!=head); //end while

    //translate box
    extra_info.box[0] += bg::get<0>(p);
    extra_info.box[1] += bg::get<0>(p);
    extra_info.box[2] += bg::get<1>(p);
    extra_info.box[3] += bg::get<1>(p);
}

///////////////////////////////////////////////////////////////////////////
//compute the Radius of the poly chain
template<typename Pt_>
typename c_ply<Pt_>::coordinate_type
c_ply<Pt_>::getRadius()
{
    if(radius<0) getCenter();

    if(radius==0){
        vertex_type * ptr=head;
        do{
            point_type v = center;
            bg::subtract_point(v, ptr->getPos());
            coordinate_type d;
            normalize(v, d);
            if(d>radius) radius=d;
            ptr=ptr->getNext();
        }while(ptr!=head); //end while
        radius=sqrt(radius);
    }

    return radius;
}


template<typename Pt_>
typename c_ply<Pt_>::coordinate_type
c_ply<Pt_>::getArea()
{
    if (area == numeric_limits<coordinate_type>::lowest())
      area = bp::area(*this);

    return area;
}


template<typename Pt_>
bool
c_ply<Pt_>::enclosed(const typename c_ply<Pt_>::point_type& p)
{
    if(triangulation.empty())
        triangulate(triangulation);

    //find the largest triangle
    for(uint i=0;i<triangulation.size();i++){
        triangle & tri=triangulation[i];
        const point_type& p1=all[tri.v[0]]->getPos();
        const point_type& p2=all[tri.v[1]]->getPos();
        const point_type& p3=all[tri.v[2]]->getPos();
        coordinate_type area1=Area(p1,p2,p);
        coordinate_type area2=Area(p2,p3,p);
        coordinate_type area3=Area(p3,p1,p);
        if(area1>=0 && area2>=0 && area3>=0) return true; //in
        if(area1<=0 && area2<=0 && area3<=0) return true; //in
    }

    return false; //out
}

template<typename Pt_>
typename c_ply<Pt_>::point_type
c_ply<Pt_>::findEnclosedPt()
{
    if(triangulation.empty())
        triangulate(triangulation);

    if(triangulation.empty()){ //too few points...
        const point_type& p1=head->getPos();
        const point_type& p2=head->getNext()->getPos();
        point_type pt;
        bg::assign_point(pt, p1);
        bg::add_point(pt, p2);
        bg::divide_value(pt, 2);
        return pt;
    }

    //find the largest triangle
    coordinate_type largest_area=-1;
    uint   largest_tri=0;
    for(uint i=0;i<triangulation.size();i++){
        triangle & tri=triangulation[i];
        const point_type& p1=all[tri.v[0]]->getPos();
        const point_type& p2=all[tri.v[1]]->getPos();
        const point_type& p3=all[tri.v[2]]->getPos();
        coordinate_type area=fabs(Area(p1,p2,p3));
        if(area>largest_area){
            largest_area=area;
            largest_tri=i;
        }
    }

    //find a node near the vertex of the triangle
    triangle & tri=triangulation[largest_tri];
    const point_type& p1=(*this)[tri.v[0]]->getPos();
    const point_type& p2=(*this)[tri.v[1]]->getPos();
    const point_type& p3=(*this)[tri.v[2]]->getPos();

    point_type pt;
    bg::assign_point(pt, p1);
    bg::add_point(pt, p2);
    bg::add_point(pt, p3);
    bg::divide_value(pt, 3);

    return pt;
}

template<typename Pt_>
void
c_ply<Pt_>::triangulate(vector<triangle>& tris)
{
     indexing();
     if(triangulation.empty()){

         const point_type& O=getHead()->getPos();

         int * ringVN=new int[1];     //number of vertices for each ring
         assert(ringVN);
         ringVN[0]=getSize();
         int vN=ringVN[0];             //total number of vertices

         if( vN<3 ){
             triangle tri;
             for(short i=0;i<3;i++) tri.v[i]=i;
             tris.push_back(tri);
             return;
         }

         //more than 3 vertices
         int tN=(vN-2);                   //# of triangles
         double * V = new double[vN*2]; //to hold vertices pos
         int *T=new int[3*tN];            //to hold resulting triangles
         assert(T&&V);

         //copy vertices
         {
             int i=0;
             vertex_type * ptr=getHead();
             do{
                 point_type pt=ptr->getPos();
                 V[i*2]=bg::get<0>(pt)-bg::get<0>(O);    // potential implicit conversion
                 V[i*2+1]=bg::get<1>(pt)-bg::get<1>(O);
                 ptr=ptr->getNext();
                 i++;
             }while( ptr!=getHead() );
         }

         FIST_PolygonalArray(1, ringVN, (double (*)[2])V, &tN, (int (*)[3])T);

         for(int i=0;i<tN;i++){
             triangle tri;
             for(int j=0;j<3;j++){
                 tri.v[j]=T[i*3+j];
             }//end j
             triangulation.push_back(tri);
         }//end i


         delete [] ringVN;
         delete [] V;
         delete [] T;
     }//end if

     tris=triangulation;
}

//check if convex
template<typename Pt_>
bool
c_ply<Pt_>::is_convex() const
{
    vertex_type * ptr=head;
    do{
        if(ptr->isReflex()) return false;
        ptr=ptr->getNext();
    }while(ptr!=head); //end while

    return true;
}

template<typename Pt_>
void
c_ply<Pt_>::delete_vertex(vertex_type * v)
{
    vertex_type *pre=v->getPre();
    vertex_type *next=v->getNext();

    pre->setNext(next);
    next->setPre(pre);
    pre->computeExtraInfo(); //recompute info

    if(head==v){
        head=next;
        tail=pre;
    }
    delete v;

    triangulation.clear(); //not valid anymore
    all.clear(); //not valid anymore
}

template<typename Pt_>
void
c_ply<Pt_>::indexing()
{
    uint vid=0;
    all.clear();
    vertex_type * ptr=head;
    do{
        ptr->setVID(vid++);
        all.push_back(ptr);
        ptr=ptr->getNext();
    }while(ptr!=head); //end while
}

template<typename Pt_>
bool
c_ply<Pt_>::identical(const c_ply<Pt_>& other)
{

    //do some basic checks
    if((*this)==other) return true; //same head ptr
    if(type!=other.type) return false;
    if(getSize()!=other.getSize()) return false;
    if( fabs(getArea()-other.getArea())>SMALLNUMBER ) return false;

    //find the first match
    const vertex_type * o_ptr=other.head;
    bool found_first_match=false;
    do{
        const point_type& pos=o_ptr->getPos();
        coordinate_type dx=fabs(bg::get<0>(pos)-bg::get<0>(head->getPos()));
        coordinate_type dy=fabs(bg::get<1>(pos)-bg::get<0>(head->getPos()));

        if( dx<SMALLNUMBER && dy<SMALLNUMBER){
            found_first_match=true;
            break;
        }

        o_ptr=o_ptr->getNext();
    }while(o_ptr!=other.head); //end while

    if(!found_first_match) return false;

    //check the rest of the match
    const vertex_type * ptr=head;

    //since we are sure that ptr matches to o_ptr, we first advance
    ptr=ptr->getNext();
    o_ptr=o_ptr->getNext();

    do
    {
        const point_type& o_pos=o_ptr->getPos();
        const point_type& pos=ptr->getPos();

        coordinate_type dx=fabs(bg::get<0>(pos)-bg::get<0>(o_pos));
        coordinate_type dy=fabs(bg::get<1>(pos)-bg::get<1>(o_pos));

        if( dx>SMALLNUMBER || dy>SMALLNUMBER){
            return false; //don't match
        }

        //advance
        ptr=ptr->getNext();
        o_ptr=o_ptr->getNext();

    }
    while(ptr!=head);

    return true;
}

//
//
// Compute the center and the box of a list of plys
//
//
template<typename Pt_>
void
c_plylist<Pt_>::buildBoxAndCenter()
{
    box[0] = box[2] = std::numeric_limits<coordinate_type>::max();
    box[0] = box[3] = std::numeric_limits<coordinate_type>::lowest();
    for(auto i=begin();i!=end();i++){

        const vertex_type * ptr=i->getHead();
        do{
            const point_type& p=ptr->getPos();
            if(bg::get<0>(p)<box[0]) box[0]= bg::get<0>(p);
            if(bg::get<0>(p)>box[1]) box[1]= bg::get<0>(p);
            if(bg::get<1>(p)<box[2]) box[2]= bg::get<1>(p);
            if(bg::get<1>(p)>box[3]) box[3]= bg::get<1>(p);
            ptr=ptr->getNext();
        }
        while(ptr!=i->getHead()); //end while
    }

    bg::set<0>(center, (box[0]+box[1])/2);
    bg::set<1>(center, (box[2]+box[3])/2);

    is_buildboxandcenter_called=true;
}

template<typename Pt_>
void
c_plylist<Pt_>::translate(const typename c_plylist<Pt_>::point_type& p)
{
    for (auto i=begin();i!=end();i++) i->translate(p);
}

template<typename Pt_>
void
c_polygon<Pt_>::reverse()
{
    for(auto i=begin();i!=end();i++) i->reverse();
    triangulation.clear();
    all.clear();
}

template<typename Pt_>
bool
c_polygon<Pt_>::valid() const
{
    if(empty()) return false;
    if(front().getType()!=ring_type::POUT) return false;
    for(auto i=++begin();i!=end();i++)
      if(i->getType()!=ring_type::PIN)
        return false;

    return true;
}

//copy from the given polygon
template<typename Pt_>
void
c_polygon<Pt_>::copy(const c_polygon& other)
{
    destroy();

    for(auto i=other.cbegin();i!=other.cend();i++){
        ring_type p(ring_type::UNKNOWN);
        p.copy(*i);
        push_back(p);
    }
    indexing();
}

//check if a point is enclosed
//the behavior is unknown if pt is on the boundary of the polygon
template<typename Pt_>
bool
c_polygon<Pt_>::enclosed(const typename c_polygon<Pt_>::point_type& p) const
{
    //find the largest triangle
    for(uint i=0;i<triangulation.size();i++){
        const triangle & tri=triangulation[i];
        const point_type& p1=(*this)[tri.v[0]]->getPos();
        const point_type& p2=(*this)[tri.v[1]]->getPos();
        const point_type& p3=(*this)[tri.v[2]]->getPos();
        coordinate_type area1=Area(p1,p2,p);
        coordinate_type area2=Area(p2,p3,p);
        coordinate_type area3=Area(p3,p1,p);
        if(area1>=0 && area2>=0 && area3>=0) return true; //in
        if(area1<=0 && area2<=0 && area3<=0) return true; //in
    }

    return false; //out
}

template<typename Pt_>
typename c_polygon<Pt_>::point_type
c_polygon<Pt_>::findEnclosedPt()
{
    if(triangulation.empty())
        triangulate(triangulation);

    if(triangulation.empty())
        return front().findEnclosedPt();

    //find the largest triangle
    coordinate_type largest_area=-1;
    uint   largest_tri=0;
    for(uint i=0;i<triangulation.size();i++){
        triangle & tri=triangulation[i];
        const point_type& p1=(*this)[tri.v[0]]->getPos();
        const point_type& p2=(*this)[tri.v[1]]->getPos();
        const point_type& p3=(*this)[tri.v[2]]->getPos();
        coordinate_type area=fabs(Area(p1,p2,p3));
        if(area>largest_area){
            largest_area=area;
            largest_tri=i;
        }
    }

    //find a node near the vertex of the triangle
    triangle & tri=triangulation[largest_tri];
    const point_type& p1=(*this)[tri.v[0]]->getPos();
    const point_type& p2=(*this)[tri.v[1]]->getPos();
    const point_type& p3=(*this)[tri.v[2]]->getPos();

    point_type pt;
    bg::assign_point(pt, p1);
    bg::add_point(pt, p2);
    bg::add_point(pt, p3);
    bg::divide_value(pt, 3);

    return pt;
}

template<typename Pt_>
void
c_polygon<Pt_>::triangulate(vector<triangle>& tris)
{
    if(triangulation.empty())
    {
         const point_type& O=front().getHead()->getPos();

         int ringN=size();             //number of rings
         int * ringVN=new int[ringN];     //number of vertices for each ring
         assert(ringVN);

         int vN=0;             //total number of vertices
         {
             int i=0;
             for(auto ip=begin();ip!=end();ip++,i++){
                 vN+=ip->getSize();
                 ringVN[i]=ip->getSize();
             }
         }

         if( vN<3 ){
//             triangle tri;
//             for(short i=0;i<3;i++) tri.v[i]=i;
//             tris.push_back(tri);
             return;
         }

         int tN=(vN-2)+2*(ringN-1);       //# of triangles, (n-2)+2*(#holes)
         double * V=new double[vN*2];     //to hold vertices pos
         int *T=new int[3*tN];            //to hold resulting triangles
         assert(T&&V);

         //copy vertices
         int i=0;
         for(auto ip=begin();ip!=end();ip++){
             vertex_type * ptr=ip->getHead();
             do{
                 point_type pt=ptr->getPos();
                 V[i*2]=bg::get<0>(pt)-bg::get<0>(O);    // potential implicit conversion to double
                 V[i*2+1]=bg::get<1>(pt)-bg::get<1>(O);
                 ptr=ptr->getNext();
                 i++;
             }while( ptr!=ip->getHead() );
         }

         FIST_PolygonalArray(ringN, ringVN,
             (double (*)[2])V, &tN, (int (*)[3])T);

         for(int i=0;i<tN;i++){
             triangle tri;
             for(int j=0;j<3;j++){
                 tri.v[j]=T[i*3+j];
             }//end j
             triangulation.push_back(tri);
         }//end i

    }

    tris=triangulation;
}

template<typename Pt_>
void
c_polygon<Pt_>::pop_front()
{
  ring_type &ply = front();
  ply.indexing(); // reset vertex IDs
  c_plylist<Pt_>::pop_front();
  /* These must be recomputed */
  all.clear();
  triangulation.clear();
  indexing();
}

template<typename Pt_>
void
c_polygon<Pt_>::destroy()
{
    for(auto i=begin();i!=end();i++){
        i->destroy();
    }
    clear(); //remove all ply from this list
    all.clear();
    triangulation.clear();
}

template<typename Pt_>
bool
c_polygon<Pt_>::is_convex() const
{
    if(size()>1) return false; //contains hole
    return front().is_convex();
}

template<typename Pt_>
void
c_polygon<Pt_>::indexing()
{
    uint vid=0;
    for(auto i=begin();i!=end();i++){
        uint vsize=i->getSize();
        for(uint j=0;j<vsize;j++){
            (*i)[j]->setVID(vid++);       //id of vertex in the polygons
            all.push_back((*i)[j]);
        }
    }//end for i
}

template<typename Pt_>
typename c_polygon<Pt_>::coordinate_type
c_polygon<Pt_>::getArea()
{
  if (area == 0) {
    area = bp::area(*this);
  }

  return area;
}


template<typename Pt_>
bool
c_polygon<Pt_>::identical(c_polygon& p)
{
    if(size()!=p.size()) return false;

    coordinate_type a1=getArea();
    coordinate_type a2=p.getArea();
    if( fabs(a1-a2)>SMALLNUMBER ) return false;



    for(auto i=begin();i!=end();i++)
    {
        bool match=false;
        for(auto j=p.begin();j!=p.end();j++){
            if( i->identical(*j) ){
                match=true;
                break;
            }
        }
        if(!match) return false;
    }

    return true;
}

template<typename Pt_>
istream& operator>>( istream& is, c_ply<Pt_>& poly)
{
    typedef typename c_ply<Pt_>::coordinate_type coordinate_type;
    typedef typename c_ply<Pt_>::point_type point_type;

    int vsize; string str_type;
    is>>vsize>>str_type;

    if( str_type.find("out")!=string::npos )
        poly.type=c_ply<Pt_>::POUT;
    else poly.type=c_ply<Pt_>::PIN;

    poly.beginPoly();
    //read in all the vertices
    int iv;
    vector<point_type> pts; pts.reserve(vsize);
    for( iv=0;iv<vsize;iv++ ){
        coordinate_type x,y;
        is>>x>>y;
        pts.emplace_back(x,y);
        //coordinate_type d=x*x+y*y;
    }
    int id;
    for( iv=0;iv<vsize;iv++ ){
        is>>id; id=id-1;
        poly.addVertex(bg::get<0>(pts[id]), bg::get<1>(pts[id]));
    }

    poly.endPoly();
    return is;
}

template<typename Pt_>
istream& operator>>( istream& is, c_plylist<Pt_>& p)
{
    //remove header commnets
    do{
        char tmp[1024];
        char c=is.peek();
        if(isspace(c)) is.get(c); //eat it
        else if(c=='#') {
            is.getline(tmp,1024);
        }
        else break;
    }while(true);

    //start reading
    uint size;
    is>>size;
    uint vid=0;
    for(uint i=0;i<size;i++){
        c_ply<Pt_> poly(c_ply<Pt_>::UNKNOWN);
        is>>poly;
        p.push_back(poly);
        uint vsize=poly.getSize();
        for(uint j=0;j<vsize;j++){
            poly[j]->setVID(vid++);       //id of vertex in the polygons
            //p.all.push_back(poly[j]);
        }
    }
    return is;
}

template<typename Pt_>
ostream&
operator<<( ostream& os, const c_ply<Pt_>& p)
{
    os<<p.getSize()<<" "<<((p.type==c_ply<Pt_>::PIN)?"PIN  ":"POUT ");
    const typename c_ply<Pt_>::vertex_type * ptr=p.head;
    os << "< " << ptr->getPos();
    ptr = ptr->getNext();
    while(ptr!=p.head)
    {
        os<< ", " << ptr->getPos();
        ptr=ptr->getNext();
    }

    return os;
}

template<typename Pt_>
ostream&
operator<<( ostream& out, const c_plylist<Pt_>& p)
{
    out<<p.size()<<"\n";
    for(auto i=p.cbegin();i!=p.cend();i++) out<<*i;
    return out;
}

template<typename Pt_>
ostream &
operator<<(ostream &os, const c_polygon<Pt_> &p)
{
  unsigned int idx = 0u;
  os << "c_polygon({" << endl;
  for (auto ply_it = p.cbegin(); ply_it != p.cend(); ++ply_it)
    os << "    [" << setw(2) << idx++ << "] " << *ply_it << endl;
  os << "})";
  return os;
}

template<typename Pt_>
ostream &
operator<<(ostream &os, const vector<c_polygon<Pt_> > &v)
{
  os << "vector<c_polygon>({" << endl;
  unsigned int idx = 0u;
  for (auto it = v.begin(); it != v.end(); ++it)
    os << "  [" << setw(2) << idx++ << "] " << *it << endl;
  os << "})";
  return os;
}

template<typename Tp_>
ostream & operator<<(ostream &os,
    boost::geometry::model::d2::point_xy<Tp_> const&p)
{
  os << "(" << bg::get<0>(p) << "," << bg::get<1>(p) << ")";
  return os;
}

PLY_INSTANTIATE(boost::geometry::model::d2::point_xy<float>, );
PLY_INSTANTIATE(boost::geometry::model::d2::point_xy<double>, );

template std::ostream &operator<< <float>(std::ostream&,
    boost::geometry::model::d2::point_xy<float> const&);
template std::ostream &operator<< <double>(std::ostream&,
    boost::geometry::model::d2::point_xy<double> const&);
