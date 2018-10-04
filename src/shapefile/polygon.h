//------------------------------------------------------------------------------
//  Copyright 2007-2011 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _POLYGON_DATA_H_
#define _POLYGON_DATA_H_

#ifdef WIN32
#pragma warning(disable : 4786)
#endif

#include <Point.h>
#include <Vector.h>
using namespace mathtool;

#include <list>
#include <vector>
#include <cassert>
#include <float.h>
using namespace std;

typedef unsigned int uint;

//
//a triangle, used in triangulating a polygon
//
struct triangle
{
    uint v[3]; // id to the vertices
};

//
// extra information about the polygon
//

struct ply_extra_info
{
    ply_extra_info()
    {
        base_height=height=0;
        has_box=false;
        box[0]=box[1]=box[2]=box[3]=0;
        flag=UINT_MAX;
    }

    bool has_box;
    double base_height;
    double height;
    double box[4];
    uint flag;
};


//
// Vertex of polygon
//
class ply_vertex
{
public:
    
    ///////////////////////////////////////////////////////////////////////////
    ply_vertex(){ init(); }
    ply_vertex( const Point2d& p ){ pos=p; init(); }
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

    ///////////////////////////////////////////////////////////////////////////
    void setPos(const Point2d& p) { pos=p; }
    virtual const Point2d& getPos() const { return pos; }

    void translate(const Vector2d& v){ pos=pos+v; }

    virtual ply_vertex * getNext() const { return next; }
    virtual ply_vertex * getPre() const { return pre; }

    const Vector2d& getNormal() const { return normal; }
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
    Point2d pos;       //position
    ply_vertex * next; //next vertex in the polygon
    ply_vertex * pre;  //previous vertex in the polygon
    Vector2d normal;   //normal, the segment normal from this v to the next.
    bool reflex;
    uint vid;
};

//
// Polygon chain
//
class c_ply{
public:

    enum POLYTYPE { UNKNOWN, PIN, POUT };

    ///////////////////////////////////////////////////////////////////////////
    c_ply(POLYTYPE t){ head=tail=NULL; type=t; radius=-1; area=-FLT_MAX; }

    ///////////////////////////////////////////////////////////////////////////
    void copy(const c_ply& ply); //copy from the other ply
    void destroy();

    ///////////////////////////////////////////////////////////////////////////
    // create c_ply
    void beginPoly();
    void addVertex( double x, double y, bool remove_duplicate=false );
    void addVertex( ply_vertex * v );
    void endPoly(bool remove_duplicate=false);

    ///////////////////////////////////////////////////////////////////////////
    void negate();
    void reverse(); //reverse vertex order

    ///////////////////////////////////////////////////////////////////////////
    void translate(const Vector2d& v);

    ///////////////////////////////////////////////////////////////////////////
    //
    Point2d findEnclosedPt(); //find a point that is enclosed by this polychain

    //triangulate the polygon
    void triangulate(vector<triangle>& tris);

    ///////////////////////////////////////////////////////////////////////////
    // Access functions
    ply_vertex * getHead() const { return head; }

    //hole or external boundary
    POLYTYPE getType() const { return type; }

    void set(POLYTYPE t,ply_vertex * h){ 
        type=t; head=h; 
        if(h!=NULL){ tail=h->getPre(); }
        else{ tail=NULL; }
    }

    int getSize() {
        if(all.empty()) indexing();
        return all.size();
    }

    ply_vertex * operator[](unsigned int id){
        if(all.empty()) indexing();
        return all[id];
    }

    //get com //center of mass
    const Point2d& getCenter();

	//compute the Radius of the poly chain
	float getRadius();
	
	//area
	float getArea();

    //
	// additional functions
	//

	//check if a point is enclosed
	//the behavior is unknown if pt is on the boundary of the polygon
	bool enclosed(const Point2d& pt);

	//check if convex
	bool is_convex() const;

	//delete a vertex
	void delete_vertex(ply_vertex * p);

	//get extra info
	ply_extra_info& extra(){ return extra_info; }
	const ply_extra_info& extra() const { return extra_info; }

    ///////////////////////////////////////////////////////////////////////////
    // Operator
    //check if give poly line is the same as this
    bool operator==( const c_ply& other ){ return other.head==head; }

    //
    //this is a more general check if other and this ply are identical
    //
    bool identical(c_ply& other);

    friend istream& operator>>( istream&, c_ply& );
    friend ostream& operator<<( ostream&, c_ply& );

protected:

    ///////////////////////////////////////////////////////////////////////////
    void doInit(); /*return # of vertice in this poly*/

    //indexing all elements and store in vector<ply_vertex*>
    void indexing();

private:

    ply_vertex * head; //the head of vertex list
    ply_vertex * tail; //end of the vertex list

	vector<ply_vertex*> all; //all vertices
	
	//additional info
	Point2d center;
	float radius;
	float area;

    //In, out or unknown.
    POLYTYPE type;

    //extrac info
    ply_extra_info extra_info;

    //triangulation
    vector<triangle> triangulation; //catched triangulation, calculated by triangulate
};


//a c_plylist is a list of c_ply
class c_plylist : public list<c_ply>
{
    friend ostream& operator<<( ostream&, c_plylist& );
    friend istream& operator>>( istream&, c_plylist& );

public:

    c_plylist()
    {
        box[0]=box[1]=box[2]=box[3]=0;
        is_buildboxandcenter_called=false;
    }

    void translate(const Vector2d& v);

    //access
    void buildBoxAndCenter();
    double * getBBox() { assert(is_buildboxandcenter_called); return box; }
    const Point2d& getCenter() { assert(is_buildboxandcenter_called); return center; }

protected:

    Point2d center;
    double box[4];

private:

    bool is_buildboxandcenter_called;
};

//
// a c_polygon is a restricted kind of c_plylist
// this defines a simple polygon so that
// the first element much be a POUT c_ply and
// the rest ply lines are holes
//
class c_polygon : public c_plylist
{
public:

    c_polygon() { area=0; }

    bool valid(); //check if this is a valid polygon

    //copy from the given polygon
    void copy(const c_polygon& other);

    //triangulate the polygon
    void triangulate(vector<triangle>& tris);

    void reverse(); //reverse the vertex order (not the list order)

    //access the vertices of the polygon as an array
    int getSize()
    {
        if(all.empty()) indexing();
        return all.size();
    }

    ply_vertex * operator[](unsigned int id){
        if(all.empty()) indexing();
        return all[id];
    }

    double getArea();

    //destroy
    void destroy();

    //check if a point is enclosed
    //the behavior is unknown if pt is on the boundary of the polygon
    bool enclosed(const Point2d& pt);

    //find a point inside the polygon
    Point2d findEnclosedPt();

    //get number of vertices
    uint getSize() const;

    bool is_convex() const;

    //check if polygons are identical (or not)
    bool identical(c_polygon& p);

private:

    //indexing the vertices and store them in vector<ply_vertex*> all
    void indexing();

    vector<ply_vertex*> all; //all vertices

    //triangulation
    vector<triangle> triangulation; //catched triangulation, calculated by triangulate

    float area;
};

//
// Given radius and resolution (# of vertices), create a circle in polygon format
//
inline void create_circle(c_polygon& p, double radius, uint res)
{
    double delta=PI2/res;
    c_ply ply(c_ply::POUT);
    ply.beginPoly();
    for(uint i=0;i<res;i++){
        double r=delta*i;
        double x=cos(r)*radius;
        double y=sin(r)*radius;
        ply.addVertex(x,y);
    }//end for i
    ply.endPoly();
    p.push_back(ply);
}


#endif //_POLYGON_H_


