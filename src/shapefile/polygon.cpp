//------------------------------------------------------------------------------
//  Copyright 2007-2011 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "polygon.h"
#include "intersection.h"
#include <vector>

using namespace cv;

#ifdef WIN32
extern "C"{
#include "triangulate.h"
}
#else
#include "triangulate.h"
#endif

ply_vertex::~ply_vertex()
{
    //doing nothing for now
}

// - compute normal
// - check if the vertex is reflex or not
void ply_vertex::computeExtraInfo()
{
    //compute normal direction
	Vec2d v=next->pos-pos;
    if( v[0]==0 ){
        if(v[1]>0){ normal[0]=1; normal[1]=0; }
        else{ normal[0]=-1; normal[1]=0; }
    }
    else if( v[0]>0 ){
        normal[1]=-1;
        normal[0]=(v[1]/v[0]);
    }
    else{//v[0]<0
        normal[1]=1;
        normal[0]=-(v[1]/v[0]);
    }
    normal=normalize(normal);

    //compute if left or right turn
    Vec2d u=pos-pre->pos;
    float z=u[0]*v[1]-u[1]*v[0];

    if(z<=0) reflex=true;
    else reflex=false;
}

void ply_vertex::negate()
{
    normal=-normal;
    pos.x=-pos.x;
    pos.y=-pos.y;
}

void ply_vertex::reverse()
{
    swap(next,pre);
    //normal=-normal;
    computeExtraInfo();
}

void ply_vertex::copy(ply_vertex * other)
{
    pos=other->pos;
    normal=other->normal;
    reflex=other->reflex;
    vid=other->vid;
}

///////////////////////////////////////////////////////////////////////////////

//copy from the given ply
void c_ply::copy(const c_ply& other)
{
    destroy();//detroy myself first

    ply_vertex* ptr=other.head;
    beginPoly();
    do{
        ply_vertex * v=new ply_vertex();
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
    center=other.center;
    radius=other.radius;
    type=other.type;
    triangulation=other.triangulation;
    extra_info=other.extra_info;
}

// clean up the space allocated
void c_ply::destroy()
{
    if( head==NULL ) return;
    ply_vertex* ptr=head;
    do{
        ply_vertex * n=ptr->getNext();
        delete ptr;
        ptr=n;
    }while( ptr!=head );
    head=tail=NULL;

    all.clear();
    triangulation.clear();
}

// Create a empty polygon
void c_ply::beginPoly()
{
    head=tail=NULL;
    all.clear();
    triangulation.clear();
}

// Add a vertex to the polygonal chian
void c_ply::addVertex( double x, double y, bool remove_duplicate )
{
    Point2d pt(x,y);

    if(tail!=NULL){
        if(tail->getPos()==pt && remove_duplicate) return; //don't add
    }

    ply_vertex * v=new ply_vertex(pt);
    if( tail!=NULL ){
        tail->setNext(v);
    }
    tail=v;
    if( head==NULL ) head=tail;
    v->setVID(all.size()); //id of the vertex in this ply
	all.push_back(v);

}

// Add a vertex to the polygonal chian
void c_ply::addVertex( ply_vertex * v )
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
void c_ply::endPoly(bool remove_duplicate)
{
    if(head!=NULL && tail!=NULL){
        if(remove_duplicate){
            if(head->getPos()==tail->getPos()){ //remove tail..
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
void c_ply::doInit()
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
    ply_vertex* ptr=head;
    do{
        ptr->computeExtraInfo();
        ptr=ptr->getNext();
    }while( ptr!=head );
}

const Point2d& c_ply::getCenter()
{
    if(radius<0){
        center = Point2d(0,0);
        ply_vertex * ptr=head;
        const Point2d& first=ptr->getPos();
        uint size=0;
        do{
            size++;
            Vec2d v=ptr->getPos()-first;
            center.x+=v[0];
            center.y+=v[1];
            ptr=ptr->getNext();
        }while(ptr!=head); //end while
        center.x=(center.x/size)+first.x;
        center.y=(center.y/size)+first.y;

        radius=0;
    }

    return center;
}


///////////////////////////////////////////////////////////////////////////
void c_ply::negate()
{
    ply_vertex * ptr=head;
    do{
        ptr->negate();
        ptr=ptr->getNext();
    }while(ptr!=head); //end while
}

//reverse the order of the vertices
void c_ply::reverse()
{
    ply_vertex * ptr=head;
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
void c_ply::translate(const Point2d& p)
{
    ply_vertex * ptr=head;
    do{
        ptr->translate(p);
        ptr=ptr->getNext();
    }while(ptr!=head); //end while

    //translate box
    extra_info.box[0]+=p.x;
    extra_info.box[1]+=p.x;
    extra_info.box[2]+=p.y;
    extra_info.box[3]+=p.y;
}

///////////////////////////////////////////////////////////////////////////
//compute the Radius of the poly chain
float c_ply::getRadius()
{
    if(radius<0) getCenter();

    if(radius==0){
        ply_vertex * ptr=head;
        do{
            float d=cv::norm(Vec2d(center-ptr->getPos()), cv::NORM_L2SQR);
            if(d>radius) radius=d;
            ptr=ptr->getNext();
        }while(ptr!=head); //end while
        radius=sqrt(radius);
    }

    return radius;
}


float c_ply::getArea()
{
    if(area==-FLT_MAX){
      area = cv::contourArea(all);
    }

    return area;
}


bool c_ply::enclosed(const Point2d& p)
{
    if(triangulation.empty())
        triangulate(triangulation);

    //find the largest triangle
    for(uint i=0;i<triangulation.size();i++){
        triangle & tri=triangulation[i];
        const Point2d& p1=all[tri.v[0]]->getPos();
        const Point2d& p2=all[tri.v[1]]->getPos();
        const Point2d& p3=all[tri.v[2]]->getPos();
        double area1=Area(p1,p2,p);
        double area2=Area(p2,p3,p);
        double area3=Area(p3,p1,p);
        if(area1>=0 && area2>=0 && area3>=0) return true; //in
        if(area1<=0 && area2<=0 && area3<=0) return true; //in
    }

    return false; //out
}

Point2d c_ply::findEnclosedPt()
{
    if(triangulation.empty())
        triangulate(triangulation);

    if(triangulation.empty()){ //too few points...
        const Point2d& p1=head->getPos();
        const Point2d& p2=head->getNext()->getPos();
        Point2d pt;
        pt.x=(p1.x+p2.x)/2;
        pt.y=(p1.y+p2.y)/2;
        return pt;
    }

    //find the largest triangle
    double largest_area=-1;
    uint   largest_tri=0;
    for(uint i=0;i<triangulation.size();i++){
        triangle & tri=triangulation[i];
        const Point2d& p1=all[tri.v[0]]->getPos();
        const Point2d& p2=all[tri.v[1]]->getPos();
        const Point2d& p3=all[tri.v[2]]->getPos();
        double area=fabs(Area(p1,p2,p3));
        if(area>largest_area){
            largest_area=area;
            largest_tri=i;
        }
    }

    //find a node near the vertex of the triangle
    triangle & tri=triangulation[largest_tri];
    const Point2d& p1=(*this)[tri.v[0]]->getPos();
    const Point2d& p2=(*this)[tri.v[1]]->getPos();
    const Point2d& p3=(*this)[tri.v[2]]->getPos();

    Point2d pt;
    pt.x=(p1.x+p2.x+p3.x)/3;
    pt.y=(p1.y+p2.y+p3.y)/3;

    return pt;
}

void c_ply::triangulate(vector<triangle>& tris)
{
     if(triangulation.empty()){

         const Point2d& O=getHead()->getPos();

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
         double * V=new double[vN*2];     //to hold vertices pos
         int *T=new int[3*tN];            //to hold resulting triangles
         assert(T&&V);

         //copy vertices
         {
             int i=0;
             ply_vertex * ptr=getHead();
             do{
                 Point2d pt=ptr->getPos();
                 V[i*2]=pt.x-O.x;
                 V[i*2+1]=pt.y-O.y;
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
bool c_ply::is_convex() const
{
    ply_vertex * ptr=head;
    do{
        if(ptr->isReflex()) return false;
        ptr=ptr->getNext();
    }while(ptr!=head); //end while

    return true;
}

void c_ply::delete_vertex(ply_vertex * v)
{
    ply_vertex *pre=v->getPre();
    ply_vertex *next=v->getNext();

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

void c_ply::indexing()
{
    uint vid=0;
    all.clear();
    ply_vertex * ptr=head;
    do{
        ptr->setVID(vid++);
        all.push_back(ptr);
        ptr=ptr->getNext();
    }while(ptr!=head); //end while
}

bool c_ply::identical(c_ply& other)
{

    //do some basic checks
    if((*this)==other) return true; //same head ptr
    if(type!=other.type) return false;
    if(getSize()!=other.getSize()) return false;
    if( fabs(getArea()-other.getArea())>SMALLNUMBER ) return false;

    //find the first match
    ply_vertex * o_ptr=other.head;
    bool found_first_match=false;
    do{
        const Point2d& pos=o_ptr->getPos();
        double dx=fabs(pos.x-head->getPos().x);
        double dy=fabs(pos.y-head->getPos().y);

        if( dx<SMALLNUMBER && dy<SMALLNUMBER){
            found_first_match=true;
            break;
        }

        o_ptr=o_ptr->getNext();
    }while(o_ptr!=other.head); //end while

    if(!found_first_match) return false;

    //check the rest of the match
    ply_vertex * ptr=head;

    //since we are sure that ptr matches to o_ptr, we first advance
    ptr=ptr->getNext();
    o_ptr=o_ptr->getNext();

    do
    {
        const Point2d& o_pos=o_ptr->getPos();
        const Point2d& pos=ptr->getPos();

        double dx=fabs(pos.x-o_pos.x);
        double dy=fabs(pos.y-o_pos.y);

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

void c_plylist::buildBoxAndCenter()
{
    //typedef list<c_ply>::iterator IT;
    box[0]=box[2]=FLT_MAX;
    box[0]=box[3]=-FLT_MAX;
    for(iterator i=begin();i!=end();i++){

        ply_vertex * ptr=i->getHead();
        do{
            const Point2d& p=ptr->getPos();
            if(p.x<box[0]) box[0]=p.x;
            if(p.x>box[1]) box[1]=p.x;
            if(p.y<box[2]) box[2]=p.y;
            if(p.y>box[3]) box[3]=p.y;
            ptr=ptr->getNext();
        }
        while(ptr!=i->getHead()); //end while
    }

    center.x=(box[0]+box[1])/2;
    center.y=(box[2]+box[3])/2;

    is_buildboxandcenter_called=true;
}

//
// Compute the center and the box of a list of plys
//
void c_plylist::translate(const Point2d& p)
{
    for(iterator i=begin();i!=end();i++) i->translate(p);
}

void c_polygon::reverse()
{
    for(iterator i=begin();i!=end();i++) i->reverse();
    triangulation.clear();
    all.clear();
}

bool c_polygon::valid() //check if this is a valid polygon
{
    //typedef list<c_ply>::iterator IT;
    if(empty()) return false;
    if(front().getType()!=c_ply::POUT) return false;
    for(iterator i=++begin();i!=end();i++) if(i->getType()!=c_ply::PIN) return false;

    return true;
}

//copy from the given polygon
void c_polygon::copy(const c_polygon& other)
{
    destroy();

    for(const_iterator i=other.begin();i!=other.end();i++){
        c_ply p(c_ply::UNKNOWN);
        p.copy(*i);
        push_back(p);
    }
}

//check if a point is enclosed
//the behavior is unknown if pt is on the boundary of the polygon
bool c_polygon::enclosed(const Point2d& p) const
{
    //find the largest triangle
    for(uint i=0;i<triangulation.size();i++){
        const triangle & tri=triangulation[i];
        const Point2d& p1=(*this)[tri.v[0]]->getPos();
        const Point2d& p2=(*this)[tri.v[1]]->getPos();
        const Point2d& p3=(*this)[tri.v[2]]->getPos();
        double area1=Area(p1,p2,p);
        double area2=Area(p2,p3,p);
        double area3=Area(p3,p1,p);
        if(area1>=0 && area2>=0 && area3>=0) return true; //in
        if(area1<=0 && area2<=0 && area3<=0) return true; //in
    }

    return false; //out
}

Point2d c_polygon::findEnclosedPt()
{
    if(triangulation.empty())
        triangulate(triangulation);

    if(triangulation.empty())
        return front().findEnclosedPt();

    //find the largest triangle
    double largest_area=-1;
    uint   largest_tri=0;
    for(uint i=0;i<triangulation.size();i++){
        triangle & tri=triangulation[i];
        const Point2d& p1=(*this)[tri.v[0]]->getPos();
        const Point2d& p2=(*this)[tri.v[1]]->getPos();
        const Point2d& p3=(*this)[tri.v[2]]->getPos();
        double area=fabs(Area(p1,p2,p3));
        if(area>largest_area){
            largest_area=area;
            largest_tri=i;
        }
    }

    //find a node near the vertex of the triangle
    triangle & tri=triangulation[largest_tri];
    const Point2d& p1=(*this)[tri.v[0]]->getPos();
    const Point2d& p2=(*this)[tri.v[1]]->getPos();
    const Point2d& p3=(*this)[tri.v[2]]->getPos();

    Point2d pt;
    pt.x=(p1.x+p2.x+p3.x)/3;
    pt.y=(p1.y+p2.y+p3.y)/3;

    return pt;
}

void c_polygon::triangulate(vector<triangle>& tris)
{
    if(triangulation.empty())
    {
         const Point2d& O=front().getHead()->getPos();
         typedef iterator   PIT;

         int ringN=size();             //number of rings
         int * ringVN=new int[ringN];     //number of vertices for each ring
         assert(ringVN);

         int vN=0;             //total number of vertices
         {
             int i=0;
             for(PIT ip=begin();ip!=end();ip++,i++){
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
         for(PIT ip=begin();ip!=end();ip++){
             ply_vertex * ptr=ip->getHead();
             do{
                 Point2d pt=ptr->getPos();
                 V[i*2]=pt.x-O.x;
                 V[i*2+1]=pt.y-O.y;
                 ptr=ptr->getNext();
                 i++;
             }while( ptr!=ip->getHead() );
         }

         FIST_PolygonalArray(ringN, ringVN, (double (*)[2])V, &tN, (int (*)[3])T);

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

void c_polygon::pop_front()
{
  c_ply &ply = front();
  ply.indexing(); // reset vertex IDs
  c_plylist::pop_front();
  /* These must be recomputed */
  all.clear();
  triangulation.clear();
  indexing();
}

void c_polygon::destroy()
{
    for(iterator i=begin();i!=end();i++){
        i->destroy();
    }
    clear(); //remove all ply from this list
    all.clear();
    triangulation.clear();
}

bool c_polygon::is_convex() const
{
    if(size()>1) return false; //contains hole
    return front().is_convex();
}

void c_polygon::indexing()
{
    uint vid=0;
    for(iterator i=begin();i!=end();i++){
        uint vsize=i->getSize();
        for(uint j=0;j<vsize;j++){
            (*i)[j]->setVID(vid++);       //id of vertex in the polygons
            all.push_back((*i)[j]);
        }
    }//end for i
}

double c_polygon::getArea()
{
    if(area==0){
        for(iterator i=begin();i!=end();i++){
            area+=i->getArea();
        }//end for i
    }

    return area;
}


bool c_polygon::identical(c_polygon& p)
{
    if(this->size()!=p.size()) return false;

    double a1=getArea();
    double a2=p.getArea();
    if( fabs(a1-a2)>SMALLNUMBER ) return false;



    for(iterator i=begin();i!=end();i++)
    {
        bool match=false;
        for(iterator j=p.begin();j!=p.end();j++){
            if( i->identical(*j) ){
                match=true;
                break;
            }
        }
        if(!match) return false;
    }

    return true;
}

istream& operator>>( istream& is, c_ply& poly)
{
    int vsize; string str_type;
    is>>vsize>>str_type;

    if( str_type.find("out")!=string::npos )
        poly.type=c_ply::POUT;
    else poly.type=c_ply::PIN;

    poly.beginPoly();
    //read in all the vertices
    int iv;
    vector< pair<double,double> > pts; pts.reserve(vsize);
    for( iv=0;iv<vsize;iv++ ){
        double x,y;
        is>>x>>y;
        pts.push_back(pair<double,double>(x,y));
        //double d=x*x+y*y;
    }
    int id;
    for( iv=0;iv<vsize;iv++ ){
        is>>id; id=id-1;
        poly.addVertex(pts[id].first,pts[id].second);
    }

    poly.endPoly();
    return is;
}

istream& operator>>( istream& is, c_plylist& p)
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
        c_ply poly(c_ply::UNKNOWN);
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

ostream& operator<<( ostream& os, const c_ply& p)
{
    os<<p.getSize()<<" "<<((p.type==c_ply::PIN)?"in":"out")<<"\n";
    const ply_vertex * ptr=p.head;
    do{
        os<<ptr->getPos().x<<" "<<ptr->getPos().y<<"\n";
        ptr=ptr->getNext();
    }while(ptr!=p.head);

    for(int i=0;i<p.getSize();i++) os<<i+1<<" ";
    os<<"\n";
    return os;
}

ostream& operator<<( ostream& out, const c_plylist& p)
{
    out<<p.size()<<"\n";
    typedef c_plylist::const_iterator PIT;
    for(PIT i=p.cbegin();i!=p.cend();i++) out<<*i;
    return out;
}

std::ostream &operator<<(std::ostream &os, const std::vector<c_polygon> &v)
{
  for (auto it = v.cbegin(); it != v.cend(); ++it)
    os << *it;
  return os;
}
