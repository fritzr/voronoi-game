#ifndef _MKSUM_INTERSECTION_H_
#define _MKSUM_INTERSECTION_H_

#include <iostream>
#include <cassert>

#include "boost_geo_poly.h"

/* range of real numbers */
#define SMALLNUMBER 1.0e-10
#define HUGENUMBER  1.0e10

/* Whether p1 -> p2 -> p3 forms a left turn.  */
template<class Pt_>
bool
leftTurn(Pt_ const& p1, Pt_ const& p2, Pt_ const& p3)
{
  auto v = p3;
  bg::subtract_point(v, p2); // normal direction
  auto u = p2;
  bg::subtract_point(u, p1);
  auto z = u.x() * v.y() - u.y() * v.x();
  return z > 0;
}

// This is copied from CG in C

template <typename Pt>
inline typename bg::coordinate_type<Pt>::type
Area(Pt a, Pt b, Pt c)
{
  // (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y)
  return ( bg::get<0>(b) - bg::get<0>(a) ) * ( bg::get<1>(c) - bg::get<1>(a) ) -
         ( bg::get<0>(c) - bg::get<0>(a) ) * ( bg::get<1>(b) - bg::get<1>(a) );
}

template <typename Pt>
inline int AreaSign(Pt a, Pt b, Pt c)
{
    auto area=Area(a,b,c);
    if      ( area >  SMALLNUMBER ) return  1;
    else if ( area < -SMALLNUMBER ) return -1;
    else     return  0;
}

template <typename Pt>
inline int Collinear(const Pt a, const Pt b, const Pt c)
{
   return AreaSign( a, b, c ) == 0;
}

/*---------------------------------------------------------------------
Returns TRUE iff point c lies on the closed segement ab.
Assumes it is already known that abc are collinear.
---------------------------------------------------------------------*/
template <typename Pt>
inline bool Between(const Pt a, const Pt b, const Pt c)
{
  // If ab not vertical, check betweenness on x; else on y.
  if ( fabs(bg::get<0>(a)-bg::get<0>(b))>fabs(bg::get<1>(a)-bg::get<1>(b)) )
    return
      ((bg::get<0>(a) <= bg::get<0>(c)) && (bg::get<0>(c) <= bg::get<0>(b)))
   || ((bg::get<0>(a) >= bg::get<0>(c)) && (bg::get<0>(c) >= bg::get<0>(b)));
  else
    return
      ((bg::get<1>(a) <= bg::get<1>(c)) && (bg::get<1>(c) <= bg::get<1>(b)))
   || ((bg::get<1>(a) >= bg::get<1>(c)) && (bg::get<1>(c) >= bg::get<1>(b)));
}

template <typename Pt>
inline bool Between_strict(const Pt a, const Pt b, const Pt c)
{
   // If ab not vertical, check betweenness on x; else on y.
    if ( fabs(bg::get<0>(a)-bg::get<0>(b))>SMALLNUMBER ){
      auto c01=bg::get<0>(c)-SMALLNUMBER;
      auto c02=bg::get<0>(c)+SMALLNUMBER;

      return ((bg::get<0>(a) < c01) && (c02 < bg::get<0>(b))) ||
             ((bg::get<0>(a) > c02) && (c01 > bg::get<0>(b)));
    }
    else{
        auto c11=bg::get<1>(c)-SMALLNUMBER;
        auto c10=bg::get<1>(c)+SMALLNUMBER;
        return ((bg::get<1>(a) < c11) && (c10 < bg::get<1>(b))) ||
               ((bg::get<1>(a) > c10) && (c11 > bg::get<1>(b)));
    }
}


template <typename C>
inline bool AlmostEqual3(const C a, const C b)
{
    return (fabs(a[0]-b[0])<SMALLNUMBER &&fabs(a[1]-b[1])<SMALLNUMBER && fabs(a[2]-b[2])<SMALLNUMBER);
}

template <typename Pt, typename T>
inline bool AlmostEqual(const Pt a, const Pt b, T tau)
{
    return (fabs(bg::get<0>(a)-bg::get<0>(b))<tau &&fabs(bg::get<1>(a)-bg::get<1>(b))<tau);
}

template <typename Pt>
inline bool AlmostEqual(const Pt a, const Pt b)
{
    return AlmostEqual(a, b, SMALLNUMBER);
}


// compute union of two colinear  segments ab and cd
// place the union in p
// return false if the union is degenerated
template <typename Pt>
inline bool Union(Pt const& a, Pt const& b, Pt const& c, Pt const& d, Pt *p)
{
    int id=0;
    if(AlmostEqual(a,c)){ p[id]=a; id++; }

    if(AlmostEqual(a,d)){ p[id]=a; id++; }
    if(id==2) return bg::get<0>(p)!=bg::get<1>(p);

    if(AlmostEqual(b,c)){ p[id]=b; id++; }
    if(id==2) return bg::get<0>(p)!=bg::get<1>(p);

    if(AlmostEqual(b,d)){ p[id]=b; id++; }
    if(id==2) return bg::get<0>(p)!=bg::get<1>(p);

    if( Between_strict(a,b,c) ){ p[id]=c; id++; }
    if(id==2) return bg::get<0>(p)!=bg::get<1>(p);

    if( Between_strict(a,b,d) ){ p[id]=d; id++; }
    if(id==2) return bg::get<0>(p)!=bg::get<1>(p);

    if( Between_strict(c,d,a) ){ p[id]=a; id++; }
    if(id==2) return bg::get<0>(p)!=bg::get<1>(p);

    if( Between_strict(c,d,b) ){ p[id]=b; id++; }
    if(id==2) return bg::get<0>(p)!=bg::get<1>(p);

    return false;
}


template <typename Pt>
inline char ParallelInt(
    Pt const& a, Pt const& b, Pt const& c, Pt const& d, Pt &p)
{
   if(!Collinear(a, b, c)) return '0';

   //check if they overlap..
   if(Between(a,b,c)) return 'e';
   if(Between(a,b,d)) return 'e';
   if(Between(c,d,a)) return 'e';
   if(Between(c,d,b)) return 'e';

   //they don't overlap but the end points may..
   //check if the end points overlap
   if(AlmostEqual(a,c)){ bg::assign_point(p, a); return 'v'; }
   if(AlmostEqual(b,c)){ bg::assign_point(p, b); return 'v'; }
   if(AlmostEqual(a,d)){ bg::assign_point(p, a); return 'v'; }
   if(AlmostEqual(b,d)){ bg::assign_point(p, b); return 'v'; }

   return '0';
}

/*---------------------------------------------------------------------
SegSegInt: Finds the point of intersection p between two closed
segments ab and cd.  Returns p and a char with the following meaning:
   'e': The segments collinearly overlap, sharing a point.
   'v': An endpoint (vertex) of one segment is on the other segment,
        but 'e' doesn't hold.
   '1': The segments intersect properly (i.e., they share a point and
        neither 'v' nor 'e' holds).
   '0': The segments do not intersect (i.e., they share no points).
Note that two collinear segments that share just one point, an endpoint
of each, returns 'e' rather than 'v' as one might expect.
---------------------------------------------------------------------*/
template <typename Pt>
inline char SegSegInt(
    Pt const& a, Pt const& b, Pt const& c, Pt const& d, Pt &p)
{
  typedef typename bg::coordinate_type<Pt>::type coord_t;
   coord_t  s, t;        // The two parameters of the parametric eqns.
   coord_t  num_s, num_t, denom;   // Numerator and denoninator of equations.
   char    code = '?';    // Return char characterizing intersection.
   const coord_t TINYNUMBER=SMALLNUMBER;

   denom = bg::get<0>(a) * ( bg::get<1>(d) - bg::get<1>(c) ) +
           bg::get<0>(b) * ( bg::get<1>(c) - bg::get<1>(d) ) +
           bg::get<0>(d) * ( bg::get<1>(b) - bg::get<1>(a) ) +
           bg::get<0>(c) * ( bg::get<1>(a) - bg::get<1>(b) );

   // If denom is zero, then segments are parallel: handle separately.
   //if (fabs(denom)<SMALLNUMBER) denom=0;
   if (fabs(denom)<TINYNUMBER) return  ParallelInt(a, b, c, d, p);


   if(AlmostEqual(a,c)){ bg::assign_point(p, a); return 'v'; }
   if(AlmostEqual(b,c)){ bg::assign_point(p, b); return 'v';}
   if(AlmostEqual(a,d)){ bg::assign_point(p, a); return 'v';}
   if(AlmostEqual(b,d)){ bg::assign_point(p, b); return 'v';}

   //compute s
   num_s =    bg::get<0>(a) * ( bg::get<1>(d) - bg::get<1>(c) ) +
              bg::get<0>(c) * ( bg::get<1>(a) - bg::get<1>(d) ) +
              bg::get<0>(d) * ( bg::get<1>(c) - bg::get<1>(a) );


   s = num_s / denom;

   if( fabs(s)<TINYNUMBER ) s=0;
   else if( fabs(1-s)<TINYNUMBER ) s=1;


   //compute t
   num_t = -( bg::get<0>(a) * ( bg::get<1>(c) - bg::get<1>(b) ) +
              bg::get<0>(b) * ( bg::get<1>(a) - bg::get<1>(c) ) +
              bg::get<0>(c) * ( bg::get<1>(b) - bg::get<1>(a) ) );



   t = num_t / denom;

   if( fabs(t)<TINYNUMBER ) t=0;
   else if( fabs(1-t)<TINYNUMBER) t=1;

   //decide the code
//   if(s==0 || s==1 || t==0 || t==1)
//        code = 'v';
   if( (0.0<s) && (s< 1.0) && (0.0< t) && (t< 1.0) )
        code = '1';
   else if ( (0.0>s) || (s>1.0) || (0.0>t) || (t>1.0) )
        return '0';
   else
       code= 'v';

   if(code!='v'){
     bg::set<0>(p, (bg::get<0>(a) + s*(bg::get<0>(b)-bg::get<0>(a))));
     bg::set<1>(p, (bg::get<1>(a) + s*(bg::get<1>(b)-bg::get<1>(a))));
   }
   else{
      if(s==0){ bg::assign_point(p, a); }
      else if(s==1){ bg::assign_point(p, b); }
      else if(t==0){ bg::assign_point(p, c); }
      else if(t==1){ bg::assign_point(p, d); }
      else{
        std::cerr<<"s="<<s<<" t="<<t<<std::endl;
        std::cerr<<"a="<<bg::get<0>(a)<<","<<bg::get<1>(a)<<std::endl;
        std::cerr<<"b="<<bg::get<0>(b)<<","<<bg::get<1>(b)<<std::endl;
        std::cerr<<"c="<<bg::get<0>(c)<<","<<bg::get<1>(c)<<std::endl;
        std::cerr<<"d="<<bg::get<0>(d)<<","<<bg::get<1>(d)<<std::endl;
        assert(false);
      }
   }

   //T check...
//   if(code=='v'){
//       T q;
//       q[0] = (c[0] + t*(d[0]-c[0]));
//       q[1] = (c[1] + t*(d[1]-c[1]));
//       if(AlmostEqual(p,q)==false){
//           //cout<<"SHOOT: "<<code<<" points:"<<p[0]<<","<<p[1]<<"  "<<q[0]<<","<<q[1]<<endl;
//           return '0';
//       }
//   }

   return code;
}



////////////////////////////////////////////////////////////////////////////////
// This is from RAPID and should be simplified


template <typename T>
inline T
VdotV(T V1[3], T V2[3])
{
  return (V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2]);
}



template <typename T>
inline void
VcrossV(T Vr[3], const T V1[3], const T V2[3])
{
  Vr[0] = V1[1]*V2[2] - V1[2]*V2[1];
  Vr[1] = V1[2]*V2[0] - V1[0]*V2[2];
  Vr[2] = V1[0]*V2[1] - V1[1]*V2[0];
}

template <typename T>
inline T
max(T a, T b, T c)
{
  T t = a;
  if (b > t) t = b;
  if (c > t) t = c;
  return t;
}

template <typename T>
inline T
min(T a, T b, T c)
{
  T t = a;
  if (b < t) t = b;
  if (c < t) t = c;
  return t;
}


template <typename T>
inline int
my_project6_2(T *ax, T *p1, T *p2, T *p3, T *q1, T *q2, T *q3) {
  T P1 = VdotV(ax, p1);
  T P2 = VdotV(ax, p2);
  T P3 = VdotV(ax, p3);
  T Q1 = VdotV(ax, q1);
  T Q2 = VdotV(ax, q2);
  T Q3 = VdotV(ax, q3);

  T mx1 = max(P1, P2, P3);
  T mn1 = min(P1, P2, P3);
  T mx2 = max(Q1, Q2, Q3);
  T mn2 = min(Q1, Q2, Q3);

  if (mn1 > mx2) return 0;
  if (mn2 > mx1) return 0;
  return 1;
}


// very robust triangle intersection test
// uses no divisions
// works on coplanar triangles

template <typename T>
inline int
my_tri_contact (const T *P1, const T *P2, const T *P3,
        const T *Q1, const T *Q2, const T *Q3)
{

  /*
     One triangle is (p1,p2,p3).  Other is (q1,q2,q3).
     Edges are (e1,e2,e3) and (f1,f2,f3).
     Normals are n1 and m1
     Outwards are (g1,g2,g3) and (h1,h2,h3).

     We assume that the triangle vertices are in the same coordinate system.

     First thing we do is establish a new c.s. so that p1 is at (0,0,0).
     */

  T p1[3], p2[3], p3[3];
  T q1[3], q2[3], q3[3];
  T e1[3], e2[3], e3[3];
  T f1[3], f2[3], f3[3];
  T g1[3], g2[3], g3[3];
  T h1[3], h2[3], h3[3];
  T n1[3], m1[3];
  T z[3];

  T ef11[3], ef12[3], ef13[3];
  T ef21[3], ef22[3], ef23[3];
  T ef31[3], ef32[3], ef33[3];

  z[0] = 0.0;  z[1] = 0.0;  z[2] = 0.0;

  p1[0] = P1[0] - P1[0];  p1[1] = P1[1] - P1[1];  p1[2] = P1[2] - P1[2];
  p2[0] = P2[0] - P1[0];  p2[1] = P2[1] - P1[1];  p2[2] = P2[2] - P1[2];
  p3[0] = P3[0] - P1[0];  p3[1] = P3[1] - P1[1];  p3[2] = P3[2] - P1[2];

  q1[0] = Q1[0] - P1[0];  q1[1] = Q1[1] - P1[1];  q1[2] = Q1[2] - P1[2];
  q2[0] = Q2[0] - P1[0];  q2[1] = Q2[1] - P1[1];  q2[2] = Q2[2] - P1[2];
  q3[0] = Q3[0] - P1[0];  q3[1] = Q3[1] - P1[1];  q3[2] = Q3[2] - P1[2];

  e1[0] = p2[0] - p1[0];  e1[1] = p2[1] - p1[1];  e1[2] = p2[2] - p1[2];
  e2[0] = p3[0] - p2[0];  e2[1] = p3[1] - p2[1];  e2[2] = p3[2] - p2[2];
  e3[0] = p1[0] - p3[0];  e3[1] = p1[1] - p3[1];  e3[2] = p1[2] - p3[2];

  f1[0] = q2[0] - q1[0];  f1[1] = q2[1] - q1[1];  f1[2] = q2[2] - q1[2];
  f2[0] = q3[0] - q2[0];  f2[1] = q3[1] - q2[1];  f2[2] = q3[2] - q2[2];
  f3[0] = q1[0] - q3[0];  f3[1] = q1[1] - q3[1];  f3[2] = q1[2] - q3[2];

  VcrossV(n1, e1, e2);
  VcrossV(m1, f1, f2);

  VcrossV(g1, e1, n1);
  VcrossV(g2, e2, n1);
  VcrossV(g3, e3, n1);
  VcrossV(h1, f1, m1);
  VcrossV(h2, f2, m1);
  VcrossV(h3, f3, m1);

  VcrossV(ef11, e1, f1);
  VcrossV(ef12, e1, f2);
  VcrossV(ef13, e1, f3);
  VcrossV(ef21, e2, f1);
  VcrossV(ef22, e2, f2);
  VcrossV(ef23, e2, f3);
  VcrossV(ef31, e3, f1);
  VcrossV(ef32, e3, f2);
  VcrossV(ef33, e3, f3);

  // now begin the series of tests

  if (!my_project6_2(n1, p1, p2, p3, q1, q2, q3)) return 0;
  if (!my_project6_2(m1, p1, p2, p3, q1, q2, q3)) return 0;

  if (!my_project6_2(ef11, p1, p2, p3, q1, q2, q3)) return 0;
  if (!my_project6_2(ef12, p1, p2, p3, q1, q2, q3)) return 0;
  if (!my_project6_2(ef13, p1, p2, p3, q1, q2, q3)) return 0;
  if (!my_project6_2(ef21, p1, p2, p3, q1, q2, q3)) return 0;
  if (!my_project6_2(ef22, p1, p2, p3, q1, q2, q3)) return 0;
  if (!my_project6_2(ef23, p1, p2, p3, q1, q2, q3)) return 0;
  if (!my_project6_2(ef31, p1, p2, p3, q1, q2, q3)) return 0;
  if (!my_project6_2(ef32, p1, p2, p3, q1, q2, q3)) return 0;
  if (!my_project6_2(ef33, p1, p2, p3, q1, q2, q3)) return 0;

  if (!my_project6_2(g1, p1, p2, p3, q1, q2, q3)) return 0;
  if (!my_project6_2(g2, p1, p2, p3, q1, q2, q3)) return 0;
  if (!my_project6_2(g3, p1, p2, p3, q1, q2, q3)) return 0;
  if (!my_project6_2(h1, p1, p2, p3, q1, q2, q3)) return 0;
  if (!my_project6_2(h2, p1, p2, p3, q1, q2, q3)) return 0;
  if (!my_project6_2(h3, p1, p2, p3, q1, q2, q3)) return 0;

  return 1;
}

//
//check if two segments ab and cd overlap
//
template <typename Pt>
inline bool is_overlapping(
    const Pt& a,const Pt& b,const Pt& c,const Pt& d)
{
    auto area = Area(a,b,c);
    if(fabs(area)>SMALLNUMBER) return false;
    area=Area(a,b,d);
    if(fabs(area)>SMALLNUMBER) return false;

    const Pt u[2];
    return Union(a,b,c,d,u);
}

#endif//_MKSUM_INTERSECTION_H_


