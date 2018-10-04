#ifndef _MKSUM_INTERSECTION_H_
#define _MKSUM_INTERSECTION_H_


#include <Basic.h>
#include <Point.h>
using namespace mathtool;

// This is copied from CG in C

inline void Assign( double p[2], double a[2])
{
   p[0] = a[0];
   p[1] = a[1];
}

inline double Area(const double a[2], const double b[2], const double c[2])
{
    return ( b[0] - a[0] ) * ( c[1] - a[1] ) -
           ( c[0] - a[0] ) * ( b[1] - a[1] );
}

inline int AreaSign(const double a[2], const double b[2], const double c[2])
{
    double area=Area(a,b,c);
    if      ( area >  SMALLNUMBER ) return  1;
    else if ( area < -SMALLNUMBER ) return -1;
    else     return  0;
}

inline int Collinear(const double a[2], const double b[2], const double c[2])
{
   return AreaSign( a, b, c ) == 0;
}

/*---------------------------------------------------------------------
Returns TRUE iff point c lies on the closed segement ab.
Assumes it is already known that abc are collinear.
---------------------------------------------------------------------*/
inline bool Between(const double a[2], const double b[2], const double c[2])
{
   // If ab not vertical, check betweenness on x; else on y.
   if ( fabs(a[0]-b[0])>fabs(a[1]-b[1]) )
      return ((a[0] <= c[0]) && (c[0] <= b[0])) ||
             ((a[0] >= c[0]) && (c[0] >= b[0]));
   else
      return ((a[1] <= c[1]) && (c[1] <= b[1])) ||
             ((a[1] >= c[1]) && (c[1] >= b[1]));
}

inline bool Between_strict(const double a[2], const double b[2], const double c[2])
{
   // If ab not vertical, check betweenness on x; else on y.
    if ( fabs(a[0]-b[0])>SMALLNUMBER ){
      double c01=c[0]-SMALLNUMBER;
      double c02=c[0]+SMALLNUMBER;

      return ((a[0] < c01) && (c02 < b[0])) ||
             ((a[0] > c02) && (c01 > b[0]));
    }
    else{
        double c11=c[1]-SMALLNUMBER;
        double c10=c[1]+SMALLNUMBER;
        return ((a[1] < c11) && (c10 < b[1])) ||
               ((a[1] > c10) && (c11 > b[1]));
    }
}


inline bool AlmostEqual3(const double a[3], const double b[3])
{
    return (fabs(a[0]-b[0])<SMALLNUMBER &&fabs(a[1]-b[1])<SMALLNUMBER && fabs(a[2]-b[2])<SMALLNUMBER);
}

inline bool AlmostEqual(const double a[2], const double b[2])
{
    return (fabs(a[0]-b[0])<SMALLNUMBER &&fabs(a[1]-b[1])<SMALLNUMBER);
}

inline bool AlmostEqual(const double a[2], const double b[2], double tau)
{
    return (fabs(a[0]-b[0])<tau &&fabs(a[1]-b[1])<tau);
}


// compute union of two colinear  segments ab and cd
// place the union in p
// return false if the union is degenerated
inline bool Union
(const double a[2], const double b[2], const double c[2], const double d[2],
 const double * p[2])
{
    int id=0;
    if(AlmostEqual(a,c)){ p[id]=a; id++; }
    
    if(AlmostEqual(a,d)){ p[id]=a; id++; }
    if(id==2) return p[0]!=p[1];
    
    if(AlmostEqual(b,c)){ p[id]=b; id++; }
    if(id==2) return p[0]!=p[1];

    if(AlmostEqual(b,d)){ p[id]=b; id++; }
    if(id==2) return p[0]!=p[1];

    if( Between_strict(a,b,c) ){ p[id]=c; id++; } 
    if(id==2) return p[0]!=p[1];
    
    if( Between_strict(a,b,d) ){ p[id]=d; id++; } 
    if(id==2) return p[0]!=p[1];

    if( Between_strict(c,d,a) ){ p[id]=a; id++; }
    if(id==2) return p[0]!=p[1];

    if( Between_strict(c,d,b) ){ p[id]=b; id++; }
    if(id==2) return p[0]!=p[1];

    return false;
}


inline char ParallelInt
( const double a[2], const double b[2], const double c[2], const double d[2], double p[2])
{
   if(!Collinear(a, b, c)) return '0';

   //check if they overlap..
   if(Between(a,b,c)) return 'e';
   if(Between(a,b,d)) return 'e';
   if(Between(c,d,a)) return 'e';
   if(Between(c,d,b)) return 'e';

   //they don't overlap but the end points may..
   //check if the end points overlap
   if(AlmostEqual(a,c)){ p[0]=a[0]; p[1]=a[1]; return 'v';}
   if(AlmostEqual(b,c)){ p[0]=b[0]; p[1]=b[1]; return 'v';}
   if(AlmostEqual(a,d)){ p[0]=a[0]; p[1]=a[1]; return 'v';}
   if(AlmostEqual(b,d)){ p[0]=b[0]; p[1]=b[1]; return 'v';}

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
inline char SegSegInt( const double a[2], const double b[2], const double c[2], const double d[2], double p[2] )
{
   double  s, t;        // The two parameters of the parametric eqns.
   double  num_s, num_t, denom;   // Numerator and denoninator of equations.
   char    code = '?';    // Return char characterizing intersection.
   const double TINYNUMBER=SMALLNUMBER;

   denom = a[0] * ( d[1] - c[1] ) +
           b[0] * ( c[1] - d[1] ) +
           d[0] * ( b[1] - a[1] ) +
           c[0] * ( a[1] - b[1] );

   // If denom is zero, then segments are parallel: handle separately.
   //if (fabs(denom)<SMALLNUMBER) denom=0;
   if (fabs(denom)<TINYNUMBER) return  ParallelInt(a, b, c, d, p);


   if(AlmostEqual(a,c)){ p[0]=a[0]; p[1]=a[1]; return 'v';}
   if(AlmostEqual(b,c)){ p[0]=b[0]; p[1]=b[1]; return 'v';}
   if(AlmostEqual(a,d)){ p[0]=a[0]; p[1]=a[1]; return 'v';}
   if(AlmostEqual(b,d)){ p[0]=b[0]; p[1]=b[1]; return 'v';}

   //compute s
   num_s =    a[0] * ( d[1] - c[1] ) +
              c[0] * ( a[1] - d[1] ) +
              d[0] * ( c[1] - a[1] );


   s = num_s / denom;

   if( fabs(s)<TINYNUMBER ) s=0;
   else if( fabs(1-s)<TINYNUMBER ) s=1;


   //compute t
   num_t = -( a[0] * ( c[1] - b[1] ) +
              b[0] * ( a[1] - c[1] ) +
              c[0] * ( b[1] - a[1] ) );



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
       p[0] = (a[0] + s*(b[0]-a[0]));
       p[1] = (a[1] + s*(b[1]-a[1]));
   }
   else{
      if(s==0){ p[0]=a[0]; p[1]=a[1]; }
      else if(s==1){p[0]=b[0]; p[1]=b[1]; }
      else if(t==0){p[0]=c[0]; p[1]=c[1]; }
      else if(t==1){p[0]=d[0]; p[1]=d[1]; }
      else{
          cout<<"s="<<s<<" t="<<t<<endl;
          cout<<"a="<<a[0]<<","<<a[1]<<endl;
          cout<<"b="<<b[0]<<","<<b[1]<<endl;
          cout<<"c="<<c[0]<<","<<c[1]<<endl;
          cout<<"d="<<d[0]<<","<<d[1]<<endl;
          assert(false);
      }
   }

   //double check...
//   if(code=='v'){
//       double q[2];
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


inline double
VdotV(double V1[3], double V2[3])
{
  return (V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2]);
}



inline void
VcrossV(double Vr[3], const double V1[3], const double V2[3])
{
  Vr[0] = V1[1]*V2[2] - V1[2]*V2[1];
  Vr[1] = V1[2]*V2[0] - V1[0]*V2[2];
  Vr[2] = V1[0]*V2[1] - V1[1]*V2[0];
}

inline double
max(double a, double b, double c)
{
  double t = a;
  if (b > t) t = b;
  if (c > t) t = c;
  return t;
}
 
inline double
min(double a, double b, double c)
{
  double t = a;
  if (b < t) t = b;
  if (c < t) t = c;
  return t;
}


inline int
my_project6_2(double *ax,
     double *p1, double *p2, double *p3,
     double *q1, double *q2, double *q3)
{
  double P1 = VdotV(ax, p1);
  double P2 = VdotV(ax, p2);
  double P3 = VdotV(ax, p3);
  double Q1 = VdotV(ax, q1);
  double Q2 = VdotV(ax, q2);
  double Q3 = VdotV(ax, q3);
  
  double mx1 = max(P1, P2, P3);
  double mn1 = min(P1, P2, P3);
  double mx2 = max(Q1, Q2, Q3);
  double mn2 = min(Q1, Q2, Q3);

  if (mn1 > mx2) return 0;
  if (mn2 > mx1) return 0;
  return 1;
}


// very robust triangle intersection test
// uses no divisions
// works on coplanar triangles

inline int
my_tri_contact (const double *P1, const double *P2, const double *P3,
        const double *Q1, const double *Q2, const double *Q3)
{

  /*
     One triangle is (p1,p2,p3).  Other is (q1,q2,q3).
     Edges are (e1,e2,e3) and (f1,f2,f3).
     Normals are n1 and m1
     Outwards are (g1,g2,g3) and (h1,h2,h3).

     We assume that the triangle vertices are in the same coordinate system.

     First thing we do is establish a new c.s. so that p1 is at (0,0,0).
     */

  double p1[3], p2[3], p3[3];
  double q1[3], q2[3], q3[3];
  double e1[3], e2[3], e3[3];
  double f1[3], f2[3], f3[3];
  double g1[3], g2[3], g3[3];
  double h1[3], h2[3], h3[3];
  double n1[3], m1[3];
  double z[3];

  double ef11[3], ef12[3], ef13[3];
  double ef21[3], ef22[3], ef23[3];
  double ef31[3], ef32[3], ef33[3];
  
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
inline bool is_overlapping(const Point2d& a,const Point2d& b,const Point2d& c,const Point2d& d)
{
    double area=Area(a.get(),b.get(),c.get());
    if(fabs(area)>SMALLNUMBER) return false;
    area=Area(a.get(),b.get(),d.get());
    if(fabs(area)>SMALLNUMBER) return false;

    const double * u[2]={NULL,NULL};
    const double * s=a.get(); const double * t=b.get();
    const double * m=c.get(); const double * n=d.get();

    return Union(s,t,m,n,u);
}

#endif//_MKSUM_INTERSECTION_H_


