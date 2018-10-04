//------------------------------------------------------------------------------
//  Copyright 2010-2017 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#pragma once
#ifndef _INTERSECTION_H_
#define _INTERSECTION_H_

#include "Basic.h"

// Most of these functions are copied/modified from CG in C

template<class real>
void Assign( real p[2], real a[2])
{
   p[0] = a[0];
   p[1] = a[1];
}

template<class real> real Area(const real a[2], const real b[2], const real c[2])
{
	return ( b[0] - a[0] ) * ( c[1] - a[1] ) -
		   ( c[0] - a[0] ) * ( b[1] - a[1] );
}

template<class real> int AreaSign(const real a[2], const real b[2], const real c[2])
{
	real area=Area(a,b,c);

	if      ( area >  SMALLNUMBER ) return  1;
	else if ( area <  -SMALLNUMBER ) return -1;
	else     return  0;
}

template<class real> int Collinear(const real a[2], const real b[2], const real c[2])
{
   return AreaSign( a, b, c ) == 0;
}

/*---------------------------------------------------------------------
Returns TRUE iff point c lies on the closed segement ab.
Assumes it is already known that abc are collinear.
---------------------------------------------------------------------*/
template<class real> bool Between(const real a[2], const real b[2], const real c[2])
{
   // If ab not vertical, check betweenness on x; else on y.
   if ( a[0]!=b[0] )
	  return ((a[0] <= c[0]) && (c[0] <= b[0])) ||
			 ((a[0] >= c[0]) && (c[0] >= b[0]));
   else
	  return ((a[1] <= c[1]) && (c[1] <= b[1])) ||
			 ((a[1] >= c[1]) && (c[1] >= b[1]));
}

template<class real> bool Between_strict(const real a[2], const real b[2], const real c[2])
{
   // If ab not vertical, check betweenness on x; else on y.
	if ( a[0]!=b[0] )
	{
	  return ((a[0] < c[0]) && (c[0] < b[0])) ||
			 ((a[0] > c[0]) && (c[0] > b[0]));
	}
	else{
		return ((a[1] < c[1]) && (c[1] < b[1])) ||
			   ((a[1] > c[1]) && (c[1] > b[1]));
	}
}


template<class real> bool AlmostEqual3(const real a[3], const real b[3])
{
	return (fabs(a[0]-b[0])<SMALLNUMBER &&fabs(a[1]-b[1])<SMALLNUMBER && fabs(a[2]-b[2])<SMALLNUMBER);
}

template<class real> bool AlmostEqual(const real a[2], const real b[2])
{
	return (fabs(a[0]-b[0])<SMALLNUMBER &&fabs(a[1]-b[1])<SMALLNUMBER);
}

template<class real> bool Equal3(const real a[3], const real b[3])
{
	return (a[0]==b[0])&& (a[1]==b[1]) && (a[2]==b[2]);
}

template<class real> bool Equal(const real a[2], const real b[2])
{
	return (a[0]==b[0])&& (a[1]==b[1]);
}

template<class real> char ParallelInt
(const real a[2], const real b[2], const real c[2], const real d[2], real p[2])
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


// compute union of two colinear  segments ab and cd
// place the union in p
// return false if the union is degenerated
template<class real> bool Union
(const real a[2], const real b[2], const real c[2], const real d[2],
 const real * p[2])
{
	int id=0;
	if(Equal(a,c)){ p[id]=a; id++; }

	if(Equal(a,d)){ p[id]=a; id++; }
	if(id==2) return p[0]!=p[1];

	if(Equal(b,c)){ p[id]=b; id++; }
	if(id==2) return p[0]!=p[1];

	if(Equal(b,d)){ p[id]=b; id++; }
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

/*---------------------------------------------------------------------
SegSegInt: Detect intersection between two closed segments ab and cd.
---------------------------------------------------------------------*/


template<class real> bool SegSegInt( const real a[2], const real b[2],
					   const real c[2], const real d[2])
{
	//check X and Y coordinates
	for(int i=0;i<2;i++){
		if(a[i]<b[i]){
			if(c[i]<d[i]){
				if( a[i]>d[i] ) return false;
				if( c[i]>b[i] ) return false;
			}
			else{ //c[i]>=d[i]
				if( a[i]>c[i] ) return false;
				if( d[i]>b[i] ) return false;
			}
		}
		else{ //a[i]>=b[i]
			if(c[i]<d[i]){
				if( b[i]>d[i] ) return false;
				if( c[i]>a[i] ) return false;
			}
			else{ //c[i]>=d[i]
				if( b[i]>c[i] ) return false;
				if( d[i]>a[i] ) return false;
			}
		}
	}//end for i

	//OK potential intersection
	int abc=AreaSign(a,b,c);
	int abd=AreaSign(a,b,d);

	if(abc==0 && abd==0 ) { //collinear
		//check if they overlap..
		if(Between(a,b,c)) return true;
		if(Between(a,b,d)) return true;
		if(Between(c,d,a)) return true;
		if(Between(c,d,b)) return true;
		return false;
	}
	else if(abc==0){
		if(Between(a,b,c)) return true;
		return false;
	}
	else if(abd==0){
		if(Between(a,b,d)) return true;
		return false;
	}
	//
	else{ // if(abc!=0 && abd!=0)
		if(abc>0 && abd>0) return false;
		if(abc<0 && abd<0) return false;
	}

	int cda=AreaSign(c,d,a);
	int cdb=AreaSign(c,d,b);

	assert(cda!=0 || cdb!=0);

	if(cda==0){
		if(Between(c,d,a)) return true;
		return false;
	}
	else if(cdb==0){
		if(Between(c,d,b)) return true;
		return false;
	}
	else{
		if(cda>0 && cdb>0) return false;
		if(cda<0 && cdb<0) return false;
	}

	return true;
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
template<class real>
char SegSegInt( const real a[2], const real b[2],
				const real c[2], const real d[2],
				real p[2] )
{
	real  s, t;                  // The two parameters of the parametric eqns.
	real  num_s, num_t, denom;   // Numerator and denoninator of equations.
	char  code = '?';            // Return char characterizing intersection.

	const real small_number=SMALLNUMBER;

	//
	if(a[0]==c[0] && a[1]==c[1]){ p[0]=a[0]; p[1]=a[1]; return 'v'; }
	if(b[0]==c[0] && b[1]==c[1]){ p[0]=b[0]; p[1]=b[1]; return 'v'; }
	if(a[0]==d[0] && a[1]==d[1]){ p[0]=a[0]; p[1]=a[1]; return 'v'; }
	if(b[0]==d[0] && b[1]==d[1]){ p[0]=b[0]; p[1]=b[1]; return 'v'; }
	//

	denom = a[0] * ( d[1] - c[1] ) +
			b[0] * ( c[1] - d[1] ) +
			d[0] * ( b[1] - a[1] ) +
			c[0] * ( a[1] - b[1] );

	// If denom is zero, then segments are parallel: handle separately.
	if (denom==0) return  ParallelInt(a, b, c, d, p);

	real denom_small=denom*small_number;

	//compute s
	num_s =    a[0] * ( d[1] - c[1] ) +
			   c[0] * ( a[1] - d[1] ) +
			   d[0] * ( c[1] - a[1] );

	//if( fabs(num_s)<SMALLNUMBER ) num_s=0;
	//else if( fabs(num_s-denom)<SMALLNUMBER ) num_s=denom;

	s = num_s / denom;

	if( fabs(s)<small_number ) s=0;
	else if( fabs(1-s)<small_number ) s=1;

	//if(num_s<0  || num_s>denom) return '0';
	if(s<0  || s>1) return '0';

	//compute t
	num_t = -( a[0] * ( c[1] - b[1] ) +
			   b[0] * ( a[1] - c[1] ) +
			   c[0] * ( b[1] - a[1] ) );


	t = num_t / denom;

	if( fabs(t)<small_number ) t=0;
	else if( fabs(1-t)<small_number ) t=1;

	if(t<0  || t>1) return '0';

	//decide the code
	if( (0.0<s) && (s<1) && (0.0<t) && (t<1) ) code = '1';
	else code= 'v';

	if(code!='v'){
		//s = num_s / denom;
		p[0] = (a[0] + s*(b[0]-a[0]));
		p[1] = (a[1] + s*(b[1]-a[1]));
	}
	else{
		if(s==0){ p[0]=a[0]; p[1]=a[1]; }
		else if(s==1){ p[0]=b[0]; p[1]=b[1]; }
		else if(t==0){ p[0]=c[0]; p[1]=c[1]; }
		else if(t==1){ p[0]=d[0]; p[1]=d[1]; }
		else
			assert(false);
	}

	return code;
}

#endif//_INTERSECTION_H_
