/********************************************************************

    Basic.h

    Jyh-Ming Lien 03/30/2002
    Computer Science.
    Texas A&M University

*********************************************************************/

#ifndef _H_UTILITY
#define _H_UTILITY

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <string>
#include <iostream>   // define C++ stream I/O routines
#include <iomanip>
using namespace std;

#include "mpreal.h"
using namespace mpfr;

namespace mathtool{

	//data type
	#define REAL double //double //mpreal

    /* range of real numbers */
    #define SMALLNUMBER 1.0e-10
    #define HUGENUMBER  1.0e10

    /* Miscellaneous Scalar Math */
    #define abs(x)      (((x) < 0) ? (-(x)) : (x))
    #define sqr(x)      ((x) * (x))

    inline int round( mpreal x, mpreal p){
        return (long int)(((long int)((x)*pow(10.0,p)+((x)<0?-0.5:0.5)).toLong())/pow(10.0,p)).toLong();
    }

    inline int round( mpreal v ){
        int integer=(int)floor((long int)(v.toLong()));
        mpreal fraction=v-integer;

        if(v>0)
            return (fraction>=0.5)?integer+1:integer;
        else
            return (fraction>=-0.5)?integer:integer+1;
    }

    #define sign(x)     ((x)>=0? 1: -1)
    //#define swap(x1, x2)  {int tmp=x1; x1=x2; x2=tmp}
    #define applysign(x, y) ((y) >= 0? abs(x): -abs(x))

    /* Angle Conversions & Constants */

    #ifndef PI
    #define PI ((REAL)3.1415926535897)
    #endif

    #ifndef PI2
    #define PI2 ((REAL)6.2831853071794)
    #endif

    #define RAD2DEG (180/PI)
    #define DEG2RAD (PI/180)

    #define DegToRad(x) ((x)*DEG2RAD)
    #define RadToDeg(x) ((x)*RAD2DEG)

    /*
      computes sqrt(a^2 + b^2) without destructive underflow or overflow
    */
    mpreal pythag(mpreal a, mpreal b);

    /*
      Utility Error message routines
    */
    // print s to stdout with trailing blank and no terminating end of line
    void prompt(char *s);

    // print s1, s2, s3 to stdout as blank separated fields with terminating eol
    void message(char *s1, char *s2 = NULL, char *s3 = NULL);

    // print Status: to stdout followed by message(s1, s2, s3)
    void status(char *s1, char *s2 = NULL, char *s3 = NULL);

    // print Error: followed by s1, s2 and s3 to stderr as blank separated fields 
    // with terminating eol
    void error(char *s1, char *s2 = NULL, char *s3 = NULL);

    // print error(s1, s2, s3) and then exit program with code 1 
    void abort(char *s1, char *s2 = NULL, char *s3 = NULL);

    ///Added by Jyh-Ming Lien
    /*
    bool getDoubleValue(char * pTag, mpreal * pValue,int size);
    bool getDoubleValue(char * pTag, mpreal * pValue);
    */


    #ifdef _WIN32

    ////////////////////////////////////////////////////////////////////////////////////////
    // Following functions define M_PI and drand48, which are not starndard c library and 
    // definitions. In addition, rint used to round off float points to int is also here.
    /////////////////////////////////////////////////////////////////////////////////////////

    #define M_PI 3.1415926 //reference PI above

    extern "C" {
        //Implementation of these functions are located in util.cpp
        mpreal drand48();
        mpreal erand48(register unsigned short *xsubi);
        long irand48(register unsigned short m);
        long krand48(register unsigned short *xsubi, unsigned short m);
        long lrand48();
        long mrand48();
        static void next();
        void srand48(long seedval);
        unsigned short * seed48(unsigned short seed16v[3]);
        void lcong48(unsigned short param[7]);
        long nrand48(register unsigned short *xsubi);
        long jrand48(register unsigned short *xsubi);

        /**Round to closest integer.
          *The rint() function rounds x to an integer value according
          *to the prevalent rounding mode.  The default rounding mode
          *is to round to the nearest integer.
          *@return The  rint() function returns the integer value as a float-
          *ing-point number.
          */
        mpreal rint(mpreal x);

    } //end extern "C"

    #endif //_WIN32

} //end of nprmlib

#endif
