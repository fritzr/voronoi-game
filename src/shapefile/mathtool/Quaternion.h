/********************************************************************

    Quaternion.h    Header File

    Jyh-Ming Lien 03/30/2002
    Computer Science.
    Texas A&M University

*********************************************************************/

#if !defined( _QUATERNION_H_NPRM_ )
#define _QUATERNION_H_NPRM_

#include "Matrix.h"
#include "Vector.h"

namespace mathtool{

    class Quaternion {
    public:
        Quaternion(){ m_s=1; m_v.set(0,0,0); }
        Quaternion( REAL s, Vector3d v ){ m_v=v; m_s=s; }
        Quaternion( const Quaternion & q){ *this=q; }
        Quaternion( const REAL r[3]){  //compute rotational quaternion
            //convert rot to quaternion
            REAL sin_r0=sin(r[0]*0.5);
            REAL cos_r0=cos(r[0]*0.5);
            REAL sin_r1=sin(r[1]*0.5);
            REAL cos_r1=cos(r[1]*0.5);
            REAL sin_r2=sin(r[2]*0.5);
            REAL cos_r2=cos(r[2]*0.5);

/*
            Quaternion qx(cos_r0,sin_r0*Vector3d(1,0,0));
            Quaternion qy(cos_r1,sin_r1*Vector3d(0,1,0));
            Quaternion qz(cos_r2,sin_r2*Vector3d(0,0,1));
            *this=(qz*qy*qx).normalize();
*/

            //Quaternion Q=(qz*qy*qx).normalize();
            m_s=(REAL)(cos_r0*cos_r1*cos_r2+sin_r0*sin_r1*sin_r2);
            m_v[0]=(REAL)(sin_r0*cos_r1*cos_r2-cos_r0*sin_r1*sin_r2);
            m_v[1]=(REAL)(cos_r0*sin_r1*cos_r2+sin_r0*cos_r1*sin_r2);
            m_v[2]=(REAL)(cos_r0*cos_r1*sin_r2-sin_r0*sin_r1*cos_r2);

            //cout<<"before="<<r[0]<<","<<r[1]<<","<<r[2]<<" after="<<getEuler()<<endl; //<<" Q="<<Q<<" this="<<*this<<endl;
        }

        ////////////////////////////////////////////////////////////////////////
        // Operations for Quaternion
        Quaternion operator*(const Quaternion & q) const {
            REAL s=q.m_s*m_s-q.m_v*m_v;
            Vector3d v=q.m_v*m_s+m_v*q.m_s+m_v%q.m_v;
            return Quaternion(s,v);
        }
        Quaternion operator*(REAL t) const {return Quaternion(m_s*t,m_v*t);}
        Quaternion operator*(const Vector3d & v) const { return *this*Quaternion(0,v); }
        Quaternion operator/(REAL s) const { return Quaternion(m_s/s,m_v/s); }
        Quaternion & operator=(const Quaternion & q){ set(q.m_s,q.m_v); return *this; }
        Quaternion operator+(const Quaternion & q) const { return Quaternion(m_s+q.m_s,m_v+q.m_v); }
        Quaternion operator-(const Quaternion & q) const { return Quaternion(m_s-q.m_s,m_v-q.m_v); }
        Quaternion operator-() const { return Quaternion(m_s,-m_v); }
        friend Quaternion operator*(const Vector3d & v, const Quaternion & q);
        friend istream& operator>>(istream & in, Quaternion & q );
        friend ostream& operator<<(ostream & out, const Quaternion & q );

        //////////////////////////////////////////////////////////////////////////
        //Normalization
        Quaternion normalize(){ 
            Quaternion q(*this);
            REAL l=q.norm();
            q=q/l;
            return q;
        }

        REAL norm(){ return std::sqrt(normsqr()); }
        REAL normsqr(){ return m_v.normsqr()+sqr(m_s); }
        
        //////////////////////////////////////////////////////////////////////////
        //Access

        Matrix3x3 getMatrix() const {
            REAL x_2=2*sqr(m_v[0]); REAL y_2=2*sqr(m_v[1]); REAL z_2=2*sqr(m_v[2]);
            REAL xy=2*m_v[0]*m_v[1]; REAL yz=2*m_v[1]*m_v[2]; REAL zx=2*m_v[2]*m_v[0]; 
            REAL sx=2*m_s*m_v[0]; REAL sy=2*m_s*m_v[1]; REAL sz=2*m_s*m_v[2]; 
            return Matrix3x3(1-y_2-z_2, xy-sz, zx+sy,
                             xy+sz, 1-x_2-z_2, yz-sx,
                             zx-sy, yz+sx, 1-x_2-y_2);
        }

        Vector3d getEuler() const { //x-y-z conversion
            Vector3d angle;
            angle[0]=atan2(2*(m_s*m_v[0]+m_v[1]*m_v[2]),1-2*(m_v[0]*m_v[0]+m_v[1]*m_v[1]));
            angle[1]=asin(2*(m_s*m_v[1]-m_v[0]*m_v[2]));
            angle[2]=atan2(2*(m_s*m_v[2]+m_v[0]*m_v[1]),1-2*(m_v[1]*m_v[1]+m_v[2]*m_v[2]));
            return angle;
        }
    

        void set(REAL s,const Vector3d & v){ m_v=v; m_s=s; }
        void set(REAL q1, REAL q2, REAL q3, REAL q4){ m_s=q1; m_v.set(q2,q3,q4); }
        const Vector3d& getComplex() const { return m_v; }
        REAL getReal() const { return m_s; }

    private:
        Vector3d m_v;
        REAL m_s;
    };


    //////////////////////////////////////////////////////////////////////////
    //Interpolation
    Quaternion slerp(const Quaternion& q0, const Quaternion& q1, REAL t);


}//emd of nprmlib

#endif //#if !defined( _QUATERNION_H_NPRM_ ) 

