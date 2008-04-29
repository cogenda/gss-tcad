#ifndef _vec6_h_
#define _vec6_h_
#include <math.h>
#include "petsc.h"

//this is the basic operate of vec and mat. speed is the most important thing
typedef struct
{
        PetscScalar v[6];
}Vec6;

typedef struct
{
        PetscScalar m[36];
}Mat6;

inline void Set_Vec6_zero(Vec6 &vec)
{
        vec.v[0]=0;
        vec.v[1]=0;
        vec.v[2]=0;
        vec.v[3]=0;
        vec.v[4]=0;
        vec.v[5]=0;
}

inline void Set_Vec6(Vec6 &vec,PetscScalar v1,PetscScalar v2,PetscScalar v3,PetscScalar v4,PetscScalar v5,PetscScalar v6)
{
        vec.v[0]=v1;
        vec.v[1]=v2;
        vec.v[2]=v3;
        vec.v[3]=v4;
        vec.v[4]=v5;
        vec.v[5]=v6;
}

inline Vec6 operator +(const Vec6 &vec1,const Vec6 &vec2)
{
        Vec6 sum;
        sum.v[0]=vec1.v[0]+vec2.v[0];
        sum.v[1]=vec1.v[1]+vec2.v[1];
        sum.v[2]=vec1.v[2]+vec2.v[2];
        sum.v[3]=vec1.v[3]+vec2.v[3];
        sum.v[4]=vec1.v[4]+vec2.v[4];
        sum.v[5]=vec1.v[5]+vec2.v[5];
        return sum;
}

inline Vec6 operator +(const Vec6 &vec1,const PetscScalar value)
{
        Vec6 sum;
        sum.v[0]=vec1.v[0]+value;
        sum.v[1]=vec1.v[1]+value;
        sum.v[2]=vec1.v[2]+value;
        sum.v[3]=vec1.v[3]+value;
        sum.v[4]=vec1.v[4]+value;
        sum.v[5]=vec1.v[5]+value;
        return sum;
}

inline Vec6 operator +(const PetscScalar value,const Vec6 &vec1)
{
        Vec6 sum;
        sum.v[0]=vec1.v[0]+value;
        sum.v[1]=vec1.v[1]+value;
        sum.v[2]=vec1.v[2]+value;
        sum.v[3]=vec1.v[3]+value;
        sum.v[4]=vec1.v[4]+value;
        sum.v[5]=vec1.v[5]+value;
        return sum;
}

inline Vec6 operator -(const Vec6 &vec1,const Vec6 &vec2)
{
        Vec6 sub;
        sub.v[0]=vec1.v[0]-vec2.v[0];
        sub.v[1]=vec1.v[1]-vec2.v[1];
        sub.v[2]=vec1.v[2]-vec2.v[2];
        sub.v[3]=vec1.v[3]-vec2.v[3];
        sub.v[4]=vec1.v[4]-vec2.v[4];
        sub.v[5]=vec1.v[5]-vec2.v[5];
        return sub;
}

inline Vec6 operator -(const Vec6 &vec1,const PetscScalar value)
{
        Vec6 sum;
        sum.v[0]=vec1.v[0]-value;
        sum.v[1]=vec1.v[1]-value;
        sum.v[2]=vec1.v[2]-value;
        sum.v[3]=vec1.v[3]-value;
        sum.v[4]=vec1.v[4]-value;
        sum.v[5]=vec1.v[5]-value;
        return sum;
}

inline Vec6 operator -(const PetscScalar value,const Vec6 &vec1)
{
        Vec6 sum;
        sum.v[0]=value-vec1.v[0];
        sum.v[1]=value-vec1.v[1];
        sum.v[2]=value-vec1.v[2];
        sum.v[3]=value-vec1.v[3];
        sum.v[4]=value-vec1.v[4];
        sum.v[5]=value-vec1.v[5];
        return sum;
}

inline Vec6 operator *(const Vec6 &vec,const PetscScalar &value)
{
        Vec6 result;
        result.v[0]=vec.v[0]*value;
        result.v[1]=vec.v[1]*value;
        result.v[2]=vec.v[2]*value;
        result.v[3]=vec.v[3]*value;
        result.v[4]=vec.v[4]*value;
        result.v[5]=vec.v[5]*value;
        return result;
}

inline Vec6 operator *(const PetscScalar &value,const Vec6 &vec)
{
        Vec6 result;
        result.v[0]=vec.v[0]*value;
        result.v[1]=vec.v[1]*value;
        result.v[2]=vec.v[2]*value;
        result.v[3]=vec.v[3]*value;
        result.v[4]=vec.v[4]*value;
        result.v[5]=vec.v[5]*value;
        return result;
}

inline Vec6 operator /(const Vec6 &vec,const PetscScalar &value)
{
        Vec6 result;
        result.v[0]=vec.v[0]/value;
        result.v[1]=vec.v[1]/value;
        result.v[2]=vec.v[2]/value;
        result.v[3]=vec.v[3]/value;
        result.v[4]=vec.v[4]/value;
        result.v[5]=vec.v[5]/value;
        return result;
}

inline Vec6 operator /(const PetscScalar &value,const Vec6 &vec)
{
        Vec6 result;
        result.v[0]=value/vec.v[0];
        result.v[1]=value/vec.v[1];
        result.v[2]=value/vec.v[2];
        result.v[3]=value/vec.v[3];
        result.v[4]=value/vec.v[4];
        result.v[5]=value/vec.v[5];
        return result;
}

inline Vec6 operator *(const Vec6 &vec1,const Vec6 &vec2)
{
        Vec6 result;
        result.v[0]=vec1.v[0]*vec2.v[0];
        result.v[1]=vec1.v[1]*vec2.v[1];
        result.v[2]=vec1.v[2]*vec2.v[2];
        result.v[3]=vec1.v[3]*vec2.v[3];
        result.v[4]=vec1.v[4]*vec2.v[4];
        result.v[5]=vec1.v[5]*vec2.v[5];
        return result;
}

inline Vec6 operator /(const Vec6 &vec1,const Vec6 &vec2)
{
        Vec6 result;
        result.v[0]=vec1.v[0]/vec2.v[0];
        result.v[1]=vec1.v[1]/vec2.v[1];
        result.v[2]=vec1.v[2]/vec2.v[2];
        result.v[3]=vec1.v[3]/vec2.v[3];
        result.v[4]=vec1.v[4]/vec2.v[4];
        result.v[5]=vec1.v[5]/vec2.v[5];
        return result;
}

//--------------------------------------------------------------------

inline PetscScalar vmaxvalue(const Vec6 &a)
{
        PetscScalar v=a.v[0];
        if(a.v[1]>v)  v=a.v[1];
        if(a.v[2]>v)  v=a.v[2];
        if(a.v[3]>v)  v=a.v[3];
        if(a.v[4]>v)  v=a.v[4];
        if(a.v[5]>v)  v=a.v[5];
        return v;
}

inline Vec6 vmax(const Vec6 &a,const Vec6 &b)
{
        Vec6 result;
        result.v[0]=a.v[0]>b.v[0] ? a.v[0]:b.v[0];
        result.v[1]=a.v[1]>b.v[1] ? a.v[1]:b.v[1];
        result.v[2]=a.v[2]>b.v[2] ? a.v[2]:b.v[2];
        result.v[3]=a.v[3]>b.v[3] ? a.v[3]:b.v[3];
        result.v[4]=a.v[4]>b.v[4] ? a.v[4]:b.v[4];
        result.v[5]=a.v[5]>b.v[5] ? a.v[5]:b.v[5];
        return result;
}

inline Vec6 vmin(const Vec6 &a,const Vec6 &b)
{
        Vec6 result;
        result.v[0]=a.v[0]<b.v[0] ? a.v[0]:b.v[0];
        result.v[1]=a.v[1]<b.v[1] ? a.v[1]:b.v[1];
        result.v[2]=a.v[2]<b.v[2] ? a.v[2]:b.v[2];
        result.v[3]=a.v[3]<b.v[3] ? a.v[3]:b.v[3];
        result.v[4]=a.v[4]<b.v[4] ? a.v[4]:b.v[4];
        result.v[5]=a.v[5]<b.v[5] ? a.v[5]:b.v[5];
        return result;
}

inline Vec6 vabs(const Vec6 &a)
{
        Vec6 result;
        result.v[0]=fabs(a.v[0]);
        result.v[1]=fabs(a.v[1]);
        result.v[2]=fabs(a.v[2]);
        result.v[3]=fabs(a.v[3]);
        result.v[4]=fabs(a.v[4]);
        result.v[5]=fabs(a.v[5]);
        return result;
}

inline Vec6 vsign(const Vec6 &a)
{
        Vec6 result;
        result.v[0]=a.v[0]>0.0? 1.0:-1.0;
        result.v[1]=a.v[1]>0.0? 1.0:-1.0;
        result.v[2]=a.v[2]>0.0? 1.0:-1.0;
        result.v[3]=a.v[3]>0.0? 1.0:-1.0;
        result.v[4]=a.v[4]>0.0? 1.0:-1.0;
        result.v[5]=a.v[5]>0.0? 1.0:-1.0;
        return result;
}

//----------------------------------------------------
inline void  Set_Mat6_zero(Mat6 &A)
{
         A.m[0]=0;       A.m[1]=0;       A.m[2]=0;       A.m[3]=0;       A.m[4]=0;       A.m[5]=0;       
         A.m[6]=0;       A.m[7]=0;       A.m[8]=0;       A.m[9]=0;       A.m[10]=0;      A.m[11]=0;      
         A.m[12]=0;      A.m[13]=0;      A.m[14]=0;      A.m[15]=0;      A.m[16]=0;      A.m[17]=0;
         A.m[18]=0;      A.m[19]=0;      A.m[20]=0;      A.m[21]=0;      A.m[22]=0;      A.m[23]=0;
         A.m[24]=0;      A.m[25]=0;      A.m[26]=0;      A.m[27]=0;      A.m[28]=0;      A.m[29]=0;
         A.m[30]=0;      A.m[31]=0;      A.m[32]=0;      A.m[33]=0;      A.m[34]=0;      A.m[35]=0;  
}

inline void Set_Mat6_I(Mat6 &A)
{
         A.m[0]=1;       A.m[1]=0;       A.m[2]=0;       A.m[3]=0;       A.m[4]=0;       A.m[5]=0;       
         A.m[6]=0;       A.m[7]=1;       A.m[8]=0;       A.m[9]=0;       A.m[10]=0;      A.m[11]=0;      
         A.m[12]=0;      A.m[13]=0;      A.m[14]=1;      A.m[15]=0;      A.m[16]=0;      A.m[17]=0;
         A.m[18]=0;      A.m[19]=0;      A.m[20]=0;      A.m[21]=1;      A.m[22]=0;      A.m[23]=0;
         A.m[24]=0;      A.m[25]=0;      A.m[26]=0;      A.m[27]=0;      A.m[28]=1;      A.m[29]=0;
         A.m[30]=0;      A.m[31]=0;      A.m[32]=0;      A.m[33]=0;      A.m[34]=0;      A.m[35]=1;  
}

inline Mat6 operator +(const Mat6 &A, const Mat6 &B)
{
        Mat6 C;
        for(int i=0;i<36;i++)
                C.m[i] = A.m[i]+B.m[i];
        return C;
}

inline Mat6 operator -(const Mat6 &A, const Mat6 &B)
{
        Mat6 C;
        for(int i=0;i<36;i++)
                C.m[i] = A.m[i]-B.m[i];
        return C;
}

inline Mat6  operator -(const Mat6 &A)
{
        Mat6 C;
        for(int i=0;i<36;i++)
                C.m[i] = -A.m[i];
        return C;
}

inline Mat6  operator *(const Mat6 &A,const PetscScalar &a)
{
        Mat6 C;
        for(int i=0;i<36;i++)
                C.m[i] = a*A.m[i];
        return C;
}

inline Mat6  operator *(const PetscScalar &a,const Mat6 &A)
{
        Mat6 C;
        for(int i=0;i<36;i++)
                C.m[i] = a*A.m[i];
        return C;
}

inline Mat6  operator /(const Mat6 &A, const PetscScalar &a)
{
        Mat6 C;
        for(int i=0;i<36;i++)
                C.m[i] = A.m[i]/a;
        return C;
}

inline Mat6  operator /(const Mat6 &A, const Vec6 &a)
{
        Mat6 C;
        for(int i=0;i<6;i++)
          for(int j=0;j<6;j++)
                C.m[6*i+j] = A.m[6*i+j]/a.v[i];
        return C;
}

inline Vec6 operator *(const Mat6 &A, const Vec6 &x)
{
        Vec6 b;
        b.v[0]=A.m[0]*x.v[0]  + A.m[1]*x.v[1]  + A.m[2]*x.v[2] + A.m[3]*x.v[3] + A.m[4]*x.v[4]  + A.m[5]*x.v[5];
        b.v[1]=A.m[6]*x.v[0]  + A.m[7]*x.v[1]  + A.m[8]*x.v[2] + A.m[9]*x.v[3] + A.m[10]*x.v[4] + A.m[11]*x.v[5];
        b.v[2]=A.m[12]*x.v[0] + A.m[13]*x.v[1] + A.m[14]*x.v[2]+ A.m[15]*x.v[3]+ A.m[16]*x.v[4] + A.m[17]*x.v[5];
        b.v[3]=A.m[18]*x.v[0] + A.m[19]*x.v[1] + A.m[20]*x.v[2]+ A.m[21]*x.v[3]+ A.m[22]*x.v[4] + A.m[23]*x.v[5];
        b.v[4]=A.m[24]*x.v[0] + A.m[25]*x.v[1] + A.m[26]*x.v[2]+ A.m[27]*x.v[3]+ A.m[28]*x.v[4] + A.m[29]*x.v[5];
        b.v[5]=A.m[30]*x.v[0] + A.m[31]*x.v[1] + A.m[32]*x.v[2]+ A.m[33]*x.v[3]+ A.m[34]*x.v[4] + A.m[35]*x.v[5];
        return b;
}


inline void mat6_transport(Mat6 &A)
{

        PetscScalar swap;
        for (int i=0;i<6;i++)
          for (int j=0;j<=i;j++)
            {
             swap=A.m[i*6+j];
             A.m[i*6+j]=A.m[j*6+i];
             A.m[j*6+i]=swap;
            }
}

inline ostream& operator << ( ostream& out, const Mat6 & A) 
{
    for (int i=0;i<6;i++)
    {
      for (int j=0;j<6;j++)
        out<<A.m[i*6+j]<<" ";
      out<<endl;
    }
    return out;
}

#endif

