#ifndef _vec5_h_
#define _vec5_h_
#include <math.h>
#include "petsc.h"

//this is the basic operate of vec and mat. speed is the most important thing
typedef struct
{
  PetscScalar v[5];
}
Vec5;

typedef struct
{
  PetscScalar m[25];
}
Mat5;

inline void Set_Vec5_zero(Vec5 &vec)
{
  vec.v[0]=0;
  vec.v[1]=0;
  vec.v[2]=0;
  vec.v[3]=0;
  vec.v[4]=0;
}

inline void Set_Vec5(Vec5 &vec,PetscScalar value)
{
  vec.v[0]=value;
  vec.v[1]=value;
  vec.v[2]=value;
  vec.v[3]=value;
  vec.v[4]=value;
}

inline void Set_Vec5(Vec5 &vec,PetscScalar v1,PetscScalar v2,PetscScalar v3,PetscScalar v4, PetscScalar v5)
{
  vec.v[0]=v1;
  vec.v[1]=v2;
  vec.v[2]=v3;
  vec.v[3]=v4;
  vec.v[4]=v5;
}

inline Vec5 operator +(const Vec5 &vec1,const Vec5 &vec2)
{
  Vec5 sum;
  sum.v[0]=vec1.v[0]+vec2.v[0];
  sum.v[1]=vec1.v[1]+vec2.v[1];
  sum.v[2]=vec1.v[2]+vec2.v[2];
  sum.v[3]=vec1.v[3]+vec2.v[3];
  sum.v[4]=vec1.v[4]+vec2.v[4];
  return sum;
}

inline Vec5 operator +(const Vec5 &vec1,const PetscScalar value)
{
  Vec5 sum;
  sum.v[0]=vec1.v[0]+value;
  sum.v[1]=vec1.v[1]+value;
  sum.v[2]=vec1.v[2]+value;
  sum.v[3]=vec1.v[3]+value;
  sum.v[4]=vec1.v[4]+value;
  return sum;
}

inline Vec5 operator +(const PetscScalar value,const Vec5 &vec1)
{
  Vec5 sum;
  sum.v[0]=vec1.v[0]+value;
  sum.v[1]=vec1.v[1]+value;
  sum.v[2]=vec1.v[2]+value;
  sum.v[3]=vec1.v[3]+value;
  sum.v[4]=vec1.v[4]+value;
  return sum;
}

inline Vec5 operator -(const Vec5 &vec1,const Vec5 &vec2)
{
  Vec5 sub;
  sub.v[0]=vec1.v[0]-vec2.v[0];
  sub.v[1]=vec1.v[1]-vec2.v[1];
  sub.v[2]=vec1.v[2]-vec2.v[2];
  sub.v[3]=vec1.v[3]-vec2.v[3];
  sub.v[4]=vec1.v[4]-vec2.v[4];
  return sub;
}

inline Vec5 operator -(const Vec5 &vec1,const PetscScalar value)
{
  Vec5 sum;
  sum.v[0]=vec1.v[0]-value;
  sum.v[1]=vec1.v[1]-value;
  sum.v[2]=vec1.v[2]-value;
  sum.v[3]=vec1.v[3]-value;
  sum.v[4]=vec1.v[4]-value;
  return sum;
}

inline Vec5 operator -(const PetscScalar value,const Vec5 &vec1)
{
  Vec5 sum;
  sum.v[0]=value-vec1.v[0];
  sum.v[1]=value-vec1.v[1];
  sum.v[2]=value-vec1.v[2];
  sum.v[3]=value-vec1.v[3];
  sum.v[4]=value-vec1.v[4];
  return sum;
}

inline Vec5 operator *(const Vec5 &vec,const PetscScalar &value)
{
  Vec5 result;
  result.v[0]=vec.v[0]*value;
  result.v[1]=vec.v[1]*value;
  result.v[2]=vec.v[2]*value;
  result.v[3]=vec.v[3]*value;
  result.v[4]=vec.v[4]*value;
  return result;
}

inline Vec5 operator *(const PetscScalar &value,const Vec5 &vec)
{
  Vec5 result;
  result.v[0]=vec.v[0]*value;
  result.v[1]=vec.v[1]*value;
  result.v[2]=vec.v[2]*value;
  result.v[3]=vec.v[3]*value;
  result.v[4]=vec.v[4]*value;
  return result;
}

inline Vec5 operator /(const Vec5 &vec,const PetscScalar &value)
{
  Vec5 result;
  result.v[0]=vec.v[0]/value;
  result.v[1]=vec.v[1]/value;
  result.v[2]=vec.v[2]/value;
  result.v[3]=vec.v[3]/value;
  result.v[4]=vec.v[4]/value;
  return result;
}

inline Vec5 operator /(const PetscScalar &value,const Vec5 &vec)
{
  Vec5 result;
  result.v[0]=value/vec.v[0];
  result.v[1]=value/vec.v[1];
  result.v[2]=value/vec.v[2];
  result.v[3]=value/vec.v[3];
  result.v[4]=value/vec.v[4];
  return result;
}

inline Vec5 operator *(const Vec5 &vec1,const Vec5 &vec2)
{
  Vec5 result;
  result.v[0]=vec1.v[0]*vec2.v[0];
  result.v[1]=vec1.v[1]*vec2.v[1];
  result.v[2]=vec1.v[2]*vec2.v[2];
  result.v[3]=vec1.v[3]*vec2.v[3];
  result.v[4]=vec1.v[4]*vec2.v[4];
  return result;
}

inline Vec5 operator /(const Vec5 &vec1,const Vec5 &vec2)
{
  Vec5 result;
  result.v[0]=vec1.v[0]/vec2.v[0];
  result.v[1]=vec1.v[1]/vec2.v[1];
  result.v[2]=vec1.v[2]/vec2.v[2];
  result.v[3]=vec1.v[3]/vec2.v[3];
  result.v[4]=vec1.v[4]/vec2.v[4];
  return result;
}

//--------------------------------------------------------------------

inline PetscScalar vmaxvalue(const Vec5 &a)
{
  PetscScalar v=a.v[0];
  if(a.v[1]>v)  v=a.v[1];
  if(a.v[2]>v)  v=a.v[2];
  if(a.v[3]>v)  v=a.v[3];
  if(a.v[4]>v)  v=a.v[4];
  return v;
}

inline Vec5 vmax(const Vec5 &a,const Vec5 &b)
{
  Vec5 result;
  result.v[0]=a.v[0]>b.v[0] ? a.v[0]:b.v[0];
  result.v[1]=a.v[1]>b.v[1] ? a.v[1]:b.v[1];
  result.v[2]=a.v[2]>b.v[2] ? a.v[2]:b.v[2];
  result.v[3]=a.v[3]>b.v[3] ? a.v[3]:b.v[3];
  result.v[4]=a.v[4]>b.v[4] ? a.v[4]:b.v[4];
  return result;
}

inline Vec5 vmin(const Vec5 &a,const Vec5 &b)
{
  Vec5 result;
  result.v[0]=a.v[0]<b.v[0] ? a.v[0]:b.v[0];
  result.v[1]=a.v[1]<b.v[1] ? a.v[1]:b.v[1];
  result.v[2]=a.v[2]<b.v[2] ? a.v[2]:b.v[2];
  result.v[3]=a.v[3]<b.v[3] ? a.v[3]:b.v[3];
  result.v[4]=a.v[4]<b.v[4] ? a.v[4]:b.v[4];
  return result;
}

inline Vec5 vabs(const Vec5 &a)
{
  Vec5 result;
  result.v[0]=fabs(a.v[0]);
  result.v[1]=fabs(a.v[1]);
  result.v[2]=fabs(a.v[2]);
  result.v[3]=fabs(a.v[3]);
  result.v[4]=fabs(a.v[4]);
  return result;
}

inline Vec5 vsign(const Vec5 &a)
{
  Vec5 result;
  result.v[0]=a.v[0]>0.0? 1.0:-1.0;
  result.v[1]=a.v[1]>0.0? 1.0:-1.0;
  result.v[2]=a.v[2]>0.0? 1.0:-1.0;
  result.v[3]=a.v[3]>0.0? 1.0:-1.0;
  result.v[4]=a.v[4]>0.0? 1.0:-1.0;
  return result;
}

//-------------------------------------------------------------------------

inline void  Set_Mat5_zero(Mat5 &A)
{
  A.m[0]=0;	 A.m[1]=0;	 A.m[2]=0;	 A.m[3]=0;	 A.m[4]=0;
  A.m[5]=0;	 A.m[6]=0;	 A.m[7]=0;	 A.m[8]=0;	 A.m[9]=0;
  A.m[10]=0;	 A.m[11]=0;	 A.m[12]=0;	 A.m[13]=0;	 A.m[14]=0;
  A.m[15]=0;	 A.m[16]=0;	 A.m[17]=0;	 A.m[18]=0;	 A.m[19]=0;
  A.m[20]=0;	 A.m[21]=0;	 A.m[22]=0;	 A.m[23]=0;	 A.m[24]=0;
}

inline void Set_Mat5_I(Mat5 &A)
{
  A.m[0]=1;	 A.m[1]=0;	 A.m[2]=0;	 A.m[3]=0;	 A.m[4]=0;
  A.m[5]=0;	 A.m[6]=1;	 A.m[7]=0;	 A.m[8]=0;	 A.m[9]=0;
  A.m[10]=0;	 A.m[11]=0;	 A.m[12]=1;	 A.m[13]=0;	 A.m[14]=0;
  A.m[15]=0;	 A.m[16]=0;	 A.m[17]=0;	 A.m[18]=1;	 A.m[19]=0;
  A.m[20]=0;	 A.m[21]=0;	 A.m[22]=0;	 A.m[23]=0;	 A.m[24]=1;

}

inline Mat5 operator +(const Mat5 &A, const Mat5 &B)
{
  Mat5 C;
  C.m[0]=A.m[0]+B.m[0];    C.m[1]=A.m[1]+B.m[1];    C.m[2]=A.m[2]+B.m[2];    C.m[3]=A.m[3]+B.m[3];    C.m[4]=A.m[4]+B.m[4];
  C.m[5]=A.m[5]+B.m[5];    C.m[6]=A.m[6]+B.m[6];    C.m[7]=A.m[7]+B.m[7];    C.m[8]=A.m[8]+B.m[8];    C.m[9]=A.m[9]+B.m[9];
  C.m[10]=A.m[10]+B.m[10]; C.m[11]=A.m[11]+B.m[11]; C.m[12]=A.m[12]+B.m[12]; C.m[13]=A.m[13]+B.m[13]; C.m[14]=A.m[14]+B.m[14];
  C.m[15]=A.m[15]+B.m[15]; C.m[16]=A.m[16]+B.m[16]; C.m[17]=A.m[17]+B.m[17]; C.m[18]=A.m[18]+B.m[18]; C.m[19]=A.m[19]+B.m[19];
  C.m[20]=A.m[20]+B.m[20]; C.m[21]=A.m[21]+B.m[21]; C.m[22]=A.m[22]+B.m[22]; C.m[23]=A.m[23]+B.m[23]; C.m[24]=A.m[24]+B.m[24];
  return C;
}

inline Mat5 operator -(const Mat5 &A, const Mat5 &B)
{
  Mat5 C;
  C.m[0]=A.m[0]-B.m[0];    C.m[1]=A.m[1]-B.m[1];    C.m[2]=A.m[2]-B.m[2];    C.m[3]=A.m[3]-B.m[3];    C.m[4]=A.m[4]-B.m[4];
  C.m[5]=A.m[5]-B.m[5];    C.m[6]=A.m[6]-B.m[6];    C.m[7]=A.m[7]-B.m[7];    C.m[8]=A.m[8]-B.m[8];    C.m[9]=A.m[9]-B.m[9];
  C.m[10]=A.m[10]-B.m[10]; C.m[11]=A.m[11]-B.m[11]; C.m[12]=A.m[12]-B.m[12]; C.m[13]=A.m[13]-B.m[13]; C.m[14]=A.m[14]-B.m[14];
  C.m[15]=A.m[15]-B.m[15]; C.m[16]=A.m[16]-B.m[16]; C.m[17]=A.m[17]-B.m[17]; C.m[18]=A.m[18]-B.m[18]; C.m[19]=A.m[19]-B.m[19];
  C.m[20]=A.m[20]-B.m[20]; C.m[21]=A.m[21]-B.m[21]; C.m[22]=A.m[22]-B.m[22]; C.m[23]=A.m[23]-B.m[23]; C.m[24]=A.m[24]-B.m[24];
  return C;
}

inline Mat5  operator *(const Mat5 &A,const PetscScalar &a)
{
  Mat5 C;
  C.m[0]=a*A.m[0];    C.m[1]=a*A.m[1];    C.m[2]=a*A.m[2];    C.m[3]=a*A.m[3];    C.m[4]=a*A.m[4]; 
  C.m[5]=a*A.m[5];    C.m[6]=a*A.m[6];    C.m[7]=a*A.m[7];    C.m[8]=a*A.m[8];    C.m[9]=a*A.m[9];
  C.m[10]=a*A.m[10];  C.m[11]=a*A.m[11];  C.m[12]=a*A.m[12];  C.m[13]=a*A.m[13];  C.m[14]=a*A.m[14];
  C.m[15]=a*A.m[15];  C.m[16]=a*A.m[16];  C.m[17]=a*A.m[17];  C.m[18]=a*A.m[18];  C.m[19]=a*A.m[19];
  C.m[20]=a*A.m[20];  C.m[21]=a*A.m[21];  C.m[22]=a*A.m[22];  C.m[23]=a*A.m[23];  C.m[24]=a*A.m[24];
  return C;
}

inline Mat5  operator *(const PetscScalar &a,const Mat5 &A)
{
  Mat5 C;
  C.m[0]=a*A.m[0];    C.m[1]=a*A.m[1];    C.m[2]=a*A.m[2];    C.m[3]=a*A.m[3];    C.m[4]=a*A.m[4]; 
  C.m[5]=a*A.m[5];    C.m[6]=a*A.m[6];    C.m[7]=a*A.m[7];    C.m[8]=a*A.m[8];    C.m[9]=a*A.m[9];
  C.m[10]=a*A.m[10];  C.m[11]=a*A.m[11];  C.m[12]=a*A.m[12];  C.m[13]=a*A.m[13];  C.m[14]=a*A.m[14];
  C.m[15]=a*A.m[15];  C.m[16]=a*A.m[16];  C.m[17]=a*A.m[17];  C.m[18]=a*A.m[18];  C.m[19]=a*A.m[19];
  C.m[20]=a*A.m[20];  C.m[21]=a*A.m[21];  C.m[22]=a*A.m[22];  C.m[23]=a*A.m[23];  C.m[24]=a*A.m[24];
  return C;
}

inline Mat5  operator /(const Mat5 &A, const PetscScalar &a)
{
  Mat5 C;
  C.m[0]=A.m[0]/a;    C.m[1]=A.m[1]/a;    C.m[2]=A.m[2]/a;    C.m[3]=A.m[3]/a;    C.m[4]=A.m[4]/a;
  C.m[5]=A.m[5]/a;    C.m[6]=A.m[6]/a;    C.m[7]=A.m[7]/a;    C.m[8]=A.m[8]/a;    C.m[9]=A.m[9]/a;
  C.m[10]=A.m[10]/a;  C.m[11]=A.m[11]/a;  C.m[12]=A.m[12]/a;  C.m[13]=A.m[13]/a;  C.m[14]=A.m[14]/a;  
  C.m[15]=A.m[15]/a;  C.m[16]=A.m[16]/a;  C.m[17]=A.m[17]/a;  C.m[18]=A.m[18]/a;  C.m[19]=A.m[19]/a;  
  C.m[20]=A.m[20]/a;  C.m[21]=A.m[21]/a;  C.m[22]=A.m[22]/a;  C.m[23]=A.m[23]/a;  C.m[24]=A.m[24]/a;  
  return C;
}

inline Mat5  operator /(const Mat5 &A, const Vec5 &a)
{
  Mat5 C;
  C.m[0]=A.m[0]/a.v[0];    C.m[1]=A.m[1]/a.v[0];    C.m[2]=A.m[2]/a.v[0];    C.m[3]=A.m[3]/a.v[0];    C.m[4]=A.m[4]/a.v[0];    
  C.m[5]=A.m[5]/a.v[1];    C.m[6]=A.m[6]/a.v[1];    C.m[7]=A.m[7]/a.v[1];    C.m[8]=A.m[8]/a.v[1];    C.m[9]=A.m[9]/a.v[1]; 
  C.m[10]=A.m[10]/a.v[2];  C.m[11]=A.m[11]/a.v[2];  C.m[12]=A.m[12]/a.v[2];  C.m[13]=A.m[13]/a.v[2];  C.m[14]=A.m[14]/a.v[2];
  C.m[15]=A.m[15]/a.v[3];  C.m[16]=A.m[16]/a.v[3];  C.m[17]=A.m[17]/a.v[3];  C.m[18]=A.m[18]/a.v[3];  C.m[19]=A.m[19]/a.v[3];
  C.m[20]=A.m[20]/a.v[4];  C.m[21]=A.m[21]/a.v[4];  C.m[22]=A.m[22]/a.v[4];  C.m[23]=A.m[23]/a.v[4];  C.m[24]=A.m[24]/a.v[4];
  return C;
}

inline Vec5 operator *(const Mat5 &A, const Vec5 &x)
{
  Vec5 b;
  b.v[0]=A.m[0]*x.v[0]+A.m[1]*x.v[1]+A.m[2]*x.v[2]+A.m[3]*x.v[3]+A.m[4]*x.v[4];
  b.v[1]=A.m[5]*x.v[0]+A.m[6]*x.v[1]+A.m[7]*x.v[2]+A.m[8]*x.v[3]+A.m[9]*x.v[4];
  b.v[2]=A.m[10]*x.v[0]+A.m[11]*x.v[1]+A.m[12]*x.v[2]+A.m[13]*x.v[3]+A.m[14]*x.v[4];
  b.v[3]=A.m[15]*x.v[0]+A.m[16]*x.v[1]+A.m[17]*x.v[2]+A.m[18]*x.v[3]+A.m[19]*x.v[4];
  b.v[4]=A.m[20]*x.v[0]+A.m[21]*x.v[1]+A.m[22]*x.v[2]+A.m[23]*x.v[3]+A.m[24]*x.v[4];
  return b;
}


inline void mat5_transport(Mat5 &A)
{
  PetscScalar swap;
  for (int i=0;i<5;i++)
    for (int j=0;j<=i;j++)
    {
      swap=A.m[i*5+j];
      A.m[i*5+j]=A.m[j*5+i];
      A.m[j*5+i]=swap;
    }
}

inline ostream& operator << ( ostream& out, const Mat5 & A) 
{
    for (int i=0;i<5;i++)
    {
      for (int j=0;j<5;j++)
        out<<A.m[i*5+j]<<" ";
      out<<endl;
    }
    return out;
}

#endif

