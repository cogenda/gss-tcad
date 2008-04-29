#ifndef _vec4_h_
#define _vec4_h_
#include <math.h>
#include "petsc.h"

//this is the basic operate of vec and mat. speed is the most important thing
typedef struct
{
	PetscScalar v[4];
}Vec4;

typedef struct
{
	PetscScalar m[16];
}Mat4;

inline void Set_Vec4_zero(Vec4 &vec)
{
	vec.v[0]=0;
	vec.v[1]=0;
	vec.v[2]=0;
	vec.v[3]=0;
}

inline void Set_Vec4(Vec4 &vec,PetscScalar value)
{
	vec.v[0]=value;
	vec.v[1]=value;
	vec.v[2]=value;
	vec.v[3]=value;
}

inline void Set_Vec4(Vec4 &vec,PetscScalar v1,PetscScalar v2,PetscScalar v3,PetscScalar v4)
{
	vec.v[0]=v1;
	vec.v[1]=v2;
	vec.v[2]=v3;
	vec.v[3]=v4;
}

inline Vec4 operator +(const Vec4 &vec1,const Vec4 &vec2)
{
	Vec4 sum;
	sum.v[0]=vec1.v[0]+vec2.v[0];
	sum.v[1]=vec1.v[1]+vec2.v[1];
	sum.v[2]=vec1.v[2]+vec2.v[2];
	sum.v[3]=vec1.v[3]+vec2.v[3];
	return sum;
}

inline Vec4 operator +(const Vec4 &vec1,const PetscScalar value)
{
	Vec4 sum;
	sum.v[0]=vec1.v[0]+value;
	sum.v[1]=vec1.v[1]+value;
	sum.v[2]=vec1.v[2]+value;
	sum.v[3]=vec1.v[3]+value;
	return sum;
}

inline Vec4 operator +(const PetscScalar value,const Vec4 &vec1)
{
	Vec4 sum;
	sum.v[0]=vec1.v[0]+value;
	sum.v[1]=vec1.v[1]+value;
	sum.v[2]=vec1.v[2]+value;
	sum.v[3]=vec1.v[3]+value;
	return sum;
}

inline Vec4 operator -(const Vec4 &vec1,const Vec4 &vec2)
{
	Vec4 sub;
	sub.v[0]=vec1.v[0]-vec2.v[0];
	sub.v[1]=vec1.v[1]-vec2.v[1];
	sub.v[2]=vec1.v[2]-vec2.v[2];
	sub.v[3]=vec1.v[3]-vec2.v[3];
	return sub;
}

inline Vec4 operator -(const Vec4 &vec1,const PetscScalar value)
{
	Vec4 sum;
	sum.v[0]=vec1.v[0]-value;
	sum.v[1]=vec1.v[1]-value;
	sum.v[2]=vec1.v[2]-value;
	sum.v[3]=vec1.v[3]-value;
	return sum;
}

inline Vec4 operator -(const PetscScalar value,const Vec4 &vec1)
{
	Vec4 sum;
	sum.v[0]=value-vec1.v[0];
	sum.v[1]=value-vec1.v[1];
	sum.v[2]=value-vec1.v[2];
	sum.v[3]=value-vec1.v[3];
	return sum;
}

inline Vec4 operator -(const Vec4 &vec1)
{
	Vec4 sub;
	sub.v[0]=-vec1.v[0];
	sub.v[1]=-vec1.v[1];
	sub.v[2]=-vec1.v[2];
	sub.v[3]=-vec1.v[3];
	return sub;
}

inline Vec4 operator *(const Vec4 &vec,const PetscScalar &value)
{
	Vec4 result;
	result.v[0]=vec.v[0]*value;
	result.v[1]=vec.v[1]*value;
	result.v[2]=vec.v[2]*value;
	result.v[3]=vec.v[3]*value;
	return result;
}

inline Vec4 operator *(const PetscScalar &value,const Vec4 &vec)
{
	Vec4 result;
	result.v[0]=vec.v[0]*value;
	result.v[1]=vec.v[1]*value;
	result.v[2]=vec.v[2]*value;
	result.v[3]=vec.v[3]*value;
	return result;
}

inline Vec4 operator /(const Vec4 &vec,const PetscScalar &value)
{
	Vec4 result;
	result.v[0]=vec.v[0]/value;
	result.v[1]=vec.v[1]/value;
	result.v[2]=vec.v[2]/value;
	result.v[3]=vec.v[3]/value;
	return result;
}

inline Vec4 operator /(const PetscScalar &value,const Vec4 &vec)
{
	Vec4 result;
	result.v[0]=value/vec.v[0];
	result.v[1]=value/vec.v[1];
	result.v[2]=value/vec.v[2];
	result.v[3]=value/vec.v[3];
	return result;
}

inline Vec4 operator *(const Vec4 &vec1,const Vec4 &vec2)
{
	Vec4 result;
	result.v[0]=vec1.v[0]*vec2.v[0];
	result.v[1]=vec1.v[1]*vec2.v[1];
	result.v[2]=vec1.v[2]*vec2.v[2];
	result.v[3]=vec1.v[3]*vec2.v[3];
	return result;
}

inline Vec4 operator /(const Vec4 &vec1,const Vec4 &vec2)
{
	Vec4 result;
	result.v[0]=vec1.v[0]/vec2.v[0];
	result.v[1]=vec1.v[1]/vec2.v[1];
	result.v[2]=vec1.v[2]/vec2.v[2];
	result.v[3]=vec1.v[3]/vec2.v[3];
	return result;
}

//--------------------------------------------------------------------
inline PetscScalar vmaxvalue(const Vec4 &a)
{
	PetscScalar v=a.v[0];
	if(a.v[1]>v)  v=a.v[1];
	if(a.v[2]>v)  v=a.v[2];
	if(a.v[3]>v)  v=a.v[3];
	return v;
}

inline Vec4 vmax(const Vec4 &a,const Vec4 &b)
{
	Vec4 result;
 	result.v[0]=a.v[0]>b.v[0] ? a.v[0]:b.v[0];
  	result.v[1]=a.v[1]>b.v[1] ? a.v[1]:b.v[1];
  	result.v[2]=a.v[2]>b.v[2] ? a.v[2]:b.v[2];
  	result.v[3]=a.v[3]>b.v[3] ? a.v[3]:b.v[3];
  	return result;
}

inline Vec4 vmin(const Vec4 &a,const Vec4 &b)
{
	Vec4 result;
  	result.v[0]=a.v[0]<b.v[0] ? a.v[0]:b.v[0];
  	result.v[1]=a.v[1]<b.v[1] ? a.v[1]:b.v[1];
  	result.v[2]=a.v[2]<b.v[2] ? a.v[2]:b.v[2];
  	result.v[3]=a.v[3]<b.v[3] ? a.v[3]:b.v[3];
  	return result;
}

inline Vec4 vabs(const Vec4 &a)
{
	Vec4 result;
	result.v[0]=fabs(a.v[0]);
	result.v[1]=fabs(a.v[1]);
	result.v[2]=fabs(a.v[2]);
	result.v[3]=fabs(a.v[3]);
        return result;
}

inline Vec4 vsign(const Vec4 &a)
{
	Vec4 result;
	result.v[0]=a.v[0]>0.0? 1.0:-1.0;
	result.v[1]=a.v[1]>0.0? 1.0:-1.0;
	result.v[2]=a.v[2]>0.0? 1.0:-1.0;
	result.v[3]=a.v[3]>0.0? 1.0:-1.0;
	return result;
}

//----------------------------------------------------
inline void  Set_Mat4_zero(Mat4 &A)
{
	 A.m[0]=0;	 A.m[1]=0;	 A.m[2]=0;	 A.m[3]=0;
	 A.m[4]=0;	 A.m[5]=0;	 A.m[6]=0;	 A.m[7]=0;
	 A.m[8]=0;	 A.m[9]=0;	 A.m[10]=0;	 A.m[11]=0;
	 A.m[12]=0;	 A.m[13]=0;	 A.m[14]=0;	 A.m[15]=0;
}

inline void Set_Mat4_I(Mat4 &A)
{
	 A.m[0]=1;	 A.m[1]=0;	 A.m[2]=0;	 A.m[3]=0;
	 A.m[4]=0;	 A.m[5]=1;	 A.m[6]=0;	 A.m[7]=0;
	 A.m[8]=0;	 A.m[9]=0;	 A.m[10]=1;	 A.m[11]=0;
	 A.m[12]=0;	 A.m[13]=0;	 A.m[14]=0;	 A.m[15]=1;

}

inline Mat4 operator +(const Mat4 &A, const Mat4 &B)
{
	Mat4 C;
	C.m[0]=A.m[0]+B.m[0];    C.m[1]=A.m[1]+B.m[1];    C.m[2]=A.m[2]+B.m[2];    C.m[3]=A.m[3]+B.m[3];
	C.m[4]=A.m[4]+B.m[4];    C.m[5]=A.m[5]+B.m[5];    C.m[6]=A.m[6]+B.m[6];    C.m[7]=A.m[7]+B.m[7];
	C.m[8]=A.m[8]+B.m[8];    C.m[9]=A.m[9]+B.m[9];    C.m[10]=A.m[10]+B.m[10]; C.m[11]=A.m[11]+B.m[11];
	C.m[12]=A.m[12]+B.m[12]; C.m[13]=A.m[13]+B.m[13]; C.m[14]=A.m[14]+B.m[14]; C.m[15]=A.m[15]+B.m[15];
	return C;
}

inline Mat4 operator -(const Mat4 &A, const Mat4 &B)
{
	Mat4 C;
	C.m[0]=A.m[0]-B.m[0];    C.m[1]=A.m[1]-B.m[1];    C.m[2]=A.m[2]-B.m[2];    C.m[3]=A.m[3]-B.m[3];
	C.m[4]=A.m[4]-B.m[4];    C.m[5]=A.m[5]-B.m[5];    C.m[6]=A.m[6]-B.m[6];    C.m[7]=A.m[7]-B.m[7];
	C.m[8]=A.m[8]-B.m[8];    C.m[9]=A.m[9]-B.m[9];    C.m[10]=A.m[10]-B.m[10]; C.m[11]=A.m[11]-B.m[11];
	C.m[12]=A.m[12]-B.m[12]; C.m[13]=A.m[13]-B.m[13]; C.m[14]=A.m[14]-B.m[14]; C.m[15]=A.m[15]-B.m[15];
	return C;
}

inline Mat4  operator -(const Mat4 &A)
{
	Mat4 C;
	C.m[0]=-A.m[0];    C.m[1]=-A.m[1];    C.m[2]=-A.m[2];    C.m[3]=-A.m[3];
	C.m[4]=-A.m[4];    C.m[5]=-A.m[5];    C.m[6]=-A.m[6];    C.m[7]=-A.m[7];
	C.m[8]=-A.m[8];    C.m[9]=-A.m[9];    C.m[10]=-A.m[10];  C.m[11]=-A.m[11];
	C.m[12]=-A.m[12];  C.m[13]=-A.m[13];  C.m[14]=-A.m[14];  C.m[15]=-A.m[15];
	return C;
}


inline Mat4  operator *(const Mat4 &A,const PetscScalar &a)
{
	Mat4 C;
	C.m[0]=a*A.m[0];    C.m[1]=a*A.m[1];    C.m[2]=a*A.m[2];    C.m[3]=a*A.m[3];
	C.m[4]=a*A.m[4];    C.m[5]=a*A.m[5];    C.m[6]=a*A.m[6];    C.m[7]=a*A.m[7];
	C.m[8]=a*A.m[8];    C.m[9]=a*A.m[9];    C.m[10]=a*A.m[10];  C.m[11]=a*A.m[11];
	C.m[12]=a*A.m[12];  C.m[13]=a*A.m[13];  C.m[14]=a*A.m[14];  C.m[15]=a*A.m[15];
	return C;
}

inline Mat4  operator *(const PetscScalar &a,const Mat4 &A)
{
	Mat4 C;
	C.m[0]=a*A.m[0];    C.m[1]=a*A.m[1];    C.m[2]=a*A.m[2];    C.m[3]=a*A.m[3];
	C.m[4]=a*A.m[4];    C.m[5]=a*A.m[5];    C.m[6]=a*A.m[6];    C.m[7]=a*A.m[7];
	C.m[8]=a*A.m[8];    C.m[9]=a*A.m[9];    C.m[10]=a*A.m[10];  C.m[11]=a*A.m[11];
	C.m[12]=a*A.m[12];  C.m[13]=a*A.m[13];  C.m[14]=a*A.m[14];  C.m[15]=a*A.m[15];
	return C;
}

inline Mat4  operator /(const Mat4 &A, const PetscScalar &a)
{
	Mat4 C;
	C.m[0]=A.m[0]/a;    C.m[1]=A.m[1]/a;    C.m[2]=A.m[2]/a;    C.m[3]=A.m[3]/a;
	C.m[4]=A.m[4]/a;    C.m[5]=A.m[5]/a;    C.m[6]=A.m[6]/a;    C.m[7]=A.m[7]/a;
	C.m[8]=A.m[8]/a;    C.m[9]=A.m[9]/a;    C.m[10]=A.m[10]/a;  C.m[11]=A.m[11]/a;
	C.m[12]=A.m[12]/a;  C.m[13]=A.m[13]/a;  C.m[14]=A.m[14]/a;  C.m[15]=A.m[15]/a;
	return C;
}

inline Mat4  operator /(const Mat4 &A, const Vec4 &a)
{
	Mat4 C;
	C.m[0]=A.m[0]/a.v[0];    C.m[1]=A.m[1]/a.v[0];    C.m[2]=A.m[2]/a.v[0];    C.m[3]=A.m[3]/a.v[0];
	C.m[4]=A.m[4]/a.v[1];    C.m[5]=A.m[5]/a.v[1];    C.m[6]=A.m[6]/a.v[1];    C.m[7]=A.m[7]/a.v[1];
	C.m[8]=A.m[8]/a.v[2];    C.m[9]=A.m[9]/a.v[2];    C.m[10]=A.m[10]/a.v[2];  C.m[11]=A.m[11]/a.v[2];
	C.m[12]=A.m[12]/a.v[3];  C.m[13]=A.m[13]/a.v[3];  C.m[14]=A.m[14]/a.v[3];  C.m[15]=A.m[15]/a.v[3];
	return C;
}

inline Vec4 operator *(const Mat4 &A, const Vec4 &x)
{
	Vec4 b;
	b.v[0]=A.m[0]*x.v[0]+A.m[1]*x.v[1]+A.m[2]*x.v[2]+A.m[3]*x.v[3];
	b.v[1]=A.m[4]*x.v[0]+A.m[5]*x.v[1]+A.m[6]*x.v[2]+A.m[7]*x.v[3];
	b.v[2]=A.m[8]*x.v[0]+A.m[9]*x.v[1]+A.m[10]*x.v[2]+A.m[11]*x.v[3];
	b.v[3]=A.m[12]*x.v[0]+A.m[13]*x.v[1]+A.m[14]*x.v[2]+A.m[15]*x.v[3];
	return b;
}

inline Mat4 operator *(const Mat4 &A, const Mat4 &B)
{
	Mat4 C;
	C.m[0]=A.m[0]*B.m[0]+A.m[1]*B.m[4]+A.m[2]*B.m[8]+A.m[3]*B.m[12];
	C.m[1]=A.m[0]*B.m[1]+A.m[1]*B.m[5]+A.m[2]*B.m[9]+A.m[3]*B.m[13];
	C.m[2]=A.m[0]*B.m[2]+A.m[1]*B.m[6]+A.m[2]*B.m[10]+A.m[3]*B.m[14];
	C.m[3]=A.m[0]*B.m[3]+A.m[1]*B.m[7]+A.m[2]*B.m[11]+A.m[3]*B.m[15];
	C.m[4]=A.m[4]*B.m[0]+A.m[5]*B.m[4]+A.m[6]*B.m[8]+A.m[7]*B.m[12];
	C.m[5]=A.m[4]*B.m[1]+A.m[5]*B.m[5]+A.m[6]*B.m[9]+A.m[7]*B.m[13];
	C.m[6]=A.m[4]*B.m[2]+A.m[5]*B.m[6]+A.m[6]*B.m[10]+A.m[7]*B.m[14];
	C.m[7]=A.m[4]*B.m[3]+A.m[5]*B.m[7]+A.m[6]*B.m[11]+A.m[7]*B.m[15];
	C.m[8]=A.m[8]*B.m[0]+A.m[9]*B.m[4]+A.m[10]*B.m[8]+A.m[11]*B.m[12];
	C.m[9]=A.m[8]*B.m[1]+A.m[9]*B.m[5]+A.m[10]*B.m[9]+A.m[11]*B.m[13];
	C.m[10]=A.m[8]*B.m[2]+A.m[9]*B.m[6]+A.m[10]*B.m[10]+A.m[11]*B.m[14];
	C.m[11]=A.m[8]*B.m[3]+A.m[9]*B.m[7]+A.m[10]*B.m[11]+A.m[11]*B.m[15];
	C.m[12]=A.m[12]*B.m[0]+A.m[13]*B.m[4]+A.m[14]*B.m[8]+A.m[15]*B.m[12];
	C.m[13]=A.m[12]*B.m[1]+A.m[13]*B.m[5]+A.m[14]*B.m[9]+A.m[15]*B.m[13];
	C.m[14]=A.m[12]*B.m[2]+A.m[13]*B.m[6]+A.m[14]*B.m[10]+A.m[15]*B.m[14];
	C.m[15]=A.m[12]*B.m[3]+A.m[13]*B.m[7]+A.m[14]*B.m[11]+A.m[15]*B.m[15];
	return C;
}

inline void mat4_transport(Mat4 &A)
{

        PetscScalar swap;
	for (int i=0;i<4;i++)
	  for (int j=0;j<=i;j++)
	    {
	     swap=A.m[i*4+j];
	     A.m[i*4+j]=A.m[j*4+i];
	     A.m[j*4+i]=swap;
	    }
}

inline ostream& operator << ( ostream& out, const Mat4 & A) 
{
    for (int i=0;i<4;i++)
    {
      for (int j=0;j<4;j++)
        out<<A.m[i*4+j]<<" ";
      out<<endl;
    }
    return out;
}
#endif

