#ifndef _vec3_h_
#define _vec3_h_
#include <math.h>
#include "petsc.h"

//--------------------------------------------------------------------
//this is the basic operate of vec and mat. speed is the most important thing
typedef struct
{
	PetscScalar v[3];
}Vec3;

typedef struct
{
	PetscScalar m[9];
}Mat3;

inline void Set_Vec3_zero(Vec3 &vec)
{
	vec.v[0]=0;
	vec.v[1]=0;
	vec.v[2]=0;
}

inline void Set_Vec3(Vec3 &vec,PetscScalar value)
{
	vec.v[0]=value;
	vec.v[1]=value;
	vec.v[2]=value;
}

inline void Set_Vec3(Vec3 &vec,PetscScalar v1, PetscScalar v2,PetscScalar v3)
{
	vec.v[0]=v1;
	vec.v[1]=v2;
	vec.v[2]=v3;
}

inline Vec3 operator +(const Vec3 &vec1,const Vec3 &vec2)
{
	Vec3 sum;
	sum.v[0]=vec1.v[0]+vec2.v[0];
	sum.v[1]=vec1.v[1]+vec2.v[1];
	sum.v[2]=vec1.v[2]+vec2.v[2];
	return sum;
}

inline Vec3 operator +(const Vec3 &vec1,const PetscScalar value)
{
	Vec3 sum;
	sum.v[0]=vec1.v[0]+value;
	sum.v[1]=vec1.v[1]+value;
	sum.v[2]=vec1.v[2]+value;
	return sum;
}

inline Vec3 operator +(const PetscScalar value,const Vec3 &vec1)
{
	Vec3 sum;
	sum.v[0]=vec1.v[0]+value;
	sum.v[1]=vec1.v[1]+value;
	sum.v[2]=vec1.v[2]+value;
	return sum;
}

inline Vec3 operator -(const Vec3 &vec1,const Vec3 &vec2)
{
	Vec3 sub;
	sub.v[0]=vec1.v[0]-vec2.v[0];
	sub.v[1]=vec1.v[1]-vec2.v[1];
	sub.v[2]=vec1.v[2]-vec2.v[2];
	return sub;
}

inline Vec3 operator -(const Vec3 &vec1,const PetscScalar value)
{
	Vec3 sum;
	sum.v[0]=vec1.v[0]-value;
	sum.v[1]=vec1.v[1]-value;
	sum.v[2]=vec1.v[2]-value;
	return sum;
}

inline Vec3 operator -(const PetscScalar value,const Vec3 &vec1)
{
	Vec3 sum;
	sum.v[0]=value-vec1.v[0];
	sum.v[1]=value-vec1.v[1];
	sum.v[2]=value-vec1.v[2];
	return sum;
}

inline Vec3 operator -(const Vec3 &vec1)
{
	Vec3 sub;
	sub.v[0]=-vec1.v[0];
	sub.v[1]=-vec1.v[1];
	sub.v[2]=-vec1.v[2];
	return sub;
}

inline Vec3 operator *(const Vec3 &vec,const PetscScalar &value)
{
	Vec3 result;
	result.v[0]=vec.v[0]*value;
	result.v[1]=vec.v[1]*value;
	result.v[2]=vec.v[2]*value;
	return result;
}

inline Vec3 operator *(const PetscScalar &value,const Vec3 &vec)
{
	Vec3 result;
	result.v[0]=vec.v[0]*value;
	result.v[1]=vec.v[1]*value;
	result.v[2]=vec.v[2]*value;
	return result;
}

inline Vec3 operator /(const Vec3 &vec,const PetscScalar &value)
{
	Vec3 result;
	result.v[0]=vec.v[0]/value;
	result.v[1]=vec.v[1]/value;
	result.v[2]=vec.v[2]/value;
	return result;
}

inline Vec3 operator /(const PetscScalar &value,const Vec3 &vec)
{
	Vec3 result;
	result.v[0]=value/vec.v[0];
	result.v[1]=value/vec.v[1];
	result.v[2]=value/vec.v[2];
	return result;
}

inline Vec3 operator *(const Vec3 &vec1,const Vec3 &vec2)
{
	Vec3 result;
	result.v[0]=vec1.v[0]*vec2.v[0];
	result.v[1]=vec1.v[1]*vec2.v[1];
	result.v[2]=vec1.v[2]*vec2.v[2];
	return result;
}

inline Vec3 operator /(const Vec3 &vec1,const Vec3 &vec2)
{
	Vec3 result;
	result.v[0]=vec1.v[0]/vec2.v[0];
	result.v[1]=vec1.v[1]/vec2.v[1];
	result.v[2]=vec1.v[2]/vec2.v[2];
	return result;
}


//----------------------------------------------------
inline void  Set_Mat3_zero(Mat3 &A)
{
	 A.m[0]=0;	 A.m[1]=0;	 A.m[2]=0;
	 A.m[3]=0;	 A.m[4]=0;	 A.m[5]=0;
	 A.m[6]=0;	 A.m[7]=0;	 A.m[8]=0;
}

inline void Set_Mat3_I(Mat3 &A)
{
	 A.m[0]=1;	 A.m[1]=0;	 A.m[2]=0;
	 A.m[3]=0;	 A.m[4]=1;	 A.m[5]=0;
	 A.m[6]=0;	 A.m[7]=0;	 A.m[8]=1;

}

inline Mat3 Set_Mat3_Scale(const Mat3 &A, const Vec3 &a)
{
	Mat3 C;
	C.m[0]=A.m[0]*a.v[0];    C.m[1]=A.m[1]*a.v[0];    C.m[2]=A.m[2]*a.v[0];
	C.m[3]=A.m[3]*a.v[1];    C.m[4]=A.m[4]*a.v[1];    C.m[5]=A.m[5]*a.v[1];
	C.m[6]=A.m[6]*a.v[2];    C.m[7]=A.m[7]*a.v[2];    C.m[8]=A.m[8]*a.v[2];
	return C;

}

inline Mat3 operator +(const Mat3 &A, const Mat3 &B)
{
	Mat3 C;
	C.m[0]=A.m[0]+B.m[0];    C.m[1]=A.m[1]+B.m[1];    C.m[2]=A.m[2]+B.m[2];
	C.m[3]=A.m[3]+B.m[3];	 C.m[4]=A.m[4]+B.m[4];    C.m[5]=A.m[5]+B.m[5];
	C.m[6]=A.m[6]+B.m[6];    C.m[7]=A.m[7]+B.m[7];	  C.m[8]=A.m[8]+B.m[8];
	return C;
}

inline Mat3 operator -(const Mat3 &A, const Mat3 &B)
{
	Mat3 C;
	C.m[0]=A.m[0]-B.m[0];    C.m[1]=A.m[1]-B.m[1];    C.m[2]=A.m[2]-B.m[2];
	C.m[3]=A.m[3]-B.m[3];	 C.m[4]=A.m[4]-B.m[4];    C.m[5]=A.m[5]-B.m[5];
	C.m[6]=A.m[6]-B.m[6];    C.m[7]=A.m[7]-B.m[7];	  C.m[8]=A.m[8]-B.m[8];
	return C;
}

inline Mat3  operator -(const Mat3 &A)
{
	Mat3 C;
	C.m[0]=-A.m[0];    C.m[1]=-A.m[1];    C.m[2]=-A.m[2];
	C.m[3]=-A.m[3];    C.m[4]=-A.m[4];    C.m[5]=-A.m[5];
	C.m[6]=-A.m[6];    C.m[7]=-A.m[7];    C.m[8]=-A.m[8];
	return C;
}

inline Mat3  operator *(const Mat3 &A,const PetscScalar &a)
{
	Mat3 C;
	C.m[0]=a*A.m[0];    C.m[1]=a*A.m[1];    C.m[2]=a*A.m[2];
	C.m[3]=a*A.m[3];    C.m[4]=a*A.m[4];    C.m[5]=a*A.m[5];
	C.m[6]=a*A.m[6];    C.m[7]=a*A.m[7];	C.m[8]=a*A.m[8];
	return C;
}

inline Mat3  operator *(const PetscScalar &a,const Mat3 &A)
{
	Mat3 C;
	C.m[0]=a*A.m[0];    C.m[1]=a*A.m[1];    C.m[2]=a*A.m[2];
	C.m[3]=a*A.m[3];    C.m[4]=a*A.m[4];    C.m[5]=a*A.m[5];
	C.m[6]=a*A.m[6];    C.m[7]=a*A.m[7];	C.m[8]=a*A.m[8];
	return C;
}

inline Mat3  operator /(const Mat3 &A, const PetscScalar &a)
{
	Mat3 C;
	C.m[0]=A.m[0]/a;    C.m[1]=A.m[1]/a;    C.m[2]=A.m[2]/a;
	C.m[3]=A.m[3]/a;    C.m[4]=A.m[4]/a;    C.m[5]=A.m[5]/a;
	C.m[6]=A.m[6]/a;    C.m[7]=A.m[7]/a;	C.m[8]=A.m[8]/a;
	return C;
}

inline Mat3  operator /(const Mat3 &A, const Vec3 &a)
{
	Mat3 C;
	C.m[0]=A.m[0]/a.v[0];    C.m[1]=A.m[1]/a.v[0];    C.m[2]=A.m[2]/a.v[0];
	C.m[3]=A.m[3]/a.v[1];    C.m[4]=A.m[4]/a.v[1];    C.m[5]=A.m[5]/a.v[1];
	C.m[6]=A.m[6]/a.v[2];    C.m[7]=A.m[7]/a.v[2];    C.m[8]=A.m[8]/a.v[2];
	return C;
}

inline Vec3 operator *(const Mat3 &A, const Vec3 &x)
{
	Vec3 b;
	b.v[0]=A.m[0]*x.v[0] + A.m[1]*x.v[1] + A.m[2]*x.v[2];
	b.v[1]=A.m[3]*x.v[0] + A.m[4]*x.v[1] + A.m[5]*x.v[2];
	b.v[2]=A.m[6]*x.v[0] + A.m[7]*x.v[1] + A.m[8]*x.v[2];
	return b;
}

inline Mat3 operator *(const Mat3 &A, const Mat3 &B)
{
	Mat3 C;
	C.m[0]=A.m[0]*B.m[0] + A.m[1]*B.m[3] + A.m[2]*B.m[6];
	C.m[1]=A.m[0]*B.m[1] + A.m[1]*B.m[4] + A.m[2]*B.m[7];
	C.m[2]=A.m[0]*B.m[2] + A.m[1]*B.m[5] + A.m[2]*B.m[8];

	C.m[3]=A.m[3]*B.m[0] + A.m[4]*B.m[3] + A.m[5]*B.m[6];
	C.m[4]=A.m[3]*B.m[1] + A.m[4]*B.m[4] + A.m[5]*B.m[7];
	C.m[5]=A.m[3]*B.m[2] + A.m[4]*B.m[5] + A.m[5]*B.m[8];

	C.m[6]=A.m[6]*B.m[0] + A.m[7]*B.m[3] + A.m[8]*B.m[6];
	C.m[7]=A.m[6]*B.m[1] + A.m[7]*B.m[4] + A.m[8]*B.m[7];
	C.m[8]=A.m[6]*B.m[2] + A.m[7]*B.m[5] + A.m[8]*B.m[8];
	return C;
}


inline void mat3_transport(Mat3 &A)
{

        PetscScalar swap;
	for (int i=0;i<3;i++)
	  for (int j=0;j<=i;j++)
	    {
	     swap=A.m[i*3+j];
	     A.m[i*3+j]=A.m[j*3+i];
	     A.m[j*3+i]=swap;
	    }
}

inline ostream& operator << ( ostream& out, const Mat3 & A) 
{
    for (int i=0;i<3;i++)
    {
      for (int j=0;j<3;j++)
        out<<A.m[i*3+j]<<" ";
      out<<endl;
    }
    return out;
}

#endif

