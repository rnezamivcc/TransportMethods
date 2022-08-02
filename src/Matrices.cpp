///////////////////////////////////////////////////////////////////////////////
// Matrice.cpp
// ===========
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>
#include "Matrices.h"

const float DEG2RAD = 3.141593f / 180.f;
const float EPSILON = 0.00001f;


///////////////////////////////////////////////////////////////////////////////
// inverse of 2x2 matrix
// If cannot find inverse, return zero matrix
///////////////////////////////////////////////////////////////////////////////
Matrix2 Matrix2::invert()
{
    float determinant = getDeterminant();
    if(fabs(determinant) <= EPSILON)
    {
        return zero();
    }

    float tmp = rows[0].x;   // copy the first element
    float invDeterminant = 1.0f / determinant;
	float2D row1(invDeterminant * rows[1].y, -invDeterminant * rows[0].y);
	float2D row2(-invDeterminant * rows[1].x, invDeterminant * rows[0].x);
	return { row1, row2 };
}

///////////////////////////////////////////////////////////////////////////////
// return determinant of 3x3 matrix
///////////////////////////////////////////////////////////////////////////////
float Matrix3::getDeterminant()
{
	return	rows[0][0] * (rows[1][1] * rows[2][2] - rows[1][2] * rows[2][1]) -
			rows[0][1] * (rows[1][0] * rows[2][2] - rows[1][2] * rows[2][0]) +
			rows[0][2] * (rows[1][0] * rows[2][1] - rows[1][1] * rows[2][0]);
}



///////////////////////////////////////////////////////////////////////////////
// inverse 3x3 matrix
// If cannot find inverse, set zero matrix
///////////////////////////////////////////////////////////////////////////////
Matrix3 Matrix3::invert()
{
	float tmp[9];

    tmp[0] = rows[1][1] * rows[2][2] - rows[1][2] * rows[2][1];
    tmp[1] = rows[1][0] * rows[2][2] - rows[1][2] * rows[2][0];
    tmp[2] = rows[1][0] * rows[2][1] - rows[1][1] * rows[2][0];
    tmp[3] = rows[1][0] * rows[2][2] - rows[1][2] * rows[2][0];
    tmp[4] = rows[0][0] * rows[2][2] - rows[0][2] * rows[2][0];
    tmp[5] = rows[0][0] * rows[2][1] - rows[1][0] * rows[2][0];
    tmp[6] = rows[0][1] * rows[1][2] - rows[1][1] * rows[0][2];
    tmp[7] = rows[0][0] * rows[1][2] - rows[0][2] * rows[1][0];
    tmp[8] = rows[0][0] * rows[1][1] - rows[0][1] * rows[1][0];

	float determinant = rows[0][0] * tmp[0] - rows[0][1] * tmp[1] + rows[0][2] * tmp[2];
	if (fabs(determinant) <= EPSILON)
	{
		return zero();
	}

	float invDeterminant = 1.0f / determinant;
	float3 row0 = invDeterminant * float3{ tmp[0], -tmp[3], tmp[6] };
	float3 row1 = invDeterminant * float3{ -tmp[1], tmp[4], -tmp[7] };
	float3 row2 = invDeterminant * float3{ tmp[2], -tmp[5], tmp[8] };
	return Matrix3{ row0, row1, row2 };
}


///////////////////////////////////////////////////////////////////////////////
// inverse 4x4 matrix
///////////////////////////////////////////////////////////////////////////////
Matrix4& Matrix4::invert()
{
    // If the 4th row is [0,0,0,1] then it is affine matrix and
    // it has no projective transformation.
  //  if(m[3] == 0 && m[7] == 0 && m[11] == 0 && m[15] == 1)
  //      this->invertAffine();
  //  else
    {
   //     this->invertGeneral();
        /*@@ invertProjective() is not optimized (slower than generic one)
        if(fabs(m[0]*m[5] - m[1]*m[4]) > EPSILON)
            this->invertProjective();   // inverse using matrix partition
        else
            this->invertGeneral();      // generalized inverse
        */
    }

    return *this;
}



///////////////////////////////////////////////////////////////////////////////
// compute the inverse of 4x4 Euclidean transformation matrix
//
// Euclidean transformation is translation, rotation, and reflection.
// With Euclidean transform, only the position and orientation of the object
// will be changed. Euclidean transform does not change the shape of an object
// (no scaling). Length and angle are reserved.
//
// Use inverseAffine() if the matrix has scale and shear transformation.
//
// M = [ R | T ]
//     [ --+-- ]    (R denotes 3x3 rotation/reflection matrix)
//     [ 0 | 1 ]    (T denotes 1x3 translation matrix)
//
// y = M*x  ->  y = R*x + T  ->  x = R^-1*(y - T)  ->  x = R^T*y - R^T*T
// (R is orthogonal,  R^-1 = R^T)
//
//  [ R | T ]-1    [ R^T | -R^T * T ]    (R denotes 3x3 rotation matrix)
//  [ --+-- ]   =  [ ----+--------- ]    (T denotes 1x3 translation)
//  [ 0 | 1 ]      [  0  |     1    ]    (R^T denotes R-transpose)
///////////////////////////////////////////////////////////////////////////////
//Matrix4& Matrix4::invertEuclidean()
//{

  //  return *this;
//}



///////////////////////////////////////////////////////////////////////////////
// compute the inverse of a 4x4 affine transformation matrix
//
// Affine transformations are generalizations of Euclidean transformations.
// Affine transformation includes translation, rotation, reflection, scaling,
// and shearing. Length and angle are NOT preserved.
// M = [ R | T ]
//     [ --+-- ]    (R denotes 3x3 rotation/scale/shear matrix)
//     [ 0 | 1 ]    (T denotes 1x3 translation matrix)
//
// y = M*x  ->  y = R*x + T  ->  x = R^-1*(y - T)  ->  x = R^-1*y - R^-1*T
//
//  [ R | T ]-1   [ R^-1 | -R^-1 * T ]
//  [ --+-- ]   = [ -----+---------- ]
//  [ 0 | 1 ]     [  0   +     1     ]
///////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////
// inverse matrix using matrix partitioning (blockwise inverse)
// It devides a 4x4 matrix into 4 of 2x2 matrices. It works in case of where
// det(A) != 0. If not, use the generic inverse method
// inverse formula.
// M = [ A | B ]    A, B, C, D are 2x2 matrix blocks
//     [ --+-- ]    det(M) = |A| * |D - ((C * A^-1) * B)|
//     [ C | D ]
//
// M^-1 = [ A' | B' ]   A' = A^-1 - (A^-1 * B) * C'
//        [ ---+--- ]   B' = (A^-1 * B) * -D'
//        [ C' | D' ]   C' = -D' * (C * A^-1)
//                      D' = (D - ((C * A^-1) * B))^-1
//
// NOTE: I wrap with () if it it used more than once.
//       The matrix is invertable even if det(A)=0, so must check det(A) before
//       calling this function, and use invertGeneric() instead.
///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
// compute the inverse of a general 4x4 matrix using Cramer's Rule
// If cannot find inverse, return indentity matrix
// M^-1 = adj(M) / det(M)
///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
// return determinant of 4x4 matrix
///////////////////////////////////////////////////////////////////////////////
float Matrix4::getDeterminant()
{
    return rows[0][0] * getCofactor(rows[1][1],rows[2][1],rows[3][1], rows[1][2],rows[2][2],rows[2][3], rows[3][1],rows[3][2],rows[3][3]) -
           rows[1][0] * getCofactor(rows[1][0],rows[1][2],rows[1][3], rows[2][0],rows[2][2],rows[2][3], rows[3][0],rows[3][2],rows[3][3]) +
           rows[2][0] * getCofactor(rows[1][0],rows[1][1],rows[1][3], rows[2][0],rows[2][1],rows[2][3], rows[3][0],rows[3][1],rows[3][3]) -
           rows[3][0] * getCofactor(rows[1][0],rows[1][1],rows[1][2], rows[2][0],rows[2][1],rows[2][2], rows[3][0],rows[3][1],rows[3][2]);
}

///////////////////////////////////////////////////////////////////////////////
// compute cofactor of 3x3 minor matrix without sign
// input params are 9 elements of the minor matrix
// NOTE: The caller must know its sign.
///////////////////////////////////////////////////////////////////////////////
float Matrix4::getCofactor(float m0, float m1, float m2,
                           float m3, float m4, float m5,
                           float m6, float m7, float m8)
{
    return m0 * (m4 * m8 - m5 * m7) -
           m1 * (m3 * m8 - m5 * m6) +
           m2 * (m3 * m7 - m4 * m6);
}



///////////////////////////////////////////////////////////////////////////////
// translate this matrix by (x, y, z)
///////////////////////////////////////////////////////////////////////////////
    // stataic transform matrix
Matrix4   Matrix4::Translate(const float3& v)           // creat translation matrix to translate by (x,y,z)
{
    Matrix4 m{ identity() };
    m.setCol(3, float4{ v.x, v.y, v.z, 1.f });
    return m;
}

Matrix4& Matrix4::translate(const float3 v)
{
    Matrix4 tr{ Translate(v) };
    (*this) = tr * (*this);
    return *this;
}

Matrix4& Matrix4::translate(float x, float y, float z)
{
    return translate({ x, y, z });
}



///////////////////////////////////////////////////////////////////////////////
// scale:  geneate and return scale matrix
///////////////////////////////////////////////////////////////////////////////
	Matrix4 Matrix4::Scale(float3 scale)  //
	{
		Matrix4 m{ identity() };
        m[0][0] = scale.x;
		m[1][1] = scale.y;
		m[2][2] = scale.z;
		return m;
	}

	Matrix4& Matrix4::scale(float s) { return scale(s, s, s); }                 // scale by uniform scale
	Matrix4& Matrix4::scale(float x, float y, float z)
	{
        *this = Scale({ x, y, z }) * *this;
		return *this;
	}


//static method: rotate angle(degree) along the given axix. Using Rodrigues formula 
Matrix4 Matrix4::Rotate(float angle, const float3& inaxis) 
{
    float3 axis = Normal(inaxis);
    float c = cosf(angle * DEG2RAD);    // cosine
    float s = sinf(angle * DEG2RAD);    // sine
    float c1 = 1.0f - c;                // 1 - c
    Matrix4 m{};
    m[0][0] = (axis.x * axis.x) * c1 + c;
    m[0][1] = (axis.y * axis.x) * c1 - (axis.z * s);
    m[0][2] = (axis.z * axis.x) * c1 + (axis.y * s);

    m[1][0] = (axis.x * axis.y) * c1 + (axis.z * s);
    m[1][1] = (axis.y * axis.y) * c1 + c;
    m[1][2] = (axis.z * axis.y) * c1 - (axis.x * s);

    m[2][0] = (axis.x * axis.z) * c1 - (axis.y * s);
    m[2][1] = (axis.y * axis.z) * c1 + (axis.x * s);
    m[2][2] = (axis.z * axis.z) * c1 + c;

    return m;
}
///////////////////////////////////////////////////////////////////////////////
// build a rotation matrix with given angle(degree) and rotation axis, then
// multiply it with this object
///////////////////////////////////////////////////////////////////////////////
Matrix4& Matrix4::rotate(float angle, const float3& axis)
{
    *this = Rotate(angle, axis) * *this;
    return *this;
}

Matrix4& Matrix4::rotateX(float angle)
{
    *this =  Rotate(angle, { 1.f, 0.f, 0.f }) * *this;
    return *this;
}

Matrix4& Matrix4::rotateY(float angle)
{
    *this = Rotate(angle, { 0.f, 1.f, 0.f }) * *this;
    return *this;

}

Matrix4& Matrix4::rotateZ(float angle)
{
    *this = Rotate(angle, { 0.f, 0.f, 1.f }) * *this;
    return *this;
}

