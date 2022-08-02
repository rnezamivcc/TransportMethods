///////////////////////////////////////////////////////////////////////////////
// Matrice.h
// =========
// NxN Matrix Math classes
// Reza Nezami: It is not row-first matrix

///////////////////////////////////////////////////////////////////////////////
#pragma once

#include <iostream>
#include <iomanip>
#include "vectors.h"
#include <assert.h>
///////////////////////////////////////////////////////////////////////////
// 2x2 matrix
///////////////////////////////////////////////////////////////////////////
class Matrix2
{
public:
    // constructors
    Matrix2();  // init with identity
	Matrix2(const float2D& row1, const float2D& row2) : rows{ row1, row2 } {}

	void setRow(size_t idx, const float2D row) {
		assert(idx < 2);
		rows[idx] = row;
	}

	float2D getRow(size_t id) const {
		assert(id < 2); 
		return rows[id];
	}

	void  setCol(size_t idx, const float2D col) {
		assert(idx < 2);
		rows[0][idx] = col[0];
		rows[1][idx] = col[1];

	}
	float2D getCol(size_t id) const {
		assert(id < 2); 
		return{ rows[0][id], rows[1][id] };
	}
	float getDeterminant() {
		return rows[0].x * rows[1].y - rows[0].y * rows[1].x;
	}

	static Matrix2 identity() { return Matrix2(	{ 1.f, 0.f }, 
												{ 0.f, 1.f }); }
	static Matrix2 zero()	{	return Matrix2({ 0.f, 0.f },
												{ 0.f, 0.f }); }
	Matrix2    transpose() {return Matrix2(rows[1], rows[0]);}

    Matrix2    invert();

    // operators
    Matrix2     operator+(const Matrix2& rhs) const;    // add rhs
    Matrix2     operator-(const Matrix2& rhs) const;    // subtract rhs
    Matrix2&    operator+=(const Matrix2& rhs);         // add rhs and update this object
    Matrix2&    operator-=(const Matrix2& rhs);         // subtract rhs and update this object
    float2D     operator*(const float2D& rhs) const;    // multiplication: v' = M * v
    Matrix2     operator*(const Matrix2& rhs) const;    // multiplication: M3 = M1 * M2
    Matrix2&    operator*=(const Matrix2& rhs);         // multiplication: M1' = M1 * M2
    bool        operator==(const Matrix2& rhs) const;   // exact compare, no epsilon
    bool        operator!=(const Matrix2& rhs) const;   // exact compare, no epsilon
    float       operator()(int idx1, int idx2) const;            // subscript operator v[0,1]
    float&      operator()(int idx1, int idx2);                  // subscript operator v[0,1]= 2.f

    friend Matrix2 operator-(const Matrix2& m);                     // unary operator (-)
    friend Matrix2 operator*(float scalar, const Matrix2& m);       // pre-multiplication
    friend float2D operator*(const float2D& vec, const Matrix2& m); // pre-multiplication
    friend std::ostream& operator<<(std::ostream& os, const Matrix2& m);

protected:

private:
	float2D rows[2];

};

///////////////////////////////////////////////////////////////////////////
// inline functions for Matrix2
///////////////////////////////////////////////////////////////////////////
inline Matrix2::Matrix2()
{
	   // initially identity matrix
	*this = identity();
}

inline Matrix2 Matrix2::operator+(const Matrix2& rhs) const
{
	return Matrix2(rows[0] + rhs.rows[0], rows[1] + rhs.rows[1]);
}

inline Matrix2 Matrix2::operator-(const Matrix2& rhs) const
{
	return Matrix2(rows[0] - rhs.rows[0], rows[1] - rhs.rows[1]);
}

inline Matrix2& Matrix2::operator+=(const Matrix2& rhs)
{
	rows[0] += rhs.rows[0];
	rows[1] += rhs.rows[1];
	return *this;
}

inline Matrix2& Matrix2::operator-=(const Matrix2& rhs)
{
	rows[0] -= rhs.rows[0];
	rows[1] -= rhs.rows[1];
	return *this;
}

inline float2D Matrix2::operator*(const float2D& rhs) const
{
	return float2D(Dot(rows[0],rhs), Dot(rows[1],rhs));
}

inline Matrix2 Matrix2::operator*(const Matrix2& rhs) const
{
	return Matrix2(	{ Dot(rows[0], rhs.getCol(0)) , Dot(rows[0], rhs.getCol(1)) },
					{ Dot(rows[1], rhs.getCol(0)) , Dot(rows[1], rhs.getCol(1)) });
}

inline Matrix2& Matrix2::operator*=(const Matrix2& rhs)
{
	*this = *this * rhs;
	return *this;
}

inline bool Matrix2::operator==(const Matrix2& rhs) const
{
	return (rows[0]==rhs.rows[0] && rows[1]==rhs.rows[1]);
}

inline bool Matrix2::operator!=(const Matrix2& rhs) const
{
	return !(*this == rhs);
}



inline float Matrix2::operator()(int idx1, int idx2) const
{
    return rows[idx1][idx2];
}



inline float& Matrix2::operator()(int idx1, int idx2)
{
    return rows[idx1][idx2];
}



inline Matrix2 operator-(const Matrix2& rhs)
{
	return Matrix2(-rhs.rows[0], -rhs.rows[1]);
}

inline Matrix2 operator*(float s, const Matrix2& rhs)
{
	return Matrix2(s*rhs.rows[0], s*rhs.rows[1]);
}

inline float2D operator*(const float2D& v, const Matrix2& rhs)
{
	return{ Dot(v, rhs.getCol(0)), Dot(v, rhs.getCol(1)) };
}

inline std::ostream& operator<<(std::ostream& os, const Matrix2& m)
{
	os << std::fixed << std::setprecision(5);
	os << "[" << std::setw(10) << m.rows[0][0] << " " << std::setw(10) << m.rows[0][1] << "]\n"
		<< "[" << std::setw(10) << m.rows[1][0] << " " << std::setw(10) << m.rows[1][1] << "]\n";
	os << std::resetiosflags(std::ios_base::fixed | std::ios_base::floatfield);
	return os;
}
// END OF MATRIX2 INLINE //////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////
// 3x3 matrix
///////////////////////////////////////////////////////////////////////////
class Matrix3
{
public:
    // constructors
    Matrix3();  // init with identity
	Matrix3(float3 row1, float3 row2, float3 row3)
		: rows{ row1, row2, row3 } {}

    void  setRow(size_t idx, float3 row)
	{
		assert(idx < 3);
		rows[idx] = row;
	}
	void  setCol(size_t idx, float3 v) 
	{
		assert(idx < 3);
		rows[0].set(idx, v[0]);
		rows[1].set(idx, v[1]);
		rows[2].set(idx, v[2]);
	}

	float3 getRow(size_t i) const 
	{
		assert(i < 3); return rows[i];
	}
	float3 getCol(size_t i) const 
	{
		assert(i < 3);
		return{ rows[0][i], rows[1][i], rows[2][i] };

	}
	float       getDeterminant();

	static Matrix3  identity() { return Matrix3({ 1.f, 0.f, 0.f }, 
												{ 0.f, 1.f, 0.f },
												{ 0.f, 0.f, 1.f }); 
	}
	static Matrix3  zero() { return Matrix3(	{ 0.f, 0.f, 0.f },
												{ 0.f, 0.f, 0.f },
												{ 0.f, 0.f, 0.f });
	}
	Matrix3    transpose() {return Matrix3(getCol(0), getCol(1), getCol(2)); }
    Matrix3    invert();

    // operators
    Matrix3     operator+(const Matrix3& rhs) const;    // add rhs
    Matrix3     operator-(const Matrix3& rhs) const;    // subtract rhs
    Matrix3&    operator+=(const Matrix3& rhs);         // add rhs and update this object
    Matrix3&    operator-=(const Matrix3& rhs);         // subtract rhs and update this object
	float3     operator*(const float3& rhs) const;    // multiplication: v' = M * v
    Matrix3     operator*(const Matrix3& rhs) const;    // multiplication: M3 = M1 * M2
    Matrix3&    operator*=(const Matrix3& rhs);         // multiplication: M1' = M1 * M2
    bool        operator==(const Matrix3& rhs) const;   // exact compare, no epsilon
    bool        operator!=(const Matrix3& rhs) const;   // exact compare, no epsilon
    float       operator[](int index) const;            // subscript operator v[0], v[1]
    float&      operator[](int index);                  // subscript operator v[0], v[1]

    friend Matrix3 operator-(const Matrix3& m);                     // unary operator (-)
    friend Matrix3 operator*(float scalar, const Matrix3& m);       // pre-multiplication
    friend float3 operator*(const float3& vec, const Matrix3& m); // pre-multiplication
    friend std::ostream& operator<<(std::ostream& os, const Matrix3& m);

protected:

private:
	float3 rows[3];
};

///////////////////////////////////////////////////////////////////////////
// inline functions for Matrix3
///////////////////////////////////////////////////////////////////////////

inline Matrix3::Matrix3()
{
	*this = identity();
}
inline Matrix3 Matrix3::operator+(const Matrix3& rhs) const
{
	return Matrix3(rows[0] + rhs.rows[0], rows[1] + rhs.rows[1], rows[2] + rhs.rows[2]);
}

inline Matrix3 Matrix3::operator-(const Matrix3& rhs) const
{
	return Matrix3(rows[0] - rhs.rows[0], rows[1] - rhs.rows[1], rows[2] - rhs.rows[2]);
}

inline Matrix3& Matrix3::operator+=(const Matrix3& rhs)
{
	rows[0] += rhs.rows[0];
	rows[1] += rhs.rows[1];
	rows[2] += rhs.rows[2];
	return *this;
}

inline Matrix3& Matrix3::operator-=(const Matrix3& rhs)
{
	rows[0] -= rhs.rows[0];
	rows[1] -= rhs.rows[1];
	rows[2] -= rhs.rows[2];
	return *this;
}

inline float3 Matrix3::operator*(const float3& rhs) const
{
	return float3(Dot(rows[0], rhs), Dot(rows[1], rhs), Dot(rows[2], rhs));
}

inline Matrix3 Matrix3::operator*(const Matrix3& rhs) const
{
	return Matrix3( { Dot(rows[0], rhs.getCol(0)) , Dot(rows[0], rhs.getCol(1)), Dot(rows[0], rhs.getCol(2)) },
					{ Dot(rows[1], rhs.getCol(0)) , Dot(rows[1], rhs.getCol(1)), Dot(rows[1], rhs.getCol(2)) },
					{ Dot(rows[2], rhs.getCol(0)) , Dot(rows[2], rhs.getCol(1)), Dot(rows[2], rhs.getCol(2)) } );
}

inline Matrix3& Matrix3::operator*=(const Matrix3& rhs)
{
    *this = *this * rhs;
    return *this;
}

inline bool Matrix3::operator==(const Matrix3& rhs) const
{
	return (rows[0] == rhs.rows[0] && rows[1] == rhs.rows[1] && rows[2] == rhs.rows[2]);
}

inline bool Matrix3::operator!=(const Matrix3& rhs) const
{
	return !(*this == rhs);
}

//////////// friends //////////////////////
inline Matrix3 operator-(const Matrix3& rhs)
{
    return Matrix3(-rhs.rows[0], -rhs.rows[1], -rhs.rows[2]);
}

inline Matrix3 operator*(float s, const Matrix3& rhs)
{
	return Matrix3(s*rhs.rows[0], s*rhs.rows[1], s*rhs.rows[2]);
}

inline float3 operator*(const float3& v, const Matrix3& rhs)
{
	return{ Dot(v, rhs.getCol(0)), Dot(v, rhs.getCol(1)), Dot(v, rhs.getCol(2)) };
}

inline std::ostream& operator<<(std::ostream& os, const Matrix3& m)
{
	os << std::fixed << std::setprecision(5);
	os << "[" << std::setw(10) << m.rows[0][0] << " " << std::setw(10) << m.rows[0][1] << " " << std::setw(10) << m.rows[0][2] << "]\n"
	   << "[" << std::setw(10) << m.rows[1][0] << " " << std::setw(10) << m.rows[1][1] << " " << std::setw(10) << m.rows[1][2] << "]\n"
	   << "[" << std::setw(10) << m.rows[2][0] << " " << std::setw(10) << m.rows[2][1] << " " << std::setw(10) << m.rows[2][2] << "]\n";
	os << std::resetiosflags(std::ios_base::fixed | std::ios_base::floatfield);
	return os;
}
// END OF MATRIX3 INLINE //////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////
// 4x4 matrix
///////////////////////////////////////////////////////////////////////////
class Matrix4
{
public:
    // constructors
	Matrix4();  // init with identity matrix\
	{ *this = identity(); }
	Matrix4(float4 row1, float4 row2, float4 row3, float4 row4)
		: rows{ row1, row2, row3, row4 } {}
		
	void  setRow(size_t idx, float4 row);
	void  setCol(size_t idx, float4 v);	
	float4 getRow(size_t i) const;
	float4 getCol(size_t i) const;
	float* get(); //{ return &rows[0][0]; }
	
	static Matrix4  identity() {
		return Matrix4({ 1.f, 0.f, 0.f, 0.f },
						{ 0.f, 1.f, 0.f, 0.f },
						{ 0.f, 0.f, 1.f, 0.f },
						{ 0.f, 0.f, 0.f, 1.f });
	}
	
	Matrix4    getTranspose() { return Matrix4(getCol(0), getCol(1), getCol(2), getCol(3)); }
	Matrix4    transpose() { *this = getTranspose(); return *this; }	
	
		
	Matrix4&    invert();
	float getDeterminant();
	 // inverse of projective matrix using partitioning
   // Matrix4&    invertGeneral();                        // inverse of generic matrix

    // transform matrix
	static Matrix4    Translate(const float3& v);
    Matrix4&    translate(float x, float y, float z);   // translation by (x,y,z)
    Matrix4&    translate(const float3 v);            //

	static Matrix4    Scale(float3 scale);                    // geneate and return uniform scale matrix
    Matrix4&    scale(float scale);                     // uniform scale
    Matrix4&    scale(float sx, float sy, float sz);    // scale by (sx, sy, sz) on each axis
 
	static Matrix4	Rotate(float angle, const float3& inaxis);
    Matrix4&    rotate(float angle, const float3& axis); // rotate angle(degree) along the given axix
    Matrix4&    rotate(float angle, float x, float y, float z);
    Matrix4&    rotateX(float angle);                   // rotate on X-axis with degree
    Matrix4&    rotateY(float angle);                   // rotate on Y-axis with degree
    Matrix4&    rotateZ(float angle);                   // rotate on Z-axis with degree
 
    // operators
    Matrix4     operator+(const Matrix4& rhs) const;    // add rhs
    Matrix4     operator-(const Matrix4& rhs) const;    // subtract rhs
    Matrix4&    operator+=(const Matrix4& rhs);         // add rhs and update this object
    Matrix4&    operator-=(const Matrix4& rhs);         // subtract rhs and update this object
    float4     operator*(float4 rhs) const;    // multiplication: v' = M * v
    float3     operator*(float3 rhs) const;    // multiplication: v' = M * v
    Matrix4     operator*(const Matrix4& rhs) const;    // multiplication: M3 = M1 * M2
    Matrix4&    operator*=(const Matrix4& rhs);         // multiplication: M1' = M1 * M2
    bool        operator==(const Matrix4& rhs) const;   // exact compare, no epsilon
    bool        operator!=(const Matrix4& rhs) const;   // exact compare, no epsilon
    float4       operator[](int index) const;            // subscript operator v[0], v[1]
    float4&      operator[](int index);                  // subscript operator v[0], v[1]

    friend Matrix4 operator-(const Matrix4& m);                     // unary operator (-)
    friend Matrix4 operator*(float scalar, const Matrix4& m);       // pre-multiplication
    friend float3 operator*(const float3& vec, const Matrix4& m); // pre-multiplication
    friend float4 operator*(const float4& vec, const Matrix4& m); // pre-multiplication
    friend std::ostream& operator<<(std::ostream& os, const Matrix4& m);

protected:

private:
    float       getCofactor(float m0, float m1, float m2,
                            float m3, float m4, float m5,
                            float m6, float m7, float m8);

	float4 rows[4];
};


///////////////////////////////////////////////////////////////////////////
// inline functions for Matrix4
///////////////////////////////////////////////////////////////////////////
inline Matrix4::Matrix4()
{
    // initially identity matrix
	rows[0] = rows[1] = rows[2] = rows[3] = float4{};
	rows[0].x = rows[1].y = rows[2].z = rows[3].w = 1.f;
}

inline Matrix4 Matrix4::operator+(const Matrix4& rhs) const
{
	return Matrix4(rows[0] + rhs.rows[0], rows[1] + rhs.rows[1], rows[2] + rhs.rows[2], rows[3] + rhs.rows[3]);
}

inline Matrix4 Matrix4::operator-(const Matrix4& rhs) const
{
	return Matrix4(rows[0] - rhs.rows[0], rows[1] - rhs.rows[1], rows[2] - rhs.rows[2], rows[3] - rhs.rows[3]);
}

inline Matrix4& Matrix4::operator+=(const Matrix4& rhs)
{
	rows[0] += rhs.rows[0];
	rows[1] += rhs.rows[1];
	rows[2] += rhs.rows[2];
	rows[3] += rhs.rows[3];
	return *this;
}

inline Matrix4& Matrix4::operator-=(const Matrix4& rhs)
{
	rows[0] -= rhs.rows[0];
	rows[1] -= rhs.rows[1];
	rows[2] -= rhs.rows[2];
	rows[3] -= rhs.rows[3];
	return *this;
}

inline float3 Matrix4::operator*(float3 rhs) const
{
	return *this * float4(rhs);
}


inline float4 Matrix4::operator*(float4 rhs) const
{
	return float4(Dot(rows[0], rhs), Dot(rows[1], rhs), Dot(rows[2], rhs), Dot(rows[3], rhs));
}

inline Matrix4 Matrix4::operator*(const Matrix4& rhs) const
{
	return Matrix4(	{ Dot(rows[0], rhs.getCol(0)), Dot(rows[0], rhs.getCol(1)), Dot(rows[0], rhs.getCol(2)), Dot(rows[0], rhs.getCol(3)) },
					{ Dot(rows[1], rhs.getCol(0)), Dot(rows[1], rhs.getCol(1)), Dot(rows[1], rhs.getCol(2)), Dot(rows[1], rhs.getCol(3)) },
					{ Dot(rows[2], rhs.getCol(0)), Dot(rows[2], rhs.getCol(1)), Dot(rows[2], rhs.getCol(2)), Dot(rows[2], rhs.getCol(3)) },
					{ Dot(rows[3], rhs.getCol(0)), Dot(rows[3], rhs.getCol(1)), Dot(rows[3], rhs.getCol(2)), Dot(rows[3], rhs.getCol(3)) });
}

inline Matrix4& Matrix4::operator*=(const Matrix4& rhs)
{
	*this = *this * rhs;
	return *this;
}

inline bool Matrix4::operator==(const Matrix4& rhs) const
{
	return (rows[0] == rhs.rows[0] && rows[1] == rhs.rows[1] && rows[2] == rhs.rows[2] && rows[3] == rhs.rows[3]);
}

inline bool Matrix4::operator!=(const Matrix4& rhs) const
{
	return !(*this == rhs);
}

inline Matrix4 operator-(const Matrix4& rhs)
{
	return -1.f * rhs;
}

inline Matrix4 operator*(float s, const Matrix4& rhs)
{
	return Matrix4(s*rhs.rows[0], s*rhs.rows[1], s*rhs.rows[2], s*rhs.rows[3]);
}

inline float4 operator*(const float4& v, const Matrix4& rhs)
{
	return{ Dot(v, rhs.getCol(0)), Dot(v, rhs.getCol(1)), Dot(v, rhs.getCol(2)), Dot(v, rhs.getCol(3)) };
}
inline float4 Matrix4::operator[](int index) const
{
    return rows[index];
}

inline float4& Matrix4::operator[](int index)
{
    return rows[index];
}

inline std::ostream& operator<<(std::ostream& os, const Matrix4& m)
{
    os << std::fixed << std::setprecision(5);
    os << "[" << std::setw(10) << m.rows[0][0] << " " << std::setw(10) << m.rows[0][1] << " " << std::setw(10) << m.rows[0][2] <<  " " << std::setw(10) << m.rows[0][3] << "]\n"
       << "[" << std::setw(10) << m.rows[1][0] << " " << std::setw(10) << m.rows[1][1] << " " << std::setw(10) << m.rows[1][2] <<  " " << std::setw(10) << m.rows[1][3] << "]\n"
       << "[" << std::setw(10) << m.rows[2][0] << " " << std::setw(10) << m.rows[2][1] << " " << std::setw(10) << m.rows[2][2] <<  " " << std::setw(10) << m.rows[2][3] << "]\n"
       << "[" << std::setw(10) << m.rows[3][0] << " " << std::setw(10) << m.rows[3][1] << " " << std::setw(10) << m.rows[3][2] <<  " " << std::setw(10) << m.rows[3][3] << "]\n";
    os << std::resetiosflags(std::ios_base::fixed | std::ios_base::floatfield);
    return os;
}


inline void Matrix4::setRow(size_t idx, float4 row)
{
	assert(idx <= 3);
	rows[idx] = row;
}
inline void  Matrix4::setCol(size_t idx, float4 v)
{
	rows[0][idx] = v[0];
	rows[1][idx] = v[1];
	rows[2][idx] = v[2];
	rows[3][idx] = v[3];
}


inline float4 Matrix4::getRow(size_t i) const
{
	assert(i <= 3); 
	return rows[i];
}

inline float4 Matrix4::getCol(size_t i) const
{
	assert(i <= 3);
	return{ rows[0][i], rows[1][i], rows[2][i], rows[3][i] };
}

inline float* Matrix4::get()
{
	return &this->transpose()[0][0];
}

// END OF MATRIX4 INLINE //////////////////////////////////////////////////////
