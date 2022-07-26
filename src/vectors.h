#pragma once


#include <cmath>
#include <cfloat>
#include <climits>
#include <tuple>
#include <assert.h>

/*The Vector2d class is an object consisting of simply an x and
y value. Certain operators are overloaded to make it easier
for vector math to be performed.*/
template<typename T>
class Vector2d {
public:
	/*The x and y values are public to give easier access for
	outside funtions. Accessors and mutators are not really
	necessary*/
	T x;
	T y;

	//Constructor assigns the inputs to x and y.
	Vector2d() : x(T(0)), y(T(0)) {}
	Vector2d(T a, T b) : x(a), y(b) {}

	void reset(void) { 
		x = T();
		y = T(); 
	}

	T operator[] (size_t i) const {
		assert(i < 2);
		T ret{ x };
		if (i == 1)
			ret = { y };
		return ret;
	}

	T& operator[] (size_t i) {
		assert(i < 2);
		if (i == 0)
			return x;
		return y;
	}
	/*The following operators simply return Vector3ds that
	have operations performed on the relative (x, y) values*/
	Vector2d& operator+=(const Vector2d& v) { x += v.x; y += v.y; return *this; }
	Vector2d& operator-=(const Vector2d& v) { x -= v.x; y -= v.y; return *this; }
	Vector2d& operator*=(const Vector2d& v) { x *= v.x; y *= v.y; return *this; }
	Vector2d& operator/=(const Vector2d& v) { x /= v.x; y /= v.y; return *this; }

	Vector2d& operator+=(const float s) { x += s; y += s;  return *this; }
	Vector2d& operator-=(const float s) { x -= s; y -= s;  return *this; }
	Vector2d& operator*=(const float s) { x *= s; y *= s;  return *this; }
	Vector2d& operator/=(const float s) { x /= s; y /= s;  return *this; }

	//Check if the Vectors have the same values (uses pairwise comparison of 
	// 'std::tuple' on the x, y values of L and R.)
	friend bool operator==(const Vector2d& L, const Vector2d& R) {
		return std::tie(L.x, L.y) == std::tie(R.x, R.y);
	}
	friend bool operator!=(const Vector2d& L, const Vector2d& R) { return !(L == R); }

	void set(T a, T b) { x = a; y = b; }
	/*Check which Vectors are closer or further from the origin.*/
	friend bool operator>(const Vector2d& L, const Vector2d& R) { return LengthSq(L) < LengthSq(R); }
	friend bool operator>=(const Vector2d& L, const Vector2d& R) { return !(L > R); }
	friend bool operator<(const Vector2d& L, const Vector2d& R) { return R < L; }
	friend bool operator<=(const Vector2d& L, const Vector2d& R) { return !(R < L); }

	//Negate both the x and y values.
	Vector2d operator-() const { return Vector2d(-x, -y); }

	//Apply scalar operations.
	Vector2d operator*(T s) { Vector2d tmp(*this); tmp.x *= s; tmp.y *= s;  return tmp; }
	Vector2d operator/(T s) { Vector2d tmp(*this); tmp.x /= s; tmp.y /= s;  return tmp; }

	//Returns the length of the vector from the origin.
	float Length() const { return sqrt(x*x + y*y); }
	float LengthSq()const { return x*x + y*y; }

};

using float2D = Vector2d<float>;
template<class T> Vector2d<T> operator*(const T& s, const Vector2d<T>& v) { return Vector2d<T>(v) *= s; }
template<class T> Vector2d<T> operator*(const Vector2d<T>& v, const T& s) { return Vector2d<T>(v) *= s; }

template<class T> Vector2d<T>  operator-(const Vector2d<T>& v1, const Vector2d<T>& v2) { return Vector2d<T>(v1.x-v2.x, v1.y-v2.y); }
template<class T> Vector2d<T>  operator+(const Vector2d<T>& v1, const Vector2d<T>& v2) { return Vector2d<T>(v1.x+v2.x, v1.y+v2.y); }

//Product functions
template<class T> T Dot(const Vector2d<T>& a, const Vector2d<T>& b) { return  ((a.x * b.x) + (a.y * b.y)); }
template<class T> T Cross(const Vector2d<T>& a, const Vector2d<T>& b) { return ((a.x * b.y) - (a.y * b.x)); }

//Return the unit vector of the input
template<class T> Vector2d<T> Normal(const Vector2d<T>& a) { double mag = a.Length(); return Vector2d<T>(a.x / mag, a.y / mag); }

//Return a vector perpendicular to the left.
template<class T> Vector2d<T> Perpendicular(const Vector2d<T>& a) { return Vector2d<T>(a.y, -a.x); }
//Return true if two line segments intersect.
template<class T> 
bool Intersect(const Vector2d<T>&aa, const Vector2d<T>& ab, const Vector2d<T>& ba, const Vector2d<T>& bb)
{
	Vector2d<T> p = aa;
	Vector2d<T> r = ab - aa;
	Vector2d<T> q = ba;
	Vector2d<T> s = bb - ba;

	double t = Cross((q - p), s) / Cross(r, s);
	double u = Cross((q - p), r) / Cross(r, s);

	return (0.0 <= t && t <= 1.0) &&
		(0.0 <= u && u <= 1.0);
}

//Return the point where two lines intersect.
template<class T>
Vector2d<T> GetIntersect(const Vector2d<T>&aa, const Vector2d<T>& ab, const Vector2d<T>& ba, const Vector2d<T>& bb)
{
	double pX = (aa.x*ab.y - aa.y*ab.x)*(ba.x - bb.x) -
		(ba.x*bb.y - ba.y*bb.x)*(aa.x - ab.x);
	double pY = (aa.x*ab.y - aa.y*ab.x)*(ba.y - bb.y) -
		(ba.x*bb.y - ba.y*bb.x)*(aa.y - ab.y);
	double denominator = (aa.x - ab.x)*(ba.y - bb.y) -
		(aa.y - ab.y)*(ba.x - bb.x);

	return Vector2d(pX / denominator, pY / denominator);
}

/////////////////// vector 3D ///////////////////////////
template<typename T>
class Vector3 {
public:
    /*The x and y values are public to give easier access for
    outside funtions. Accessors and mutators are not really
    necessary*/
    T x;
    T y;
    T z;

    //Constructor assigns the inputs to x and y.
    Vector3<T>() : x(T()), y(T()), z(T()) {}
    Vector3<T>(T a, T b, T c = T()) : x(a), y(b), z(c) {}
    Vector3<T>(T* p) : x(p[0]), y(p[1]), z(p[2]) {}
  //  Vector3<T>(const Vector4d<T>& p4);

    void reset(void) {
        x = T();
        y = T();
        z = T();
    }
    T operator[] (size_t i) const {
        assert(i < 3);
        T ret{ x };
        if (i == 1)
            ret = { y };
        else
            ret = { z };
        return ret;
    }

    T& operator[] (size_t i) {
        assert(i < 2);
        if (i == 0)
            return x;
        else if (i == 1)
            return y;
        return z;
    }

    void set(size_t i, T val) {
        assert(i < 3);
        (i == 0) ? x = val : (i == 1) ? y = val : z = val;
    }
    /*The following operators simply return Vector3ds that
    have operations performed on the relative (x, y) values*/
    Vector3& operator+=(const Vector3& v) { x += v.x; y += v.y; z += v.z; return *this; }
    Vector3& operator-=(const Vector3& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
    Vector3& operator*=(const Vector3& v) { x *= v.x; y *= v.y; z *= v.z; return *this; }
    Vector3& operator/=(const Vector3& v) { x /= v.x; y /= v.y; z /= v.z; return *this; }

    Vector3& operator+=(const float s) { x += s; y += s; z += s; return *this; }
    Vector3& operator-=(const float s) { x -= s; y -= s; z -= s; return *this; }
    Vector3& operator*=(const float s) { x *= s; y *= s; z *= s; return *this; }
    Vector3& operator/=(const float s) { x /= s; y /= s; z /= s; return *this; }

    //Check if the Vectors have the same values (uses pairwise comparison of 
    // 'std::tuple' on the x, y values of L and R.)
    friend bool operator==(const Vector3& L, const Vector3& R) {
        return std::tie(L.x, L.y, L.z) == std::tie(R.x, R.y, R.z);
    }
    friend bool operator!=(const Vector3& L, const Vector3& R) { return !(L == R); }

    void set(T a, T b, T c) { x = a; y = b; z = c; }
    /*Check which Vectors are closer or further from the origin.*/
    friend bool operator>(const Vector3& L, const Vector3& R) { return LengthSq(L) < LengthSq(R); }
    friend bool operator>=(const Vector3& L, const Vector3& R) { return !(L > R); }
    friend bool operator<(const Vector3& L, const Vector3& R) { return R < L; }
    friend bool operator<=(const Vector3& L, const Vector3& R) { return !(R < L); }

    //Negate both the x and y values.
    Vector3 operator-() const { return Vector3(-x, -y, -z); }

    //Apply scalar operations.
    Vector3 operator*(T s) { Vector3 tmp(*this); tmp.x *= s; tmp.y *= s; tmp.z *= s;  return tmp; }
    Vector3 operator/(T s) { Vector3 tmp(*this); tmp.x /= s; tmp.y /= s; tmp.z /= s;  return tmp; }

    //Returns the length of the vector from the origin.
    double Length() const { return sqrt(x * x + y * y + z * z); }
    double LengthSq()const { return x * x + y * y + z * z; }

};
using float3D = Vector3<float>;
template<class T> Vector3<T> operator*(const T& s, const Vector3<T>& v) { return Vector3<T>(v) *= s; }
template<class T> Vector3<T> operator*(const Vector3<T>& v, const T& s) { return Vector3<T>(v) *= s; }

template<class T> Vector3<T>  operator-(const Vector3<T>& v1, const Vector3<T>& v2) { return Vector3<T>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z); }
template<class T> Vector3<T>  operator+(const Vector3<T>& v1, const Vector3<T>& v2) { return Vector3<T>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z); }

//Product functions
template<class T> T Dot(const Vector3<T>& a, const Vector3<T>& b) { return  ((a.x * b.x) + (a.y * b.y) + (a.z * b.z)); }
template<class T>  Vector3<T>  Cross(const Vector3<T>& a, const Vector3<T>& b)
{
    Vector3<T> ret{ a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
    return ret;
}

//Return the unit vector of the input
template<class T> Vector3<T> Normal(const Vector3<T>& a) { auto mag = a.Length(); return Vector3<T>(a.x / (T)mag, a.y / (T)mag, a.z / (T)mag); }
template<class T>
float SquarDist(const Vector3<T>& a, const Vector3<T>& b)
{
    Vector3<T> dist{ a - b };
    return dist.LengthSq();
}
//Return a vector perpendicular to the left.
//template<class T> Vector3d<T> Perpendicular(const Vector3d<T>& a) { return Vector3d<T>(a.y, -a.x); }
//Return true if two line segments intersect.
template<class T>
bool Intersect(const Vector3<T>& aa, const Vector3<T>& ab, const Vector3<T>& ba, const Vector3<T>& bb)
{
    Vector3<T> p = aa;
    Vector3<T> r = ab - aa;
    Vector3<T> q = ba;
    Vector3<T> s = bb - ba;

    T t = Cross((q - p), s) / Cross(r, s);
    T u = Cross((q - p), r) / Cross(r, s);

    return (0.0 <= t && t <= 1.0) && (0.0 <= u && u <= 1.0);
}

//Return the point where two lines intersect.
template<class T>
Vector3<T> GetIntersect(const Vector3<T>& aa, const Vector3<T>& ab, const Vector3<T>& ba, const Vector3<T>& bb)
{
    T pX = (aa.x * ab.y - aa.y * ab.x) * (ba.x - bb.x) -
        (ba.x * bb.y - ba.y * bb.x) * (aa.x - ab.x);
    T pY = (aa.x * ab.y - aa.y * ab.x) * (ba.y - bb.y) -
        (ba.x * bb.y - ba.y * bb.x) * (aa.y - ab.y);
    T denominator = (aa.x - ab.x) * (ba.y - bb.y) -
        (aa.y - ab.y) * (ba.x - bb.x);

    return Vector3d(pX / denominator, pY / denominator);
}


///////////////////////////////////////////////////////////////////////////////
// 4D vector
///////////////////////////////////////////////////////////////////////////////
struct Vector4
{
    float x;
    float y;
    float z;
    float w;

    // ctors
    Vector4() : x(0), y(0), z(0), w(0) {};
    Vector4(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {};

    // utils functions
    void        set(float x, float y, float z, float w);
    float       length() const;                         //
    float       distance(const Vector4& vec) const;     // distance between two vectors
    Vector4&    normalize();                            //
    float       dot(const Vector4& vec) const;          // dot product
    bool        equal(const Vector4& vec, float e) const; // compare with epsilon

    // operators
    Vector4     operator-() const;                      // unary operator (negate)
    Vector4     operator+(const Vector4& rhs) const;    // add rhs
    Vector4     operator-(const Vector4& rhs) const;    // subtract rhs
    Vector4&    operator+=(const Vector4& rhs);         // add rhs and update this object
    Vector4&    operator-=(const Vector4& rhs);         // subtract rhs and update this object
    Vector4     operator*(const float scale) const;     // scale
    Vector4     operator*(const Vector4& rhs) const;    // multiply each element
    Vector4&    operator*=(const float scale);          // scale and update this object
    Vector4&    operator*=(const Vector4& rhs);         // multiply each element and update this object
    Vector4     operator/(const float scale) const;     // inverse scale
    Vector4&    operator/=(const float scale);          // scale and update this object
    bool        operator==(const Vector4& rhs) const;   // exact compare, no epsilon
    bool        operator!=(const Vector4& rhs) const;   // exact compare, no epsilon
    bool        operator<(const Vector4& rhs) const;    // comparison for sort
    float       operator[](int index) const;            // subscript operator v[0], v[1]
    float&      operator[](int index);                  // subscript operator v[0], v[1]

    friend Vector4 operator*(const float a, const Vector4 vec);
    friend std::ostream& operator<<(std::ostream& os, const Vector4& vec);
};



// fast math routines from Doom3 SDK
inline float invSqrt(float x)
{
    float xhalf = 0.5f * x;
    int i = *(int*)&x;          // get bits for floating value
    i = 0x5f3759df - (i>>1);    // gives initial guess
    x = *(float*)&i;            // convert bits back to float
    x = x * (1.5f - xhalf*x*x); // Newton step
    return x;
}



///////////////////////////////////////////////////////////////////////////////
// inline functions for Vector2
///////////////////////////////////////////////////////////////////////////////
/*inline Vector2 Vector2::operator-() const {
    return Vector2(-x, -y);
}

inline Vector2 Vector2::operator+(const Vector2& rhs) const {
    return Vector2(x+rhs.x, y+rhs.y);
}

inline Vector2 Vector2::operator-(const Vector2& rhs) const {
    return Vector2(x-rhs.x, y-rhs.y);
}

inline Vector2& Vector2::operator+=(const Vector2& rhs) {
    x += rhs.x; y += rhs.y; return *this;
}

inline Vector2& Vector2::operator-=(const Vector2& rhs) {
    x -= rhs.x; y -= rhs.y; return *this;
}

inline Vector2 Vector2::operator*(const float a) const {
    return Vector2(x*a, y*a);
}

inline Vector2 Vector2::operator*(const Vector2& rhs) const {
    return Vector2(x*rhs.x, y*rhs.y);
}

inline Vector2& Vector2::operator*=(const float a) {
    x *= a; y *= a; return *this;
}

inline Vector2& Vector2::operator*=(const Vector2& rhs) {
    x *= rhs.x; y *= rhs.y; return *this;
}

inline Vector2 Vector2::operator/(const float a) const {
    return Vector2(x/a, y/a);
}

inline Vector2& Vector2::operator/=(const float a) {
    x /= a; y /= a; return *this;
}

inline bool Vector2::operator==(const Vector2& rhs) const {
    return (x == rhs.x) && (y == rhs.y);
}

inline bool Vector2::operator!=(const Vector2& rhs) const {
    return (x != rhs.x) || (y != rhs.y);
}

inline bool Vector2::operator<(const Vector2& rhs) const {
    if(x < rhs.x) return true;
    if(x > rhs.x) return false;
    if(y < rhs.y) return true;
    if(y > rhs.y) return false;
    return false;
}

inline float Vector2::operator[](int index) const {
    return (&x)[index];
}

inline float& Vector2::operator[](int index) {
    return (&x)[index];
}

inline void Vector2::set(float x, float y) {
    this->x = x; this->y = y;
}

inline float Vector2::length() const {
    return sqrtf(x*x + y*y);
}

inline float Vector2::distance(const Vector2& vec) const {
    return sqrtf((vec.x-x)*(vec.x-x) + (vec.y-y)*(vec.y-y));
}

inline Vector2& Vector2::normalize() {
    //@@const float EPSILON = 0.000001f;
    float xxyy = x*x + y*y;
    //@@if(xxyy < EPSILON)
    //@@    return *this;

    //float invLength = invSqrt(xxyy);
    float invLength = 1.0f / sqrtf(xxyy);
    x *= invLength;
    y *= invLength;
    return *this;
}

inline float Vector2::dot(const Vector2& rhs) const {
    return (x*rhs.x + y*rhs.y);
}

inline bool Vector2::equal(const Vector2& rhs, float epsilon) const {
    return fabs(x - rhs.x) < epsilon && fabs(y - rhs.y) < epsilon;
}

inline Vector2 operator*(const float a, const Vector2 vec) {
    return Vector2(a*vec.x, a*vec.y);
}

inline std::ostream& operator<<(std::ostream& os, const Vector2& vec) {
    os << "(" << vec.x << ", " << vec.y << ")";
    return os;
}
// END OF VECTOR2 /////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////
// inline functions for Vector3
///////////////////////////////////////////////////////////////////////////////
inline Vector3 Vector3::operator-() const {
    return Vector3(-x, -y, -z);
}

inline Vector3 Vector3::operator+(const Vector3& rhs) const {
    return Vector3(x+rhs.x, y+rhs.y, z+rhs.z);
}

inline Vector3 Vector3::operator-(const Vector3& rhs) const {
    return Vector3(x-rhs.x, y-rhs.y, z-rhs.z);
}

inline Vector3& Vector3::operator+=(const Vector3& rhs) {
    x += rhs.x; y += rhs.y; z += rhs.z; return *this;
}

inline Vector3& Vector3::operator-=(const Vector3& rhs) {
    x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this;
}

inline Vector3 Vector3::operator*(const float a) const {
    return Vector3(x*a, y*a, z*a);
}

inline Vector3 Vector3::operator*(const Vector3& rhs) const {
    return Vector3(x*rhs.x, y*rhs.y, z*rhs.z);
}

inline Vector3& Vector3::operator*=(const float a) {
    x *= a; y *= a; z *= a; return *this;
}

inline Vector3& Vector3::operator*=(const Vector3& rhs) {
    x *= rhs.x; y *= rhs.y; z *= rhs.z; return *this;
}

inline Vector3 Vector3::operator/(const float a) const {
    return Vector3(x/a, y/a, z/a);
}

inline Vector3& Vector3::operator/=(const float a) {
    x /= a; y /= a; z /= a; return *this;
}

inline bool Vector3::operator==(const Vector3& rhs) const {
    return (x == rhs.x) && (y == rhs.y) && (z == rhs.z);
}

inline bool Vector3::operator!=(const Vector3& rhs) const {
    return (x != rhs.x) || (y != rhs.y) || (z != rhs.z);
}

inline bool Vector3::operator<(const Vector3& rhs) const {
    if(x < rhs.x) return true;
    if(x > rhs.x) return false;
    if(y < rhs.y) return true;
    if(y > rhs.y) return false;
    if(z < rhs.z) return true;
    if(z > rhs.z) return false;
    return false;
}

inline float Vector3::operator[](int index) const {
    return (&x)[index];
}

inline float& Vector3::operator[](int index) {
    return (&x)[index];
}

inline void Vector3::set(float x, float y, float z) {
    this->x = x; this->y = y; this->z = z;
}

inline float Vector3::length() const {
    return sqrtf(x*x + y*y + z*z);
}

inline float Vector3::distance(const Vector3& vec) const {
    return sqrtf((vec.x-x)*(vec.x-x) + (vec.y-y)*(vec.y-y) + (vec.z-z)*(vec.z-z));
}

inline Vector3& Vector3::normalize() {
    //@@const float EPSILON = 0.000001f;
    float xxyyzz = x*x + y*y + z*z;
    //@@if(xxyyzz < EPSILON)
    //@@    return *this; // do nothing if it is ~zero vector

    //float invLength = invSqrt(xxyyzz);
    float invLength = 1.0f / sqrtf(xxyyzz);
    x *= invLength;
    y *= invLength;
    z *= invLength;
    return *this;
}

inline float Vector3::dot(const Vector3& rhs) const {
    return (x*rhs.x + y*rhs.y + z*rhs.z);
}

inline Vector3 Vector3::cross(const Vector3& rhs) const {
    return Vector3(y*rhs.z - z*rhs.y, z*rhs.x - x*rhs.z, x*rhs.y - y*rhs.x);
}

inline bool Vector3::equal(const Vector3& rhs, float epsilon) const {
    return fabs(x - rhs.x) < epsilon && fabs(y - rhs.y) < epsilon && fabs(z - rhs.z) < epsilon;
}

inline Vector3 operator*(const float a, const Vector3 vec) {
    return Vector3(a*vec.x, a*vec.y, a*vec.z);
}

inline std::ostream& operator<<(std::ostream& os, const Vector3& vec) {
    os << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")";
    return os;
}
// END OF VECTOR3 /////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
// inline functions for Vector4
///////////////////////////////////////////////////////////////////////////////
inline Vector4 Vector4::operator-() const {
    return Vector4(-x, -y, -z, -w);
}

inline Vector4 Vector4::operator+(const Vector4& rhs) const {
    return Vector4(x+rhs.x, y+rhs.y, z+rhs.z, w+rhs.w);
}

inline Vector4 Vector4::operator-(const Vector4& rhs) const {
    return Vector4(x-rhs.x, y-rhs.y, z-rhs.z, w-rhs.w);
}

inline Vector4& Vector4::operator+=(const Vector4& rhs) {
    x += rhs.x; y += rhs.y; z += rhs.z; w += rhs.w; return *this;
}

inline Vector4& Vector4::operator-=(const Vector4& rhs) {
    x -= rhs.x; y -= rhs.y; z -= rhs.z; w -= rhs.w; return *this;
}

inline Vector4 Vector4::operator*(const float a) const {
    return Vector4(x*a, y*a, z*a, w*a);
}

inline Vector4 Vector4::operator*(const Vector4& rhs) const {
    return Vector4(x*rhs.x, y*rhs.y, z*rhs.z, w*rhs.w);
}

inline Vector4& Vector4::operator*=(const float a) {
    x *= a; y *= a; z *= a; w *= a; return *this;
}

inline Vector4& Vector4::operator*=(const Vector4& rhs) {
    x *= rhs.x; y *= rhs.y; z *= rhs.z; w *= rhs.w; return *this;
}

inline Vector4 Vector4::operator/(const float a) const {
    return Vector4(x/a, y/a, z/a, w/a);
}

inline Vector4& Vector4::operator/=(const float a) {
    x /= a; y /= a; z /= a; w /= a; return *this;
}

inline bool Vector4::operator==(const Vector4& rhs) const {
    return (x == rhs.x) && (y == rhs.y) && (z == rhs.z) && (w == rhs.w);
}

inline bool Vector4::operator!=(const Vector4& rhs) const {
    return (x != rhs.x) || (y != rhs.y) || (z != rhs.z) || (w != rhs.w);
}

inline bool Vector4::operator<(const Vector4& rhs) const {
    if(x < rhs.x) return true;
    if(x > rhs.x) return false;
    if(y < rhs.y) return true;
    if(y > rhs.y) return false;
    if(z < rhs.z) return true;
    if(z > rhs.z) return false;
    if(w < rhs.w) return true;
    if(w > rhs.w) return false;
    return false;
}

inline float Vector4::operator[](int index) const {
    return (&x)[index];
}

inline float& Vector4::operator[](int index) {
    return (&x)[index];
}

inline void Vector4::set(float x, float y, float z, float w) {
    this->x = x; this->y = y; this->z = z; this->w = w;
}

inline float Vector4::length() const {
    return sqrtf(x*x + y*y + z*z + w*w);
}

inline float Vector4::distance(const Vector4& vec) const {
    return sqrtf((vec.x-x)*(vec.x-x) + (vec.y-y)*(vec.y-y) + (vec.z-z)*(vec.z-z) + (vec.w-w)*(vec.w-w));
}

inline Vector4& Vector4::normalize() {
    //NOTE: leave w-component untouched
    //@@const float EPSILON = 0.000001f;
    float xxyyzz = x*x + y*y + z*z;
    //@@if(xxyyzz < EPSILON)
    //@@    return *this; // do nothing if it is zero vector

    //float invLength = invSqrt(xxyyzz);
    float invLength = 1.0f / sqrtf(xxyyzz);
    x *= invLength;
    y *= invLength;
    z *= invLength;
    return *this;
}

inline float Vector4::dot(const Vector4& rhs) const {
    return (x*rhs.x + y*rhs.y + z*rhs.z + w*rhs.w);
}

inline bool Vector4::equal(const Vector4& rhs, float epsilon) const {
    return fabs(x - rhs.x) < epsilon && fabs(y - rhs.y) < epsilon &&
           fabs(z - rhs.z) < epsilon && fabs(w - rhs.w) < epsilon;
}

inline Vector4 operator*(const float a, const Vector4 vec) {
    return Vector4(a*vec.x, a*vec.y, a*vec.z, a*vec.w);
}

inline std::ostream& operator<<(std::ostream& os, const Vector4& vec) {
    os << "(" << vec.x << ", " << vec.y << ", " << vec.z << ", " << vec.w << ")";
    return os;
}
*/
// END OF VECTOR4 /////////////////////////////////////////////////////////////


