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

template<typename T>
class Vector4;

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
    Vector3<T>(const Vector4<T>& p4);

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
using float3 = Vector3<float>;
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


/////////////////// vector 4D ///////////////////////////
template<typename T>
class Vector4 {
public:
	/*These values are public to give easier access for
	outside funtions. Accessors and mutators are not really
	necessary*/
	T x;
	T y;
	T z;
	T w;

	//Constructor assigns the inputs to x and y.
	Vector4(T a=T(), T b=T(), T c=T(), T d = T()) : x(a), y(b), z(c), w(d) {}
	Vector4(T *p) :x(p[0]), y(p[1]), z(p[2]), w(p[3]) {}
	Vector4(Vector3<T> p) : x(p.x), y(p.y), z(p.z), w(0.f) {}

	T operator[] (size_t i) const {
		assert(i <= 3);
		if (i == 0)
			return x;
		if (i == 1)
			return y;
		if (i == 2)
			return z;
			
		return w;
	}
	T& operator[] (size_t i) {
		assert(i <= 3);
		if (i == 0)
			return x;
		if (i == 1)
			return y;
		if (i == 2)
			return z;
		return w;
	}	

	std::string toString() {
		std::stringstream s;
		s << "(" << x << ", " << y << "," << z << ","<<w<<")";
		return s.str();
	}
	/*The following operators simply return Vector4ds that
	have operations performed on the relative (x, y, z, w) values*/
	Vector4& operator+=(const Vector4<T>& v) { x += v.x; y += v.y; z += v.z; w += v.w; return *this; }
	Vector4& operator-=(const Vector4<T>& v) { x -= v.x; y -= v.y; z -= v.z; w -= v.w; return *this; }
	Vector4& operator*=(const Vector4<T>& v) { x *= v.x; y *= v.y; z *= v.z; w *= v.w; return *this; }
	Vector4& operator/=(const Vector4<T>& v) { x /= v.x; y /= v.y; z /= v.z; w /= v.w; return *this; }

	Vector4& operator+=(const float s) { x += s; y += s; z += s; return *this; }
	Vector4& operator-=(const float s) { x -= s; y -= s; z -= s; return *this; }
	Vector4& operator*=(const float s) { x *= s; y *= s; z *= s; return *this; }
	Vector4& operator/=(const float s) { x /= s; y /= s; z /= s; return *this; }

	//Check if the Vectors have the same values (uses pairwise comparison of 
	// 'std::tuple' on the x, y values of L and R.)
	friend bool operator==(const Vector4& L, const Vector4& R) {
		return std::tie(L.x, L.y, L.z, L.w) == std::tie(R.x, R.y, R.z, R.w);
	}
	friend bool operator!=(const Vector4& L, const Vector4& R) { return !(L == R); }

	void set(T a, T b, T c) { x = a; y = b; z = c; }
	/*Check which Vectors are closer or further from the origin.*/
	friend bool operator>(const Vector4& L, const Vector4& R) { return LengthSq(L) < LengthSq(R); }
	friend bool operator>=(const Vector4& L, const Vector4& R) { return !(L > R); }
	friend bool operator<(const Vector4& L, const Vector4& R) { return R < L; }
	friend bool operator<=(const Vector4& L, const Vector4& R) { return !(R < L); }

	//Negate both the x and y values.
	Vector3<T> operator-() const { return Vector3d<T>(-x, -y, -z, -w); }

	//Apply scalar operations.
	Vector4 operator*(T s) { Vector4 tmp(*this); tmp.x *= s; tmp.y *= s; tmp.z *= s; tmp.w *= s;  return tmp; }
	Vector4 operator/(T s) { Vector4 tmp(*this); tmp.x /= s; tmp.y /= s; tmp.z /= s tmp.w /= s;  return tmp; }

	//Returns the length of the vector from the origin.
	double Length() const { return sqrt(x*x + y * y + z * z + w*w); }
	double LengthSq()const { return x * x + y * y + z * z + w*w; }

};
using float4 = Vector4<float>;
template<class T> Vector4<T> operator*(const T& s, const Vector4<T>& v) { return Vector4<T>(v) *= s; }
template<class T> Vector4<T> operator*(const Vector4<T>& v, const T& s) { return Vector4<T>(v) *= s; }

template<class T> Vector4<T>  operator-(const Vector4<T>& v1, const Vector4<T>& v2) { return Vector4<T>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z, v1.w-v2.w); }
template<class T> Vector4<T>  operator+(const Vector4<T>& v1, const Vector4<T>& v2) { return Vector4<T>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z, v1.w+ v2.w); }

//Product functions
template<class T> T Dot(const Vector4<T>& a, const Vector4<T>& b) { return  ((a.x * b.x) + (a.y * b.y) + (a.z * b.z) + a.w*b.w); }

template<typename T>
inline Vector3<T>::Vector3(const Vector4<T>& p4) : x(p4.x), y(p4.y), z(p4.z)
{
}



