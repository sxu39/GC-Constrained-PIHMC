#ifndef MD_MC_VECTOR3_H
#define MD_MC_VECTOR3_H

#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

template <class T>
class Vector3
{
public:
    T x;
    T y;
    T z;

    explicit Vector3(const T &x1 = 0,const T &y1 = 0,const T &z1 = 0) :x(x1),y(y1),z(z1){};
    void set(const T &x1, const T &y1,const T &z1) { x = x1; y = y1; z = z1; }

    Vector3<T>& operator =(const Vector3<T> &u) { x=u.x; y=u.y; z=u.z;     return *this; }
    Vector3<T>& operator+=(const Vector3<T> &u) { x+=u.x; y+=u.y; z+=u.z;  return *this; }
    Vector3<T>& operator-=(const Vector3<T> &u) { x-=u.x; y-=u.y; z-=u.z;  return *this; }
    Vector3<T>& operator*=(const Vector3<T> &u);
    Vector3<T>& operator*=(const T &s)          { x*=s; y*=s; z*=s;        return *this; }
    Vector3<T>& operator/=(const Vector3<T> &u);
    Vector3<T>& operator/=(const T &s)          { x/=s; y/=s; z/=s;        return *this; }
    Vector3<T>& operator%=(const T &s); // modulo operation of the vector
    Vector3<T>  operator -() const              { return Vector3<T>(-x,-y,-z); }			// Peize Lin add 2017-01-10

    T norm2() const	{ return x*x + y*y + z*z; }
    T norm() const	{ return sqrt(norm2()); }
    void normalize(){ const T m=norm(); x/=m; y/=m; z/=m; }
    void reverse()	{ x=-x; y=-y; z=-z; }

    void print()const ;		// mohan add 2009-11-29
};

template <class T>
Vector3<T>& Vector3<T>::operator%=(const T &s){
    x -= floor(x/s)*s;
    y -= floor(y/s)*s;
    z -= floor(z/s)*s;
    return *this;
}

template <class T> Vector3<T> operator%( const Vector3<T> &u, const T &s ) { return Vector3<T>( u.x-floor(u.x/s)*s, u.y-floor(u.y/s)*s, u.z-floor(u.z/s)*s ); }

template <class T> Vector3<T> operator+( const Vector3<T> &u, const Vector3<T> &v ) { return Vector3<T>( u.x+v.x, u.y+v.y, u.z+v.z ); }
template <class T> Vector3<T> operator-( const Vector3<T> &u, const Vector3<T> &v ) { return Vector3<T>( u.x-v.x, u.y-v.y, u.z-v.z ); }
//u.v=(ux*vx)+(uy*vy)+(uz*vz)
template <class T> T          operator*( const Vector3<T> &u, const Vector3<T> &v ) { return ( u.x*v.x + u.y*v.y + u.z*v.z ); }
template <class T> Vector3<T> operator*( const T &s,          const Vector3<T> &u ) { return Vector3<T>( u.x*s, u.y*s, u.z*s ); }
template <class T> Vector3<T> operator*( const Vector3<T> &u, const T &s          ) { return Vector3<T>( u.x*s, u.y*s, u.z*s ); } // mohan add 2009-5-10
template <class T> Vector3<T> operator/( const Vector3<T> &u, const T &s          ) { return Vector3<T>( u.x/s, u.y/s, u.z/s ); }
//u.v=(ux*vx)+(uy*vy)+(uz*vz)
template <class T> T          dot      ( const Vector3<T> &u, const Vector3<T> &v ) { return ( u.x*v.x + u.y*v.y + u.z*v.z ); }
// | i  j  k  |
// | ux uy uz |
// | vx vy vz |
// u.v=(uy*vz-uz*vy)i+(-ux*vz+uz*vx)j+(ux*vy-uy*vx)k
template <class T> Vector3<T> operator^(const Vector3<T> &u,const Vector3<T> &v)
{
    return Vector3<T> ( u.y * v.z - u.z * v.y,
                        -u.x * v.z + u.z * v.x,
                        u.x * v.y - u.y * v.x);
}
// | i  j  k  |
// | ux uy uz |
// | vx vy vz |
// u.v=(uy*vz-uz*vy)i+(-ux*vz+uz*vx)j+(ux*vy-uy*vzx)k
template <class T> Vector3<T> cross(const Vector3<T> &u,const Vector3<T> &v)
{
    return Vector3<T> ( u.y * v.z - u.z * v.y,
                        -u.x * v.z + u.z * v.x,
                        u.x * v.y - u.y * v.x);
}
//s = u.(v x w)
//template <class T> T TripleScalarProduct(Vector3<T> u, Vector3<T> v, Vector3<T> w)
//{
//	return T((u.x * (v.y * w.z - v.z * w.y)) +
//	         (u.y * (-v.x * w.z + v.z * w.x)) +
//	         (u.z * (v.x * w.y - v.y * w.x)));
//}

//whether m1 != m2
template <class T> bool operator !=(const Vector3<T> &u, const Vector3<T> &v){ return !(u == v); }
//whether u == v
template <class T> bool operator ==(const Vector3<T> &u, const Vector3<T> &v)
{
    if(u.x == v.x && u.y == v.y && u.z == v.z)
        return true;
    return false;
}

template <class T> void Vector3<T>::print()const
{
    cout.precision(5) ;
    cout << "(" << setw(10) << x << "," << setw(10) << y << ","
         << setw(10) << z  << ")"  << endl;
}
template<class T> static std::ostream & operator << ( std::ostream &os, const Vector3<T> &u )
{
    os << setprecision(12);
    os << u.x << " " << u.y << " " << u.z;
    return os;
}



#endif