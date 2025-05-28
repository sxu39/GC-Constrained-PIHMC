#ifndef CPIHMC_VEC3_H
#define CPIHMC_VEC3_H

namespace cpihmc
{
    template <class T>
    class vec3
    {
        public:
            T x, y, z;
        public:
            explicit vec3(const T &X = 0, const T &Y = 0, const T &Z = 0):x(X), y(Y), z(Z){}
            void set_zero(){x = y = z = 0;}
            virtual const vec3<T> &operator=(const vec3<T> &u){x = u.x; y = u.y; z = u.z; return *this;}
            virtual vec3<T> &operator+=(const vec3<T> &u){x += u.x; y += u.y; z += u.z; return *this;}
            virtual vec3<T> &operator-=(const vec3<T> &u){x -= u.x; y -= u.y; z -= u.z; return *this;}
            vec3<T> &operator*=(const T &s){x *= s; y *= s; z *= s; return *this;}
            vec3<T> &operator/=(const T &s){x /= s; y /= s; z /= s; return *this;}
            const vec3<T> operator-(){return vec3<T>(-x, -y, -z);}
            const T norm2() const {return x*x + y*y + z*z;}
            const T norm() const {return sqrt(norm2());}
    };

    template <class T> vec3<T> operator+(const vec3<T> &u, const vec3<T> &v){return vec3<T>(u.x+v.x, u.y+v.y, u.z+v.z);}
    template <class T> vec3<T> operator-(const vec3<T> &u, const vec3<T> &v){return vec3<T>(u.x-v.x, u.y-v.y, u.z-v.z);}
    template <class T> T operator*(const vec3<T> &u, const vec3<T> &v){return u.x*v.x + u.y*v.y + u.z*v.z;}
    template <class T> vec3<T> operator*(const vec3<T> &u, const T &s){return vec3<T>(u.x*s, u.y*s, u.z*s);}
    template <class T> vec3<T> operator*(const T &s, const vec3<T> &u){return vec3<T>(s*u.x, s*u.y, s*u.z);}
    template <class T> vec3<T> operator/(const vec3<T> &u, const T &s){return vec3<T>(u.x/s, u.y/s, u.z/s);}
}

#endif