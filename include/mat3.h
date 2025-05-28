#ifndef CPIHMC_MAT3_H
#define CPIHMC_MAT3_H

#include "vec3.h"

namespace cpihmc
{
    template <class T>
    class mat3
    {
        public:
            T e11, e12, e13, e21, e22, e23, e31, e32, e33;
        public:
            explicit mat3(const T &E11 = 1, const T &E12 = 0, const T &E13 = 0,
                          const T &E21 = 0, const T &E22 = 1, const T &E23 = 0,
                          const T &E31 = 0, const T &E32 = 0, const T &E33 = 1):
                          e11(E11), e12(E12), e13(E13),
                          e21(E21), e22(E22), e23(E23),
                          e31(E31), e32(E32), e33(E33){}
            const mat3<T> &operator+=(const mat3<T> &);
            const mat3<T> &operator-=(const mat3<T> &);
            const T det() const;
            const mat3<T> inverse() const;
    };

    template <class T> const mat3<T> operator*(const mat3<T> &Mat, const prec_t S); // Mat * S
    template <class T> const vec3<T> operator*(const vec3<T> &Vec, const mat3<T> &Mat); // Vec * Mat
}

template <class T>
const cpihmc::mat3<T> &cpihmc::mat3<T>::operator+=(const mat3<T> &Mat)
{
    e11 += Mat.e11; e12 += Mat.e12; e13 += Mat.e13;
    e21 += Mat.e21; e22 += Mat.e22; e23 += Mat.e23;
    e31 += Mat.e31; e32 += Mat.e32; e33 += Mat.e33;
    return *this;
}

template <class T>
const cpihmc::mat3<T> &cpihmc::mat3<T>::operator-=(const mat3<T> &Mat)
{
    e11 -= Mat.e11; e12 -= Mat.e12; e13 -= Mat.e13;
    e21 -= Mat.e21; e22 -= Mat.e22; e23 -= Mat.e23;
    e31 -= Mat.e31; e32 -= Mat.e32; e33 -= Mat.e33;
    return *this;
}

template <class T>
const T cpihmc::mat3<T>::det() const
{
    return e11 * e22 * e33 - e11 * e32 * e23 + e21 * e32 * e13 - e21 * e12 * e33 + e31 * e12 * e23 - e31 * e22 * e13;
}

template <class T>
const cpihmc::mat3<T> cpihmc::mat3<T>::inverse() const
{
    prec_t Det = det();
    if (Det == 0) Det = 1;
    return mat3<T>((e22 * e33 - e23 * e32) / Det,
                   -(e12 * e33 - e13 * e32) / Det,
                   (e12 * e23 - e13 * e22) / Det,
                   -(e21 * e33 - e23 * e31) / Det,
                   (e11 * e33 - e13 * e31) / Det,
                   -(e11 * e23 - e13 * e21) / Det,
                   (e21 * e32 - e22 * e31) / Det,
                   -(e11 * e32 - e12 * e31) / Det,
                   (e11 * e22 - e12 * e21) / Det);
}

template <class T>
const cpihmc::mat3<T> cpihmc::operator*(const cpihmc::mat3<T> &Mat, const cpihmc::prec_t S)
{
    return cpihmc::mat3<T>(Mat.e11 * S, Mat.e12 * S, Mat.e13 * S,
                           Mat.e21 * S, Mat.e22 * S, Mat.e23 * S,
                           Mat.e31 * S, Mat.e32 * S, Mat.e33 * S);
}

// Vec * Mat
template <class T>
const cpihmc::vec3<T> cpihmc::operator*(const cpihmc::vec3<T> &Vec, const cpihmc::mat3<T> &Mat)
{
    return cpihmc::vec3<T>(Vec.x * Mat.e11 + Vec.y * Mat.e21 + Vec.z * Mat.e31,
                           Vec.x * Mat.e12 + Vec.y * Mat.e22 + Vec.z * Mat.e32,
                           Vec.x * Mat.e13 + Vec.y * Mat.e23 + Vec.z * Mat.e33);
}

#endif