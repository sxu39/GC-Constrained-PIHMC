//
// Created by Jin Bin on 2021/08/16.
//

#include "NSphere.h"

RanNSphere::RanNSphere(double step, double sig):step(step), sigma(sig){}

Vector3<double> RanNSphere::ran_vector() const{
    Vector3<double> vector;
    double Norm;
    do {
        vector.x = RanGaussian();
        vector.y = RanGaussian();
        vector.z = RanGaussian();
        Norm = vector.norm();
    }while (Norm == 0 or isinf(Norm)); // if the vector is 0 or has infinity number, produce vector again.

    // change the vector length
    double length = Random() * step;
        vector *= length / Norm;
    return vector;
}