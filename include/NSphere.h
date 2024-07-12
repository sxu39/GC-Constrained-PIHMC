//
// Created by Jin Bin on 2021/08/16.
//

#ifndef MD_MC_NSPHERE_H
#define MD_MC_NSPHERE_H

#include <cmath>
#include <iostream>
#include <ctime>
#include "vector3.h"
using namespace std;

#define pi acos(-1)

class RanNSphere{
public:
    int N_dim = 3; // system DOF(degree of freedom) dimension
    double step; // maximum step length of the trial moves
    double sigma; // sigma of one dimension Gaussian distribution
    RanNSphere(double step, double sig = 1);// initialize the class, main parameter is dimension
    Vector3<double> ran_vector() const; // main function, provide a vector on N dimension sphere homogeneously
private:
    double RanGaussian()
    // provide a number with Gaussian distribution
    const{
        return sigma * std::sqrt(-2*std::log(Random()))*cos(2*pi*Random());
    }
    static double Random(){
        return rand()/(double)RAND_MAX;
    }
};


#endif //MD_MC_NSPHERE_H
