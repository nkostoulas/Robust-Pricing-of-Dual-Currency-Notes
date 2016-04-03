//
//  BSB.hpp
//  uncertainVolatilityModel
//
//  Created by Nick Kostoulas on 01/04/16.
//  Copyright © 2016 Nick Kostoulas. All rights reserved.
//

#ifndef BSB_hpp
#define BSB_hpp

#include <stdio.h>
#include <cmath>
#include <iostream>

#define MRT_STEPS 50

class BSB{
    double dt, sigmaMax, sigmaMin, r;
    double** F;
    int n;
    const double lowerP = sigmaMin*sigmaMin/(2*sigmaMax*sigmaMax);
    const double upperP = 0.5;
    const double step = (upperP - lowerP)/MRT_STEPS;
    
public:
    BSB(int n, double dt, double smax, double smin, double r, double** F);
    double upperBound();
    double lowerBound();
    double supremum(double U, double M, double D);
    double infinum(double U, double M, double D);
    double probU(double p);
    double probM(double p);
    double probD(double p);
};



#endif /* BSB_hpp */
