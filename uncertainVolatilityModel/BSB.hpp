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

#define MRT_STEPS 2


class BSB{
private:
    double dt, sigmaMax, sigmaMin, r;
    double** F; //payoff 2D array based on trinomial tree implementation
    int n;  //number of periods
    const double lowerP = sigmaMin*sigmaMin/(2*sigmaMax*sigmaMax);
    const double upperP = 0.5;
    const double step = (upperP - lowerP)/MRT_STEPS;
    
public:    
    BSB(int n, double dt, double smax, double smin, double r, double** F);
    BSB(const BSB &source);
    virtual ~BSB();
    
    double upperBound() const;
    double lowerBound() const;
    double supremum(double U, double M, double D) const;
    double infinum(double U, double M, double D) const;
    inline double probU(double p) const{return p*(1-sigmaMax*sqrt(dt)/2);};
    inline double probM(double p) const{return 1 - 2*p;};
    inline double probD(double p) const{return p*(1+sigmaMax*sqrt(dt)/2);};
};



#endif /* BSB_hpp */
