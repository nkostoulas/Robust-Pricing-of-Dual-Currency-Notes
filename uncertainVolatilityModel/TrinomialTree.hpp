//
//  TrinomialTree.hpp
//  uncertainVolatilityModel
//
//  Created by Nick Kostoulas on 01/04/16.
//  Copyright © 2016 Nick Kostoulas. All rights reserved.
//

#ifndef TrinomialTree_hpp
#define TrinomialTree_hpp

#include <stdio.h>
#include <cmath>
#include <iostream>

class TrinomialTree {
    double T, dt, N, sigmaMax, sigmaMin, r, So;
    double p_lower, p_upper;
    int lower, upper;
    
public:
    TrinomialTree(double T, double N, double smax, double smin, double r, double S);
    double nodePrice(int n, int j);
    void printTree();
    double* nodeProbability();
};



#endif /* TrinomialTree_hpp */

