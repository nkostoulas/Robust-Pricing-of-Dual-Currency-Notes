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

class BSB{
    double dt, sigmaMax, sigmaMin, r;
    double** F;
    int n;
    
public:
    BSB(int n, double dt, double smax, double smin, double r, double** F);
    double upperBound();
    double lowerBound();
};



#endif /* BSB_hpp */
