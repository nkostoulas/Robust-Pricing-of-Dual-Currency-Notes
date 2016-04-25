//
//  BS.hpp
//  uncertainVolatilityModel
//
//  Created by Nick Kostoulas on 02/04/16.
//  Copyright © 2016 Nick Kostoulas. All rights reserved.
//

#ifndef BS_hpp
#define BS_hpp

#include <stdio.h>
#include <cmath>

class BS{
    double E, S, T, t, r, sigma;
    double d1, d2;
    
public:
    BS(double strP, double underlP, double expT, double currT, double r, double sigma);
    BS(const BS &source);
    virtual ~BS();
    
    double callOptionPrice() const;
    double putOptionPrice() const;
    double normalCDF(double x) const;
    double normalCDF2(double x) const;
    
};

#endif /* BS_hpp */
