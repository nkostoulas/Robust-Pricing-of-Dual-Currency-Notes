//
//  toms462.hpp
//  uncertainVolatilityModel
//
//  Created by Nick Kostoulas on 16/06/16.
//  Copyright © 2016 Nick Kostoulas. All rights reserved.
//

#ifndef toms462_hpp
#define toms462_hpp

#include <stdio.h>

void bivariate_normal_cdf_values ( int &n_data, double &x, double &y,
                                  double &r, double &fxy );
double bivnor ( double ah, double ak, double r );
double gauss ( double t );
double r8_abs ( double x );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
void timestamp ( );

#endif /* toms462_hpp */
