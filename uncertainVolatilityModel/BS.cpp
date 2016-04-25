//
//  BS.cpp
//  uncertainVolatilityModel
//
//  Created by Nick Kostoulas on 02/04/16.
//  Copyright © 2016 Nick Kostoulas. All rights reserved.
//

#include "BS.hpp"

BS::BS(double strP, double underlP, double expT, double currT, double r, double sigma):E(strP),S(underlP),T(expT),t(currT),r(r),sigma(sigma){
    d1 = (log(S/E)+(r+sigma*sigma/2)*(T-t))/(sigma*sqrt(T-t));
    d2 = (log(S/E)+(r-sigma*sigma/2)*(T-t))/(sigma*sqrt(T-t));
    
}

BS::~BS(){
    
}

double BS::normalCDF(double x) const{
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
    
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);
    
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
    
    return 0.5*(1.0 + sign*y);
}

double BS::normalCDF2(double x) const{
    int neg = (x<0);
    if(neg) x *= -1;
    double k(1/(1+0.2316419*x));
    double y=((((1.330274429*k-1.821255978)*k+1.781477937)*k-0.356563782)*k+0.319381530)*k;
    y = 1.0 - 0.398942280401*exp(-0.5*x*x)*y;
    return (1-neg)*y + neg*(1-y);
    
}

double BS::callOptionPrice() const{
    
    double cdfn1 = normalCDF(d1);
    double cdfn2 = normalCDF(d2);
    
    return S*cdfn1 - E*exp(-r*(T-t))*cdfn2;
}

double BS::putOptionPrice() const{
    
    double cdfn1 = normalCDF(-d1);
    double cdfn2 = normalCDF(-d2);
    
    return E*exp(-r*(T-t))*cdfn2 - S*cdfn1;
}













