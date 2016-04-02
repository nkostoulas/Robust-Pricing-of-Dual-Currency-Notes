//
//  BSB.cpp
//  uncertainVolatilityModel
//
//  Created by Nick Kostoulas on 02/04/16.
//  Copyright © 2016 Nick Kostoulas. All rights reserved.
//

#include "BSB.hpp"

BSB::BSB(int n_, double dt_, double smax, double smin, double r_, double** F_):n(n_),dt(dt_),sigmaMax(smax),sigmaMin(smin),r(r_),F(F_){
    
}

double BSB::upperBound(){
    double Wp[2*n+1][n+1];
    double Lp = 0;
    for (int i=n; i>=0; i--){
        for (int j=0; j<2*i+1; j++){
            if (i==n){
                Wp[j][i] = F[j][i];
            }else{
                Lp = (1-sigmaMax*sqrt(dt)/2)*Wp[j+2][i+1] + (1+sigmaMax*sqrt(dt)/2)*Wp[j][i+1] - Wp[j+1][i+1];
                //std::cout<<Lp<<"    ";
                if (Lp>=0){
                    Wp[j][i] = F[j][i] + exp(-r*dt)*(Wp[j+1][i+1] + Lp/2);
                }else{
                    Wp[j][i] = F[j][i] + exp(-r*dt)*(Wp[j+1][i+1] + Lp*pow(sigmaMin,2)/(2*pow(sigmaMax,2)));
                }
                
            }
            //std::cout << Wp[j][i] << " ";
        }
        //std::cout << "\n";
    }
    return Wp[0][0];
    
}

double BSB::lowerBound(){
    double Wm[2*n+1][n+1];
    double Lm = 0;
    for (int i=n; i>=0; i--){
        for (int j=0; j<2*i+1; j++){
            if (i==n){
                Wm[j][i] = F[j][i];
            }else{
                Lm = (1-sigmaMax*sqrt(dt)/2)*Wm[j+2][i+1] + (1+sigmaMax*sqrt(dt)/2)*Wm[j][i+1] - Wm[j+1][i+1];
                if (Lm<0){
                    Wm[j][i] = F[j][i] + exp(-r*dt)*(Wm[j+1][i+1] + Lm/2);
                }else{
                    Wm[j][i] = F[j][i] + exp(-r*dt)*(Wm[j+1][i+1] + pow(sigmaMin,2)/(2*pow(sigmaMax,2))*Lm);
                }
            }
        }
    }
    return Wm[0][0];
}
