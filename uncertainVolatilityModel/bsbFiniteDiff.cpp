//
//  bsbFiniteDiff.cpp
//  uncertainVolatilityModel
//
//  Created by Nick Kostoulas on 08/06/16.
//  Copyright © 2016 Nick Kostoulas. All rights reserved.
//

#include "bsbFiniteDiff.hpp"

void bsbFiniteDiff(){
    // Variable definitions
    double r = 0.05;
    double smax = 0.4;
    double smin = 0.1;
    
    double timeToExpiry = 0.5;
    double buyStrike = 90;
    double sellStrike = 100;
    
    double T = timeToExpiry;
    int N = 2000;
    int NS = 200;
    int S0 = 0;
    
    double dt = T/N;
    
    double timesExpiry = 4;
    double dS = buyStrike/((NS-1)/timesExpiry);
    dS = 1;
    
    double** F = new double*[N];
        for(int i = 0; i < N; ++i)
        F[i] = new double[NS];
    
    // Initial and boudnary conditions
    for (int i=0; i<N; i++){
        for (int j=0; j<NS; j++){
            
            if(i==0){
                if(j==NS-1){
                    F[i][j] = (sellStrike - buyStrike)*exp(-r*i*dt);
                }else if(j==0){
                    F[i][j] = 0;
                }else{
                    if(S0+j*dS>buyStrike && S0+j*dS<=sellStrike){
                        F[i][j] = S0 + j*dS - buyStrike;
                    }else if(S0+j*dS>sellStrike){
                        F[i][j] = (sellStrike - buyStrike)*exp(-r*i*dt);
                    }else{
                        F[i][j] = 0;
                    }
                }
            }else{
                if(j==NS-1){
                    F[i][j] = (sellStrike - buyStrike)*exp(-r*i*dt);
                }else{
                    F[i][j] = 0;
                }
            }
        }
    }
    
    // Explicit finite difference to calculate upper and lower bounds
    double b = r*dt;
    
    for(int i=0; i<N-1; i++){
        for(int j=1; j<NS-1; j++){
            double gamma = (F[i][j+1] - 2*F[i][j] + F[i][j-1])/(dS*dS);
            double s;
            if(gamma>=0){
                s = smax;
            }else{
                s = smin;
            }
            double a = s*s*dt;
            
            F[i+1][j] = 0.5*(a*j*j - b*j)*F[i][j-1] + (1-a*j*j-b)*F[i][j] + 0.5*(a*j*j + b*j)*F[i][j+1];
        }
    }

    for(int j=0; j<NS; j++){
        std::cout<<S0+j*dS<<"\t"<<F[N-1][j]<<"\n";
    } 
}