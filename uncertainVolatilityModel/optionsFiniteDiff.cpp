//
//  optionsFiniteDiff.cpp
//  uncertainVolatilityModel
//
//  Created by Nick Kostoulas on 08/06/16.
//  Copyright © 2016 Nick Kostoulas. All rights reserved.
//

#include "optionsFiniteDiff.hpp"
#include "BS.hpp"

void optionsFiniteDiff(){
    //*********** Basic Definitions *************************//
    double r = 0.1;    //risk free interest rate
    double s = 0.4;
    
    double timeToExpiry = 0.25;
    double buyStrike = 10;
    
    double T = timeToExpiry;
    int N = 5;//7000;
    int NS = 10;//201;
    int S0 = 0;
    
    double dt = T/N;        //percentage of year for each period

    double timesExpiry = 4; //how many times expiry is Smax
    double dS = buyStrike/((NS-1)/timesExpiry);
    
    double** F = new double*[N];
    
    for(int i = 0; i < N; ++i)
        F[i] = new double[NS];
    
    // INITIALISE
    for (int i=0; i<N; i++){
        for (int j=0; j<NS; j++){
            
            if(i==0){
                if(j==NS-1){
                    F[i][j] = S0 + (NS-1)*dS - buyStrike*exp(-r*i*dt);
                }else if(j==0){
                    F[i][j] = 0;
                }else{
                    if(S0+j*dS>buyStrike){
                        F[i][j] = S0 + j*dS - buyStrike;
                    }else{
                        F[i][j] = 0;
                    }
                }
            }else{
                if(j==NS-1){
                    F[i][j] = S0 + (NS-1)*dS - buyStrike*exp(-r*i*dt);
                }else{
                    F[i][j] = 0;
                }
            }
            std::cout<<F[i][j]<<"\t";
        }
        std::cout<<"\n";
    }
    
    
    // FINITE DIFFERENCE
    for(int i=0; i<N-1; i++){
        for(int j=1; j<NS-1; j++){
            
            F[i+1][j] = 0.5*(s*s*j*j*dt - r*j*dt)*F[i][j-1] + (1-s*s*j*j*dt-r*dt)*F[i][j] + 0.5*(s*s*j*j*dt + r*j*dt)*F[i][j+1];
        }
    }
    
    for(int j=0; j<NS; j++){
        BS upperBuy(buyStrike, S0+j*dS, timeToExpiry, 0, r, s);
        std::cout<<S0+j*dS<<"\t"<<upperBuy.callOptionPrice()<<"\t"<<F[N-1][j]<<"\n";
    }
    
    
    
    
    
}