//
//  bsbFiniteDiff.cpp
//  uncertainVolatilityModel
//
//  Created by Nick Kostoulas on 08/06/16.
//  Copyright © 2016 Nick Kostoulas. All rights reserved.
//

#include "bsbFiniteDiff.hpp"

void bsbFiniteDiff(){
    //*********** Basic Definitions *************************//
    double r = 0.05;    //risk free interest rate
    double smax = 0.4;
    double smin = 0.1;
    
    double timeToExpiry = 0.5;
    double buyStrike = 90;
    double sellStrike = 100;
    
    double T = timeToExpiry;
    int N = 2000;
    int NS = 200;
    int S0 = 0;
    
    double dt = T/N;        //percentage of year for each period
    
    double timesExpiry = 4; //how many times expiry is Smax
    double dS = buyStrike/((NS-1)/timesExpiry);
    dS = 1;
    
    
    double** F = new double*[N];
        for(int i = 0; i < N; ++i)
        F[i] = new double[NS];
    
    // INITIALISE call spread example
    for (int i=0; i<N; i++){
        for (int j=0; j<NS; j++){
            
            if(i==0){       //at time 0 or time T equivalently calculate payoff
                if(j==NS-1){    //upper bound
                    F[i][j] = (sellStrike - buyStrike)*exp(-r*i*dt);
                }else if(j==0){ //lower bound
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
                if(j==NS-1){    //upper boundary
                    F[i][j] = (sellStrike - buyStrike)*exp(-r*i*dt);
                }else{          //lower boundary and all the rest initialised to 0
                    F[i][j] = 0;
                }
            }
        }
    }
    
    int plus = 0;
    int minus = 0;
    
    // FINITE DIFFERENCE
    double b = r*dt;
    
    for(int i=0; i<N-1; i++){
        for(int j=1; j<NS-1; j++){
            double gamma = (F[i][j+1] - 2*F[i][j] + F[i][j-1])/(dS*dS);
            double s;
            if(gamma>=0){
                s = smax;
                plus += 1;
            }else{
                s = smin;
                minus += 1;
            }
            double a = s*s*dt;
            
            F[i+1][j] = 0.5*(a*j*j - b*j)*F[i][j-1] + (1-a*j*j-b)*F[i][j] + 0.5*(a*j*j + b*j)*F[i][j+1];
        }
    }
    
    std::cout<<plus<<" "<<minus<<"\n";

    for(int j=0; j<NS; j++){
        std::cout<<S0+j*dS<<"\t"<<F[N-1][j]<<"\n";
    }
    
    
    
    
}