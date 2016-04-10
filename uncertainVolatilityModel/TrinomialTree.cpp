//
//  TrinomialTree.cpp
//  uncertainVolatilityModel
//
//  Created by Nick Kostoulas on 02/04/16.
//  Copyright © 2016 Nick Kostoulas. All rights reserved.
//

#include "TrinomialTree.hpp"

TrinomialTree::TrinomialTree(double T_, double N_, double smax, double smin, double r_, double So_):T(T_),N(N_),sigmaMax(smax),sigmaMin(smin),r(r_),So(So_){
    dt = T/N;
    p_lower = pow(sigmaMin,2)/(2*pow(sigmaMax, 2));
    p_upper = 0.5;
    upper = (int)(p_upper*1e6);
    lower = (int)(p_lower*1e6);
}

double TrinomialTree::nodePrice(int n, int j){
    return So*exp(j*sigmaMax*sqrt(dt)+n*r*dt);
}

void TrinomialTree::printTree(){
    
     std::cout<<"Printing trinomial tree \n";
     for (int i=1; i<=N; i++) {
         for (int j=-i; j<=i; j++){
            std::cout << this->nodePrice(i,j) << " ";
         }
         std::cout<<"\n";
     }
}

double* TrinomialTree::nodeProbability(){
    /*
    double p = (lower + rand() % (upper - lower))*1e-6;
    double* prob = new double[3];
    prob[0] = p*(1-sigmaMax*sqrt(dt)/2);
    prob[1] = 1 - 2*p;
    prob[2] = p*(1+sigmaMax*sqrt(dt)/2);
    return prob;
     */
    return 0;
}


