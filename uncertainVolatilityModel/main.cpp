//
//  main.cpp
//  uncertainVolatilityModel
//
//  Created by Nick Kostoulas on 01/04/16.
//  Copyright © 2016 Nick Kostoulas. All rights reserved.
//

#include <iostream>
#include <cmath>
#include "TrinomialTree.hpp"
#include "BSB.hpp"
#include "BS.hpp"

int main(int argc, const char * argv[]) {
    
    //*********** basic definitions
    double r = 0.05;    //risk free interest rate
    double S = 75;      //current underlying asset price
    double smax = 0.4;  //maximum volatility
    double smin = 0.1;  //minimum volatility
    double smid = (smax+smin)/2;
    
    double timeToExpiry = 0.5;  // time to expiry of call spread
    double buyStrike = 90;      // strike price of call option bought
    double sellStrike = 100;    // strike price of call option sold
    
    //********************  BSB Spread pricing *************//
    
    double T = 1; //years
    double N = 2; //periods per year
    TrinomialTree tree(T, N, smax, smin, r, S);
    //tree.printTree();
    
    int n = N*timeToExpiry; //number of periods
    double dt = T/N;        //percentage of year for each period
    
    double** F = new double*[2*n+1];    //payoff structure matrix
    double price = 0;
    for(int i = 0; i < 2*n+1; ++i)
        F[i] = new double[n+1];
    
    // call spread payoff structure
    for (int i=n; i>=0; i--){
        for (int j=0; j<2*i+1; j++){
            if (i==n){
                price = tree.nodePrice(i, j-i);
                if(price>buyStrike && price<=sellStrike){
                    F[j][i] = price - buyStrike;
                }
                else if(price>sellStrike){
                    F[j][i] = sellStrike - buyStrike;
                }else{
                    F[j][i] = 0;
                }
            }
            else{
                F[j][i] = 0;
          }
            std::cout << F[j][i] << " ";
        }
        std::cout<<"\n";
        
    }
    
    //upper and lower value of example using BSB
    BSB bsb(n, dt, smax, smin, r, F);

    std::cout<<"For an initial stock price "<<S<<":"<<"\n";
    std::cout<<"W+ value is "<<bsb.upperBound()<<"\n";
    std::cout<<"W- value is "<<bsb.lowerBound()<<"\n";
    
    
    //****************** Black Scholes spread pricing ********//
    
    double timeToEval = 0.0;
    
    BS upperBuy(buyStrike, S, timeToExpiry, timeToEval, r, smax);
    BS lowerBuy(buyStrike, S, timeToExpiry, timeToEval, r, smin);
    BS midBuy(buyStrike, S, timeToExpiry, timeToEval, r, smid);
    
    BS upperSell(sellStrike, S, timeToExpiry, timeToEval, r, smax);
    BS lowerSell(sellStrike, S, timeToExpiry, timeToEval, r, smin);
    BS midSell(sellStrike, S, timeToExpiry, timeToEval, r, smid);
    
    std::cout<<"Upper buy - lower sell is "<<upperBuy.callOptionPrice()-lowerSell.callOptionPrice()<<"\n";
    std::cout<<"Lower buy - upper sell is "<<lowerBuy.callOptionPrice()-upperSell.callOptionPrice()<<"\n";
    std::cout<<"Mid buy - mid sell is "<<midBuy.callOptionPrice()-midSell.callOptionPrice()<<"\n";
    
}






