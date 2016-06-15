//
//  calendarSpreadExample.cpp
//  uncertainVolatilityModel
//
//  Created by Nick Kostoulas on 02/04/16.
//  Copyright © 2016 Nick Kostoulas. All rights reserved.
//
#include "calendarSpreadExample.hpp"

void calendarSpreadExample(){
    ofstream bsbAsk("./data/bsbAskClnd.txt");
    ofstream bsbBid("./data/bsbBidClnd.txt");
    ofstream bsAsk("./data/bsAskClnd.txt");
    ofstream bsBid("./data/bsBidClnd.txt");
    ofstream bsMidl("./data/bsMidClnd.txt");
    ofstream prices("./data/pricesClnd.txt");
    
    cout<<"This program compares BS and BSB on a Calendar Spread of 2 options. \n";
    cout<<"Price\tBSUpper\tBSMid\tBSLower\tBSBLower\tBSBUpper\n";
    
    //*********** Basic Definitions *************************//
    double r = 0.05;    //risk free interest rate
    double S = 90;      //current underlying asset price
    double smax = 0.4;  //maximum volatility
    double smin = 0.1;  //minimum volatility
    double smid = (smax+smin)/2;
    
    double bsUpper = 0.0, bsMid=0.0, bsLower=0.0, bsbUpper=0.0, bsbLower=0.0;
    
    double timeToExpiry1 = 1;
    double timeToExpiry2 = 0.5;
    double buyStrike = 90;      // strike price of call option bought
    double sellStrike = 100;    // strike price of call option sold
    
    double T = 1; //years
    double N = 252; //periods per year

    int n = N*timeToExpiry1; //number of periods
    double dt = T/N;        //percentage of year for each period
    double timeToEval = 0.0;
    
    double** F = new double*[n+1];
    double price = 0;
    for(int i = 0; i < n+1; ++i)
        F[i] = new double[2*i+1];
    for (int i=n; i>=0; i--){
        for (int j=0; j<2*i+1; j++){
            F[i][j] = 0;
        }
    }
    
    array<int,40> smatr= {0, 10, 20, 30, 40, 50, 55, 60, 65, 70, 72, 75, 78, 80, 82, 85,
        88, 90, 92, 95, 98, 100, 105, 110, 115, 120, 125, 130, 135, 140,
        145, 150, 155, 160, 165, 170, 175, 180, 190, 200 };
    
    for (int sind=0; sind<smatr.size(); sind++){
        
        S = smatr[sind];
        
        //********************  BSB Spread pricing *************//

        TrinomialTree tree(T, N, smax, smin, r, S);
        //tree.printTree();
        
        //***************** Calendar Spread Payoff Structure ******//
        for (int i=n; i>=0; i--){
            for (int j=0; j<2*i+1; j++){
                if (i==n){
                    price = tree.nodePrice(i, j-i);
                    F[i][j] = max(price - buyStrike, 0.0);
                }else if(i==n/2){
                    price = tree.nodePrice(i, j-i);
                    F[i][j] = -max(price - sellStrike, 0.0);
                }
            }
        }
        
        //****************** BSB spread pricing ********//
        BSB bsb(n, dt, smax, smin, r, F);
        bsbUpper = bsb.upperBound();
        bsbLower = bsb.lowerBound();
        
        //****************** Black Scholes spread pricing ********//
        BS upperBuy(buyStrike, S, timeToExpiry1, timeToEval, r, smax);  //ASK
        BS lowerBuy(buyStrike, S, timeToExpiry1, timeToEval, r, smin);  //BID
        BS midBuy(buyStrike, S, timeToExpiry1, timeToEval, r, smid);
        
        BS upperSell(sellStrike, S, timeToExpiry2, timeToEval, r, smax); //BID
        BS lowerSell(sellStrike, S, timeToExpiry2, timeToEval, r, smin); //ASK
        BS midSell(sellStrike, S, timeToExpiry2, timeToEval, r, smid);
        
        bsUpper = upperBuy.callOptionPrice() - lowerSell.callOptionPrice();
        bsMid = lowerBuy.callOptionPrice() - upperSell.callOptionPrice();
        bsLower = midBuy.callOptionPrice() - midSell.callOptionPrice();
        
        prices << S << "\n";
        bsbAsk << bsbUpper << "\n";
        bsbBid << bsbLower << "\n";
        bsAsk <<  bsUpper<< "\n";
        bsBid << bsMid << "\n";
        bsMidl << bsLower << "\n";
        
        cout<<S<<"\t"<<bsUpper<<"\t"<<bsMid<<"\t"<<bsLower<<"\t"<<bsbLower<<"\t"<<bsbUpper<<"\n";
        
    }
    for(int i = 0; i < n+1; ++i)
        delete [] F[i];
    delete [] F;
    
    bsbAsk.close();
    bsbBid.close();
    bsAsk.close();
    bsBid.close();
    bsMidl.close();
    prices.close();
}