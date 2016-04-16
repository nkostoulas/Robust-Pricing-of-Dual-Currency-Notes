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
    ofstream bsMid("./data/bsMidClnd.txt");
    ofstream prices("./data/pricesClnd.txt");
    
    int smatr[40] = {0, 10, 20, 30, 40, 50, 55, 60, 65, 70, 72, 75, 78, 80, 82, 85,
        88, 90, 92, 95, 98, 100, 105, 110, 115, 120, 125, 130, 135, 140,
        145, 150, 155, 160, 165, 170, 175, 180, 190, 200 };
    
    for (int sind=0; sind<40; sind++){
        
        //*********** basic definitions
        double r = 0.05;    //risk free interest rate
        double S = smatr[sind];      //current underlying asset price
        double smax = 0.4;  //maximum volatility
        double smin = 0.1;  //minimum volatility
        double smid = (smax+smin)/2;
        
        double timeToExpiry1 = 1;
        double timeToExpiry2 = 0.5;
        double buyStrike = 90;      // strike price of call option bought
        double sellStrike = 100;    // strike price of call option sold
        
        //********************  BSB Spread pricing *************//
        
        double T = 1; //years
        double N = 252; //periods per year
        TrinomialTree tree(T, N, smax, smin, r, S);
        //tree.printTree();
        
        int n = N*timeToExpiry1; //number of periods
        double dt = T/N;        //percentage of year for each period
        
        double** F = new double*[2*n+1];    //payoff structure matrix
        double price = 0;
        for(int i = 0; i < 2*n+1; ++i)
            F[i] = new double[n+1];
        
        // calendar spread payoff structure
        for (int i=n; i>=0; i--){
            for (int j=0; j<2*i+1; j++){
                if (i==n){
                    price = tree.nodePrice(i, j-i);
                    
                    if(price>buyStrike){
                        F[j][i] = price - buyStrike;
                    }else{
                        F[j][i] = 0;
                    }
                }else if(i==n/2){
                    price = tree.nodePrice(i, j-i);
                    
                    if(price>sellStrike){
                        F[j][i] = -(price - sellStrike);
                    }else{
                        F[j][i] = 0;
                    }
                }else{
                    F[j][i] = 0;
                }
            }
        }
        
        //upper and lower value of example using BSB
        BSB bsb(n, dt, smax, smin, r, F);
        
        //****************** Black Scholes spread pricing ********//
        
        double timeToEval = 0.0;
        
        BS upperBuy(buyStrike, S, timeToExpiry1, timeToEval, r, smax);
        BS lowerBuy(buyStrike, S, timeToExpiry1, timeToEval, r, smin);
        BS midBuy(buyStrike, S, timeToExpiry1, timeToEval, r, smid);
        
        BS upperSell(sellStrike, S, timeToExpiry2, timeToEval, r, smax);
        BS lowerSell(sellStrike, S, timeToExpiry2, timeToEval, r, smin);
        BS midSell(sellStrike, S, timeToExpiry2, timeToEval, r, smid);
        
        prices << S << "\n";
        bsbAsk << bsb.upperBound() << "\n";
        bsbBid << bsb.lowerBound() << "\n";
        bsAsk << upperBuy.callOptionPrice() - lowerSell.callOptionPrice() << "\n";
        bsBid << lowerBuy.callOptionPrice() - upperSell.callOptionPrice() << "\n";
        bsMid << midBuy.callOptionPrice() - midSell.callOptionPrice() << "\n";
        /*
        cout<<S<<"   ";
        cout<<bsb.upperBound()<<" ";
        cout<<bsb.lowerBound()<<"  ";
        cout<<upperBuy.callOptionPrice()-lowerSell.callOptionPrice()<<" ";
        cout<<lowerBuy.callOptionPrice()-upperSell.callOptionPrice()<<"  ";
        cout<<midBuy.callOptionPrice()-midSell.callOptionPrice();
        cout<<"\n";
         */
    }
    bsbAsk.close();
    bsbBid.close();
    bsAsk.close();
    bsBid.close();
    bsMid.close();
    prices.close();
}