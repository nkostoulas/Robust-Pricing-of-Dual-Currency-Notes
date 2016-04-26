//
//  callSpreadExample.cpp
//  uncertainVolatilityModel
//
//  Created by Nick Kostoulas on 02/04/16.
//  Copyright © 2016 Nick Kostoulas. All rights reserved.
//
#include "callSpreadExample.hpp"

void callSpreadExample(){
    ofstream bsbAsk("./data/bsbAskCall.txt");
    ofstream bsbBid("./data/bsbBidCall.txt");
    ofstream bsAsk("./data/bsAskCall.txt");
    ofstream bsBid("./data/bsBidCall.txt");
    ofstream bsMid("./data/bsMidCall.txt");
    ofstream bsLow("./data/bsLowCall.txt");
    ofstream bsHigh("./data/bsHighCall.txt");
    ofstream prices("./data/pricesCall.txt");
    
    //*********** Basic Definitions *************************//
    double r = 0.05;    //risk free interest rate
    double S = 90;      //current underlying asset price
    double smax = 0.4;  //maximum volatility
    double smin = 0.1;  //minimum volatility
    double smid = (smax+smin)/2;
    
    double timeToExpiry = 0.5;  // time to expiry of call spread
    double buyStrike = 90;      // strike price of call option bought
    double sellStrike = 100;    // strike price of call option sold
    
    double T = 1; //years for trinomial tree
    double N = 252; //periods per year for trinomial tree
    
    int n = N*timeToExpiry; //number of periods
    double dt = T/N;        //percentage of year for each period
    double timeToEval = 0.0;//time at which black scholes is evaluated
    
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
    
    // Loop for all initial prices of underlying in smatr
    for (int sind=0; sind<smatr.size(); sind++){
        
        S = smatr[sind];

        //********************  BSB Spread pricing *************//
        
        TrinomialTree tree(T, N, smax, smin, r, S);
        //tree.printTree();
        
        //***************** Call Spread Payoff Structure ******//
        for (int j=0; j<2*n+1; j++){
            price = tree.nodePrice(n, j-n);
            if(price>buyStrike && price<=sellStrike){
                F[n][j] = price - buyStrike;
            }
            else if(price>sellStrike){
                F[n][j] = sellStrike - buyStrike;
            }else{
                F[n][j] = 0;
            }
        }
        //****************** BSB spread pricing ********//
        BSB bsb(n, dt, smax, smin, r, F);

        //****************** Black Scholes spread pricing ********//

        BS upperBuy(buyStrike, S, timeToExpiry, timeToEval, r, smax);   //ASK
        BS lowerBuy(buyStrike, S, timeToExpiry, timeToEval, r, smin);   //BID
        BS midBuy(buyStrike, S, timeToExpiry, timeToEval, r, smid);

        BS upperSell(sellStrike, S, timeToExpiry, timeToEval, r, smax); //BID
        BS lowerSell(sellStrike, S, timeToExpiry, timeToEval, r, smin); //ASK
        BS midSell(sellStrike, S, timeToExpiry, timeToEval, r, smid);
        
        prices << S << "\n";
        bsbAsk << bsb.upperBound() << "\n";
        bsbBid << bsb.lowerBound() << "\n";
        bsAsk << upperBuy.callOptionPrice() - lowerSell.callOptionPrice() << "\n";
        bsBid << lowerBuy.callOptionPrice() - upperSell.callOptionPrice() << "\n";
        bsMid << midBuy.callOptionPrice() - midSell.callOptionPrice() << "\n";
        bsLow << lowerBuy.callOptionPrice() - lowerSell.callOptionPrice()<< "\n";
        bsHigh << upperBuy.callOptionPrice() - upperSell.callOptionPrice() << "\n";
        
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
    for(int i = 0; i < n+1; ++i)
        delete [] F[i];
    delete [] F;
    
    bsbAsk.close();
    bsbBid.close();
    bsAsk.close();
    bsBid.close();
    bsMid.close();
    bsLow.close();
    bsHigh.close();
    prices.close();
}