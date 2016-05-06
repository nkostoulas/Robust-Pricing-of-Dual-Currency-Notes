//
//  dualCurrencyNoteExample.cpp
//  uncertainVolatilityModel
//
//  Created by Nick Kostoulas on 04/05/16.
//  Copyright © 2016 Nick Kostoulas. All rights reserved.
//

#include "dualCurrencyNoteExample.hpp"

void dualCurrencyNoteExample(){
    ofstream bsbAsk("./data/bsbAskDual.txt");
    ofstream bsbBid("./data/bsbBidDual.txt");
    ofstream bsUp("./data/bsAskDual.txt");
    ofstream bsMd("./data/bsMidDual.txt");
    ofstream bsLow("./data/bsBidDual.txt");
    ofstream prices("./data/pricesDual.txt");
    
    //*********** Basic Definitions *************************//
    double rf = 0.15;    //foreign interest rate
    double rd = 0.10;   //domestic interest rate
    
    double FX = 0.0;      //current underlying asset price
    double FX0 = 100;
    double FXstrike = FX0*rd/rf;
    double L = rf/FX0;
    
    double smax = 0.05;  //maximum volatility
    double smin = 0.001;  //minimum volatility
    double smid = 0.01;
    smax = 0.4;
    smin = 0.2;
    smid = 0.3;
    
    double bsUpper = 0.0;
    double bsMid = 0.0;
    double bsLower = 0.0;
    
    
    double T = 10; //years for trinomial tree
    double N = 65; //periods per year for trinomial tree IF I INCREASE EXC_BAD_ACCESS
    
    
    int n = N*T; //number of periods
    double dt = T/n;
    
    double** F = new double*[n+1];
    double price = 0;
    for(int i = 0; i < n+1; ++i)
        F[i] = new double[2*i+1];
    for (int i=n; i>=0; i--){
        for (int j=0; j<2*i+1; j++){
            F[i][j] = 0;
        }
    }
 
    array<int,27> smatr= {50, 55, 60, 65, 70, 72, 75, 78, 80, 82, 85, 88, 90, 92,
                            95, 98, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150 };

    // Loop for all initial prices of underlying in smatr
    for (int sind=0; sind<smatr.size(); sind++){
        
        FX = smatr[sind];
        
        //********************  FX price dynamics *************//
        
        TrinomialTree tree(T, n, smax, smin, rd-rf, FX);
        //tree.printTree();
        
        //***************** Dual Currency Note Payoff  ******//
        for (int i=n; i>=0; i--){
            if(i%(int)N == 0 && i!=0){
                for (int j=0; j<2*i+1; j++){
                    price = tree.nodePrice(i, j-i);

                    if(price>FXstrike){
                        F[i][j] = L*(price - FXstrike);
                    }else{
                        F[i][j] = 0;
                    }
                }
            }
        }
        
        
        //****************** BSB note pricing ********//
        BSB bsb(n, dt, smax, smin, rd-rf, F);
        
        //****************** Black Scholes note pricing ********//
        
        for (int opts=1; opts<=(int)T; opts++){
            
            //BS upperBuy(FXstrike, FX, opts, 0, rd, rf, smax);
            //bsUpper += upperBuy.callFXOptionPrice();
            BS upperBuy(FXstrike, FX, opts, 0, rd-rf, smax);
            bsUpper += L*upperBuy.callOptionPrice();
            
            //BS midBuy(FXstrike, FX, opts, 0, rd, rf, smid);
            //bsMid += midBuy.callFXOptionPrice();
            BS midBuy(FXstrike, FX, opts, 0, rd-rf, smid);
            bsMid += L*midBuy.callOptionPrice();
            
            //BS lowerBuy(FXstrike, FX, opts, 0, rd, rf, smin);
            //bsLower += lowerBuy.callFXOptionPrice();
            BS lowerBuy(FXstrike, FX, opts, 0, rd-rf, smin);
            bsLower += L*lowerBuy.callOptionPrice();
        }
        
        cout<<FX<<"   ";
        cout<<bsUpper<<" ";
        cout<<bsMid<<" ";
        cout<<bsLower<<" ";
        cout<< bsb.lowerBound()<<" ";
        cout<<bsb.upperBound()<<"  ";
        cout<<"\n";
        
        bsbAsk<<bsb.upperBound()<<"\n";
        bsbBid<<bsb.lowerBound()<<"\n";
        bsUp<<bsUpper<<"\n";
        bsMd<<bsMid<<"\n";
        bsLow<<bsLower<<"\n";
        prices<<FX<<"\n";
        
        bsUpper = 0.0;
        bsMid = 0.0;
        bsLower = 0.0;
        
    }
    
    for(int i = 0; i < n+1; ++i)
        delete [] F[i];
    delete [] F;
    
    bsbAsk.close();
    bsbBid.close();
    bsUp.close();
    bsMd.close();
    bsLow.close();
    prices.close();

    
    
    
    
}