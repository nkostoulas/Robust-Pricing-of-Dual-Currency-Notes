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
    
    int SELECTOR;
    cout<<"This program compares BS and BSB on pricing dual currency notes. \n";
    cout<<"Enter 1 for a Bull Call Spread of PRDC notes or 2 for a Calendar Spread of PRDC notes: ";
    cin >> SELECTOR;
    cout<<"\n";
    cout<<"Price\tBSUpper\tBSMid\tBSLower\tBSBLower\tBSBUpper\n";
    
    //*********** Basic Definitions *************************//
    double rf = 0.15;    //foreign interest rate
    double rd = 0.10;   //domestic interest rate
    
    double FX = 0.0;      //current underlying asset price
    double FX0 = 90;
    double FX1 = 110;
    
    double FXstrike0 = FX0*rd/rf;
    double L0 = rf/FX0;
    double FXstrike1 = FX1*rd/rf;
    double L1 = rf/FX1;
    
    //FXstrike0 = FX0;
    //L0 = 1;
    //FXstrike1 = FX1;
    //L1 = 1;

    double smax = 0.25;  //maximum volatility
    double smin = 0.05;  //minimum volatility
    double smid = 0.15;  //mid volatility
    
    double bsUpper = 0.0, bsMid=0.0, bsLower=0.0, bsbUpper=0.0, bsbLower=0.0;
    
    double T = 10; //years for trinomial tree
    double N = 66; //periods per year for trinomial tree IF I INCREASE EXC_BAD_ACCESS
    double callOffset = 0.5;    //time of short call option on calendar spread
    
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
 
    array<int,40> smatr= {1, 10, 20, 30, 35, 40, 45, 50, 55, 60, 65, 70, 72, 75, 78, 80, 82, 85,
        88, 90, 92, 95, 98, 100, 105, 110, 115, 120, 125, 130, 135, 140,
        145, 150, 155, 160, 170, 180, 190, 200};
    
    // Loop for all initial prices of underlying in smatr
    for (int sind=0; sind<smatr.size(); sind++){
        
        FX = smatr[sind];
        
        //********************  FX price dynamics *************//
        
        TrinomialTree tree(T, n, smax, smin, rd-rf, FX);
        //tree.printTree();
        
        if(SELECTOR == 1){
            //***************** CALL SPREAD Dual Currency Note Payoff  ******//
            callOffset = 0.0;
            
            for (int i=n; i>=0; i--){
                if(i%(int)N == 0 && i!=0){
                    for (int j=0; j<2*i+1; j++){
                        price = tree.nodePrice(i, j-i);

                        if(price>FXstrike0 && price<=FXstrike1){
                            F[i][j] = L0*(price - FXstrike0);
                        }else if(price>FXstrike1){
                            F[i][j] = L0*(price - FXstrike0) - L1*(price - FXstrike1);
                        }else{
                            F[i][j] = 0;
                        }
                    }
                }
            }
        }else{
            //***************** Calendar SPREAD Dual Currency Note Payoff  ******//
            
            for (int i=n; i>=0; i--){
                if(i%(int)N == 0 && i!=0){
                    
                    for (int j=0; j<2*i+1; j++){
                        price = tree.nodePrice(i, j-i);
                        
                        if(price>FXstrike0){
                            F[i][j] = L0*(price - FXstrike0);
                        }else{
                            F[i][j] = 0;
                        }
                    }
                    
                    int callIter = N*callOffset + i - N;
                    
                    for (int j=0; j<2*callIter+1; j++){
                        price = tree.nodePrice(callIter, j-callIter);
                        if(price>FXstrike1){
                            F[callIter][j] = L1*(FXstrike1 - price);
                        }else{
                            F[callIter][j] = 0;
                        }
                    }
                    
                }
            }
        }

        //****************** BSB note pricing ********//
        BSB bsb(n, dt, smax, smin, rd-rf, F);
        bsbUpper = bsb.upperBound();
        bsbLower = bsb.lowerBound();
        
        //****************** Black Scholes note pricing ********//
        
        for (int opts=1; opts<=(int)T; opts++){
            
            BS upperBuy(FXstrike0, FX, opts, 0, rd-rf, smax);
            BS midBuy(FXstrike0, FX, opts, 0, rd-rf, smid);
            BS lowerBuy(FXstrike0, FX, opts, 0, rd-rf, smin);

            double k  = opts - callOffset;
            
            BS upperSell(FXstrike1, FX, k, 0, rd-rf, smax);
            BS midSell(FXstrike1, FX, k, 0, rd-rf, smid);
            BS lowerSell(FXstrike1, FX, k, 0, rd-rf, smin);
            
            bsUpper += L0*upperBuy.callOptionPrice() - L1*lowerSell.callOptionPrice();
            bsMid += L0*midBuy.callOptionPrice() - L1*midSell.callOptionPrice();
            bsLower += L0*lowerBuy.callOptionPrice() - L1*upperSell.callOptionPrice();
            
        }
        
        cout<<FX<<"\t"<<bsUpper<<"\t"<<bsMid<<"\t"<<bsLower<<"\t"<<bsbLower<<"\t"<<bsbUpper<<"\n";
        
        bsbAsk<<bsbUpper<<"\n";
        bsbBid<<bsbLower<<"\n";
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
    
    system("sh scripts/dualCurr_script.sh");  
}