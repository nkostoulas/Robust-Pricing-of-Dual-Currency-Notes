//
//  stochasticVolatilityExample.cpp
//  uncertainVolatilityModel
//
//  Created by Nick Kostoulas on 08/04/16.
//  Copyright © 2016 Nick Kostoulas. All rights reserved.
//

#include "stochasticVolatilityExample.hpp"

void stochasticVolatilityExample(){
    //General Parameters
    const int Npath = 10000;
    const int Nsim = 10000;
    const int Nhedge = 101;
    
    double dtSIM = pow(Nsim,-1);
    
    double T = 0.5;             //time horizon
    int NsimYear = Nhedge/T;        //sim per year
    double dtHEDGE = pow(NsimYear,-1);
    int step = Nsim/(Nhedge-1);
    
    double randNorm = 0.0, randNorm2 = 0.0;
    std::mt19937 generator((unsigned int)time(NULL));
    std::normal_distribution<double> distribution(0, 1);
    
    //Stock Simulation Parameters
    array<double,Nsim> stockSim;
    array<double,Nsim> stockSimConstVol;
    stockSim[0] = 95;
    stockSimConstVol[0] = 95;
    double miu = 0.1;
    double sigmaAvg = 0.2214;
    
    //Volatility Simulation Parameters
    array<double,Nsim> volSim;
    volSim[0] = 0.2;
    double a = log(2);
    double g = log(0.2);
    double Z95 = 3;  //or 1.64 / 1.96
    double ro = 2*sqrt(log(2))/Z95;

    int count = 0;
    int k = 0;
    
    //////////
    ofstream myfile;
    myfile.open ("/Users/nkostoulas/Documents/xcode_/uncertainVolatilityModel/uncertainVolatilityModel/example.txt");
    
    //Pricing parameters
    double r = 0.05;    //risk free interest rate
    double smax = 0.4;  //maximum volatility
    double smin = 0.1;  //minimum volatility
    
    double buyStrike = 90;      // strike price of call option bought
    double sellStrike = 100;    // strike price of call option sold
    int per = Nhedge;
    double timeToExp = T;
    double bsbLower_ = 0.0;
    double bsbLower = 0.0;
    int d = 1; //decider
    
    //**** Monte Carlo Simulation of Stock price and Volatility
    for(int PATH = 0; PATH < Npath; PATH++){
        per = Nhedge;
        timeToExp = T;
        for(int i=1; i<Nsim; i++){
            randNorm = distribution(generator);
            randNorm2 = distribution(generator);
            volSim[i] = exp(log(volSim[i-1]) + a*(g - log(volSim[i-1]))*dtSIM + ro*randNorm*sqrt(dtSIM));
            stockSim[i] = stockSim[i-1] + miu*dtSIM*stockSim[i-1] + volSim[i-1]*randNorm2*sqrt(dtSIM)*stockSim[i-1];
            //stockSimConstVol[i] = stockSimConstVol[i-1] + miu*dtSIM*stockSimConstVol[i-1] + sigmaAvg*randNorm2*sqrt(dtSIM)*stockSimConstVol[i-1];
            
            if(volSim[i]>0.4 || volSim[i]<0.1){
                count++;
            }
            
            //NO HEDGE AT TIME 0
            if((i+1)%step==0 && i!=Nsim-1){ //
                if(d==1){
                    TrinomialTree tree(1, per, smax, smin, r, stockSim[i]);
                    
                    double** F = new double*[2*per+1];    //payoff structure matrix
                    double price = 0;
                    for(int i = 0; i < 2*per+1; ++i)
                        F[i] = new double[per+1];
                    for (int i=per; i>=0; i--){
                        for (int j=0; j<2*i+1; j++){
                            if (i==per){
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
                        }
                    }
                    
                    BSB bsb(per, dtHEDGE, smax, smin, r, F);
                    bsbLower = bsb.lowerBound();
                    
                    if(volSim[i]<smax){
                        BSB bsb2(per, dtHEDGE, smax, volSim[i], r, F);
                        bsbLower_ = bsb2.lowerBound();
                    }else{
                        BSB bsb2(per, dtHEDGE, volSim[i], smax, r, F);
                        bsbLower_ = bsb2.lowerBound();
                    }

                    myfile<<bsbLower_ - bsbLower<<"\n";

                    
                    per--;
                    timeToExp -= dtHEDGE;
                }else{
                    // comparison of BS with constant and stoch volatility
                    BS currBuy(buyStrike, stockSim[i], timeToExp, 0, r, volSim[i]);
                    BS currSell(sellStrike, stockSim[i], timeToExp, 0, r, volSim[i]);
                    BS avgBuy(buyStrike, stockSim[i], timeToExp, 0, r, sigmaAvg);
                    BS avgSell(sellStrike, stockSim[i], timeToExp, 0, r, sigmaAvg);
                    per--;
                    timeToExp -= dtHEDGE;
                    
                    double currPrice = currBuy.callOptionPrice() - currSell.callOptionPrice();
                    double avgPrice = avgBuy.callOptionPrice() - avgSell.callOptionPrice();
                    //myfile<<avgPrice - currPrice<<"\n";
                    
                    if((avgPrice-currPrice)>=-0.1 && (avgPrice-currPrice)<=0.1){
                        if(k<0){
                            myfile<<avgPrice - currPrice<<"\n";
                        }
                        k--;
                    }else{
                        myfile<<avgPrice - currPrice<<"\n";
                    }
                }
            }
        }
    }
    cout<<(double)100*count/(Npath*Nsim)<<"%"<<"\t";
    
    /*
    ofstream myfile;
    myfile.open ("/Users/nkostoulas/Documents/xcode_/uncertainVolatilityModel/uncertainVolatilityModel/example.txt");

    //Pricing parameters
    double r = 0.05;    //risk free interest rate
    double smax = 0.4;  //maximum volatility
    double smin = 0.1;  //minimum volatility
    
    double buyStrike = 90;      // strike price of call option bought
    double sellStrike = 100;    // strike price of call option sold
    int per = Nhedge;
    double dt = T/Nhedge;
    double timeToExp = 0.5;
    
    double bsbUpper = 0.0;
    double bsbLower = 0.0;
    
    
    //** Hedging
    
    for(int k=0; k<1; k++){
        if(k%Nhedge==0){
            per = Nhedge;
            timeToExp = T;
        }
        switch (4){
            case 1:
            {
                TrinomialTree tree(1, per, smax, smin, r, stockSim[k]);
                
                double** F = new double*[2*per+1];    //payoff structure matrix
                double price = 0;
                for(int i = 0; i < 2*per+1; ++i)
                    F[i] = new double[per+1];
                for (int i=per; i>=0; i--){
                    for (int j=0; j<2*i+1; j++){
                        if (i==per){
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
                    }
                }
                
                
                BSB bsb(per, Dt, smax, smin, r, F);
                bsbUpper = bsb.upperBound();
                bsbLower = bsb.lowerBound();
                BS currBuy(buyStrike, stockSim[k], timeToExp, 0, r, volSim[k]);
                BS currSell(sellStrike, stockSim[k], timeToExp, 0, r, volSim[k]);
                double currPrice = currBuy.callOptionPrice() - currSell.callOptionPrice();
                if(volSim[k]>(smax-smin)/2){
                    myfile<<bsbUpper - currPrice<<"\n";
                }else{
                    myfile<<currPrice - bsbLower<<"\n";
                }
                per--;
                timeToExp -= dt;
     
                break;
            }
            case 2:
            {
                // comparison of BS with constant and stoch volatility
                BS currBuy(buyStrike, stockSim[k], timeToExp, 0, r, volSim[k]);
                BS currSell(sellStrike, stockSim[k], timeToExp, 0, r, volSim[k]);
                BS avgBuy(buyStrike, stockSim[k], timeToExp, 0, r, sigmaAvg);
                BS avgSell(sellStrike, stockSim[k], timeToExp, 0, r, sigmaAvg);
                per--;
                timeToExp -= dt;

                double currPrice = currBuy.callOptionPrice() - currSell.callOptionPrice();
                double avgPrice = avgBuy.callOptionPrice() - avgSell.callOptionPrice();
                myfile<<avgPrice - currPrice<<"\n";
                break;
            }
            case 3:
            {
                BS currBuy(buyStrike, stockSimConstVol[k], timeToExp, 0, r, volSim[k]);
                BS currSell(sellStrike, stockSimConstVol[k], timeToExp, 0, r, volSim[k]);
                BS avgBuy(buyStrike, stockSimConstVol[k], timeToExp, 0, r, sigmaAvg);
                BS avgSell(sellStrike, stockSimConstVol[k], timeToExp, 0, r, sigmaAvg);
                per--;
                timeToExp -= dt;
                
                double currPrice = currBuy.callOptionPrice() - currSell.callOptionPrice();
                double avgPrice = avgBuy.callOptionPrice() - avgSell.callOptionPrice();
                myfile<<avgPrice - currPrice<<"\n";
                break;
                
                
            }
            case 4:
                
            break;
        }
    }
 
    myfile.close();
     */
    
}

