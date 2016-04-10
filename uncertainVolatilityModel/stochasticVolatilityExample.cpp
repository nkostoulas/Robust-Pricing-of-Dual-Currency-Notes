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
    const int Nsim = 100000;
    const int Nhedge = 100;
    double T = 0.5;             //time horizon
    int NsimYear = 200;        //sim per year
    double Dt = pow(NsimYear,-1);
    
    
    double randNorm = 0.0, randNorm2 = 0.0;
    //std::default_random_engine generator((unsigned int)time(NULL));
    std::mt19937 generator((unsigned int)time(NULL));
    std::normal_distribution<double> distribution(0, 1);
    
    //Stock Simulation Parameters
    array<double,Nsim> stockSim;
    stockSim[0] = 90;
    double miu = 0.1;
    
    //Volatility Simulation Parameters
    array<double,Nsim> volSim;
    volSim[0] = 0.2;
    double a = log(2);
    double g = log(0.2);
    double Z95 = 1.64;
    double ro = 2*sqrt(log(2))/Z95;

    int count = 0;
    double outsideConfInt = 0.0;
    
    //**** Monte Carlo Simulation of Stock price and Volatility
    for(int i=1; i<Nsim; i++){
        
        if(i%Nhedge==0){
            stockSim[i] = stockSim[0];
            volSim[i] = volSim[0];
        }else{
            randNorm = distribution(generator);
            randNorm2 = distribution(generator);
            volSim[i] = exp(log(volSim[i-1]) + a*(g - log(volSim[i-1]))*Dt + ro*randNorm*sqrt(Dt));
            stockSim[i] = stockSim[i-1] + miu*Dt*stockSim[i-1] + volSim[i-1]*randNorm2*sqrt(Dt)*stockSim[i-1];
        }
        
        if(volSim[i]>0.4 || volSim[i]<0.1){
            count++;
        }
    }
    
    outsideConfInt = (double)100*count/Nsim;
    cout<<count<<"%"<<"\t";
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
    //double bsbUpper2 = 0.0;
    //double bsbLower2 = 0.0;
    double sigmaAvg = 0.2214;
    
    //** Hedging
    
    for(int k=0; k<Nsim; k++){
        if(k%Nhedge==0){
            per = Nhedge;
            timeToExp = T;
        }
        switch (1){
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
                myfile<<bsbUpper - currPrice<<"\n";
                
     
                //BSB bsb2(per, Dt, volSim[k], smin, r, F);
                //bsbUpper2 = bsb2.upperBound();
                //bsbLower2 = bsb2.lowerBound();
                //myfile<<bsbUpper - bsbUpper2<<"\n";
     
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
        }
    }
 
    myfile.close();
    */
}

