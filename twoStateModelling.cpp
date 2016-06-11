//
//  twoStateModelling.cpp
//  uncertainVolatilityModel
//
//  Created by Nick Kostoulas on 10/06/16.
//  Copyright © 2016 Nick Kostoulas. All rights reserved.
//

#include "twoStateModelling.hpp"
struct priceArray{
    double upperValue;
    double lowerValue;
};
pair<double, double> findMaxValue(double d1, double d2, double d12, double smax1, double smax2, double smin1, double smin2)
{
    return std::make_pair(smax1, smax2);
}
pair<double, double> findMinValue(double d1, double d2, double d12, double smax1, double smax2, double smin1, double smin2)
{
    return std::make_pair(smin1, smin2);
}

void twoStateModelling(){
    
    double r = 0.05;
    double smax1 = 0.5, smax2 = 0.5, smin1 = 0.3, smin2 = 0.3;
    double s1 = smax1;
    double s2 = smax2;
    double corr = 0.3;
    
    double T = 0.5;
    double E = 40;

    int NT = 2500;  //need HIGH FOR STABILITY comment
    int NS1 = 100;
    int NS2 = 100;
    
    double dt = T/NT;
    double dS = 1;
    double df = r*dt;
    
    // INITIALISE PRICE ARRAY
    vector<vector<vector<double>>> Fu, Fl;
    Fu.resize(NT);
    Fl.resize(NT);
    for (int i = 0; i < NT; ++i) {
        Fu[i].resize(NS1);
        Fl[i].resize(NS1);
        for (int j = 0; j < NS1; ++j){
            Fu[i][j].resize(NS2);
            Fl[i][j].resize(NS2);
        }
    }
    
    // *******************  Set Boundary and Initial Conditions for the payoff **********************************
    //Initial
    for(int j=0; j<NS1; j++){
        for(int k=0; k<NS2; k++){
            double p1 = j*dS;
            double p2 = k*dS;
            Fu[0][j][k] = max(max(p1,p2)-E,0.0);
        }
    }
    //Boundary
    for(int i=1; i<NT; i++){
        Fu[i][0][0] = 0;
        Fu[i][NS1-1][NS2-1] = max((NS2-1)*dS,(NS1-1)*dS) - E*exp(-r*i*dt);
        int j = 0;  //left boundary
        for(int k=1; k<NS2-1;k++){
            Fu[i][0][k] = Fu[i-1][0][k] + (0.5*s2*s2*k*k*dt)*(Fu[i-1][0][k+1] - 2*Fu[i-1][0][k] + Fu[i-1][0][k-1])
                     -df*Fu[i-1][0][k] + (0.5*r*k*dt)*(Fu[i-1][0][k+1] - Fu[i-1][0][k-1]);
        }
        int k = 0;  //right boundary
        for(int j=1; j<NS1-1; j++){
            Fu[i][j][0] = Fu[i-1][j][0] + (0.5*s1*s1*j*j*dt)*(Fu[i-1][j+1][0] - 2*Fu[i-1][j][0] + Fu[i-1][j-1][0])
                     -df*Fu[i-1][j][0] + (0.5*r*j*dt)*(Fu[i-1][j+1][0] - Fu[i-1][j-1][0]);
        }
        j = NS1-1;  //s1 = max boundary
        double p = j*dS - E*exp(-r*i*dt);
        for(int k=0; k<NS2; k++){
            Fu[i][j][k] = p;
            
        }
        k = NS2-1;  //s2 = max boundary
        p = k*dS - E*exp(-r*i*dt);
        for(int j=0; j<NS1; j++){
            Fu[i][j][k] = p;
        }
    }
    Fl = Fu;

    // *******************  Calculate derivative price using finite difference ************************************
    int cnt = 0;
    
    vector<vector<priceArray>> price(NS1, vector<priceArray>(NS2));
    
    for(int i=0; i<NT-1; i++){
        for(int j=1; j<NS1-1; j++){
            for(int k=1; k<NS2-1; k++){
                
                double d1u = (Fu[i][j+1][k] - 2*Fu[i][j][k] + Fu[i][j-1][k]);
                double d2u = (Fu[i][j][k+1] - 2*Fu[i][j][k] + Fu[i][j][k-1]);
                double d12u = (Fu[i][j+1][k+1] - Fu[i][j+1][k-1] - Fu[i][j-1][k+1] + Fu[i][j-1][k-1])/(4);
                
        
                tie(s1,s2) = findMaxValue(d1u,d2u,d12u,smax1,smax2,smin1,smin2);
                corr = 0.3;
                
                Fu[i+1][j][k] = (1-df)*Fu[i][j][k] + (0.5*s1*s1*j*j*dt)*d1u
                                                 + (0.5*s2*s2*k*k*dt)*d2u
                                                 + (0.5*r*k*dt)*(Fu[i][j][k+1] - Fu[i][j][k-1])
                                                 + (0.5*r*j*dt)*(Fu[i][j+1][k] - Fu[i][j-1][k])
                                                 + (0.25*corr*s1*s2*j*k*dt)*4*d12u;
                
                double d1l = (Fl[i][j+1][k] - 2*Fl[i][j][k] + Fl[i][j-1][k]);
                double d2l = (Fl[i][j][k+1] - 2*Fl[i][j][k] + Fl[i][j][k-1]);
                double d12l = (Fl[i][j+1][k+1] - Fl[i][j+1][k-1] - Fl[i][j-1][k+1] + Fl[i][j-1][k-1])/(4);
                
                
                tie(s1,s2) = findMinValue(d1l,d2l,d12l,smax1,smax2,smin1,smin2);
                corr = 0.5;
                
                Fl[i+1][j][k] = (1-df)*Fl[i][j][k] + (0.5*s1*s1*j*j*dt)*d1l
                                                    + (0.5*s2*s2*k*k*dt)*d2l
                                                    + (0.5*r*k*dt)*(Fl[i][j][k+1] - Fl[i][j][k-1])
                                                    + (0.5*r*j*dt)*(Fl[i][j+1][k] - Fl[i][j-1][k])
                                                    + (0.25*corr*s1*s2*j*k*dt)*4*d12l;
                
                double det = d1u*d2u - d12u*d12u;
                if(det==0 && d1u!=0 && d2u!=0 && d12u!=0){
                    cnt+=1;
                    // so apparently this is very rare and only happens when second derivatives are super small or zero
                    //cout<<F[i+1][j][k]<<"\t"<<d1<<"\t"<<d2<<"\t"<<d12<<"\t"<<det<<"\n";
                }
                if(i==NT-2){
                    price[j][k].upperValue = Fu[i+1][j][k];
                }
            }
        }
    }
    cout << cnt << "\n";
    //for(int i=1; i<NS1-1; i++){
      //  for(int j=1; j<NS2-2; j++){
        //    cout<<price[i][j].upperValue<<"\n";
        //}
    //}
    cout<<"40 40"<<" "<<Fl[NT-1][40][40]<<"\t"<<Fu[NT-1][40][40]<<"\n";
    
    
    
        
        
        
}