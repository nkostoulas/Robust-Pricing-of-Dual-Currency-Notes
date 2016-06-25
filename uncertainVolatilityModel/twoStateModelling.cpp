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
double normalCDF(double x){
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
    
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);
    
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
    
    return 0.5*(1.0 + sign*y);
}

//Calculate value of Vanilla European option on the minimum of two assets
double minOptionValue(double S1, double S2, double E, double s1, double s2, double corr, double r, double T){
    double g1 = (log(S2/E)+(r-0.5*s2*s2)*T)/(s2*sqrt(T));
    double g2 = (log(S1/E)+(r-0.5*s1*s1)*T)/(s1*sqrt(T));
    double s = sqrt(s1*s1 + s2*s2 - 2*corr*s1*s2);
    double rc = (corr*s1 - s2)/s;
    
    double a1 = g1 + s2*sqrt(T);
    double a2 = (log(S1/S2)-0.5*s*s*sqrt(T))/(s*sqrt(T));
    double b1 = g2 + s1*sqrt(T);
    double b2 = (log(S2/S1)-0.5*s*s*sqrt(T))/(s*sqrt(T));
    
    return S2*bivnor(-a1, -a2, rc) + S1*bivnor(-b1, -b2, rc) - E*exp(-r*T)*bivnor(-g1, -g2, corr);
}
//Calculate value of Vanilla European option on the maximum of two assets
double maxOptionValue(double S1, double S2, double E, double s1, double s2, double corr, double r, double T){
    BS bs1(E, S1, T, 0, r, s1);
    BS bs2(E, S2, T, 0, r, s2);
    return bs1.callOptionPrice() + bs2.callOptionPrice() - minOptionValue(S1, S2, E, s1, s2, corr, r, T);
}
//Calculate determinant of second derivative matrix
double valueToOpt(double d1, double d2, double d12, double s1, double s2){
    return d1*s1*s1 + 2*d12*s1*s2 + d2*s2*s2;
}
//Find the volatility pair that maximises the 2-asset option spread
pair<double, double> findMaxValue(double d1, double d2, double d12, double smax1, double smax2, double smin1, double smin2, double corr){
    d12 = d12*corr;
    double det = d1*d2 - d12*d12;
    double max = -INT_MAX, val = 0.0;
    double s1=smax1, s2=smax2, s1_return = smax1, s2_return=smax2;

    val = valueToOpt(d1, d2, d12, s1=smax1, s2=smax2);
    if(val > max){
        if((det==0 && (d1*s1==-d12*s2)) || det!=0){
                max = val;
                s1_return = s1;
                s2_return = s2;
        }
    }
    val = valueToOpt(d1, d2, d12, s1=smin1, s2=smin2);
    if(val > max){
        if((det==0 && (d1*s1==-d12*s2)) || det!=0){
            max = val;
            s1_return = s1;
            s2_return = s2;
        }
    }
    val = valueToOpt(d1, d2, d12, s1=smax1, s2=smin2);
    if(val > max){
        if((det==0 && (d1*s1==-d12*s2)) || det!=0){
            max = val;
            s1_return = s1;
            s2_return = s2;
        }
    }
    val = valueToOpt(d1, d2, d12, s1=smin1, s2=smax2);
    if(val > max){
        if((det==0 && (d1*s1==-d12*s2)) || det!=0){
            max = val;
            s1_return = s1;
            s2_return = s2;
        }
    }
    val = valueToOpt(d1, d2, d12, s1=smax1, s2=-d12*smax1/d2);
    if(val > max && s2>=smin2 && s2<=smax2){
        if((det==0 && (d1*s1==-d12*s2)) || det!=0){
            max = val;
            s1_return = s1;
            s2_return = s2;
        }
    }
    val = valueToOpt(d1, d2, d12, s1=smin1, s2=-d12*smin1/d2);
    if(val > max && s2>=smin2 && s2<=smax2){
        if((det==0 && (d1*s1==-d12*s2)) || det!=0){
            max = val;
            s1_return = s1;
            s2_return = s2;
        }
    }
    val = valueToOpt(d1, d2, d12, s1=-d12*smax2/d1, s2=smax2);
    if(val > max && s1>=smin1 && s1<=smax1){
        if((det==0 && (d1*s1==-d12*s2)) || det!=0){
            max = val;
            s1_return = s1;
            s2_return = s2;
        }
    }
    val = valueToOpt(d1, d2, d12, s1=-d12*smin2/d1, s2=smin2);
    if(val > max && s1>=smin1 && s1<=smax1){
        if((det==0 && (d1*s1==-d12*s2)) || det!=0){
            max = val;
            s1_return = s1;
            s2_return = s2;
        }
    }
    
    return std::make_pair(s1_return, s2_return);
}
//Find the volatility pair that minimises the 2-asset option spread
pair<double, double> findMinValue(double d1, double d2, double d12, double smax1, double smax2, double smin1, double smin2, double corr){
    d12 = d12*corr;
    double det = d1*d2 - d12*d12;
    double min = INT_MAX, val = 0.0;
    double s1=smin1, s2=smin2, s1_return = smin1, s2_return=smin2;

    val = valueToOpt(d1, d2, d12, s1=smax1, s2=smax2);
    if(val < min){
        if((det==0 && (d1*s1==-d12*s2)) || det!=0){
            min = val;
            s1_return = s1;
            s2_return = s2;
        }
    }
    val = valueToOpt(d1, d2, d12, s1=smin1, s2=smin2);
    if(val < min){
        if((det==0 && (d1*s1==-d12*s2)) || det!=0){
            min = val;
            s1_return = s1;
            s2_return = s2;
        }
    }
    val = valueToOpt(d1, d2, d12, s1=smax1, s2=smin2);
    if(val < min){
        if((det==0 && (d1*s1==-d12*s2)) || det!=0){
            min = val;
            s1_return = s1;
            s2_return = s2;
        }
    }
    val = valueToOpt(d1, d2, d12, s1=smin1, s2=smax2);
    if(val < min){
        if((det==0 && (d1*s1==-d12*s2)) || det!=0){
            min = val;
            s1_return = s1;
            s2_return = s2;
        }
    }
    val = valueToOpt(d1, d2, d12, s1=smax1, s2=-d12*smax1/d2);
    if(val < min && s2>=smin2 && s2<=smax2){
        if((det==0 && (d1*s1==-d12*s2)) || det!=0){
            min = val;
            s1_return = s1;
            s2_return = s2;
        }
    }
    val = valueToOpt(d1, d2, d12, s1=smin1, s2=-d12*smin1/d2);
    if(val < min && s2>=smin2 && s2<=smax2){
        if((det==0 && (d1*s1==-d12*s2)) || det!=0){
            min = val;
            s1_return = s1;
            s2_return = s2;
        }
    }
    val = valueToOpt(d1, d2, d12, s1=-d12*smax2/d1, s2=smax2);
    if(val < min && s1>=smin1 && s1<=smax1){
        if((det==0 && (d1*s1==-d12*s2)) || det!=0){
            min = val;
            s1_return = s1;
            s2_return = s2;
        }
    }
    val = valueToOpt(d1, d2, d12, s1=-d12*smin2/d1, s2=smin2);
    if(val < min && s1>=smin1 && s1<=smax1){
        if((det==0 && (d1*s1==-d12*s2)) || det!=0){
            min = val;
            s1_return = s1;
            s2_return = s2;
        }
    }

    return std::make_pair(s1_return, s2_return);
}

//Butterfly spread under the Uncertain Volatility Model using Finite Difference
void Bspread_BSB_bounds(double r, double smax1, double smax2, double smin1, double smin2, double corr, double T, double E1, double E2, double NT, double NS, double dS){
    ofstream bsbUpper("./data/2Dbsb_spread_U.txt");
    ofstream bsbLower("./data/2Dbsb_spread_L.txt");
    ofstream prices("./data/2Dbsb_spread_prices.txt");

    double s1 = smax1, s2 = smax2;
    double Eavg = 0.5*(E1+E2);
    int NS1 = NS, NS2 = NS;
    
    double dt = T/NT;
    double df = r*dt;
    
    // Initialise Price Array
    vector<vector<double>> Fu, Fu_n, Fl, Fl_n;
    Fu.resize(NS1);
    Fu_n.resize(NS1);
    Fl.resize(NS1);
    Fl_n.resize(NS1);
    for (int j = 0; j < NS1; ++j){
        Fu[j].resize(NS2);
        Fl[j].resize(NS2);
        Fu_n[j].resize(NS2);
        Fl_n[j].resize(NS2);
    }
    
    // *******************  Set Initial Conditions for the payoff **********************************
    for(int j=0; j<NS1; j++){
        for(int k=0; k<NS2; k++){
            double p1 = j*dS;
            double p2 = k*dS;
            Fu[j][k] = max(max(p1,p2)-E1,0.0) + max(max(p1,p2)-E2,0.0) - 2*max(max(p1,p2)-Eavg,0.0);
        }
    }
    Fl = Fu;

    // *******************  Calculate derivative price using finite difference ************************************
    for(int i=0; i<NT-1; i++){
        //Boundary Conditions
        s1 = smax1;
        s2 = smax2;
        Fu_n[0][0] = 0;
        Fu_n[NS1-1][NS2-1] = 0;
        for(int k=1; k<NS2-1;k++){
            Fu_n[0][k] = Fu[0][k] + (0.5*s2*s2*k*k*dt)*(Fu[0][k+1] - 2*Fu[0][k] + Fu[0][k-1])
            -df*Fu[0][k] + (0.5*r*k*dt)*(Fu[0][k+1] - Fu[0][k-1]);
        }
        for(int j=1; j<NS1-1; j++){
            Fu_n[j][0] = Fu[j][0] + (0.5*s1*s1*j*j*dt)*(Fu[j+1][0] - 2*Fu[j][0] + Fu[j-1][0])
            -df*Fu[j][0] + (0.5*r*j*dt)*(Fu[j+1][0] - Fu[j-1][0]);
        }
        Fl_n = Fu_n;
        
        //Finite difference
        for(int j=1; j<NS1-1; j++){
            for(int k=1; k<NS2-1; k++){
                
                double d1u = (Fu[j+1][k] - 2*Fu[j][k] + Fu[j-1][k]);
                double d2u = (Fu[j][k+1] - 2*Fu[j][k] + Fu[j][k-1]);
                double d12u = (Fu[j+1][k+1] - Fu[j+1][k-1] - Fu[j-1][k+1] + Fu[j-1][k-1])/(4);

                tie(s1,s2) = findMaxValue(d1u,d2u,d12u,smax1,smax2,smin1,smin2,corr);
                
                Fu_n[j][k] = (1-df)*Fu[j][k] + (0.5*s1*s1*j*j*dt)*d1u
                + (0.5*s2*s2*k*k*dt)*d2u
                + (0.5*r*k*dt)*(Fu[j][k+1] - Fu[j][k-1])
                + (0.5*r*j*dt)*(Fu[j+1][k] - Fu[j-1][k])
                + (0.25*corr*s1*s2*j*k*dt)*4*d12u;
                
                double d1l = (Fl[j+1][k] - 2*Fl[j][k] + Fl[j-1][k]);
                double d2l = (Fl[j][k+1] - 2*Fl[j][k] + Fl[j][k-1]);
                double d12l = (Fl[j+1][k+1] - Fl[j+1][k-1] - Fl[j-1][k+1] + Fl[j-1][k-1])/(4);

                tie(s1,s2) = findMinValue(d1l,d2l,d12l,smax1,smax2,smin1,smin2,corr);
                
                Fl_n[j][k] = (1-df)*Fl[j][k] + (0.5*s1*s1*j*j*dt)*d1l
                + (0.5*s2*s2*k*k*dt)*d2l
                + (0.5*r*k*dt)*(Fl[j][k+1] - Fl[j][k-1])
                + (0.5*r*j*dt)*(Fl[j+1][k] - Fl[j-1][k])
                + (0.25*corr*s1*s2*j*k*dt)*4*d12l;
            }
        }
        Fu = Fu_n;
        Fl = Fl_n;
    }
    
    for(int j = 1; j<NS1-1; j++){
        for(int k=1; k<NS2-1; k++){
            bsbUpper << Fu[j][k] << "\n";
            bsbLower << Fl[j][k] << "\n";
        }
        prices << j*dS << "\n";
    }
    prices.close();
    bsbUpper.close();
    bsbLower.close();
}
//Butterfly spread under Black Scholes using Finite Difference
void Bspread_BS_bounds(double r, double smax1, double smax2, double smin1, double smin2, double corr, double T, double E1, double E2, double NT, double NS, double dS){
    ofstream bsUpper("./data/2Dbs_spread_U.txt");
    ofstream bsLower("./data/2Dbs_spread_L.txt");
    ofstream prices("./data/2Dbs_spread_prices.txt");
    
    double s1 = smax1, s2 = smax2;
    double Eavg = 0.5*(E1+E2);
    int NS1 = NS, NS2 = NS;
    
    double dt = T/NT;
    double df = r*dt;

    // Initialise Price Array
    vector<vector<vector<double>>> Fu, Fu_n, Fl, Fl_n;
    Fu.resize(3);
    Fu_n.resize(3);
    Fl.resize(3);
    Fl_n.resize(3);
    for(int i=0; i<3; i++){
        Fu[i].resize(NS1);
        Fu_n[i].resize(NS1);
        Fl[i].resize(NS1);
        Fl_n[i].resize(NS1);
        for (int j = 0; j < NS1; ++j){
            Fu[i][j].resize(NS2);
            Fl[i][j].resize(NS2);
            Fu_n[i][j].resize(NS2);
            Fl_n[i][j].resize(NS2);
        }
    }
    
    // *******************  Set Initial Conditions for the payoff **********************************
    for(int j=0; j<NS1; j++){
        for(int k=0; k<NS2; k++){
            double p1 = j*dS;
            double p2 = k*dS;
            Fu[0][j][k] = max(max(p1,p2)-E1,0.0);
            Fu[1][j][k] = max(max(p1,p2)-E2,0.0);
            Fu[2][j][k] = -2*max(max(p1,p2)-Eavg,0.0);
        }
    }
    Fl = Fu;
    
    // *******************  Calculate derivative price using finite difference ************************************
    for(int i=0; i<NT-1; i++){
        //Boundary Conditions
        Fu_n[0][0][0] = 0;Fu_n[1][0][0] = 0;Fu_n[2][0][0] = 0;
        Fu_n[0][NS1-1][NS2-1] = (NS-1)*dS - E1*exp(-r*i*dt);
        Fu_n[1][NS1-1][NS2-1] = (NS-1)*dS - E2*exp(-r*i*dt);
        Fu_n[2][NS1-1][NS2-1] = -2*((NS-1)*dS - Eavg*exp(-r*i*dt));
        for(int k=1; k<NS2-1;k++){
            Fu_n[0][0][k] = Fu[0][0][k] + (0.5*s2*s2*k*k*dt)*(Fu[0][0][k+1] - 2*Fu[0][0][k] + Fu[0][0][k-1])
            -df*Fu[0][0][k] + (0.5*r*k*dt)*(Fu[0][0][k+1] - Fu[0][0][k-1]);
            Fu_n[1][0][k] = Fu[1][0][k] + (0.5*s2*s2*k*k*dt)*(Fu[1][0][k+1] - 2*Fu[1][0][k] + Fu[1][0][k-1])
            -df*Fu[1][0][k] + (0.5*r*k*dt)*(Fu[1][0][k+1] - Fu[1][0][k-1]);
            Fu_n[2][0][k] = Fu[2][0][k] + (0.5*s2*s2*k*k*dt)*(Fu[2][0][k+1] - 2*Fu[2][0][k] + Fu[2][0][k-1])
            -df*Fu[2][0][k] + (0.5*r*k*dt)*(Fu[2][0][k+1] - Fu[2][0][k-1]);
        }
        for(int j=1; j<NS1-1; j++){
            Fu_n[0][j][0] = Fu[0][j][0] + (0.5*s1*s1*j*j*dt)*(Fu[0][j+1][0] - 2*Fu[0][j][0] + Fu[0][j-1][0])
            -df*Fu[0][j][0] + (0.5*r*j*dt)*(Fu[0][j+1][0] - Fu[0][j-1][0]);
            Fu_n[1][j][0] = Fu[1][j][0] + (0.5*s1*s1*j*j*dt)*(Fu[1][j+1][0] - 2*Fu[1][j][0] + Fu[1][j-1][0])
            -df*Fu[1][j][0] + (0.5*r*j*dt)*(Fu[1][j+1][0] - Fu[1][j-1][0]);
            Fu_n[2][j][0] = Fu[2][j][0] + (0.5*s1*s1*j*j*dt)*(Fu[2][j+1][0] - 2*Fu[2][j][0] + Fu[2][j-1][0])
            -df*Fu[2][j][0] + (0.5*r*j*dt)*(Fu[2][j+1][0] - Fu[2][j-1][0]);
        }
        for(int k=0; k<NS2; k++){
            Fu_n[0][NS-1][k] = (NS-1)*dS - E1*exp(-r*i*dt);
            Fu_n[1][NS-1][k] = (NS-1)*dS - E2*exp(-r*i*dt);
            Fu_n[2][NS-1][k] = -2*((NS-1)*dS - Eavg*exp(-r*i*dt));
        }
        for(int j=0; j<NS1; j++){
            Fu_n[0][j][NS-1] = (NS-1)*dS - E1*exp(-r*i*dt);
            Fu_n[1][j][NS-1] = (NS-1)*dS - E2*exp(-r*i*dt);
            Fu_n[2][j][NS-1] = -2*((NS-1)*dS - Eavg*exp(-r*i*dt));
        }
        Fl_n = Fu_n;
        
        //Finite difference
        for(int j=1; j<NS1-1; j++){
            for(int k=1; k<NS2-1; k++){
                Fu_n[0][j][k] = (1-df)*Fu[0][j][k] + (0.5*smax1*smax1*j*j*dt)*(Fu[0][j+1][k] - 2*Fu[0][j][k] + Fu[0][j-1][k])
                + (0.5*smax2*smax2*k*k*dt)*(Fu[0][j][k+1] - 2*Fu[0][j][k] + Fu[0][j][k-1])
                + (0.5*r*k*dt)*(Fu[0][j][k+1] - Fu[0][j][k-1])
                + (0.5*r*j*dt)*(Fu[0][j+1][k] - Fu[0][j-1][k])
                + (0.25*corr*smax1*smax2*j*k*dt)*(Fu[0][j+1][k+1] - Fu[0][j+1][k-1] - Fu[0][j-1][k+1] + Fu[0][j-1][k-1]);
                
                Fu_n[1][j][k] = (1-df)*Fu[1][j][k] + (0.5*smax1*smax1*j*j*dt)*(Fu[1][j+1][k] - 2*Fu[1][j][k] + Fu[1][j-1][k])
                + (0.5*smax2*smax2*k*k*dt)*(Fu[1][j][k+1] - 2*Fu[1][j][k] + Fu[1][j][k-1])
                + (0.5*r*k*dt)*(Fu[1][j][k+1] - Fu[1][j][k-1])
                + (0.5*r*j*dt)*(Fu[1][j+1][k] - Fu[1][j-1][k])
                + (0.25*corr*smax1*smax2*j*k*dt)*(Fu[1][j+1][k+1] - Fu[1][j+1][k-1] - Fu[1][j-1][k+1] + Fu[1][j-1][k-1]);
                
                Fu_n[2][j][k] = (1-df)*Fu[2][j][k] + (0.5*smin1*smin1*j*j*dt)*(Fu[2][j+1][k] - 2*Fu[2][j][k] + Fu[2][j-1][k])
                + (0.5*smin2*smin2*k*k*dt)*(Fu[2][j][k+1] - 2*Fu[2][j][k] + Fu[2][j][k-1])
                + (0.5*r*k*dt)*(Fu[2][j][k+1] - Fu[2][j][k-1])
                + (0.5*r*j*dt)*(Fu[2][j+1][k] - Fu[2][j-1][k])
                + (0.25*corr*smin1*smin2*j*k*dt)*(Fu[2][j+1][k+1] - Fu[2][j+1][k-1] - Fu[2][j-1][k+1] + Fu[2][j-1][k-1]);
                
                Fl_n[0][j][k] = (1-df)*Fl[0][j][k] + (0.5*smin1*smin1*j*j*dt)*(Fl[0][j+1][k] - 2*Fl[0][j][k] + Fl[0][j-1][k])
                + (0.5*smin2*smin2*k*k*dt)*(Fl[0][j][k+1] - 2*Fl[0][j][k] + Fl[0][j][k-1])
                + (0.5*r*k*dt)*(Fl[0][j][k+1] - Fl[0][j][k-1])
                + (0.5*r*j*dt)*(Fl[0][j+1][k] - Fl[0][j-1][k])
                + (0.25*corr*smin1*smin2*j*k*dt)*(Fl[0][j+1][k+1] - Fl[0][j+1][k-1] - Fl[0][j-1][k+1] + Fl[0][j-1][k-1]);
                
                Fl_n[1][j][k] = (1-df)*Fl[1][j][k] + (0.5*smin1*smin1*j*j*dt)*(Fl[1][j+1][k] - 2*Fl[1][j][k] + Fl[1][j-1][k])
                + (0.5*smin2*smin2*k*k*dt)*(Fl[1][j][k+1] - 2*Fl[1][j][k] + Fl[1][j][k-1])
                + (0.5*r*k*dt)*(Fl[1][j][k+1] - Fl[1][j][k-1])
                + (0.5*r*j*dt)*(Fl[1][j+1][k] - Fl[1][j-1][k])
                + (0.25*corr*smin1*smin2*j*k*dt)*(Fl[1][j+1][k+1] - Fl[1][j+1][k-1] - Fl[1][j-1][k+1] + Fl[1][j-1][k-1]);
                
                Fl_n[2][j][k] = (1-df)*Fl[2][j][k] + (0.5*smax1*smax1*j*j*dt)*(Fl[2][j+1][k] - 2*Fl[2][j][k] + Fl[2][j-1][k])
                + (0.5*smax2*smax2*k*k*dt)*(Fl[2][j][k+1] - 2*Fl[2][j][k] + Fl[2][j][k-1])
                + (0.5*r*k*dt)*(Fl[2][j][k+1] - Fl[2][j][k-1])
                + (0.5*r*j*dt)*(Fl[2][j+1][k] - Fl[2][j-1][k])
                + (0.25*corr*smax1*smax2*j*k*dt)*(Fl[2][j+1][k+1] - Fl[2][j+1][k-1] - Fl[2][j-1][k+1] + Fl[2][j-1][k-1]);
            }
        }
        Fu = Fu_n;
        Fl = Fl_n;
    }
    
    for(int j = 1; j<NS1-1; j++){
        for(int k=1; k<NS2-1; k++){
            bsUpper << Fu[0][j][k] + Fu[1][j][k] + Fu[2][j][k] << "\n";
            bsLower << Fl[0][j][k] + Fl[1][j][k] + Fl[2][j][k]<< "\n";
        }
        prices << j*dS << "\n";
    }
    prices.close();
    bsUpper.close();
    bsLower.close();
    
}
//Butterfly spread under Black Scholes closed form solution
void Bspread_BScf_bounds(double r, double smax1, double smax2, double smin1, double smin2, double corr, double T, double E1, double E2, double NT, double NS, double dS){
    ofstream bsUpper("./data/2Dbs_spread_U.txt");
    ofstream bsLower("./data/2Dbs_spread_L.txt");
    ofstream prices("./data/2Dbs_spread_prices.txt");
    
    double Eavg = 0.5*(E1+E2);
    int NS1 = NS, NS2 = NS;
    
    for(int j = 1; j<NS1-1; j++){
        for(int k=1; k<NS2-1; k++){
            double upper = maxOptionValue(j*dS, k*dS, E1, smax1, smax2, corr, r, T) + maxOptionValue(j*dS, k*dS, E2, smax1, smax2, corr, r, T) - 2*maxOptionValue(j*dS, k*dS, Eavg, smin1, smin2, corr, r, T);
            double lower = maxOptionValue(j*dS, k*dS, E1, smin1, smin2, corr, r, T) + maxOptionValue(j*dS, k*dS, E2, smin1, smin2, corr, r, T) - 2*maxOptionValue(j*dS, k*dS, Eavg, smax1, smax2, corr, r, T);
            bsUpper << upper << "\n";
            bsLower << lower << "\n";
        }
        prices << j*dS << "\n";
    }
    prices.close();
    bsUpper.close();
    bsLower.close();
}

//Bull call spread under the Uncertain Volatility Model using Finite Difference
void Cspread_BSB_bounds(double r, double smax1, double smax2, double smin1, double smin2, double corr, double T, double E1, double E2, double NT, double NS, double dS){
    ofstream bsbUpper("./data/2Dbsb_spread_U.txt");
    ofstream bsbLower("./data/2Dbsb_spread_L.txt");
    ofstream prices("./data/2Dbsb_spread_prices.txt");
    
    double s1 = smax1, s2 = smax2;
    int NS1 = NS, NS2 = NS;
    
    double dt = T/NT;
    double df = r*dt;
    
    // Initialise Price Array
    vector<vector<double>> Fu, Fu_n, Fl, Fl_n;
    Fu.resize(NS1);
    Fu_n.resize(NS1);
    Fl.resize(NS1);
    Fl_n.resize(NS1);
    for (int j = 0; j < NS1; ++j){
        Fu[j].resize(NS2);
        Fl[j].resize(NS2);
        Fu_n[j].resize(NS2);
        Fl_n[j].resize(NS2);
    }
    
    // *******************  Set Initial Conditions for the payoff **********************************
    for(int j=0; j<NS1; j++){
        for(int k=0; k<NS2; k++){
            double p1 = j*dS;
            double p2 = k*dS;
            Fu[j][k] = max(max(p1,p2)-E1,0.0) - max(max(p1,p2)-E2,0.0);
        }
    }
    Fl = Fu;
    
    // *******************  Calculate derivative price using finite difference ************************************
    for(int i=0; i<NT-1; i++){
        //Boundary Conditions
        s1 = smax1;
        s2 = smax2;
        Fu_n[0][0] = 0;
        Fu_n[NS1-1][NS2-1] = (E2-E1)*exp(-r*(i+1)*dt);
        for(int k=1; k<NS2-1;k++){
            Fu_n[0][k] = Fu[0][k] + (0.5*s2*s2*k*k*dt)*(Fu[0][k+1] - 2*Fu[0][k] + Fu[0][k-1])
            -df*Fu[0][k] + (0.5*r*k*dt)*(Fu[0][k+1] - Fu[0][k-1]);
        }
        for(int j=1; j<NS1-1; j++){
            Fu_n[j][0] = Fu[j][0] + (0.5*s1*s1*j*j*dt)*(Fu[j+1][0] - 2*Fu[j][0] + Fu[j-1][0])
            -df*Fu[j][0] + (0.5*r*j*dt)*(Fu[j+1][0] - Fu[j-1][0]);
        }
        for(int k=0; k<NS2; k++){
            Fu_n[NS-1][k] = (E2-E1)*exp(-r*(i+1)*dt);
        }
        for(int j=0; j<NS1; j++){
            Fu_n[j][NS-1] = (E2-E1)*exp(-r*(i+1)*dt);
         }
        Fl_n = Fu_n;
        
        //Finite difference
        for(int j=1; j<NS1-1; j++){
            for(int k=1; k<NS2-1; k++){
                
                double d1u = (Fu[j+1][k] - 2*Fu[j][k] + Fu[j-1][k]);
                double d2u = (Fu[j][k+1] - 2*Fu[j][k] + Fu[j][k-1]);
                double d12u = (Fu[j+1][k+1] - Fu[j+1][k-1] - Fu[j-1][k+1] + Fu[j-1][k-1])/(4);
                
                tie(s1,s2) = findMaxValue(d1u,d2u,d12u,smax1,smax2,smin1,smin2,corr);
                
                Fu_n[j][k] = (1-df)*Fu[j][k] + (0.5*s1*s1*j*j*dt)*d1u
                + (0.5*s2*s2*k*k*dt)*d2u
                + (0.5*r*k*dt)*(Fu[j][k+1] - Fu[j][k-1])
                + (0.5*r*j*dt)*(Fu[j+1][k] - Fu[j-1][k])
                + (0.25*corr*s1*s2*j*k*dt)*4*d12u;
                
                double d1l = (Fl[j+1][k] - 2*Fl[j][k] + Fl[j-1][k]);
                double d2l = (Fl[j][k+1] - 2*Fl[j][k] + Fl[j][k-1]);
                double d12l = (Fl[j+1][k+1] - Fl[j+1][k-1] - Fl[j-1][k+1] + Fl[j-1][k-1])/(4);
                
                tie(s1,s2) = findMinValue(d1l,d2l,d12l,smax1,smax2,smin1,smin2,corr);
                
                Fl_n[j][k] = (1-df)*Fl[j][k] + (0.5*s1*s1*j*j*dt)*d1l
                + (0.5*s2*s2*k*k*dt)*d2l
                + (0.5*r*k*dt)*(Fl[j][k+1] - Fl[j][k-1])
                + (0.5*r*j*dt)*(Fl[j+1][k] - Fl[j-1][k])
                + (0.25*corr*s1*s2*j*k*dt)*4*d12l;
            }
        }
        Fu = Fu_n;
        Fl = Fl_n;
    }
    
    cout<<"40 40"<<" "<<Fl[30][30]<<"\t"<<Fu[30][30]<<"\n";
    
    for(int j = 1; j<NS1-1; j++){
        for(int k=1; k<NS2-1; k++){
            bsbUpper << Fu[j][k] << "\n";
            bsbLower << Fl[j][k] << "\n";
        }
        prices << j*dS << "\n";
    }
    prices.close();
    bsbUpper.close();
    bsbLower.close();
}
//Bull call spread under Black Scholes using Finite Difference
void Cspread_BS_bounds(double r, double smax1, double smax2, double smin1, double smin2, double corr, double T, double E1, double E2, double NT, double NS, double dS){
    ofstream bsUpper("./data/2Dbs_spread_U.txt");
    ofstream bsLower("./data/2Dbs_spread_L.txt");
    ofstream prices("./data/2Dbs_spread_prices.txt");
    
    double s1 = smax1, s2 = smax2;
    int NS1 = NS, NS2 = NS;
    
    double dt = T/NT;
    double df = r*dt;
    
    // Initialise Price Array
    vector<vector<vector<double>>> Fu, Fu_n, Fl, Fl_n;
    Fu.resize(2);
    Fu_n.resize(2);
    Fl.resize(2);
    Fl_n.resize(2);
    for(int i=0; i<2; i++){
        Fu[i].resize(NS1);
        Fu_n[i].resize(NS1);
        Fl[i].resize(NS1);
        Fl_n[i].resize(NS1);
        for (int j = 0; j < NS1; ++j){
            Fu[i][j].resize(NS2);
            Fl[i][j].resize(NS2);
            Fu_n[i][j].resize(NS2);
            Fl_n[i][j].resize(NS2);
        }
    }
    
    // *******************  Set Initial Conditions for the payoff **********************************
    for(int j=0; j<NS1; j++){
        for(int k=0; k<NS2; k++){
            double p1 = j*dS;
            double p2 = k*dS;
            Fu[0][j][k] = max(max(p1,p2)-E1,0.0);
            Fu[1][j][k] = -max(max(p1,p2)-E2,0.0);
        }
    }
    Fl = Fu;
    
    // *******************  Calculate derivative price using finite difference ************************************
    for(int i=0; i<NT-1; i++){
        //Boundary Conditions
        Fu_n[0][0][0] = 0;Fu_n[1][0][0] = 0;
        Fu_n[0][NS1-1][NS2-1] = (NS-1)*dS - E1*exp(-r*(i+1)*dt);
        Fu_n[1][NS1-1][NS2-1] = -(NS-1)*dS + E2*exp(-r*(i+1)*dt);
        for(int k=1; k<NS2-1;k++){
            Fu_n[0][0][k] = Fu[0][0][k] + (0.5*s2*s2*k*k*dt)*(Fu[0][0][k+1] - 2*Fu[0][0][k] + Fu[0][0][k-1])
            -df*Fu[0][0][k] + (0.5*r*k*dt)*(Fu[0][0][k+1] - Fu[0][0][k-1]);
            Fu_n[1][0][k] = Fu[1][0][k] + (0.5*s2*s2*k*k*dt)*(Fu[1][0][k+1] - 2*Fu[1][0][k] + Fu[1][0][k-1])
            -df*Fu[1][0][k] + (0.5*r*k*dt)*(Fu[1][0][k+1] - Fu[1][0][k-1]);
        }
        for(int j=1; j<NS1-1; j++){
            Fu_n[0][j][0] = Fu[0][j][0] + (0.5*s1*s1*j*j*dt)*(Fu[0][j+1][0] - 2*Fu[0][j][0] + Fu[0][j-1][0])
            -df*Fu[0][j][0] + (0.5*r*j*dt)*(Fu[0][j+1][0] - Fu[0][j-1][0]);
            Fu_n[1][j][0] = Fu[1][j][0] + (0.5*s1*s1*j*j*dt)*(Fu[1][j+1][0] - 2*Fu[1][j][0] + Fu[1][j-1][0])
            -df*Fu[1][j][0] + (0.5*r*j*dt)*(Fu[1][j+1][0] - Fu[1][j-1][0]);
        }
        for(int k=0; k<NS2; k++){
            Fu_n[0][NS-1][k] = (NS-1)*dS - E1*exp(-r*(i+1)*dt);
            Fu_n[1][NS-1][k] = -(NS-1)*dS + E2*exp(-r*(i+1)*dt);
        }
        for(int j=0; j<NS1; j++){
            Fu_n[0][j][NS-1] = (NS-1)*dS - E1*exp(-r*(i+1)*dt);
            Fu_n[1][j][NS-1] = -(NS-1)*dS + E2*exp(-r*(i+1)*dt);
        }
        Fl_n = Fu_n;
        
        //Finite difference
        for(int j=1; j<NS1-1; j++){
            for(int k=1; k<NS2-1; k++){
                Fu_n[0][j][k] = (1-df)*Fu[0][j][k] + (0.5*smax1*smax1*j*j*dt)*(Fu[0][j+1][k] - 2*Fu[0][j][k] + Fu[0][j-1][k])
                + (0.5*smax2*smax2*k*k*dt)*(Fu[0][j][k+1] - 2*Fu[0][j][k] + Fu[0][j][k-1])
                + (0.5*r*k*dt)*(Fu[0][j][k+1] - Fu[0][j][k-1])
                + (0.5*r*j*dt)*(Fu[0][j+1][k] - Fu[0][j-1][k])
                + (0.25*corr*smax1*smax2*j*k*dt)*(Fu[0][j+1][k+1] - Fu[0][j+1][k-1] - Fu[0][j-1][k+1] + Fu[0][j-1][k-1]);
                
                Fu_n[1][j][k] = (1-df)*Fu[1][j][k] + (0.5*smin1*smin1*j*j*dt)*(Fu[1][j+1][k] - 2*Fu[1][j][k] + Fu[1][j-1][k])
                + (0.5*smin2*smin2*k*k*dt)*(Fu[1][j][k+1] - 2*Fu[1][j][k] + Fu[1][j][k-1])
                + (0.5*r*k*dt)*(Fu[1][j][k+1] - Fu[1][j][k-1])
                + (0.5*r*j*dt)*(Fu[1][j+1][k] - Fu[1][j-1][k])
                + (0.25*corr*smin1*smin2*j*k*dt)*(Fu[1][j+1][k+1] - Fu[1][j+1][k-1] - Fu[1][j-1][k+1] + Fu[1][j-1][k-1]);
                
                Fl_n[0][j][k] = (1-df)*Fl[0][j][k] + (0.5*smin1*smin1*j*j*dt)*(Fl[0][j+1][k] - 2*Fl[0][j][k] + Fl[0][j-1][k])
                + (0.5*smin2*smin2*k*k*dt)*(Fl[0][j][k+1] - 2*Fl[0][j][k] + Fl[0][j][k-1])
                + (0.5*r*k*dt)*(Fl[0][j][k+1] - Fl[0][j][k-1])
                + (0.5*r*j*dt)*(Fl[0][j+1][k] - Fl[0][j-1][k])
                + (0.25*corr*smin1*smin2*j*k*dt)*(Fl[0][j+1][k+1] - Fl[0][j+1][k-1] - Fl[0][j-1][k+1] + Fl[0][j-1][k-1]);
                
                Fl_n[1][j][k] = (1-df)*Fl[1][j][k] + (0.5*smax1*smax1*j*j*dt)*(Fl[1][j+1][k] - 2*Fl[1][j][k] + Fl[1][j-1][k])
                + (0.5*smax2*smax2*k*k*dt)*(Fl[1][j][k+1] - 2*Fl[1][j][k] + Fl[1][j][k-1])
                + (0.5*r*k*dt)*(Fl[1][j][k+1] - Fl[1][j][k-1])
                + (0.5*r*j*dt)*(Fl[1][j+1][k] - Fl[1][j-1][k])
                + (0.25*corr*smax1*smax2*j*k*dt)*(Fl[1][j+1][k+1] - Fl[1][j+1][k-1] - Fl[1][j-1][k+1] + Fl[1][j-1][k-1]);
            }
        }
        Fu = Fu_n;
        Fl = Fl_n;
    }
    
    for(int j = 1; j<NS1-1; j++){
        for(int k=1; k<NS2-1; k++){
            bsUpper << Fu[0][j][k] + Fu[1][j][k] << "\n";
            bsLower << Fl[0][j][k] + Fl[1][j][k]<< "\n";
        }
        prices << j*dS << "\n";
    }
    prices.close();
    bsUpper.close();
    bsLower.close();
}
//Bull call spread under Black Scholes closed form solution
void Cspread_BScf_bounds(double r, double smax1, double smax2, double smin1, double smin2, double corr, double T, double E1, double E2, double NT, double NS, double dS){
    ofstream bsUpper("./data/2Dbs_spread_U.txt");
    ofstream bsLower("./data/2Dbs_spread_L.txt");
    ofstream prices("./data/2Dbs_spread_prices.txt");
    
    int NS1 = NS, NS2 = NS;
    
    for(int j = 1; j<NS1-1; j++){
        for(int k=1; k<NS2-1; k++){
            double upper = maxOptionValue(j*dS, k*dS, E1, smax1, smax2, corr, r, T) - maxOptionValue(j*dS, k*dS, E2, smin1, smin2, corr, r, T);
            double lower = maxOptionValue(j*dS, k*dS, E1, smin1, smin2, corr, r, T) - maxOptionValue(j*dS, k*dS, E2, smax1, smax2, corr, r, T);
            bsUpper << upper << "\n";
            bsLower << lower << "\n";
        }
        prices << j*dS << "\n";
    }
    prices.close();
    bsUpper.close();
    bsLower.close();
}

//Dual currency spread under Uncertain Volatility Model using Finite Difference
void vanilla_BSB_DC(double rf, double rd, double smax, double smin, double corr, double T, double E1, double E2, int NT, int NS, int nOfOpts, int DCorPRDC, double dS){
    ofstream bsbUpper("./data/2Dbsb_spread_U.txt");
    ofstream bsbLower("./data/2Dbsb_spread_L.txt");
    ofstream prices("./data/2Dbsb_spread_prices.txt");
    
    //*********** Variable definitions *************************//
    double r = rf-rd;
    double s1 = smax, s2 = smax;
    double smax1 = smax, smax2 = smax, smin1 = smin, smin2 = smin;
    
    double FX0 = E1, FX1 = E2;
    
    double FXstrike0 = FX0*rd/rf;
    double L0 = rf/FX0;
    double FXstrike1 = FX1*rd/rf;
    double L1 = rf/FX1;
    if(DCorPRDC==1){
        L0=1; L1=1;
    }
    
    int NTstep = NT/nOfOpts;
    int NS1 = NS, NS2 = NS;
    
    double dt = T/NT;
    double df = r*dt;
    
    double FXmax = dS * (NS-1);
    
    double boundTime = 0.0, findifTime = 0.0;
    
    clock_t begin1 = clock();
    // Initialise Price Array
    vector<vector<double>> Fu, Fu_n, Fl, Fl_n;
    Fu.resize(NS1);
    Fu_n.resize(NS1);
    Fl.resize(NS1);
    Fl_n.resize(NS1);
    for (int j = 0; j < NS1; ++j){
        Fu[j].resize(NS2);
        Fl[j].resize(NS2);
        Fu_n[j].resize(NS2);
        Fl_n[j].resize(NS2);
    }
    
    // *******************  Set Initial Conditions for the payoff **********************************
    for(int j=0; j<NS1; j++){
        for(int k=0; k<NS2; k++){
            double p1 = j*dS;
            double p2 = k*dS;
            Fu[j][k] = L0*max(max(p1,p2)-FXstrike0,0.0) - L1*max(max(p1,p2)-FXstrike1,0.0);
        }
    }
    Fl = Fu;
    clock_t end1 = clock();
    double initTime = (double)(end1 - begin1);
    // *******************  Calculate derivative price using finite difference ************************************
    for(int i=0; i<NT-1; i++){
        clock_t begin2 = clock();
        //Boundary Conditions
        if((i+1)%NTstep==0){
            for(int j=0; j<NS1; j++){
                for(int k=0; k<NS2; k++){
                    double p1 = j*dS;
                    double p2 = k*dS;
                    Fu_n[j][k] = L0*max(max(p1,p2)-FXstrike0,0.0) - L1*max(max(p1,p2)-FXstrike1,0.0);
                }
            }
        }else{
            s1 = smax;
            s2 = smax;
            Fu_n[0][0] = 0;
            Fu_n[NS1-1][NS2-1] = (L0-L1)*FXmax + (FXstrike1-FXstrike0)*exp(-r*(i+1)*dt);
            for(int k=1; k<NS2-1;k++){
                Fu_n[0][k] = Fu[0][k] + (0.5*s2*s2*k*k*dt)*(Fu[0][k+1] - 2*Fu[0][k] + Fu[0][k-1])
                -df*Fu[0][k] + (0.5*r*k*dt)*(Fu[0][k+1] - Fu[0][k-1]);
            }
            for(int j=1; j<NS1-1; j++){
                Fu_n[j][0] = Fu[j][0] + (0.5*s1*s1*j*j*dt)*(Fu[j+1][0] - 2*Fu[j][0] + Fu[j-1][0])
                -df*Fu[j][0] + (0.5*r*j*dt)*(Fu[j+1][0] - Fu[j-1][0]);
            }
            for(int k=0; k<NS2; k++){
                Fu_n[NS-1][k] = (L0-L1)*FXmax + (FXstrike1-FXstrike0)*exp(-r*(i+1)*dt);
            }
            for(int j=0; j<NS1; j++){
                Fu_n[j][NS-1] = (L0-L1)*FXmax + (FXstrike1-FXstrike0)*exp(-r*(i+1)*dt);
            }
        }
        Fl_n = Fu_n;
        clock_t end2 = clock();
        boundTime += (double)(end2 - begin2);
        
        clock_t begin3 = clock();
        //Finite difference
        for(int j=1; j<NS1-1; j++){
            for(int k=1; k<NS2-1; k++){
                
                double d1u = (Fu[j+1][k] - 2*Fu[j][k] + Fu[j-1][k]);
                double d2u = (Fu[j][k+1] - 2*Fu[j][k] + Fu[j][k-1]);
                double d12u = (Fu[j+1][k+1] - Fu[j+1][k-1] - Fu[j-1][k+1] + Fu[j-1][k-1])/(4);
                
                tie(s1,s2) = findMaxValue(d1u,d2u,d12u,smax1,smax2,smin1,smin2,corr);
                
                Fu_n[j][k] = (1-df)*Fu[j][k] + (0.5*s1*s1*j*j*dt)*d1u
                + (0.5*s2*s2*k*k*dt)*d2u
                + (0.5*r*k*dt)*(Fu[j][k+1] - Fu[j][k-1])
                + (0.5*r*j*dt)*(Fu[j+1][k] - Fu[j-1][k])
                + (0.25*corr*s1*s2*j*k*dt)*4*d12u;
                
                double d1l = (Fl[j+1][k] - 2*Fl[j][k] + Fl[j-1][k]);
                double d2l = (Fl[j][k+1] - 2*Fl[j][k] + Fl[j][k-1]);
                double d12l = (Fl[j+1][k+1] - Fl[j+1][k-1] - Fl[j-1][k+1] + Fl[j-1][k-1])/(4);
                
                tie(s1,s2) = findMinValue(d1l,d2l,d12l,smax1,smax2,smin1,smin2,corr);
                
                Fl_n[j][k] = (1-df)*Fl[j][k] + (0.5*s1*s1*j*j*dt)*d1l
                + (0.5*s2*s2*k*k*dt)*d2l
                + (0.5*r*k*dt)*(Fl[j][k+1] - Fl[j][k-1])
                + (0.5*r*j*dt)*(Fl[j+1][k] - Fl[j-1][k])
                + (0.25*corr*s1*s2*j*k*dt)*4*d12l;
            }
        }
        Fu = Fu_n;
        Fl = Fl_n;
        clock_t end3 = clock();
        findifTime += (double)(end3 - begin3);
    }
    
    cout << "Init time: "<< initTime/ CLOCKS_PER_SEC << "\n";
    cout << "Boundary time: "<< boundTime/ CLOCKS_PER_SEC << "\n";
    cout << "Finite diff time: "<< findifTime/ CLOCKS_PER_SEC << "\n";

    for(int j = 1; j<NS1-1; j++){
        for(int k=1; k<NS2-1; k++){
            bsbUpper << Fu[j][k] << "\n";
            bsbLower << Fl[j][k] << "\n";
        }
        prices << j*dS << "\n";
    }
    prices.close();
    bsbUpper.close();
    bsbLower.close();
}
//Duall currency spread under Black Scholes using Finite Difference
void vanilla_BS_DC(double rf, double rd, double smax, double smin, double corr, double T, double E1, double E2, int NT, int NS, int nOfOpts, int DCorPRDC, double dS){
    ofstream bsUpper("./data/2Dbs_spread_U.txt");
    ofstream bsLower("./data/2Dbs_spread_L.txt");
    ofstream prices("./data/2Dbs_spread_prices.txt");
    //*********** Variable definitions *************************/
    double r = rf-rd;
    double s1 = smax, s2 = smax;
    double smax1 = smax, smax2 = smax, smin1 = smin, smin2 = smin;
    
    double FX0 = E1, FX1 = E2;
    
    double FXstrike0 = FX0*rd/rf;
    double L0 = rf/FX0;
    double FXstrike1 = FX1*rd/rf;
    double L1 = rf/FX1;
    if(DCorPRDC==1){
        L0=1; L1=1;
    }
    int NTstep = NT/nOfOpts;
    int NS1 = NS, NS2 = NS;
    
    double dt = T/NT;
    double df = r*dt;
    
    double FXmax = dS * (NS-1);
    
    // Initialise Price Array
    vector<vector<vector<double>>> Fu, Fu_n, Fl, Fl_n;
    Fu.resize(2);
    Fu_n.resize(2);
    Fl.resize(2);
    Fl_n.resize(2);
    for(int i=0; i<2; i++){
        Fu[i].resize(NS1);
        Fu_n[i].resize(NS1);
        Fl[i].resize(NS1);
        Fl_n[i].resize(NS1);
        for (int j = 0; j < NS1; ++j){
            Fu[i][j].resize(NS2);
            Fl[i][j].resize(NS2);
            Fu_n[i][j].resize(NS2);
            Fl_n[i][j].resize(NS2);
        }
    }
    // *******************  Set Initial Conditions for the payoff **********************************
    for(int j=0; j<NS1; j++){
        for(int k=0; k<NS2; k++){
            double p1 = j*dS;
            double p2 = k*dS;
            Fu[0][j][k] = L0*max(max(p1,p2)-FXstrike0,0.0);
            Fu[1][j][k] = -L1*max(max(p1,p2)-FXstrike1,0.0);
        }
    }
    Fl = Fu;
    
    // *******************  Calculate derivative price using finite difference ************************************
    for(int i=0; i<NT-1; i++){
        //Boundary Conditions
        if((i+1)%NTstep==0){
            for(int j=0; j<NS1; j++){
                for(int k=0; k<NS2; k++){
                    double p1 = j*dS;
                    double p2 = k*dS;
                    Fu_n[0][j][k] = L0*max(max(p1,p2)-FXstrike0,0.0);
                    Fu_n[1][j][k] = -L1*max(max(p1,p2)-FXstrike1,0.0);
                }
            }
        }else{
            s1 = smax;
            s2 = smax;
            Fu_n[0][0][0] = 0;Fu_n[1][0][0] = 0;
            Fu_n[0][NS1-1][NS2-1] = L0*FXmax - FXstrike0*exp(-r*(i+1)*dt);
            Fu_n[1][NS1-1][NS2-1] = -L1*FXmax + FXstrike1*exp(-r*(i+1)*dt);
            for(int k=1; k<NS2-1;k++){
                Fu_n[0][0][k] = Fu[0][0][k] + (0.5*s2*s2*k*k*dt)*(Fu[0][0][k+1] - 2*Fu[0][0][k] + Fu[0][0][k-1])
                -df*Fu[0][0][k] + (0.5*r*k*dt)*(Fu[0][0][k+1] - Fu[0][0][k-1]);
                Fu_n[1][0][k] = Fu[1][0][k] + (0.5*s2*s2*k*k*dt)*(Fu[1][0][k+1] - 2*Fu[1][0][k] + Fu[1][0][k-1])
                -df*Fu[1][0][k] + (0.5*r*k*dt)*(Fu[1][0][k+1] - Fu[1][0][k-1]);
            }
            for(int j=1; j<NS1-1; j++){
                Fu_n[0][j][0] = Fu[0][j][0] + (0.5*s1*s1*j*j*dt)*(Fu[0][j+1][0] - 2*Fu[0][j][0] + Fu[0][j-1][0])
                -df*Fu[0][j][0] + (0.5*r*j*dt)*(Fu[0][j+1][0] - Fu[0][j-1][0]);
                Fu_n[1][j][0] = Fu[1][j][0] + (0.5*s1*s1*j*j*dt)*(Fu[1][j+1][0] - 2*Fu[1][j][0] + Fu[1][j-1][0])
                -df*Fu[1][j][0] + (0.5*r*j*dt)*(Fu[1][j+1][0] - Fu[1][j-1][0]);
            }
            for(int k=0; k<NS2; k++){
                Fu_n[0][NS-1][k] = L0*FXmax - FXstrike0*exp(-r*(i+1)*dt);
                Fu_n[1][NS-1][k] = -L1*FXmax + FXstrike1*exp(-r*(i+1)*dt);
            }
            for(int j=0; j<NS1; j++){
                Fu_n[0][j][NS-1] = L0*FXmax - FXstrike0*exp(-r*(i+1)*dt);
                Fu_n[1][j][NS-1] = -L1*FXmax + FXstrike1*exp(-r*(i+1)*dt);
            }
        }
        Fl_n = Fu_n;
        
        //Finite difference
        for(int j=1; j<NS1-1; j++){
            for(int k=1; k<NS2-1; k++){
                Fu_n[0][j][k] = (1-df)*Fu[0][j][k] + (0.5*smax1*smax1*j*j*dt)*(Fu[0][j+1][k] - 2*Fu[0][j][k] + Fu[0][j-1][k])
                + (0.5*smax2*smax2*k*k*dt)*(Fu[0][j][k+1] - 2*Fu[0][j][k] + Fu[0][j][k-1])
                + (0.5*r*k*dt)*(Fu[0][j][k+1] - Fu[0][j][k-1])
                + (0.5*r*j*dt)*(Fu[0][j+1][k] - Fu[0][j-1][k])
                + (0.25*corr*smax1*smax2*j*k*dt)*(Fu[0][j+1][k+1] - Fu[0][j+1][k-1] - Fu[0][j-1][k+1] + Fu[0][j-1][k-1]);
                
                Fu_n[1][j][k] = (1-df)*Fu[1][j][k] + (0.5*smin1*smin1*j*j*dt)*(Fu[1][j+1][k] - 2*Fu[1][j][k] + Fu[1][j-1][k])
                + (0.5*smin2*smin2*k*k*dt)*(Fu[1][j][k+1] - 2*Fu[1][j][k] + Fu[1][j][k-1])
                + (0.5*r*k*dt)*(Fu[1][j][k+1] - Fu[1][j][k-1])
                + (0.5*r*j*dt)*(Fu[1][j+1][k] - Fu[1][j-1][k])
                + (0.25*corr*smin1*smin2*j*k*dt)*(Fu[1][j+1][k+1] - Fu[1][j+1][k-1] - Fu[1][j-1][k+1] + Fu[1][j-1][k-1]);
                
                Fl_n[0][j][k] = (1-df)*Fl[0][j][k] + (0.5*smin1*smin1*j*j*dt)*(Fl[0][j+1][k] - 2*Fl[0][j][k] + Fl[0][j-1][k])
                + (0.5*smin2*smin2*k*k*dt)*(Fl[0][j][k+1] - 2*Fl[0][j][k] + Fl[0][j][k-1])
                + (0.5*r*k*dt)*(Fl[0][j][k+1] - Fl[0][j][k-1])
                + (0.5*r*j*dt)*(Fl[0][j+1][k] - Fl[0][j-1][k])
                + (0.25*corr*smin1*smin2*j*k*dt)*(Fl[0][j+1][k+1] - Fl[0][j+1][k-1] - Fl[0][j-1][k+1] + Fl[0][j-1][k-1]);
                
                Fl_n[1][j][k] = (1-df)*Fl[1][j][k] + (0.5*smax1*smax1*j*j*dt)*(Fl[1][j+1][k] - 2*Fl[1][j][k] + Fl[1][j-1][k])
                + (0.5*smax2*smax2*k*k*dt)*(Fl[1][j][k+1] - 2*Fl[1][j][k] + Fl[1][j][k-1])
                + (0.5*r*k*dt)*(Fl[1][j][k+1] - Fl[1][j][k-1])
                + (0.5*r*j*dt)*(Fl[1][j+1][k] - Fl[1][j-1][k])
                + (0.25*corr*smax1*smax2*j*k*dt)*(Fl[1][j+1][k+1] - Fl[1][j+1][k-1] - Fl[1][j-1][k+1] + Fl[1][j-1][k-1]);
            }
        }
        Fu = Fu_n;
        Fl = Fl_n;
    }
    for(int j = 1; j<NS1-1; j++){
        for(int k=1; k<NS2-1; k++){
            bsUpper << Fu[0][j][k] + Fu[1][j][k] << "\n";
            bsLower << Fl[0][j][k] + Fl[1][j][k]<< "\n";
        }
        prices << j*dS << "\n";
    }
    prices.close();
    bsUpper.close();
    bsLower.close();
}

void twoStateModelling(){
    
    double r = 0.05, corr = 0.5;
    double smax1 = 0.5, smax2 = 0.5, smin1 = 0.1, smin2 = 0.1;
    double T = 0.5;
    double E1 = 35, E2 = 45;
    int NT = 1000, NS = 50; double dS = 2;  //need NT = 20*NS

    
    int i = 3; //1 for bull call vanilla, 2 for butterfly vanilla, 3 for bull call prdc
    if(i == 1){
        Cspread_BSB_bounds(r, smax1, smax2, smin1, smin2, corr, T, E1, E2, NT, NS, dS);
        Cspread_BS_bounds(r, smax1, smax2, smin1, smin2, corr, T, E1, E2, NT, NS, dS);
    }else if(i==2){
        Bspread_BSB_bounds(r, smax1, smax2, smin1, smin2, corr, T, E1, E2, NT, NS, dS);
        Bspread_BS_bounds(r, smax1, smax2, smin1, smin2, corr, T, E1, E2, NT, NS, dS);
    }else if(i==3){
        int nOfOpts = 5;
        double rf = 0.15, rd = 0.10;
        E1 = 45, E2 = 55;
        NT = NT*nOfOpts;
        T = T*nOfOpts;
        int DCorPRDC = 2; //1fordc 2forprdc
        vanilla_BS_DC(rf, rd, smax1, smin1, corr, T, E1, E2, NT, NS, nOfOpts, DCorPRDC, dS);
        vanilla_BSB_DC(rf, rd, smax1, smin1, corr, T, E1, E2, NT, NS, nOfOpts, DCorPRDC, dS);
    }
    //system("sh scripts/uncVol2D_script.sh");
}