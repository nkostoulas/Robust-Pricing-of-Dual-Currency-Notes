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
double valueToOpt(double d1, double d2, double d12, double s1, double s2){
    return d1*s1*s1 + 2*d12*s1*s2 + d2*s2*s2;
}
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
void twoStateBSBVanilla(){
    ofstream bsbUpper("./data/2Dbsb_spread_U.txt");
    ofstream bsbLower("./data/2Dbsb_spread_L.txt");
    ofstream prices("./data/2Dbsb_spread_prices.txt");
    double r = 0.05;
    double smax1 = 0.5, smax2 = 0.5, smin1 = 0.1, smin2 = 0.1;
    double s1 = smax1, s2 = smax2;
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
    
    for(int i=0; i<NT-1; i++){
        for(int j=1; j<NS1-1; j++){
            for(int k=1; k<NS2-1; k++){
                
                double d1u = (Fu[i][j+1][k] - 2*Fu[i][j][k] + Fu[i][j-1][k]);
                double d2u = (Fu[i][j][k+1] - 2*Fu[i][j][k] + Fu[i][j][k-1]);
                double d12u = (Fu[i][j+1][k+1] - Fu[i][j+1][k-1] - Fu[i][j-1][k+1] + Fu[i][j-1][k-1])/(4);
                if(d12u>0){
                    corr = 0.5;
                }else{
                    corr = 0.3;
                }
                tie(s1,s2) = findMaxValue(d1u,d2u,d12u,smax1,smax2,smin1,smin2,corr);
                
                Fu[i+1][j][k] = (1-df)*Fu[i][j][k] + (0.5*s1*s1*j*j*dt)*d1u
                + (0.5*s2*s2*k*k*dt)*d2u
                + (0.5*r*k*dt)*(Fu[i][j][k+1] - Fu[i][j][k-1])
                + (0.5*r*j*dt)*(Fu[i][j+1][k] - Fu[i][j-1][k])
                + (0.25*corr*s1*s2*j*k*dt)*4*d12u;
                
                double d1l = (Fl[i][j+1][k] - 2*Fl[i][j][k] + Fl[i][j-1][k]);
                double d2l = (Fl[i][j][k+1] - 2*Fl[i][j][k] + Fl[i][j][k-1]);
                double d12l = (Fl[i][j+1][k+1] - Fl[i][j+1][k-1] - Fl[i][j-1][k+1] + Fl[i][j-1][k-1])/(4);
                if(d12u>0){
                    corr = 0.3;
                }else{
                    corr = 0.5;
                }
                tie(s1,s2) = findMinValue(d1l,d2l,d12l,smax1,smax2,smin1,smin2,corr);
                
                Fl[i+1][j][k] = (1-df)*Fl[i][j][k] + (0.5*s1*s1*j*j*dt)*d1l
                + (0.5*s2*s2*k*k*dt)*d2l
                + (0.5*r*k*dt)*(Fl[i][j][k+1] - Fl[i][j][k-1])
                + (0.5*r*j*dt)*(Fl[i][j+1][k] - Fl[i][j-1][k])
                + (0.25*corr*s1*s2*j*k*dt)*4*d12l;
                
                double det = d1u*d2u - d12u*d12u;
                if(det==0 && d1u!=0 && d2u!=0 && d12u!=0){
                    cnt+=1;
                }
            }
        }
    }
    cout << cnt << "\n";

    cout<<"40 40"<<" "<<Fl[NT-1][40][40]<<"\t"<<Fu[NT-1][40][40]<<"\n";
    
    for(int j = 1; j<NS1-1; j++){
        for(int k=1; k<NS2-1; k++){
            bsbUpper << Fu[NT-1][j][k] << "\n";
            bsbLower << Fl[NT-1][j][k] << "\n";
        }
        prices << j*dS << "\n";
    }
    prices.close();
    bsbUpper.close();
    bsbLower.close();
}
void twoStateBSBSpread(){
    ofstream bsbUpper("./data/2Dbsb_spread_U.txt");
    ofstream bsbLower("./data/2Dbsb_spread_L.txt");
    ofstream prices("./data/2Dbsb_spread_prices.txt");
    double r = 0.05;
    double smax1 = 0.5, smax2 = 0.5, smin1 = 0.3, smin2 = 0.3;
    double s1 = smax1, s2 = smax2;
    double corr = 0.3;
    
    double T = 0.5;
    double E1 = 35, E2 = 45, Eavg = 40;
    
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
            Fu[0][j][k] = max(max(p1,p2)-E1,0.0) + max(max(p1,p2)-E2,0.0) - 2*max(max(p1,p2)-Eavg,0.0);
        }
    }
    //Boundary
    for(int i=1; i<NT; i++){
        Fu[i][0][0] = 0;
        Fu[i][NS1-1][NS2-1] = 0;
        for(int k=1; k<NS2-1;k++){
            Fu[i][0][k] = Fu[i-1][0][k] + (0.5*s2*s2*k*k*dt)*(Fu[i-1][0][k+1] - 2*Fu[i-1][0][k] + Fu[i-1][0][k-1])
            -df*Fu[i-1][0][k] + (0.5*r*k*dt)*(Fu[i-1][0][k+1] - Fu[i-1][0][k-1]);
        }
        for(int j=1; j<NS1-1; j++){
            Fu[i][j][0] = Fu[i-1][j][0] + (0.5*s1*s1*j*j*dt)*(Fu[i-1][j+1][0] - 2*Fu[i-1][j][0] + Fu[i-1][j-1][0])
            -df*Fu[i-1][j][0] + (0.5*r*j*dt)*(Fu[i-1][j+1][0] - Fu[i-1][j-1][0]);
        }

    }
    Fl = Fu;/*
    for(int i=0; i<5; i++){
    for(int j =0; j<NS1; j++){
        for(int k=0; k<NS2; k++){
            cout<<Fu[i][j][k]<<"\t";
        }
        cout<<"\n";
    }
    cout<<"\n";
    }*/
    // *******************  Calculate derivative price using finite difference ************************************
    int cnt = 0;
    
    for(int i=0; i<NT-1; i++){
        for(int j=1; j<NS1-1; j++){
            for(int k=1; k<NS2-1; k++){
                
                double d1u = (Fu[i][j+1][k] - 2*Fu[i][j][k] + Fu[i][j-1][k]);
                double d2u = (Fu[i][j][k+1] - 2*Fu[i][j][k] + Fu[i][j][k-1]);
                double d12u = (Fu[i][j+1][k+1] - Fu[i][j+1][k-1] - Fu[i][j-1][k+1] + Fu[i][j-1][k-1])/(4);
                if(d12u>0){
                    corr = 0.5;
                }else{
                    corr = 0.3;
                }
                tie(s1,s2) = findMaxValue(d1u,d2u,d12u,smax1,smax2,smin1,smin2,corr);
                
                Fu[i+1][j][k] = (1-df)*Fu[i][j][k] + (0.5*s1*s1*j*j*dt)*d1u
                + (0.5*s2*s2*k*k*dt)*d2u
                + (0.5*r*k*dt)*(Fu[i][j][k+1] - Fu[i][j][k-1])
                + (0.5*r*j*dt)*(Fu[i][j+1][k] - Fu[i][j-1][k])
                + (0.25*corr*s1*s2*j*k*dt)*4*d12u;
                
                double d1l = (Fl[i][j+1][k] - 2*Fl[i][j][k] + Fl[i][j-1][k]);
                double d2l = (Fl[i][j][k+1] - 2*Fl[i][j][k] + Fl[i][j][k-1]);
                double d12l = (Fl[i][j+1][k+1] - Fl[i][j+1][k-1] - Fl[i][j-1][k+1] + Fl[i][j-1][k-1])/(4);
                if(d12l>0){
                    corr = 0.3;
                }else{
                    corr = 0.5;
                }
                tie(s1,s2) = findMinValue(d1l,d2l,d12l,smax1,smax2,smin1,smin2,corr);
                
                Fl[i+1][j][k] = (1-df)*Fl[i][j][k] + (0.5*s1*s1*j*j*dt)*d1l
                + (0.5*s2*s2*k*k*dt)*d2l
                + (0.5*r*k*dt)*(Fl[i][j][k+1] - Fl[i][j][k-1])
                + (0.5*r*j*dt)*(Fl[i][j+1][k] - Fl[i][j-1][k])
                + (0.25*corr*s1*s2*j*k*dt)*4*d12l;
                
                double det = d1l*d2l - d12l*d12l;
                if(det==0 && d1l!=0 && d2l!=0 && d12l!=0){
                    cnt+=1;
                    // so apparently this is very rare and only happens when second derivatives are super small or zero
                }
            }
        }
    }
    cout << cnt << "\n";
    
    cout<<"40 40"<<" "<<Fl[NT-1][40][40]<<"\t"<<Fu[NT-1][40][40]<<"\n";
    
    for(int j = 1; j<NS1-1; j++){
        for(int k=1; k<NS2-1; k++){
            bsbUpper << Fu[NT-1][j][k] << "\n";
            bsbLower << Fl[NT-1][j][k] << "\n";
        }
        prices << j*dS << "\n";
    }
    prices.close();
    bsbUpper.close();
    bsbLower.close();
}

void Bspread_BSB_bounds(double r, double smax1, double smax2, double smin1, double smin2, double corr, double T, double E1, double E2, double NT, double NS){
    ofstream bsbUpper("./data/2Dbsb_spread_U.txt");
    ofstream bsbLower("./data/2Dbsb_spread_L.txt");
    ofstream prices("./data/2Dbsb_spread_prices.txt");

    double s1 = smax1, s2 = smax2;
    double Eavg = 0.5*(E1+E2);
    int NS1 = NS, NS2 = NS;
    
    double dt = T/NT;
    double dS = 1;
    double df = r*dt;
    
    // INITIALISE PRICE ARRAY
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
    
    // *******************  Set Boundary and Initial Conditions for the payoff **********************************
    //Initial
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
        //Boundary
        
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
    
    cout<<"40 40"<<" "<<Fl[40][40]<<"\t"<<Fu[40][40]<<"\n";
    
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
void Bspread_BS_bounds(double r, double smax1, double smax2, double smin1, double smin2, double corr, double T, double E1, double E2, double NT, double NS){
    ofstream bsUpper("./data/2Dbs_spread_U.txt");
    ofstream bsLower("./data/2Dbs_spread_L.txt");
    ofstream prices("./data/2Dbs_spread_prices.txt");
    
    double s1 = smax1, s2 = smax2;
    double Eavg = 0.5*(E1+E2);
    int NS1 = NS, NS2 = NS;
    
    double dt = T/NT;
    double dS = 1;
    double df = r*dt;

    // INITIALISE PRICE ARRAY
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
    
    // *******************  Set Boundary and Initial Conditions for the payoff **********************************
    //Initial
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
        //Boundary
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
    cout << "How long?";
    prices.close();
    bsUpper.close();
    bsLower.close();
    
}

void Cspread_BSB_bounds(double r, double smax1, double smax2, double smin1, double smin2, double corr, double T, double E1, double E2, double NT, double NS){
    ofstream bsbUpper("./data/2Dbsb_spread_U.txt");
    ofstream bsbLower("./data/2Dbsb_spread_L.txt");
    ofstream prices("./data/2Dbsb_spread_prices.txt");
    
    double s1 = smax1, s2 = smax2;
    int NS1 = NS, NS2 = NS;
    
    double dt = T/NT;
    double dS = 1;
    double df = r*dt;
    
    // INITIALISE PRICE ARRAY
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
    
    // *******************  Set Boundary and Initial Conditions for the payoff **********************************
    //Initial
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
        //Boundary
        
        s1 = smax1;
        s2 = smax2;
        Fu_n[0][0] = 0;
        Fu_n[NS1-1][NS2-1] = (E2-E1)*exp(-r*i*dt);
        for(int k=1; k<NS2-1;k++){
            Fu_n[0][k] = Fu[0][k] + (0.5*s2*s2*k*k*dt)*(Fu[0][k+1] - 2*Fu[0][k] + Fu[0][k-1])
            -df*Fu[0][k] + (0.5*r*k*dt)*(Fu[0][k+1] - Fu[0][k-1]);
        }
        for(int j=1; j<NS1-1; j++){
            Fu_n[j][0] = Fu[j][0] + (0.5*s1*s1*j*j*dt)*(Fu[j+1][0] - 2*Fu[j][0] + Fu[j-1][0])
            -df*Fu[j][0] + (0.5*r*j*dt)*(Fu[j+1][0] - Fu[j-1][0]);
        }
        for(int k=0; k<NS2; k++){
            Fu_n[NS-1][k] = (E2-E1)*exp(-r*i*dt);
        }
        for(int j=0; j<NS1; j++){
            Fu_n[j][NS-1] = (E2-E1)*exp(-r*i*dt);
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
void Cspread_BS_bounds(double r, double smax1, double smax2, double smin1, double smin2, double corr, double T, double E1, double E2, double NT, double NS){
    ofstream bsUpper("./data/2Dbs_spread_U.txt");
    ofstream bsLower("./data/2Dbs_spread_L.txt");
    ofstream prices("./data/2Dbs_spread_prices.txt");
    
    double s1 = smax1, s2 = smax2;
    int NS1 = NS, NS2 = NS;
    
    double dt = T/NT;
    double dS = 1;
    double df = r*dt;
    
    // INITIALISE PRICE ARRAY
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
    
    // *******************  Set Boundary and Initial Conditions for the payoff **********************************
    //Initial
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
        //Boundary
        Fu_n[0][0][0] = 0;Fu_n[1][0][0] = 0;Fu_n[2][0][0] = 0;
        Fu_n[0][NS1-1][NS2-1] = (NS-1)*dS - E1*exp(-r*i*dt);
        Fu_n[1][NS1-1][NS2-1] = -(NS-1)*dS + E2*exp(-r*i*dt);
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
            Fu_n[0][NS-1][k] = (NS-1)*dS - E1*exp(-r*i*dt);
            Fu_n[1][NS-1][k] = -(NS-1)*dS + E2*exp(-r*i*dt);
        }
        for(int j=0; j<NS1; j++){
            Fu_n[0][j][NS-1] = (NS-1)*dS - E1*exp(-r*i*dt);
            Fu_n[1][j][NS-1] = -(NS-1)*dS + E2*exp(-r*i*dt);
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
    int NT = 2500, NS = 100;
    /*
    double r = 0.05, corr = 1;
    double smax1 = 0.4, smax2 = 0.4, smin1 = 0.1, smin2 = 0.1;
    double T = 0.5;
    double E1 = 90, E2 = 100;
    int NT = 2500, NS = 100;
    */
    
    int i = 1; //1 for calendar, else for butterfly
    if(i == 1){
        Cspread_BSB_bounds(r, smax1, smax2, smin1, smin2, corr, T, E1, E2, NT, NS);
        //Cspread_BS_bounds(r, smax1, smax2, smin1, smin2, corr, T, E1, E2, NT, NS);
    }else{
        Bspread_BSB_bounds(r, smax1, smax2, smin1, smin2, corr, T, E1, E2, NT, NS);
        Bspread_BS_bounds(r, smax1, smax2, smin1, smin2, corr, T, E1, E2, NT, NS);
    }
    
    //system("sh scripts/uncVol2D_script.sh");
}