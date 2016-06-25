//
//  3AssetModelling.cpp
//  uncertainVolatilityModel
//
//  Created by Nick Kostoulas on 23/06/16.
//  Copyright © 2016 Nick Kostoulas. All rights reserved.
//

#include "3AssetModelling.hpp"
#include "twoStateModelling.hpp"

double optValue(double d1, double d2, double d3, double d12, double d13, double d23, double s1, double s2, double s3){
    return d1*s1*s1 + d2*s2*s2 + d3*s3*s3 + 2*d12*s1*s2 + 2*d13*s1*s3 + 2*d23*s2*s3;
}
tuple<double, double, double> findMax(double d1, double d2, double d3, double d12, double d13, double d23, double smax, double smin, double corr){
    d12 = d12*corr; d13 = d13*corr; d23 = d23*corr;
    
    double det = d1*(d2*d3 - d23*d23) - d12*(d12*d3 - d13*d23) + d13*(d12*d23 - d13*d2);
    double max = -INT_MAX, val = 0.0;
    
    double s1_return = smin, s2_return = smin, s3_return = smin;
    double s1 = smin, s2 = smin, s3 = smin;
    
    val = optValue(d1,d2,d3,d12,d13,d23,s1=smax,s2=smax,s3=smax);
    if(val > max){
        if(det!=0){
            max = val;
            s1_return = s1;
            s2_return = s2;
            s3_return = s3;
        }
    }
    val = optValue(d1,d2,d3,d12,d13,d23,s1=smax,s2=smax,s3=smin);
    if(val > max){
        if(det!=0){
            max = val;
            s1_return = s1;
            s2_return = s2;
            s3_return = s3;
        }
    }
    val = optValue(d1,d2,d3,d12,d13,d23,s1=smax,s2=smin,s3=smax);
    if(val > max){
        if(det!=0){
            max = val;
            s1_return = s1;
            s2_return = s2;
            s3_return = s3;
        }
    }
    val = optValue(d1,d2,d3,d12,d13,d23,s1=smax,s2=smin,s3=smin);
    if(val > max){
        if(det!=0){
            max = val;
            s1_return = s1;
            s2_return = s2;
            s3_return = s3;
        }
    }
    val = optValue(d1,d2,d3,d12,d13,d23,s1=smin,s2=smax,s3=smax);
    if(val > max){
        if(det!=0){
            max = val;
            s1_return = s1;
            s2_return = s2;
            s3_return = s3;
        }
    }
    val = optValue(d1,d2,d3,d12,d13,d23,s1=smin,s2=smax,s3=smin);
    if(val > max){
        if(det!=0){
            max = val;
            s1_return = s1;
            s2_return = s2;
            s3_return = s3;
        }
    }
    val = optValue(d1,d2,d3,d12,d13,d23,s1=smin,s2=smin,s3=smax);
    if(val > max){
        if(det!=0){
            max = val;
            s1_return = s1;
            s2_return = s2;
            s3_return = s3;
        }
    }
    val = optValue(d1,d2,d3,d12,d13,d23,s1=smin,s2=smin,s3=smin);
    if(val > max){
        if(det!=0){
            max = val;
            s1_return = s1;
            s2_return = s2;
            s3_return = s3;
        }
    }
    return std::make_tuple(s1_return, s2_return, s3_return);
}
tuple<double, double, double> findMin(double d1, double d2, double d3, double d12, double d13, double d23, double smax, double smin, double corr){
    d12 = d12*corr; d13 = d13*corr; d23 = d23*corr;
    
    double det = d1*(d2*d3 - d23*d23) - d12*(d12*d3 - d13*d23) + d13*(d12*d23 - d13*d2);
    
    double min = INT_MAX, val = 0.0;
    
    double s1_return = smin, s2_return = smin, s3_return = smin;
    double s1 = smin, s2 = smin, s3 = smin;
    
    val = optValue(d1,d2,d3,d12,d13,d23,s1=smax,s2=smax,s3=smax);
    if(val < min){
        if(det!=0){
            min = val;
            s1_return = s1;
            s2_return = s2;
            s3_return = s3;
        }
    }
    val = optValue(d1,d2,d3,d12,d13,d23,s1=smax,s2=smax,s3=smin);
    if(val < min){
        if(det!=0){
            min = val;
            s1_return = s1;
            s2_return = s2;
            s3_return = s3;
        }
    }
    val = optValue(d1,d2,d3,d12,d13,d23,s1=smax,s2=smin,s3=smax);
    if(val < min){
        if(det!=0){
            min = val;
            s1_return = s1;
            s2_return = s2;
            s3_return = s3;
        }
    }
    val = optValue(d1,d2,d3,d12,d13,d23,s1=smax,s2=smin,s3=smin);
    if(val < min){
        if(det!=0){
            min = val;
            s1_return = s1;
            s2_return = s2;
            s3_return = s3;
        }
    }
    val = optValue(d1,d2,d3,d12,d13,d23,s1=smin,s2=smax,s3=smax);
    if(val < min){
        if(det!=0){
            min = val;
            s1_return = s1;
            s2_return = s2;
            s3_return = s3;
        }
    }
    val = optValue(d1,d2,d3,d12,d13,d23,s1=smin,s2=smax,s3=smin);
    if(val < min){
        if(det!=0){
            min = val;
            s1_return = s1;
            s2_return = s2;
            s3_return = s3;
        }
    }
    val = optValue(d1,d2,d3,d12,d13,d23,s1=smin,s2=smin,s3=smax);
    if(val < min){
        if(det!=0){
            min = val;
            s1_return = s1;
            s2_return = s2;
            s3_return = s3;
        }
    }
    val = optValue(d1,d2,d3,d12,d13,d23,s1=smin,s2=smin,s3=smin);
    if(val < min){
        if(det!=0){
            min = val;
            s1_return = s1;
            s2_return = s2;
            s3_return = s3;
        }
    }
    return std::make_tuple(s1_return, s2_return, s3_return);
}
 
void threeAssetModelling(){
    //*********** Variable Definitions *************************//
    double corr = 0.5, smax = 0.5, smin = 0.1;
    double T = 0.5;
    int NT = 1000, NS = 50; double dS = 2;
    
    int nOfOpts = 1;
    double rf = 0.15, rd = 0.10;
    double E1 = 45, E2 = 55;
    
    NT = NT*nOfOpts;
    T = T*nOfOpts;
    int DCorPRDC = 2; //1fordc 2forprdc
    
    // Dual Currency Note Definitions
    double r = rf-rd;
    double s1, s2, s3;
    
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
    vector<vector<vector<double>>> Fu, Fu_n, Fl, Fl_n;
    Fu.resize(NS);
    Fu_n.resize(NS);
    Fl.resize(NS);
    Fl_n.resize(NS);
    for (int j = 0; j < NS; ++j){
        Fu[j].resize(NS);
        Fl[j].resize(NS);
        Fu_n[j].resize(NS);
        Fl_n[j].resize(NS);
        for (int k = 0; k < NS; k++){
            Fu[j][k].resize(NS);
            Fl[j][k].resize(NS);
            Fu_n[j][k].resize(NS);
            Fl_n[j][k].resize(NS);
        }
    }
    
    // *******************  Set Initial Conditions for the payoff **********************************
    for(int i=0; i<NS; i++){
        for(int j=0; j<NS1; j++){
            for(int k=0; k<NS2; k++){
                double p1 = i*dS, p2 = j*dS, p3 = k*dS;
                Fu[i][j][k] = L0*max(max(max(p1,p2),p3)-FXstrike0,0.0) - L1*max(max(max(p1,p2),p3)-FXstrike1,0.0);
            }
        }
    }
    Fl = Fu;
    clock_t end1 = clock();
    double initTime = (double)(end1 - begin1);
    // *******************  Calculate derivative price using finite difference ************************************
    for(int t=0; t<NT-1; t++){
        clock_t begin2 = clock();
        //Boundary Conditions
        if((t+1)%NTstep==0){
            for(int i=0; i<NS; i++){
                for(int j=0; j<NS; j++){
                    for(int k=0; k<NS; k++){
                        double p1 = i*dS, p2 = j*dS, p3 = k*dS;
                        Fu[i][j][k] = L0*max(max(max(p1,p2),p3)-FXstrike0,0.0) - L1*max(max(max(p1,p2),p3)-FXstrike1,0.0);
                    }
                }
            }
        }else{
            Fu_n[0][0][0] = 0;
            for(int i=1; i<NS-1; i++){
                Fu_n[0][0][i] = Fu[0][0][i] + (0.5*smax*smax*i*i*dt)*(Fu[0][0][i+1] - 2*Fu[0][0][i] + Fu[0][0][i-1])
                -df*Fu[0][0][i] + (0.5*r*i*dt)*(Fu[0][0][i+1] - Fu[0][0][i-1]);
                
                Fu_n[0][i][0] = Fu[0][i][0] + (0.5*smax*smax*i*i*dt)*(Fu[0][i+1][0] - 2*Fu[0][i][0] + Fu[0][i-1][0])
                -df*Fu[0][i][0] + (0.5*r*i*dt)*(Fu[0][i+1][0] - Fu[0][i-1][0]);
                
                Fu_n[i][0][0] = Fu[i][0][0] + (0.5*smax*smax*i*i*dt)*(Fu[i+1][0][0] - 2*Fu[i][0][0] + Fu[i-1][0][0])
                -df*Fu[i][0][0] + (0.5*r*i*dt)*(Fu[i+1][0][0] - Fu[i-1][0][0]);
                
                for(int j=1; j<NS-1; j++){
                    Fu_n[i][j][0] = Fu[i][j][0] -df*Fu[i][j][0];
                    Fu_n[i][0][j] = Fu[i][0][j] -df*Fu[i][0][j];
                    Fu_n[0][i][j] = Fu[0][i][j] -df*Fu[0][i][j];
                }
            }
            for(int i=0; i<NS; i++){
                double val = (L0-L1)*FXmax + (FXstrike1-FXstrike0)*exp(-r*(i+1)*dt);
                Fu_n[NS-1][i][0] = val;
                Fu_n[0][NS-1][i] = val;
                Fu_n[i][NS-1][0] = val;
                Fu_n[0][i][NS-1] = val;
                Fu_n[i][0][NS-1] = val;
                Fu_n[NS-1][0][i] = val;
                Fu_n[i][NS-1][NS-1] = val;
                Fu_n[NS-1][NS-1][i] = val;
                Fu_n[NS-1][i][NS-1] = val;
            }
        }
        Fl_n = Fu_n;
        clock_t end2 = clock();
        boundTime += (double)(end2 - begin2);
        
        clock_t begin3 = clock();
        //Finite difference
        for(int i=1; i<NS-1; i++){
            for(int j=1; j<NS-1; j++){
                for(int k=1; k<NS-1; k++){
                    
                    double d1u = (Fu[i+1][j][k] - 2*Fu[i][j][k] + Fu[i-1][j][k]);
                    double d2u = (Fu[i][j+1][k] - 2*Fu[i][j][k] + Fu[i][j-1][k]);
                    double d3u = (Fu[i][j][k+1] - 2*Fu[i][j][k] + Fu[i][j][k-1]);
                    double d12u = (Fu[i+1][j+1][k] - Fu[i+1][j-1][k] - Fu[i-1][j+1][k] + Fu[i-1][j-1][k])/(4);
                    double d13u = (Fu[i+1][j][k+1] - Fu[i+1][j][k-1] - Fu[i-1][j][k+1] + Fu[i-1][j][k-1])/(4);
                    double d23u = (Fu[i][j+1][k+1] - Fu[i][j-1][k+1] - Fu[i][j+1][k-1] + Fu[i][j-1][k-1])/(4);
                
                    //s1 = smax; s2 = smin; s3 = smin;
                    tie(s1,s2,s3) = findMax(d1u,d2u,d3u,d12u,d13u,d23u,smax,smin,corr);
                    
                    Fu_n[i][j][k] = (1-df)*Fu[i][j][k] + (0.5*s1*s1*i*i*dt)*d1u
                                                       + (0.5*s2*s2*j*j*dt)*d2u
                                                       + (0.5*s3*s3*k*k*dt)*d3u
                                    + (0.5*r*i*dt)*(Fu[i+1][j][k] - Fu[i-1][j][k])
                                    + (0.5*r*j*dt)*(Fu[i][j+1][k] - Fu[i][j-1][k])
                                    + (0.5*r*k*dt)*(Fu[i][j][k+1] - Fu[i][j][k-1])
                                    + (0.25*corr*s1*s2*i*j*dt)*4*d12u
                                    + (0.25*corr*s1*s3*i*k*dt)*4*d13u
                                    + (0.25*corr*s2*s3*j*k*dt)*4*d23u;
                    
                    
                    double d1l = (Fl[i+1][j][k] - 2*Fl[i][j][k] + Fl[i-1][j][k]);
                    double d2l = (Fl[i][j+1][k] - 2*Fl[i][j][k] + Fl[i][j-1][k]);
                    double d3l = (Fl[i][j][k+1] - 2*Fl[i][j][k] + Fl[i][j][k-1]);
                    double d12l = (Fl[i+1][j+1][k] - Fl[i+1][j-1][k] - Fl[i-1][j+1][k] + Fl[i-1][j-1][k])/(4);
                    double d13l = (Fl[i+1][j][k+1] - Fl[i+1][j][k-1] - Fl[i-1][j][k+1] + Fl[i-1][j][k-1])/(4);
                    double d23l = (Fl[i][j+1][k+1] - Fl[i][j-1][k+1] - Fl[i][j+1][k-1] + Fl[i][j-1][k-1])/(4);
                    
                    //s1 = smin; s2 = smax; s3 = smax;
                    tie(s1,s2,s3) = findMin(d1l,d2l,d3l,d12l,d13l,d23l,smax,smin,corr);
                    
                    Fl_n[i][j][k] = (1-df)*Fl[i][j][k] + (0.5*s1*s1*i*i*dt)*d1l
                                                       + (0.5*s2*s2*j*j*dt)*d2l
                                                       + (0.5*s3*s3*k*k*dt)*d3l
                                    + (0.5*r*i*dt)*(Fl[i+1][j][k] - Fl[i-1][j][k])
                                    + (0.5*r*j*dt)*(Fl[i][j+1][k] - Fl[i][j-1][k])
                                    + (0.5*r*k*dt)*(Fl[i][j][k+1] - Fl[i][j][k-1])
                                    + (0.25*corr*s1*s2*i*j*dt)*4*d12l
                                    + (0.25*corr*s1*s3*i*k*dt)*4*d13l
                                    + (0.25*corr*s2*s3*j*k*dt)*4*d23l;
                }
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
    cout<<Fu[21][21][21]<<" "<<Fl[21][21][21]<<"\n";

}