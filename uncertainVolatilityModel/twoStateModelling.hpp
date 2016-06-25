//
//  twoStateModelling.hpp
//  uncertainVolatilityModel
//
//  Created by Nick Kostoulas on 10/06/16.
//  Copyright © 2016 Nick Kostoulas. All rights reserved.
//

#ifndef twoStateModelling_hpp
#define twoStateModelling_hpp

#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "toms462.hpp"
#include "BS.hpp"

using namespace std;

double normalCDF(double x);
double minOptionValue(double S1, double S2, double E, double s1, double s2, double corr, double r, double T);
double maxOptionValue(double S1, double S2, double E, double s1, double s2, double corr, double r, double T);
void Bspread_BSB_bounds(double r, double smax1, double smax2, double smin1, double smin2, double corr, double T, double E1, double E2, double NT, double NS, double dS);
void Bspread_BS_bounds(double r, double smax1, double smax2, double smin1, double smin2, double corr, double T, double E1, double E2, double NT, double NS, double dS);
void Bspread_BScf_bounds(double r, double smax1, double smax2, double smin1, double smin2, double corr, double T, double E1, double E2, double NT, double NS, double dS);
void Cspread_BSB_bounds(double r, double smax1, double smax2, double smin1, double smin2, double corr, double T, double E1, double E2, double NT, double NS, double dS);
void Cspread_BS_bounds(double r, double smax1, double smax2, double smin1, double smin2, double corr, double T, double E1, double E2, double NT, double NS, double dS);
void Cspread_BScf_bounds(double r, double smax1, double smax2, double smin1, double smin2, double corr, double T, double E1, double E2, double NT, double NS, double dS);

void vanilla_BSB_DC(double rf, double rd, double smax, double smin, double corr, double T, double E1, double E2, int NT, int NS, int nOfOpts, int DCorPRDC, double dS);
void vanilla_BS_DC(double rf, double rd, double smax, double smin, double corr, double T, double E1, double E2, int NT, int NS, int nOfOpts, int DCorPRDC, double dS);

pair<double, double> findMinValue(double d1, double d2, double d12, double smax1, double smax2, double smin1, double smin2, double corr);
pair<double, double> findMaxValue(double d1, double d2, double d12, double smax1, double smax2, double smin1, double smin2, double corr);
void twoStateModelling();

#endif /* twoStateModelling_hpp */
