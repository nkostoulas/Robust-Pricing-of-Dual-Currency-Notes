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

using namespace std;

void twoStateBSBVanilla();
void twoStateBSBSpread();
void twoStateBSBSpread_2();

void Bspread_BSB_bounds(double r, double smax1, double smax2, double smin1, double smin2, double corr, double T, double E1, double E2, double NT, double NS);
void Bspread_BS_bounds(double r, double smax1, double smax2, double smin1, double smin2, double corr, double T, double E1, double E2, double NT, double NS);
void Cspread_BSB_bounds(double r, double smax1, double smax2, double smin1, double smin2, double corr, double T, double E1, double E2, double NT, double NS);
void Cspread_BS_bounds(double r, double smax1, double smax2, double smin1, double smin2, double corr, double T, double E1, double E2, double NT, double NS);

void twoStateModelling();

#endif /* twoStateModelling_hpp */
