//
//  main.cpp
//  uncertainVolatilityModel
//
//  Created by Nick Kostoulas on 01/04/16.
//  Copyright © 2016 Nick Kostoulas. All rights reserved.
//

#include "callSpreadExample.hpp"
#include "calendarSpreadExample.hpp"
#include "stochasticVolatilityExample.hpp"
#include "dualCurrencyNoteExample.hpp"
#include <iostream>
#include <ctime>

using namespace std;

int main(int argc, const char * argv[]) {
    clock_t begin = clock();
    
    // 1. Bull Spread and Calendar Spread on 2 options
    callSpreadExample();
    calendarSpreadExample();
    system("sh scripts/uncVol_script.sh");
    
    // 2. Stoch vol hedging BSB vs BS
    stochasticVolatilityExample();
    
    // 3. Pricing of dual currency notes
    dualCurrencyNoteExample();
    
    clock_t end = clock();
    double elapsed_secs = (double)(end - begin);
    cout <<  elapsed_secs/ CLOCKS_PER_SEC << "\t";
}





