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
#include <iostream>
#include <ctime>

using namespace std;

int main(int argc, const char * argv[]) {
    clock_t begin = clock();
    //callSpreadExample();
    //calendarSpreadExample();
    stochasticVolatilityExample();
    clock_t end = clock();
    double elapsed_secs = (double)(end - begin);
    cout <<  elapsed_secs/ CLOCKS_PER_SEC << "\t";
}





