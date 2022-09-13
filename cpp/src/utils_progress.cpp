#ifndef CLASS_PROGRESS_HPP
#define CLASS_PROGRESS_HPP


#define PBSTR "============================================================"
#define PBWIDTH 60

#include <iostream>
#include <fstream>
#include <string>

using namespace std;


void printProgress(double percentage) {
    percentage = max(percentage,0.);
    int val = (int) (percentage);
    int lpad = (int) (percentage * PBWIDTH / 100);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
};

#endif