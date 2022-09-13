#ifndef UTILS_INT
#define UTILS_INT

#include <string>
#include <tuple>
#include <vector>
#include <map>
#include <iostream>

#include "utils.hpp"

using namespace std;

// --------------------------------------------------------------------------
double integ(const vector<double>& f, double a, double b, string method){
	double n = f.size();
	double h = (b - a)/(n-1);

	double val = 0.;
	if (method == "S"){
		for (int i = 1; i < f.size()/2; ++i){
			val += 2*f[2*i];
		}
		for (int i = 1; i <= f.size()/2; ++i){
			val += 4*f[2*i-1];
		}
		val += f[0] + f[f.size()-1];
		val *= h/3;
	}
	else if (method == "R"){
		for (int i = 0; i < f.size()-1; ++i){
			val += f[i];
		}
		val *= h;
	}
	return val;
};

// --------------------------------------------------------------
double integ(const vector<double>& x, const vector<double>& f){
	double val = 0.;
	for (int i = 0; i < f.size()-1; ++i){
		val += 0.5 * (f[i] + f[i+1]) * (x[i+1] - x[i]);
	}
	return val;
};


#endif