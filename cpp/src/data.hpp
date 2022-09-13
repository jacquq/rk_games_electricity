#ifndef CLASS_PARAMS_HPP
#define CLASS_PARAMS_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "data.hpp"
#include "utils.hpp"

using namespace std;

struct OptimData {

	int _nb_k;

	vector<double> _sigma;
	double _T;
	double _p;
	vector<double> _c;
	vector<double> _c_der;
	vector<double> _X0;

	double _alpha;
	vector<double> _cr;

	double _x_min;
	double _x_max;
	int _nb_x;
	vector<double> _range_x;


	OptimData(){};
	OptimData(Config* config, char delim = '\t');
	~OptimData(){};

	double get_Reward_pure_ranked(double rank);
	double get_Reward(int k, const vector<double>& cumprob, double x);
	void get_u_and_ux(int k, double t, double xt, const vector<double>& cumprob, double& u, double& u_x);
	double get_optimal_effort(int k, double t, double xt, const vector<double>& cumprob);
	vector<double> get_new_iterate_cdf(int k, const vector<double>& R_mu);

	double get_Wass_dist(const vector<double>& cumprob_1, const vector<double>& cumprob_2);


	// ----- retailer
	double retailer_Beta(int k, const vector<double>& b_rangex);
	double retailer_intReward_xpart(const vector<double>& x_cma, const vector<double>& B_x, const vector<double>& cumprob);	
    double retailer_intReward_rankpart(const vector<double>& r_cma, const vector<double>& reward);   
	double retailer_gain(int k, double mean);
	pair<double,double> retailer_constraintV0(int k, double int_Beta);

};

#endif
