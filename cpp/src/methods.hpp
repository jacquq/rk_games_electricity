#ifndef METHODS_HPP
#define METHODS_HPP

#include <string>
#include <tuple>
#include <vector>
#include <map>

#include "gnuplot_i.hpp"

#include "data.hpp"
#include "utils.hpp"

// ======================= equilibrium.cpp ==============================

void plot_its_equilibrium(vector<string> legend, vector<vector<double>> list_its);

// ======================= retailer.cpp ==============================
struct OptimalRewardFinder {
	OptimData* _data;
	Config* _config;
	vector<double> _r_cma;
	vector<double> _x_cma;

	vector<double> _opt_cma;
	vector<double> _opt_B_x;
	vector<double> _opt_B_r;
	vector<vector<double>> _cumprob;
	vector<vector<double>> _R_mu;
	
	bool _fixed_price;

	double _bigM;
	double _bigM_x;

	OptimalRewardFinder(){};
	~OptimalRewardFinder(){};
	OptimalRewardFinder(OptimData* data, Config* config);

	void run(int verbose = 1);
	double fitfun(const double *v, const int N, int verbose = 1, bool store_cumprob = false);
	void plot_reward(bool plot_opt_reward = false, bool usetex = false, bool show = true);

	void construct_reward_from_CMA(const double* v, vector<double>& B_r, vector<double>& B_x);
	void construct_reward_from_cumprob(const vector<double>& cumprob,  
                const vector<double>& B_r, const vector<double>& B_x, 
                vector<double>& R_mu);

	vector<double> compute_equilibrium(int k, const vector<double>& B_r, 
		const vector<double>& B_x, double& mean, vector<double>& cumprob, int verbose = 1);

	void get_equilibrium(int k, const vector<double>& reward,
		double& mean, vector<double>& cumprob);

};
tuple<double,vector<double>,vector<double>> get_optimal_reward(OptimData* data, int k);     // only for homogeneous pop


// ======================= monte_carlo.cpp ==============================
struct MonteCarlo {
	OptimData* _data;
	double _delta_t;
	int _nb_simu;

	vector<vector<double>> _store_trajectory;
	vector<vector<double>> _store_effortcost;
	vector<vector<double>> _store_brownian;
	vector<double> _terminal_values;

	MonteCarlo(){};
	MonteCarlo(	OptimData* data, double delta_t, int nb_simu) {
		_data   = data;
		_delta_t = delta_t;
		_nb_simu = nb_simu;
	}

	void run(int k, const vector<double>& cumprob, int num_threads, int seed = 0, bool store = false);
	void fill_new_cumprob(vector<double>& new_cumprob);
	void plot_trajectories(int nb_traj, bool show = true);
	void plot_effortcost(int nb_traj, bool show = true);
	void plot_brownian(int nb_traj, bool show = true);
	void plot_mean(map<string, string> keywords, bool show = true);
};

void plot_cumprob(const vector<double>& range_x, const vector<double>& cumprob, const map<string,string>& keywords, bool show = true);

#endif