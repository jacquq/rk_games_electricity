#ifndef UTILS_SDE
#define UTILS_SDE

#include <string>
#include <tuple>
#include <vector>
#include <map>
#include "math.h"
#include <random>
#include <utility>
#include <algorithm>
#include <omp.h>

#include "methods.hpp"
#include "data.hpp"

using namespace std;

/* ----------------------------

This class provides a simulator of trajectories (Monte-Carlo) for a given cumulative distribution,
and so a given effort. A first function computes the trajectories using a Euler-Maruyama discretization
scheme, and a second fuction computes the new cumulative distribution once the trajectories are drawn.

---------------------------- */

// ---------------------------------------------------------------------------------
void MonteCarlo::run(int k, const vector<double>& R_mu, int num_threads, int seed, bool store){
	/*
	This function computes the trajectories for a process Xt defined by

		dXt = a*(t,Xt)dt + sigma dBt ,

	where a* is the optimal effort associated with the given cumulative distribution and
	Bt a standard Brownian motion, see [2].
	It uses a Euler-Maruyama discretization scheme, which convergence is guarenteed when the 
	number of simulations goes to infinity, see [*].
	*/

	// clear vectors
	_store_trajectory = vector<vector<double>>(_nb_simu);
	_store_effortcost = vector<vector<double>>(_nb_simu);
	_store_brownian   = vector<vector<double>>(_nb_simu);
	_terminal_values  = vector<double>(_nb_simu);

	// local variables
	int nb_points = (int) floor(_data->_T / _delta_t);
	double h = _data->_T / nb_points;				// caution: do not touch t = T, effort = +oo, outch !
	double sqrt_h = sqrt(h);

	// standard normal distribution
	random_device rd;
	mt19937 gen(seed); // tmp, otherwise rd()
	normal_distribution<double> d(0.,1.);

	// create trajectory according to Euler-Maruyama scheme
	int iter = 0;
	//_best_utility = 0.;
	for (int i = 0; i < _nb_simu; ++i) {
		double effortcost = 0.;
		double brownian = 0.;
		if (store) {
			_store_trajectory[i] = vector<double>(nb_points);
			_store_trajectory[i][0] = _data->_X0[k];
			_store_effortcost[i] = vector<double>(nb_points);
			_store_effortcost[i][0] = 0.;
			_store_brownian[i]   = vector<double>(nb_points);
			_store_brownian[i][0]   = 0.;
		}
		double x = _data->_X0[k];
		double t, a, w;
		for (int n = 1; n < nb_points; ++n){
			t = h*n;
			a = _data->get_optimal_effort(k,t,x,R_mu);
			w = sqrt_h * d(gen);
			brownian += w;
			x += a * h + _data->_sigma[k] * w; 
			effortcost += pow(a,2) * (_data->_c[k] + _data->_c_der[k]*t)  * h;
			if (store) {
				_store_trajectory[i][n] = x;
				_store_effortcost[i][n] = effortcost;
				_store_brownian[i][n]   = brownian;
			}
		}
		int j_x = 0;
		//while (_data->_range_x[j_x + 1] < x) { ++ j_x; }
        //_best_utility += R_mu[j_x] + (R_mu[j_x+1] - R_mu[j_x]) / (_data->_range_x[j_x+1] - _data->_range_x[j_x]) 
        //					*(x - _data->_range_x[j_x]);
		_terminal_values[i] = x;

		++iter;
	}
	//_best_utility /= _nb_simu;

};

// ---------------------------------------------------------------------------------
void MonteCarlo::fill_new_cumprob(vector<double>& new_cumprob){
	/*
	This function computes the new cumulative distribution once the trajectories are drawn.
	*/

	// reinit new_cumprob 
	new_cumprob = vector<double>(_data->_nb_x, 0.);

	// local variables
	double h = (_data->_x_max - _data->_x_min) / (_data->_nb_x - 1);
	int idx_cumprob;
	double loss = 0.;

	// sort terminal values (ascending sort)
	sort(_terminal_values.begin(), _terminal_values.end());

	// fill cumprob
	for (vector<double>::iterator it = _terminal_values.begin(); it != _terminal_values.end(); ++it) {

		if (*it > _data->_x_min && *it < _data->_x_max){
			idx_cumprob = (int) floor( (*it - _data->_x_min) / h); // rouding of index

			for (int i = idx_cumprob + 1; i < _data->_nb_x; ++i) {
				new_cumprob[i] += 1./_terminal_values.size();
			}
		}
		else {
			loss += 1./_terminal_values.size();
		}
	}
	if (loss > 1e-5) cout << "[WARNING] loss = " << loss << endl; // as we are not in \R but on
															      //   restricted interval, we need to 
																  //   take care of the loss

};

#endif