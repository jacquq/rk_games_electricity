#ifndef RETAILER
#define RETAILER

#include <string>
#include <tuple>
#include <vector>
#include <map>
#include <math.h>
#include <chrono>

#include "cmaes.h"

#include "data.hpp"
#include "utils.hpp"
#include "methods.hpp"
#include <omp.h>

using namespace std;
using namespace chrono;

/* ----------------------------

This class provides functions for the search of the optimal reward.
In the general case, no analytic solution is known. We therefore use the solver CMA-ES through
the library:
    https://github.com/AlexanderFabisch/CMA-ESpp
In the simpler case where the population is homogeneous, analytic solution is known.


---------------------------- */


bool is_feasible(double const*pop_i, double *lb, double *ub, int dim) {
    /*
    Intern function for bounds in CMA-ES
    */

    for (int i = 0; i < dim; ++i){
        if ((pop_i[i] < lb[i]) || (pop_i[i] > ub[i])) return false; 
    }
    return true;
}

// --------------------------------------
double OptimalRewardFinder::fitfun(const double *v, const int N, int verbose, bool store_cumprob)  {

    double obj = 0.;

    // construct reward B(r) from variable of CMA
    vector<double> B_r, B_x;
    construct_reward_from_CMA(v, B_r, B_x);

    for (int k = 0; k < _data->_nb_k; ++k){
        
        double mean;
        vector<double> cumprob_k, R_mu;

        // use the analytic solution if possible
        if (_fixed_price) {
            // Analytic result in this case
            get_equilibrium(k, B_r, mean, cumprob_k); // return mean, quantile
        }
        // or compute the fixed point
        else {
            compute_equilibrium(k, B_r, B_x, mean, cumprob_k,0);
        }

        construct_reward_from_cumprob(cumprob_k, B_r, B_x, R_mu);
        // compute integrale of reward + integrale of beta (Value for consumers)
        double int_add_reward = _data->retailer_intReward_rankpart(_r_cma, B_r);
        double beta           = _data->retailer_Beta(k,R_mu);
        if (!_fixed_price)  int_add_reward += _data->retailer_intReward_xpart(_x_cma, B_x,cumprob_k);
        else                int_add_reward -= _data->_p * mean;
        double gain           = _data->retailer_gain(k,mean);
        pair<double,double> V = _data->retailer_constraintV0(k,beta);

        obj += int_add_reward - gain + 1e3*max(V.second - V.first,0.); // penalization of the constraint

        if (verbose >= 2) { 
            cout << "cluster " << k << ", cma Beta : " << beta << endl;
            cout << "cluster " << k << ", cma gain : " << gain << endl;
            cout << "cluster " << k << ", cma mean : " << mean << endl;
            cout << "cluster " << k << ", cma V :" << V.first << ", V^pi : " << V.second << endl;
            cout << "cluster " << k << ", cma obj : " << obj << endl;
        }

        if (store_cumprob) {
            _R_mu[k] = R_mu;
            _cumprob[k] = cumprob_k;
        }
    }
    if (store_cumprob){
        _opt_B_r = B_r;
        _opt_B_x = B_x;
        _opt_cma = vector<double>(N);
        for (int i = 0; i < N; ++i) { _opt_cma[i] = v[i]; } 
    }

    return obj;
};

// ---------------------------------------------------------------------
OptimalRewardFinder::OptimalRewardFinder(OptimData* data, Config* config){
    _data   = data;
    _config = config;
    int N_r = config->cma_nbvar();
    int N_x = config->cma_nbvar();
    _r_cma = vector<double>(N_r);
    _x_cma = vector<double>(N_x);
    _bigM = config->cma_bigM();
    _bigM_x = 5.;

    _cumprob = vector<vector<double>>(_data->_nb_k);
    _R_mu = vector<vector<double>>(_data->_nb_k);

    // ----- discretization grid for optimization
    for (int i = 0; i < N_r; ++i) {
        _r_cma[i] = i * 1./(N_r - 1); 
    }
    double h = (_data->_x_max - _data->_x_min) / (N_x - 1); 
    for (int i = 0; i < N_x; ++i) {
        _x_cma[i] = _data->_x_min + i *h; 
    }
};

// ----------------------------------------------------------------------
void OptimalRewardFinder::run(int verbose){
    /*
    This function returns the best solution found by the black-box optimization algorithm
    CMA-ES.
    For homogeneous population, we can successfully recover the optimal reward.

    The optimization is performed in a box [-1,1]^N, and we have constructed a bijection 
    between this set and the set of decreasing (discretized) reward function.
    */
    
    cout << "[INFO] " << "\033[0;32mOptimal Reward resolution in progress ... \033[0m" << endl;
    auto start_solve = chrono::high_resolution_clock::now();

    int N_r = _r_cma.size();
    int N_x = _x_cma.size();
    if (_fixed_price) N_x = 0;
    int N   = N_r + N_x;
    
    // ------ Initialization CMA-ES
    double* v_init = new double[N];
    double* lb     = new double[N];
    double* ub     = new double[N];
    double* stddev = new double[N];
    for (int i = 0; i < N_r; ++i) {  
        v_init[i] = _config->cma_init(); 
        lb[i]     = _config->cma_lb(); 
        ub[i]     = _config->cma_ub(); 
        stddev[i] = _config->cma_sigma();
    }
    lb[0] = 0.; // special lb for first variable

    for (int i = N_r; i < N; ++i) {  
        v_init[i] = -1./_bigM_x; 
        lb[i]     = -2./_bigM_x;        // TMP: to generalize
        ub[i]     = -0.8/_bigM_x; 
        stddev[i] = 0.4*_config->cma_sigma();
    }

    CMAES<double> evo;
    double *arFunvals, *const*pop, *xfinal;
    Parameters<double> parameters;
    parameters.init(N, v_init, stddev);
    parameters.stopMaxIter = _config->cma_maxiter();
    arFunvals = evo.init(parameters);

    // Iterate until stop criterion holds
    double best_obj = 1e8;
    int iter = 0;
    while(!evo.testForTermination() && iter < _config->cma_maxiter())
    {
        // Generate lambda new search points, sample population
        pop = evo.samplePopulation(); // Do not change content of pop

        for (int i = 0; i < evo.get(CMAES<double>::PopSize); ++i) {
             while (!is_feasible(pop[i],lb,ub,N)) evo.reSampleSingle(i);
        }

        // evaluate the new search points using fitfun from above
        omp_set_num_threads(_config->get_num_threads()); // parallelization
        int nb_i = evo.get(CMAES<double>::Lambda);
        #pragma omp parallel for
        for (int i = 0; i < nb_i; ++i){
          arFunvals[i] = fitfun(pop[i], N, verbose);
        }
        for (int i = 0; i < nb_i; ++i){
          if (best_obj > arFunvals[i] + 1e-5) {
            if (verbose >= 1) {
                cout << "Iter " << std::setfill(' ') << std::setw(4) << iter  << " | Obj : " << arFunvals[i] << endl;
            }
            best_obj = arFunvals[i];
          }
        }

        // update the search distribution used for sampleDistribution()
        evo.updateDistribution(arFunvals);

        ++ iter;
    }

    int cma_lambda  = (int) evo.get(CMAES<double>::Lambda);     // nb eval per generation
    int cma_nbiter  = (int) evo.get(CMAES<double>::Generation); // nb generation
    int cma_nbeval  = (int) evo.get(CMAES<double>::Eval);       // nb eval = nb generation * lambda
    double cma_time = (double) evo.eigenTimings.totaltime;      // temps dans CMA
    double cma_obj  = (double) evo.get(CMAES<double>::Fitness); // optimal objective

    if (verbose >= 1) {
        cout << "----- Final solution -----" << endl;
        cout << "CMA lambda    : " << cma_lambda << endl;
        cout << "CMA nb iter   : " << cma_nbiter << endl;
        cout << "CMA nb eval   : " << cma_nbeval << endl;
        cout << "CMA time      : " << cma_time << " sec" << endl;
        cout << "CMA objective : " << cma_obj << endl;
        cout << "--------------------------" << endl;
    }

    // get best estimator for the optimum, xmean
    xfinal = evo.getNew(CMAES<double>::XMean); // "XBestEver" might be used as well

    verbose = 2;
    fitfun(xfinal, N, verbose, true); // to store 
    delete[] xfinal;

    double time = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now()
            - start_solve).count()/1000.0;
    cout << "[INFO] " << "\033[0;32mOptimal Reward resolution added in " + to_string(time) + " s\033[0m" << endl;
};

// ------------------------------------------------------------------------------
void OptimalRewardFinder::construct_reward_from_CMA(const double* v, 
                vector<double>& B_r, vector<double>& B_x){
    /* 
    This function constructs reward B(r) using the variables of CMA 
    */

    int N_r = _r_cma.size();
    B_r = vector<double>(N_r);
    B_r[0] = _bigM * v[0];
    for (int i = 1; i < N_r; ++i) B_r[i] = 0.5*(B_r[i-1]-_bigM) + 0.5*(B_r[i-1]+_bigM) * v[i] ; 

    int N_x = _x_cma.size();
    B_x = vector<double>(N_x);
    if (!_fixed_price){
        for (int i = 0; i < N_x; ++i) B_x[i] = _bigM_x * v[N_r + i] * _x_cma[i]; 
        // - _data->_p * _x_cma[i] - 0.3 * (_x_cma[i] - 14.) * (_x_cma[i] - 20.); // tmp: to test non linear  
    }
    else {
        for (int i = 0; i < N_x; ++i) B_x[i] = - _data->_p * _x_cma[i];
    }
};

// -----------------------------------------------------------------------------------------
void OptimalRewardFinder::construct_reward_from_cumprob(const vector<double>& cumprob, 
                const vector<double>& B_r, const vector<double>& B_x, 
                vector<double>& R_mu){
    /* 
    This function constructs the reward as a function of the value x, 
    knowing the distribution (equilibrium) 
    The computation is exact, as the reward B is piecewise affine
    */
 
    R_mu = vector<double>(_data->_nb_x, 0.);

    int j_r = 0;
    int j_x = 0;
    for (int i = 0; i < _data->_nb_x; ++i){
        double F = cumprob[i];
        double x = _data->_range_x[i];
        while (_r_cma[j_r + 1] < F) { ++ j_r; } // it ends since "_r_cma[-1] = 1"
        while (_x_cma[j_x + 1] < x) { ++ j_x; }
        R_mu[i] += B_r[j_r] + (B_r[j_r+1] - B_r[j_r]) / (_r_cma[j_r+1] - _r_cma[j_r]) *(F - _r_cma[j_r]);
        R_mu[i] += B_x[j_x] + (B_x[j_x+1] - B_x[j_x]) / (_x_cma[j_x+1] - _x_cma[j_x]) *(x - _x_cma[j_x]);
    }
};


// ----------------------------------------------------------------------
tuple<double,vector<double>,vector<double>> get_optimal_reward(OptimData* data, int k){
    /* 
    In the case of homogeneous population, the optimal reward is analytically known,
    see theory 
    */

    double m = ( data->_X0[k] - (data->_cr[k] * data->_T / data->_c[k] / 2.) ) 
                    / ( 1 + data->_alpha*data->_T / data->_c[k] );
    vector<double> pts_cumprob;
    vector<double> pts_reward;

    double constant_B = data->_c[k]/data->_T*(pow(data->_X0[k],2) - pow(m,2))
                        - data->_p * data->_X0[k]
                        + pow(data->_p,2) * data->_T / data->_c[k] / 4.;

    for (int i = 0; i < data->_nb_x; ++i) {
        double x = data->_range_x[i];
        pts_cumprob.push_back(0.5 * (1. + erf( (x - m) / (data->_sigma[k] * sqrt(data->_T * 2)) )) );
        pts_reward.push_back(constant_B + x * (data->_p - data->_cr[k] - 2.*data->_alpha*m) );
    }

    //double int_reward = constant_B + m * (data->_p - data->_cr - 2.*data->_alpha*m);
    double total_gain = data->_alpha * (pow(data->_X0[k],2) + pow(m,2)) - constant_B;

    return make_tuple(total_gain,pts_cumprob,pts_reward);

};

#endif