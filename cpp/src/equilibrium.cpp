#ifndef EQUILIBRIUM
#define EQUILIBRIUM

#include <string>
#include <tuple>
#include <vector>
#include <map>
#include <math.h>
#include <chrono>
#include <numeric>

#include "asa241.hpp"

#include "data.hpp"
#include "methods.hpp"

using namespace std;
using namespace chrono;

/* ----------------------------

This class provides a fixed-point solver (with an iterative scheme).
Moreover, in the specific case 'R = B - px', functions are provided to computes
the analytically known equilibrium, see [1] and [*].

---------------------------- */


// ---------------------------------------------------------------------------------
vector<double> OptimalRewardFinder::compute_equilibrium(int k, const vector<double>& B_r, const vector<double>& B_x, 
    double& mean, vector<double>& cumprob, int verbose) {
    /*
    This function a fixed-point solver (with an iterative scheme).
    
        ---> cum_prob -|-> optimal effort ---> MC trajectories -|
                       '----------------------------------------'

    Terminal condition is a sufficiently small Wasserstein distance between two iterates.
    Krasnoselskii-Mann damping is used, see [*].
    */

    if (verbose >= 1){
        cout << "[INFO] " << "\033[0;32mFixed-point resolution in progress ... \033[0m" << endl;
    }
    auto start_solve = chrono::high_resolution_clock::now();

    vector<double> cumprob_next;

    // Create Monte-Carlo simulator
    MonteCarlo mc(_data, _config->get_delta_t(), _config->mc_nb_simu());

    vector<double> its;

    // Initial _cumprob (todo: find best one, for now linear one)
    cumprob.clear();
    for (int i = 0; i < _data->_nb_x; ++i) {
        cumprob.push_back( i*1./(_data->_nb_x - 1) );
    }

    double Wass_dist = _data->_x_max - _data->_x_min;
    int nb_iter = 0;
    double eps = _config->fxdpt_eps();
    while (Wass_dist > eps && nb_iter < _config->fxdpt_max_iter()) {

        if (verbose >= 1){
            printProgress(ceil( log(Wass_dist/(_data->_x_max - _data->_x_min))
                                /log(eps/(_data->_x_max - _data->_x_min))*1e3
                            )/10); // to display
        }

        vector<double> R_mu;
        construct_reward_from_cumprob(cumprob, B_r, B_x, R_mu);
        if (_config->fxdpt_method() == "MC") {
            // run montecarlo with best effort
            mc.run(k, R_mu, _config->get_num_threads(),0,true);
            mc.fill_new_cumprob(cumprob_next);
        }
        else {
            // We directly known the next cum probs
            cumprob_next = _data->get_new_iterate_cdf(k, R_mu);
        }

        // update
        ++ nb_iter;
        double w;
        string weights = _config->fxdpt_damping();
        if (weights == "KM")            w = 0.5;
        else if (weights == "inv_prop") w = 1./(1+nb_iter);
        else if (weights == "WOD")      w = 1.;
        for (int i = 0; i < _data->_nb_x; ++i){
            cumprob_next[i] = w*cumprob_next[i] + (1. - w)*cumprob[i]; // stabilization
        }
         // compute distance between the two iterates
        Wass_dist = _data->get_Wass_dist(cumprob, cumprob_next);
        if (verbose >= 2)  cout << "Cluster " << k << ", Wass dist : " << Wass_dist << endl;
        cumprob = cumprob_next;

        its.push_back(Wass_dist);
    }
    if (verbose >= 2) cout << endl;

    // compute the mean value of the process as the integral of (1-F(x)) between 0 in +oo
    mean = 0.;
    double h = (_data->_x_max - _data->_x_min) / (_data->_nb_x - 1); 
    double x;

    for (int i = 1; i < _data->_nb_x; ++i){
        x = _data->_x_min + i*h;
        mean += x * (cumprob[i] - cumprob[i-1]);
    }

    double time = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now()
            - start_solve).count()/1000.0;
    if (verbose >= 1){
        cout << "[INFO] " << "\033[0;32mFixed-point resolution added in " + to_string(time) + " s\033[0m" << endl;
    }

    return its;
};


// -------------------------------------------------------------
void OptimalRewardFinder::get_equilibrium(int k, const vector<double>& reward,
    double& mean, vector<double>& cumprob){
    
    /*
    This function computes the equilibrium in the specific case 'B - px'.
    This is a shifted quantile function of the pure-ranked case described in 
    Theorem 3.2 of [2]. See [*] for more details.

    Return the mean and the cumultaive distribution
    */

    // ----- precompute integrals on coarse grid (exact calculation)
    int size_r = _r_cma.size();
    double r_i, B_i, r_ii, B_ii, h_i;
    vector<double> integs(size_r - 1);
    for (int i = 0; i < size_r-1; ++i) {
        B_i = reward[i];
        r_i = _r_cma[i];
        B_ii = reward[i+1];
        r_ii = _r_cma[i+1];
        h_i = r_ii - r_i;
        if (abs(B_ii - B_i) > 1e-5){
            integs[i] = 2.*_data->_c[k]*pow(_data->_sigma[k],2)*
                        h_i / (B_i - B_ii) * (  exp(- B_ii/2./_data->_c[k]/pow(_data->_sigma[k],2))
                                                - exp(- B_i/2./_data->_c[k]/pow(_data->_sigma[k],2))
                                                );
        }
        else integs[i] = h_i * exp(- B_i/2./_data->_c[k]/pow(_data->_sigma[k],2));
    }

    double tot_integ = accumulate(integs.begin(), integs.end(), 0.);

    // ----- compute cumprob on range_x (exact calculation)
    cumprob    = vector<double>(_data->_nb_x);
    vector<double> cumprob_bar, range_pos_x;
    double cum_r, q, m;
    int idx_coarse_r;

    double r = 0.;
    for (int i = 0; i < _data->_nb_x; ++i) {

        // achieve a binary search on rank to find r such that ''q(r) = _range_x[i]''
        double m_r = r;
        double M_r = 1.;
        while (M_r - m_r > 1e-8 ){
            r = 0.5*(m_r + M_r);

            // compute integrale from r_0 = 0 to r, using the precomputed integrales (Chasles relation) 
            cum_r = 0.;
            idx_coarse_r = 1;
            r_i  = _r_cma[0];
            B_i  = reward[0];
            B_ii = reward[1];
            r_ii = _r_cma[1];
            h_i  = r_ii - r_i;
            while (_r_cma[idx_coarse_r] < r) {
                B_i    = B_ii;
                r_i    = _r_cma[idx_coarse_r];
                cum_r += integs[idx_coarse_r-1];
                ++idx_coarse_r;
                r_ii = _r_cma[idx_coarse_r];
                B_ii = reward[idx_coarse_r];
                h_i  = r_ii - r_i;
            }

            // compute the last part (integrale from r_i to r)
            if (abs(B_ii - B_i) > 1e-5){
                cum_r += 2.*_data->_c[k]*pow(_data->_sigma[k],2)*
                        h_i / (B_i - B_ii) * (  exp(- (B_i + (B_ii - B_i)/h_i*(r - r_i))
                                                        /2./_data->_c[k]/pow(_data->_sigma[k],2))
                                                   - exp(- B_i/2./_data->_c[k]/pow(_data->_sigma[k],2))
                                                   );
            }
            else cum_r += (r - r_i)* exp(- B_i/2./_data->_c[k]/pow(_data->_sigma[k],2));

            // compute the quantile value
            // using equation (3.2) of [2], the inverse of the cumulative distribution of N(0,1) is 
            //   computed using 'r8_normal_01_cdf_inverse' of 'asa241.cpp'. 
            q = _data->_X0[k] - _data->_p*_data->_T/2./_data->_c[k] 
                + _data->_sigma[k] * sqrt(_data->_T) 
                    * r8_normal_01_cdf_inverse (cum_r / tot_integ);
            
            if (q <= _data->_range_x[i]+1e-6)  m_r = r;
            else                        M_r = r;
        }
        r = min(max(r, 1e-8),1. - 1e-8); // pour la stabilité
        cumprob[i]     = r;
        if (q >= 0) {
            range_pos_x.push_back(q);
            cumprob_bar.push_back(1. - r);
        }
        //else cumprob_bar.push_back(-r);
    }

    //map<string,string> keywords;
    //plot_cumprob(data->_range_x, cumprob, keywords,false);
    //plot_cumprob(range_pos_x, cumprob_bar, keywords);

    // compute the mean value, using the relation between mean and cumprob, see wiki
    mean = 0.;
    mean += integ(range_pos_x, cumprob_bar);
    mean += range_pos_x[0]; // caution, this is for the case where the first point is not 0

};

#endif