#ifndef PLOT_FUNCTION
#define PLOT_FUNCTION

#include <string>
#include <tuple>
#include <vector>
#include <map>
#include <math.h>

#include "matplotlibcpp.h"

#include "data.hpp"
#include "methods.hpp"

using namespace std;
namespace plt = matplotlibcpp;

map<string,string> get_rcparams(bool usetex){
    map<string,string> rc_params;
    rc_params["font.size"]   = "12";
    rc_params["font.family"] = "serif";
    rc_params["font.serif"]  = "Helvetica"; // Helvetica or Times;
    if (usetex) rc_params["text.usetex"] = "True";
    return rc_params;
};


// -------------------------------------------------------------------------------
void OptimalRewardFinder::plot_reward(bool plot_opt_reward, bool usetex, bool show){
    cout << usetex << endl;
    plt::rcparams(get_rcparams(usetex));
    if (usetex) plt::named_plot("Reward CMA $B^*$",_r_cma,_opt_B_r, "dg--");
    else        plt::named_plot("Reward CMA B*",_r_cma,_opt_B_r, "dg--");
    
    // ----- Compair with optimal one
    if (plot_opt_reward) {
        vector<string> colors = {"r","purple"};
        for (int k = 0; k < _data->_nb_k; ++k){
            auto t = get_optimal_reward(_data,k);
            plt::named_plot("Optimal reward cluster "+to_string(k+1), get<1>(t), get<2>(t),colors[k]);
        }
    }
    plt::ylim(-_bigM,_bigM);
    plt::xlim(0.,1.);
    if (usetex) plt::xlabel("Ranking $r$");
    else        plt::xlabel("Ranking");
    plt::ylabel("Reward (1e2€)");
    plt::legend();
    plt::grid(true);
    if (show) plt::show();

    // ----- Deuxième graph
    if (!_fixed_price){
        plt::rcparams(get_rcparams(usetex));
        if (usetex) plt::named_plot("C(x)",_x_cma,_opt_B_x, "db--");
        else        plt::named_plot("C(x)",_x_cma,_opt_B_x, "db--");
        
        //plt::ylim(_bigM_x*_data->_x_min,_bigM_x*_data->_x_max);
        if (usetex) plt::xlabel("Terminal value $x$");
        else        plt::xlabel("Terminal value");
        plt::ylabel("Reward (1e2€)");
        plt::legend();
        plt::grid(true);
        if (show) plt::show();
    }


};

// ----------------------------------------------
void plot_its_equilibrium(vector<string> legend, vector<vector<double>> list_its){
    for (int i = 0; i < list_its.size(); ++i){
        vector<double> log_y(list_its[i].size());
        for (int j = 0; j < list_its[i].size(); ++j){
            log_y[j] = log10(list_its[i][j]);
        }
        plt::named_plot(legend[i],log_y);
    }
    plt::legend();
    plt::grid(true);
    plt::show();
};

// ----------------------------------------------
void MonteCarlo::plot_trajectories(int nb_traj, bool show){

    int nb_points = (int) floor(_data->_T / _delta_t);
    double h = _data->_T / nb_points;
    
    // ----- times
    vector<double> pts_t;
    for (int n = 0; n < nb_points; ++n) {
        pts_t.push_back(h*n);
    }

    // ----- trajectories
    for (int i = 0; i < nb_traj; ++i){
        plt::plot(pts_t, _store_trajectory[i]);
    }
    
    if (show) {
        plt::legend();
        plt::grid(true);
        plt::xlabel("Time");
        plt::ylabel("Forecasted consumption (MWh)");
        plt::show();
    }
};

// ----------------------------------------------
void MonteCarlo::plot_effortcost(int nb_traj, bool show){

    int nb_points = (int) floor(_data->_T / _delta_t);
    double h = _data->_T / nb_points;
    
    // ----- times
    vector<double> pts_t;
    for (int n = 0; n < nb_points; ++n) {
        pts_t.push_back(h*n);
    }

    // ----- cost of effort
    for (int i = 0; i < nb_traj; ++i){
        plt::plot(pts_t, _store_effortcost[i]);
    }
    
    if (show) {
        plt::grid(true);
        plt::xlabel("Time");
        plt::ylabel("Cost of effort (1e2€)");
        plt::show();
    }
};

// ----------------------------------------------
void MonteCarlo::plot_brownian(int nb_traj, bool show){

    int nb_points = (int) floor(_data->_T / _delta_t);
    double h = _data->_T / nb_points;
    
    // ----- times
    vector<double> pts_t;
    for (int n = 0; n < nb_points; ++n) {
        pts_t.push_back(h*n);
    }

    // ----- cost of effort
    for (int i = 0; i < nb_traj; ++i){
        plt::plot(pts_t, _store_brownian[i]);
    }
    
    if (show) {
        plt::grid(true);
        plt::xlabel("Time");
        plt::ylabel("Brownian deviation (MWh");
        plt::show();
    }
};

// -----------------------------------------------------------------
void MonteCarlo::plot_mean(map<string, string> keywords, bool show){

    int nb_points = (int) floor(_data->_T / _delta_t);
    double h = _data->_T / nb_points;
    
    // ----- times
    vector<double> pts_t;
    for (int n = 0; n < nb_points; ++n) {
        pts_t.push_back(h*n);
    }

    // ----- mean value
    vector<double> means(nb_points);
    for (int n = 0; n < nb_points; ++n){
        double m = 0.;
        for (int i = 0; i < _nb_simu; ++i) {
            m += _store_trajectory[i][n];
        }
        means[n] = m/_nb_simu;
    }
    plt::plot(pts_t, means, keywords);
    
    if (show) {
        plt::legend();
        plt::grid(true);
        plt::xlabel("Time");
        plt::ylabel("Forecasted consumption (MWh)");
        plt::show();
    }
};

// ------------------------------------------------------------------------------------------
void plot_cumprob(const vector<double>& range_x, const vector<double>& cumprob, const map<string,string>& keywords, bool show){
    plt::plot(range_x, cumprob, keywords);
    if (show) {
        plt::legend();
        plt::grid(true);
        plt::xlabel("Forecasted consumption (MWh)");
        plt::ylabel("Cum. density function");
        plt::show();
    }
};

#endif