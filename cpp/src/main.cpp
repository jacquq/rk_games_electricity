#include <vector>
#include "math.h"

#include "data.hpp"
#include "utils.hpp"
#include "methods.hpp"

#include "math.h"
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

using namespace std;

/* =====================================

This work is based on two papers:

    [1] Erhan Bayraktar and Yuchong Zhang. “A rank-based mean field game in the strong formulation”. 
        In: Electronic Communications in Probability 21 (2016), pp. 1–12.
    [2] Erhan Bayraktar and Yuchong Zhang. “Terminal Ranking Games”. 
        In: Mathematics of Operations Research Forthcoming (2019).

References of these two papers are made in the code. When we refer to results from the new papers,
we use '[*]'.

@Quentin Jacquet

======================================= */

int main(int argc, char *argv[]) {

    // ======================== Read config & data =================================
    Config* config  = new Config(argc, argv); // a file 'config.ini' can be provided
    OptimData* data = new OptimData(config);  // a file 'data.dat' is expected, see 'data.cpp'
    // =============================================================================


    // ======================== Optimization of the reward =========================
    OptimalRewardFinder* orf = new OptimalRewardFinder(data,config);    // declare the dedicated class, see 'retailer.cpp'
    orf->_fixed_price = config->cma_fixedprice();                       // have a fixed dependance on price "-px"
    orf->run();                                                         // run the search of the optimal reward (using CMAES)
    orf->plot_reward(true, config->use_tex());                          // plot (in matplotlib) the solution
    // =============================================================================

    // ============================= Plots =========================================
    map<string,string> keywords;
    MonteCarlo mc(data, config->get_delta_t(),  config->mc_nb_simu());
    for (int k = 0; k < data->_nb_k; ++k) {
        
        // ----- plot optimal cdf
        keywords["color"] = "green";
        keywords["label"] = "PDF with B*";
        plot_cumprob(data->_range_x, orf->_cumprob[k],keywords,false);
        
        // ----- plot pdf with B = 0
        vector<double> B_r(orf->_r_cma.size(), 0.);
        vector<double> cumprob;
        double mean;
        orf->get_equilibrium(k, B_r, mean, cumprob);
        keywords["color"] = "blue";
        keywords["label"] = "PDF with B=0";
        plot_cumprob(data->_range_x, cumprob, keywords,false);
        
        // ----- plot pdf with B=0 and p=0
        double mem_p = data->_p;
        data->_p = 0.;
        orf->get_equilibrium(k, B_r, mean ,cumprob);
        keywords["color"] = "brown";
        keywords["label"] = "PDF with B=0 and p = 0";
        plot_cumprob(data->_range_x, cumprob, keywords);
        data->_p = mem_p;

        // ----- compute numerical fixed point
        vector<double> its;
        vector<vector<double>> list_its;
        vector<string> dampings = {"inv_prop", "KM", "WOD"};
        config->labels["fxdpt_method"] = "SA";
        for (vector<string>::iterator it = dampings.begin(); it != dampings.end(); ++it){
            config->labels["fxdpt_damping"] = *it;
            its = orf->compute_equilibrium(k, orf->_opt_B_r, orf->_opt_B_x, mean, cumprob); // solve the fixed-point equation with semi-analytic formula
            list_its.push_back(its);
        }
        vector<string> legend = {"Damping 1/2", "Damping 1/(k+1)", "No damping"};
        plot_its_equilibrium(legend, list_its);

        // ---- Plot trajectories with B*
        mc.run(k,orf->_R_mu[k], 1, k, true);
        mc.plot_trajectories(20,false);
        keywords.clear();
        keywords["linewidth"] = "4";
        keywords["label"]     = "Mean of consumption with B*";
        keywords["color"]     = "green";
        mc.plot_mean(keywords);

        mc.plot_effortcost(20);

        mc.plot_brownian(20);

        // ---- Plot cost of effort
    }
    // trajectories with B=0
    int size_b = config->cma_nbvar();
    double* b = new double[size_b];
    b[0] = 0.;
    for (int i = 1; i < size_b; ++i) b[i] = 1.; // cosntant reward
    orf->fitfun(b,size_b,1, true);
    for (int k = 0; k < data->_nb_k; ++k) {
        mc.run(k,orf->_R_mu[k], 1, k, true);
        mc.plot_trajectories(20,false);
        keywords["label"] = "Mean of consumption B=0";
        keywords["color"] = "blue";
        mc.plot_mean(keywords);    
    }
    // =============================================================================

    return 0;
};
