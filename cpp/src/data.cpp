#ifndef CLASS_PARAMS
#define CLASS_PARAMS

#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <algorithm>
#include "math.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#include "data.hpp"
#include "utils.hpp"

using namespace std;

/* ----------------------------

This class provides a parser for the data input file,
and provides functions that precomputes interesting 
quantities for 'equilibrium.cpp' and 'retailer.cpp'

---------------------------- */

// --------------------------------------------------------
OptimData::OptimData(Config* config, char delim){
    /*
    Parse the input file at the construction of the class
    */

    // ----- Load data file
    string datapath = config->get_data_path();
    vector<vector<double>> data;
    string m_line;
    fstream file;
    file.open(datapath,ios::in);
    if (file.is_open()){
        while (getline(file, m_line)) 
        {
            if (m_line[0] != '#') {
                vector<double> row = vector<double>();
                string item = "";
                int pos = 0;
                while (pos < m_line.size())
                {
                    char c = m_line[pos];
                    if (c != delim){ item.push_back(c); }
                    else {                            
                        replace(item.begin(), item.end(), ',', '.');
                        row.push_back(stod(item));
                        item = "";
                    }
                    ++ pos;
                }
                if (item.size() > 0){
                    replace(item.begin(), item.end(), ',', '.');
                    row.push_back(stod(item));
                }
                if (row.size() > 0){
                    data.push_back(row);
                }
            }
        }
        cout << "[INFO] " << "\033[0;32mFile '"+datapath.substr(datapath.find_last_of("/\\")+1)+"' loaded\033[0m" << endl;
        file.close();       

        // ============== parsing
        vector<vector<double>>::iterator it = data.begin();
        _nb_k  = (*it)[0]; ++it;
        for (int i=0; i < _nb_k; ++i){ _sigma.push_back( (*it)[0] ); ++it; }
        _T     = (*it)[0]; ++it;
        _p     = (*it)[0]; ++it;
        for (int i=0; i < _nb_k; ++i){ 
            _c.push_back( (*it)[0] ); 
            if (it->size() > 1) _c_der.push_back( (*it)[1] );
            else                _c_der.push_back(     0.   );
            ++it; 
        }
        for (int i=0; i < _nb_k; ++i){ _X0.push_back( (*it)[0] ); ++it; }
        _alpha = (*it)[0]; ++it;
        for (int i=0; i < _nb_k; ++i){ _cr.push_back( (*it)[0] ); ++it; }
        

        // -----construction of grids on x
        _x_min = (*it)[0]; 
        _x_max = (*it)[1]; ++it; 
        _nb_x = (int) floor((_x_max - _x_min) / 2./config->get_delta_x());
        _nb_x = 2 * _nb_x + 1; // impair pour Simpson
        _range_x = vector<double>(_nb_x);
        double h = (_x_max - _x_min) / (_nb_x - 1); 
        for (int j = 0; j < _nb_x; ++j){ _range_x[j] = _x_min + j*h ; } 

        // ----- reward given in config
    	//while (it != data.end()) { _R_mu.push_back( (*it)[0] ); ++ it; }
        //_nb_mu = (int) _R_mu.size();	
    }
    else {
        cout << "[INFO] " <<  "\033[0;33mFile '"+datapath.substr(datapath.find_last_of("/\\")+1)+"' not open\033[0m" << endl;
    }
};

// ---------------------------------------------------------------------------------
/*double OptimData::get_Reward_pure_ranked(double rank){
    
    //This function returns the value of the additionnal reward B at 'rank', where
    //B is defined as a linear interpolation of the points given in the input file.
    
    
    // Interpolate R(rank)
    double h = 1. / (_nb_mu - 1);
    double idx = rank / h;
    int idx_lw = (int) idx;
    //cout << idx_lw << ", " << _R_mu.size() << endl;
    return _R_mu[idx_lw] * (idx - idx_lw) + _R_mu[idx_lw+1] * (idx_lw + 1 - idx);
};*/

// ---------------------------------------------------------------------------------
double OptimData::get_Reward(int k, const vector<double>& R_mu, double x){
    /*
    This function computes the value of the total reward ('B - px') for the cluster 'k' at
    the terminal value 'x' of the processus 'X_t'. 
    This corresponds the quantity 'R_\mu(x) = B(q_\mu(x)) - p x', see [*].
    */
    
    // first, interpolate rank(x), the cumulative probability is known through 'cumprob',
    //    and the grid in 'x' is supposed to be regular between '_x_min' and '_x_max'.
    double h = (_x_max - _x_min) / (_nb_x - 1); 
    double idx = (x - _x_min) / h;
    int idx_lw = (int) idx;
    if (idx_lw < 0) return R_mu[0];
    if (idx_lw > R_mu.size() - 2) return R_mu[R_mu.size() - 1];
    
    double r = R_mu[idx_lw] * (idx - idx_lw) + R_mu[idx_lw+1] * (idx_lw + 1 - idx);
    return r;
};

// ---------------------------------------------------------------------------------
void OptimData::get_u_and_ux(int k, double t, double xt, const vector<double>& R_mu, double& u, double& u_x){
    /*
    This function is the implementation of equations (3.3) and (3.4) of [1].
    We respect as much as we can the notations of the paper.
    */
    double min_z = -7.;
    double max_z = 7.;
    int size_z = _nb_x;
    double h = (max_z - min_z) / (size_z - 1);
    double sqrt_T_t = sqrt(_T - t);
    double z_i, u_i;

    u = 0.;
    u_x = 0.;

    // discretization of the integral as a Riemann sum
    for (int i = 0; i < size_z-1; ++i) {
        z_i = min_z + i*h;
        u_i = exp(   1./(2 * _sigma[k]*_sigma[k] * (_c[k]+_c_der[k]*t)) * 
                            get_Reward(k, R_mu, xt + _sigma[k] * sqrt_T_t * z_i)
                            - z_i * z_i / 2
                        );
        u += u_i;
        u_x += u_i * z_i; 
    }
    u   /= sqrt(2*M_PI);
    u_x /= (sqrt(2*M_PI) * _sigma[k] * sqrt_T_t);
};

// ---------------------------------------------------------------------------------
double OptimData::get_optimal_effort(int k, double t, double xt, const vector<double>& R_mu){
    /* 
    This function is a direct implementation of the optimal effort 'a* = v_x / (2c)'. 
    */
    double u, u_x;
    get_u_and_ux(k, t, xt, R_mu, u, u_x); 
    return u_x / u * pow(_sigma[k], 2);
};

// ---------------------------------------------------------------------------------
double OptimData::get_Wass_dist(const vector<double>& cumprob_1, const vector<double>& cumprob_2) {
    /*
    This function computes the Wasserstein distance between two distributions, given by their
    densition function.
    As we are in R, The distance W_1 can be computed using cumulative distribution, see [*].
    */

    double h = (_x_max - _x_min) / (_nb_x - 1);
    double W1 = 0.;
    // We discretize the integrale as a Riemann sum
    for (int i = 0; i < _nb_x; ++i) {
        W1 += abs(cumprob_1[i] - cumprob_2[i]);
    }
    return W1 * h;
    
};

// --------------------------------------------------------------------------------------

vector<double> OptimData::get_new_iterate_cdf(int k, const vector<double>& R_mu){
    double x,f;
    double int_f = 0.;
    vector<double> new_cumprob(_nb_x);
    for (int i = 0; i < _nb_x; ++i){

        int_f += exp(-0.5*pow(_range_x[i] - _X0[k],2) / (_T * pow(_sigma[k],2))
                            + R_mu[i]/(2.*_c[k] * pow(_sigma[k],2)));
        new_cumprob[i] = int_f;
    }
    for (int i = 0; i < _nb_x; ++i){
        new_cumprob[i] /= int_f;
    }

    return new_cumprob;    
};

// =====================================================================================

double OptimData::retailer_Beta(int k, const vector<double>& b_rangex){
    /*
    This function computes the quantity '\Beta' defines in Proposition 2.1 of [2], in the
    specific case where the reward is 'R = B - px', see [*].
    */

    vector<double> integs;
    for (int i = 0; i < _nb_x; ++i){
        integs.push_back( exp(    -0.5*pow(_range_x[i] - _X0[k],2) / (_T * pow(_sigma[k],2))
                            + (b_rangex[i])/(2.*_c[k] * pow(_sigma[k],2))
                        ));
    }
    double int_Beta = integ(integs, _x_min, _x_max) / sqrt(2.*M_PI*_T) / _sigma[k];

    return int_Beta;
};

// --------------------------------------------------------------------------------------
double OptimData::retailer_intReward_xpart(const vector<double>& x_cma, const vector<double>& B_x, const vector<double>& cumprob){
    /*
    This function computes the integrale of Reward, providing a range of ranks 'rank' 
    and the correspondingvalues of reward 'b'.

    TODO: use analytic formula for piecewise affine B_x.
    */

    double val = 0.;
    int j_x = 0;
    for (int i = 0; i < _nb_x - 1; ++i){
        double x = _range_x[i];
        while (x_cma[j_x+1] < x){ ++ j_x; } 
        double B = B_x[j_x]
                + (B_x[j_x+1] - B_x[j_x]) / (x_cma[j_x+1] - x_cma[j_x]) * (x - x_cma[j_x]);
        //cout << x << ", " << B << endl;
        val += B * (cumprob[i+1]-cumprob[i]); 
    }

    return val;
};

// --------------------------------------------------------------------------------------
double OptimData::retailer_intReward_rankpart(const vector<double>& r_cma, const vector<double>& reward){
    /*
    This function computes the integrale of Reward, providing a range of ranks 'rank' 
    and the correspondingvalues of reward 'b'.
    Trapezian method is exact for piecewise affine reward.
    */

    double int_b = 0.;
    for (int i = 0; i < r_cma.size()-1; ++i) {
        int_b += 0.5*(r_cma[i+1] - r_cma[i])*(reward[i] + reward[i+1]);
    }
    return int_b;
};

// --------------------------------------------------------------------------------------
double OptimData::retailer_gain(int k, double mean){
    /*
    This function computes the gain for the retailer defines as 'g(m) = s(m) - c_r*m'.
    */

    return _alpha * (pow(_X0[k],2) - pow(mean,2) ) - _cr[k] * mean;
};

// --------------------------------------------------------------------------------------
pair<double,double> OptimData::retailer_constraintV0(int k, double int_Beta){
    /*
    This function computes the utility constraint of customers from cluster k, where the utility is
    defined in Proposition 2.1 in [2].
    */
    
    return make_pair(2.*_c[k] * pow(_sigma[k],2)*log(int_Beta), // -V(R)
                - _p * _X0[k] + pow(_p,2) / _c[k] / pow(_sigma[k],2) / 4.); // + V_0
};

#endif
