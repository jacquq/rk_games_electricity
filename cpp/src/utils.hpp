#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <tuple>
#include <vector>
#include <map>

using namespace std;

/* 
========================================================================

Outils pour la lecture du fichier de config:
Le fichier de config .ini est au format. Celui est parsé manuellement,
et seuls les champs pré-définis.

========================================================================
*/


struct Config 
{
	map<string, string> paths;
	map<string, string> labels;
	map<string, bool> options;
	map<string, double> params;
	string path_folder, path_config;

	Config(int argc, char *argv []);
	~Config(){};

	string get_path();
	string PrettyHeader();
	void parse_configfile(string filepath);
	void fill_default_values();

	// ----- Getters
	string get_data_path(){return paths["data"];};
	string get_data_name(){return paths["data"].substr(paths["data"].find_last_of("/")+1);};
	string get_solution_path(){return paths["solution"];};
	string get_solution_name(){return paths["solution"].substr(paths["solution"].find_last_of("/")+1);};

	string fxdpt_method(){return labels["fxdpt_method"];};
	string fxdpt_damping(){return labels["fxdpt_damping"];};
	
	bool use_tex(){return options["use_tex"];};
	bool cma_fixedprice(){return options["cma_fixedprice"];};

	int get_verbose(){return (int) params["verbose"];};
	int get_num_threads(){return (int) params["num_threads"];};
	double get_delta_x(){return params["delta_x"];};
	double get_delta_t(){return params["delta_t"];};
	
	int fxdpt_max_iter(){return (int) params["fxdpt_max_iter"];};
	double fxdpt_eps(){return params["fxdpt_eps"];};
	int mc_nb_simu(){return (int) params["mc_nb_simu"];};

	double cma_sigma(){return params["cma_sigma"];};
	int cma_nbvar(){return (int) params["cma_nbvar"];};
	int cma_maxiter(){return (int) params["cma_maxiter"];};
	double cma_bigM(){return params["cma_bigM"];};
	double cma_init(){return params["cma_init"];};
	double cma_lb(){return params["cma_lb"];};
	double cma_ub(){return params["cma_ub"];};
};


/* 
========================================================================

Outils pour récupérer le path de l'exe:
Code récupéré sur Internet tel quel

========================================================================
*/
string getExecutablePath();
string getExecutableDir();

/* 
========================================================================

Outils pour faire une jolie barre de progression

========================================================================
*/
void printProgress(double percentage);


float my_logf (float);
float erfinv (float);


/* 
========================================================================

Outils pour l'intégration (Simpson)

========================================================================
*/
double integ(const vector<double>& f, double a, double b, string method = "S");
double integ(const vector<double>& x, const vector<double>& f);

#endif