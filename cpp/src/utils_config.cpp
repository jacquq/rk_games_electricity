#ifndef UTILS_CONFIG_CPP
#define UTILS_CONFIG_CPP

#ifdef _WIN64
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include <map>

#include "utils.hpp"

using namespace std;

/* 
========================================================================

Reading tools

========================================================================
*/

void Config::fill_default_values(){
	/*Valeurs par défaut si non renseigné*/

	// ----- Paths par défaut
	paths["data"] 	  = get_path();
	paths["solution"] = get_path();

	// ----- Labels
	labels["fxdpt_method"]  = "SA";
	labels["fxdpt_damping"] = "inv_prop";

	// ----- Options d'affichages
	options["use_tex"]        = false;
	options["cma_fixedprice"] = true;

	// ----- Paramètres d'optimisation
	params["verbose"]     = 0;
	params["num_threads"] = 1;
	params["delta_x"]     = 0.01;
	params["delta_t"]     = 0.01;

	params["mc_nb_simu"]     = 1000;
	params["fxdpt_max_iter"] = 1000;
	params["fxdpt_eps"]      = 0.01;

	params["cma_sigma"]   = 0.1;
	params["cma_nbvar"]   = 20;
	params["cma_maxiter"] = 1000;
	params["cma_bigM"]    = 30.0;
	params["cma_init"]    = 0.8;
	params["cma_lb"]      = 0.5;
	params["cma_ub"]      = 1.0;

};

// ------------------------------------------
string Config::get_path()
{
	/*Pour récuperer le path courant*/
	char buff[FILENAME_MAX];
	GetCurrentDir( buff, FILENAME_MAX );
	string current_working_dir(buff);
	return current_working_dir;
};

// ------------------------------------------
string remove_char(string input, char a){
	string output;
	for (int i = 0; i < input.size(); ++i){
		if (input[i] != a){ output += input[i]; }
	}
	return output;
};

// ------------------------------------------
void Config::parse_configfile(string filepath)
{
	/*Lecture du fichier de config*/
	string m_line;
	fstream file;
	file.open(filepath,ios::in);
	if (file.is_open()){
		while (getline(file, m_line)) 
		{
			if (m_line[0] != '#' && m_line[0] != '[' && m_line.find('=') != string::npos) 
			{
				// argument
				string arg = m_line.substr(0,m_line.find_first_of("="));
				arg = remove_char(arg, ' ');
				arg = remove_char(arg, '\t');

				// valeur associée
				string val = m_line.substr(m_line.find_first_of("=")+1);
				val = val.substr(val.find_first_not_of(" "));
				val = val.substr(0, val.find_last_not_of(" ")+1);
				val = remove_char(val, '"');

				// On remplace la valeur par défaut
				if (labels.find(arg) != labels.end()) 			labels[arg]  = val;
				else if (paths.find(arg) != paths.end()){
					if (val[0] == '~') 							paths[arg]   = get_path() + val.substr(1);
					else 										paths[arg]   = val;
				}
				else if (options.find(arg) != options.end()) 	options[arg] = (val=="true");
				else if (params.find(arg) != params.end())		params[arg]  = stof(val);
			}
		}
		cout << "[INFO] " <<  "\033[0;32mConfig file loaded\033[0m" << endl;
	}
	else {
		cout << "[INFO] " <<  "\033[0;33mConfig file not open\033[0m" << endl;
	}
};

// ------------------------------------------
Config::Config(int argc, char * argv[]){
	fill_default_values();
	vector<string> string_argv;
    string_argv.assign(argv, argv + argc);
    path_folder = getExecutableDir();
	if(argc > 1) { path_config = string_argv[1]; }
	else { path_config = path_folder + "/config.ini"; }// default path

	parse_configfile(path_config);
};


#endif