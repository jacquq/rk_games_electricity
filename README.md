Ranking Games for Energy Savings
============

Folders
------------

The code contains:

  * `src` contains all the C++ source files
  * `include` contains some dependencies (mainly the blackbox solver and the plotting libraries)


Installation
------------

To compile the source code, Cmake is recommended since we
provide a `CMakeLists.txt` file along with a bash script to facilitate the job.

### Dependencies: ###
The code depends on 
  * `Eigen3` (mendatory),
  * `OpenMP` (recommended)
  * `Matplotlib` (recommended)

Config and Data
------------

The code contains a config.ini file and the option can be modified directly using this file.
In particular, it contains the path to the data file. Examples of data (used in the article [[1]](#1) )
are also provided.


## References
<a id="1">[1]</a> 
Clémence Alasseur, Erhan Bayraktar, Roxana Dumitrescu, Quentin Jacquet,
A Rank-Based Reward between a Principal and a Field of Agents: Application to Energy Savings, arXiv preprint, 2022

