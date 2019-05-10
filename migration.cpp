// Collective migration when resources vary
// Bram Kuijper & Simon Evans
// 2019
//
#include <iostream>

// random number generation
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


// various functions, such as unique filename creation
#include "auxiliary.h"

#define DEBUG

// standard namespace
using namespace std;

// number of individuals in population
const int N = 5000;

// number of generations
long int number_generations = 10;

struct Individual {
    
    double resources;

    // elevation (baseline leaving rate)
    double theta_a[2];

    // reaction norm, slope on the number of individuals in the pool
    double theta_b[2];
};

// wintering ground
Individual WinterPop[N];
Individual SummerPop[N];

string filename("sim_migration")
string filename_new(create_filename(filename));





#ifdef DEBUG
    // write stuff
#endif


