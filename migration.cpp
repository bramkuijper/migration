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

// parameters & variables:

// number of individuals in population
const int N = 5000;

// number of generations
long int number_generations = 10;

// initial values for a and b
double init_a = 0.0;
double init_b = 0.0;


// keep track of the current number of 
// individuals signaling to disperse
int NStaging = 0;

string filename("sim_migration");
string filename_new(create_filename(filename));



struct Individual {
    
    double resources;

    // resource reaction norm (determines entry into staging pool)
    //
    // elevation (baseline leaving rate) 
    double theta_a[2];

    // reaction norm, slope on the amount of resources
    double theta_b[2];

    // collective dispersal reaction norm 
    // determines migration dependent on number of individuals
    double phi_a[2];
    double theta_a[2];
};

// wintering ground
Individual WinterPop[N];
Individual StagingPool[N];
Individual FlockPool[N];

Individual SummerPop[N];


// get parameters from the command line when 
// running the executable file
void init_arguments()
{
    init_a = // some command line argument
    init_b = // some command line argument
}

// initialize the population
void init_population()
{
    // loop through all individuals in the wintering ground
    // and give them values
    for (int i = 0; i < N; ++i)
    {
        WinterPop[i].resources = 0.0;

        for (int j = 0; j < 2; ++j)
        {
            WinterPop[i].theta_a[j] = init_a;
            WinterPop[i].theta_b[j] = init_b;
        }
    }
}

// the dynamics of the population at the wintering ground
void winter_dynamics()
{
    // individuals forage
    // individuals accumulate resources
    // individuals make dispersal decisions

    // probability of encountering a good resource
    double pgood = pgood_init - decay_good * t;

    // set lower boundary to the probability
    if (pgood <= 0)
    {
        pgood = 0;
    }

    for (int i = 0; i < NWinter; ++i)
    {
        if (gsl_rng_uniform(rng_r) < pgood) // good resource chosen
        {
            WinterPop[i].resources += rgood;
        }
        else
        {
            WinterPop[i].resources += rbad;
        }
    
    } // ok resource dynamic done

    for (int i = 0; i < NStaging; ++i)
    {
        if (gsl_rng_uniform(rng_r) < pgood) // good resource chosen
        {
            StagingPool[i].resources += rgood;
        }
        else
        {
            StagingPool[i].resources += rbad;
        }
    }

    assert(NWinter <= N);
    assert(NWinter >= 0);

    // signal to disperse
    for (int i = 0; i < NWinter; ++i) 
    {
        // reaction norm dependent on resources
        // resulting in signaling a willingness to disperse
        // => go to the staging level
        psignal = 0.5 * (WinterPop[i].theta_a[0] + WinterPop[i].theta_a[1])
            + 0.5 * (WinterPop[i].theta_b[0] + WinterPop[i].theta_b[1]) * WinterPop[i].resources;

        // does individual want to signal?
        if (gsl_rng_uniform(rng_r) < psignal)
        {
            // add individual to the staging pool
            StagingPool[NStaging] = WinterPop[i];
            ++NStaging; // increment the number of individuals in the staging pool

            assert(NStaging <= N);
            assert(NStaging >= 0);

            // delete this individual from the winter population
            WinterPop[i] = WinterPop[NWinter - 1];

            // decrement the number of individuals in the winter population
            --NWinter;
            --i;
        }
    }

    // store current number of individuals at the breeding ground
    // so that we know which individuals have just arrived there
    // (we need to update their resources dependent on their migration
    // group size)
    NSummer_old = NSummer;

    NFlock = 0;

    // actual dispersal
    for (int i = 0; i < NStaging; ++i)
    {
        // later we will consider collective dispersal decisions
        // for now, individuals leave dependent on the current amount of individuals
        // within the staging pool

        pdisperse = 0.5 * (StagingPool[i].phi_a[0] + StagingPool[i].phi_a[1])
            + 0.5 * (StagingPool[i].phi_b[0] + StagingPool[i].phi_b[1]) * NStaging;

        // yes individual goes
        if (gsl_rng_uniform(rng_r) < pdisperse)
        {
            SummerPool[NSummer] = StagingPool[i];
            ++NSummer;
            
            assert(NSummer <= N);


            // delete this individual from the staging population
            StagingPool[i] = StagingPool[NStaging - 1];

            // decrement the number of individuals in the staging population
            --NStaging;
            --i;

            assert(NStaging <= N);
            assert(NStaging >= 0);

            // increase flock size
            ++NFlock;
            
            assert(NFlock <= N);
        }
    }

    // update resource levels for all new individuals that have just
    // been added to the summer pool dependent on their flock size
    // TODO

}

int main()
{
    for (int generation = 0; generation < number_generations; ++generation)
    {
        // time during winter (i.e., days)
        for (int t = 0; t < tmax; ++t)
        {
            winter_dynamics();
        }

        summer_dynamics();

        write_stats();
    }
}



