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
int NWinter = 0;

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
    //
    // collective dispersal elevation
    double phi_a[2];

    // collective dispersal slope
    double phi_b[2];
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

// initialize the population at the start of the simulation
void init_population()
{
    // loop through all individuals in the wintering ground
    // and give them values
    for (int i = 0; i < N; ++i)
    {
        WinterPop[i].resources = 0.0;

        for (int j = 0; j < 2; ++j)
        {
            // initialize allelic values for theta elevation and slope
            WinterPop[i].theta_a[j] = init_theta_a;
            WinterPop[i].theta_b[j] = init_theta_b;
            
            // initialize allelic values for phi elevation and slope
            WinterPop[i].phi_a[j] = init_phi_a;
            WinterPop[i].phi_b[j] = init_phi_b;
        }
    }

    NWinter = N;
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

    // foraging of individuals who are just at the wintering site
    // and who have yet to decide to go to the staging site
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


    // foraging of individuals who are already at the staging site
    for (int i = 0; i < NStaging; ++i)
    { 
        // indivuals who are already at the staging site
        // continue to forage at the staging site
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

    // individuals decide whether to go to staging site
    // i.e., prepare for dispersal
    // signal to disperse
    for (int i = 0; i < NWinter; ++i) 
    {
        // reaction norm dependent on resources
        // resulting in signaling a willingness to disperse
        // => go to the staging level
        psignal = 0.5 * (WinterPop[i].theta_a[0] + WinterPop[i].theta_a[1])
            + 0.5 * (WinterPop[i].theta_b[0] + WinterPop[i].theta_b[1]) * WinterPop[i].resources;

        // does individual want to signal to others to be ready for departure?
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
    for (int i = NSummer_old; i < NSummer; ++i)
    {
        // TODO: think about relationship between flock size and
        // resource reduction
        SummerPool[i].resources -= 1.0 / NFlock;

        // and reduce it by time of arrival
        // TODO think more about this function
        SummerPool[i].resources -= arrival_resource_decay * t;
    }
}

// mutation of a certain allele with value val
// given mutation rate mu and mutational distribution stdev sdmu
double mutation(double val, double mu, double sdmu)
{
    if (gsl_rng_uniform(rng_r) < mu)
    {
        val += gsl_ran_gaussian(rng_r, sdmu);
    }

    return(val);
}

// create a new offspring
void create_offspring(Individual &mother, Individual &offspring)
{
    offspring.resources = 0.0;

    // inherit theta loci

    // each parental allele has probability 0.5 to make it into offspring
    offspring.theta_a[0] = mutation(mother.theta_a[gsl_rng_uniform_int(rng_r,2)], mu_theta, sdmu_theta);

    offspring.theta_a[1] = mutation(father.theta_a[gsl_rng_uniform_int(rng_r,2)], mu_theta, sdmu_theta);

    // put boundaries on the elevation between 0 and 1
    // to help with the interpretation of the evolved values of the slope
    if (offspring.theta_a[0] < 0.0)
    {
        offspring.theta_a[0] = 0.0;
    }
    else if (offspring.theta[0] > 1.0)
    {
        offspring.theta[0] = 1.0;
    }
    
    if (offspring.theta_a[1] < 0.0)
    {
        offspring.theta_a[1] = 0.0;
    }
    else if (offspring.theta[1] > 1.0)
    {
        offspring.theta[1] = 1.0;
    }

    
    offspring.theta_b[0] = mutation(mother.theta_b[gsl_rng_uniform_int(rng_r,2)], mu_theta, sdmu_theta);
    offspring.theta_b[1] = mutation(father.theta_b[gsl_rng_uniform_int(rng_r,2)], mu_theta, sdmu_theta);

    // inherit phi loci
    offspring.phi_a[0] = mutation(mother.phi_a[gsl_rng_uniform_int(rng_r,2)], mu_phi, sdmu_phi);
    offspring.phi_a[1] = mutation(father.phi_a[gsl_rng_uniform_int(rng_r,2)], mu_phi, sdmu_phi);
    
    offspring.phi_b[0] = mutation(mother.phi_b[gsl_rng_uniform_int(rng_r,2)], mu_phi, sdmu_phi);
    offspring.phi_b[1] = mutation(father.phi_b[gsl_rng_uniform_int(rng_r,2)], mu_phi, sdmu_phi);

    // put boundaries on the elevation between 0 and 1
    // to help with the interpretation of the evolved values of the slope
    if (offspring.phi_a[0] < 0.0)
    {
        offspring.phi_a[0] = 0.0;
    }
    else if (offspring.phi[0] > 1.0)
    {
        offspring.phi[0] = 1.0;
    }
    
    if (offspring.phi_a[1] < 0.0)
    {
        offspring.phi_a[1] = 0.0;
    }
    else if (offspring.phi[1] > 1.0)
    {
        offspring.phi[1] = 1.0;
    }
}



// in summery they reproduce dependent on 
// resources and arrival time
void summer_reproduction()
{
    // mix the list of individuals so no statistical associations
    // can build up between individuals in the same flocks
    // with all sorts of hidden consequences to polymorphism
    // later models may relax this and look at mating dynamics.
    // This will affect things.
//    random_shuffle(SummerPool, SummerPool + NSummer);

    Individual mother, father;

    // auxiliary variable specifing the rounded amount of a mother's
    // resources
    int resource_integer;

    /// auxilary variable specifying the id of the randomly sampled
    // fater
    int father_id;

    // see if population is extinct
    if (NSummer == 1)
    {
        // quit if extinct 
        write_parameters();
        exit(1);
    }

    // mating dynamic. Presumes that there an even 
    // number of individuals
    // so we just discard the last individual
    for (int i = 0; i < NSummer; ++i)
    {
        // get the mother
        mother = SummerPool[i];

        // if mom does not meet minimum standards
        // no reproduction through female function
        if (mother.resources < resource_reproduce_threshold)
        {
            continue;
        }

        // now randomly select a father
        do {
            // sample integer uniformly between 0 and NSummer
            // (not including NSummer itself)
            father_id = gsl_rng_uniform_int(rng_r, NSummer);
        }
        while (father_id == i);

        father = SummerPool[father_id];

        // translate maternal resources to numbers of offspring
        //
        // first round to lowest integer
        resource_integer = floor(mother.resources);

        // TODO (slightly digressing): can we come up with an analytical 
        // description of this rounding process of w into integer values?
        if (gsl_rng_uniform(rng_r) < mother.resources - resource_integer)
        {
            // make an additional offspring
            ++resource_integer;
        }
        
        // for each parent create the offspring
        for (int kid_i = 0; kid_i < resource_integer; ++kid_i)
        {
            Individual kid;

            create_offspring(mother, father, kid);

            // add kid to the stack
            Kids[NKids++] = kid;
        }

        Ndead = N - NSummer;

        // recruit new individuals to the summer pool
        for (int i = 0; i < Ndead; ++i)
        {
            // no kids left to recruit
            if (Nkids == 0)
            {
                break;
            }

            random_kid = gsl_rng_uniform_int(rng_r, NKids);

            // add random kid to population
            SummerPool[NSummer++] = Kids[random_kid];

            //  delete random kid as it has been sampled
            Kids[random_kid] = Kids[NKids - 1];
            --NKids;

        }
    }
}

// gaining resources at breeding ground
// & fly back
void summer_dynamics()
{
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
            summer_reproduction();
        }


        write_stats();
    }
}



