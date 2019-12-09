// Collective migration when resources vary
// Bram Kuijper & Simon Evans
// 2019
//
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <vector>
#include <cmath>
#include <random>

// various functions, such as unique filename creation
#include "auxiliary.h"

#define DEBUG

// standard namespace
using namespace std;

// set random seed etc
unsigned int seed = get_nanoseconds();
mt19937 rng_r{seed};
uniform_real_distribution<> uniform(0.0,1.0);


// parameters & variables:
// values of the most of these are overridden in the init_arguments()
// function

// number of individuals in population
const int N = 2000; 

// number of generations
long int number_generations = 100000; // simulation time doubled to see if evolving traits will stabilise

// initial values for phi (social dependency) and theta (resource dependency)
// a is an intercept, b is a gradient
double init_theta_a = 0.0; 
double init_theta_b = 0.0;
double init_phi_a = 0.0;
double init_phi_b = 0.0;

// mortality probability 
double pmort = 0.0;

// initial probability per season to encounter a good resource patch
double pgood_init = 0.0;

// the time point at which 
// the probability of encountering a good environment becomes 0
int t_good_ends = 0.0;

// how much resource individuals obtain on a good vs bad patch
double rgood = 0.0;
double rbad = 0.0;

// minimum resource level necessary to reproduce 
double resource_reproduce_threshold = 0.0; 

// how quickly resources decay per day arriving later than 0
double arrival_resource_decay = 0.0;

// mutation rates
double mu_theta = 0.0;
double mu_phi = 0.0;
double sdmu_theta = 0.0;
double sdmu_phi = 0.0;

// migration cost function
double max_migration_cost = 0.0;
double min_migration_cost = 0.0;
double migration_cost_decay = 0.0;
double migration_cost_power = 0.0;

// max number of intervals per season (two seasons: summer, winter)
int tmax = 5000;

int skip = 10;

// stats of flock size and staging
double mean_spring_flock_size = 0.0;
double mean_spring_staging_size = 0.0;
double mean_autumn_flock_size = 0.0;
double mean_autumn_staging_size = 0.0;
double var_spring_flock_size = 0.0;
double var_spring_staging_size = 0.0;
double var_autumn_flock_size = 0.0;
double var_autumn_staging_size = 0.0;
double ss_spring_flock_size = 0.0;
double ss_autumn_flock_size = 0.0;

// keep track of the current number of 
// individuals in various seasons/demographics
int staging_pop = 0;
int winter_pop = 0;
int spring_migrant_pop = 0;
int summer_pop = 0;
int breeder_pop = 0;
int offspring_pop = 0;
int autumn_migrant_pop = 0;
int n_spring_flocks = 0;  // recording the number of spring flocks (tmax - n(unusued departure intervals))
int n_autumn_flocks = 0;

double ss_spring_migrant_pop = 0.0;
double ss_autumn_migrant_pop = 0.0;
double ss_spring_staging_size = 0.0;
double ss_autumn_staging_size = 0.0;

struct Individual {
    
    double resources;

    // RESOURCE SENSITIVITY reaction norm (determines entry into staging pool)
    //
    // elevation (baseline leaving rate) 
    double theta_a[2];

    // reaction norm, dependency on the amount of resources
    double theta_b[2];

    // COLLECTIVE DISPERSAL reaction norm 
    // determines migration dependent on number of individuals
    //
    // collective dispersal elevation
    double phi_a[2];

    // collective dispersal slope, dependency on number of individuals
    double phi_b[2];
};

// wintering ground
Individual WinterPop[N];
Individual StagingPool[N];
Individual SummerPop[N];

// get parameters from the command line when 
// running the executable file
void init_arguments(int argc, char **argv)
{
    init_phi_a = atof(argv[1]);
    init_phi_b = atof(argv[2]);
    init_theta_a = atof(argv[3]);
    init_theta_b = atof(argv[4]);
    pmort = atof(argv[5]);
    pgood_init = atof(argv[6]);
    t_good_ends = atoi(argv[7]);
    rgood = atof(argv[8]);
    rbad = atof(argv[9]);
    arrival_resource_decay = atof(argv[10]);
    resource_reproduce_threshold = atof(argv[11]);
    mu_theta = atof(argv[12]);
    mu_phi = atof(argv[13]);
    sdmu_theta = atof(argv[14]);
    sdmu_phi = atof(argv[15]);
    max_migration_cost = atof(argv[16]);
    min_migration_cost = atof(argv[16])/2;
    migration_cost_decay = atof(argv[18]);
    migration_cost_power = atof(argv[19]);
    tmax = atoi(argv[20]);
}

// write down all parameters in the file
void write_parameters(ofstream &DataFile)
{
    DataFile << endl << endl
            << "init_theta_a;" << init_theta_a << endl
            << "init_theta_b;" << init_theta_b << endl
            << "init_phi_a;" << init_phi_a << endl
            << "init_phi_b;" << init_phi_b << endl
            << "pmort;" << pmort << endl
            << "pgood_init;" << pgood_init << endl
            << "t_good_ends;" << t_good_ends << endl
            << "rgood;" << rgood << endl
            << "rbad;" << rbad << endl
            << "arrival_resource_decay;" << arrival_resource_decay << endl
            << "resource_reproduce_threshold;" << resource_reproduce_threshold << endl
            << "mu_theta;" << mu_theta << endl
            << "mu_phi;" << mu_phi << endl
            << "sdmu_theta;" << sdmu_theta << endl
            << "sdmu_phi;" << sdmu_phi << endl
            << "tmax;" << tmax << endl
            << "N;" << N << endl
            << "max_migration_cost;" << max_migration_cost << endl
            << "min_migration_cost;" << min_migration_cost << endl
            << "migration_cost_decay;" << migration_cost_decay << endl
            << "migration_cost_power;" << migration_cost_power << endl
            << "seed;" << seed << endl;
}

// list of the data headers at the start of the file
void write_data_headers(ofstream &DataFile)
{
    DataFile << "generation;time_interval;mean_theta_a;mean_theta_b;mean_phi_a;mean_phi_b;var_theta_a;var_theta_b;var_phi_a;var_phi_b;mean_resources;var_resources;winter_pop;mean_spring_staging_size;var_spring_staging_size;spring_migrant_pop;n_spring_flocks;mean_spring_flock_size;var_spring_flock_size;breeder_pop;offspring_pop;mean_autumn_staging_size;var_autumn_staging_size;autumn_migrant_pop;n_autumn_flocks;mean_autumn_flock_size;var_autumn_flock_size;" << endl;
}

// write data both for winter and summer populations
void write_stats(ofstream &DataFile, int generation, int timestep)
{
    double mean_theta_a = 0.0;
    double ss_theta_a = 0.0;
    double mean_theta_b = 0.0;
    double ss_theta_b = 0.0;

    double mean_phi_a = 0.0;
    double ss_phi_a = 0.0;
    double mean_phi_b = 0.0;
    double ss_phi_b = 0.0;

    double mean_resources = 0.0;
    double ss_resources = 0.0;
    
    for (int i = 0; i < winter_pop; ++i)  // So here we are cycling one by one through the winter population
    {
        // each character (elevation and slope of the 
		// two reaction norms) is genetically controlled
		// by a single gene (diploid) exhibiting incomplete dominance
		// (hence *0.5)
		mean_theta_a += 0.5 * (WinterPop[i].theta_a[0] + WinterPop[i].theta_a[1]);
        mean_theta_b += 0.5 * (WinterPop[i].theta_b[0] + WinterPop[i].theta_b[1]);
        
        ss_theta_a += 0.5 * (WinterPop[i].theta_a[0] + WinterPop[i].theta_a[1]) * 0.5 * (WinterPop[i].theta_a[0] + WinterPop[i].theta_a[1]);
        ss_theta_b += 0.5 * (WinterPop[i].theta_b[0] + WinterPop[i].theta_b[1]) * 0.5 * (WinterPop[i].theta_b[0] + WinterPop[i].theta_b[1]);

        mean_phi_a += 0.5 * (WinterPop[i].phi_a[0] + WinterPop[i].phi_a[1]);
        mean_phi_b += 0.5 * (WinterPop[i].phi_b[0] + WinterPop[i].phi_b[1]);
        
        ss_phi_a += 0.5 * (WinterPop[i].phi_a[0] + WinterPop[i].phi_a[1]) * 0.5 * (WinterPop[i].phi_a[0] + WinterPop[i].phi_a[1]);
        ss_phi_b += 0.5 * (WinterPop[i].phi_b[0] + WinterPop[i].phi_b[1]) * 0.5 * (WinterPop[i].phi_b[0] + WinterPop[i].phi_b[1]);

        mean_resources += WinterPop[i].resources;  // the resource level of individual i 
        ss_resources += WinterPop[i].resources * WinterPop[i].resources;
    }

    for (int i = 0; i < summer_pop; ++i)  // for each individual in the summer population:
    {
        mean_theta_a += 0.5 * (SummerPop[i].theta_a[0] + SummerPop[i].theta_a[1]);
        mean_theta_b += 0.5 * (SummerPop[i].theta_b[0] + SummerPop[i].theta_b[1]);
        
        ss_theta_a += 0.5 * (SummerPop[i].theta_a[0] + SummerPop[i].theta_a[1]) * 0.5 * (SummerPop[i].theta_a[0] + SummerPop[i].theta_a[1]);
        ss_theta_b += 0.5 * (SummerPop[i].theta_b[0] + SummerPop[i].theta_b[1]) * 0.5 * (SummerPop[i].theta_b[0] + SummerPop[i].theta_b[1]);

        mean_phi_a += 0.5 * (SummerPop[i].phi_a[0] + SummerPop[i].phi_a[1]);
        mean_phi_b += 0.5 * (SummerPop[i].phi_b[0] + SummerPop[i].phi_b[1]);
        
        ss_phi_a += 0.5 * (SummerPop[i].phi_a[0] + SummerPop[i].phi_a[1]) * 0.5 * (SummerPop[i].phi_a[0] + SummerPop[i].phi_a[1]);
        ss_phi_b += 0.5 * (SummerPop[i].phi_b[0] + SummerPop[i].phi_b[1]) * 0.5 * (SummerPop[i].phi_b[0] + SummerPop[i].phi_b[1]);

        mean_resources += SummerPop[i].resources;
        ss_resources += SummerPop[i].resources * SummerPop[i].resources;
    }

    // Here we are gathering the average over the two populations (i.e., summer AND winter) WHY DO WE DO THIS? Maybe to include individuals that are spending the summer at the overwintering grounds?
	mean_theta_a /= summer_pop + winter_pop;
    mean_theta_b /= summer_pop + winter_pop;
    mean_phi_a /= summer_pop + winter_pop;
    mean_phi_b /= summer_pop + winter_pop;
	mean_resources /= summer_pop + winter_pop;

    
    ss_theta_a /= summer_pop + winter_pop;
    ss_theta_b /= summer_pop + winter_pop;
    ss_phi_a /= summer_pop + winter_pop;
    ss_phi_b /= summer_pop + winter_pop;
    ss_resources /= summer_pop + winter_pop;

    double var_theta_a = ss_theta_a - mean_theta_a * mean_theta_a;
    double var_theta_b = ss_theta_b - mean_theta_b * mean_theta_b;
    double var_phi_a = ss_phi_a - mean_phi_a * mean_phi_a;
    double var_phi_b = ss_phi_b - mean_phi_b * mean_phi_b;
	double var_resources = ss_resources - mean_resources * mean_resources;

    DataFile 
        << generation << ";"
        << timestep << ";"
        << mean_theta_a << ";"
        << mean_theta_b << ";"
        << mean_phi_a << ";"
        << mean_phi_b << ";"
        << var_theta_a << ";"
        << var_theta_b << ";"
        << var_phi_a << ";"
        << var_phi_b << ";"
	    << mean_resources << ";"
        << var_resources << ";"
        << winter_pop << ";"
	    << mean_spring_staging_size << ";" 
		<< var_spring_staging_size << ";"
		<< spring_migrant_pop << ";"
		<< n_spring_flocks << ";"
		<< mean_spring_flock_size << ";" 
		<< var_spring_flock_size << ";"
		<< breeder_pop << ";"
        << offspring_pop << ";" 
		<< mean_autumn_staging_size << ";"
		<< var_autumn_staging_size << ";"			
        << autumn_migrant_pop << ";"
		<< n_autumn_flocks << ";"
		<< mean_autumn_flock_size << ";" 
		<< var_autumn_flock_size << ";"
		<< endl;
// ENDS: write data both for winter and summer populations
}


// initialize the population at the start of the simulation
void init_population()
{
    // loop through all individuals in the wintering ground
    // and give them values
    for (int i = 0; i < N; ++i)
    {
        WinterPop[i].resources = 0.0;  // at start of simulation, initial resource level for all individuals is 0

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

    winter_pop = N;
} // ENDS: initialize the population at the start of the simulation

// individuals in both summer and winter populations
// die at a certain rate. If this function is called in winter
// the summer pool will be empty so no individuals die there
// If this function is called in summer, however, both the summer
// ground individuals die, as well as the individuals who have
// stayed at the wintering ground
void mortality()
{
    for (int i = 0; i < winter_pop;++i)
    {
        // individual dies; replace with end of the stack individual
        if (uniform(rng_r) < pmort)
        {
            WinterPop[i] = WinterPop[winter_pop - 1];
            --winter_pop;
            --i;
        }
    }

    for (int i = 0; i < summer_pop;++i)
    {
        // individual dies; replace with end of the stack individual
        if (uniform(rng_r) < pmort)
        {
            SummerPop[i] = SummerPop[summer_pop - 1];
            --summer_pop;
            --i;
        }
    }
}

// remove individuals from the staging pool and put them
// back in the original population
void clear_staging_pool()
{
	
	// put individuals from staging pool (which haven't migrated) 
    // back in the original population
    for (int i = 0; i < staging_pop; ++i)
    {
        WinterPop[winter_pop++] = StagingPool[i];

    }

    // just double check that winter_pop does not exceed max population size
	assert(winter_pop <= N);

    staging_pop = 0;
}  // ENDS STAGING POOL CLEARANCE

// MIGRATION COST
double get_migration_cost(int const flock_size)
{
    double total_migration_cost = max_migration_cost - migration_cost_decay * 
            pow((double) flock_size / N, migration_cost_power);

    if (total_migration_cost < min_migration_cost)
    {
        total_migration_cost = min_migration_cost;
    }

    assert(total_migration_cost >= 0.0);
    assert(total_migration_cost <= max_migration_cost);

    return(total_migration_cost);
}

// the dynamics of the population at the wintering ground
void winter_dynamics(int t)
{
	// Ensure that at least two individuals are still alive (to allow sexual reproduction)
	
	// individuals forage
    // individuals accumulate resources
    // individuals make dispersal decisions

    // probability of encountering a good resource
    // decays with time until hits time t = t_good_ends 
    // when pgood = 0
    double pgood = pgood_init - pgood_init / t_good_ends * t;

    // set lower boundary to the probability
    if (pgood <= 0)
    {
        pgood = 0;
    }

    // foraging of individuals who are just at the wintering site
    // and who have yet to decide to go to the staging site
    for (int i = 0; i < winter_pop; ++i)
    {
        if (uniform(rng_r) < pgood) // good resource chosen
        {
            WinterPop[i].resources += rgood;
        }
        else
        {
            WinterPop[i].resources += rbad;
        } // ENDS winter foraging loop
    
    } // ok, resource dynamic done


    // foraging of individuals who are already at the staging site
    for (int i = 0; i < staging_pop; ++i)
    { 
        // indivuals who are already at the staging site
        // continue to forage at the staging site
        if (uniform(rng_r) < pgood) // good resource chosen
        {
            StagingPool[i].resources += rgood;
        }
        else
        {
            StagingPool[i].resources += rbad;
        }
    } // ENDS staging site foraging loop

    assert(winter_pop <= N);
    assert(winter_pop >= 0);  

    double psignal = 0.0;

    // individuals decide whether to go to staging site
    // i.e., prepare for dispersal
    // signal to disperse
    for (int i = 0; i < winter_pop; ++i) 
    {
        // reaction norm dependent on resources
        // resulting in signaling a willingness to disperse
        // => go to the staging level
        psignal = 0.5 * (WinterPop[i].theta_a[0] + WinterPop[i].theta_a[1])
            + 0.5 * (WinterPop[i].theta_b[0] + WinterPop[i].theta_b[1]) * WinterPop[i].resources;

        // does individual want to signal to others to be ready for departure?
        if (uniform(rng_r) < psignal)
        {
            // add individual to the staging pool
            StagingPool[staging_pop] = WinterPop[i];
            ++staging_pop; // increment the number of individuals in the staging pool

            assert(staging_pop <= N);
            assert(staging_pop >= 0);

            // delete this individual from the winter population
            WinterPop[i] = WinterPop[winter_pop - 1];

            // decrement the number of individuals in the winter population
            --winter_pop;
            --i;
        }
    }

    // store current number of individuals at the breeding ground
    // so that we know which individuals have just arrived there
    // (we need to update their resources dependent on their migration
    // group size)
    int summer_pop_old = summer_pop;

    int NFlock = 0;

    double pdisperse = 0.0;

    int staging_pop_start = staging_pop;

    // actual spring dispersal
    for (int i = 0; i < staging_pop; ++i)
    {
        // later we will consider collective dispersal decisions
        // for now, individuals leave dependent on the current amount of individuals
        // within the staging pool

        pdisperse = 0.5 * (StagingPool[i].phi_a[0] + StagingPool[i].phi_a[1])
            + 0.5 * (StagingPool[i].phi_b[0] + StagingPool[i].phi_b[1]) * staging_pop_start;

        // yes individual goes
        if (uniform(rng_r) < pdisperse)
        {
            SummerPop[summer_pop] = StagingPool[i];
            ++summer_pop;
            
            assert(summer_pop <= N);

			// increment the number of individuals recorded as spring migrants
			//spring_migrant_pop = spring_migrant_pop +1;  // Same result by sum below
			
            // delete this individual from the staging population
            StagingPool[i] = StagingPool[staging_pop - 1];

            // decrement the number of individuals in the staging population
            --staging_pop;
            --i;

            assert(staging_pop <= N);
            assert(staging_pop >= 0);

            // increase flock size
            ++NFlock;
            
            assert(NFlock <= N);
        

		} // ENDS: yes individual goes
		
    } // ENDS ACTUAL SPRING DISPERSAL
	
	double total_migration_cost;

	mean_spring_staging_size += staging_pop_start;
	ss_spring_staging_size += staging_pop_start * staging_pop_start;
	
    if (NFlock > 0){
		n_spring_flocks = n_spring_flocks+1;
		mean_spring_flock_size += NFlock;
		ss_spring_flock_size += NFlock * NFlock;  // Also serves as sum of squares of spring migrant population size
		

	}
	
	// update resource levels for all new individuals that have just
    // been added to the summer pool dependent on their flock size
    for (int i = summer_pop_old; i < summer_pop; ++i)
    {
		total_migration_cost = get_migration_cost(NFlock);

        // resources are reduced due to migration,
        // yet this depends on group size in a curvilinear fashion
        SummerPop[i].resources = SummerPop[i].resources - total_migration_cost;

        // and reduce it by time of arrival
        // TODO think more about this function
        SummerPop[i].resources = arrival_resource_decay * t;

		// Surviving spring migrants are added to the count of breeders
		if (SummerPop[i].resources >= resource_reproduce_threshold)
		{
			breeder_pop = breeder_pop+1;
		}
        
		// death due to starvation
        if (SummerPop[i].resources < 0)
        {
            SummerPop[i] = SummerPop[summer_pop - 1];
            --summer_pop;
            --i;
        } // ends: death due to starvation
		
    } // ENDS: updating resources of migrants
	
} // ENDS WINTER DYNAMICS (looping through t)

// mutation of a certain allele with value val
// given mutation rate mu and mutational distribution stdev sdmu
double mutation(double val, double mu, double sdmu)
{
    if (uniform(rng_r) < mu)
    {
        normal_distribution<> allelic_dist(0,sdmu);
        val += allelic_dist(rng_r);
    }

    return(val);
// ENDS ALLELE MUTATION
}

// create a new offspring
void create_offspring(Individual &mother, Individual &father, Individual &offspring)
{
    bernoulli_distribution allele_sample(0.5);

    offspring.resources = 0.0;

    // inherit theta loci

    // each parental allele has probability 0.5 to make it into offspring
    offspring.theta_a[0] = mutation(mother.theta_a[allele_sample(rng_r)], mu_theta, sdmu_theta);
    offspring.theta_a[1] = mutation(father.theta_a[allele_sample(rng_r)], mu_theta, sdmu_theta);

    offspring.theta_b[0] = mutation(mother.theta_b[allele_sample(rng_r)], mu_theta, sdmu_theta);
    offspring.theta_b[1] = mutation(father.theta_b[allele_sample(rng_r)], mu_theta, sdmu_theta);

    // inherit phi loci
    offspring.phi_a[0] = mutation(mother.phi_a[allele_sample(rng_r)], mu_phi, sdmu_phi);
    offspring.phi_a[1] = mutation(father.phi_a[allele_sample(rng_r)], mu_phi, sdmu_phi);
    
    offspring.phi_b[0] = mutation(mother.phi_b[allele_sample(rng_r)], mu_phi, sdmu_phi);
    offspring.phi_b[1] = mutation(father.phi_b[allele_sample(rng_r)], mu_phi, sdmu_phi);

    for (int allele_i = 0; allele_i < 2; ++allele_i)
    {
        // put boundaries on the elevation between 0 and 1
        // to help with the interpretation of the evolved values of the slope
        if (offspring.theta_a[allele_i] < 0.0)
        {
            offspring.theta_a[allele_i] = 0.0;
        }
        else if (offspring.theta_a[allele_i] > 1.0)
        {
            offspring.theta_a[allele_i] = 1.0;
        }
        
        if (offspring.phi_a[allele_i] < 0.0)
        {
            offspring.phi_a[allele_i] = 0.0;
        }
        else if (offspring.phi_a[allele_i] > 1.0)
        {
            offspring.phi_a[allele_i] = 1.0;
        }
    }   
// ENDS OFFSPRING PRODUCTION
}



// in summary, they reproduce dependent on 
// resources and arrival time
void summer_reproduction(ofstream &DataFile)
{
    // auxiliary variables storing current mom and dad
    Individual mother, father;

    // auxiliary variable specifing the rounded amount of a mother's
    // resources
    int resource_integer;

    /// auxilary variable specifying the id of the randomly sampled
    // father
    int father_id;

    // see if population is extinct
    if (summer_pop == 1)
    {
        // quit if extinct 
        write_parameters(DataFile);
        exit(1);
    }

    // use a flexible array for the kids
    vector<Individual> Kids;

    uniform_int_distribution<> summer_sample(0, summer_pop - 1);
    // mating dynamic. Presumes that there an even 
    // number of individuals
    // so we just discard the last individual
    for (int i = 0; i < summer_pop; ++i)
    {
        // get the mother
        mother = SummerPop[i];

        // if mom does not meet minimum standards
        // no reproduction through female function
        if (mother.resources < resource_reproduce_threshold)
        {
            continue;
        }

        // now randomly select a father
        do {
            // sample integer uniformly between 0 and summer_pop
            // (not including summer_pop itself)
            father_id = summer_sample(rng_r);
        }
        while (father_id == i);

        father = SummerPop[father_id];

        // translate maternal resources to numbers of offspring
        //
        // first round to lowest integer
        resource_integer = floor(mother.resources);

        // TODO (slightly digressing): can we come up with an analytical 
        // description of this rounding process of w into integer values?
        if (uniform(rng_r) < mother.resources - resource_integer)
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
            Kids.push_back(kid);
        }
    }

    offspring_pop = Kids.size();

    // number of dead individuals is the max population
    // minus the current individuals in the summer population
    // minus the current individuals who stayed at the wintering ground
    int Ndead = N - summer_pop - winter_pop;

    int random_kid = 0;

    // recruit new individuals to the summer pool
    for (int i = 0; i < Ndead; ++i)
    {
        // no kids left to recruit
        if (Kids.size() == 0)
        {
            break;
        }
        
        uniform_int_distribution<> kids_sample(0, Kids.size() - 1);
        
        random_kid = kids_sample(rng_r);

        // add random kid to population
        SummerPop[summer_pop++] = Kids[random_kid];

        //  delete random kid as it has been sampled
        Kids[random_kid] = Kids[Kids.size() - 1];

        Kids.pop_back();
    }
// ENDS SUMMER REPRODUCTION
}

// gaining resources at breeding ground
// & fly back
void summer_dynamics(int t)
{
    // probability of encountering a good resource
    double pgood = pgood_init - pgood_init / t_good_ends * t;

    // set lower boundary to the probability
    if (pgood <= 0)
    {
        pgood = 0;
    }

    // foraging of individuals who are just at the breeding (Bram: am I correct that 
	// this should be breeding rather than winter?) site
    // and who have yet to decide to go to the staging site
    for (int i = 0; i < summer_pop; ++i)
    {
        if (uniform(rng_r) < pgood) // good resource chosen
        {
            SummerPop[i].resources += rgood;
        }
        else
        {
            SummerPop[i].resources += rbad;
        }
    
    } // ok resource dynamic done


    // foraging of individuals who are already at the staging site
    for (int i = 0; i < staging_pop; ++i)
    { 
        // indivuals who are already at the staging site
        // continue to forage at the staging site
        if (uniform(rng_r) < pgood) // good resource chosen
        {
            StagingPool[i].resources += rgood;
        }
        else
        {
            StagingPool[i].resources += rbad;
        }
    }

    assert(summer_pop<= N);
    assert(summer_pop>= 0);

    double psignal = 0.0;


    // individuals decide whether to go to staging site
    // i.e., prepare for dispersal
    // signal to disperse
    for (int i = 0; i < summer_pop; ++i) 
    {
        // reaction norm dependent on resources
        // resulting in signaling a willingness to disperse
        // => go to the staging level
        psignal = 0.5 * (SummerPop[i].theta_a[0] + SummerPop[i].theta_a[1])
            + 0.5 * (SummerPop[i].theta_b[0] + SummerPop[i].theta_b[1]) * SummerPop[i].resources;

        // does individual want to signal to others to be ready for departure?
        if (uniform(rng_r) < psignal)
        {
            // add individual to the staging pool
            StagingPool[staging_pop] = SummerPop[i];
            ++staging_pop; // increment the number of individuals in the staging pool

            assert(staging_pop <= N);
            assert(staging_pop >= 0);

            // delete this individual from the summer population
            SummerPop[i] = SummerPop[summer_pop - 1];

            // decrement the number of individuals in the summer population
            --summer_pop;
            --i;
        }
    }

    // store current number of individuals at the breeding ground
    // so that we know which individuals have just arrived there
    // (we need to update their resources dependent on their migration
    // group size)
    int winter_pop_old = winter_pop;

    int NFlock = 0;

    double pdisperse = 0.0;

    int staging_pop_start = staging_pop;

    // actual autumn dispersal
    for (int i = 0; i < staging_pop; ++i)
    {
        // later we will consider collective dispersal decisions
        // for now, individuals leave dependent on the current amount of individuals
        // within the staging pool

        pdisperse = 0.5 * (StagingPool[i].phi_a[0] + StagingPool[i].phi_a[1])
            + 0.5 * (StagingPool[i].phi_b[0] + StagingPool[i].phi_b[1]) * staging_pop_start;

        // yes individual goes
        if (uniform(rng_r) < pdisperse)
        {
            WinterPop[winter_pop] = StagingPool[i];
            ++winter_pop;
            
            assert(winter_pop <= N);


            // delete this individual from the staging population
            StagingPool[i] = StagingPool[staging_pop - 1];

            // decrement the number of individuals in the staging population
            --staging_pop;
            --i;
			
			// but increment the number of individuals recorded as autumn migrants
			++autumn_migrant_pop;
			//++i;  // deleted this because I think it's a mistake (30 Oct 2019)

            assert(staging_pop <= N);
            assert(staging_pop >= 0);

            // increase flock size
            ++NFlock;
            
            assert(NFlock <= N);
        } // Ends: individual goes
	    
    } // ENDS: Autumn dispersal
	
    double total_migration_cost = 0.0;
	
	mean_autumn_flock_size += NFlock;
	ss_autumn_flock_size += NFlock * NFlock;
	mean_autumn_staging_size += staging_pop_start;
	ss_autumn_staging_size += staging_pop_start * staging_pop_start;
	
    if (winter_pop_old < winter_pop){
		++ n_autumn_flocks;
	}

    // update resource levels for all new individuals that have just
    // been added to the pool dependent on their flock size
    for (int i = winter_pop_old; i < winter_pop; ++i)
    {
        total_migration_cost = get_migration_cost(NFlock);

        if (total_migration_cost < min_migration_cost)
        {
            total_migration_cost = min_migration_cost;
        }

        assert(total_migration_cost >= 0.0);
        assert(total_migration_cost <= max_migration_cost);

        // resources are reduced due to migration,
        // yet this depends on group size in a curvilinear fashion
        SummerPop[i].resources = SummerPop[i].resources - total_migration_cost;

        // and reduce it by time of arrival
        // TODO think more about this function
        WinterPop[i].resources -= arrival_resource_decay * t;

        // death due to starvation
        if (WinterPop[i].resources < 0)
        {
            WinterPop[i] = WinterPop[winter_pop - 1];
            --winter_pop;
            --i;
        }

    } // Ends: update resource levels of winter arrivals

} // ENDS: SUMMER DYNAMICS 


// the key part of the code
// accepting command line arguments
int main(int argc, char **argv)
{
    string filename = "sim_migration";
    create_filename(filename);
    ofstream DataFile(filename.c_str());  // output file 

    init_arguments(argc, argv);

    write_data_headers(DataFile);

    init_population();

    for (int generation = 0; generation < number_generations; ++generation)
    {
        mean_spring_flock_size = 0.0;
		mean_spring_staging_size = 0.0;
		var_spring_flock_size = 0.0;
		var_spring_staging_size = 0.0;
		ss_spring_flock_size = 0.0;
		n_spring_flocks = 0;
		spring_migrant_pop = 0.0;
		ss_spring_staging_size = 0.0;

        staging_pop = 0.0;  // Set staging population count to zero before winter dynamics
		
		// time during winter (i.e., days)
        // during which individuals forage
        for (int t = 0; t < tmax; ++t)
        {
            winter_dynamics(t);
			
        }
		
		spring_migrant_pop = mean_spring_flock_size;
		
        // now take averages over all timesteps that individuals did (can) join groups
        mean_spring_flock_size /= n_spring_flocks;
		mean_spring_staging_size /= tmax;
		
		// now record variance in flock size and staging size over the season
		var_spring_flock_size = (ss_spring_flock_size / n_spring_flocks) - (mean_spring_flock_size * mean_spring_flock_size);
		var_spring_staging_size = (ss_spring_staging_size / tmax) - (mean_spring_staging_size * mean_spring_staging_size);	
        
        // all individuals that wanted to migrate have migrated now
        // all remainers are going to stay at wintering ground
        clear_staging_pool();

        // let individuals die with a certain probability 
        mortality();
		
		// summer_pop = 0.0;  // So summer population is counted after mortality event (because population size is otherwise the same as spring migrant population size)

        // have individuals reproduce after they migrated to the summer spot
        summer_reproduction(DataFile);
        
        // set flock size stats to 0 before summer dynamics starts
        mean_autumn_flock_size = 0.0;
        mean_autumn_staging_size = 0.0;
		var_autumn_flock_size = 0.0;
		var_autumn_staging_size = 0.0;
		ss_autumn_flock_size = 0.0;
		n_autumn_flocks = 0.0;
		autumn_migrant_pop = 0.0;
		ss_autumn_migrant_pop = 0.0;
		ss_autumn_staging_size = 0.0;

        // time during summer (i.e., days)
        // during which individuals forage
        for (int t = 0; t < tmax; ++t)
        {
            summer_dynamics(t);
			
        }
        
        // now take averages over all timesteps that individuals did (can) join groups
        mean_autumn_flock_size /= n_autumn_flocks;
        mean_autumn_staging_size /= tmax;
		
		// now record variance in autumn flock size and staging size over the season
		var_autumn_flock_size = (ss_autumn_flock_size / n_autumn_flocks) - (mean_autumn_flock_size * mean_autumn_flock_size);
		var_autumn_staging_size = (ss_autumn_staging_size / tmax) - (mean_autumn_staging_size * mean_autumn_staging_size);
        
		if (generation % skip == 0)
        {
            write_stats(DataFile, generation, 1000);  // Surely we want to have let the seasons play out? (so timestep =/= 2, which was the previous setting [07 November 2019])
        }

        // all individuals who remain at the summer grounds die
        summer_pop = 0;
        staging_pop = 0;
		breeder_pop = 0;

        // let individuals die with a certain probability 
        mortality();
				
    } // ENDS: GENERATION

    write_parameters(DataFile);
}
