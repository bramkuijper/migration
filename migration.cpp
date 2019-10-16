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
int seed = get_nanoseconds();
mt19937 rng_r{static_cast<long unsigned int>(seed)};
uniform_real_distribution<> uniform(0.0,1.0);


// parameters & variables:
// values of the most of these are overridden in the init_arguments()
// function

// number of individuals in population
const int N = 500;

// number of generations
long int number_generations = 50000;

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

// how much resources individuals obtain on a good vs bad patch
double rgood = 0.0;
double rbad = 0.0;

// minimum level of resources necessary to reproduce 
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

// max number of intervals / season (two seasons: summer, winter)
int tmax = 1000;

int skip = 10;

// stats of flock size and staging
double mean_flock_size_winter = 0.0;
double mean_staging_size_winter = 0.0;
double mean_flock_size_summer = 0.0;
double mean_staging_size_summer = 0.0;
double var_flock_size_winter = 0.0;
double var_staging_size_winter = 0.0;
double var_flock_size_summer = 0.0;
double var_staging_size_summer = 0.0;

// keep track of the current number of 
// individuals signaling to disperse
int NStaging = 0;
int NWinter = 0;
int NSpring_migrant = 0;
int NSummer = 0;
int NBreeders = 0;
int NKids = 0;
int NAutumn_migrant = 0;

struct Individual {
    
    double resources;

    // resource reaction norm (determines entry into staging pool)
    //
    // elevation (baseline leaving rate) 
    double theta_a[2];

    // reaction norm, dependency on the amount of resources
    double theta_b[2];

    // collective dispersal reaction norm 
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
    init_theta_a = atof(argv[1]);
    init_theta_b = atof(argv[2]);
    init_phi_a = atof(argv[3]);
    init_phi_b = atof(argv[4]);
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
    min_migration_cost = atof(argv[17]);
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
    DataFile << "generation;time;mean_theta_a;mean_theta_b;mean_phi_a;mean_phi_b;mean_resources;var_theta_a;var_theta_b;var_phi_a;var_phi_b;var_resources;nwinter;nstaging;nspring_migrant;nsummer;nkids;nautumn_migrant;mean_flock_size_summer;mean_flock_size_winter;mean_staging_size_winter;mean_staging_size_summer;var_flock_size_summer;var_flock_size_winter;var_staging_size_winter;var_staging_size_summer;" << endl;
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
    
    for (int i = 0; i < NWinter; ++i)
    {
        mean_theta_a += 0.5 * (WinterPop[i].theta_a[0] + WinterPop[i].theta_a[1]);
        mean_theta_b += 0.5 * (WinterPop[i].theta_b[0] + WinterPop[i].theta_b[1]);
        
        ss_theta_a += 0.5 * (WinterPop[i].theta_a[0] + WinterPop[i].theta_a[1]) * 0.5 * (WinterPop[i].theta_a[0] + WinterPop[i].theta_a[1]);
        ss_theta_b += 0.5 * (WinterPop[i].theta_b[0] + WinterPop[i].theta_b[1]) * 0.5 * (WinterPop[i].theta_b[0] + WinterPop[i].theta_b[1]);

        mean_phi_a += 0.5 * (WinterPop[i].phi_a[0] + WinterPop[i].phi_a[1]);
        mean_phi_b += 0.5 * (WinterPop[i].phi_b[0] + WinterPop[i].phi_b[1]);
        
        ss_phi_a += 0.5 * (WinterPop[i].phi_a[0] + WinterPop[i].phi_a[1]) * 0.5 * (WinterPop[i].phi_a[0] + WinterPop[i].phi_a[1]);
        ss_phi_b += 0.5 * (WinterPop[i].phi_b[0] + WinterPop[i].phi_b[1]) * 0.5 * (WinterPop[i].phi_b[0] + WinterPop[i].phi_b[1]);

        mean_resources += WinterPop[i].resources;
        ss_resources += WinterPop[i].resources * WinterPop[i].resources;
    }

    for (int i = 0; i < NSummer; ++i)
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

    mean_theta_a /= NSummer + NWinter;
    mean_theta_b /= NSummer + NWinter;
    mean_phi_a /= NSummer + NWinter;
    mean_phi_b /= NSummer + NWinter;
    mean_resources /= NSummer + NWinter;
    
    ss_theta_a /= NSummer + NWinter;
    ss_theta_b /= NSummer + NWinter;
    ss_phi_a /= NSummer + NWinter;
    ss_phi_b /= NSummer + NWinter;
    ss_resources /= NSummer + NWinter;

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
        << mean_resources << ";"
        << var_theta_a << ";"
        << var_theta_b << ";"
        << var_phi_a << ";"
        << var_phi_b << ";"
        << var_resources << ";"
        << NWinter << ";"
        << NStaging << ";"
		<< NSpring_migrant << ";"
        << NSummer << ";" 
		<< NBreeders << ";"
        << NKids << ";" 
		<< NAutumn_migrant << ";"
		<< mean_flock_size_winter << ";" 
        << mean_flock_size_summer << ";" 
        << mean_staging_size_winter << ";" 
        << mean_staging_size_summer << ";"
		<< var_flock_size_winter << ";"
		<< var_flock_size_summer << ";"
		<< var_staging_size_winter << ";"
		<< var_staging_size_summer << ";"			
        << endl;

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

// individuals in both summer and winter populations
// die at a certain rate. If this function is called in winter
// the summer pool will be empty so no individuals die there
// If this function is called in summer, however, both the summer
// ground individuals die, as well as the individuals who have
// stayed at the wintering ground
void mortality()
{
    for (int i = 0; i < NWinter;++i)
    {
        // individual dies; replace with end of the stack individual
        if (uniform(rng_r) < pmort)
        {
            WinterPop[i] = WinterPop[NWinter - 1];
            --NWinter;
            --i;
        }
    }

    for (int i = 0; i < NSummer;++i)
    {
        // individual dies; replace with end of the stack individual
        if (uniform(rng_r) < pmort)
        {
            SummerPop[i] = SummerPop[NSummer - 1];
            --NSummer;
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
    for (int i = 0; i < NStaging; ++i)
    {
        WinterPop[NWinter++] = StagingPool[i];

    }

    // just double check that NWinter does not exceed max population size
    assert(NWinter <= N);

    NStaging = 0;  // Is this the problem? NStaging is defined to be 0
}

double get_migration_cost(int const flock_size)
{
    double total_migration_cost = max_migration_cost - migration_cost_decay * 
            pow((double) flock_size / N, migration_cost_power);

    if (total_migration_cost < min_migration_cost)
    {
        total_migration_cost = min_migration_cost;
    }

    assert(total_migration_cost >= 0.0);
    assert(total_migration_cost <= 1.0);

    return(total_migration_cost);
}


// the dynamics of the population at the wintering ground
void winter_dynamics(int t)
{
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
    for (int i = 0; i < NWinter; ++i)
    {
        if (uniform(rng_r) < pgood) // good resource chosen
        {
            WinterPop[i].resources += rgood;
        }
        else
        {
            WinterPop[i].resources += rbad;
        }
    
    } // ok, resource dynamic done


    // foraging of individuals who are already at the staging site
    for (int i = 0; i < NStaging; ++i)
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

    assert(NWinter <= N);
    assert(NWinter >= 0);

    double psignal = 0.0;

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
        if (uniform(rng_r) < psignal)
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
    int NSummer_old = NSummer;

    int NFlock = 0;

    double pdisperse = 0.0;

    int NStaging_start = NStaging;

    // actual dispersal
    for (int i = 0; i < NStaging; ++i)
    {
        // later we will consider collective dispersal decisions
        // for now, individuals leave dependent on the current amount of individuals
        // within the staging pool

        pdisperse = 0.5 * (StagingPool[i].phi_a[0] + StagingPool[i].phi_a[1])
            + 0.5 * (StagingPool[i].phi_b[0] + StagingPool[i].phi_b[1]) * NStaging_start;

        // yes individual goes
        if (uniform(rng_r) < pdisperse)
        {
            SummerPop[NSummer] = StagingPool[i];
            ++NSummer;
            
            assert(NSummer <= N);


            // delete this individual from the staging population
            StagingPool[i] = StagingPool[NStaging - 1];

            // decrement the number of individuals in the staging population
            --NStaging;
            --i;
			
			// but increment the number of individuals recorded as spring migrants
			++NSpring_migrant;

            assert(NStaging <= N);
            assert(NStaging >= 0);

            // increase flock size
            ++NFlock;
            
            assert(NFlock <= N);
        }
    }

    double total_migration_cost;

    // update resource levels for all new individuals that have just
    // been added to the summer pool dependent on their flock size
    for (int i = NSummer_old; i < NSummer; ++i)
    {
		total_migration_cost = get_migration_cost(NFlock);



        // resources are reduced due to migration,
        // yet this depends on group size in a curvilinear fashion
        SummerPop[i].resources = SummerPop[i].resources * total_migration_cost;

        // and reduce it by time of arrival
        // TODO think more about this function
        SummerPop[i].resources = arrival_resource_decay * t;

		// Surviving spring migrants are added to the count of breeders
		if (SummerPop[i].resources >= resource_reproduce_threshold)
		{
			++NBreeders;
		}
        
		// death due to starvation
        if (SummerPop[i].resources < 0)
        {
            SummerPop[i] = SummerPop[NSummer - 1];
            --NSummer;
            --i;
        }
		
		
    }

    // add current dispersal flock size to stats...IF =/= 0
    if (NFlock > 0)
	{
		mean_flock_size_winter += NFlock;
	    mean_staging_size_winter += NStaging_start;
	
		var_flock_size_winter += NFlock * NFlock;
		var_staging_size_winter += NStaging_start * NStaging_start;
	}		
} // end winter_dynamics

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
    if (NSummer == 1)
    {
        // quit if extinct 
        write_parameters(DataFile);
        exit(1);
    }

    // use a flexible array for the kids
    vector<Individual> Kids;

    uniform_int_distribution<> summer_sample(0, NSummer - 1);
    // mating dynamic. Presumes that there an even 
    // number of individuals
    // so we just discard the last individual
    for (int i = 0; i < NSummer; ++i)
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
            // sample integer uniformly between 0 and NSummer
            // (not including NSummer itself)
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

    NKids = Kids.size();

    // number of dead individuals is the max population
    // minus the current individuals in the summer population
    // minus the current individuals who stayed at the wintering ground
    int Ndead = N - NSummer - NWinter;

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
        SummerPop[NSummer++] = Kids[random_kid];

        //  delete random kid as it has been sampled
        Kids[random_kid] = Kids[Kids.size() - 1];

        Kids.pop_back();
    }
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
    for (int i = 0; i < NSummer; ++i)
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
    for (int i = 0; i < NStaging; ++i)
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

    assert(NSummer<= N);
    assert(NSummer>= 0);

    double psignal = 0.0;


    // individuals decide whether to go to staging site
    // i.e., prepare for dispersal
    // signal to disperse
    for (int i = 0; i < NSummer; ++i) 
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
            StagingPool[NStaging] = SummerPop[i];
            ++NStaging; // increment the number of individuals in the staging pool

            assert(NStaging <= N);
            assert(NStaging >= 0);

            // delete this individual from the summer population
            SummerPop[i] = SummerPop[NSummer - 1];

            // decrement the number of individuals in the summer population
            --NSummer;
            --i;
        }
    }

    // store current number of individuals at the breeding ground
    // so that we know which individuals have just arrived there
    // (we need to update their resources dependent on their migration
    // group size)
    int NWinter_old = NWinter;

    int NFlock = 0;

    double pdisperse = 0.0;

    int NStaging_start = NStaging;

    // actual dispersal
    for (int i = 0; i < NStaging; ++i)
    {
        // later we will consider collective dispersal decisions
        // for now, individuals leave dependent on the current amount of individuals
        // within the staging pool

        pdisperse = 0.5 * (StagingPool[i].phi_a[0] + StagingPool[i].phi_a[1])
            + 0.5 * (StagingPool[i].phi_b[0] + StagingPool[i].phi_b[1]) * NStaging_start;

        // yes individual goes
        if (uniform(rng_r) < pdisperse)
        {
            WinterPop[NWinter] = StagingPool[i];
            ++NWinter;
            
            assert(NWinter <= N);


            // delete this individual from the staging population
            StagingPool[i] = StagingPool[NStaging - 1];

            // decrement the number of individuals in the staging population
            --NStaging;
            --i;
			
			// but increment the number of individuals recorded as autumn migrants
			++NSpring_migrant;
			++i;

            assert(NStaging <= N);
            assert(NStaging >= 0);

            // increase flock size
            ++NFlock;
            
            assert(NFlock <= N);
        }
    }

    double total_migration_cost = 0.0;

    // update resource levels for all new individuals that have just
    // been added to the pool dependent on their flock size
    for (int i = NWinter_old; i < NWinter; ++i)
    {
        total_migration_cost = get_migration_cost(NFlock);

        if (total_migration_cost < min_migration_cost)
        {
            total_migration_cost = min_migration_cost;
        }

        assert(total_migration_cost >= 0.0);
        assert(total_migration_cost <= 1.0);

        // resources are reduced due to migration,
        // yet this depends on group size in a curvilinear fashion
        SummerPop[i].resources = SummerPop[i].resources * total_migration_cost;

        // and reduce it by time of arrival
        // TODO think more about this function
        WinterPop[i].resources -= arrival_resource_decay * t;

        // death due to starvation
        if (WinterPop[i].resources < 0)
        {
            WinterPop[i] = WinterPop[NWinter - 1];
            --NWinter;
            --i;
        }

    }
    
    // add current dispersal flock size to stats
    mean_flock_size_summer += NFlock;
    mean_staging_size_summer += NStaging_start;
	
	var_flock_size_summer += NFlock * NFlock;
	var_staging_size_summer += NStaging_start * NStaging_start;
}


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
        mean_flock_size_winter = 0.0;
        mean_staging_size_winter = 0.0;
		var_flock_size_winter = 0.0;
		var_staging_size_winter = 0.0;

        // time during winter (i.e., days)
        // during which individuals forage
        for (int t = 0; t < tmax; ++t)
        {
            winter_dynamics(t);
        }

        // now take averages over all timesteps that individuals can join groups
        mean_flock_size_winter /= tmax;
        mean_staging_size_winter /= tmax;
		
		// now record variance in flock size and staging size over the season
		var_flock_size_winter = (var_flock_size_winter / tmax) - (mean_flock_size_winter * mean_flock_size_winter);
		var_staging_size_winter = (var_staging_size_winter / tmax) - (mean_staging_size_winter * mean_staging_size_winter);	
        
        // all individuals that wanted to migrate have migrated now
        // all remainers are going to stay at wintering ground
        clear_staging_pool();

        // let individuals die with a certain probability 
        mortality();


        // have individuals reproduce after they migrated to the summer spot
        summer_reproduction(DataFile);
        
        // set flock size stats to 0 before summer dynamics starts
        mean_flock_size_summer = 0.0;
        mean_staging_size_summer = 0.0;

        // time during summer (i.e., days)
        // during which individuals forage
        for (int t = 0; t < tmax; ++t)
        {
            summer_dynamics(t);
        }
        
        // now take averages over all timesteps that individuals can join groups
        mean_flock_size_summer /= tmax;
        mean_staging_size_summer /= tmax;
		
		// now record variance in summer flock size and staging size over the season
		var_flock_size_summer = (var_flock_size_summer / tmax) - (mean_flock_size_summer * mean_flock_size_summer);
		var_staging_size_summer = (var_staging_size_summer / tmax) - (mean_staging_size_summer * mean_staging_size_summer);


        // all individuals who remain at the summer ground die
        NSummer = 0;
        NStaging = 0;

        // let individuals die with a certain probability 
        mortality();

        if (generation % skip == 0)
        {
            write_stats(DataFile, generation, 2);
        }
    }

    write_parameters(DataFile);
}
