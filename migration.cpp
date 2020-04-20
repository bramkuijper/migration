// Collective migration when resources vary
// Bram Kuijper & Simon Evans
// 2019
//
#define DEBUG

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
const int N = 1500;

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

// how much resource individuals obtain on a good vs bad patch
double rgood_init = 0.0;
double rbad_init = 0.0;
double rgood = 0.0;
double rbad = 0.0;

double resource_reproduction_threshold = 0.0;  // minimum resource level necessary to reproduce 
double resource_starvation_threshold = 0.0;  // minimum resource level necessary to survive
double resource_max = 0.0;  // maximum resource value an individual can achieve

// how quickly resources decay per day arriving later than 0
// TODO: do we really need this?
double arrival_resource_decay = 0.0;

// mutation rates
double mu_theta = 0.0;
double mu_phi = 0.0;
double sdmu_theta = 0.0;
double sdmu_phi = 0.0;

// migration cost function
double migration_cost_decay = 0.0;
double migration_cost_power = 0.0;
double max_migration_cost = 0.0;

// max number of intervals per season (two seasons: summer, winter)
int tmax = 5000;
int twinter = 0.0;

int skip = 100;

// stats of flock size and staging
double mean_spring_flock_size = 0.0;
double mean_spring_staging_size = 0.0;
double mean_autumn_flock_size = 0.0;
double mean_autumn_staging_size = 0.0;
double var_spring_flock_size = 0.0;
double var_spring_staging_size = 0.0;
double mean_spring_cost = 0.0;
double var_spring_cost = 0.0;
double var_autumn_flock_size = 0.0;
double var_autumn_staging_size = 0.0;
double ss_spring_flock_size = 0.0;
double ss_autumn_flock_size = 0.0;
double mean_autumn_cost = 0.0;
double var_autumn_cost = 0.0;
double mean_latency = 0.0;
double var_latency = 0.0;
double ss_latency = 0.0;
double mean_cost = 0.0;
double var_cost = 0.0;
double ss_cost = 0.0;
double mean_departure = 0.0;
double ss_departure = 0.0;
double mean_signal_timing = 0.0;
double var_signal_timing = 0.0;
double ss_signal_timing = 0.0;
double mean_age = 0.0;
double var_age = 0.0;
double ss_age = 0.0;

// keep track of the current number of 
// individuals in various seasons/demographics
int staging_pop = 0;
int winter_pop = 0;
int spring_pop_start = 0;
int spring_nonmigrant_pop = 0;
int spring_migrant_pop = 0;
int summer_pop = 0;
int autumn_pop_start = 0;
int breeder_pop = 0;
int offspring_pop = 0;
int autumn_nonmigrant_pop = 0;
int autumn_migrant_pop = 0;
int n_spring_flocks = 0;  // recording the number of spring flocks (tmax - n(unusued departure intervals))
int n_autumn_flocks = 0;
int summer_pop_old = 0;  // 06/02/20: So that I can track summer_pop old

double ss_spring_migrant_pop = 0.0;
double ss_autumn_migrant_pop = 0.0;
double ss_spring_staging_size = 0.0;
double ss_autumn_staging_size = 0.0;

struct Individual 
{
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
	
	// individual departure latency
	int latency;
	
	// individual departure timing
	int timing;
	
	// resource cost of migration for individual
	double cost;
	
	// maximum flock size they could have formed part of
	int potential;
	
	// phenology of signalling
	int signal_timing;
	
	// individual age (start out at 1 as they are reproductively mature)
	int age;
};


// wintering ground
Individual WinterPop[N];
Individual StagingPool[N];
Individual SummerPop[N];

// bounds value val between (min and max)
double clamp(double const val, double const min, double const max) 
{
    return(val > max ? max : val < min ? min : val);
}


// get parameters from the command line when 
// running the executable file
void init_arguments(int argc, char **argv)
{
    init_phi_a = atof(argv[1]);  // Elevation of the reaction norm for group size-dependent miratory departure
    init_phi_b = atof(argv[2]);  // Slope of the reaction norm for group size-dependent miratory departure
    init_theta_a = atof(argv[3]);  // // Elevation of the reaction norm for resource-dependent entry to staging pool
    init_theta_b = atof(argv[4]);  // Slope of the reaction norm for resource-dependent entry to staging pool
    pmort = atof(argv[5]);
    pgood_init = atof(argv[6]);
    t_good_ends = atoi(argv[7]);
    rgood_init = atof(argv[8]);
    rbad_init = atof(argv[9]);
    arrival_resource_decay = atof(argv[10]);
    resource_reproduction_threshold = atof(argv[11]);
    resource_starvation_threshold = atof(argv[12]);
    mu_theta = atof(argv[13]);
    mu_phi = atof(argv[14]);
    sdmu_theta = atof(argv[15]);
    sdmu_phi = atof(argv[16]);
    max_migration_cost = atof(argv[17]);
    migration_cost_decay = atof(argv[18]);
    migration_cost_power = atof(argv[19]);
    tmax = atoi(argv[20]);
	twinter = atoi(argv[21]);
	resource_max = atoi(argv[22]);
		
    // some bounds checking on parameters
    // probability of encountering a good environment
    // initially should be 0 <= pgood <= 1
    assert(pgood_init >= 0.0);
    assert(pgood_init <= 1.0);
    
    //max number of days per season > 0
    assert(tmax > 0);
	assert(twinter >= 0);

    // probability that encountering good resource should be
    // set to 0 after t_good_ends timesteps
    assert(t_good_ends >= 0);
    assert(t_good_ends <= tmax);

    // resource increments
    assert(rgood_init > 0);
    assert(rbad_init > 0);
    assert(rbad_init < rgood_init);
	
	assert(resource_max > 0);

}

// write down all parameters in the file
void write_parameters(ofstream &DataFile)  // at end of outputted file
{
    DataFile << endl << endl
            << "init_theta_a;" << init_theta_a << endl
            << "init_theta_b;" << init_theta_b << endl
            << "init_phi_a;" << init_phi_a << endl
            << "init_phi_b;" << init_phi_b << endl
            << "pmort;" << pmort << endl
            << "pgood_init;" << pgood_init << endl
            << "t_good_ends;" << t_good_ends << endl
            << "rgood_init;" << rgood_init << endl
            << "rbad_init;" << rbad_init << endl
            << "arrival_resource_decay;" << arrival_resource_decay << endl
			<< "resource_max;"  << resource_max << endl
            << "resource_reproduction_threshold;" << resource_reproduction_threshold << endl
            << "resource_starvation_threshold;" << resource_starvation_threshold << endl
            << "mu_theta;" << mu_theta << endl
            << "mu_phi;" << mu_phi << endl
            << "sdmu_theta;" << sdmu_theta << endl
            << "sdmu_phi;" << sdmu_phi << endl
            << "tmax;" << tmax << endl
			<< "twinter;" << twinter << endl
            << "N;" << N << endl
            << "migration_cost_decay;" << migration_cost_decay << endl
            << "migration_cost_power;" << migration_cost_power << endl
            << "max_migration_cost;" << max_migration_cost << endl
            << "seed;" << seed << endl;
}

// list of the data headers at the start of the file
void write_data_headers(ofstream &DataFile)
{
    DataFile << "generation;"
        << "time_interval;"
				
		// SPRING MIGRATION STATS:
        << "mean_spring_staging_size;"
        << "var_spring_staging_size;"
        << "spring_migrant_pop;"
		<< "spring_nonmigrant_pop;"
		<< "mean_spring_signal_timing;"
		<< "var_spring_signal_timing;"
		<< "mean_spring_latency;"
		<< "var_spring_latency;"
		<< "mean_spring_departure;"
		<< "var_spring_departure;"
        << "n_spring_flocks;"
        << "mean_spring_flock_size;"
        << "var_spring_flock_size;"
		<< "mean_spring_cost;"
		<< "var_spring_cost;"
		
		// SUMMER STATS:	
		<< "summer_pop;"
	    << "mean_resources_summer;"
	    << "var_resources_summer;"
	    << "mean_theta_a_summer;"
	    << "var_theta_a_summer;"
	    << "mean_theta_b_summer;"
	    << "var_theta_b_summer;"
	    << "mean_phi_a_summer;"
	    << "var_phi_a_summer;"
	    << "mean_phi_b_summer;"
	    << "var_phi_b_summer;"
        << "breeder_pop;"
        << "offspring_pop;"
		
		// AUTUMN MIGRATION STATS
        << "mean_autumn_staging_size;"
        << "var_autumn_staging_size;"
        << "autumn_migrant_pop;"
		<< "autumn_nonmigrant_pop;"
		<< "mean_autumn_signal_timing;"
		<< "var_autumn_signal_timing;"
		<< "mean_autumn_latency;"
		<< "var_autumn_latency;"
		<< "mean_autumn_departure;"
		<< "var_autumn_departure;"
        << "n_autumn_flocks;"
        << "mean_autumn_flock_size;"
        << "var_autumn_flock_size;"
		<< "mean_autumn_cost;"
		<< "var_autumn_cost;"
	
		// WINTER STATS:
		<< "winter_pop;"
		<< "mean_resources_winter;"
	    << "var_resources_winter;"
		<< "mean_theta_a_winter;"
	    << "var_theta_a_winter;"
	    << "mean_theta_b_winter;"
	    << "var_theta_b_winter;"
	    << "mean_phi_a_winter;"
	    << "var_phi_a_winter;"
	    << "mean_phi_b_winter;"
	    << "var_phi_b_winter;"
		<< "mean_age;"
		<< "var_age;"<< endl;
}


// write data for winter population (post mortality)
void write_winter_stats(ofstream &DataFile, int generation, int timestep)
{
    double mean_theta_a[2] = { 0.0, 0.0 };
    double ss_theta_a[2] = { 0.0, 0.0 };
    double mean_theta_b[2] = { 0.0, 0.0 };
    double ss_theta_b[2] = { 0.0, 0.0 };

    double mean_phi_a[2] = { 0.0, 0.0 };
    double ss_phi_a[2] = { 0.0, 0.0 };
    double mean_phi_b[2] = { 0.0, 0.0 };
    double ss_phi_b[2] = { 0.0, 0.0 };

    double mean_resources[2] = { 0.0, 0.0 };
    double ss_resources[2] = { 0.0, 0.0 };
	
	double val;
	
	mean_age = 0.0;
	ss_age = 0.0;
	int age;
    
    for (int i = 0; i < winter_pop; ++i)  // So here we are cycling one by one through the winter population
    {
        // each character (elevation and slope of the 
		// two reaction norms) is genetically controlled
		// by a single gene (diploid) exhibiting incomplete dominance
		// (hence *0.5)
		val = 0.5 * (WinterPop[i].theta_a[0] + WinterPop[i].theta_a[1]);
        mean_theta_a[0] += val;
        ss_theta_a[0] += val * val;

        val = 0.5 * (WinterPop[i].theta_b[0] + WinterPop[i].theta_b[1]);
        mean_theta_b[0] += val;
        ss_theta_b[0] += val * val;
        
        val = 0.5 * (WinterPop[i].phi_a[0] + WinterPop[i].phi_a[1]);
        mean_phi_a[0] += val;
        ss_phi_a[0] += val * val;

        val = 0.5 * (WinterPop[i].phi_b[0] + WinterPop[i].phi_b[1]);
        mean_phi_b[0] += val;
        ss_phi_b[0] += val * val;
        
        val = WinterPop[i].resources;  // the resource level of individual i  // 18 Feb 2020: I've changed this from StagingPool[i].resources
        mean_resources[0] += val;
        ss_resources[0] += val * val;
		
		age = WinterPop[i].age;
		mean_age += age;
		ss_age += age * age;
		
    }

    if (winter_pop > 0)
    {
        // calculate means and variances of the winter population
        mean_theta_a[0] /=  winter_pop;
        mean_theta_b[0] /=  winter_pop;
        mean_phi_a[0] /=  winter_pop;
        mean_phi_b[0] /=  winter_pop;
        mean_resources[0] /=  winter_pop;
		mean_age /= winter_pop;
         
        ss_theta_a[0] /= winter_pop; 
        ss_theta_b[0] /= winter_pop; 
        ss_phi_a[0] /= winter_pop; 
        ss_phi_b[0] /= winter_pop;
        ss_resources[0] /= winter_pop; 
		ss_age /= winter_pop;
    }
	
    // write statistics to a file
    DataFile 
		<< winter_pop << ";"
        << mean_resources[0] << ";"
        << (ss_resources[0] - mean_resources[0] * mean_resources[0]) << ";"
        << mean_theta_a[0] << ";"
        << (ss_theta_a[0] - mean_theta_a[0] * mean_theta_a[0]) << ";"
        << mean_theta_b[0] << ";"
        << (ss_theta_b[0] - mean_theta_b[0] * mean_theta_b[0]) << ";"
        << mean_phi_a[0] << ";"
        << (ss_phi_a[0] - mean_phi_a[0] * mean_phi_a[0]) << ";"
        << mean_phi_b[0] << ";"
        << (ss_phi_b[0] - mean_phi_b[0] * mean_phi_b[0]) << ";"
		<< mean_age << ";"
		<< (ss_age - mean_age * mean_age) << ";"
		<< endl;
// ENDS: write data both for winter population and ends line entry in DataFile
}

void write_summer_stats(ofstream &DataFile, int generation, int timestep)
{
    double mean_theta_a[2] = { 0.0, 0.0 };
    double ss_theta_a[2] = { 0.0, 0.0 };
    double mean_theta_b[2] = { 0.0, 0.0 };
    double ss_theta_b[2] = { 0.0, 0.0 };

    double mean_phi_a[2] = { 0.0, 0.0 };
    double ss_phi_a[2] = { 0.0, 0.0 };
    double mean_phi_b[2] = { 0.0, 0.0 };
    double ss_phi_b[2] = { 0.0, 0.0 };

    double mean_resources[2] = { 0.0, 0.0 };
    double ss_resources[2] = { 0.0, 0.0 };

    double val;
   	
    for (int i = 0; i < summer_pop; ++i)  // for each individual in the summer population:
    {
		val = 0.5 * (SummerPop[i].theta_a[0] + SummerPop[i].theta_a[1]);
        mean_theta_a[1] += val;
        ss_theta_a[1] += val * val;

        val = 0.5 * (SummerPop[i].theta_b[0] + SummerPop[i].theta_b[1]);
        mean_theta_b[1] += val;
        ss_theta_b[1] += val * val;
        
        val = 0.5 * (SummerPop[i].phi_a[0] + SummerPop[i].phi_a[1]);
        mean_phi_a[1] += val;
        ss_phi_a[1] += val * val;

        val = 0.5 * (SummerPop[i].phi_b[0] + SummerPop[i].phi_b[1]);
        mean_phi_b[1] += val;
        ss_phi_b[1] += val * val;
        
        val = SummerPop[i].resources;  // the resource level of individual i 
        mean_resources[1] += val;
        ss_resources[1] += val * val;

	}

    if (summer_pop > 0)
    {
        // calculate means and variances of the summer population
        mean_theta_a[1] /= summer_pop;
        mean_theta_b[1] /= summer_pop;
        mean_phi_a[1] /= summer_pop;
        mean_phi_b[1] /= summer_pop;
        mean_resources[1] /= summer_pop;
        
        ss_theta_a[1] /= summer_pop; 
        ss_theta_b[1] /= summer_pop; 
        ss_phi_a[1] /= summer_pop; 
        ss_phi_b[1] /= summer_pop;
        ss_resources[1] /= summer_pop; 
    }

	
    // write statistics to a file
    DataFile 
        << summer_pop << ";"
	    << mean_resources[1] << ";"
	    << (ss_resources[1] - mean_resources[1] * mean_resources[1]) << ";"
		<< mean_theta_a[1] << ";"
        << (ss_theta_a[1] - mean_theta_a[1] * mean_theta_a[1]) << ";"
        << mean_theta_b[1] << ";"
        << (ss_theta_b[1] - mean_theta_b[1] * mean_theta_b[1]) << ";"
        << mean_phi_a[1] << ";"
        << (ss_phi_a[1] - mean_phi_a[1] * mean_phi_a[1]) << ";"
        << mean_phi_b[1] << ";"
        << (ss_phi_b[1] - mean_phi_b[1] * mean_phi_b[1]) << ";"
		<< breeder_pop << ";"
        << offspring_pop << ";";
// ENDS: write summer stats
}

void write_spring_stats(ofstream &DataFile, int generation, int timestep)
{
	mean_latency = 0.0;
	ss_latency = 0.0;
	mean_departure = 0.0;
	ss_departure = 0.0;
	mean_signal_timing = 0.0;
	ss_signal_timing = 0.0;
	int lat;
	int ticktock;
	int calendar;
   	
    for (int i = 0; i < summer_pop; ++i)  // for each individual in the population of migrants:
    {
		lat = SummerPop[i].latency;  // the migratory latency of individual i
		mean_latency += lat;
		ss_latency += lat * lat;
		
		ticktock = SummerPop[i].timing;
		mean_departure += ticktock;
		ss_departure += ticktock * ticktock;
		
		calendar = SummerPop[i].signal_timing;
		mean_signal_timing += calendar;
		ss_signal_timing += calendar * calendar;
	}

    if (summer_pop > 0)
    {
        mean_latency /= summer_pop;
        ss_latency /= summer_pop;
        mean_departure /= summer_pop;
        ss_departure /= summer_pop;
		mean_signal_timing /= summer_pop;
		ss_signal_timing /= summer_pop;
    }

    // write statistics to a file
    DataFile
        << generation << ";"
        << timestep << ";" 
        << mean_spring_staging_size << ";" 
		<< var_spring_staging_size << ";"
		<< spring_migrant_pop << ";"
		<< spring_nonmigrant_pop << ";"
		<< mean_signal_timing << ";"
		<< (ss_signal_timing - mean_signal_timing * mean_signal_timing) << ";"
		<< mean_latency << ";"
		<< (ss_latency - mean_latency * mean_latency) << ";"
		<< mean_departure << ";"
		<< (ss_departure - mean_departure * mean_departure) << ";"
		<< n_spring_flocks << ";"
		<< mean_spring_flock_size << ";" 
		<< var_spring_flock_size << ";"
		<< mean_spring_cost << ";"
		<< var_spring_cost << ";";

}  // ENDS: write data for spring migrants

void write_autumn_stats(ofstream &DataFile, int generation, int timestep)
{
	mean_latency = 0.0;
	ss_latency = 0.0;
	mean_departure = 0.0;
	ss_departure = 0.0;
	mean_signal_timing = 0.0;
	ss_signal_timing = 0.0;
	int lat;
	int ticktock;
	int calendar;
	
	for (int i = autumn_nonmigrant_pop; i < winter_pop; ++i)	
	{
		lat = WinterPop[i].latency;  // the migratory latency of individual i
		mean_latency += lat;
		ss_latency += lat * lat;
		
		ticktock = WinterPop[i].timing;  // the migratory departure timing of individual i
		mean_departure += ticktock;
		ss_departure += ticktock * ticktock;
		
		calendar = WinterPop[i].signal_timing;  // the signalling phenology of individual i
		mean_signal_timing += calendar;
		ss_signal_timing += calendar * calendar;
	}
	
    if (autumn_migrant_pop > 0)
    {
        mean_latency /= autumn_migrant_pop;
        ss_latency /= autumn_migrant_pop;
        mean_departure /= autumn_migrant_pop;
        ss_departure /= autumn_migrant_pop;
		mean_signal_timing /= autumn_migrant_pop;
		ss_signal_timing /= autumn_migrant_pop;
    }
	
    // write statistics to a file
    DataFile 
        << mean_autumn_staging_size << ";"
		<< var_autumn_staging_size << ";"			
        << autumn_migrant_pop << ";"
		<< autumn_nonmigrant_pop << ";"
		<< mean_signal_timing << ";"
		<< (ss_signal_timing - mean_signal_timing * mean_signal_timing) << ";"
		<< mean_latency << ";"
		<< (ss_latency - mean_latency * mean_latency) << ";"
		<< mean_departure << ";"
		<< (ss_departure - mean_departure * mean_departure) << ";"
		<< n_autumn_flocks << ";"
		<< mean_autumn_flock_size << ";" 
		<< var_autumn_flock_size << ";"
		<< mean_autumn_cost << ";"
		<< var_autumn_cost << ";";
// ENDS: write data both for autumn migrants
}



// initialize the population at the start of the simulation
void init_population()
{
    // loop through all individuals in the wintering ground
    // and give them values
    for (int i = 0; i < N; ++i)
    {
        WinterPop[i].resources = 0.0;  // at start of simulation, initial resource level for all individuals is 0
		
		// set individual latency to 0
		WinterPop[i].latency = 0;
		
		// set individual signal phenology to 1 (Signalling on the first day of the season will give a value of 1)
		WinterPop[i].signal_timing = 1;
		
		// set individual timing value to 1 (Departure on the first day will give a timing value of 1)
		WinterPop[i].timing = 1;
		
        for (int j = 0; j < 2; ++j)
        {
            // initialize allelic values for theta elevation and slope
            WinterPop[i].theta_a[j] = init_theta_a;
            WinterPop[i].theta_b[j] = init_theta_b;
            
            // initialize allelic values for phi elevation and slope
            WinterPop[i].phi_a[j] = init_phi_a;
            WinterPop[i].phi_b[j] = init_phi_b;
			
        }
		
		WinterPop[i].age = 0;
		
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
		else
		{
			WinterPop[i].timing = 1;  // individual survives: timing is reset to 1 for time t+1
			WinterPop[i].signal_timing = 1;
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
		else
		{
			SummerPop[i].timing = 1;  // individual survives: timing is reset to 1 for autumn migration
			SummerPop[i].signal_timing = 1;  // signal phenology is also reset to 1 for autumn
			SummerPop[i].age += 1;
		}
		
    }

    assert((winter_pop > 0 || staging_pop > 0) || summer_pop > 0);
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

// migration cost as a fraction of resources
// wasted on migration
// if flock size is 1, function returns c = 1 (i.e., max cost)
// if flock size is >1, function returns 0 <= c < 1  
double get_migration_cost_proportion(int const flock_size, int const max_flock_size)
{
    // flock size = 1: max cost
    // flock size = N: no cost
    // as this is a proportional cost, it should be bounded between 0 and 1
    // hence costs decay as a fraction of the maximum flock size N
    // if it would just decay with flock size it may well be that 
    double total_migration_cost = max_migration_cost * (1.0 - migration_cost_decay * 
            pow((double) flock_size / max_flock_size, migration_cost_power));

    // make sure function is bounded between 0 and 1
    total_migration_cost = clamp(total_migration_cost, 0.0, 1.0);

    return(total_migration_cost);
}


// the foraging dynamics of the population during winter
void winter_dynamics(int t)
{
	// individuals forage
    // individuals accumulate resources
	// rgood does not run out during winter
	// TODO introduce maximum resourcce value for individuals
	
    for (int i = 0; i < winter_pop; ++i)
    {
		if (uniform(rng_r) < pgood_init) // good resource chosen
        {
            WinterPop[i].resources += rgood;
        }
        else
        {
            WinterPop[i].resources += rbad;
        }
		
	WinterPop[i].resources = clamp(WinterPop[i].resources, 0.0, resource_max); // individuals can reach a (uniform) maximum resource value
    
	} // ok, resource dynamic done

    assert(winter_pop <= N);
    assert(winter_pop >= 0);    

} // ENDS WINTER DYNAMICS (looping through t)


// the dynamics of the population at the wintering ground in spring, time t
void spring_dynamics(int t)
{
	// individuals can continue to forage
    // individuals can continue to accumulate resources
    // individuals make dispersal decisions
	
	// determine probability of encountering a good resource:
    //  if the time is later than t_good_ends
    //  one can only encounter bad resources, hence p_good = 0
    double pgood = t < t_good_ends ? pgood_init : 0.0;
   
    // foraging of individuals who are just at the wintering site
    // and who have yet to decide to go to the staging site
    for (int i = 0; i < winter_pop; ++i)
    {
        WinterPop[i].potential = 0;
		
		if (uniform(rng_r) < pgood) // good resource chosen
        {
            WinterPop[i].resources += rgood;
        }
        else
        {
            WinterPop[i].resources += rbad;
        }
		
	WinterPop[i].resources = clamp(WinterPop[i].resources, 0.0, resource_max);
		 
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
	
	StagingPool[i].resources = clamp(StagingPool[i].resources, 0.0, resource_max);	
		
    } // ENDS staging site foraging loop

    assert(winter_pop <= N);
    assert(winter_pop >= 0);  
    assert((winter_pop > 0 || staging_pop > 0) || summer_pop > 0);  

    double psignal = 0.0;

    // individuals decide whether to go to staging site
    // i.e., prepare for dispersal
    // signal to disperse
    for (int i = 0; i < winter_pop; ++i) 
    {
        
		// individuals must acquire the resources necessary for migration and reproduction before
		// they can consider migrating. 
		if (WinterPop[i].resources < (resource_reproduction_threshold + min_migration_cost))
		{
			continue;
		}
		
		// reaction norm dependent on resources
        // resulting in signaling a willingness to disperse
        // => go to the staging level
        psignal = 0.5 * (WinterPop[i].theta_a[0] + WinterPop[i].theta_a[1])
            + 0.5 * (WinterPop[i].theta_b[0] + WinterPop[i].theta_b[1]) * WinterPop[i].resources;

        // bound the probability
        psignal = clamp(psignal, 0.0, 1.0);

        // does individual want to signal to others to be ready for departure?
        if (uniform(rng_r) < psignal)
        {
            // add individual to the staging pool
            StagingPool[staging_pop] = WinterPop[i];
			StagingPool[staging_pop].latency = 0;
			StagingPool[staging_pop].potential = winter_pop + staging_pop;
            ++staging_pop; // increment the number of individuals in the staging pool

            assert(staging_pop <= N);
            assert(staging_pop > 0);

            // delete this individual from the winter population
            WinterPop[i] = WinterPop[winter_pop - 1];

            // decrement the number of individuals in the winter population
            --winter_pop;
            --i;
        }
		else
		{
			WinterPop[i].signal_timing += 1;  // Individual has not signalled readiness to leave (equivalent to not entering 'staging pool')
			
			WinterPop[i].timing +=1;  // Individuals that do not enter the staging population will not be departing at time t
		}
    } // end for: move dispersers to staging

    // store current number of individuals at the breeding ground
    // so that we know which individuals have just arrived there
    // (we need to update their resources dependent on their migration
    // group size)
    int summer_pop_old = summer_pop;

    // keep track of flock size of individuals who will disperse
    // at this timestep. This should be initialized at 0 at each
    // timestep t
    int NFlock = 0;

    double pdisperse = 0.0;

    // remember numbers of individuals in the staging pop at the start
    int staging_pop_start = staging_pop;
	
    // actual spring dispersal from winter to summer population
    for (int i = 0; i < staging_pop; ++i)
    {
        // later we will consider collective dispersal decisions
        // for now, individuals leave dependent on the current amount of individuals
        // within the staging pool
        pdisperse = 0.5 * (StagingPool[i].phi_a[0] + StagingPool[i].phi_a[1])
            + 0.5 * (StagingPool[i].phi_b[0] + StagingPool[i].phi_b[1]) * (double) staging_pop_start / (staging_pop_start + winter_pop);

        // yes individual goes
        if (uniform(rng_r) < pdisperse)
        {
            SummerPop[summer_pop] = StagingPool[i]; // Individual transfers to SummerPop
            ++summer_pop;
            
            assert(summer_pop <= N);

			// increment the number of individuals recorded as spring migrants
			//spring_migrant_pop = spring_migrant_pop +1;  // Same result by sum below
			
            // delete this individual from the staging population
            StagingPool[i] = StagingPool[staging_pop - 1];

            // decrement the number of individuals in the staging population
            --staging_pop;
            --i;

            assert(staging_pop < N);
            assert(staging_pop >= 0);

            // increase flock size
            ++NFlock;
            
            assert(NFlock <= N);			
        
		} else {
			StagingPool[i].latency += 1;  // If individual does not depart, increment its latency score
			StagingPool[i].timing += 1;  // Also increment its timing score
		}
		
    } // ENDS ACTUAL SPRING DISPERSAL and making of flocks
	
	double total_migration_scalar;

    // keep track of mean and variance in flock sizes
	mean_spring_staging_size += staging_pop_start;

	ss_spring_staging_size += staging_pop_start * staging_pop_start;
	
    if (NFlock > 0)
    {
		n_spring_flocks = n_spring_flocks+1;
		mean_spring_flock_size += NFlock;
		ss_spring_flock_size += NFlock * NFlock;  // Also serves as sum of squares of spring migrant population size
	}
	
	double cost = 0.0;
	
	// update resource levels for all new individuals that have just
    // been added to the summer pool dependent on their flock size
    for (int i = summer_pop_old; i < summer_pop; ++i)  // Selecting individuals that have been added to the summer pop this timepoint
		
    {
		
		total_migration_scalar = (1.0 - get_migration_cost_proportion(NFlock, spring_pop_start)) * 
            (1.0 - arrival_resource_decay * (double)t/tmax);

        assert(total_migration_scalar >= 0);
        assert(total_migration_scalar <= 1.0);
		
		// Resource cost of migration to the individual
		// must be calculated before migration-induced mortality or non-survivors will be excluded
		// THIS IS WHERE WE MAKE A CHANGE (20/04/20), because we have now decided we want to calculate group size based on survivors
		cost = 1 - total_migration_scalar;  // Individual's resource cost scalar
		mean_cost += cost;
		ss_cost += cost * cost;

        // resources are reduced due to migration,
        // yet this depends on group size in a curvilinear fashion
        SummerPop[i].resources = SummerPop[i].resources * total_migration_scalar;
		SummerPop[i].resources = clamp(SummerPop[i].resources, 0.0, resource_max);

		// death due to starvation
        if (SummerPop[i].resources < resource_starvation_threshold)
        {
            SummerPop[i] = SummerPop[summer_pop - 1];
            --summer_pop;
            --i;
        } // ends: death due to starvation
		
    } // ENDS: updating resources of migrants


//    cout << 
//        "mean_spring_flock_size: " << (n_spring_flocks == 0 ? 0 : (double) mean_spring_flock_size / n_spring_flocks) << " ";
//        
//    cout << "summer cost: " << (summer_pop == 0 ? 0 : mean_cost / summer_pop) << " ";
//    cout << "total mean cost: " << mean_cost << " " << endl;
	
} // ENDS SPRING DYNAMICS (looping through t)

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
void create_offspring(
        Individual &mother
        ,Individual &father
        ,Individual &offspring)
{
    bernoulli_distribution allele_sample(0.5);

    offspring.resources = 0.0;
	offspring.latency = 0;
	offspring.timing = 1;
	offspring.signal_timing = 1;
	offspring.age = 0;

    // inherit theta loci

    // each parental allele has probability 0.5 to make it into offspring
    offspring.theta_a[0] = mutation(mother.theta_a[allele_sample(rng_r)], mu_theta, sdmu_theta);
    offspring.theta_a[0] = clamp(offspring.theta_a[0], 0.0, 1.0);

    offspring.theta_a[1] = mutation(father.theta_a[allele_sample(rng_r)], mu_theta, sdmu_theta);
    offspring.theta_a[1] = clamp(offspring.theta_a[1], 0.0, 1.0);

    offspring.theta_b[0] = mutation(mother.theta_b[allele_sample(rng_r)], mu_theta, sdmu_theta);
    offspring.theta_b[0] = clamp(offspring.theta_b[0], 0.0, 1.0);

    offspring.theta_b[1] = mutation(father.theta_b[allele_sample(rng_r)], mu_theta, sdmu_theta);
    offspring.theta_b[1] = clamp(offspring.theta_b[1], 0.0, 1.0);

    // inherit phi loci
    offspring.phi_a[0] = mutation(mother.phi_a[allele_sample(rng_r)], mu_phi, sdmu_phi);
    offspring.phi_a[0] = clamp(offspring.phi_a[0], 0.0, 1.0);

    offspring.phi_a[1] = mutation(father.phi_a[allele_sample(rng_r)], mu_phi, sdmu_phi);
    offspring.phi_a[1] = clamp(offspring.phi_a[1], 0.0, 1.0);
    
    offspring.phi_b[0] = mutation(mother.phi_b[allele_sample(rng_r)], mu_phi, sdmu_phi);
    offspring.phi_b[0] = clamp(offspring.phi_b[0], 0.0, 1.0);

    offspring.phi_b[1] = mutation(father.phi_b[allele_sample(rng_r)], mu_phi, sdmu_phi);
    offspring.phi_b[1] = clamp(offspring.phi_b[1], 0.0, 1.0);

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

    // see if population is extinct. 13 Feb: TODO Need to change this. Population could be alive but no migrants that particular year.
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

        assert(mother.theta_b[1] >= 0.0);
        assert(mother.theta_b[1] <= 1.0);

        // if mom does not meet minimum standards
        // no reproduction through female function
        if (mother.resources < resource_reproduction_threshold)
        {
            continue;  // breaks current iteration in the loop and proceeds to the next one
        }

        // update stats
        ++breeder_pop;

        // now randomly select a father
        do {
            // sample integer uniformly between 0 and summer_pop
            // (not including summer_pop itself)
            father_id = summer_sample(rng_r);
        }
        while (father_id == i);

        father = SummerPop[father_id];
        
        assert(father.theta_b[1] >= 0.0);
        assert(father.theta_b[1] <= 1.0);

        // translate maternal resources to numbers of offspring
        //
        // first round to lowest integer
        resource_integer = floor((mother.resources - resource_reproduction_threshold) / 5);  // 17/04/20: Prior to resetting mothers' resource values, mother.resources was divided by five to reduce family size

        // TODO (slightly digressing): can we come up with an analytical 
        // description of this rounding process of w into integer values?
        if (uniform(rng_r) < mother.resources - resource_integer)
        {
            // make an additional offspring (adding stochasticity to fecundity)
            ++resource_integer;
        }
        
		// reset mother's resource value to zero
		mother.resources = 0;
		
        // for each parent create the number of offspring prescribed by their resource value + noise (i.e, resource_integer)
        for (int kid_i = 0; kid_i < resource_integer; ++kid_i)
        {
            Individual kid;

            create_offspring(mother, father, kid);

            // add kid to the stack
            Kids.push_back(kid);
        }
    } // end for (int i = 0; i < summer_pop; ++i)

    offspring_pop = Kids.size();

    // number of dead individuals is the max population
    // minus the current individuals in the summer population
    // minus the current individuals who stayed at the wintering ground
    int Ndead = N - summer_pop - winter_pop;

    assert(Ndead >= 0);

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

} // ENDS SUMMER REPRODUCTION

// gaining resources at breeding ground
// & fly back
void postbreeding_dynamics(int t)
{
    // determine probability of encountering a good resource:
    //  if the time is later than t_good_ends
    //  one can only encounter bad resources, hence p_good = 0
    double pgood = t < t_good_ends ? pgood_init : 0.0;

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
        SummerPop[i].potential = 0;
		
		if (uniform(rng_r) < pgood) // good resource chosen
        {
            SummerPop[i].resources += rgood;
        }
        else
        {
            SummerPop[i].resources += rbad;
        }
    
	SummerPop[i].resources = clamp(SummerPop[i].resources, 0.0, resource_max);
	
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
    
	StagingPool[i].resources = clamp(StagingPool[i].resources, 0.0, resource_max);
	
	}

    assert(summer_pop<= N);
    assert(summer_pop >= 0);
    assert((summer_pop > 0 || winter_pop > 0) || staging_pop > 0);  // Might have no breeders in a given year, but non-zero population (is that right, Bram?)

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

        // bound the probability
        psignal = clamp(psignal, 0.0, 1.0);

        // does individual want to signal to others to be ready for departure?
        if (uniform(rng_r) < psignal)
        {
            // add individual to the staging pool
            StagingPool[staging_pop] = SummerPop[i];
			StagingPool[staging_pop].latency = 0.0;
			StagingPool[staging_pop].potential = summer_pop + staging_pop;
            ++staging_pop; // increment the number of individuals in the staging pool

            assert(staging_pop <= N);
            assert(staging_pop >= 0);

            // delete this individual from the summer population
            SummerPop[i] = SummerPop[summer_pop - 1];

            // decrement the number of individuals in the summer population
            --summer_pop;
            --i;
        }
		else
		{
			SummerPop[i].signal_timing +=1;  // Individual did not signal at time t
			
			SummerPop[i].timing +=1;  // Individuals that do not enter the staging population will not be departing at time t
		}
    } // end for (int i = 0; i < summer_pop; ++i)

    // store current number of individuals at the wintering ground
    // so that we know which individuals have just arrived there
    // (we need to update their resources dependent on their migration
    // group size)
    int winter_pop_old = winter_pop;

    int NFlock = 0;

    double pdisperse = 0.0;
	mean_latency = 0.0;
	ss_latency = 0.0;
	int lat = 0;

    int staging_pop_start = staging_pop;

    // actual autumn dispersal
    for (int i = 0; i < staging_pop; ++i)
    {
        // later we will consider collective dispersal decisions
        // for now, individuals leave dependent on the current amount of individuals
        // within the staging pool

        pdisperse = 0.5 * (StagingPool[i].phi_a[0] + StagingPool[i].phi_a[1])
            + 0.5 * (StagingPool[i].phi_b[0] + StagingPool[i].phi_b[1]) * (double) staging_pop_start / (staging_pop_start + summer_pop);

        // yes individual goes
        if (uniform(rng_r) < pdisperse)
        {
			
			WinterPop[winter_pop] = StagingPool[i];  // Individual moves from staging pool to first empty position in WinterPop
			lat = WinterPop[winter_pop].latency;
			mean_latency += lat;
			ss_latency += (lat * lat);			
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
		else 
		{
			StagingPool[i].latency += 1;  // Individual does not depart at time t
			StagingPool[i].timing += 1;  // Ditto
		}
	    
    } // ENDS: Autumn dispersal at time t
	
    double total_migration_scalar = 0.0;
	
	mean_autumn_flock_size += NFlock;
	ss_autumn_flock_size += NFlock * NFlock;
	mean_autumn_staging_size += staging_pop_start;
	ss_autumn_staging_size += staging_pop_start * staging_pop_start;
	
    if (winter_pop_old < winter_pop){
		++n_autumn_flocks;
	}

    double cost = 0.0;
	
	// update resource levels for all new individuals that have just
    // been added to the pool dependent on their flock size
    for (int i = winter_pop_old; i < winter_pop; ++i)
    {
		total_migration_scalar = (1.0 - get_migration_cost_proportion(NFlock, autumn_pop_start)) * 
            1.0; //(1.0 - arrival_resource_decay * (double)t/tmax);  13 Feb, SRE: Removed arrival resource decay from wintering ground
		
		// Resource cost of migration to the individual
		cost = 1 - total_migration_scalar;  // Individual's resource cost scalar
		mean_cost += cost;
		ss_cost += cost * cost;

        // resources are reduced due to migration,
        // yet this depends on group size in a curvilinear fashion
        SummerPop[i].resources = SummerPop[i].resources * total_migration_scalar;
		SummerPop[i].resources = clamp(SummerPop[i].resources, 0.0, resource_max);

//        cout << get_migration_cost_proportion(NFlock, winter_pop) << ";" << SummerPop[i].resources << ";" << total_migration_scalar << ";" << NFlock << ";" << endl;

        // death due to starvation
        if (WinterPop[i].resources < resource_starvation_threshold)
        {
            WinterPop[i] = WinterPop[winter_pop - 1];
            --winter_pop;
            --i;
        }

    } // Ends: update resource levels of winter arrivals

//    cout << 
//        "mean_spring_flock_size: " << (n_autumn_flocks == 0 ? 0 : (double) mean_autumn_flock_size / n_autumn_flocks) << " ";
//        
//    cout << "winter cost: " << (winter_pop == 0 ? 0 : mean_cost / winter_pop) << " ";
//    cout << "total mean cost: " << mean_cost << " " << endl;

} // ENDS: POST-BREEDING DYNAMICS 


// THE KEY PART OF THE CODE
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
		mean_spring_cost = 0.0;
		var_spring_cost = 0.0;
		mean_cost = 0.0;
		ss_cost = 0.0;
		spring_pop_start = 0.0;
		autumn_pop_start = 0.0;

        staging_pop = 0.0;  // Set staging population count to zero before winter dynamics
		
		rgood = rgood_init;
		rbad = rbad_init;  
		
		// CHANGING RESOURCE AVAILABILITY MID-SIMULATION (if desired)
		if(generation > number_generations*1)  //EXPERIMENTAL SWITCH
			{	
				rgood = rgood_init/1;
				rbad = rbad_init/1;
			}
		
		if(generation > number_generations*1)  //EXPERIMENTAL SWITCH
			{	
				rgood = rgood_init/1;
				rbad = rbad_init/1;
			}	
		
		spring_pop_start = winter_pop;
		
		
		// Winter foraging (migration is not an option)
		for (int t = 0; t < twinter; ++t)
		{
			winter_dynamics(t);
		}
		
		
		// time during spring
        // during which individuals can migrate (or carry on foraging)
		for (int t = 0; t < tmax; ++t)
        {
            spring_dynamics(t);
			
        }
				
		spring_migrant_pop = mean_spring_flock_size;
			
        // now take averages over all timesteps that individuals did (can) join groups
        mean_spring_flock_size = n_spring_flocks > 0 ? mean_spring_flock_size / n_spring_flocks : 0;
		mean_spring_staging_size /= tmax;
		
		// now record variance in flock size and staging size over the season
		var_spring_flock_size = n_spring_flocks > 0 ? (ss_spring_flock_size / n_spring_flocks) - (mean_spring_flock_size * mean_spring_flock_size) : 0;
		var_spring_staging_size = (ss_spring_staging_size / tmax) - (mean_spring_staging_size * mean_spring_staging_size);	
        
		// migration cost statistics
		mean_spring_cost = spring_migrant_pop > 0 ? mean_cost / spring_migrant_pop : 0;
		var_spring_cost = spring_migrant_pop > 0 ? (ss_cost / spring_migrant_pop) - (mean_spring_cost * mean_spring_cost) : 0;
		
		
		if (generation % skip == 0)
		 {
			 write_spring_stats(DataFile, generation, 1000);
		  }
		
		// all individuals that wanted to migrate have migrated now
        // all remainers are going to stay at wintering ground
		spring_nonmigrant_pop = winter_pop + staging_pop;  
		  
        clear_staging_pool();

        // let individuals die with a certain probability 
        mortality();
		
        // Individuals reproduce after they migrated to the summer spot
        summer_reproduction(DataFile);
        
		if (generation % skip == 0)
		 {
			 write_summer_stats(DataFile, generation, 1000);
		  }
		
        // set autumn migration stats to 0 before postbreeding_dynamics starts
        mean_autumn_flock_size = 0.0;
        mean_autumn_staging_size = 0.0;
		var_autumn_flock_size = 0.0;
		var_autumn_staging_size = 0.0;
		ss_autumn_flock_size = 0.0;
		n_autumn_flocks = 0.0;
		autumn_migrant_pop = 0.0;
		autumn_nonmigrant_pop = 0.0;
		ss_autumn_migrant_pop = 0.0;
		ss_autumn_staging_size = 0.0;
		mean_cost = 0.0;
		ss_cost = 0.0;
		mean_autumn_cost = 0.0;
		var_autumn_cost = 0.0;
		
        autumn_pop_start = summer_pop;
		
		// time during summer during which individuals forage
        for (int t = 0; t < tmax; ++t)
        {
            postbreeding_dynamics(t);
			
        }

		autumn_migrant_pop = mean_autumn_flock_size;
		
		//mean_signal_timing /= autumn_migrant_pop;
		//var_signal_timing = (ss_signal_timing / autumn_migrant_pop) - (autumn_migrant_pop * autumn_migrant_pop);
		
		//mean_latency /= autumn_migrant_pop;  // Denominator is the migrant population size
		//var_latency = (ss_latency / autumn_migrant_pop) - (autumn_migrant_pop * autumn_migrant_pop);
		
        // now take averages over all timesteps that individuals did (can) join groups
        mean_autumn_flock_size = n_autumn_flocks > 0 ?  mean_autumn_flock_size / n_autumn_flocks : 0;
        mean_autumn_staging_size /= tmax;
		
		// now record variance in autumn flock size and staging size over the season
		var_autumn_flock_size = n_autumn_flocks > 0 ? (ss_autumn_flock_size / n_autumn_flocks) - (mean_autumn_flock_size * mean_autumn_flock_size) : 0;
		var_autumn_staging_size = (ss_autumn_staging_size / tmax) - (mean_autumn_staging_size * mean_autumn_staging_size);
		
		// migration cost statistics
		mean_autumn_cost =autumn_migrant_pop > 0 ? mean_cost / autumn_migrant_pop : 0;
		var_autumn_cost = autumn_migrant_pop > 0 ? (ss_cost / autumn_migrant_pop) - (mean_autumn_cost * mean_autumn_cost) : 0;
		
		autumn_nonmigrant_pop = summer_pop + staging_pop;
	  
		if (generation % skip == 0)
		{
			write_autumn_stats(DataFile, generation, 1000);
		}
		
		// all individuals who remain at the summer grounds die
        summer_pop = 0;
        staging_pop = 0;
		breeder_pop = 0;
		summer_pop_old = 0;  // 06/02/20: Again, to track summer_pop_old
		spring_nonmigrant_pop = 0;
		spring_migrant_pop = 0;

        // let individuals die with a certain probability 
        mortality();
		
		if (generation % skip == 0)
        {
            write_winter_stats(DataFile, generation, 1000); 
        }
				
    } // ENDS: GENERATION

    write_parameters(DataFile);
}
