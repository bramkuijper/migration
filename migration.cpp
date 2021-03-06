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
//unsigned int seed =550476143; 
mt19937 rng_r{seed};
uniform_real_distribution<> uniform(0.0,1.0);


// parameters & variables:
// values of the most of these are overridden in the init_arguments()
// function

// number of individuals in population
const int N = 200;  // DEAFULT: 2000

// number of generations
long int number_generations = 200;  // DEFAULT: 1000000

// sampling interval
//int skip = ceil(number_generations / 500);
int skip = 5;


// initial values for phi (social dependency) and theta (resource dependency)
// a is an intercept, b is a gradient
double init_theta_a = 0.0; 
double init_theta_b = 0.0;
double init_phi_a = 0.0;
double init_phi_b = 0.0;

// mortality probability 
double pmort = 0.0;
double relative_mortality_risk_of_migration = 0.0;

// initial probability per season to encounter a good resource patch
double pgood = 0.0;

// 
double patch_consistency_factor = 0.0;

// how much resource individuals obtain on a good vs bad patch
double rgood_init = 0.0;
double rbad_init = 0.0;
double rgood = 0.0;
double rbad = 0.0;
double preparation_penalty = 0.0;

double resource_reproduction_threshold = 0.0;  // minimum resource level necessary to reproduce 
double resource_starvation_threshold = 0.0;  // minimum resource level necessary to survive
double resource_max = 0.0;  // maximum resource value an individual can achieve
double breeding_threshold = 0.0;

// mutation rates
double mu_theta = 0.0;
double mu_phi = 0.0;
double sdmu_theta = 0.0;
double sdmu_phi = 0.0;

// migration cost function
double migration_cost_power = 0.0;
double max_migration_cost = 0.0;
double min_migration_cost = 0.0;
double capacity = 0.0;  // group size at which minimum migration cost is reached. Made this a double rather than an integer so as to avoid problems with division in the migration cost function.
double cost = 0.0;

// reproductive cost function
double min_offspring_cost = 0.0;
double offspring_cost_magnifier = 0.0;

double carryover_proportion = 0.0;  // proportion of an individual's resource value that can be carried over to the following year

// max number of intervals per season (two seasons: summer, winter)
int twinter = 0;
int tspring = 0;

// stats of flock size and staging
double population_mean_spring_flock_size = 0.0;
double individual_mean_spring_flock_size = 0.0;
double mean_spring_staging_size = 0.0;
double population_mean_autumn_flock_size = 0.0;
double individual_mean_autumn_flock_size = 0.0;
double mean_autumn_staging_size = 0.0;
double population_var_spring_flock_size = 0.0;
double individual_var_spring_flock_size = 0.0;
double var_spring_staging_size = 0.0;
double mean_spring_cost = 0.0;
double ss_spring_cost = 0.0;
double population_var_autumn_flock_size = 0.0;
double individual_var_autumn_flock_size = 0.0;
double var_autumn_staging_size = 0.0;
double population_ss_spring_flock_size = 0.0;
double individual_ss_spring_flock_size = 0.0;
double population_ss_autumn_flock_size = 0.0;
double individual_ss_autumn_flock_size = 0.0;
double mean_autumn_cost = 0.0;
double ss_autumn_cost = 0.0;
double mean_latency = 0.0;
double var_latency = 0.0;
double ss_latency = 0.0;
double mean_cost = 0.0;
double var_cost = 0.0;
double ss_cost = 0.0;
double mean_departure_timing = 0.0;
double ss_departure_timing = 0.0;
double mean_signal_timing = 0.0;
double var_signal_timing = 0.0;
double ss_signal_timing = 0.0;
double mean_signal_resources = 0.0;
double var_signal_resources = 0.0;
double ss_signal_resources = 0.0;
double mean_age = 0.0;
double var_age = 0.0;
double ss_age = 0.0;
double mean_summer_cost_breederpop = 0.0;
double ss_summer_cost_breederpop = 0.0;
double mean_fecundity_breederpop = 0.0;
double ss_fecundity_breederpop = 0.0;
double mean_resources = 0.0;
double ss_resources = 0.0;
double rv =0.0;

// keep track of the current number of 
// individuals in various seasons/demographics
int staging_pop = 0;
int winter_pop = 0;
int spring_pop_start = 0;
int spring_nonmigrant_pop = 0;
double spring_migrant_pop = 0;
int spring_migrants_resource_cap = 0;
int spring_migrant_deaths = 0;
int summer_pop = 0;
int autumn_pop_start = 0;
int breeder_pop = 0;
int nonreproductive_pop = 0;
int offspring_pop = 0;
int postbreeding_pop = 0;
int autumn_nonmigrant_pop = 0;
double autumn_migrant_pop = 0;
int autumn_migrants_resource_cap = 0;
int n_spring_flocks = 0;  // recording the number of spring flocks (tspring - n(unusued departure intervals))
int n_autumn_flocks = 0;
int summer_pop_old = 0;  // 06/02/20: So that I can track summer_pop old
int Nvacancies = 0;
int Nsurplus = 0;
int autumn_migrant_deaths = 0;
double autumn_migrant_mortality_rate = 0.0;
int remainer_pop = 0;

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
	
	// size of flock individual was in
	int flock_size;
	
	// resource cost of migration for individual
	double cost;
	
	// phenology of signalling
	int signal_timing;
	
	// resource value when signalling begins
	double signal_resources;
	
	double signalling_proportion;
	
	// individual age (start out at 1 as they are reproductively mature)
	int age;
	
	// number of offspring produced
	int fecundity;
	
	int patch_quality;
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
    pgood = atof(argv[6]);
	patch_consistency_factor = atof(argv[7]);
    rgood_init = atof(argv[8]);
    rbad_init = atof(argv[9]);
	preparation_penalty = atof(argv[10]);
    resource_reproduction_threshold = atof(argv[11]);
    resource_starvation_threshold = atof(argv[12]);
    mu_theta = atof(argv[13]);
    mu_phi = atof(argv[14]);
    sdmu_theta = atof(argv[15]);
    sdmu_phi = atof(argv[16]);
    max_migration_cost = atof(argv[17]);
	min_migration_cost = atof(argv[18]);
    migration_cost_power = atof(argv[19]);
	twinter = atoi(argv[20]);
    tspring = atoi(argv[21]);
	resource_max = atof(argv[22]);
	min_offspring_cost = atof(argv[23]);
	offspring_cost_magnifier = atof(argv[24]);
	carryover_proportion = atof(argv[25]);
	relative_mortality_risk_of_migration = atof(argv[26]);
	capacity = atof(argv[27 ]);
		
    // some bounds checking on parameters
    // probability of encountering a good environment
    // initially should be 0 <= pgood <= 1
    assert(pgood >= 0.0);
    assert(pgood <= 1.0);
    
    //max number of days per season > 0
    assert(tspring > 0);
	assert(twinter >= 0);

    // resource increments
    assert(rgood_init > 0);
    assert(rbad_init > 0);
    assert(rbad_init <= rgood_init);
	
	assert(resource_max > 0);
	
	assert(max_migration_cost >= min_migration_cost);
	assert(capacity <= N);

}

// write down all parameters in the file
void write_parameters(ofstream &DataFile)  // at top of outputted file
{
    DataFile
            << "init_theta_a;" << init_theta_a << endl
            << "init_theta_b;" << init_theta_b << endl
            << "init_phi_a;" << init_phi_a << endl
            << "init_phi_b;" << init_phi_b << endl
            << "pmort;" << pmort << endl
            << "pgood;" << pgood << endl
            << "rgood_init;" << rgood_init << endl
            << "rbad_init;" << rbad_init << endl
			<< "preparation_penalty;" << preparation_penalty << endl
            << "resource_max;"  << resource_max << endl
            << "resource_reproduction_threshold;" << resource_reproduction_threshold << endl
            << "resource_starvation_threshold;" << resource_starvation_threshold << endl
            << "mu_theta;" << mu_theta << endl
            << "mu_phi;" << mu_phi << endl
            << "sdmu_theta;" << sdmu_theta << endl
            << "sdmu_phi;" << sdmu_phi << endl
            << "tspring;" << tspring << endl
			<< "twinter;" << twinter << endl
            << "N;" << N << endl
			<< "number_generations;" << number_generations << endl
            << "migration_cost_power;" << migration_cost_power << endl
            << "max_migration_cost;" << max_migration_cost << endl
			<< "min_migration_cost;" << min_migration_cost << endl
			<< "capacity;" << capacity << endl
			<< "offspring_cost_magnifier;" << offspring_cost_magnifier << endl
			<< "carryover_proportion;" << carryover_proportion << endl
			<< "relative_mortality_risk_of_migration;" << relative_mortality_risk_of_migration << endl
            << "seed;" << seed << endl
			<< endl;
}

// write the distributin of all individuals
// ofstream &DataFile: the distribution file to write it to
// int const generation: the particular generation in which the function is called
// int const factor: a particular number allowing you to distinguish
// between different writes of the distribution within the same generation
void write_dist(ofstream &DataFile, 
        int const generation,
        int const factor)

{
	
    for (int summer_idx = 0; summer_idx < summer_pop; ++summer_idx)
    {
        DataFile << generation << ";"
			<< "summer;"
			<< SummerPop[summer_idx].signal_timing << ";"
			<< SummerPop[summer_idx].signal_resources << ";"
            << SummerPop[summer_idx].timing << ";"
			<< SummerPop[summer_idx].signalling_proportion << ";"
            << SummerPop[summer_idx].latency << ";"
			<< SummerPop[summer_idx].flock_size << ";"
           	<< SummerPop[summer_idx].cost << ";"
			<< SummerPop[summer_idx].resources << ";"  // At arrival
	        << SummerPop[summer_idx].theta_a[0]*0.5 + SummerPop[summer_idx].theta_a[1]*0.5 << ";"
	        << SummerPop[summer_idx].theta_b[0]*0.5 + SummerPop[summer_idx].theta_b[1]*0.5 << ";"
	        << SummerPop[summer_idx].phi_a[0]*0.5 + SummerPop[summer_idx].phi_a[1]*0.5 << ";"
	        << SummerPop[summer_idx].phi_b[0]*0.5 + SummerPop[summer_idx].phi_b[1]*0.5 << ";"
            << SummerPop[summer_idx].age << ";"
			<< SummerPop[summer_idx].patch_quality << ";" << std::endl;
	    }
} // ENDS: ()

void write_dist_data_headers(ofstream &DataFile)
{
    DataFile << "generation;"
//        << "factor;" // allows you to distinguish between multiple calls of ()
        << "season;"
		<< "signal_timing;"	
		<< "signal_resources;"
		<< "departure_timing;"
		<< "signalling_proportion;"
		<< "latency;"
		<< "flock_size;"
		<< "cost;"
		<< "arrival_resources;"
        << "theta_a;"
        << "theta_b;"
        << "phi_a;"
        << "phi_b;"
        << "age;"
		<< "patch_quality;" << std::endl;
}

// list of the data headers 
void write_data_headers(ofstream &DataFile)
{
    // SPRING MIGRATION STATS (n = 24):
	DataFile << "generation;"  // 1		
		<< "spring_pop;"  // 2
		<< "mean_spring_staging_size;"  // 3
        << "var_spring_staging_size;"  // 4
        << "spring_migrant_pop;"  // 5
		<< "spring_migrants_resource_cap;"  // 6
		<< "spring_nonmigrant_pop;"  // 7
		<< "mean_spring_signal_resources;"  // 8
		<< "var_spring_signal_resources;"  // 9
		<< "mean_spring_signal_timing;"  // 10
		<< "var_spring_signal_timing;"  // 11
		<< "mean_spring_latency;"
		<< "var_spring_latency;"  // 13
		<< "mean_spring_departure_timing;"
		<< "var_spring_departure_timing;"  // 15
		<< "mean_spring_departure_resources;"  // 16
		<< "var_spring_departure_resources;"  // 17
        << "n_spring_flocks;"
        << "population_mean_spring_flock_size;"
        << "population_var_spring_flock_size;"  // 20
	    << "individual_mean_spring_flock_size;"  //21
	    << "individual_var_spring_flock_size;"  // 22
		<< "mean_spring_cost;"  // 23
		<< "var_spring_cost;"  // 24
		
		// SUMMER STATS (10):	
		<< "spring_migrant_mortality_rate;"
		<< "summer_pop;"
	    << "breeder_pop;"
		<< "nonreproductive_pop;"
        << "offspring_pop;"
		<< "mean_reproductive_cost_breederpop;"  // 6
		<< "var_reproductive_cost_breederpop;"
		<< "mean_fecundity_breederpop;"
		<< "var_fecundity_breederpop;"
		<< "postbreeding_pop;"
		
		// AUTUMN MIGRATION STATS (22)
        << "mean_autumn_staging_size;"
        << "var_autumn_staging_size;"
        << "autumn_migrant_pop;"
		<< "autumn_migrants_resource_cap;"
		<< "autumn_nonmigrant_pop;"
		<< "mean_autumn_signal_resources;"
		<< "var_autumn_signal_resources;"
		<< "mean_autumn_signal_timing;"
		<< "var_autumn_signal_timing;"
		<< "mean_autumn_latency;"
		<< "var_autumn_latency;"
		<< "mean_autumn_departure_timing;"
		<< "var_autumn_departure_timing;"
		<< "mean_autumn_departure_resources;"
		<< "var_autumn_departure_resources;"
        << "n_autumn_flocks;"
        << "population_mean_autumn_flock_size;"
        << "population_var_autumn_flock_size;"
	    << "individual_mean_autumn_flock_size;"
	    << "individual_var_autumn_flock_size;"
		<< "mean_autumn_cost;"
		<< "var_autumn_cost;"
	
		// WINTER STATS (17):
		<< "autumn_migrant_deaths;"	
		<< "autumn_migrant_mortality_rate;"
		<< "total_autumn_mortality_rate;"	
		<< "remainer_pop;"
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
		<< "var_age;"
		<< endl;
}


// write data for winter population (post mortality)
void write_winter_stats(ofstream &DataFile)
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
        mean_theta_a += val;
        ss_theta_a += val * val;

        val = 0.5 * (WinterPop[i].theta_b[0] + WinterPop[i].theta_b[1]);
        mean_theta_b += val;
        ss_theta_b += val * val;
        
        val = 0.5 * (WinterPop[i].phi_a[0] + WinterPop[i].phi_a[1]);
        mean_phi_a += val;
        ss_phi_a += val * val;

        val = 0.5 * (WinterPop[i].phi_b[0] + WinterPop[i].phi_b[1]);
        mean_phi_b += val;
        ss_phi_b += val * val;
        
        val = WinterPop[i].resources;  // the resource level of individual i  // 18 Feb 2020: I've changed this from StagingPool[i].resources
        mean_resources += val;
        ss_resources += val * val;
		
		age = WinterPop[i].age;
		mean_age += age;
		ss_age += age * age;
		
    }

    if (winter_pop > 0)
    {
        // calculate means and variances of the winter population
        mean_theta_a /=  winter_pop;
        mean_theta_b /=  winter_pop;
        mean_phi_a /=  winter_pop;
        mean_phi_b /=  winter_pop;
        mean_resources /=  winter_pop;
		mean_age /= winter_pop;
         
        ss_theta_a /= winter_pop; 
        ss_theta_b /= winter_pop; 
        ss_phi_a /= winter_pop; 
        ss_phi_b /= winter_pop;
        ss_resources /= winter_pop; 
		ss_age /= winter_pop;
    }
	
    // write statistics to a file
    DataFile 
		<< autumn_migrant_deaths << ";"
		<< (autumn_migrant_deaths / autumn_migrant_pop) << ";"
		<< (autumn_nonmigrant_pop + autumn_migrant_deaths) / (autumn_nonmigrant_pop + autumn_migrant_pop ) << ";"
		<< remainer_pop << ";"
		<< winter_pop << ";"
        << mean_resources << ";"
        << (ss_resources - mean_resources * mean_resources) << ";"
        << mean_theta_a << ";"
        << (ss_theta_a - mean_theta_a * mean_theta_a) << ";"
        << mean_theta_b << ";"
        << (ss_theta_b - mean_theta_b * mean_theta_b) << ";"
        << mean_phi_a << ";"
        << (ss_phi_a - mean_phi_a * mean_phi_a) << ";"
        << mean_phi_b << ";"
        << (ss_phi_b - mean_phi_b * mean_phi_b) << ";"
		<< mean_age << ";"
		<< (ss_age - mean_age * mean_age) << ";"
		<< endl;
// ENDS: write data both for winter population and ends line entry in DataFile
}

void write_summer_stats(ofstream &DataFile)
{
    double val;
		
	double mean_summer_cost = 0.0;
	double ss_summer_cost = 0.0;
	
	double mean_fecundity = 0.0;
	double ss_fecundity = 0.0;
		
    for (int i = 0; i < (breeder_pop + nonreproductive_pop); ++i)  // for each individual in the summer population:
    {
		
		val = SummerPop[i].cost;
		mean_summer_cost += val;
		ss_summer_cost += val * val;
		
		val = SummerPop[i].fecundity;
		mean_fecundity += val;
		ss_fecundity += val * val;
	}

    if (summer_pop > 0)
    {
        // calculate means and variances of the summer population
        mean_summer_cost /= breeder_pop;
		mean_fecundity /= breeder_pop;
        ss_summer_cost /= breeder_pop;
		ss_fecundity /= breeder_pop;
    }
	
    // write statistics to a file
    DataFile 
        << (spring_migrant_deaths / spring_migrant_pop) << ";"		
		<< (breeder_pop + nonreproductive_pop) << ";"
	    << breeder_pop << ";"
		<< nonreproductive_pop << ";"
        << offspring_pop << ";"
		<< mean_summer_cost << ";"
		<< (ss_summer_cost - mean_summer_cost * mean_summer_cost) << ";"
		<< mean_fecundity << ";"
		<< (ss_fecundity - mean_fecundity * mean_fecundity) << ";"
		<< postbreeding_pop << ";";

}  // ENDS: write summer stats


void write_spring_stats(ofstream &DataFile, int generation)
{
	mean_latency = 0.0;
	ss_latency = 0.0;
	mean_departure_timing = 0.0;
	ss_departure_timing = 0.0;
	mean_signal_timing = 0.0;
	ss_signal_timing = 0.0;
	individual_mean_spring_flock_size = 0.0;
	individual_ss_spring_flock_size = 0.0;
	mean_spring_cost = 0.0;
	ss_spring_cost = 0.0;
	mean_signal_resources = 0.0;
	ss_signal_resources = 0.0;
	int lat;
	int ticktock;
	int calendar;
	int fatness;
	int group_size;
	double spring_cost;
   	
    for (int i = 0; i < summer_pop; ++i)  // for each individual in the population of migrants:
    {
		lat = SummerPop[i].latency;  // the migratory latency of individual i
		mean_latency += lat;
		ss_latency += lat * lat;
		
		ticktock = SummerPop[i].timing;
		mean_departure_timing += ticktock;
		ss_departure_timing += ticktock * ticktock;
		
		calendar = SummerPop[i].signal_timing;
		mean_signal_timing += calendar;
		ss_signal_timing += calendar * calendar;
		
		fatness = SummerPop[i].signal_resources;
		mean_signal_resources += fatness;
		ss_signal_resources += fatness * fatness;
		
		group_size = SummerPop[i].flock_size;
		individual_mean_spring_flock_size += group_size;
		individual_ss_spring_flock_size += group_size * group_size;
		
		spring_cost = SummerPop[i].cost;
		mean_spring_cost += spring_cost;
		ss_spring_cost += spring_cost * spring_cost;
	}

    if (summer_pop > 0)
    {
        mean_latency /= summer_pop;
        ss_latency /= summer_pop;
        mean_departure_timing /= summer_pop;
        ss_departure_timing /= summer_pop;
		mean_signal_timing /= summer_pop;
		ss_signal_timing /= summer_pop;
		individual_mean_spring_flock_size /= summer_pop;
		individual_ss_spring_flock_size /= summer_pop;
		mean_spring_cost /= summer_pop;
		ss_spring_cost /= summer_pop;
		mean_resources /= summer_pop;
		ss_resources /= summer_pop;
		mean_signal_resources /= summer_pop;
		ss_signal_resources /= summer_pop;
    }

    // write statistics to a file
    DataFile
        << generation << ";"  // 1
		<< spring_pop_start << ";"
		<< mean_spring_staging_size << ";"  // 3
		<< var_spring_staging_size << ";"
		<< spring_migrant_pop << ";"  // 5
		<< spring_migrants_resource_cap << ";"
		<< spring_nonmigrant_pop << ";"  // 7
		<< mean_signal_resources << ";"
		<< (ss_signal_resources - mean_signal_resources * mean_signal_resources) << ";"  // 9
		<< mean_signal_timing << ";"
		<< (ss_signal_timing - mean_signal_timing * mean_signal_timing) << ";"  // 11
		<< mean_latency << ";"
		<< (ss_latency - mean_latency * mean_latency) << ";"  // 13
		<< mean_departure_timing << ";"
		<< (ss_departure_timing - mean_departure_timing * mean_departure_timing) << ";"  // 15
		<< mean_resources << ";"
		<< (ss_resources - mean_resources * mean_resources) << ";"  // 17
		<< n_spring_flocks << ";"
		<< population_mean_spring_flock_size << ";"  // 19
		<< population_var_spring_flock_size << ";"
		<< individual_mean_spring_flock_size << ";"  // 21
		<< (individual_ss_spring_flock_size - individual_mean_spring_flock_size * individual_mean_spring_flock_size) << ";"	
		<< mean_spring_cost << ";"  // 23
		<< (ss_spring_cost - mean_spring_cost * mean_spring_cost) << ";";  // 24

}  // ENDS: write data for spring migrants

void write_autumn_stats(ofstream &DataFile)
{
	mean_latency = 0.0;
	ss_latency = 0.0;
	mean_departure_timing = 0.0;
	ss_departure_timing = 0.0;
	mean_signal_timing = 0.0;
	ss_signal_timing = 0.0;
	individual_mean_autumn_flock_size = 0.0;
	individual_ss_autumn_flock_size = 0.0;
	mean_autumn_cost = 0.0;
	ss_autumn_cost = 0.0;
	mean_signal_resources = 0.0;
	ss_signal_resources = 0.0;
	int lat;
	int ticktock;
	int calendar;
	int fatness;
	int group_size;
	double autumn_cost;
	
	for (int i = remainer_pop; i < winter_pop; ++i)	
	{
		lat = WinterPop[i].latency;  // the migratory latency of individual i
		mean_latency += lat;
		ss_latency += lat * lat;
		
		ticktock = WinterPop[i].timing;  // the migratory departure timing of individual i
		mean_departure_timing += ticktock;
		ss_departure_timing += ticktock * ticktock;
		
		calendar = WinterPop[i].signal_timing;  // the signalling phenology of individual i
		mean_signal_timing += calendar;
		ss_signal_timing += calendar * calendar;
		
		fatness = WinterPop[i].signal_resources;
		mean_signal_resources += fatness;
		ss_signal_resources += fatness * fatness;
		
		group_size = WinterPop[i].flock_size;
		individual_mean_autumn_flock_size += group_size;
		individual_ss_autumn_flock_size += group_size * group_size;
		
		autumn_cost = WinterPop[i].cost;
		mean_autumn_cost += autumn_cost;
		ss_autumn_cost += autumn_cost * autumn_cost;
	}
	
    if (autumn_migrant_pop > 0)
    {
        mean_latency /= autumn_migrant_pop;
        ss_latency /= autumn_migrant_pop;
        mean_departure_timing /= autumn_migrant_pop;
        ss_departure_timing /= autumn_migrant_pop;
		mean_resources /= autumn_migrant_pop;
		ss_resources /= autumn_migrant_pop;
		mean_signal_timing /= autumn_migrant_pop;
		ss_signal_timing /= autumn_migrant_pop;
		mean_signal_resources /= autumn_migrant_pop;
		ss_signal_resources /= autumn_migrant_pop;
		individual_mean_autumn_flock_size /= autumn_migrant_pop;
		individual_ss_autumn_flock_size /= autumn_migrant_pop;
		mean_autumn_cost /= autumn_migrant_pop;
		ss_autumn_cost /= autumn_migrant_pop;
    }
	
    // write statistics to a file
    DataFile 
        << mean_autumn_staging_size << ";"
		<< var_autumn_staging_size << ";"			
        << autumn_migrant_pop << ";"
		<< autumn_migrants_resource_cap << ";"
		<< autumn_nonmigrant_pop << ";"
		<< mean_signal_resources << ";"
		<< (ss_signal_resources - mean_signal_resources * mean_signal_resources) << ";"
		<< mean_signal_timing << ";"
		<< (ss_signal_timing - mean_signal_timing * mean_signal_timing) << ";"
		<< mean_latency << ";"
		<< (ss_latency - mean_latency * mean_latency) << ";"
		<< mean_departure_timing << ";"
		<< (ss_departure_timing - mean_departure_timing * mean_departure_timing) << ";"
		<< mean_resources << ";"
		<< (ss_resources - mean_resources * mean_resources) << ";"
		<< n_autumn_flocks << ";"
		<< population_mean_autumn_flock_size << ";" 
		<< population_var_autumn_flock_size << ";"
		<< individual_mean_autumn_flock_size << ";" 
		<< (individual_ss_autumn_flock_size - individual_mean_autumn_flock_size * individual_mean_autumn_flock_size) << ";"
		<< mean_autumn_cost << ";"
		<< (ss_autumn_cost - mean_autumn_cost * mean_autumn_cost) << ";";
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
		
		WinterPop[i].signal_resources = 0.0;
		
		// set individual latency to 0
		WinterPop[i].latency = 0;
		
		// set individual signal phenology to 1 (Signalling on the first day of the season will give a value of 1)
		WinterPop[i].signal_timing = 1;
		
		// set individual timing value to 1 (Departure on the first day will give a timing value of 1)
		WinterPop[i].timing = 1;
		
		WinterPop[i].signalling_proportion = 0.0;
		
		WinterPop[i].flock_size = 0.0;
		
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
		
		WinterPop[i].cost = 0.0; 
		
		WinterPop[i].fecundity = 0; 
		
		if (uniform(rng_r) < pgood)
		{
			WinterPop[i].patch_quality = 1;  // Will secure high quality foraging
		}
		else
		{
			WinterPop[i].patch_quality = 0;  // Will access low quality foraging
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

double migration_cost(int NFlock)
{
	//cost = min_migration_cost + ((max_migration_cost - min_migration_cost) / pow(NFlock, migration_cost_power));
	// Flock size of 1 (i.e., travelling alone) yields the maximum cost
	
	// 17 Nov 2020: Revising the cost function to allow it to take an accelerating / decelarating function
	if (NFlock >= capacity)
	{
		cost = min_migration_cost;
	}
	else{
		//cost = min_migration_cost + (capacity - NFlock)*pow(capacity, -1);	
		cost = min_migration_cost + ((max_migration_cost - min_migration_cost) * pow(1 - ((NFlock-1) / (capacity-1)), exp(migration_cost_power)));
	}
	return(cost);
}  // ENDS: migration cost function

void spring_mortality()
{
    // NON-MIGRANTS
	for (int i = 0; i < winter_pop;++i)
    {
        // random mortality of non-migrants
        if (uniform(rng_r) < pmort / relative_mortality_risk_of_migration)
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
	
	// MIGRANTS
    for (int i = 0; i < summer_pop;++i)
    {
		// migration-induced starvation
        if (SummerPop[i].resources < resource_starvation_threshold)
        {
            SummerPop[i] = SummerPop[summer_pop - 1];
            --summer_pop;
            --i;
			
			++spring_migrant_deaths;
        } // ends: death due to starvation
		
		// random mortality of migrants
        else if (uniform(rng_r) < 1-sqrt(1 - pmort))
        {
            SummerPop[i] = SummerPop[summer_pop - 1];
            --summer_pop;
            --i;
			
			++spring_migrant_deaths;
        }
		
		else
		{
			SummerPop[i].timing = 1;  // individual survives: timing is reset to 1 for autumn migration
			SummerPop[i].signal_timing = 1;  // signal phenology is also reset to 1 for autumn
		}
    }
}

void autumn_mortality()
{
    
	// MIGRANTS
	for (int i = remainer_pop; i < winter_pop; ++i)
		
		// migration-induced mortality
		if (WinterPop[i].resources < resource_starvation_threshold)        
			{            
				WinterPop[i] = WinterPop[winter_pop - 1];           
				--winter_pop;            
				--i;
				
				++autumn_migrant_deaths;
			}
				
		//random mortality
	    else if (uniform(rng_r) < 1 - sqrt(1 - pmort))
	        {
	            WinterPop[i] = WinterPop[winter_pop - 1];
	            --winter_pop;
	            --i;
				
				++autumn_migrant_deaths;
	        } 
			
		else
			{
				WinterPop[i].timing = 1;  // individual survives: timing is reset to 1, ready for spring phenology monitoring
				WinterPop[i].signal_timing = 1;
				WinterPop[i].age += 1;  // individuals become a year older during winter, so that first-year spring migrants are scored as one year-olds rather than 0 year-olds
			} 
	
	// NON-MIGRANTS MORTALITY WAS DEALT WITH IN SPRING_MORTALITY()
	// (they effectively have two bouts of mortality in one go
	// ('pmort' being the annual mortality probability, 
	// modified by 'relative_mortality_risk')
	
}

// remove individuals from the staging pool and put them
// back in the original population
void clear_staging_pool()
{
	// put individuals still in the staging pool (i.e., those that signalled but didn't depart) back in the original population
    for (int i = 0; i < staging_pop; ++i)
    {
        WinterPop[winter_pop] = StagingPool[i];
		winter_pop++;
    }

    // just double check that winter_pop does not exceed max population size
	assert(winter_pop <= N);
	assert(winter_pop >= staging_pop);
	
    staging_pop = 0;
	spring_nonmigrant_pop = winter_pop;
	
}  // ENDS STAGING POOL CLEARANCE

// the foraging dynamics of the population during winter
void winter_dynamics(int t)
{
	// individuals forage
    // individuals accumulate resources
	// rgood does not run out during winter
	
    for (int i = 0; i < winter_pop; ++i)
    {
		if (WinterPop[i].patch_quality == 1)
		{
			WinterPop[i].resources += rgood;
		}
		else
		{
			WinterPop[i].resources += rbad;
		}
		
		// possible change of foraging success in next timestep:	
		if (uniform(rng_r) < 1/pow(10, patch_consistency_factor))
        {
            WinterPop[i].patch_quality = 1 - WinterPop[i].patch_quality;  // switch of patch quality (good to poor; poor to good)
        }
        
		
	WinterPop[i].resources = min(WinterPop[i].resources, resource_max); // individuals can reach a (uniform) maximum resource value
    
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
	
	// foraging of individuals who are just at the wintering site
    // and who have yet to decide to go to the staging site
    for (int i = 0; i < winter_pop; ++i)
    {		
		if (WinterPop[i].patch_quality == 1)
		{
			WinterPop[i].resources += rgood;
		}
		else
		{
			WinterPop[i].resources += rbad;
		}
		
		// possible change of foraging success in next timestep:	
		if (uniform(rng_r) < 1/pow(10, patch_consistency_factor))
        {
            WinterPop[i].patch_quality = 1 - WinterPop[i].patch_quality;  // switch of patch quality (good to poor; poor to good)
        }
		
	WinterPop[i].resources = min(WinterPop[i].resources, resource_max);
		 
    } // ok, resource dynamic done



    // foraging of individuals who are already at the staging site
    for (int i = 0; i < staging_pop; ++i)
    { 
        // indivuals who are already at the staging site
        // continue to forage at the staging site
		if (StagingPool[i].patch_quality == 1)
		{
			StagingPool[i].resources += rgood * preparation_penalty;
		}
		else
		{
			StagingPool[i].resources += rbad * preparation_penalty;
		}
		
		// possible change of foraging success in next timestep:	
		if (uniform(rng_r) < 1/pow(10, patch_consistency_factor))
        {
            StagingPool[i].patch_quality = 1 - StagingPool[i].patch_quality;  // switch of patch quality (good to poor; poor to good)
        }
	
	StagingPool[i].resources = min(StagingPool[i].resources, resource_max);	
		
    } // ENDS staging site foraging loop

    assert(winter_pop <= N);
    assert(winter_pop >= 0);  
    assert(winter_pop > 0 || staging_pop > 0 || summer_pop > 0);  
	assert(winter_pop + staging_pop + summer_pop <= N);  
    double psignal = 0.0;

    // individuals decide whether to go to staging site
    // i.e., prepare for dispersal
    // signal to disperse
    for (int i = 0; i < winter_pop; ++i) 
    {	
		// reaction norm dependent on resources
        // resulting in signaling a willingness to disperse
        // => go to the staging level
		
		// SIGMOIDAL MODEL
		// of the form psignal = (1 + e^-theta_b.(resources - theta_a))^-1
		psignal = pow(1 + exp(-0.5 * (WinterPop[i].theta_b[0] + WinterPop[i].theta_b[1]) 
			* (WinterPop[i].resources - 0.5 * (WinterPop[i].theta_a[0] + WinterPop[i].theta_a[1]))), -1);

        // bound the probability
        psignal = clamp(psignal, 0, 1);

        // does individual want to signal to others to be ready for departure?
        if (uniform(rng_r) < psignal)
        {
            // add individual to the staging pool
            StagingPool[staging_pop] = WinterPop[i];
			StagingPool[staging_pop].latency = 0;
			StagingPool[staging_pop].signal_resources = StagingPool[staging_pop].resources;
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
        assert(staging_pop < N);
		
		// later we will consider collective dispersal decisions
        // for now, individuals leave dependent on the current amount of individuals
        // within the staging pool
        
		// SIGMOIDAL MODEL
		pdisperse = pow(1 + exp(-0.5 * (StagingPool[i].phi_b[0] + StagingPool[i].phi_b[1]) 
			* (((double) staging_pop_start / (staging_pop_start + winter_pop)) - 0.5 * (StagingPool[i].phi_a[0] + StagingPool[i].phi_a[1]))), -1);
		
		// bound the probability (not really necessary)
        pdisperse = clamp(pdisperse, 0, 1);

        // yes individual goes
        if (uniform(rng_r) < pdisperse)
        {
			SummerPop[summer_pop] = StagingPool[i]; // Individual transfers to SummerPop
			SummerPop[summer_pop].signalling_proportion = (double) staging_pop_start / (staging_pop_start + winter_pop);
			++summer_pop;
            
            assert(summer_pop <= N);

			if (StagingPool[i].resources == resource_max)
			{
				++spring_migrants_resource_cap;
			}

			rv = StagingPool[i].resources;
			mean_resources += rv;
			ss_resources += rv * rv;
			
            // delete this individual from the staging population
            StagingPool[i] = StagingPool[staging_pop - 1];

            // decrement the number of individuals in the staging population
            --staging_pop;
            --i;

            assert(staging_pop >= 0);

            // increase flock size
            ++NFlock;
            
            assert(NFlock <= N);			
        
		} 
		
		else {
			StagingPool[i].latency += 1;  // If individual does not depart, increment its latency score
			StagingPool[i].timing += 1;  // Also increment its timing score
		}
		
    } // ENDS ACTUAL SPRING DISPERSAL and making of flocks

    // keep track of mean and variance in flock sizes
	mean_spring_staging_size += staging_pop_start;

	ss_spring_staging_size += staging_pop_start * staging_pop_start;
	
    if (NFlock > 0)
    {
		n_spring_flocks = n_spring_flocks+1;
		population_mean_spring_flock_size += NFlock;
		population_ss_spring_flock_size += NFlock * NFlock;  // Also serves as sum of squares of spring migrant population size
	}
	
	// update resource levels for all new individuals that have just
    // been added to the summer pool dependent on their flock size
    for (int i = summer_pop_old; i < summer_pop; ++i)  // Selecting individuals that have been added to the summer pop this timepoint
		
    {		
		// Resource cost of migration to the individual
		SummerPop[i].flock_size = NFlock;
		SummerPop[i].cost = migration_cost(NFlock);
		SummerPop[i].resources -= migration_cost(NFlock);
		SummerPop[i].resources = min(SummerPop[i].resources, resource_max);
		SummerPop[i].fecundity = 0;
		
    } // ENDS: updating resources of migrants

} // ENDS SPRING DYNAMICS (looping through t)

// mutation of a certain allele with value val
// given mutation rate mu and mutational distribution stdev sdmu
double mutation(double val, double mu, double sdmu)
{
    if (uniform(rng_r) < mu)
    {
        normal_distribution<> allelic_dist(0,sdmu*val);
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
	offspring.signal_resources = 0.0;
	offspring.latency = 0;
	offspring.timing = 1;
	offspring.signal_timing = 1;
	offspring.age = 0;
	offspring.fecundity = 0;
	offspring.signalling_proportion = 0.0;
	offspring.flock_size = 0;

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
	
}  // ENDS OFFSPRING PRODUCTION



// in summary, they reproduce dependent on 
// resources and arrival time
void summer_reproduction(ofstream &DataFile)
{
	// the next three lines are now duplicated in the main text of the model at the bottom of the script, so I think this is redundant
	mean_resources = 0.0;
	ss_resources = 0.0;
	rv = 0.0;
	
	// auxiliary variables storing current mom and dad
    Individual mother, father;

    // auxiliary variable specifing the rounded amount of a mother's
    // resources
    int resource_integer = 0;
	double offspring_equivalence = 0.0;

    /// auxilary variable specifying the identity of the randomly sampled
    // father
    int father_id;

    // use a flexible array for the kids
    vector<Individual> Kids;

    uniform_int_distribution<> summer_sample(0, summer_pop - 1);

    // Mating dynamic. Presumes that there are an even number of individuals so we  
    // just discard the last individual
    for (int i = 0; i < summer_pop; ++i)
    {

        // if individual does not meet minimum standards then no reproduction through female function
        if (SummerPop[i].resources < breeding_threshold * (tspring + SummerPop[i].timing - 1) / tspring)  // Cost of clutch size of one. Further offspring incur a smaller, incremental cost that also increases through the season (below)
        {
            SummerPop[i].fecundity = 0.0;
			SummerPop[i].cost = 0.0;
			
			++nonreproductive_pop;  // Tally of non-reproductive adults.
			
			//continue;  // breaks current iteration in the loop and proceeds to the next one
        }  // Closes IF loop
		
		else 
		{
			
	        // update stats
	        ++breeder_pop;
		
			rv = SummerPop[i].resources;
			mean_resources += rv;
			ss_resources += rv * rv;
		
		    // get the mother
	        mother = SummerPop[i];

	        //assert(mother.theta_b[1] >= 0.0);
	        //assert(mother.theta_b[1] <= 1.0);

	        // now randomly select a father
	        do {
	            // sample integer uniformly between 0 and summer_pop
	            // (not including summer_pop itself)
	            father_id = summer_sample(rng_r);
	        }
	        while (father_id == i);

	        father = SummerPop[father_id];
        
	        //assert(father.theta_b[1] >= 0.0);
	        //assert(father.theta_b[1] <= 1.0);

	        // translate maternal resources to numbers of offspring
	        //
	        // first round to lowest integer (+ 1 to account for first offspring, cost of which is represented by the phenologically-sensitive breeding_threshold)
	        //resource_integer = 1 + floor((SummerPop[i].resources - (breeding_threshold * (tspring + SummerPop[i].timing - 1) / tspring)) / (min_offspring_cost + ((SummerPop[i].timing - 1) * offspring_cost_magnifier * min_offspring_cost / tspring)));
	        offspring_equivalence = 1 + (SummerPop[i].resources - (breeding_threshold * (tspring + SummerPop[i].timing - 1) / tspring)) / (min_offspring_cost + (min_offspring_cost * offspring_cost_magnifier * (SummerPop[i].timing - 1) / tspring));
			resource_integer = floor(offspring_equivalence);
		
	        // TODO (slightly digressing): can we come up with an analytical 
	        // description of this rounding process of w into integer values?
	        if (uniform(rng_r) < offspring_equivalence - resource_integer)
	        {
	            // make an additional offspring (adding stochasticity to fecundity)
	            ++resource_integer;
	        }
        
			// record resources invested into reproduction
			SummerPop[i].cost = SummerPop[i].resources;
		
			// reset mother's resource value to zero
			SummerPop[i].resources = 0;
		
			// record mother's fecundity
			SummerPop[i].fecundity = resource_integer;
		
	        // for each parent create the number of offspring prescribed by their resource value + noise (i.e, resource_integer)
	        for (int kid_i = 0; kid_i < resource_integer; ++kid_i)
	        {
	            Individual kid;

	            create_offspring(mother, father, kid);

	            // add kid to the stack
	            Kids.push_back(kid);
	        }
		}  // Close ELSE loop
	} // end for (int i = 0; i < summer_pop; ++i)

	assert((breeder_pop + nonreproductive_pop) <= summer_pop);

	offspring_pop = Kids.size();

    // number of dead individuals is the max population
    // minus the current individuals in the summer population
    Nvacancies = N - summer_pop - winter_pop;

    assert(staging_pop == 0);
    assert(Nvacancies >= 0);
	assert(Nvacancies <= N);

    int random_kid_id = 0;

    // recruit new individuals to the summer pool
    for (int i = 0; i < Nvacancies; ++i)
    {
        // no kids left to recruit
        if (Kids.size() == 0)
        {
            break;
        }
    
        uniform_int_distribution<> kids_sample(0, Kids.size() - 1);
    
        random_kid_id = kids_sample(rng_r);

        // add random kid to population
        SummerPop[summer_pop] = Kids[random_kid_id];
		++summer_pop;

        //  delete randomly selected kid as it has been sampled
		// and replace with kid from end of stack
        Kids[random_kid_id] = Kids[Kids.size() - 1];
		
		// delete the now-replicated record
        Kids.pop_back();

    }  // Ends recruitment of offspring (kids)
	
	postbreeding_pop = summer_pop;
	assert(summer_pop + winter_pop <= N);
	
	// Randomly assign all individuals an initial habitat quality for postbreeding feeding
	
   // set lower boundary to the probability	
	for (int i = 0; i < summer_pop; ++i)
	    {			
			if (uniform(rng_r) < pgood)
			{
				SummerPop[i].patch_quality = 1;  // Will secure high quality foraging
			}
			else
			{
				SummerPop[i].patch_quality = 0;  // Will access low quality foraging
			}
	}
		
} // ENDS SUMMER REPRODUCTION


// gaining resources at breeding ground
// & fly back
void postbreeding_dynamics(int t)
{
    // foraging of individuals who are just at the breeding site
    // and who have yet to decide to go to the staging site
    for (int i = 0; i < summer_pop; ++i)
    {
		if (SummerPop[i].patch_quality == 1)
		{
			SummerPop[i].resources += rgood;
		}
		else
		{
			SummerPop[i].resources += rbad;
		}
		
		// possible change of foraging success in next timestep:	
		if (uniform(rng_r) < 1/pow(10, patch_consistency_factor))
        {
            SummerPop[i].patch_quality = 1 - SummerPop[i].patch_quality;  // switch of patch quality (good to poor; poor to good)
        }
		
		SummerPop[i].resources = min(SummerPop[i].resources, resource_max);
	
    } // ok resource dynamic done


    // foraging of individuals who are already at the staging site
    for (int i = 0; i < staging_pop; ++i)
    { 
        // indivuals who are already at the staging site
        // continue to forage at the staging site
		if (StagingPool[i].patch_quality == 1)
		{
			StagingPool[i].resources += rgood * preparation_penalty;
		}
		else
		{
			StagingPool[i].resources += rbad * preparation_penalty;
		}
		
		// possible change of foraging success in next timestep:	
		if (uniform(rng_r) < 1/pow(10, patch_consistency_factor))
        {
            StagingPool[i].patch_quality = 1 - StagingPool[i].patch_quality;  // switch of patch quality (good to poor; poor to good)
        }
    
	StagingPool[i].resources = min(StagingPool[i].resources, resource_max);
	
	}

    assert(summer_pop + winter_pop + staging_pop <= N);
    assert(summer_pop >= 0);

    double psignal = 0.0;

    // individuals decide whether to go to staging site
    // i.e., prepare for dispersal
    // signal to disperse
    for (int i = 0; i < summer_pop; ++i) 
    {
        // reaction norm dependent on resources
        // resulting in signaling a willingness to disperse
        // => go to the staging level
        //psignal = 0.5 * (SummerPop[i].theta_a[0] + SummerPop[i].theta_a[1])
        //    + 0.5 * (SummerPop[i].theta_b[0] + SummerPop[i].theta_b[1]) * SummerPop[i].resources; // resource_max;

		psignal = pow(1 + exp(-0.5 * (SummerPop[i].theta_b[0] + SummerPop[i].theta_b[1]) 
			* (SummerPop[i].resources - 0.5 * (SummerPop[i].theta_a[0] + SummerPop[i].theta_a[1]))), -1);
		
		// bound the probability
        psignal = clamp(psignal, 0.0, 1.0);

        // does individual want to signal to others to be ready for departure?
        if (uniform(rng_r) < psignal)
        {
            // add individual to the staging pool
            StagingPool[staging_pop] = SummerPop[i];
			StagingPool[staging_pop].latency = 0.0;
			StagingPool[staging_pop].signal_resources = StagingPool[staging_pop].resources;
            ++staging_pop; // increment the number of individuals in the staging pool

            // delete this individual from the summer population
			// and replace with the individual from the end of the summer_pop stack
            SummerPop[i] = SummerPop[summer_pop - 1];

            // decrement the number of individuals in the summer population
            --summer_pop;
            --i;
         
            // bounds checking   
            assert(staging_pop + winter_pop + summer_pop <= N);
            assert(staging_pop >= 0);
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

    // actual autumn dispersal at time t
    for (int i = 0; i < staging_pop; ++i)
    {
        // later we will consider collective dispersal decisions
        // for now, individuals leave dependent on the current amount of individuals
        // within the staging pool

        //pdisperse = 0.5 * (StagingPool[i].phi_a[0] + StagingPool[i].phi_a[1])
        //    + 0.5 * (StagingPool[i].phi_b[0] + StagingPool[i].phi_b[1]) * (double) staging_pop_start / (staging_pop_start + summer_pop);

		// SIGMOIDAL MODEL
		pdisperse = pow(1 + exp(-0.5 * (StagingPool[i].phi_b[0] + StagingPool[i].phi_b[1]) 
			* (((double) staging_pop_start / (staging_pop_start + summer_pop)) - 0.5 * (StagingPool[i].phi_a[0] + StagingPool[i].phi_a[1]))), -1);
		
		pdisperse = clamp(pdisperse, 0, 1);
		
		// yes individual goes
        if (uniform(rng_r) < pdisperse)
        {
			
			rv = StagingPool[i].resources;  // Record individual's resource value at departure
			mean_resources += rv;
			ss_resources += rv * rv;
			
			if (StagingPool[i].resources == resource_max)
			{
				++autumn_migrants_resource_cap;
			}
			
			WinterPop[winter_pop] = StagingPool[i];  // Individual moves from staging pool to first empty position in WinterPop
			WinterPop[winter_pop].signalling_proportion = staging_pop_start / (staging_pop_start + summer_pop);  // record the proportion of the population that was signalling at the point of departure
			lat = WinterPop[winter_pop].latency;
			mean_latency += lat;
			ss_latency += (lat * lat);			
			++winter_pop;
			
            // delete this individual from the staging population
			// replace with the last individual in the staging_pop stack
            StagingPool[i] = StagingPool[staging_pop - 1];

            // decrement the number of individuals in the staging population
            --staging_pop;
            --i;
			
			// increment the number of individuals recorded as autumn migrants
			++autumn_migrant_pop;
			
            assert(staging_pop + summer_pop + winter_pop <= N);
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
	
	population_mean_autumn_flock_size += NFlock;
	population_ss_autumn_flock_size += NFlock * NFlock;
	mean_autumn_staging_size += staging_pop_start;
	ss_autumn_staging_size += staging_pop_start * staging_pop_start;
	
    if (winter_pop_old < winter_pop){
		++n_autumn_flocks;
	}
	
	// update resource levels for all new individuals that have just
    // been added to the pool dependent on their flock size
    for (int i = winter_pop_old; i < winter_pop; ++i)
    {	
        // resources are reduced due to migration,
		WinterPop[i].cost = migration_cost(NFlock);
		WinterPop[i].resources -= migration_cost(NFlock);
		WinterPop[i].resources = min(WinterPop[i].resources, resource_max);  // Individual resource values cannot exceed resource max
		WinterPop[i].resources *= carryover_proportion;
		WinterPop[i].flock_size = NFlock;

    } // Ends: update resource levels of winter arrivals

} // ENDS: POST-BREEDING DYNAMICS 





// THE KEY PART OF THE CODE
// accepting command line arguments
int main(int argc, char **argv)
{
    string filename = "sim_migration";
    create_filename(filename);
    ofstream DataFile(filename.c_str());  // output file 

    // setting up filename for distribution file
    string dist_filename{filename + "_dist"};
    ofstream DistFile(dist_filename.c_str());

    init_arguments(argc, argv);

    write_parameters(DataFile);
	
	write_data_headers(DataFile);

    // write the data headers for the distribution file
    write_dist_data_headers(DistFile);

    init_population();

    for (int generation = 0; generation < number_generations; ++generation)
    {
        population_mean_spring_flock_size = 0.0;
		mean_spring_staging_size = 0.0;
		population_var_spring_flock_size = 0.0;
		var_spring_staging_size = 0.0;
		population_ss_spring_flock_size = 0.0;
		n_spring_flocks = 0;
		ss_spring_staging_size = 0.0;
		mean_spring_cost = 0.0;
		ss_spring_cost = 0.0;
		mean_signal_resources = 0.0;
		var_signal_resources = 0.0;
		ss_signal_resources = 0.0;
		mean_cost = 0.0;
		ss_cost = 0.0;
		staging_pop = 0;
		spring_pop_start = 0;
		autumn_pop_start = 0;
		summer_pop = 0;
		breeder_pop = 0;
		nonreproductive_pop = 0;
		offspring_pop = 0;
		summer_pop_old = 0;
		postbreeding_pop = 0;
		staging_pop = 0;  // Set staging population count to zero before winter dynamics
		remainer_pop = 0;
		
		rgood = rgood_init;
		rbad = rbad_init;  

		spring_nonmigrant_pop = 0;
		spring_migrants_resource_cap = 0;
		spring_migrant_pop = 0;
		spring_migrant_deaths = 0;
		autumn_nonmigrant_pop = 0;
		autumn_migrant_pop = 0;
		autumn_migrants_resource_cap = 0;
		
		
		// Allow populations to become established (individuals must acquire resources)
		if(generation < 0)  
			{	
				breeding_threshold = resource_reproduction_threshold * generation / 10000;
			}
		else
			{
				breeding_threshold = resource_reproduction_threshold;
			}

		// If individuals did not migrate to the breeding ground then the population size has potentially swollen above N. Random mortality returns population to N.
		// uniform_int_distribution <> pop_sample(0, winter_pop - 1);
		
		// while(winter_pop > N)
		//{
		//	WinterPop[pop_sample(rng_r)] = WinterPop[winter_pop - 1];
		//	--winter_pop;
			//} 
			
			// TODO Check that it is phenologically non-random to simply take the first N individuals as the population for year t+1
	
		assert(winter_pop <= N);
		
		spring_pop_start = winter_pop;
		
	    if (winter_pop <= 1)
	    { 
	        exit(1);
	    }
		
		// Winter foraging (migration is not an option)
		for (int t = 0; t < twinter; ++t)
		{
			winter_dynamics(t);
		}
		
		// time during spring during which individuals can migrate (or carry on foraging)
		mean_resources = 0.0;
		ss_resources = 0.0;
		rv = 0;
		
		for (int t = 0; t < tspring; ++t)
        {
            spring_dynamics(t);
			
			if (generation == number_generations - 1 && t >= tspring - 1)
				{
					write_dist(DistFile, generation, tspring);
				}
        }
				
		spring_migrant_pop = population_mean_spring_flock_size;
			
        // now take averages over all timesteps that individuals did (can) join groups
        population_mean_spring_flock_size = n_spring_flocks > 0 ? population_mean_spring_flock_size / n_spring_flocks : 0;
		mean_spring_staging_size /= tspring;
		
		// now record variance in flock size and staging size over the season
		population_var_spring_flock_size = n_spring_flocks > 0 ? (population_ss_spring_flock_size / n_spring_flocks) - (population_mean_spring_flock_size * population_mean_spring_flock_size) : 0;
		var_spring_staging_size = (ss_spring_staging_size / tspring) - (mean_spring_staging_size * mean_spring_staging_size);	
		
		clear_staging_pool();
		
		if (generation == 0)
		{
			write_spring_stats(DataFile, generation);
		}
		
		if ((generation + 1) % skip == 0)
		 {
			 write_spring_stats(DataFile, generation);
		  }  

        // let individuals die with a certain probability 
        spring_mortality();
		
	    if ((winter_pop + summer_pop) <= 1)
	    { 
	        exit(1);
	    }
		
		remainer_pop = winter_pop;
		
		// Individuals reproduce after they migrated to the summer spot
		Nvacancies = 0;
		rv = 0;
				
		mean_resources = 0.0;
		ss_resources = 0.0;
		rv = 0.0;
		
		if (summer_pop > 1)
		{
			summer_reproduction(DataFile);	
		}
		
		if (generation == 0)
		{
			write_summer_stats(DataFile);
		}
		
		if ((generation + 1) % skip == 0)
		{
			write_summer_stats(DataFile);
		}
		
		// set autumn migration stats to 0 before postbreeding_dynamics starts
        population_mean_autumn_flock_size = 0.0;
        mean_autumn_staging_size = 0.0;
		population_var_autumn_flock_size = 0.0;
		var_autumn_staging_size = 0.0;
		population_ss_autumn_flock_size = 0.0;
		n_autumn_flocks = 0.0;
		autumn_migrant_pop = 0.0;
		autumn_nonmigrant_pop = 0.0;
		ss_autumn_migrant_pop = 0.0;
		ss_autumn_staging_size = 0.0;
		mean_autumn_cost = 0.0;
		ss_autumn_cost = 0.0;
		mean_autumn_cost = 0.0;
		ss_autumn_cost = 0.0;
		mean_resources = 0.0;
		ss_resources = 0.0;
		rv = 0;
		autumn_migrant_deaths = 0;
		autumn_migrant_mortality_rate = 0.0;
		
        autumn_pop_start = summer_pop;
		
		// time during summer during which individuals forage
        for (int t = 0; t < tspring; ++t)
        {
            postbreeding_dynamics(t);
			
        }

		autumn_migrant_pop = population_mean_autumn_flock_size;
		
        // now take averages over all timesteps that individuals did (can) join groups
        population_mean_autumn_flock_size = n_autumn_flocks > 0 ?  population_mean_autumn_flock_size / n_autumn_flocks : 0;
        mean_autumn_staging_size /= tspring;
		
		// now record variance in autumn flock size and staging size over the season
		population_var_autumn_flock_size = n_autumn_flocks > 0 ? (population_ss_autumn_flock_size / n_autumn_flocks) - (population_mean_autumn_flock_size * population_mean_autumn_flock_size) : 0;
		var_autumn_staging_size = (ss_autumn_staging_size / tspring) - (mean_autumn_staging_size * mean_autumn_staging_size);
		
		autumn_nonmigrant_pop = summer_pop + staging_pop;
	  
		if (generation == 0)
		{
			write_autumn_stats(DataFile);
		}
		
		if ((generation+1) % skip == 0)
		{
			write_autumn_stats(DataFile);
		}
		
	    if (winter_pop <= 1)
	    { 
	        exit(1);
	    }
		
        // let individuals die with a certain probability 
        autumn_mortality();
		
		assert(winter_pop = remainer_pop + autumn_migrant_pop - autumn_migrant_deaths);
		
		if (generation == 0)
		{
			write_winter_stats(DataFile);
		}
		
		if ((generation+1) % skip == 0)
        {
            write_winter_stats(DataFile); 
        }
		
	    if (winter_pop <= 1)
	    { 
	        exit(1);
	    }
		
		// all individuals who remain at the summer grounds die
        summer_pop = 0;
        postbreeding_pop = 0;
		staging_pop = 0;
		breeder_pop = 0;
		nonreproductive_pop = 0;
		offspring_pop = 0;
		summer_pop_old = 0;
		spring_nonmigrant_pop = 0;
		spring_migrant_pop = 0;
		spring_migrants_resource_cap = 0;
		autumn_migrants_resource_cap = 0;
		autumn_nonmigrant_pop = 0;
		autumn_migrant_pop = 0;
				
    } // ENDS: GENERATION
}
