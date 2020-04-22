#!/usr/bin/env python3

# generate all parameter combinations to run the migration simulation

init_theta_a = [0.005] # [0.001, 0.005, 0.01, 0.05, 0.1, 0.5]
init_theta_b = [0] #[ 0.0, 0.25, 0.5, 0.75, 1 ]
init_phi_a = [0.005] #[0.001, 0.005, 0.01, 0.05, 0.1, 0.5]
init_phi_b = [0] #[ 0.0, 0.25, 0.5, 0.75, 1 ]

tmax = 5000
twinter = 5000

pmort = [ 0.05 ]  # Short simulations on 27 November suggested pmort of 0.2 was too high. After further simulations, I decided to fix it at 0.05
pgood_init = [ 1.0 ] # Deleted 0.5 option on 27 November 2019
t_good_ends = [ 1500 ]

rgood = [ 0.01 ] # Changed from 1 on 27 November 2019
rbad = [ 0.005 ] # Changed from 0.5 on 27 November 2019

arrival_resource_decay = [0.1] # [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9] changed 19 Nov 2019. Outcome was that for value of 0.01, the winter population size was < 25. On 27th Nov, I fixed it to 0.1 as values of 0.2 were associated with very low population sizes.
resource_reproduction_threshold = [100, 200]  # [1]
resource_starvation_threshold = [0.0]
resource_max = 500

# mutation rates
mu_theta = 0.01
mu_phi = 0.01
sdmu_theta = 0.01
sdmu_phi = 0.01
    
# migration cost parameters
max_migration_cost = [25, 50, 75, 100] # [0.75, 0.625, 0.5, 0.375, 0.25, 0.2, 0.15, 0.125, 0.1, 0.075, 0.05, 0.025, 0.01]
min_migration_cost = [ 10 ]
migration_cost_decay = [ 1 ]
migration_cost_power = [ 0.5, 1, 2 ]

# reproductive cost parameters
min_offspring_cost = [ 5 ]
offspring_cost_magnifier = [ 3 ] # The relative difference in resource cost per offspring having migrated at the earliest opportunity versus the last


executable = "./xmigration"

counter = 1

# should jobs run on the background
# only do this if the number of jobs is smaller
# than the number of cores you have 
# do not do this on the carson cluster
background =  False

backgroundstr = ""

if background:
    backgroundstr = "&"


number_replicates = 5

for rep_i in range(0, number_replicates):
    for init_phi_a_i in init_phi_a:
        for init_phi_b_i in init_phi_b:
            for init_theta_a_i in init_phi_a:
                for init_theta_b_i in init_theta_b:
                    for pmort_i in pmort:
                        for pgood_init_i in pgood_init:
                            for t_good_ends_i in t_good_ends:
                                for rgood_i in rgood:
                                    for rbad_i in rbad:
                                        for arrival_resource_decay_i in arrival_resource_decay:
                                            for resource_reproduction_threshold_i in resource_reproduction_threshold:
                                                for resource_starvation_threshold_i in resource_starvation_threshold:
                                                    for max_migration_cost_i in max_migration_cost:
                                                        for min_migration_cost_i in min_migration_cost:
                                                            for migration_cost_decay_i in migration_cost_decay:
                                                                for migration_cost_power_i in migration_cost_power:
                                                                    for min_offspring_cost_i in min_offspring_cost:
                                                                        for offspring_cost_magnifier_i in offspring_cost_magnifier:

                                                                            # increment the counter for the number of 
                                                                            # runs
                                                                            counter += 1

            #                                                                print("echo " + str(counter))


                                                                            print(executable + " " 
                                                                                    + str(init_phi_a_i) + " "
                                                                                    + str(init_phi_b_i) + " "
                                                                                    + str(init_theta_a_i) + " "
                                                                                    + str(init_theta_b_i) + " "
                                                                                    + str(pmort_i) + " "
                                                                                    + str(pgood_init_i) + " "
                                                                                    + str(t_good_ends_i) + " "
                                                                                    + str(rgood_i) + " "
                                                                                    + str(rbad_i) + " "
                                                                                    + str(arrival_resource_decay_i) + " "
                                                                                    + str(resource_reproduction_threshold_i) + " "
                                                                                    + str(resource_starvation_threshold_i) + " "
                                                                                    + str(mu_theta) + " "
                                                                                    + str(mu_phi) + " "
                                                                                    + str(sdmu_theta) + " "
                                                                                    + str(sdmu_phi) + " "
                                                                                    + str(max_migration_cost_i) + " "
                                                                                    + str(min_migration_cost_i) + " "
                                                                                    + str(migration_cost_decay_i) + " "
                                                                                    + str(migration_cost_power_i) + " "
                                                                                    + str(tmax) + " " 
                                                                                    + str(twinter) + " "
                                                                                    + str(resource_max) + " "
                                                                                    + str(min_offspring_cost) + " "
                                                                                    + str(offspring_cost_magnifier) + " "
                                                                                    + backgroundstr)
