#!/usr/bin/env python3

# generate all parameter combinations to run the migration simulation

init_theta_a = [0.0005]
init_theta_b = [0]
init_phi_a = [0.0005]
init_phi_b = [0]

tmax = 5000
twinter = 5000

pmort = [0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22]
pgood_init = [ 0.5 ] # Deleted 0.5 option on 27 November 2019
t_good_ends = [ 5000 ]

rgood = [ 0.01 ]
rbad = [ 0.005 ]
preparation_penalty = [0.0] # The reduction in resource acquisition, relative to the normal feeding phase.

resource_reproduction_threshold = [50]
resource_starvation_threshold = [0.0]
resource_max = [100]

# mutation rates
mu_theta = 0.01
mu_phi = 0.01
sdmu_theta = [0.01]
sdmu_phi = sdmu_theta
    
# migration cost parameters
max_migration_cost = [20] 
min_migration_cost = [10]
migration_cost_power = [2]

# reproductive cost parameters
min_offspring_cost = [ 10 ]
offspring_cost_magnifier = [ 1.5 ] # The relative difference in resource cost per offspring having migrated at the earliest opportunity versus the last

number_replicates = 8

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
                                        for preparation_penalty_i in preparation_penalty:
                                            for resource_reproduction_threshold_i in resource_reproduction_threshold:
                                                for resource_starvation_threshold_i in resource_starvation_threshold:
                                                    for resource_max_i in resource_max:
                                                        for sdmu_theta_i in sdmu_theta:
                                                            for max_migration_cost_i in max_migration_cost:
                                                                for min_migration_cost_i in min_migration_cost:
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
                                                                                        + str(preparation_penalty_i) + " "
                                                                                        + str(resource_reproduction_threshold_i) + " "
                                                                                        + str(resource_starvation_threshold_i) + " "
                                                                                        + str(mu_theta) + " "
                                                                                        + str(mu_phi) + " "
                                                                                        + str(sdmu_theta_i) + " "
                                                                                        + str(sdmu_theta_i) + " "
                                                                                        + str(max_migration_cost_i) + " "
                                                                                        + str(min_migration_cost_i) + " "
                                                                                        + str(migration_cost_power_i) + " "
                                                                                        + str(tmax) + " " 
                                                                                        + str(twinter) + " "
                                                                                        + str(resource_max_i) + " "
                                                                                        + str(min_offspring_cost_i) + " "
                                                                                        + str(offspring_cost_magnifier_i) + " "
                                                                                        + backgroundstr)
