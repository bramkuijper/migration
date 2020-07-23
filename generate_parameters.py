#!/usr/bin/env python3

# generate all parameter combinations to run the migration simulation

init_theta_a = [0.0005]  # Default is 0.0005
init_theta_b = [0]  # Default is 0
init_phi_a = [0.0005]  # Default is 0.0005
init_phi_b = [0]  # Default is 0

twinter = 0
tspring = 10000

pmort = [0.1] #[0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25]  # Default is 0.05
pgood_init = 0.5 # Deleted 0.5 option on 27 November 2019
t_good_ends = 500000

rgood = [ 0.005, 0.01, 0.015, 0.02, 0.025, 0.03 ]  # Default 0.01
rbad = [ 0.0075 ]  # Default 0.005
preparation_penalty = [0] # The reduction in resource acquisition, relative to the normal feeding phase.

resource_reproduction_threshold = [30]
resource_starvation_threshold = 0
resource_max = [100]

# mutation rates
mu_theta = 0.01
mu_phi = mu_theta
sdmu_theta = 0.002
sdmu_phi = sdmu_theta
    
# migration cost parameters
max_migration_cost = 20  # Default is 20 
min_migration_cost = [20] #[20, 18, 16, 14, 12, 10, 8, 6, 4]
migration_cost_power = [2]

# reproductive cost parameters
min_offspring_cost = [ 5 ]
offspring_cost_magnifier = [ 1 ] # The relative difference in resource cost per offspring having migrated at the earliest opportunity versus the last
relative_mortality_risk_of_migration = [5]

carryover_proportion = [1]

number_replicates = 5

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
                        for rgood_i in rgood:
                            for rbad_i in rbad:
                                for preparation_penalty_i in preparation_penalty:
                                    for resource_reproduction_threshold_i in resource_reproduction_threshold:
                                        for resource_max_i in resource_max:
                                            for sdmu_theta_i in sdmu_theta:
                                                for min_migration_cost_i in min_migration_cost:
                                                    for migration_cost_power_i in migration_cost_power:
                                                        for min_offspring_cost_i in min_offspring_cost:
                                                            for offspring_cost_magnifier_i in offspring_cost_magnifier:
                                                                for carryover_proportion_i in carryover_proportion:
                                                                    for relative_mortality_risk_of_migration_i in relative_mortality_risk_of_migration:

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
                                                                                + str(pgood_init) + " "
                                                                                + str(t_good_ends) + " "
                                                                                + str(rgood_i) + " "
                                                                                + str(rgood_i) + " "
                                                                                + str(preparation_penalty_i) + " "
                                                                                + str(resource_reproduction_threshold_i) + " "
                                                                                + str(resource_starvation_threshold) + " "
                                                                                + str(mu_theta) + " "
                                                                                + str(mu_theta) + " "
                                                                                + str(sdmu_theta_i) + " "
                                                                                + str(sdmu_theta_i) + " "
                                                                                + str(max_migration_cost) + " "
                                                                                + str(min_migration_cost_i) + " "
                                                                                + str(migration_cost_power_i) + " "
                                                                                + str(twinter) + " "
                                                                                + str(tspring) + " " 
                                                                                + str(resource_max_i) + " "
                                                                                + str(min_offspring_cost_i) + " "
                                                                                + str(offspring_cost_magnifier_i) + " "
                                                                                + str(carryover_proportion_i) + " "
                                                                                + str(relative_mortality_risk_of_migration_i) + " "
                                                                                + backgroundstr)
