#!/usr/bin/env python3

# generate all parameter combinations to run the migration simulation

init_theta_a = [60]
init_theta_b = [1]
init_phi_a = [0.3]
init_phi_b = [25]

twinter = 0
tspring = 10000  # twinter in manuscript

pmort = [0.1]
pgood = 0.5
patch_consistency_factor = [0.5, 1, 1.5, 2, 2.5, 3]  # Varies from 0 (switching patch type EVERY timestep, so uniformity in resource value) and upwards. Above 3 (equating to a switch probability of 0.001) the distribution of resource values starts to become bimodal.

rgood = [ 0.04 ]
rbad = [ 0.02 ]
preparation_penalty = [0.5] # The reduction in resource acquisition, relative to the normal feeding phase.

resource_reproduction_threshold = [30]
resource_starvation_threshold = 0
resource_max = [60]

# mutation rates
mu_theta = 0.01
mu_phi = 0.01
sdmu_theta = 0.05
sdmu_phi = [0.05]
    
# migration cost parameters
max_migration_cost = 20  # Default is 20 
min_migration_cost = [10] # [20, 18, 16, 14, 12, 10, 8, 6, 4, 2]
migration_cost_power = [2] #[3, 2, 1, 0, -1, -2, -3]
capacity = 8

# reproductive cost parameters
min_offspring_cost = [ 5 ]
offspring_cost_magnifier = [ 1 ] # The relative difference in resource cost per offspring having migrated at the earliest opportunity versus the last
relative_mortality_risk_of_migration = [5]

carryover_proportion = [0]

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
            for init_theta_a_i in init_theta_a:
                for init_theta_b_i in init_theta_b:
                    for pmort_i in pmort:
                        for patch_consistency_factor_i in patch_consistency_factor:
                            for rgood_i in rgood:
                                for rbad_i in rbad:
                                    for preparation_penalty_i in preparation_penalty:
                                        for resource_reproduction_threshold_i in resource_reproduction_threshold:
                                            for resource_max_i in resource_max:
                                                for sdmu_phi_i in sdmu_phi:
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
                                                                                    + str(pmort_i) + " " #5
                                                                                    + str(pgood) + " "
                                                                                    + str(patch_consistency_factor_i) + " "
                                                                                    + str(rgood_i) + " "
                                                                                    + str(rbad_i) + " "  #9
                                                                                    + str(preparation_penalty_i) + " "  #10
                                                                                    + str(resource_reproduction_threshold_i) + " "
                                                                                    + str(resource_starvation_threshold) + " "
                                                                                    + str(mu_theta) + " "
                                                                                    + str(mu_phi) + " "  #14
                                                                                    + str(sdmu_theta) + " "  #15
                                                                                    + str(sdmu_phi_i) + " "  #16
                                                                                    + str(max_migration_cost) + " "
                                                                                    + str(min_migration_cost_i) + " "
                                                                                    + str(migration_cost_power_i) + " "  #19
                                                                                    + str(twinter) + " "  #20
                                                                                    + str(tspring) + " " 
                                                                                    + str(resource_max_i) + " "
                                                                                    + str(min_offspring_cost_i) + " "
                                                                                    + str(offspring_cost_magnifier_i) + " "  #24
                                                                                    + str(carryover_proportion_i) + " "  #25
                                                                                    + str(relative_mortality_risk_of_migration_i) + " "
                                                                                    + str(capacity) + " " #27
                                                                                    + backgroundstr)
