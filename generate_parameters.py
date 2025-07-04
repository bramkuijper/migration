#!/usr/bin/env python3
import datetime
# generate all parameter combinations to run the migration simulation

init_theta_a = [60]  # Default is 60
init_theta_b = [1]  # Default is 1
init_phi_a = [0.3]  # Default is 0.3
init_phi_b = [25]  # Default is 25
init_psi_a = 1
init_psi_b = 0

twinter = 0
tspring = 10000  # twinter in manuscript (default is 10K)

pmort = 0.1
pgood = 0.5
patch_consistency_factor = [0.30103] #[0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3] # Varies from 0 (switching patch type EVERY timestep, so uniformity in resource value) and upwards. Above 3 (equating to a switch probability of 0.001) the distribution of resource values starts to become bimodal. log(2) gives the previous behaviour, which was a 50:50 likelihood of switching.

rgood = 0.04
rbad = 0.02
preparation_penalty = [0.5] # [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0] # The reduction in resource acquisition, relative to the normal feeding phase. Default is 0.5

resource_reproduction_threshold = 30
resource_starvation_threshold = 0
resource_max = 60

# mutation rates
mu = 0.05
sdmu = 0.2
    
# migration cost parameters
max_migration_cost = 20  # Default is 20 
min_migration_cost = [20] # [20, 16, 12, 8, 4]  # Default is 20
cost_power = [2]  # Default is 2
capacity = [8]
# reproductive cost parameters
min_offspring_cost = 5
offspring_cost_magnifier = [1] # The relative difference in resource cost per offspring having migrated at the earliest opportunity versus the last. 1 represents seasonal invariability
relative_mortality_risk_of_migration = 5
socially_sensitive_mortality = [0, 0.2, 0.4, 0.6, 0.8, 1] # [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]  # Default would be 0, where the mortality rate is indepdendent of flock size
postequilibrialisation_experimental_runtime = 0
K_decline_factor = [1] # 1 represents no decline in carrying capacity
autumn_harvest = [0.75] # Proportion of the population to be harvested: 0 represents no harvest
carryover_proportion = 0 # carryover of resources from one year to the next

#costs_sourcefile = ["/Users/simonevans/Library/CloudStorage/OneDrive-UniversityofExeter/Research/Modelling\ migration/hpcbatch_30_08_2021_174744/sim_migration_14_9_2021_133518_1302909500_dist"] #["/nobackup/beegfs/home/ISAD/sre206/hpcbatch_29_08_2021_172728/sim_migration_30_8_2021_075422_139529932_dist", "/nobackup/beegfs/home/ISAD/sre206/hpcbatch_29_08_2021_172728/sim_migration_30_8_2021_075227_1496541961_dist", "/nobackup/beegfs/home/ISAD/sre206/hpcbatch_29_08_2021_172728/sim_migration_30_8_2021_075428_1418409188_dist", "/nobackup/beegfs/home/ISAD/sre206/hpcbatch_29_08_2021_172728/sim_migration_30_8_2021_075227_1191646591_dist", "/nobackup/beegfs/home/ISAD/sre206/hpcbatch_29_08_2021_172728/sim_migration_30_8_2021_075227_1140373795_dist", "/nobackup/beegfs/home/ISAD/sre206/hpcbatch_29_08_2021_172948/sim_migration_30_8_2021_075428_1605237909_dist", "/nobackup/beegfs/home/ISAD/sre206/hpcbatch_29_08_2021_172948/sim_migration_30_8_2021_075432_1732123732_dist", "/nobackup/beegfs/home/ISAD/sre206/hpcbatch_29_08_2021_172948/sim_migration_30_8_2021_075431_665834704_dist", "/nobackup/beegfs/home/ISAD/sre206/hpcbatch_29_08_2021_172948/sim_migration_30_8_2021_075434_333892663_dist", "/nobackup/beegfs/home/ISAD/sre206/hpcbatch_29_08_2021_172948/sim_migration_30_8_2021_075429_1503304252_dist", "/nobackup/beegfs/home/ISAD/sre206/hpcbatch_29_08_2021_172948/sim_migration_30_8_2021_075431_660393031_dist", "/nobackup/beegfs/home/ISAD/sre206/hpcbatch_29_08_2021_172948/sim_migration_30_8_2021_075429_1499475050_dist", "/nobackup/beegfs/home/ISAD/sre206/hpcbatch_29_08_2021_172948/sim_migration_30_8_2021_075428_1777764160_dist"]
costs_sourcefile = ["none"]  # If none, enter "none"

risks_sourcefile = ["none"]  # If none, enter "none"

#risks_sourcefile = ["/Users/simonevans/Library/CloudStorage/OneDrive-UniversityofExeter/Research/Modelling\ migration/hpcbatch_29_08_2021_172728/sim_migration_30_8_2021_075422_139529932_dist"] 

## Control simulations for main results:
#risks_sourcefile = ["/nobackup/beegfs/home/ISAD/sre206/hpcbatch_22_11_2023_104936/sim_migration_22_11_2023_104929_11_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_22_11_2023_104936/sim_migration_22_11_2023_104929_22_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_22_11_2023_104936/sim_migration_22_11_2023_104929_33_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_22_11_2023_104936/sim_migration_22_11_2023_104929_44_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_22_11_2023_104936/sim_migration_22_11_2023_104929_55_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_22_11_2023_104936/sim_migration_22_11_2023_104929_66_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_22_11_2023_104936/sim_migration_22_11_2023_104929_77_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_22_11_2023_104936/sim_migration_22_11_2023_104929_88_dist"]

## Control simulations for autumnal cull (taking distributions from the 25 simulations with 100% social sensitivity)
## Remember to set number_replicates = 1
#risks_sourcefile = ["/nobackup/beegfs/home/ISAD/sre206/hpcbatch_24_11_2023_104810/sim_migration_24_11_2023_104801_42_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_24_11_2023_104810/sim_migration_24_11_2023_104801_10_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_24_11_2023_104810/sim_migration_24_11_2023_104801_16_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_24_11_2023_104810/sim_migration_24_11_2023_104801_12_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_24_11_2023_104810/sim_migration_24_11_2023_104801_48_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_24_11_2023_104810/sim_migration_24_11_2023_104801_44_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_24_11_2023_104810/sim_migration_24_11_2023_104801_30_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_24_11_2023_104810/sim_migration_24_11_2023_104801_40_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_24_11_2023_104810/sim_migration_24_11_2023_104801_26_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_24_11_2023_104810/sim_migration_24_11_2023_104801_32_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_24_11_2023_104810/sim_migration_24_11_2023_104801_34_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_24_11_2023_104810/sim_migration_24_11_2023_104801_46_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_24_11_2023_104810/sim_migration_24_11_2023_104801_50_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_24_11_2023_104810/sim_migration_24_11_2023_104801_18_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_24_11_2023_104810/sim_migration_24_11_2023_104801_2_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_24_11_2023_104810/sim_migration_24_11_2023_104801_20_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_24_11_2023_104810/sim_migration_24_11_2023_104801_24_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_24_11_2023_104810/sim_migration_24_11_2023_104801_28_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_24_11_2023_104810/sim_migration_24_11_2023_104801_36_dist",
#"/nobackup/beegfs/home/ISAD/sre206/hpcbatch_24_11_2023_104810/sim_migration_24_11_2023_104801_4_dist"]

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


date = datetime.datetime.now()
base_name = "sim_migration_" +\
        f"{date:%d}_{date:%m}_{date:%Y}_{date:%H}{date:%M}{date:%S}"


for rep_i in range(0, number_replicates):
    for patch_consistency_factor_i in patch_consistency_factor:
        for preparation_penalty_i in preparation_penalty:
            for min_migration_cost_i in min_migration_cost:
                for cost_power_i in cost_power:
                   for capacity_i in capacity:
                        for offspring_cost_magnifier_i in offspring_cost_magnifier:
                            for socially_sensitive_mortality_i in socially_sensitive_mortality:
                                for K_decline_factor_i in K_decline_factor:
                                    for autumn_harvest_i in autumn_harvest:
                                        for costs_sourcefile_i in costs_sourcefile:
                                            for risks_sourcefile_i in risks_sourcefile:
                                                for init_phi_a_i in init_phi_a:
                                                    for init_phi_b_i in init_phi_b:
                                                        for init_theta_a_i in init_theta_a:
                                                            for init_theta_b_i in init_theta_b:

                # print("echo " + str(counter))
                                                                print(executable + " " 
                                                                    + str(init_theta_a_i) + " "
                                                                    + str(init_theta_b_i) + " "
                                                                    + str(init_phi_a_i) + " "
                                                                    + str(init_phi_b_i) + " "
                                                                    + str(init_psi_a) + " "
                                                                    + str(init_psi_b) + " "  #6
                                                                    + str(pmort) + " "
                                                                    + str(pgood) + " "
                                                                    + str(patch_consistency_factor_i) + " "
                                                                    + str(rgood) + " "
                                                                    + str(rbad) + " "
                                                                    + str(preparation_penalty_i) + " "         #12
                                                                    + str(resource_reproduction_threshold) + " "
                                                                    + str(resource_starvation_threshold) + " "
                                                                    + str(mu) + " "
                                                                    + str(sdmu) + " "                          #16
                                                                    + str(max_migration_cost) + " "
                                                                    + str(min_migration_cost_i) + " "
                                                                    + str(cost_power_i) + " "
                                                                    + str(twinter) + " "
                                                                    + str(tspring) + " " 
                                                                    + str(resource_max) + " "
                                                                    + str(min_offspring_cost) + " "
                                                                    + str(offspring_cost_magnifier_i) + " "    #24
                                                                    + str(carryover_proportion) + " "
                                                                    + str(relative_mortality_risk_of_migration) + " "
                                                                    + str(socially_sensitive_mortality_i) + " "
                                                                    + str(capacity_i) + " "
                                                                    + str(postequilibrialisation_experimental_runtime) + " "
                                                                    + str(K_decline_factor_i) + " "             #30
                                                                    + str(autumn_harvest_i) + " "
                                                                    + str(costs_sourcefile_i) + " "
                                                                    + str(risks_sourcefile_i) + " "
                                                                    + str(base_name) + "_" + str(counter) + " "
                                                                    + backgroundstr)
                                                                # increment the counter for the number of 
                                                                # runs
                                                                counter += 1
