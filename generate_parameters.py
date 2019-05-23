#!/usr/bin/env python3

# generate all parameter combinations to run the migration simulation

init_theta_a = 0.05
init_theta_b = 0.0
init_phi_a = 0.05
init_phi_b = 0.0

tmax = 20

pmort = [ 0.05, 0.1, 0.2 ]
pgood_init = [ 1.0, 0.5 ]
decay_good = [ 1.0/tmax ]

rgood = [ 1 ]
rbad = [ 0.5 ]

arrival_resource_decay = [0.1, 0.5]
resource_reproduction_threshold = [ 1 ]

mu_theta = 0.01
mu_phi = 0.01
sdmu_theta = 0.01
sdmu_phi = 0.01

executable = "./xmigration"

counter = 0

for pmort_i in pmort:
    for pgood_init_i in pgood_init:
        for decay_good_i in decay_good:
            for rgood_i in rgood:
                for rbad_i in rbad:
                    for arrival_resource_decay_i in arrival_resource_decay:
                        for resource_reproduction_threshold_i in resource_reproduction_threshold:

                            # increment the counter for the number of 
                            # runs
                            counter += 1

                            print("echo " + str(counter))


                            print(executable + " " 
                                    + str(init_phi_a) + " "
                                    + str(init_phi_b) + " "
                                    + str(init_theta_a) + " "
                                    + str(init_theta_b) + " "
                                    + str(pmort_i) + " "
                                    + str(pgood_init_i) + " "
                                    + str(decay_good_i) + " "
                                    + str(rgood_i) + " "
                                    + str(rbad_i) + " "
                                    + str(arrival_resource_decay_i) + " "
                                    + str(resource_reproduction_threshold_i) + " "
                                    + str(mu_theta) + " "
                                    + str(mu_phi) + " "
                                    + str(sdmu_theta) + " "
                                    + str(sdmu_phi) + " "
                                    + str(tmax))




