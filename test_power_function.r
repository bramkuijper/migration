library("lattice")

# just to test our function
# relating flock size with 
# survival cost of migration

N <- seq(0,5000,10)

pow <- 3.2

max_migration_cost <- 1.0
migration_cost_decay <- 0
migration_cost_nonlinear_decay <- 0.0000000000012

power.shizzle <- max_migration_cost - migration_cost_decay * N - migration_cost_nonlinear_decay * N^pow

print(
        xyplot(power.shizzle ~ N
                ,type="l"
                ,ylim=c(0,1)))
                

pow <- 1

max_migration_cost <- 1.0
migration_cost_decay <- 0.0002
migration_cost_nonlinear_decay <- 0

power.shizzle <- max_migration_cost - migration_cost_decay * N - migration_cost_nonlinear_decay * N^pow

print(
        xyplot(power.shizzle ~ N
                ,type="l"
                ,ylim=c(0,1)))                  
                

pow <- 0.62

max_migration_cost <- 1.0
migration_cost_decay <- -0.0003
migration_cost_nonlinear_decay <- 0.012

power.shizzle <- max_migration_cost - migration_cost_decay * N - migration_cost_nonlinear_decay * N^pow

print(
        xyplot(power.shizzle ~ N
                ,type="l"
                ,ylim=c(0,1)))                
