library("lattice")

# just to test our function
# relating flock size with 
# survival cost of migration

N <- seq(0,1000,1)

pow <- 1.2

max_migration_cost <- 1.0
migration_cost_decay <- 0
migration_cost_nonlinear_decay <- 0.001

power.shizzle <- max_migration_cost - migration_cost_decay * N - migration_cost_nonlinear_decay *   N^pow

print(
        xyplot(power.shizzle ~ N
                ,type="l"
                ,ylim=c(0,1)))
