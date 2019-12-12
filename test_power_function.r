library("lattice")

# just to test our function
# relating flock size with 
# survival cost of migration

parameters <- as.data.frame(expand.grid(
        rel_pop=seq(0,1,0.01)
        ,slope=c(2.0,1.0,0.5,0.25)
        ,power=c(1.0,2.0,0.5)))

max_migration_cost <- 1.0

calc_mort_prob <- function(x) {
    val <- with(data=as.list(x) 
            ,expr= 1.0 
                    - slope * (rel_pop)^power)

    if (val > 1.0)
    {
        val <- 1.0
    } else if (val < 0.0)
    {
        val <- 0.0
    }

    return(val)
}

# calculate mortality probability for all parameter values
parameters[,"mortality_prob"] <- apply(
        X=parameters
        ,MARGIN=1
        ,FUN=calc_mort_prob)

# sort on population size to draw straightish lines
# rather than lots of tangled up lines
parameters <- parameters[order(parameters$rel_pop),]

pdf("example_lines.pdf")
print(
        xyplot(mortality_prob ~ rel_pop | power
                ,type="l"
                ,data=parameters
                ,groups=slope
                ,auto.key=T
                ,xlab="Ngroup/Ngroupmax"
                ,ylab="mortality prob"
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ,ylim=c(-0.1,1.1)))

dev.off()


