library("lattice")

# a plot to assess how costs of 
# migration should decay with group size
# using a polynomial function

# different group sizes
group.size <- seq(0,200)

# the polynomial cost function
# evaluated against a row of a dataframe containing
# necessary parameters
cost <- function(row)
{
    total_cost <- 0
    
    with(as.list(row),{
        # evaluate the function
        total_cost <<- with(as.list(row),
                a - b * N - c * N^(powerN)
        )

        # if total cost is lower than minimum, 
        # give the minimum cost
        if (total_cost < min_cost)
        {
            total_cost <<- min_cost
        }
    })
    return(total_cost)
}

# grid with parameters for a 
# concave up decreasing function
# (i.e., where power of the nonlinear term above
# is <1.)
par.grid <- expand.grid(
        a=c(1.0)
        ,b=seq(0,0.01,0.002)
        ,c=seq(0,0.05,0.01)
        ,powerN=c(0.5)
        ,min_cost=c(0.05)
        ,N=group.size)

# calculate the cost using apply
par.grid["total_cost"] = apply(par.grid,FUN=cost,MARGIN=1)

# legend_labels
labels = as.character(sort(unique(par.grid$c)))
labels = sapply(X=labels,FUN=function(x) { return(paste("c = ",x,sep="",collapse=""))})

# now print a multivariate plot to see what happens
pdf("cost_functions_concave_up_decreasing.pdf",width=10,height=5)
print(xyplot(total_cost ~ N | b 
                ,data=par.grid
                ,groups=par.grid$c
                ,type="l"
                ,key=simpleKey(text=labels,lines=T,points=F)
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ))

dev.off()


# grid with parameters for concave down
# (i.e., where power of the nonlinear term above
# is >1.)
par.grid <- expand.grid(
        a=c(1.0)
        ,b=seq(0,0.01,0.002)
        ,c=seq(0,0.001,0.0002)
        ,powerN=c(2.0)
        ,min_cost=c(0.05)
        ,N=group.size)

par.grid["total_cost"] = apply(par.grid,FUN=cost,MARGIN=1)

# legend_labels
labels = as.character(sort(unique(par.grid$c)))
labels = sapply(X=labels,FUN=function(x) { return(paste("c = ",x,sep="",collapse=""))})

# now print a multivariate plot to see what happens
pdf("cost_functions_concave_down_decreasing.pdf",width=10,height=5)
print(xyplot(total_cost ~ N | b 
                ,data=par.grid
                ,groups=par.grid$c
                ,type="l"
                ,key=simpleKey(text=labels,lines=T,points=F)
                ,strip=function(strip.levels,...) { strip.default(strip.levels=T,...) }
                ))

dev.off()
