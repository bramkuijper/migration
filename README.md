# The evolution of collective migration 

## Model variables 

### Reaction norm loci (TODO: reaction norm on what?, possible values)
- `theta_a` reaction norm on resources: elevation 
- `theta_b` reaction norm on resources: slope
- `phi_a` reaction norm on resources: elevation
- `phi_b` reaction norm on resources: slope

## Model variables 

### Life-history parameters
- `pmort` mortality probability (TODO: per what, for whom)

### Environmental parameters
- `pgood_init` probability at the start of each season that individuals encounter a good resource environment 
- `t_good_ends` time in each season after which it is no longer possible to find good resources
- `rgood` the value of good resources
- `rbad` the value of poor resources
- `resource_reproduce_threshold` the amount of resources necessary to reproduce

## The model

### Cost of migration
The cost of migration is calculated in terms of energy resources, where the fraction of resources `z` usurped by migration depends on flock size `n`, where

`z = max * (1.0 - decay * (n/N)^power`)

where `0 <= z <= 1`. Here, `max` reflects the maximum to the proportional cost, say if `max = 0.6` then only 60% of all resources will be removed by migrating alone. 

Next, `decay` is the strength with which costs decay with increasing flock size `n`, while `power` reflects whether costs rapidly decrease with even small group sizes (convex downward: `power > 1`), or whether costs start to decrease only with large group sizes (concave downward, `power < 1`).

Note that we take group size as a fraction `n/N` relative to the total population size `N`, as it ismathematically more convenient to scale relative to the maximum achievable flock size, as this implies that when `0 <= decay <= 1.0` we conveniently have `0 <= z <= 1`. Moreover, it also means we do not have to choose different values of `decay` and `power` when the expected group size `E[n]` (proportional to `N` changes. 
