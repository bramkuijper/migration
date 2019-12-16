# The evolution of collective migration 

## Model variables 

### Reaction norm loci (TODO: reaction norm on what?, possible values)
- `theta_a` reaction norm of sigresources: elevation 
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

### Winter dynamics

#### Foraging
During each timestep `t`, individuals obtain with probability `pgood` a good resource (of value `rgood`) and with probability `1 - pgood` a bad resource (of value `rbad`). Here `pgood` is larger than 0 only for the first `t_good_ends` timesteps of the season, after which all individuals receive `rbad` resources. 

#### Migration decisions
There are several ways of modeling collective dispersal. One way is to have individuals join another individual, this other individual can either accept or reject. I am not too sure this applies here. Rather it may well be that individuals signal to one another that they are willing to migrate. The way we implement this is by having individuals first enter a staging pool, which indicates their willingness to migrate. Individuals can then decide to disperse dependent on the size of the staging pool (i.e., the potential migrating group)

##### Signaling willingness to migrate 
Individuals signal their willingness to disperse (i.e., enter the staging pool) dependent on resources. The probability to enter the staging pool is given by

`Prob(enter staging pool) = theta_a + theta_b * resources`

which is just a simple reaction norm. If the reaction norm is zero, individuals always signal their willingness to disperse. In the staging pool individuals still acquire resources each timestep. 

##### Signaling willingness to migrate 
In the staging pool, individuals can then decide to migrate dependent on the size of the staging pool, according to:

`Prob(migrate) = phi_a + phi_b * size_staging_pool`


### Cost of migration
The cost of migration is calculated in terms of energy resources, where the fraction of resources `z` usurped by migration depends on flock size `n`, where

`z = max * (1.0 - decay * (n/N)^power`)

where `0 <= z <= 1`. Here, `max` reflects the maximum to the proportional cost, say if `max = 0.6` then only 60% of all resources will be removed by migrating alone. 

Next, `decay` is the strength with which costs decay with increasing flock size `n`, while `power` reflects whether costs rapidly decrease with even small group sizes (convex downward: `power > 1`), or whether costs start to decrease only with large group sizes (concave downward, `power < 1`).

Note that we take group size as a fraction `n/N` relative to the total population size `N`, because it is mathematically more convenient to scale relative to the maximum achievable flock size, as this implies that when `0 <= decay <= 1.0` we conveniently have `0 <= z <= 1`. Most importantly, however, it means that the functional form of `z` is the same for populations that can have different maximum values of the migratory flock size `max(n)`, because `max(n) == N`. In other words: `z` has the same shape regardless of the population size considered.
