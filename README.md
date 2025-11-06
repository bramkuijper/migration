# The evolution of collective migration 

## Running the simulation

First, clone this git repository to the computer/server on which you want to run the simulations.

Our simulation is written in C++. To compile it, open Terminal within the repository clone and enter `make clean` (not necessary if this is the first time you are compiling, but probably good practice to get into the habit). On the next command line, enter `make`. You should see confirmation that the C++ code has been compiled. If there is an issue, details should be provided, allowing you to make corrections and then try again (once again, `make clean` and then `make`).

Next, use `generate_parameters.py` to generate a list of individual simulations, defined by the combination of parameter values you have specified and the number of replicates for each unique combination.

e.g., `./generate_parameters.py > batch_list`

You can then use another Python file to generate a new folder with a title starting `hpcbatch_` followed by today's date and the current time. 

`./jobarray.py -i batch_list`

If you move to this directory (i.e., `cd ~/hpcbatch_`...) it should contain a document (filename starting `hpcjob_`; use, e.g., `ls -lh` to obtain a list of the folder contents).
Entering the command`sbatch hpcjob_`... will set these jobs running (or enter them in a queue). If you have more jobs to run than available cores on your computer then we reccomend using a cluster to run batches of simulations.

## Model variables 

### Decision rule loci
- `theta_a` condition-dependent signalling: sigmoidal midpoint
- `theta_b` condition-dependent signalling: logarithmic growth rate (i.e., sensitivity)
- `phi_a` socially-dependent departure: sigmoidal midpoint
- `phi_b` socially-dependent departure: logarithmic growth rate (i.e., sensitivity)

There are an additional two loci facilitating condition-dependent departure that are, by default, fixed to zero
– `psi_a` condition-dependent departure: sigmoidal midpoint
- `psi_b` condition-dependent departure: logarithmic growth rate (i.e., sensitivity)

## Model variables 

### Life-history parameters
- `pmort` mortality probability (TODO: per what, for whom)

### Environmental parameters
- `pgood_init` probability at the start of each season that individuals encounter a good resource environment 
- `t_good_ends` time in each season after which it is no longer possible to find good resources
- `rgood` the value of good resources
- `rbad` the value of poor resources
- `resource_reproduce_threshold` the minimum amount of resources necessary to reproduce

## The model

### Winter dynamics

#### Foraging
During each timestep `t`, individuals obtain with probability `pgood` a good resource (of value `rgood`) and with probability `1 - pgood` a bad resource (of value `rbad`). Here `pgood` is larger than 0 only for the first `t_good_ends` timesteps of the season, after which all individuals receive `rbad` resources. 

#### Migration decisions
There are several ways of modelling collective dispersal. One way is to have individuals join another individual, this other individual can either accept or reject. I am not too sure this applies here. Rather it may well be that individuals signal to one another that they are willing to migrate. The way we implement this is by having individuals first enter a staging pool, which indicates their willingness to migrate. Individuals can then decide to disperse dependent on the size of the staging pool (i.e., the potential migrating group)

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
