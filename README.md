# Running the simulation

First, clone this git repository to the computer/server on which you want to run the simulations.

Our simulation is written in C++. To compile it, open Terminal within the repository clone and enter `make clean` (not necessary if this is the first time you are compiling, but probably good practice to get into the habit). On the next command line, enter `make`. You should see confirmation that the C++ code has been compiled. If there is an issue, details should be provided, allowing you to make corrections and then try again (once again, `make clean` and then `make`).

Next, use `generate_parameters.py` to generate a list of individual simulations, defined by the combination of parameter values you have specified and the number of replicates for each unique combination.

e.g., `./generate_parameters.py > batch_list`

You can then use another Python file to generate a new folder with a title starting `hpcbatch_` followed by today's date and the current time. 

`./jobarray.py -i batch_list`

If you move to this directory (i.e., `cd ~/hpcbatch_`...) it should contain a document (filename starting `hpcjob_`; use, e.g., `ls -lh` to obtain a list of the folder contents).
Entering the command`sbatch hpcjob_`... will set these jobs running (or enter them in a queue). If you have more jobs to run than available cores on your computer then we reccomend using a cluster to run batches of simulations.

## Model details
The model is detailed in the paper reporting the results of our study of coordinated actions, including a description of its variables.
