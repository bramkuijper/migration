#!/usr/bin/env python3

# plot the division of labor simulations

import pandas as pd
import itertools
import subprocess 
import math
import argparse
import numpy as np
import sys, re, os.path
import matplotlib
matplotlib.use('Agg')

from pathlib import Path
from pathlib import PurePath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import rcParams
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from matplotlib import cm


#########################################
#           check where data ends
#########################################
def find_parameter_linenum(filename):

    # get all the lines
    f = open(filename)
    fl = f.readlines()
    f.close()

    # find the first line from below that
    # starts with a number. That is the last
    # line of the data
    frange = list(range(0,len(fl)))
    frange.reverse()

    # ok loop over all lines (backwards)
    for line_num in frange:
        if re.match("^\d.*",fl[line_num]) is not None:
            return(line_num+1)

    return(len(fl))


#########################################
#           read in the data
#########################################

if len(sys.argv) < 2:
    print("plz provide a filename")
    sys.exit()

line_num_params = find_parameter_linenum(sys.argv[1])

dat = pd.read_csv(sys.argv[1],
        nrows=line_num_params-1,
        sep=';',
        index_col=False)

#########################################

#           make figure

#########################################

# initialize the figure
fig = plt.figure(figsize=(10,50))

# generate the grid of the graph
# see: 
widths = [ 1 ]
heights = [ 1 for x in range(10)]
numrows = len(heights)
numcols  = len(widths)

rowctr= 0

# make the grid
gs = gridspec.GridSpec(
        nrows=numrows,
        ncols=numcols,
        width_ratios=widths,
        height_ratios=heights)


# plot the resulting joining probability 
ax = plt.subplot(gs[rowctr,0])

# make a range of resources
resources = list(np.arange(0,100,1))

# now calculate the corresponding signal probability for each resource level
# based on the evolved values of the reaction norm for signaling
psignal = [ float(dat["mean_theta_a"][-1:]) + float(dat["mean_theta_b"][-1:]) * x for x in resources]

ax.plot(
        resources
        ,psignal
        ,label=r"reaction norm")

ax.set_ylabel(r"Willingness 2 disperse")
ax.set_xlabel(r"Resource level")
ax.set_ylim(0,1)

rowctr +=1

# plot the resulting migration probz
ax = plt.subplot(gs[rowctr,0])

# make a range of resources
group_size = list(np.arange(0,5000,1))


# now calculate the corresponding signal probability for each resource level
# based on the evolved values of the reaction norm for signaling
pdisperse = [ float(dat["mean_phi_a"][-1:]) + float(dat["mean_phi_b"][-1:]) * x for x in group_size]

ax.plot(
        group_size 
        ,pdisperse
        ,label=r"reaction norm")

ax.set_ylabel(r"Group size dispersal")
ax.set_xlabel(r"Group size")
ax.set_ylim(0,1)

rowctr +=1

# start next entry of the graph
# the number of acts * their efficiency
ax = plt.subplot(gs[rowctr,0])

ax.plot(
        dat["generation"],
        dat["mean_theta_a"],
        label=r"entry staging pool")

ax.plot(
        dat["generation"],
        dat["mean_phi_a"],
        label=r"migrate")

ax.legend()

ax.set_ylabel(r"Elevation, $a$")

ax.tick_params(
        axis="x",
        which="both",
        labelbottom=False)

rowctr +=1
ax = plt.subplot(gs[rowctr,0])


ax.plot(
        dat["generation"],
        dat["mean_theta_b"],
        label=r"entry staging pool")

ax.plot(
        dat["generation"],
        dat["mean_phi_b"],
        label=r"migrate")

ax.legend()

ax.set_ylabel(r"Slope, $b$")

rowctr +=1

ax = plt.subplot(gs[rowctr,0])


ax.plot(
        dat["generation"],
        dat["nwinter"],
        label=r"$N_{\mathrm{winter}}$")

ax.plot(
        dat["generation"],
        dat["nsummer"],
        label=r"$N_{\mathrm{summer}}$")

ax.plot(
        dat["generation"],
        dat["nstaging"],
        label=r"$N_{\mathrm{staging}}$")

ax.plot(
        dat["generation"],
        dat["nkids"],
        label=r"$N_{\mathrm{kids}}$")

ax.legend()

ax.set_ylabel(r"Popsize")


rowctr +=1

ax = plt.subplot(gs[rowctr,0])


ax.plot(
        dat["generation"],
        dat["mean_resources"],
        label=r"resources")

ax.set_ylabel(r"Resources")
ax.set_xlabel(r"Generation")

rowctr +=1

ax = plt.subplot(gs[rowctr,0])


ax.plot(
        dat["generation"],
        dat["mean_staging_size_summer"],
        label=r"$N_{\mathrm{staging,summer}}i$")

ax.plot(
        dat["generation"],
        dat["mean_staging_size_winter"],
        label=r"$N_{\mathrm{staging,winter}}i$")

ax.plot(
        dat["generation"],
        dat["mean_flock_size_winter"],
        label=r"$N_{\mathrm{flock,winter}}i$")

ax.plot(
        dat["generation"],
        dat["mean_flock_size_summer"],
        label=r"$N_{\mathrm{flock,summer}}i$")

ax.set_ylabel(r"Flock, Staging size")
ax.set_xlabel(r"Generation")
ax.legend()

format = "pdf"

filename = os.path.join(
        os.path.dirname(sys.argv[1]),
        "graph_" + os.path.basename(sys.argv[1]) + "." + format
        )

plt.savefig(
        fname=filename,
        format=format, 
        bbox_inches="tight")
