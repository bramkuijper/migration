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

plt.style.use('base')

# some stuff to render fonts in graphs
rcParams['axes.labelsize'] = 15
rcParams['text.usetex'] = True
rcParams['font.family'] = 'sans-serif'

# some stuff to render fonts in graphs
# see http://stackoverflow.com/questions/2537868/sans-serif-math-with-latex-in-matplotlib 
rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]

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
line_num_params = find_parameter_linenum(sys.argv[1])

dat = pd.read_csv(sys.argv[1],
        nrows=line_num_params-1,
        sep=';',
        index_col=False)

#########################################

#           make figure

#########################################

# initialize the figure
fig = plt.figure(figsize=(10,15))

# generate the grid of the graph
# see: 
widths = [ 1 ]
heights = [ 1, 1, 1, 1, 1]
numrows = len(heights)
numcols  = len(widths)

# make the grid
gs = gridspec.GridSpec(
        nrows=numrows,
        ncols=numcols,
        width_ratios=widths,
        height_ratios=heights)


# start next entry of the graph
# the number of acts * their efficiency
ax = plt.subplot(gs[0,0])

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

ax = plt.subplot(gs[1,0])


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

ax = plt.subplot(gs[2,0])


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


ax = plt.subplot(gs[3,0])


ax.plot(
        dat["generation"],
        dat["mean_resources"],
        label=r"resources")

ax.set_ylabel(r"Resources")
ax.set_xlabel(r"Generation")

format = "pdf"

filename = os.path.join(
        os.path.dirname(sys.argv[1]),
        "graph_" + os.path.basename(sys.argv[1]) + "." + format
        )

plt.savefig(
        fname=filename,
        format=format, 
        bbox_inches="tight")
