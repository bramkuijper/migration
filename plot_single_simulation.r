#!/usr/bin/env Rscript 

#--vanilla

library("ggplot2")

# get command line arguments
args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1)
{
    print("provide a simulation file name")
    stop()
}

find_out_param_line <- function(filename) {

    f <- readLines(filename)

    # make a reverse sequence
    seqq <- seq(length(f),1,-1)

    # go through each line in the data file and find first line
    # where data is printed (i.e., a line which starts with a digit)
    for (line_i in seqq)
    {
        if (length(grep("^\\d",f[[line_i]])) > 0)
        {
            return(line_i)
        }
    }

    return(NA)
}

parameter_row <- find_out_param_line(args[1])

if (is.na(parameter_row))
{
    print("cannot find data...")
    stop()
}

# read in data frame of corresponding simulation
the.data <- read.table(args[1], header=T, nrow=parameter_row)


# now plot stuff
