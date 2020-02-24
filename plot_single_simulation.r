#!/usr/bin/env Rscript 

#--vanilla

library("ggplot2")
library("gridExtra")

# get command line arguments
args = commandArgs(trailingOnly=TRUE)

# give an error message if you do not provide it with a simulation file name
if (length(args) < 1)
{
    print("provide a simulation file name")
    stop()
}

# find out where the parameter listing starts
# so that we can read in the data part of the file 
# without having it messed up by the subsequent parameter listing
find_out_param_line <- function(filename) {

    f <- readLines(filename)

    # make a reverse sequence
    seqq <- seq(length(f),1,-1)

    # go through each line in the data file and find first line
    # where data is printed (i.e., a line which starts with a digit)
    for (line_i in seqq)
    {
        print(f[[line_i]])
        print(line_i)
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
the.data <- read.table(args[1], header=T, nrow=parameter_row - 1, sep=";")

# now use ggplot2 to plot stuff

str(the.data)

p1.a <- ggplot(data=the.data
        ,aes(x=generation)) +
            geom_line(aes(y = winter_pop, colour="Winter")) +
            theme_classic() + 
            xlab("Generation") + 
            ylab("Population size")

p1.b <- ggplot(data=the.data
        ,aes(x=generation)) +
            geom_line(aes(y = summer_pop, colour="Summer")) +
            theme_classic() + 
            xlab("Generation") + 
            ylab("Population size")

p2 <- ggplot(data=the.data
        ,aes(x=generation)) +
            geom_line(aes(y = mean_spring_staging_size, colour="Spring")) +
            geom_line(aes(y = mean_autumn_staging_size, colour = "Autumn")) + 
            theme_classic() + 
            xlab("Generation") + 
            ylab("Mean staging size")

p3 <- ggplot(data=the.data
        ,aes(x=generation)) +
            geom_line(aes(y = mean_spring_flock_size, colour="Spring")) +
            geom_line(aes(y = mean_autumn_flock_size, colour = "Autumn")) + 
            theme_classic() + 
            xlab("Generation") + 
            #ylim(c(0,10)) +
            ylab("Mean flock size")

p3b <- ggplot(data=the.data
        ,aes(x=generation)) +
            geom_line(aes(y = var_spring_flock_size, colour="Spring")) +
            geom_line(aes(y = var_autumn_flock_size, colour = "Autumn")) + 
            theme_classic() + 
            xlab("Generation") + 
            #ylim(c(0,10)) +
            ylab("Var flock size")

p3c <- ggplot(data=the.data
        ,aes(x=generation)) +
            geom_line(aes(y = n_spring_flocks, colour="Spring")) +
            geom_line(aes(y = n_autumn_flocks, colour = "Autumn")) + 
            theme_classic() + 
            xlab("Generation") + 
            #ylim(c(0,10)) +
            ylab("Number flocks")

p3d <- ggplot(data=the.data
		,aes(x=generation)) +
			geom_line(aes(y = mean_spring_latency, colour="Spring")) +
			geom_line(aes(y = mean_autumn_latency, colour="Autumn")) +
			theme_classic() +
			xlab("Generation") +
			ylab("Mean latency")
			
p3e <- ggplot(data=the.data
		,aes(x=generation)) +
			geom_line(aes(y = mean_spring_departure, colour="Spring")) +
			geom_line(aes(y = mean_autumn_departure, colour="Autumn")) +
			theme_classic() +
			xlab("Generation") +
			ylab("Mean timing")

p3f <- ggplot(data=the.data
        ,aes(x=generation)) +
            geom_line(aes(y = mean_spring_cost, colour="Spring")) +
            geom_line(aes(y = mean_autumn_cost, colour = "Autumn")) + 
            theme_classic() + 
            xlab("Generation") + 
            #ylim(c(0,10)) +
            ylab("Mean cost")

p4 <- ggplot(data=the.data
        ,aes(x=generation)) +
            geom_line(aes(y = breeder_pop, colour="N breeders")) +
            theme_classic() + 
            xlab("Generation") + 
            ylab(expression("N"[breeder]))

p5 <- ggplot(data=the.data
        ,aes(x=generation)) +
            geom_line(aes(y = offspring_pop, colour="N offspring")) +
            theme_classic() + 
            xlab("Generation") + 
            ylab(expression("N"[offspring]))

p6 <- ggplot(data=the.data
        ,aes(x=generation)) +
            geom_line(aes(y = mean_resources_summer, colour="Winter")) +
            geom_line(aes(y = mean_resources_winter, colour="Summer")) +
            theme_classic() + 
            xlab("Generation") + 
            ylab("Mean resources")


p7 <- ggplot(data=the.data
        ,aes(x=generation)) +
            geom_line(aes(y = mean_theta_a_winter, colour="Psignal (resources), theta winter")) +
            #geom_line(aes(y = mean_theta_a_summer, colour="Stage (resources), theta summer")) +
            geom_line(aes(y = mean_phi_a_winter, colour="Pdisperse (group size), phi winter")) +
            #geom_line(aes(y = mean_phi_a_summer, colour="Disperse (group size), phi summer")) +
            theme_classic() + 
            xlab("Generation") + 
            ylab("Elevation")

p8 <- ggplot(data=the.data
        ,aes(x=generation)) +
            geom_line(aes(y = mean_theta_b_winter, colour="Psignal (resources), theta winter")) +
            #geom_line(aes(y = mean_theta_b_summer, colour="Stage (resources), theta summer")) +
            geom_line(aes(y = mean_phi_b_winter, colour="Pdisperse (group size), phi winter")) +
            #geom_line(aes(y = mean_phi_b_summer, colour="Disperse (group size), phi summer")) +
            theme_classic() + 
            xlab("Generation") + 
            ylab("Slope")

big_plot <- arrangeGrob(p1.a, p1.b, p2, p3, p3b, p3c, p3d, p3e, p3f, p4, p5, p6, p7, p8, nrow=14,ncol=1)
the.base.name <- basename(args[1])

output_file_name <- paste(
        "graph_"
        ,the.base.name
        ,".pdf"
        ,sep="")

ggsave(output_file_name, big_plot, height = 25)

