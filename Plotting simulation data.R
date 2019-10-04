rm(list=ls())

setwd("~/migration")

prefix <- "sim_migration_18_6_2019"  # ENTER FILE PREFIX

myFiles <- list.files(pattern = prefix)
myFiles <- myFiles[2:length(myFiles)]

for (j in 1:length(myFiles)){
	label <- paste("sim_", j, sep = "")
	sim.data <- read.csv(myFiles[j], sep = ";", header = T, row.names = NULL, stringsAsFactors = F)	
	sim.data <- sim.data[1:5000, 1:ncol(sim.data)]
	assign(label, sim.data)
	}
	
#data_1 <- read.csv("sim_migration_18_6_2019_113552_-483780392", sep = ";", header = T, row.names = NULL, stringsAsFactors = F)
#dim(simulation_1)

par(mar = c(5, 4, 1, 1), mfrow = c(6, 2))
plot(y = sim_1$mean_theta_a,
	 x = as.numeric(sim_1$generation),
	 type = "l",
	 main = "Baseline departure rate",
	 ylab = "",
	 xlab = "Generation",
	 bty = "l",
	 las = 1
	 )
	 
plot(y = sim_1$mean_theta_b,
	 x = as.numeric(sim_1$generation),
	 type = "l",
	 main = "Resource-dependent reaction norm for departure",
	 ylab = "",
	 xlab = "Generation",
	 bty = "l",
	 las = 1
	 )
	 
plot(y = sim_1$mean_phi_a,
	 x = as.numeric(sim_1$generation),
	 type = "l",
	 main = "Collective dispersal baseline",
	 ylab = "",
	 xlab = "Generation",
	 bty = "l",
	 las = 1
	 )
	 
plot(y = sim_1$mean_phi_b,
	 x = as.numeric(sim_1$generation),
	 type = "l",
	 main = "Group size-dependent reaction norm for departure",
	 ylab = "",
	 xlab = "Generation",
	 bty = "l",
	 las = 1
	 )
	 
plot(y = sim_1$var_theta_a,
	 x = as.numeric(sim_1$generation),
	 type = "l",
	 main = "Variance in baseline departure rate",
	 ylab = "",
	 xlab = "Generation",
	 bty = "l",
	 las = 1
	 )
	 
plot(y = sim_1$var_theta_b,
	 x = as.numeric(sim_1$generation),
	 type = "l",
	 main = "Variance in resource-dependent reaction norm for departure",
	 ylab = "",
	 xlab = "Generation",
	 bty = "l",
	 las = 1
	 )
	 
plot(y = sim_1$var_phi_a,
	 x = as.numeric(sim_1$generation),
	 type = "l",
	 main = "Variance in collective dispersal baseline",
	 ylab = "",
	 xlab = "Generation",
	 bty = "l",
	 las = 1
	 )
	 
plot(y = sim_1$var_phi_b,
	 x = as.numeric(sim_1$generation),
	 type = "l",
	 main = "Variance in group size-dependent reaction norm for departure",
	 ylab = "",
	 xlab = "Generation",
	 bty = "l",
	 las = 1
	 )

#############################################

####################
# POPULATION SIZES #
####################
 
par(mfcol = c(3,length(myFiles)/3), mar = c(3.5, 3.7, 1, 0.7))

ymax <- max(sim_1[3:5000,17:20], sim_2[3:5000,17:20], sim_3[3:5000,17:20], sim_4[3:5000,17:20], sim_5[3:5000,17:20], sim_6[3:5000,17:20], sim_7[3:5000,17:20], sim_8[3:5000,17:20], sim_9[3:5000,17:20], sim_10[3:5000,17:20], sim_11[3:5000,17:20], sim_12[3:5000,17:20]) # Define y-axis limits by finding maximum population size
	
plot(y = sim_1$mean_staging_size_winter,
	 x = as.numeric(sim_1$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 1",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax)
	 )
	 
mtext(text = expression(italic("n")), side = 2, line = 2.8, cex = 0.7)	 
lines(y = sim_1$mean_flock_size_summer,
	 x = as.numeric(sim_1$generation),
	 col = "red",
	 )	 
lines(y = sim_1$mean_staging_size_summer,
	 x = as.numeric(sim_1$generation),
	 col = "brown",
	 )	 
lines(y = sim_1$mean_flock_size_winter,
	 x = as.numeric(sim_1$generation),
	 col = "lightblue",
	 )	

plot(y = sim_2$mean_staging_size_winter,
	 x = as.numeric(sim_2$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 2",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax)
	 ) 
mtext(text = expression(italic("n")), side = 2, line = 2.8, cex = 0.7)	 
lines(y = sim_2$mean_flock_size_summer,
	 x = as.numeric(sim_2$generation),
	 col = "red",
	 )	 
lines(y = sim_2$mean_staging_size_summer,
	 x = as.numeric(sim_2$generation),
	 col = "brown",
	 )	 
lines(y = sim_2$mean_flock_size_winter,
	 x = as.numeric(sim_2$generation),
	 col = "lightblue"
	 )
	
plot(y = sim_3$mean_staging_size_winter,
	 x = as.numeric(sim_3$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 3",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax)
	 ) 
mtext(text = expression(italic("n")), side = 2, line = 2.8, cex = 0.7)	 
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 
lines(y = sim_3$mean_flock_size_summer,
	 x = as.numeric(sim_3$generation),
	 col = "red",
	 )	 
lines(y = sim_3$mean_staging_size_summer,
	 x = as.numeric(sim_3$generation),
	 col = "brown",
	 )	 
lines(y = sim_3$mean_flock_size_winter,
	 x = as.numeric(sim_3$generation),
	 col = "lightblue"
	 )
	
plot(y = sim_4$mean_staging_size_winter,
	 x = as.numeric(sim_4$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 4",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax)
	 ) 
lines(y = sim_4$mean_flock_size_summer,
	 x = as.numeric(sim_4$generation),
	 col = "red",
	 )	 
lines(y = sim_4$mean_staging_size_summer,
	 x = as.numeric(sim_4$generation),
	 col = "brown",
	 )	 
lines(y = sim_4$mean_flock_size_winter,
	 x = as.numeric(sim_4$generation),
	 col = "lightblue"
	 )
	
plot(y = sim_5$mean_staging_size_winter,
	 x = as.numeric(sim_5$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 5",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax)
	 ) 
lines(y = sim_5$mean_flock_size_summer,
	 x = as.numeric(sim_5$generation),
	 col = "red",
	 )	 
lines(y = sim_5$mean_staging_size_summer,
	 x = as.numeric(sim_5$generation),
	 col = "brown",
	 )	 
lines(y = sim_5$mean_flock_size_winter,
	 x = as.numeric(sim_5$generation),
	 col = "lightblue"
	 )

plot(y = sim_6$mean_staging_size_winter,
	 x = as.numeric(sim_6$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 6",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax)
	 ) 
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 
lines(y = sim_6$mean_flock_size_summer,
	 x = as.numeric(sim_6$generation),
	 col = "red",
	 )	 
lines(y = sim_6$mean_staging_size_summer,
	 x = as.numeric(sim_6$generation),
	 col = "brown",
	 )	 
lines(y = sim_6$mean_flock_size_winter,
	 x = as.numeric(sim_6$generation),
	 col = "lightblue"
	 )
	
plot(y = sim_7$mean_staging_size_winter,
	 x = as.numeric(sim_7$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 7",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax)
	 ) 
lines(y = sim_7$mean_flock_size_summer,
	 x = as.numeric(sim_7$generation),
	 col = "red",
	 )	 
lines(y = sim_7$mean_staging_size_summer,
	 x = as.numeric(sim_7$generation),
	 col = "brown",
	 )	 
lines(y = sim_7$mean_flock_size_winter,
	 x = as.numeric(sim_7$generation),
	 col = "lightblue"
	 )
	
plot(y = sim_8$mean_staging_size_winter,
	 x = as.numeric(sim_8$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 8",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax)
	 ) 
lines(y = sim_8$mean_flock_size_summer,
	 x = as.numeric(sim_8$generation),
	 col = "red",
	 )	 
lines(y = sim_8$mean_staging_size_summer,
	 x = as.numeric(sim_8$generation),
	 col = "brown",
	 )	 
lines(y = sim_8$mean_flock_size_winter,
	 x = as.numeric(sim_8$generation),
	 col = "lightblue"
	 )
	
plot(y = sim_9$mean_staging_size_winter,
	 x = as.numeric(sim_9$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 9",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax)
	 ) 
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 
lines(y = sim_9$mean_flock_size_summer,
	 x = as.numeric(sim_9$generation),
	 col = "red",
	 )	 
lines(y = sim_9$mean_staging_size_summer,
	 x = as.numeric(sim_9$generation),
	 col = "brown",
	 )	 
lines(y = sim_9$mean_flock_size_winter,
	 x = as.numeric(sim_9$generation),
	 col = "lightblue"
	 )
	 
plot(y = sim_10$mean_staging_size_winter,
	 x = as.numeric(sim_10$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 10",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax)
	 ) 
lines(y = sim_10$mean_flock_size_summer,
	 x = as.numeric(sim_10$generation),
	 col = "red",
	 )	 
lines(y = sim_10$mean_staging_size_summer,
	 x = as.numeric(sim_10$generation),
	 col = "brown",
	 )	 
lines(y = sim_10$mean_flock_size_winter,
	 x = as.numeric(sim_10$generation),
	 col = "lightblue"
	 )
	
plot(y = sim_11$mean_staging_size_winter,
	 x = as.numeric(sim_11$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 11",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax)
	 ) 
lines(y = sim_11$mean_flock_size_summer,
	 x = as.numeric(sim_11$generation),
	 col = "red",
	 )	 
lines(y = sim_11$mean_staging_size_summer,
	 x = as.numeric(sim_11$generation),
	 col = "brown",
	 )	 
lines(y = sim_11$mean_flock_size_winter,
	 x = as.numeric(sim_11$generation),
	 col = "lightblue"
	 )
	
plot(y = sim_12$mean_staging_size_winter,
	 x = as.numeric(sim_12$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 12",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax)
	 ) 
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 	 
lines(y = sim_12$mean_flock_size_summer,
	 x = as.numeric(sim_12$generation),
	 col = "red",
	 )	 
lines(y = sim_12$mean_staging_size_summer,
	 x = as.numeric(sim_12$generation),
	 col = "brown",
	 )	 
lines(y = sim_12$mean_flock_size_winter,
	 x = as.numeric(sim_12$generation),
	 col = "lightblue"
	 )
legend("topright", legend = c("Winter", "Spring migration", "Summer", "Autumn migration"), lty = c(1,1,1,1), col = c("lightblue", "green", "red", "brown"), bty = "n", cex = 0.75)



################################
# VARIANCE IN POPULATION SIZES #
################################

dev.new()

ymax.var <- max(sim_1[3:5000,21:24], sim_2[3:5000,21:24], sim_3[3:5000,21:24], sim_4[3:5000,21:24], sim_5[3:5000,21:24], sim_6[3:5000,21:24], sim_7[3:5000,21:24], sim_8[3:5000,21:24], sim_9[3:5000,21:24], sim_10[3:5000,21:24], sim_11[3:5000,21:24], sim_12[3:5000,21:24]) # Define y-axis limits by finding maximum population size

par(mfcol = c(3,length(myFiles)/3), mar = c(3.5, 4.5, 1, 0.8))

plot(y = sim_1$var_staging_size_winter,
	 x = as.numeric(sim_1$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 1",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.var)
	 )
	 
mtext(text = "Variance", side = 2, line = 3.5, cex = 0.7)	 
lines(y = sim_1$var_flock_size_summer,
	 x = as.numeric(sim_1$generation),
	 col = "red",
	 )	 
lines(y = sim_1$var_staging_size_summer,
	 x = as.numeric(sim_1$generation),
	 col = "brown",
	 )	 
lines(y = sim_1$var_flock_size_winter,
	 x = as.numeric(sim_1$generation),
	 col = "lightblue",
	 )	

plot(y = sim_2$var_staging_size_winter,
	 x = as.numeric(sim_2$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 2",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.var)
	 ) 
mtext(text = "Variance", side = 2, line = 3.5, cex = 0.7)	 
lines(y = sim_2$var_flock_size_summer,
	 x = as.numeric(sim_2$generation),
	 col = "red",
	 )	 
lines(y = sim_2$var_staging_size_summer,
	 x = as.numeric(sim_2$generation),
	 col = "brown",
	 )	 
lines(y = sim_2$var_flock_size_winter,
	 x = as.numeric(sim_2$generation),
	 col = "lightblue"
	 )
	
plot(y = sim_3$var_staging_size_winter,
	 x = as.numeric(sim_3$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 3",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.var)
	 ) 
mtext(text = "Variance", side = 2, line = 3.5, cex = 0.7)	 
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 
lines(y = sim_3$var_flock_size_summer,
	 x = as.numeric(sim_3$generation),
	 col = "red",
	 )	 
lines(y = sim_3$var_staging_size_summer,
	 x = as.numeric(sim_3$generation),
	 col = "brown",
	 )	 
lines(y = sim_3$var_flock_size_winter,
	 x = as.numeric(sim_3$generation),
	 col = "lightblue"
	 )
	
plot(y = sim_4$var_staging_size_winter,
	 x = as.numeric(sim_4$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 4",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.var)
	 ) 
lines(y = sim_4$var_flock_size_summer,
	 x = as.numeric(sim_4$generation),
	 col = "red",
	 )	 
lines(y = sim_4$var_staging_size_summer,
	 x = as.numeric(sim_4$generation),
	 col = "brown",
	 )	 
lines(y = sim_4$var_flock_size_winter,
	 x = as.numeric(sim_4$generation),
	 col = "lightblue"
	 )
	
plot(y = sim_5$var_staging_size_winter,
	 x = as.numeric(sim_5$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 5",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.var)
	 ) 
lines(y = sim_5$var_flock_size_summer,
	 x = as.numeric(sim_5$generation),
	 col = "red",
	 )	 
lines(y = sim_5$var_staging_size_summer,
	 x = as.numeric(sim_5$generation),
	 col = "brown",
	 )	 
lines(y = sim_5$var_flock_size_winter,
	 x = as.numeric(sim_5$generation),
	 col = "lightblue"
	 )

plot(y = sim_6$var_staging_size_winter,
	 x = as.numeric(sim_6$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 6",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.var)
	 ) 
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 
lines(y = sim_6$var_flock_size_summer,
	 x = as.numeric(sim_6$generation),
	 col = "red",
	 )	 
lines(y = sim_6$var_staging_size_summer,
	 x = as.numeric(sim_6$generation),
	 col = "brown",
	 )	 
lines(y = sim_6$var_flock_size_winter,
	 x = as.numeric(sim_6$generation),
	 col = "lightblue"
	 )
	
plot(y = sim_7$var_staging_size_winter,
	 x = as.numeric(sim_7$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 7",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.var)
	 ) 
lines(y = sim_7$var_flock_size_summer,
	 x = as.numeric(sim_7$generation),
	 col = "red",
	 )	 
lines(y = sim_7$var_staging_size_summer,
	 x = as.numeric(sim_7$generation),
	 col = "brown",
	 )	 
lines(y = sim_7$var_flock_size_winter,
	 x = as.numeric(sim_7$generation),
	 col = "lightblue"
	 )
	
plot(y = sim_8$var_staging_size_winter,
	 x = as.numeric(sim_8$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 8",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.var)
	 ) 
lines(y = sim_8$var_flock_size_summer,
	 x = as.numeric(sim_8$generation),
	 col = "red",
	 )	 
lines(y = sim_8$var_staging_size_summer,
	 x = as.numeric(sim_8$generation),
	 col = "brown",
	 )	 
lines(y = sim_8$var_flock_size_winter,
	 x = as.numeric(sim_8$generation),
	 col = "lightblue"
	 )
	
plot(y = sim_9$var_staging_size_winter,
	 x = as.numeric(sim_9$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 9",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.var)
	 ) 
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 
lines(y = sim_9$var_flock_size_summer,
	 x = as.numeric(sim_9$generation),
	 col = "red",
	 )	 
lines(y = sim_9$var_staging_size_summer,
	 x = as.numeric(sim_9$generation),
	 col = "brown",
	 )	 
lines(y = sim_9$var_flock_size_winter,
	 x = as.numeric(sim_9$generation),
	 col = "lightblue"
	 )
	 
plot(y = sim_10$var_staging_size_winter,
	 x = as.numeric(sim_10$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 10",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.var)
	 ) 
lines(y = sim_10$var_flock_size_summer,
	 x = as.numeric(sim_10$generation),
	 col = "red",
	 )	 
lines(y = sim_10$var_staging_size_summer,
	 x = as.numeric(sim_10$generation),
	 col = "brown",
	 )	 
lines(y = sim_10$var_flock_size_winter,
	 x = as.numeric(sim_10$generation),
	 col = "lightblue"
	 )
	
plot(y = sim_11$var_staging_size_winter,
	 x = as.numeric(sim_11$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 11",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.var)
	 ) 
lines(y = sim_11$var_flock_size_summer,
	 x = as.numeric(sim_11$generation),
	 col = "red",
	 )	 
lines(y = sim_11$var_staging_size_summer,
	 x = as.numeric(sim_11$generation),
	 col = "brown",
	 )	 
lines(y = sim_11$var_flock_size_winter,
	 x = as.numeric(sim_11$generation),
	 col = "lightblue"
	 )
	
plot(y = sim_12$var_staging_size_winter,
	 x = as.numeric(sim_12$generation),
	 type = "l",
	 col = "green",
	 main = "Simulation 12",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.var)
	 ) 
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 	 
lines(y = sim_12$var_flock_size_summer,
	 x = as.numeric(sim_12$generation),
	 col = "red",
	 )	 
lines(y = sim_12$var_staging_size_summer,
	 x = as.numeric(sim_12$generation),
	 col = "brown",
	 )	 
lines(y = sim_12$var_flock_size_winter,
	 x = as.numeric(sim_12$generation),
	 col = "lightblue"
	 )
legend("topright", legend = c("Winter", "Spring migration", "Summer", "Autumn migration"), lty = c(1,1,1,1), col = c("lightblue", "green", "red", "brown"), bty = "n", cex = 0.75)


##################
# RESOURCE LEVEL #
##################

dev.new()

ymax.mean_resources <- max(sim_1[4:5000,7], sim_2[4:5000,7], sim_3[4:5000,7], sim_4[4:5000,7], sim_5[4:5000,7], sim_6[4:5000,7], sim_7[4:5000,7], sim_8[4:5000,7], sim_9[4:5000,7], sim_10[4:5000,7], sim_11[4:5000,7], sim_12[4:5000,7]) # Define y-axis limits by finding maximum resource level (ignoring first 20 generations)

ymax.var_resources <- max(sim_1[3:5000,12], sim_2[3:5000,12], sim_3[3:5000,12], sim_4[3:5000,12], sim_5[3:5000,12], sim_6[3:5000,12], sim_7[3:5000,12], sim_8[3:5000,12], sim_9[3:5000,12], sim_10[3:5000,12], sim_11[3:5000,12], sim_12[3:5000,12]) # Define y-axis limits by finding maximum variance in resource level

par(mfcol = c(3,length(myFiles)/3), mar = c(3.5, 3.7, 1, 3.7))

# SIMULATION 1
plot(y = sim_1$mean_resources,
	 x = as.numeric(sim_1$generation),
	 type = "l",
	 col = "red",
	 main = "Simulation 1",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_resources)
	 ) 
mtext(text = "Mean", side = 2, line = 2.2, cex = 0.7)	 
par(new = T)
plot(y = sim_1$var_resources,
	 x = as.numeric(sim_1$generation),
	 type = "l",
	 col = "green",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_resources)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

# SIMULTATION 2
plot(y = sim_2$mean_resources,
	 x = as.numeric(sim_2$generation),
	 type = "l",
	 col = "red",
	 main = "Simulation 2",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_resources)
	 ) 
mtext(text = "Mean", side = 2, line = 2.2, cex = 0.7)	 
par(new = T)
plot(y = sim_2$var_resources,
	 x = as.numeric(sim_2$generation),
	 type = "l",
	 col = "green",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_resources)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 	 

# SIMULATION 3
plot(y = sim_3$mean_resources,
	 x = as.numeric(sim_3$generation),
	 type = "l",
	 col = "red",
	 main = "Simulation 3",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_resources)
	 ) 
mtext(text = "Mean", side = 2, line = 2.2, cex = 0.7)	 
par(new = T)
plot(y = sim_3$var_resources,
	 x = as.numeric(sim_3$generation),
	 type = "l",
	 col = "green",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_resources)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 	 

# SIMULTATION 4
plot(y = sim_4$mean_resources,
	 x = as.numeric(sim_4$generation),
	 type = "l",
	 col = "red",
	 main = "Simulation 4",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_resources)
	 ) 
mtext(text = "Mean", side = 2, line = 2.2, cex = 0.7)	 
par(new = T)
plot(y = sim_4$var_resources,
	 x = as.numeric(sim_4$generation),
	 type = "l",
	 col = "green",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_resources)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7) 

# SIMULATION 5
plot(y = sim_5$mean_resources,
	 x = as.numeric(sim_5$generation),
	 type = "l",
	 col = "red",
	 main = "Simulation 5",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_resources)
	 ) 
mtext(text = "Mean", side = 2, line = 2.2, cex = 0.7)	 
par(new = T)
plot(y = sim_5$var_resources,
	 x = as.numeric(sim_5$generation),
	 type = "l",
	 col = "green",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_resources)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

# SIMULTATION 6
plot(y = sim_6$mean_resources,
	 x = as.numeric(sim_6$generation),
	 type = "l",
	 col = "red",
	 main = "Simulation 6",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_resources)
	 ) 
mtext(text = "Mean", side = 2, line = 2.2, cex = 0.7)	 
par(new = T)
plot(y = sim_6$var_resources,
	 x = as.numeric(sim_6$generation),
	 type = "l",
	 col = "green",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_resources)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 

# SIMULATION 7
plot(y = sim_7$mean_resources,
	 x = as.numeric(sim_7$generation),
	 type = "l",
	 col = "red",
	 main = "Simulation 7",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_resources)
	 ) 
mtext(text = "Mean", side = 2, line = 2.2, cex = 0.7)	 
par(new = T)
plot(y = sim_7$var_resources,
	 x = as.numeric(sim_7$generation),
	 type = "l",
	 col = "green",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_resources)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

# SIMULTATION 8
plot(y = sim_8$mean_resources,
	 x = as.numeric(sim_8$generation),
	 type = "l",
	 col = "red",
	 main = "Simulation 8",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_resources)
	 ) 
mtext(text = "Mean", side = 2, line = 2.2, cex = 0.7)	 
par(new = T)
plot(y = sim_8$var_resources,
	 x = as.numeric(sim_8$generation),
	 type = "l",
	 col = "green",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_resources)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 	 

# SIMULATION 9
plot(y = sim_9$mean_resources,
	 x = as.numeric(sim_9$generation),
	 type = "l",
	 col = "red",
	 main = "Simulation 9",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_resources)
	 ) 
mtext(text = "Mean", side = 2, line = 2.2, cex = 0.7)	 
par(new = T)
plot(y = sim_9$var_resources,
	 x = as.numeric(sim_9$generation),
	 type = "l",
	 col = "green",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_resources)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 

# SIMULTATION 10
plot(y = sim_10$mean_resources,
	 x = as.numeric(sim_10$generation),
	 type = "l",
	 col = "red",
	 main = "Simulation 10",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_resources)
	 ) 
mtext(text = "Mean", side = 2, line = 2.2, cex = 0.7)	 
par(new = T)
plot(y = sim_10$var_resources,
	 x = as.numeric(sim_10$generation),
	 type = "l",
	 col = "green",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_resources)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7) 

# SIMULATION 11
plot(y = sim_11$mean_resources,
	 x = as.numeric(sim_11$generation),
	 type = "l",
	 col = "red",
	 main = "Simulation 11",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_resources)
	 ) 
mtext(text = "Mean", side = 2, line = 2.2, cex = 0.7)	 
par(new = T)
plot(y = sim_11$var_resources,
	 x = as.numeric(sim_11$generation),
	 type = "l",
	 col = "green",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_resources)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

# SIMULTATION 12
plot(y = sim_12$mean_resources,
	 x = as.numeric(sim_12$generation),
	 type = "l",
	 col = "red",
	 main = "Simulation 12",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_resources)
	 ) 
mtext(text = "Mean", side = 2, line = 2.2, cex = 0.7)
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 	 
par(new = T)
plot(y = sim_12$var_resources,
	 x = as.numeric(sim_12$generation),
	 type = "l",
	 col = "green",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_resources)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

legend("topright", legend = c("Resources mean", "Resources variance"), lty = c(1,1,1,1), col = c("red", "green"), bty = "n", cex = 0.75)



#######################################
# THETA: RESOURCE-DEPENDENT DEPARTURE #
#######################################

# THETA_A

dev.new()
						   
ymax.mean_theta_a <- max(c(sim_1[3:5000,3], sim_2[3:5000,3], sim_3[3:5000,3], sim_4[3:5000,3], sim_5[3:5000,3], sim_6[3:5000,3], sim_7[3:5000,3], sim_8[3:5000,3], sim_9[3:5000,3], sim_10[3:5000,3], sim_11[3:5000,3], sim_12[3:5000,3]))

ymax.var_theta_a <- max(c(sim_1[3:5000,8], sim_2[3:5000,8], sim_3[3:5000,8], sim_4[3:5000,8], sim_5[3:5000,8], sim_6[3:5000,8], sim_7[3:5000,8], sim_8[3:5000,8], sim_9[3:5000,8], sim_10[3:5000,8], sim_11[3:5000,8], sim_12[3:5000,8]))

par(mfcol = c(3,length(myFiles)/3), mar = c(4, 3.7, 1, 4))

# SIMULATION 1
plot(y = sim_1$mean_theta_a,
	 x = as.numeric(sim_1$generation),
	 type = "l",
	 col = "blue",
	 main = "Simulation 1",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_1$var_theta_a,
	 x = as.numeric(sim_1$generation),
	 type = "l",
	 col = "yellow",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_a)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

# SIMULTATION 2
plot(y = sim_2$mean_theta_a,
	 x = as.numeric(sim_2$generation),
	 type = "l",
	 col = "blue",
	 main = "Simulation 2",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_2$var_theta_a,
	 x = as.numeric(sim_2$generation),
	 type = "l",
	 col = "yellow",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_a)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 	 

# SIMULATION 3
plot(y = sim_3$mean_theta_a,
	 x = as.numeric(sim_3$generation),
	 type = "l",
	 col = "blue",
	 main = "Simulation 3",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_3$var_theta_a,
	 x = as.numeric(sim_3$generation),
	 type = "l",
	 col = "yellow",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_a)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 	 

# SIMULTATION 4
plot(y = sim_4$mean_theta_a,
	 x = as.numeric(sim_4$generation),
	 type = "l",
	 col = "blue",
	 main = "Simulation 4",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_4$var_theta_a,
	 x = as.numeric(sim_4$generation),
	 type = "l",
	 col = "yellow",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_a)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7) 

# SIMULATION 5
plot(y = sim_5$mean_theta_a,
	 x = as.numeric(sim_5$generation),
	 type = "l",
	 col = "blue",
	 main = "Simulation 5",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_5$var_theta_a,
	 x = as.numeric(sim_5$generation),
	 type = "l",
	 col = "yellow",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_a)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

# SIMULTATION 6
plot(y = sim_6$mean_theta_a,
	 x = as.numeric(sim_6$generation),
	 type = "l",
	 col = "blue",
	 main = "Simulation 6",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_6$var_theta_a,
	 x = as.numeric(sim_6$generation),
	 type = "l",
	 col = "yellow",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_a)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 

# SIMULATION 7
plot(y = sim_7$mean_theta_a,
	 x = as.numeric(sim_7$generation),
	 type = "l",
	 col = "blue",
	 main = "Simulation 7",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_7$var_theta_a,
	 x = as.numeric(sim_7$generation),
	 type = "l",
	 col = "yellow",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_a)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

# SIMULTATION 8
plot(y = sim_8$mean_theta_a,
	 x = as.numeric(sim_8$generation),
	 type = "l",
	 col = "blue",
	 main = "Simulation 8",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_8$var_theta_a,
	 x = as.numeric(sim_8$generation),
	 type = "l",
	 col = "yellow",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_a)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 	 

# SIMULATION 9
plot(y = sim_9$mean_theta_a,
	 x = as.numeric(sim_9$generation),
	 type = "l",
	 col = "blue",
	 main = "Simulation 9",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_9$var_theta_a,
	 x = as.numeric(sim_9$generation),
	 type = "l",
	 col = "yellow",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_a),
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 

# SIMULTATION 10
plot(y = sim_10$mean_theta_a,
	 x = as.numeric(sim_10$generation),
	 type = "l",
	 col = "blue",
	 main = "Simulation 10",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_10$var_theta_a,
	 x = as.numeric(sim_10$generation),
	 type = "l",
	 col = "yellow",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_a),
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7) 

# SIMULATION 11
plot(y = sim_11$mean_theta_a,
	 x = as.numeric(sim_11$generation),
	 type = "l",
	 col = "blue",
	 main = "Simulation 11",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_11$var_theta_a,
	 x = as.numeric(sim_11$generation),
	 type = "l",
	 col = "yellow",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_a),
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

# SIMULTATION 12
plot(y = sim_12$mean_theta_a,
	 x = as.numeric(sim_12$generation),
	 type = "l",
	 col = "blue",
	 main = "Simulation 12",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 	 
par(new = T)
plot(y = sim_12$var_theta_a,
	 x = as.numeric(sim_12$generation),
	 type = "l",
	 col = "yellow",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_a),
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

legend("topright", legend = c("Mean(theta_a)", "Var(theta_a)"), lty = c(1,1,1,1), col = c("blue", "yellow"), bty = "n", cex = 0.75)


# THETA_B
dev.new()
						   
ymax.mean_theta_b <- max(c(sim_1[3:5000,4], sim_2[3:5000,4], sim_3[3:5000,4], sim_4[3:5000,4], sim_5[3:5000,4], sim_6[3:5000,4], sim_7[3:5000,4], sim_8[3:5000,4], sim_9[3:5000,4], sim_10[3:5000,4], sim_11[3:5000,4], sim_12[3:5000,4]))

ymax.var_theta_b <- max(c(sim_1[3:5000,9], sim_2[3:5000,9], sim_3[3:5000,9], sim_4[3:5000,9], sim_5[3:5000,9], sim_6[3:5000,9], sim_7[3:5000,9], sim_8[3:5000,9], sim_9[3:5000,9], sim_10[3:5000,9], sim_11[3:5000,9], sim_12[3:5000,9]))

par(mfcol = c(3,length(myFiles)/3), mar = c(4, 3.7, 1, 4))

# SIMULATION 1
plot(y = sim_1$mean_theta_b,
	 x = as.numeric(sim_1$generation),
	 type = "l",
	 col = "purple",
	 main = "Simulation 1",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_1$var_theta_b,
	 x = as.numeric(sim_1$generation),
	 type = "l",
	 col = "orange",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_b)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

# SIMULTATION 2
plot(y = sim_2$mean_theta_b,
	 x = as.numeric(sim_2$generation),
	 type = "l",
	 col = "purple",
	 main = "Simulation 2",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_2$var_theta_b,
	 x = as.numeric(sim_2$generation),
	 type = "l",
	 col = "orange",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_b)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 	 

# SIMULATION 3
plot(y = sim_3$mean_theta_b,
	 x = as.numeric(sim_3$generation),
	 type = "l",
	 col = "purple",
	 main = "Simulation 3",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_3$var_theta_b,
	 x = as.numeric(sim_3$generation),
	 type = "l",
	 col = "orange",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_b)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 	 

# SIMULTATION 4
plot(y = sim_4$mean_theta_b,
	 x = as.numeric(sim_4$generation),
	 type = "l",
	 col = "purple",
	 main = "Simulation 4",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_4$var_theta_b,
	 x = as.numeric(sim_4$generation),
	 type = "l",
	 col = "orange",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_b)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7) 

# SIMULATION 5
plot(y = sim_5$mean_theta_b,
	 x = as.numeric(sim_5$generation),
	 type = "l",
	 col = "purple",
	 main = "Simulation 5",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_5$var_theta_b,
	 x = as.numeric(sim_5$generation),
	 type = "l",
	 col = "orange",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_b)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

# SIMULTATION 6
plot(y = sim_6$mean_theta_b,
	 x = as.numeric(sim_6$generation),
	 type = "l",
	 col = "purple",
	 main = "Simulation 6",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_6$var_theta_b,
	 x = as.numeric(sim_6$generation),
	 type = "l",
	 col = "orange",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_b)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 

# SIMULATION 7
plot(y = sim_7$mean_theta_b,
	 x = as.numeric(sim_7$generation),
	 type = "l",
	 col = "purple",
	 main = "Simulation 7",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_7$var_theta_b,
	 x = as.numeric(sim_7$generation),
	 type = "l",
	 col = "orange",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_b)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

# SIMULTATION 8
plot(y = sim_8$mean_theta_b,
	 x = as.numeric(sim_8$generation),
	 type = "l",
	 col = "purple",
	 main = "Simulation 8",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_8$var_theta_b,
	 x = as.numeric(sim_8$generation),
	 type = "l",
	 col = "orange",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_b)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 	 

# SIMULATION 9
plot(y = sim_9$mean_theta_b,
	 x = as.numeric(sim_9$generation),
	 type = "l",
	 col = "purple",
	 main = "Simulation 9",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_9$var_theta_b,
	 x = as.numeric(sim_9$generation),
	 type = "l",
	 col = "orange",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_b),
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 

# SIMULTATION 10
plot(y = sim_10$mean_theta_b,
	 x = as.numeric(sim_10$generation),
	 type = "l",
	 col = "purple",
	 main = "Simulation 10",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_10$var_theta_b,
	 x = as.numeric(sim_10$generation),
	 type = "l",
	 col = "orange",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_b),
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7) 

# SIMULATION 11
plot(y = sim_11$mean_theta_b,
	 x = as.numeric(sim_11$generation),
	 type = "l",
	 col = "purple",
	 main = "Simulation 11",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_11$var_theta_b,
	 x = as.numeric(sim_11$generation),
	 type = "l",
	 col = "orange",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_b),
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

# SIMULTATION 12
plot(y = sim_12$mean_theta_b,
	 x = as.numeric(sim_12$generation),
	 type = "l",
	 col = "purple",
	 main = "Simulation 12",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_theta_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 	 
par(new = T)
plot(y = sim_12$var_theta_b,
	 x = as.numeric(sim_12$generation),
	 type = "l",
	 col = "orange",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_theta_b),
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

legend("topright", legend = c("Mean(theta_b)", "Var(theta_b)"), lty = c(1,1,1,1), col = c("purple", "orange"), bty = "n", cex = 0.75)


#####################################
# PHI: RESOURCE-DEPENDENT DEPARTURE #
#####################################

# PHI_A

dev.new()
						   
ymax.mean_phi_a <- max(c(sim_1[3:5000,5], sim_2[3:5000,5], sim_3[3:5000,5], sim_4[3:5000,5], sim_5[3:5000,5], sim_6[3:5000,5], sim_7[3:5000,5], sim_8[3:5000,5], sim_9[3:5000,5], sim_10[3:5000,5], sim_11[3:5000,5], sim_12[3:5000,5]))

ymax.var_phi_a <- max(c(sim_1[3:5000,10], sim_2[3:5000,10], sim_3[3:5000,10], sim_4[3:5000,10], sim_5[3:5000,10], sim_6[3:5000,10], sim_7[3:5000,10], sim_8[3:5000,10], sim_9[3:5000,10], sim_10[3:5000,10], sim_11[3:5000,10], sim_12[3:5000,10]))

par(mfcol = c(3,length(myFiles)/3), mar = c(4, 3.7, 1, 4))

# SIMULATION 1
plot(y = sim_1$mean_phi_a,
	 x = as.numeric(sim_1$generation),
	 type = "l",
	 col = "darkblue",
	 main = "Simulation 1",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_1$var_phi_a,
	 x = as.numeric(sim_1$generation),
	 type = "l",
	 col = "pink",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_a)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

# SIMULTATION 2
plot(y = sim_2$mean_phi_a,
	 x = as.numeric(sim_2$generation),
	 type = "l",
	 col = "darkblue",
	 main = "Simulation 2",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_2$var_phi_a,
	 x = as.numeric(sim_2$generation),
	 type = "l",
	 col = "pink",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_a)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 	 

# SIMULATION 3
plot(y = sim_3$mean_phi_a,
	 x = as.numeric(sim_3$generation),
	 type = "l",
	 col = "darkblue",
	 main = "Simulation 3",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_3$var_phi_a,
	 x = as.numeric(sim_3$generation),
	 type = "l",
	 col = "pink",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_a)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 	 

# SIMULTATION 4
plot(y = sim_4$mean_phi_a,
	 x = as.numeric(sim_4$generation),
	 type = "l",
	 col = "darkblue",
	 main = "Simulation 4",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_4$var_phi_a,
	 x = as.numeric(sim_4$generation),
	 type = "l",
	 col = "pink",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_a)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7) 

# SIMULATION 5
plot(y = sim_5$mean_phi_a,
	 x = as.numeric(sim_5$generation),
	 type = "l",
	 col = "darkblue",
	 main = "Simulation 5",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_5$var_phi_a,
	 x = as.numeric(sim_5$generation),
	 type = "l",
	 col = "pink",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_a)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

# SIMULTATION 6
plot(y = sim_6$mean_phi_a,
	 x = as.numeric(sim_6$generation),
	 type = "l",
	 col = "darkblue",
	 main = "Simulation 6",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_6$var_phi_a,
	 x = as.numeric(sim_6$generation),
	 type = "l",
	 col = "pink",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_a)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 

# SIMULATION 7
plot(y = sim_7$mean_phi_a,
	 x = as.numeric(sim_7$generation),
	 type = "l",
	 col = "darkblue",
	 main = "Simulation 7",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_7$var_phi_a,
	 x = as.numeric(sim_7$generation),
	 type = "l",
	 col = "pink",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_a)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

# SIMULTATION 8
plot(y = sim_8$mean_phi_a,
	 x = as.numeric(sim_8$generation),
	 type = "l",
	 col = "darkblue",
	 main = "Simulation 8",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_8$var_phi_a,
	 x = as.numeric(sim_8$generation),
	 type = "l",
	 col = "pink",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_a)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 	 

# SIMULATION 9
plot(y = sim_9$mean_phi_a,
	 x = as.numeric(sim_9$generation),
	 type = "l",
	 col = "darkblue",
	 main = "Simulation 9",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_9$var_phi_a,
	 x = as.numeric(sim_9$generation),
	 type = "l",
	 col = "pink",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_a),
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 

# SIMULTATION 10
plot(y = sim_10$mean_phi_a,
	 x = as.numeric(sim_10$generation),
	 type = "l",
	 col = "darkblue",
	 main = "Simulation 10",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_10$var_phi_a,
	 x = as.numeric(sim_10$generation),
	 type = "l",
	 col = "pink",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_a),
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7) 

# SIMULATION 11
plot(y = sim_11$mean_phi_a,
	 x = as.numeric(sim_11$generation),
	 type = "l",
	 col = "darkblue",
	 main = "Simulation 11",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_11$var_phi_a,
	 x = as.numeric(sim_11$generation),
	 type = "l",
	 col = "pink",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_a),
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

# SIMULTATION 12
plot(y = sim_12$mean_phi_a,
	 x = as.numeric(sim_12$generation),
	 type = "l",
	 col = "darkblue",
	 main = "Simulation 12",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_a)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 	 
par(new = T)
plot(y = sim_12$var_phi_a,
	 x = as.numeric(sim_12$generation),
	 type = "l",
	 col = "pink",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_a),
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

legend("topleft", legend = c("Mean(phi_a)", "Var(phi_a)"), lty = c(1,1,1,1), col = c("darkblue", "pink"), bty = "n", cex = 0.75)


# PHI_B
dev.new()
						   
ymax.mean_phi_b <- max(c(sim_1[3:5000,6], sim_2[3:5000,6], sim_3[3:5000,6], sim_4[3:5000,6], sim_5[3:5000,6], sim_6[3:5000,6], sim_7[3:5000,6], sim_8[3:5000,6], sim_9[3:5000,6], sim_10[3:5000,6], sim_11[3:5000,6], sim_12[3:5000,6]))

ymax.var_phi_b <- max(c(sim_1[3:5000,11], sim_2[3:5000,11], sim_3[3:5000,11], sim_4[3:5000,11], sim_5[3:5000,11], sim_6[3:5000,11], sim_7[3:5000,11], sim_8[3:5000,11], sim_9[3:5000,11], sim_10[3:5000,11], sim_11[3:5000,11], sim_12[3:5000,11]))

par(mfcol = c(3,length(myFiles)/3), mar = c(4, 3.7, 1, 4))

# SIMULATION 1
plot(y = sim_1$mean_phi_b,
	 x = as.numeric(sim_1$generation),
	 type = "l",
	 col = "aquamarine",
	 main = "Simulation 1",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_1$var_phi_b,
	 x = as.numeric(sim_1$generation),
	 type = "l",
	 col = "coral",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_b)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

# SIMULTATION 2
plot(y = sim_2$mean_phi_b,
	 x = as.numeric(sim_2$generation),
	 type = "l",
	 col = "aquamarine",
	 main = "Simulation 2",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_2$var_phi_b,
	 x = as.numeric(sim_2$generation),
	 type = "l",
	 col = "coral",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_b)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 	 

# SIMULATION 3
plot(y = sim_3$mean_phi_b,
	 x = as.numeric(sim_3$generation),
	 type = "l",
	 col = "aquamarine",
	 main = "Simulation 3",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_3$var_phi_b,
	 x = as.numeric(sim_3$generation),
	 type = "l",
	 col = "coral",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_b)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 	 

# SIMULTATION 4
plot(y = sim_4$mean_phi_b,
	 x = as.numeric(sim_4$generation),
	 type = "l",
	 col = "aquamarine",
	 main = "Simulation 4",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_4$var_phi_b,
	 x = as.numeric(sim_4$generation),
	 type = "l",
	 col = "coral",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_b)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7) 

# SIMULATION 5
plot(y = sim_5$mean_phi_b,
	 x = as.numeric(sim_5$generation),
	 type = "l",
	 col = "aquamarine",
	 main = "Simulation 5",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_5$var_phi_b,
	 x = as.numeric(sim_5$generation),
	 type = "l",
	 col = "coral",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_b)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

# SIMULTATION 6
plot(y = sim_6$mean_phi_b,
	 x = as.numeric(sim_6$generation),
	 type = "l",
	 col = "aquamarine",
	 main = "Simulation 6",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_6$var_phi_b,
	 x = as.numeric(sim_6$generation),
	 type = "l",
	 col = "coral",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_b)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 

# SIMULATION 7
plot(y = sim_7$mean_phi_b,
	 x = as.numeric(sim_7$generation),
	 type = "l",
	 col = "aquamarine",
	 main = "Simulation 7",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_7$var_phi_b,
	 x = as.numeric(sim_7$generation),
	 type = "l",
	 col = "coral",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_b)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

# SIMULTATION 8
plot(y = sim_8$mean_phi_b,
	 x = as.numeric(sim_8$generation),
	 type = "l",
	 col = "aquamarine",
	 main = "Simulation 8",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_8$var_phi_b,
	 x = as.numeric(sim_8$generation),
	 type = "l",
	 col = "coral",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_b)	,
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 	 

# SIMULATION 9
plot(y = sim_9$mean_phi_b,
	 x = as.numeric(sim_9$generation),
	 type = "l",
	 col = "aquamarine",
	 main = "Simulation 9",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_9$var_phi_b,
	 x = as.numeric(sim_9$generation),
	 type = "l",
	 col = "coral",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_b),
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 

# SIMULTATION 10
plot(y = sim_10$mean_phi_b,
	 x = as.numeric(sim_10$generation),
	 type = "l",
	 col = "aquamarine",
	 main = "Simulation 10",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_10$var_phi_b,
	 x = as.numeric(sim_10$generation),
	 type = "l",
	 col = "coral",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_b),
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7) 

# SIMULATION 11
plot(y = sim_11$mean_phi_b,
	 x = as.numeric(sim_11$generation),
	 type = "l",
	 col = "aquamarine",
	 main = "Simulation 11",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)	 
par(new = T)
plot(y = sim_11$var_phi_b,
	 x = as.numeric(sim_11$generation),
	 type = "l",
	 col = "coral",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_b),
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

# SIMULTATION 12
plot(y = sim_12$mean_phi_b,
	 x = as.numeric(sim_12$generation),
	 type = "l",
	 col = "aquamarine",
	 main = "Simulation 12",
	 ylab = "",
	 xlab = "",
	 bty = "l",
	 las = 1,
	 ylim = c(0, ymax.mean_phi_b)
	 ) 
mtext(text = "Mean", side = 2, line = 2.8, cex = 0.7)
mtext(text = "Generation", side = 1, line = 2.3, cex = 0.7)	 	 
par(new = T)
plot(y = sim_12$var_phi_b,
	 x = as.numeric(sim_12$generation),
	 type = "l",
	 col = "coral",
	 ylab = "",
	 xlab = "",
	 bty = "n",
	 ylim = c(0, ymax.var_phi_b),
	 axes = F
	 )
axis(side = 4)
axis(side = 4, at = c(-1000, 9000))
mtext(text = "Variance", side = 4, line = 2.2, cex = 0.7)	 

legend("topleft", legend = c("Mean(phi_b)", "Var(phi_b)"), lty = c(1,1,1,1), col = c("aquamarine", "coral"), bty = "n", cex = 0.75)


# PLOTTING EVOLUTION OF GROUP SIZE DEPENDENCY FOR ALL TWELVE SIMULATIONS	
ymax.mean_phi_b <- max(c(sim_1[3:5000,6], sim_2[3:5000,6], sim_3[3:5000,6], sim_4[3:5000,6], sim_5[3:5000,6], sim_6[3:5000,6], sim_7[3:5000,6], sim_8[3:5000,6], sim_9[3:5000,6], sim_10[3:5000,6], sim_11[3:5000,6], sim_12[3:5000,6]))

par(mar = c(4, 4, 1, 1))

cl = rainbow(12)

plot(y = sim_1$mean_phi_b,
	 x = as.numeric(sim_1$generation),
	 type = "l",
	 col = cl[1],
	 main = "",
	 ylab = "Social sensitivity of departure timing",
	 xlab = "Time (generations)",
	 bty = "l",
	 las = 1,
	 lwd = 1.5,
	 xaxs = "i",
	 yaxs = "i",
	 ylim = c(0, ymax.mean_phi_b)
	 ) 

lines(y = sim_2$mean_phi_b, x = as.numeric(sim_2$generation), col = cl[2], lwd = 1.5)
lines(y = sim_3$mean_phi_b, x = as.numeric(sim_3$generation), col = cl[3], lwd = 1.5)
lines(y = sim_4$mean_phi_b, x = as.numeric(sim_4$generation), col = cl[4], lwd = 1.5)
lines(y = sim_5$mean_phi_b, x = as.numeric(sim_5$generation), col = cl[5], lwd = 1.5)
lines(y = sim_6$mean_phi_b, x = as.numeric(sim_6$generation), col = cl[6], lwd = 1.5)
lines(y = sim_7$mean_phi_b, x = as.numeric(sim_7$generation), col = cl[7], lwd = 1.5)
lines(y = sim_8$mean_phi_b, x = as.numeric(sim_8$generation), col = cl[8], lwd = 1.5)
lines(y = sim_9$mean_phi_b, x = as.numeric(sim_9$generation), col = cl[9], lwd = 1.5)
lines(y = sim_10$mean_phi_b, x = as.numeric(sim_10$generation), col = cl[10], lwd = 1.5)
lines(y = sim_11$mean_phi_b, x = as.numeric(sim_11$generation), col = cl[11], lwd = 1.5)
lines(y = sim_12$mean_phi_b, x = as.numeric(sim_12$generation), col = cl[12], lwd = 1.5)