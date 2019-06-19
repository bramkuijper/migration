data_1 <- read.csv("sim_migration_18_6_2019_113552_-483780392", sep = ";", header = T, row.names = NULL)
dim(data_1)

simulation_data_1 <- data_1[1:5000,]
