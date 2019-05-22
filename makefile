xmigration : migration.cpp auxiliary.h
	g++ -Wall -O3 migration.cpp -o xmigration -lm -lrt -lgsl -lgslcblas
