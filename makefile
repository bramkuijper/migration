xmigration : migration.cpp auxiliary.h
	g++ -Wall -O3 -ggdb migration.cpp -o xmigration -lm -lrt -lgsl -lgslcblas

.PHONY: clean

clean:
	rm -rf xmigration
