xmigration : migration.cpp auxiliary.h
	g++ -std=c++11 -Wall -O3 -ggdb migration.cpp -o xmigration -lgsl -lgslcblas

.PHONY: clean

clean:
	rm -rf xmigration
