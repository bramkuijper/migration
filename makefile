xmigration : migration.cpp auxiliary.h
	g++ -Wall -O3 -ggdb migration.cpp -o xmigration -lgsl -lgslcblas

.PHONY: clean

clean:
	rm -rf xmigration
