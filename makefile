xmigration : migration.cpp auxiliary.h
	g++ -std=c++11 -Wall -O3 -ggdb migration.cpp -o xmigration

.PHONY: clean

clean:
	rm -rf xmigration
