#!/bin/bash
# Checks that the new migration.cpp, with all new switches off, gives exactly
# the same output as the old one when both are run with the same seed.
#
# Usage (from the repository folder):
#     git show revision:migration.cpp > migration_old.cpp
#     bash regression_test.sh migration_old.cpp migration.cpp
#
# Both versions are compiled for a short run (40 years). The old version is
# given the same optional seed argument as the new one, so both start from
# the same random numbers.

set -e
OLD=$1
NEW=$2
YEARS=40
SEED=12345

mkdir -p regtest
cp auxiliary.h regtest/

perl -pe "s/long int number_years = 250000;/long int number_years = $YEARS;/; s/^(\s*filename_output = argv\[33\];)/\$1\n    if (argc > 34) { seed = std::strtoul(argv[34], nullptr, 10); rng_r.seed(seed); }/" "$OLD" > regtest/old.cpp
perl -pe "s/long int number_years = 250000;/long int number_years = $YEARS;/" "$NEW" > regtest/new.cpp

cd regtest
g++ -std=c++11 -O3 old.cpp -o old
g++ -std=c++11 -O3 new.cpp -o new

# Table A1 parameters, followed by S, F, post-equilibrium run time, K decline, harvest, files
P="60 1 0.3 25 0 0 0.1 0.5 0.30103 0.04 0.02 0.5 30 0 0.05 0.1 20 20 2 10000 60 5 1 0 5"

for S in 0 1; do
    ./old $P $S 8 0 1 0 none none old_S$S $SEED > /dev/null
    ./new $P $S 8 0 1 0 none none new_S$S $SEED > /dev/null
    if cmp -s old_S$S new_S$S && cmp -s old_S${S}_dist new_S${S}_dist; then
        echo "S = $S: outputs identical"
    else
        echo "S = $S: OUTPUTS DIFFER"
    fi
done
