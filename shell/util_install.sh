#!/bin/bash

cd utils/load_balancing
make clean
make

mkdir -p ./../bin

mv merge_graph ./../bin
mv elem2graph ./../bin
mv merge_dist_val ./../bin

cd ./../..
