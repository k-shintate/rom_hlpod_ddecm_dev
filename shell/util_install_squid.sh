#!/bin/bash

cd utils/load_balancing
make -f Makefile_squid clean
make -f Makefile_squid

mkdir -p ./../bin

cp merge_graph ./../bin
cp merge_graph_bc ./../bin
rm merge_graph
rm merge_graph_bc

cd ./../..
