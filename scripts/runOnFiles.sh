#!/bin/bash
DATADIR="../../data/bipartite_sparse"
for f in $DATADIR/*.gr;
do
	../Release/NetworkFlows -r 1 -i $f -a 2 -p 0 -l ../results/LDA_sparse_bipart.csv
	../Release/NetworkFlows -r 1 -i $f -a 0 -p 0 -l ../results/SIA_sparse_bipart.csv
	../Release/NetworkFlows -r 1 -i $f -a 5 -p 0 -l ../results/OLemon_sparse_bipart.csv
done

echo "Done."
