#!/bin/bash
DATADIR="../../data/bipartite_sparse"
for f in $DATADIR/*.gr;
do
	../Release/NetworkFlows -r 1 -i $f -a 2 -p 0 -l LDA_sparse_bipart.csv
	../Release/NetworkFlows -r 1 -i $f -a 0 -p 0 -l SIA_sparse_bipart.csv
done

echo "Done."
