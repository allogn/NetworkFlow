#!/bin/bash
DATADIR="../../data/bipartite_sparse"
for f in $DATADIR/*.gr;
do
	../Release/NetworkFlows -r 1 -i $f -a 5 -p 0 -l OLemon_sparse_bipart.csv
done

echo "Done."
