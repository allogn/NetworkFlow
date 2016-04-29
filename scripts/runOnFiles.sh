#!/bin/bash
DATADIR="../../data/bipartite"
for f in $DATADIR/*.gr;
do
	../Release/NetworkFlows -r 1 -i $f -a 0 -p 0 -l ../results/SIAvsSIAnew_on_bipart.csv
	../Release/NetworkFlows -r 1 -i $f -a 0 -p 1 -l ../results/SIAvsSIAnew_on_bipart.csv
done

echo "Done."
