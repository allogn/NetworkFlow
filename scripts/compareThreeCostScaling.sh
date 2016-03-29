#!/bin/bash
DATADIR="../../data/bipartite"

for f in $DATADIR/*;
do
	../Debug/NetworkFlows -r 1 -i $f -a 3 -l dijkCSonBipart.csv
	../Debug/NetworkFlows -r 1 -i $f -a 5 -l dijkCSonBipart.csv
	../Debug/NetworkFlows -r 1 -i $f -a 4 -l dijkCSonBipart.csv
done

echo "Done."
