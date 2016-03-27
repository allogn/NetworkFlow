#!/bin/bash
DATADIR="../../data/bipartite"

for f in $DATADIR/*;
do
	../Debug/NetworkFlows -r 3 -i $f -a 3 -l CSonBipart.csv
	../Debug/NetworkFlows -r 3 -i $f -a 5 -l CSonBipart.csv
	../Debug/NetworkFlows -r 3 -i $f -a 4 -l CSonBipart.csv
done

echo "Done."
