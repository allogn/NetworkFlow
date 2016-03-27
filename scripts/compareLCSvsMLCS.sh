#!/bin/bash
DATADIR="../../data/bipartite"

for f in $DATADIR/*;
do
	../Debug/NetworkFlows -r 3 -i $f -a 3 -l LCSvsMLCSonBipart.csv
	../Debug/NetworkFlows -r 3 -i $f -a 5 -l LCSvsMLCSonBipart.csv
done

echo "Done."
