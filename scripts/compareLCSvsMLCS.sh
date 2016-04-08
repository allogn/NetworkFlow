#!/bin/bash
DATADIR="../../data/bipartite"

for f in $DATADIR/*.gr;
do
	../Release/NetworkFlows -r 3 -i $f -a 2 -l LDfullBipart.csv
done

echo "Done."
