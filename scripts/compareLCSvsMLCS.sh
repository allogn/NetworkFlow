#!/bin/bash
DATADIR="../../data/bipartite"

for f in $DATADIR/*.gr;
do
	../Debug/NetworkFlows -r 1 -i $f -a 0 -l SIAfullBipart.csv
done

echo "Done."
