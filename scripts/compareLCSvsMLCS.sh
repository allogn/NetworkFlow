#!/bin/bash
DATADIR="../../data/bipartite"

for f in $DATADIR/*.gr;
do
	../Debug/NetworkFlows -r 1 -i $f -a 7 -l ASIAfullBipart.csv
done

echo "Done."
