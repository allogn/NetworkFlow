#!/bin/bash
DATADIR="../../../data/bipartite"

for f in $DATADIR/*;
do
	../../Debug/NetworkFlows -r 3 -i $f -a 0 -l SIAsimple.csv
done

echo "Done."
