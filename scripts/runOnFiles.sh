#!/bin/bash
DATADIR="../../data/spatial"

for f in $DATADIR/*;
do
	../Debug/NetworkFlows -r 1 -i $f -a 5 -l testSpatial.csv
done

echo "Done."
