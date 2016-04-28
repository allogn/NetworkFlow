#!/bin/bash
DATADIR="../../data/spatial"
for f in $DATADIR/*.sgr;
do
	../Release/NetworkFlows -r 1 -i $f -a 4 -p 0 -l uniSpatial.csv
	../Release/NetworkFlows -r 1 -i $f -a 4 -p 1 -l uniSpatial.csv
done

echo "Done."
