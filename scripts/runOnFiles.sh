#!/bin/bash
DATADIR="../../data/graphs"

for f in $DATADIR/*.gr;
do
	../Release/NetworkFlows -r 10 -i $f -a 3 -p 0 -l rmat.csv
	../Release/NetworkFlows -r 10 -i $f -a 3 -p 1 -l rmat.csv
	../Release/NetworkFlows -r 10 -i $f -a 3 -p 2 -l rmat.csv
done

echo "Done."
