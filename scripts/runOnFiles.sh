#!/bin/bash
DATADIR="../../data/spatial/clustManyTargets"
TH=9
for f in $DATADIR/*.sgr;
do
	../Release/NetworkFlows -r 1 -i $f -a 5 -p 0 -d 0.$TH -l rmat2.csv
	../Release/NetworkFlows -r 1 -i $f -a 3 -p 1 -d 0.$TH -l rmat.csv
	../Release/NetworkFlows -r 1 -i $f -a 3 -p 0 -d 0.$TH -l rmat.csv
	../Release/NetworkFlows -r 1 -i $f -a 3 -p 2 -d 0.$TH -l rmat.csv
done

echo "Done."
