#!/bin/bash
DATADIR="../../data/randSSST"

for f in $DATADIR/*;
do
	../Debug/NetworkFlows -i $f -a 3 -l LCSvsMLCSonRandSSST.csv
	../Debug/NetworkFlows -i $f -a 5 -l LCSvsMLCSonRandSSST.csv
done

echo "Done."
