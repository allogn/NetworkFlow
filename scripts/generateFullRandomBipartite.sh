#!/bin/bash
DATADIR="../../data/bipartite/"

#generate uniform weights from 0 to 1000, 10 files per size
for SIZE in {6,1000,2000};
do
	for i in {1,1,10};
	do
		NAME="${DATADIR}bi_full_${SIZE}_0_1000_uni.gr"
		../Debug/Generator -g 0 -n $SIZE -l 0 -u 1000 -d 0 -o $NAME
	done
done

echo "Generation done."
