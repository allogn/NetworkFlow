#!/bin/bash
DATADIR="../data/bipartite/"

#generate uniform weights from 0 to 1000, 10 files per size
for SIZE in {100,1000,5000};
do
	for i in {1,1,10};
	do
		../Debug/
		CURSIZE=$(($SIZE*$i))
		$UTILDIR/randLocalGraph $CURSIZE $TESTDIR/rand_${CURSIZE}_nw.adj
		$UTILDIR/adjGraphAddWeights $TESTDIR/rand_${CURSIZE}_nw.adj $TESTDIR/rand_${CURSIZE}.adj
	done
done

echo "Generation done."
