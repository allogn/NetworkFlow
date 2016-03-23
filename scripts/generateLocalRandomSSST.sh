#!/bin/bash
TESTDIR="../../data/randSSST"
UTILDIR="../../MinCostFlowLigra/utils"
rm $TESTDIR/*
for SIZE in {100..1000..200};
do
	for i in $(seq 1 5);
	do
	    echo $TESTDIR/rand_${SIZE}_${i}.adj
		$UTILDIR/randLocalGraph $SIZE $TESTDIR/rand_${SIZE}_${i}.adj
		../Debug/AdjToGr $TESTDIR/rand_${SIZE}_${i}.adj $TESTDIR/rand_${SIZE}_${i}.gr
	done
done
rm $TESTDIR/*.adj
echo "Generation done."
