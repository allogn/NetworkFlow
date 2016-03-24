#!/bin/bash
TESTDIR="../../data/randSSST"
UTILDIR="../../MinCostFlowLigra/utils"

rm $TESTDIR/*

SIZE=16
#7 : 1024 (2^10)
#10: 8192
#11: 16384
#14: 131 072
#17: 1M (2^20)

for mult in {1..10..1};
do
	for i in $(seq 1 5);
	do
	    echo $TESTDIR/rand_${SIZE}_${i}.adj
		$UTILDIR/randLocalGraph $SIZE $TESTDIR/rand_${SIZE}_${i}.adj
		../Debug/AdjToGr $TESTDIR/rand_${SIZE}_${i}.adj $TESTDIR/rand_${SIZE}_${i}.gr
	done
	SIZE=$(($SIZE * 2))
done
rm $TESTDIR/*.adj
echo "Generation done."
