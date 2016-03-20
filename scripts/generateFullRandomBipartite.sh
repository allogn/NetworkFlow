#!/bin/bash
DATADIR="../../data/NetworkFlowTests/bipartite/"

#generate uniform weights from 0 to 1000, 10 files per size
for SIZE in {100..1000..100};
do
	for i in {1..3..1};
	do
		NAME="${DATADIR}bi_full_${SIZE}_1_10_uni_${i}.gr"
		NAME2="${DATADIR}bi_full_${SIZE}_1_10_uni_${i}.bl"
		../Debug/Generator -g 0 -n $SIZE -l 0 -u 100 -d 0 -o ${NAME} --blossom ${NAME2}
	done
	for i in {1..3..1};
	do
		NAME="${DATADIR}bi_full_${SIZE}_30_10_gauss_${i}.gr"
		NAME2="${DATADIR}bi_full_${SIZE}_30_10_gauss_${i}.bl"
		../Debug/Generator -g 0 -n $SIZE -l 30 -u 10 -d 1 -o $NAME --blossom $NAME2
	done
	for i in {1..3..1};
	do
		NAME="${DATADIR}bi_full_${SIZE}_1_exp_${i}.gr"
		NAME2="${DATADIR}bi_full_${SIZE}_1_exp_${i}.bl"
		../Debug/Generator -g 0 -n $SIZE -l 1 -u 1 -d 2 -o $NAME --blossom $NAME2
	done
done

echo "Generation done."
