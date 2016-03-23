#!/bin/bash
DATADIR="../../data/bipartite/"

SIZE=16
#7 : 1024 (2^10)
#10: 8192
#11: 16384
#14: 131 072
#17: 1M (2^20)
for mult in {1..7..1};
do
	for i in {1..3..1};
	do
		NAME="${DATADIR}bi_full_${SIZE}_1_100_uni_${i}.gr"
		NAME2="${DATADIR}bi_full_${SIZE}_1_100_uni_${i}.bl"
		../Debug/Generator -g 0 -n $SIZE -l 1 -u 100 -d 0 -o ${NAME} --blossom ${NAME2}
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
	SIZE=$(($SIZE*2))
done

echo "Generation done."
