#!/bin/bash
DATADIR="../../data/bipartite_sparse/"

SIZE=0
#7 : 1024 (2^10)
#10: 8192
#11: 16384
#14: 131 072
#17: 1M (2^20)
for mult in {1..8..1};
do
	for i in {1..2..1};
	do
		NAME="${DATADIR}bi_full_4096_m${SIZE}_1_100_uni_${i}.gr"
		NAME2="${DATADIR}bi_full_4096_m${SIZE}_1_100_uni_${i}.bl"
		../Release/Generator -g 0 -n 4096 -l 1 -u 100 --missed $SIZE -d 0 -o ${NAME} --blossom ${NAME2}
	done
	for i in {1..2..1};
	do
		NAME="${DATADIR}bi_full_4096_m${SIZE}_30_10_gauss_${i}.gr"
		NAME2="${DATADIR}bi_full_4096_m${SIZE}_30_10_gauss_${i}.bl"
		../Release/Generator -g 0 -n 4096 -l 30 -u 10 -d 1 -o $NAME --blossom $NAME2
	done
	for i in {1..2..1};
	do
		NAME="${DATADIR}bi_full_4096_m${SIZE}_1_exp_${i}.gr"
		NAME2="${DATADIR}bi_full_4096_m${SIZE}_1_exp_${i}.bl"
		../Release/Generator -g 0 -n 4096 -l 1 -u 1 -d 2 -o $NAME --blossom $NAME2
	done
	SIZE=$(($SIZE+256))
done

echo "Generation done."
