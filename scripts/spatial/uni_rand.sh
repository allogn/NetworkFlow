#!/bin/bash
DATADIR="../../data/spatial"
BINDIR="../Debug"

SIZE=128
for mult in {1..4..1};
do
	#fully random
	SOURCENUM=1
	while [ $SOURCENUM -le $(($SIZE/2)) ];
	do
		TARGETNUM=$SOURCENUM
		while [ $TARGETNUM -le $(($SIZE/2)) ];
		do
			echo $BINDIR/SpatialGenerator -r 2 -n $SIZE -p 0 -d 0 --density 1 -s $SOURCENUM -t $TARGETNUM -o $DATADIR/${SIZE}_rand_uni_s${SOURCENUM}_t${TARGETNUM}_d1
			$BINDIR/SpatialGenerator -r 2 -n $SIZE -p 0 -d 0 --density 1 -s $SOURCENUM -t $TARGETNUM -o $DATADIR/${SIZE}_rand_uni_s${SOURCENUM}_t${TARGETNUM}_d1
			TARGETNUM=$(($TARGETNUM*2)) 
		done
		SOURCENUM=$(($SOURCENUM*2))
	done
	SIZE=$(($SIZE*2))
done
