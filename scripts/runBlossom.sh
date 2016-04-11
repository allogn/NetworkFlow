#!/bin/bash
DATADIR="../../data/bipartite"

for f in $DATADIR/*.bl;
do
	../utils/blossom5-v2.05.src/blossom5 -e $f >> blossomRes.txt
done

echo "Done."
