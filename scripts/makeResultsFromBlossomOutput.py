#!/usr/bin/python3

import subprocess
import re
import os
from os import listdir
from os.path import isfile, join
import datetime
import random

dirname = "../../data/bipartite_sparse"
onlyfiles = [f for f in listdir(dirname) if isfile(join(dirname, f))]
r1 = re.compile('^.*\.bl$')
graphs = list(filter(r1.match, onlyfiles))
print("total ", len(graphs))
fres = open("blossomResultSparse.csv",'w')
for f in graphs:
    random.seed()
    expId = random.randint(10000000,99999999)
    with open(join(dirname,f)) as graph:
        for line in graph:
            if line[0] == '#':
                fres.write(str(expId) + ',' + line[1:line.find(' ')] + ',' + line[line.find(' '):])
            else:
                break
    out = subprocess.check_output(["../utils/blossom5-v2.05.src/blossom5", "-e", join(dirname,f)])
    print("Output:",out.decode())
    r = re.compile('^(.* with )([0-9]+)( nodes and )([0-9]+)( edges.*done \[)([0-9\.]+)( secs.*= )([0-9]+)(.*)$', re.DOTALL)
    r_cost = r.search(out.decode())
    nodes = r_cost.group(2)
    edges = r_cost.group(4)
    time = r_cost.group(6)
    cost = r_cost.group(8)
    fres.write(str(expId) + ',Algorithm,BlossomV\n');
    fres.write(str(expId) + ',Input file,' + f + '\n')
    fres.write(str(expId) + ',Nodes,' + str(nodes) + '\n')
    fres.write(str(expId) + ',Edges,' + str(edges) + '\n')
    fres.write(str(expId) + ',Total time,' + str(time) + '\n')
    fres.write(str(expId) + ',Total cost,' + str(cost) + '\n')
print("Done.")
