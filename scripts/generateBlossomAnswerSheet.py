#!/usr/bin/python3

import subprocess
import re
import os
from os import listdir
from os.path import isfile, join
import datetime

dirname = "../../data/tests/bipartite"
onlyfiles = [f for f in listdir(dirname) if isfile(join(dirname, f))]
r1 = re.compile('^.*\.bl$')
graphs = list(filter(r1.match, onlyfiles))
print("total ", len(graphs))
for f in graphs:
    out = subprocess.check_output(["../../AlgorithmComparisonExperiments/blossom5-v2.05.src/blossom5", "-e", join(dirname,f)])
    r = re.compile('^(.*= )([0-9]+)(.*)$', re.DOTALL)
    r_cost = r.search(out.decode())
    cost = r_cost.group(2)
    
    with open("AnswerSheet.txt","a") as answers:
        cont = f + " " + str(cost) + "\n"
        answers.write(cont)
print("Done.")
