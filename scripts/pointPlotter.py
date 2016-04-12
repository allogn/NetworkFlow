#!/usr/bin/python3

__author__ = 'alvis'
from os.path import isfile
import argparse
import numpy as np
from matplotlib import pyplot as plt

colors = ['b', 'c', 'r']

parser = argparse.ArgumentParser();
parser.add_argument("input",help="input file")
# parser.add_argument("output",help="plot file")
args = parser.parse_args()
infile = args.input
# outfile = args.output
outfile = infile + ".png"

if not isfile(infile):
    print("Input file does not exist")
    exit()

x = []
y = []
color = []
nodes = 0
with open(infile, "r") as f:
    for line in f:
        if (line[0] == "#"):
            continue
        else:
            nodes = int(line)
            break
    print("Nodes:",nodes)
    for line in f:
        data = line.split()
        x.append(float(data[0]))
        y.append(float(data[1]))
        if int(data[2]) > 0:
               color.append(colors[2])
        else:
               if (int(data[2]) < 0):
                       color.append(colors[0])
               else:
                       color.append(colors[1])
plt.scatter(x, y, c=color, s=5, linewidths=0)
plt.title("blue: neutral, green: targets, red: sources")
plt.savefig(outfile, dpi=300)
print("Saved to", outfile)
