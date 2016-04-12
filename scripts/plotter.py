#!/usr/bin/python3

__author__ = 'alvis'
from os.path import isfile
import argparse
import numpy as np
from matplotlib import pyplot as plt

styles = ['rs-', 'g^-', "bo-", "co-", "kx-", "y1-"]
colors = ['r', 'g', 'b', 'c', 'k', 'y']

parser = argparse.ArgumentParser();
parser.add_argument("input",help="input file")
parser.add_argument("output",help="plot file")
args = parser.parse_args()
infile = args.input
outfile = args.output

if not isfile(infile):
    print("Input file does not exist")
    exit()

print("Input parameters >>")
title = input("Title:")
y_name = input("Y:")
x_name = input("X:")
param = input("PARAM:")
where_list = []
while True:
    new_name = input("WHERE:")
    if new_name=="":
        break;
    new_value = input("EQ:")
    where_list.append([new_name,new_value])

if x_name == "":
    x_name = "Nodes"
if y_name == "":
    y_name = "Total time"
if param == "":
    param = "Algorithm"

print("Your query: PLOT",y_name,"(",x_name,") FOR EACH",param,"WHERE")
for cond in where_list:
    print(cond[0],"=",cond[1])

db = np.loadtxt(infile, delimiter=',', dtype={'names': ('expID', 'param', 'value'),
                                         'formats': ('S30', 'S30', 'S30')})

params = np.array(db['param'],dtype=str)
values = np.array(db['value'],dtype=str)
exp = np.array(db['expID'],dtype=int)
print(len(params),"objects in DB")
for cond in where_list:
    q = np.logical_and(params == cond[0],values == cond[1])
    expIDs = np.unique(exp[q])
    q = np.in1d(exp,expIDs)
    params = params[q]
    values = values[q]
    exp = exp[q]
print(len(params),"objects after filtering")

# getting all parameter values
q = params == param
vars = np.unique(values[q])
print(len(vars),"params")

plt.title(title)
plt.xlabel(x_name)
plt.ylabel(y_name)
ax = plt.gca()
ax.set_xscale("log")#, nonposx='clip')
ax.set_yscale("log")#, nonposy='clip')

# plot for each param
i = 0
for var in vars:
    q = np.logical_and(params == param, values == var)
    expIDs = np.unique(exp[q])
    if len(expIDs) == 0:
        print("Nothing found for", var)
    else:
        print(len(expIDs),"objects for", var)
        q = np.in1d(exp,expIDs)
        # assuming they go together
        # np.set_printoptions(threshold=np.nan)
        x = np.array(values[np.logical_and(params == x_name, q)], dtype=float)
        y = np.array(values[np.logical_and(params == y_name, q)], dtype=float)
        if len(x) != len(y):
            print("Wrong input data: x and y sizes differ for",var,":",len(x),"and",len(y))
            exit()
        argx = np.argsort(x)
        x = x[argx]
        y = y[argx]

        # prune zero data at the beginning
        ind = 0
        while y[ind] == 0:
            ind += 1
        x = x[ind:]
        y = y[ind:]

        y_split = np.split(y, np.where(np.diff(x))[0]+1)
        x_split = np.unique(x)
        y_mean = np.mean(y_split, 1)
        y_std = np.std(y_split, 1)
        print(y_mean)
        print(y_std)
        if i >= len(styles):
            print("Provide more styles")
            exit()
        plt.errorbar(x_split, y_mean, yerr=y_std, fmt=styles[i], label=var)
        plt.scatter(x, y, c=colors[i], s=5, linewidths=0)
        i += 1
plt.legend(loc='upper center', fontsize='8')
plt.savefig(outfile, dpi=300)
print("Saved to", outfile)
