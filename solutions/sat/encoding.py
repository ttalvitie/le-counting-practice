#!/usr/bin/env python3

import tempfile
import os
import os.path
import sys
import subprocess
import time

pj = os.path.join

scriptdir = os.path.dirname(os.path.abspath(__file__))

filename = sys.argv[1]
with open(filename) as fp:
	data = fp.read()
data = data.strip().split()

n = 0
while n**2 < len(data):
	n += 1
assert n**2 == len(data)

poset = [[int(data[n * i + j]) for j in range(n)] for i in range(n)]

for k in range(n):
	for i in range(n):
		for j in range(n):
			if poset[i][k] and poset[k][j]:
				poset[i][j] = 1

formula = []
var_count = None

encoding = sys.argv[2]

if encoding == "1":
	var_count = 0
	var = {}
	for i in range(n):
		for j in range(i + 1, n):
			var_count += 1
			var[(i, j)] = var_count
			var[(j, i)] = -var_count
	
	for k in range(n):
		for i in range(n):
			for j in range(n):
				if i == k or j == k or i == j:
					continue
				formula.append([-var[(i, k)], -var[(k, j)], var[(i, j)]])
	
	for i in range(n):
		for j in range(n):
			if poset[i][j]:
				formula.append([var[(i, j)]])
elif encoding == "2":
	var_count = n**2
	for i in range(n):
		clause = []
		for j in range(n):
			clause.append(1 + (n * i + j))
		formula.append(clause)

	for i in range(n):
		for j1 in range(n):
			for j2 in range(j1 + 1, n):
				formula.append([-1 - (n * i + j1), -1 - (n * i + j2)])

	for j in range(n):
		clause = []
		for i in range(n):
			clause.append(1 + (n * i + j))
		formula.append(clause)

	for j in range(n):
		clause = []
		for i1 in range(n):
			for i2 in range(i1 + 1, n):
				formula.append([-1 - (n * i1 + j), -1 - (n * i2 + j)])

	for a in range(n):
		for b in range(n):
			if not poset[a][b]:
				continue
			
			ok = True
			for k in range(n):
				if poset[a][k] and poset[k][b]:
					ok = False
					break
			if not ok:
				continue
			
			for p1 in range(n):
				for p2 in range(p1 + 1, n):
					formula.append([-1 - (n * a + p2), -1 - (n * b + p1)])
else:
	sys.stderr.write("Unknown encoding\n")
	exit(1)

print("p cnf {} {}".format(var_count, len(formula)))
for clause in formula:
    for x in clause:
        assert x != 0 and abs(x) <= var_count
    print(" ".join(list(map(str, clause)) + ["0"]))
