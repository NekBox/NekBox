#!/usr/bin/env python
from sys import argv, exit
import numpy as np

files = []
files.append(open(argv[1], 'r'))
files.append(open(argv[2], 'r'))
vals = [[1]]
failed = False
while len(vals[0]) > 0:
    vals = [np.array([float(x) for x in fp.readline().split()[2:]], dtype=float) for fp in files]
    diffs = np.maximum(*vals) - np.minimum(*vals)
    print(diffs)
    if np.any(diffs > 0):
        failed = True

if failed:
    exit(6)
else:
    exit(0)
