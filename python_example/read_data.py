#!/usr/bin/python

import numpy as np

## Works fine:
solun_1 = np.loadtxt("solution_single.dat", unpack=True)
print(solun_1)

## Works not so fine:
solun_2 = np.loadtxt("solution_example.dat", delimiter="-----", unpack=True)
print(solun_2)

