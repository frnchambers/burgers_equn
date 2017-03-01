#!/usr/bin/python3

import numpy as np
import matplotlib

import sys as sys


# def read_blocks(input_file):
#     empty_lines = 0
#     blocks = []
#     for line in open(input_file):
#         # Check for empty/commented lines
#         if not line or line.startswith('#'):
#             # If 1st one: new block
#             if empty_lines == 0:
#                 blocks.append([])
#             empty_lines += 1
#         # Non empty line: add line in current(last) block
#         else:
#             print line
#             empty_lines = 0
#             b=np.array(line.split()).astype(float)
#             blocks[-1].append(b)
#     return blocks[i:j + 1]



# s="data/solution_example.dat"

# for block in read_blocks(s, 0, 1):
#     print '-> block'
#     for line in block:
#         print line




        





solun = np.loadtxt("data/solution_example.dat", delimiter="\n\n", unpack=True)
print(solun)


# import itertools

# def block_separator(line):
#     return line=="e\n"

# with open("data/solution.dat") as f:
#     for block in itertools.groupby(f,block_separator):
#         # print(key,list(group))  # uncomment to see what itertools.groupby does.
#         u = block
#         print u
#         # if not key:
#         #     data={}
#         #     for item in group:
#         #         field,value=item.split(':')
#         #         value=value.strip()
#         #         data[field]=value
#         #     print('{FamilyN} {Name} {Age}'.format(**data))
