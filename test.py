import numpy as np
import random
import csv
import sys
import math
from itertools import product

# current set of tuning parameters and global variables is below
step = 20
# options currently just include "EM" but will be extended later on
method = "EM"
outfile = ""
haplo_temp = []
hetero_cnt = geno.count('1')
haps = set()
permute_list = ['1' for ii in range(hetero_cnt)] + ['0' for ii in range(hetero_cnt)]
for p in product(permute_list, repeat=(2*hetero_cnt)):
	if p 
print()
print(geno)
print(hetero_cnt)

# expects a haplotype in the first of a list, i.e. ['1', '0', '1']
# returns the complementary haplotype, ['0', '1', '0']
def compl(haplotype):
	return ['0' if h == '1' else '1' for h in haplotype]
