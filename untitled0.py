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

def reorient_genos(genos_t):
	gtlen = len(genos_t)
	slen = len(genos_t[0])
	genos = [[genos_t[jj][ii] for jj in range(gtlen)] for ii in range(slen)]
	print(len(genos))
	print(len(genos[0]))
	return genos

def print_haplotypes(hapl):
	h_t = reorient_genos(hapl)
	with(open(outfile, "a+")) as f:
		rlen = len(h_t[0])
		for row in h_t:
			for ii in range(rlen - 1):
				f.write(row[ii] + ' ')
			f.write(row[rlen - 1] + '\n')
	return

def gen_all_compatible(genos):
	haplotype_possibilities = []
	for geno in genos:
		haplo_temp = []
		hetero_cnt = geno.count('1')
		haps = set()

		print(geno)
		print(hetero_cnt)
		haplotype_possibilities.append(haplo_temp)

def process_EM(genos):
	gen_all_compatible(genos)

def process_geno(genos):
	if method == "EM":
		process_EM(genos)
	return

def main(fname):
	print(fname)
	outfile = fname[:-3] + "_sol.txt"
	input_mat = []
	with(open(fname, 'rb')) as f:
		reader = csv.reader(f, delimiter=' ')
		for row in reader:
			input_mat.append(row) #[int(x) for x in row])
	idx = 0
	ilen = len(input_mat)

	for ii in range(step, ilen, step):
		# uncomment this in a bit!!!
		#process_geno(reorient_genos(input_mat[ii - step:ii]))
		idx = ii
	# need to finish off the remainder here
	process_geno(reorient_genos(input_mat[idx:]))
	print(len(input_mat))

if __name__=="__main__":
	if len(sys.argv) != 2:
		print("Incorrect number of arguments. Provide one filename.")
	else:
		main(sys.argv[1])