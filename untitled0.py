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

# return the transposed matrix
def reorient_genos(genos_t):
	gtlen = len(genos_t)
	slen = len(genos_t[0])
	genos = [[genos_t[jj][ii] for jj in range(gtlen)] for ii in range(slen)]
	print(len(genos))
	print(len(genos[0]))
	return genos

# TODO: THIS ISN'T ACCEPTING THE RIGHT SORT OF INPUT...
# expects hapl to be a matrix consisting of pairs of haplotypes
#
def print_haplotypes(hapl):
	h_t = reorient_genos(hapl)
	with(open(outfile, "a+")) as f:
		rlen = len(h_t[0])
		for row in h_t:
			for ii in range(rlen - 1):
				f.write(row[ii] + ' ')
			f.write(row[rlen - 1] + '\n')
	return

# expects a haplotype in the first of a list, i.e. ['1', '0', '1']
# returns the complementary haplotype, ['0', '1', '0']
def compl(haplotype):
	return ['0' if h == '1' else '1' for h in haplotype]

# takes as input the genotype, a list, like ['1', '0', '0', '2', '1', '1'] and a
# permutation like ['0', '0', '1'], of length equal to the number of heterozygous
# alleles in the genotype and returns the corresponding haplotype in which
# 0's remain, 1's are mapped to 0's or 1's depending on the corresponding entry
# in the permutation and 2's are mapped to 1's
def gen_hapl(geno, perm):
	pidx = 0
	glen = len(geno)
	hapl = []
	for ii in range(glen):
		if geno[ii] == '0':
			hapl.append('0')
		elif geno[ii] == '2':
			hapl.append('1')
		else:
			hapl.append(perm[pidx])
			pidx += 1
	return hapl

# for each genotype in genos, generates the set of all compatible haplotype
def gen_all_compatible(genos):
	haplotype_possibilities = []
	#print(genos)
	for geno in genos:
		haplo_temp = []
		hetero_cnt = geno.count('1')
		#print(geno)
		#print(hetero_cnt)
		if hetero_cnt == 0:
			homo_hapl = tuple(gen_hapl(geno, []))
			haplotype_possibilities.append((homo_hapl, homo_hapl))
			continue
		for p in product(['0','1'], repeat=hetero_cnt):
			# because we've been constructing the complement as well, by the time this
			# condition is satisfied, we will have examined all unique pairs of haplotypes
			if p[0] == "1":
				break
			hapl = tuple(gen_hapl(geno, p))
			c_hapl = tuple(gen_hapl(geno, compl(p)))
			haplo_temp.append((hapl, c_hapl))
			haplotype_possibilities.append(haplo_temp)
	return haplotype_possibilities

# updates the probability of each pair of haplotypes
def e-step(compat_probs, hapl_probs):
	for hapl_list in compat_probs:
		pass

# updates the probability of each unique haplotype
def m-step(compat_probs, hapl_probs):


# function orchestrating the meat of the EM algorithm
def process_EM(genos):
	compat_hapls = gen_all_compatible(genos)
	hapl_set = set()
	# generate the set of all unique haplotypes
	for hapl_list in compat_hapls:
		for hapl_tuple in hapl_list:
			for hapl in hapl_tuple:
				hapl_set.add(hapl)
	hapl_probs = {}
	num_hapls = len(hapl_set)
	# intialize the probabilities for each unique haplotype
	for hapl in hapl_set:
		hapl_probs[hapl] = 1.0/num_hapls
	compat_probs = []
	# initialize probabilities for each haplotype in C(g)
	for hapl_list in compat_hapls:
		aug_hapls = []
		l_len = len(hapl_list)
		for hapl_tuple in hapl_list:
			aug_hapls.append([hapl_tuple, 1/l_len])
		compat_probs.append(aug_hapls)
	print(num_hapls)
	# run the EM algorithm some number of times
	for ii in range(10):
		m-step(compat_probs, hapl_probs)
		e-step(compat_probs, hapl_probs)
	# select the haplotype pair for each genotype that maximizes the probability
	glen = len(compat_hapls)
	for ii in range()
	return


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