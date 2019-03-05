import numpy as np
import random
import csv
import sys
import math
from itertools import product

# current set of tuning parameters and global variables is below
step = 15
# options currently just include "EM" but will be extended later on
method = "EM"
outfile = ""
end_bits = []

# return the transposed matrix
def reorient_genos(genos_t):
	gtlen = len(genos_t)
	slen = len(genos_t[0])
	genos = [[genos_t[jj][ii] for jj in range(gtlen)] for ii in range(slen)]
	#print(len(genos))
	#print(len(genos[0]))
	return genos

# expects h_t to be the transpose of the matrix of haplotypes (i.e. h_t[0]
# is the first allele? of each haplotype
def print_haplotypes(h_t):
	#h_t = reorient_genos(hapl)
	#print(outfile)
	with(open(outfile, "a+")) as f:
		rlen = len(h_t[0])
		row_cnt = len(h_t)
		#if row_cnt == step + 1:
		#	row_cnt = step
		'''for jj in range(row_cnt):
			for ii in range(rlen - 1):
				f.write(h_t[jj][ii] + ' ')
			f.write(h_t[jj][rlen - 1] + '\n')'''
		for row in h_t:
			for ii in range(rlen - 1):
				f.write(row[ii] + ' ')
			f.write(row[rlen - 1] + '\n')

		'''
		hapl_len = len(hapl[0])
		num_hapls = len(hapl)
		for ii in range(hapl_len):
			for jj in range(num_hapls - 1):
				f.write(hapl[jj][ii] + ' ')
			f.write(hapl[num_hapls - 1][ii] + '\n')'''
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
		#print(geno)
		#print(len(haplotype_possibilities))
		haplo_temp = []
		hetero_cnt = geno.count('1')
		#print(geno)
		#print(hetero_cnt)
		if hetero_cnt == 0:
			homo_hapl = tuple(gen_hapl(geno, []))
			haplotype_possibilities.append([(homo_hapl, homo_hapl)])
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

# checks for parity between the end bits of each inferred haplotype and flips
# if necessary
def check_for_parity(hapl):
	num_hapls = len(hapl)
	elen = len(end_bits)
	#print(end_bits)
	#print([hapl[ii][0] for ii in range(num_hapls)])
	if num_hapls != elen:
		print('Houston, we have a problem. Mismatched hapl and end_bits. {} {}'.format(hapl, elen))
	else:
		for ii in range(0, num_hapls, 2):
			if end_bits[ii] != hapl[ii][0]:
				temp_lst = hapl[ii + 1]
				hapl[ii + 1] = hapl[ii]
				hapl[ii] = temp_lst
	return

# updates the probability of each pair of haplotypes
def e_step(compat_probs, hapl_probs):
	for hapl_list in compat_probs:
		# two for loops; one just computes all P_{h_1}P_{h_2} and the other
		# updates the probabilities in compat_probs with P_{h_1}P_{h_2} / \sigma
		ph_sum = 0.0
		hlen = len(hapl_list)
		for ii in range(hlen):
			h1 = hapl_list[ii][0][0]
			h2 = hapl_list[ii][0][1]
			ph1_ph2 = hapl_probs[h1] * hapl_probs[h2]
			hapl_list[ii][1] = ph1_ph2
			ph_sum += ph1_ph2
		for ii in range(hlen):
			hapl_list[ii][1] /= ph_sum
	return

# updates the probability of each unique haplotype
def m_step(compat_probs, hapl_probs):
	two_n = float(2 * len(compat_probs))
	num_uniq = len(hapl_probs)
	prob_dict = {}
	for uhapl in hapl_probs:
		prob_dict[uhapl] = 0
	# sum all the entries in compat_probs corresponding to a given haplotype
	for hapl_list in compat_probs:
		for hapl_tuple in hapl_list:
			for hapl in hapl_tuple[0]:
				prob_dict[hapl] += hapl_tuple[1]
	# normalize the probabilities (divide by 2n) and assign to hapl_probs
	for uhapl in hapl_probs:
		hapl_probs[uhapl] = prob_dict[uhapl] / two_n
	return

def set_init_hapl_probs(hapl_set, hapl_probs, em_init):
        if em_init == "rand":
                print('hello')
                ps = np.random.uniform(.0001, .9999, len(hapl_set))
                ps /= np.sum(ps)
                ii = 0
                for hapl in hapl_set:
                        hapl_probs[hapl] = ps[ii]
                        ii += 1
        # Change below line if we add more EM initialization methods
        #if em_init == "unif":
        else:
                num_hapls = len(hapl_set)
	        # intialize the probabilities for each unique haplotype
	        for hapl in hapl_set:
		        hapl_probs[hapl] = 1.0/num_hapls
        return

def q_metric(compat_probs, hapl_probs):
        q_sum = 0.0
	for hlist in compat_probs:
		max_prob = hlist[0][1]
		cur_pair = hlist[0][0]
		cnt_pairs = len(hlist)
		for ii in range(1, cnt_pairs):
			pair_prob = hlist[ii][1]
			if pair_prob > max_prob:
				max_prob = pair_prob
				cur_pair = hlist[ii][0]
		q_sum += np.log(hapl_probs[cur_pair[0]])
		q_sum += np.log(hapl_probs[cur_pair[1]])
        return q_sum
                        
# function orchestrating the meat of the EM algorithm
def process_EM(genos):
	compat_hapls = gen_all_compatible(genos)
	#print('len(compat_hapls): {}'.format(len(compat_hapls)))
	#print('len(compat_hapls[0]): {}'.format(len(compat_hapls[0])))
	hapl_set = set()
	# generate the set of all unique haplotypes
	for hapl_list in compat_hapls:
		for hapl_tuple in hapl_list:
			for hapl in hapl_tuple:
				hapl_set.add(hapl)
	hapl_probs = {}
        #hapl_probs_0 = {}
        set_init_hapl_probs(hapl_set, hapl_probs, "unif")
        #set_init_hapl_probs(hapl_set, hapl_probs_0, "rand")
        '''for h in hapl_probs:
                print(hapl_probs[h])
                break
        for h in hapl_probs_0:
                print(hapl_probs_0[h])
                break'''
        compat_probs = []
        #compat_probs_0 = []
	# initialize probabilities for each haplotype in C(g)
	for hapl_list in compat_hapls:
		aug_hapls = []
		l_len = float(len(hapl_list))
		for hapl_tuple in hapl_list:
			aug_hapls.append([hapl_tuple, 1/l_len])
		compat_probs.append(aug_hapls)
                #compat_probs_0.append(aug_hapls)
        #e_step(compat_probs_0, hapl_probs_0)
        #print(compat_probs_0[0])
        #print(compat_probs[0])
	#print(num_hapls)
	# run the EM algorithm some number of times
	for ii in range(6):
		m_step(compat_probs, hapl_probs)
		e_step(compat_probs, hapl_probs)
  		'''q0 = q_metric(compat_probs, hapl_probs)
        q1 = q_metric(compat_probs_0, hapl_probs_0)
        print(q0)
        print(q1)'''
        
        #m_step(compat_probs_0, hapl_probs_0)
		#e_step(compat_probs_0, hapl_probs_0)
		# select the haplotype pair for each genotype that maximizes the probability
        # goodness metric: sum of log of hapl probs
        #m_step(compat_probs, hapl_probs)
        #m_step(compat_probs_0, hapl_probs_0)
        #q0 = q_metric(compat_probs, hapl_probs)
        #q1 = q_metric(compat_probs_0, hapl_probs_0)
        #print(q0)
        #print(q1)
        #if q1 > q0:
        #        compat_probs = compat_probs_0
	inferred_hapls = []
	for hlist in compat_probs:
		max_prob = hlist[0][1]
		cur_pair = hlist[0][0]
		cnt_pairs = len(hlist)
		for ii in range(1, cnt_pairs):
			pair_prob = hlist[ii][1]
			if pair_prob > max_prob:
				max_prob = pair_prob
				cur_pair = hlist[ii][0]
		inferred_hapls.append(cur_pair[0])
		inferred_hapls.append(cur_pair[1])
	global end_bits
	if len(end_bits) != 0:
		check_for_parity(inferred_hapls)
	h_t = reorient_genos(inferred_hapls)
	'''if len(end_bits) == 0:
		print_haplotypes(h_t)
	else:
		print_haplotypes(h_t[1:])'''
	print_haplotypes(h_t)
	end_bits = h_t[-1]
	#print_haplotypes(h_t)
	return


def process_geno(genos):
	if method == "EM":
		process_EM(genos)
	return

def main(fname):
	print(fname)
	global outfile
	# CHANGE THIS LINE LATER!!! CURRENTLY TRYING TO AVOID OVERWRITING THE ACTUAL FILE
	outfile = fname[:-4] + "_sol_wip_0.txt"
	print(outfile)
	input_mat = []
	with(open(fname, 'rb')) as f:
		reader = csv.reader(f, delimiter=' ')
		for row in reader:
			input_mat.append(row) #[int(x) for x in row])
	idx = 0
	ilen = len(input_mat)
	row_len = len(input_mat[0])
	num_hetero = [0 for x in range(row_len)]
	hetero_hyper = 16
	skip_next = False
	for ii in range(ilen):
		if skip_next:
			skip_next = False
			continue
		max_h = num_hetero[0]
		for jj in range(row_len):
			if input_mat[ii][jj] == '1':
				nh = num_hetero[jj] + 1
				num_hetero[jj] = nh
				if nh > max_h:
					max_h = nh
		if max_h >= hetero_hyper:
			process_geno(reorient_genos(input_mat[idx : ii + 1]))
			skip_next = True
			idx = ii + 1
			num_hetero = [0 for x in range(row_len)]
			print(ii)
		elif ii == ilen - 1:
			process_geno(reorient_genos(input_mat[idx:]))
			print(ii)
	'''for ii in range(step, ilen, step):
		# orig
		#process_geno(reorient_genos(input_mat[ii - step:ii]))
		idx = ii
		print(ii)
		# checking end agreement
		if ii + 1 < ilen:
			process_geno(reorient_genos(input_mat[ii - step:ii + 1]))
		else:
			process_geno(reorient_genos(input_mat[ii - step:ii]))
	process_geno(reorient_genos(input_mat[idx:]))'''
        '''for ii in range(step, 1000, step):
                idx = ii
                print(ii)
                if ii + 1 < ilen:
                        process_geno(reorient_genos(input_mat[ii - step:ii + 1]))
		else:
			process_geno(reorient_genos(input_mat[ii - step:ii]))
	process_geno(reorient_genos(input_mat[idx:1001]))'''
	print(len(input_mat))

if __name__=="__main__":
	if len(sys.argv) != 2:
		print("Incorrect number of arguments. Provide one filename.")
	else:
		main(sys.argv[1])
