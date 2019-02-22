import numpy as np
import random
import csv
import sys
import math
from itertools import product

outfile = ""

# return the transposed matrix
def reorient_genos(genos_t):
	gtlen = len(genos_t)
	slen = len(genos_t[0])
	genos = [[genos_t[jj][ii] for jj in range(gtlen)] for ii in range(slen)]
	return genos

def rand_row(input_row, repeat_cnt):
	r_row = []
	for allele in input_row:
		if allele == '1':
			r_row += ['0' if ch == 0 else '1' for ch in np.random.binomial(1, 0.5, repeat_cnt)]
		elif allele == '0':
			r_row += ['0' for ii in range(repeat_cnt)]
		else:
			r_row += ['1' for ii in range(repeat_cnt)]
	return r_row

def rand_compat_hapls(input_mat, repeat_cnt):
	ilen = len(input_mat)
	return [rand_row(input_mat[ii], repeat_cnt) for ii in range(ilen)]

# generates the next level of the localized haplotype cluster model
def split_l(l_graph_row, hapls_row):
	node_map = {}
	n_cnt = 0
	for node in l_graph_row:
		for index in node[3]:
			hval = hapls_row[index]
			if (n_cnt, hval) not in node_map:
				node_map[(n_cnt, hval)] = {index}
			else:
				node_map[(n_cnt, hval)].add(index)
		n_cnt += 1
	new_row = []
	for k in node_map:
		new_row.append([[k[0]], [k[1]], [len(node_map[k])], node_map[k]])
	return new_row

# combines all fields of two nodes
def merge_two(n1, n2):
	merged = [[n1[ii] + n2[ii]] for ii in range(2)]
	merged.append(n1[2] + n2[2])
	merged.append(n1[3].union(n2[3]))
	return merged

# merges all the nodes and returns a list containing that single node
def merge_all(hapl_dict):
	while len(hapl_dict) != 1:
		indices = hapl_dict.keys()
		hapl_dict[indices[0]] = merge_two(hapl_dict[indices[0]], hapl_dict[indices[1]])
		hapl_dict.remove(indices[1])
	return [hapl_dict[hapl_dict.keys()[0]]]

# compute pairwise threshold
def thresh(n1, n2):
	return math.sqrt((1.0/sum(n1[2])) + (1.0/sum(n2[2])))

# get all pairs and compute the threshold. find the minimum and we're chilling
def min_thresh(hapl_dict):
	indices = hapl_dict.keys()
	ilen = len(indices)
	min_t = thresh(hapl_dict[indices[0]], hapl_dict[indices[1]])
	for ii in range(ilen):
		for jj in range(ii + 1, ilen):
			t = thresh(hapl_dict[indices[ii]], hapl_dict[indices[jj]])
			if t < min_t:
				min_t = t
	return min_t

def merge_score(l1, idx2, hapls, ii):

# attempts to merge nodes that contain similar information
# idx: the next haplotype to look at. If no more, we can just merge all of them!
def merge_l(l_row, hapls, idx):
	l_len = len(l_row)
	# nothing to merge, so return
	if l_len == 1:
		return l_row
	# create a dictionary mapping indices to nodes. easier to combine.
	hapl_dict = {}
	for ii in range(l_len):
		hapl_dict[ii] = l_row[ii]
	# merge all of them if no more rows of haplotype data
	if idx >= len(hapls):
		return merge_all(hapl_dict)
	# compute minimum threshold
	min_t = min_thresh(hapl_dict)
	# compute merging score for each pair of nodes
	# find min and if it is less than corresponding threshold, accept
	# combine the two nodes and repeat the process until no satisfying pair
	# convert back from dictionary to list (i.e. dict.values())

# generates the localized haplotype cluster model
def localized_hapl_cluster(hapls):
	# list of lists of lists corresp. to each level
	# each low-level list contains a list of adjacent nodes, corresponding edge
	# labels, edge counts, and the set of indices of haplotypes that arrive at 
	# a given node
	l_graph = []
	hlen = len(hapls)
	r = range(0, len(hapls[0]))
	index_set = set()
	for ii in r:
		index_set.add(ii)
	l_graph.append([[[-1], ['-1'], [], index_set]])
	for ii in range(hlen):
		# use idx rather than ii because l_graph has one level before 1st iter
		idx = ii + 1
		# split (i.e., create next level)
		l_row = split_l(l_graph[-1], hapls[ii])
		#for ele in l_graph:
		#	print(ele)
		# merge (combine nodes whose futures are sufficiently similar)
		l_graph.append(merge_l(l_row, hapls, idx))
	# afterwards, just add a node that all of the previous label maps to
	return l_graph

# also have function to sample the haplotypes
def main(fname):
	print(fname)
	global outfile
	# CHANGE THIS LINE LATER!!! CURRENTLY TRYING TO AVOID OVERWRITING THE ACTUAL FILE
	outfile = fname[:-4] + "_sol_wip.txt"
	print(outfile)
	input_mat = []
	with(open(fname, 'rb')) as f:
		reader = csv.reader(f, delimiter=' ')
		for row in reader:
			input_mat.append(row)
	# for each genotype, randomly generate some number (currently 10) compatible haplotypes
	hapls = rand_compat_hapls(input_mat, 10)
	local_graph = []
	dipl_hmm = []
	print('Constructed random haplotypes')
	for ii in range(10):
		local_graph = localized_hapl_cluster(hapls)
		# dipl_hmm = construct_dipl(local_graph)
		# hapls = generate_new_samples(dipl_hmm, input_mat)
		print('Completed iteration {}'.format(ii))
	print('Running Viterbi algorithm on final HMM')
	#hapls = dipl_viterbi(dipl_hmm, input_mat)
	#then, print hapls
	#process_geno(reorient_genos(input_mat[idx:]))

if __name__=="__main__":
	if len(sys.argv) != 2:
		print("Incorrect number of arguments. Provide one filename.")
	else:
		main(sys.argv[1])
