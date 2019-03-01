import numpy as np
import random
import csv
import sys
import math
import itertools
from collections import defaultdict

dups = 3
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
	merged = [n1[ii] + n2[ii] for ii in range(2)]
	merged.append(n1[2] + n2[2])
	merged.append(n1[3].union(n2[3]))
	return merged

# merges all the nodes and returns a list containing that single node
def merge_all(hapl_dict):
	while len(hapl_dict) != 1:
		indices = hapl_dict.keys()
		hapl_dict[indices[0]] = merge_two(hapl_dict[indices[0]], hapl_dict[indices[1]])
		del hapl_dict[indices[1]]
		#hapl_dict.remove(indices[1])
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

# compute the merge score between a pair of nodes, n1 and n2, stopping 
def merge_score(n1, n2, hapls, idx2, min_diff, hapl_len):
	# compute the threshold. if ever exceeded, return 100.0 or something
	max_diff = 0.0
	thr_diff = thresh(n1, n2)
	#print(thr_diff)
	o1 = float(len(n1[3]))
	o2 = float(len(n2[3]))
	# for every level before the end of hapls, check if the difference is
	# too large: i.e. create a dictionary, which initially consists of a
	# blank string. then, take the allele at the given index (and incr)
	# and construct the next level. once all have been added to the next level,
	# purge the ones that involve too few haplotypes to possibly reach the
	# minimum threshold. at that same time, compare the difference and check
	# if it exceeds the threshold. if it doesn't move to the next pair of
	# haplotypes. terminate when all haplotypes have been exhausted. return
	# the maximum observed difference over the courseo of this.
	merge_dict_1 = {"": n1[3]}
	merge_dict_2 = {"": n2[3]}
	# windowing hyperparameter: this might really be a mistake. check later.
	depth = 0
	#print('Before while loop')
	while idx2 < hapl_len and (len(merge_dict_1) != 0 or len(merge_dict_2) != 0) and depth < 50:
		depth += 1
		post_merge_1 = {}
		post_merge_2 = {}
		#print('Construct next level')
		# construct the next level
		for h_str in merge_dict_1:
			for ii in merge_dict_1[h_str]:
				ch = hapls[idx2][ii]
				new_str = h_str + ch
				if new_str not in post_merge_1:
					post_merge_1[new_str] = {ii}
				else:
					post_merge_1[new_str].add(ii)
		for h_str in merge_dict_2:
			for ii in merge_dict_2[h_str]:
				ch = hapls[idx2][ii]
				new_str = h_str + ch
				if new_str not in post_merge_2:
					post_merge_2[new_str] = {ii}
				else:
					post_merge_2[new_str].add(ii)
		merge_dict_1 = post_merge_1
		merge_dict_2 = post_merge_2
		poss_hapls = set()
		for h_str in merge_dict_1:
			poss_hapls.add(h_str)
		for h_str in merge_dict_2:
			poss_hapls.add(h_str)
		'''print('Possible haplotypes:')
		print(poss_hapls)'''
		# check each edge 
		for h in poss_hapls:
			h1_cnt = 0
			h2_cnt = 0
			if h in merge_dict_1:
				h1_cnt = len(merge_dict_1[h])
			if h in merge_dict_2:
				h2_cnt = len(merge_dict_2[h])
			# should also just compute the difference first and also check the
			# below condition
			'''print(o2)
			print(h2_cnt)'''
			c = [h1_cnt/o1, h2_cnt/o2]
			#print(c)
			diff = abs(c[0] - c[1])
			#print(diff)
			#print('c and diff above')
			if diff >= thr_diff:
				return 100.0
			if diff >= max_diff:
				max_diff = diff
			if max(c) <= min_diff:
				if h in merge_dict_1:
					del merge_dict_1[h]
				if h in merge_dict_2:
					del merge_dict_2[h]
		idx2 += 1
	return max_diff


# attempts to merge nodes that contain similar information
# idx: the next haplotype to look at. If no more, we can just merge all of them!
def merge_l(l_row, hapls, idx):
	l_len = len(l_row)
	hapl_len = len(hapls)
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
	pmap = {}
	# compute merging score for each pair of nodes
	# find min and if it is less than corresponding threshold, accept
	# combine the two nodes and repeat the process until no satisfying pair
	min_merge_score = 99.0
	while min_merge_score < 100.0 and len(hapl_dict) > 1:
		min_merge_score = 100.0
		pairs = list(itertools.combinations(hapl_dict.keys(), 2))
		#print(pairs)
		mpair = pairs[0]
		for p in pairs:
			# no need to recalculate for pairs that didn't change (i.e. if both are unchanged)
			sc = 100.0
			if p in pmap and pmap[p][0] == 0:
				sc = pmap[p][1]
			else:
				sc = merge_score(hapl_dict[p[0]], hapl_dict[p[1]], hapls, idx, min_t, hapl_len)
				pmap[p] = [0, sc]
			if sc < min_merge_score:
				min_merge_score = sc
				mpair = p
		'''print(min_merge_score)
		print(mpair)'''
		if min_merge_score < 100.0:
			hapl_dict[mpair[0]] = merge_two(hapl_dict[mpair[0]], hapl_dict[mpair[1]])
			del hapl_dict[mpair[1]]
			pmap[mpair] = [-1, 100.0]
			for p in pmap:
				if p[0] == mpair[0] or p[0] == mpair[1] or p[1] == mpair[0] or p[1] == mpair[1]:
					pmap[p] = [-1, 100.0]
	# convert back from dictionary to list
	return hapl_dict.values()

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
		if ii % 100 == 0:
			print(ii)
		# use idx rather than ii because l_graph has one level before 1st iter
		idx = ii + 1
		# split (i.e., create next level)
		l_row = split_l(l_graph[-1], hapls[ii])
		#for ele in l_graph:
		#	print(ele)
		# merge (combine nodes whose futures are sufficiently similar)
		rows = merge_l(l_row, hapls, idx)
		print(len(rows))
		l_graph.append(rows)
		'''for ele in l_graph:
			print(ele)
			print(len(ele))'''
		#sys.exit(-1)
	# afterwards, just add a node that all of the previous label maps to
	return l_graph

# so first construct a haplotype HMM and then see if a diplotype one is necessary
def construct_dipl(local_graph):
    l_len = len(local_graph)
    tot = 0
    r1 = 0
    pis = []
    transitions = []
    t_prev_map = []
    lnodes = len(local_graph[1])
    for ii in range(lnodes):
        incoming = len(local_graph[1][ii][0])
        for jj in range(incoming):
            ecnt = local_graph[1][ii][2][jj]
            pis.append(float(ecnt))
    sump = sum(pis)
    pis = [pis[0]/sump, pis[1]/sump]
    for ii in range(2, l_len):
        rlen = len(local_graph[ii])
        emap = {}
        edges = defaultdict(list)
        for jj in range(rlen):
            incoming = len(local_graph[ii][jj][0])
            for kk in range(incoming):
                prev_node = local_graph[ii][jj][0][kk]
                sym = local_graph[ii][jj][1][kk]
                ecnt = local_graph[ii][jj][2][kk]
                emap[(prev_node, jj, sym)] = float(ecnt)
                edges[prev_node].append((prev_node, jj, sym))
        t_prev_map.append(edges)
        for e in edges:
            outgoing = len(edges[e])
            tsum = sum([emap[edges[e][ii]] for ii in range(outgoing)])
            for ii in range(outgoing):
                emap[edges[e][ii]] /= tsum
        transitions.append(emap)
    return pis, transitions, t_prev_map

# generates new phased haplotypes given the haplotype hmm and the input matrix
def generate_new_samples(hapl_hmm, pis, t_prev, input_mat):
	# construct pairs first and then split them up to match the format
	hlen = len(hapl_hmm)
	ilen = len(input_mat[0])
	# initialize based on the pis and compatibility
	# then, randomly choose the next state based on the probability of the
	# transition (which incorporates both the haplotypes)
	# map pairs of current vertices to next vertices (that satisfy genotype)
	# (cur_0, cur_1, next_geno) -> [[next_0, next_1, t_prob]
	poss_next_probs = {}
	for ii in 
	for ii in range()

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
	hapls = rand_compat_hapls(input_mat, dups) #10)
	local_graph = []
	hapl_hmm = []
	print('Constructed random haplotypes')
	for ii in range(10):
		local_graph = localized_hapl_cluster(hapls)
		#lens = []
		#for ele in local_graph:
		#	lens.append(len(ele))
		#print(lens)
		'''print('-------')
		print(len(local_graph))'''
		# I'm thinking I might be able to get away with just constructing a hapl_hmm
		# and then just keep track of pairs?
		pis, hapl_hmm, t_prev = construct_dipl(local_graph)
		'''for row in hapl_hmm:
			print(row)
		print(len(hapl_hmm))'''
		sys.exit(-1)
		# hapls = generate_new_samples(hapl_hmm, pis, t_prev, input_mat)
		# remember to reverse the order with every other iteration
		print('Completed iteration {}'.format(ii))
	print('Running Viterbi algorithm on final HMM')
	#hapls = dipl_viterbi(hapl_hmm, input_mat)
	#then, print hapls
	#process_geno(reorient_genos(input_mat[idx:]))

if __name__=="__main__":
	if len(sys.argv) != 2:
		print("Incorrect number of arguments. Provide one filename.")
	else:
		main(sys.argv[1])
