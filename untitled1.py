import numpy as np
import random
import csv
import sys
import math
import itertools
from collections import defaultdict

# Number of duplicated haplotype pairs per genotype
dups = 2
outfile = ""

# Returns the transposed matrix
def reorient_genos(genos_t):
	gtlen = len(genos_t)
	slen = len(genos_t[0])
	genos = [[genos_t[jj][ii] for jj in range(gtlen)] for ii in range(slen)]
	return genos

# Generates a random pair of haplotypes compatible with the genotype at a SNP
def rand_row(input_row, repeat_cnt):
	r_row = []
	for allele in input_row:
		if allele == '1':
			trow = [['0', '1'] if ch == 0 else ['1', '0'] for ch in np.random.binomial(1, 0.5, repeat_cnt)]
			r_row += list(itertools.chain.from_iterable(trow))
		elif allele == '0':
			r_row += ['0' for ii in range(2 * repeat_cnt)]
		else:
			r_row += ['1' for ii in range(2 * repeat_cnt)]
	return r_row

# Generates the entire set of compatible haplotypes
def rand_compat_hapls(input_mat, repeat_cnt):
	ilen = len(input_mat)
	return [rand_row(input_mat[ii], repeat_cnt) for ii in range(ilen)]

# Generates the next level of the localized haplotype cluster model
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

# Combines all fields of two nodes
def merge_two(n1, n2):
	merged = [n1[ii] + n2[ii] for ii in range(2)]
	merged.append(n1[2] + n2[2])
	merged.append(n1[3].union(n2[3]))
	return merged

# Merges all the nodes and returns a list containing that single node
def merge_all(hapl_dict):
	while len(hapl_dict) != 1:
		indices = hapl_dict.keys()
		hapl_dict[indices[0]] = merge_two(hapl_dict[indices[0]], hapl_dict[indices[1]])
		del hapl_dict[indices[1]]
		#hapl_dict.remove(indices[1])
	return [hapl_dict[hapl_dict.keys()[0]]]

# Computes pairwise threshold
def thresh(n1, n2):
	# Constant at the end may be tuned. Browning recommended .75 but .9 seems to
	# yield a more parsimonious representation without losing too much info
	return math.sqrt((1.0/sum(n1[2])) + (1.0/sum(n2[2]))) * math.sqrt(dups) * .9

# Finds the minimum pairwise threshold (used later as terminating condition of traversal)
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

# Computes the merge score between a pair of nodes, n1 and n2
def merge_score(n1, n2, hapls, idx2, min_diff, hapl_len):
	# Compute the threshold. If ever exceeded, return 100.0
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
	# the maximum observed difference over the course of this.
	merge_dict_1 = {"": n1[3]}
	merge_dict_2 = {"": n2[3]}
	# Windowing hyperparameter added for efficiency reasons
	depth = 0
	while idx2 < hapl_len and (len(merge_dict_1) != 0 or len(merge_dict_2) != 0) and depth < 40:
		depth += 1
		post_merge_1 = {}
		post_merge_2 = {}
		# Construct the next level
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
		# Check each edge 
		for h in poss_hapls:
			h1_cnt = 0
			h2_cnt = 0
			if h in merge_dict_1:
				h1_cnt = len(merge_dict_1[h])
			if h in merge_dict_2:
				h2_cnt = len(merge_dict_2[h])
			c = [h1_cnt/o1, h2_cnt/o2]
			diff = abs(c[0] - c[1])
			if diff >= thr_diff:
				return 100.0
			if diff >= max_diff:
				max_diff = diff
			# If maximal possible difference between two vertices is too
			# small, stop exploring
			if max(c) <= min_diff:
				if h in merge_dict_1:
					del merge_dict_1[h]
				if h in merge_dict_2:
					del merge_dict_2[h]
		idx2 += 1
	return max_diff


# Attempts to merge nodes that contain similar information
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
		#if ii % 100 == 0:
		#	print(ii)
		idx = ii + 1
		# split (i.e., create next level)
		l_row = split_l(l_graph[-1], hapls[ii])
		# merge (combine nodes whose futures are sufficiently similar)
		rows = merge_l(l_row, hapls, idx)
		l_graph.append(rows)
	return l_graph

# actually just constructs a haplotype HMM. diplotype unnecessary
def construct_dipl(local_graph):
    l_len = len(local_graph)
    tot = 0
    r1 = 0
    pis = []
    psym = []
    transitions = []
    t_prev_map = []
    lnodes = len(local_graph[1])
    for ii in range(lnodes):
        incoming = len(local_graph[1][ii][0])
        for jj in range(incoming):
            ecnt = local_graph[1][ii][2][jj]
            pis.append(float(ecnt))
            # store pair of node index and symbol
            psym.append((jj, local_graph[1][ii][1][jj]))
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
    return pis, transitions, t_prev_map, psym

# returns next_node
def find_next(cur_node, sym, t_prev_row):
	clen = len(t_prev_row[cur_node])
	for ii in range(clen):
		if t_prev_row[cur_node][ii][2] == sym:
			return t_prev_row[cur_node][ii][1]
	return -1

# returns (prob, next_node)
def find_next_with_prob(cur_node, sym, t_prev_row, hapl_hmm_row):
	clen = len(t_prev_row[cur_node])
	for ii in range(clen):
		if t_prev_row[cur_node][ii][2] == sym:
			prob = hapl_hmm_row[t_prev_row[cur_node][ii]]
			return prob, t_prev_row[cur_node][ii][1]
	return 0, -1

# generates new phased haplotypes given the haplotype hmm and the input matrix
def generate_new_samples(hapl_hmm, pis, t_prev, psym, input_mat):
	# construct pairs first and then split them up to match the format
	hlen = len(hapl_hmm)
	ilen = len(input_mat[0])
	# initialize based on the pis and compatibility
	# then, randomly choose the next state based on the probability of the
	# transition (which incorporates both the haplotypes)
	# map pairs of current vertices to next vertices (that satisfy genotype)
	# (cur_0, cur_1, next_geno) -> [[next_0, next_1, t_prob]
	next_hapls = []
	thapls = []
	node_list = []
	indices = [0, 1]
	# determine the start indices corresponding to the '1' and '0'
	plen = len(psym)
	for ii in range(plen):
		if psym[ii][1] == '0':
			indices[0] = ii
		else:
			indices[1] = ii
	# choose the first positions
	for jj in range(ilen):
		if input_mat[0][jj] == '2':
			for qq in range(dups):
				for lolz in range(2):
					thapls.append('1')
					node_list.append(indices[1])
		elif input_mat[0][jj] == '0':
			for qq in range(dups):
				for lolz in range(2):
					thapls.append('0')
					node_list.append(indices[0])
		else:
			hsyms = np.random.binomial(1, pis[indices[1]], dups)
			for qq in range(dups):
				if hsyms[qq] == 1:
					thapls.append('1')
					thapls.append('0')
					node_list.append(indices[1])
					node_list.append(indices[0])
				else:
					thapls.append('0')
					thapls.append('1')
					node_list.append(indices[0])
					node_list.append(indices[1])
	next_hapls.append(thapls)
	# for each following position, randomly sample a pair of compatible haplotypes
	for ii in range(hlen):
		thapls = []
		next_nodes = []
		t_prev_row = t_prev[ii]
		hapl_row = hapl_hmm[ii]
		dipls = {}
		nl_idx = 0
		for jj in range(ilen):
			if input_mat[ii + 1][jj] == '2':
				for qq in range(dups):
					for lolz in range(2):
						thapls.append('1')
						nn = find_next(node_list[nl_idx], '1', t_prev_row)
						if nn == -1:
							next_nodes.append(find_next(node_list[nl_idx], '0', t_prev_row))
						else:
							next_nodes.append(nn)
						nl_idx += 1
			elif input_mat[ii + 1][jj] == '0':
				for qq in range(dups):
					for lolz in range(2):
						thapls.append('0')
						nn = find_next(node_list[nl_idx], '0', t_prev_row)
						if nn == -1:
							# prior approach to generating next pair
							'''for ele in hapl_row.keys():
								if ele[2] == '0':
									next_nodes.append(ele[1])
									break'''
							next_nodes.append(find_next(node_list[nl_idx], '1', t_prev_row))
						else:
							next_nodes.append(nn)
						nl_idx += 1
			else:
				for qq in range(dups):
					zero_zero, n00 = find_next_with_prob(node_list[nl_idx], '0', t_prev_row, hapl_row)
					zero_one, n01 = find_next_with_prob(node_list[nl_idx], '1', t_prev_row, hapl_row)
					one_zero, n10 = find_next_with_prob(node_list[nl_idx + 1], '0', t_prev_row, hapl_row)
					one_one, n11 = find_next_with_prob(node_list[nl_idx + 1], '1', t_prev_row, hapl_row)
					first = zero_one * one_zero
					second = zero_zero * one_one
					# if this happens, randomly move over to some valid pair of nodes
					# not a great solution, but we'll see if it does at least ok...
					if first == 0.0 and second == 0.0:
						if np.random.binomial(1, 0.5, 1) == 1:
							if zero_zero != 0.0:
								thapls.append('0')
								thapls.append('1')
								next_nodes.append(n00)
								next_nodes.append(n10)
							else:
								thapls.append('1')
								thapls.append('0')
								next_nodes.append(n01)
								next_nodes.append(n11)
						else:
							if one_one != 0.0:
								thapls.append('0')
								thapls.append('1')
								next_nodes.append(n01)
								next_nodes.append(n11)
							else:
								thapls.append('1')
								thapls.append('0')
								next_nodes.append(n00)
								next_nodes.append(n10)
						nl_idx += 2
						#sys.exit(-1)
					else:
						one_z_prob = first / (first + second)
						hsym = np.random.binomial(1, one_z_prob, 1)
						if hsym == 1:
							thapls.append('1')
							thapls.append('0')
							next_nodes.append(n01)
							next_nodes.append(n10)
						else:
							thapls.append('0')
							thapls.append('1')
							next_nodes.append(n00)
							next_nodes.append(n11)
						nl_idx += 2
		next_hapls.append(thapls)
		node_list = next_nodes
	return next_hapls

# checks if moving from a given pair onwards improves likelihood. if so, update next_hapl_pairs 
def update_next_hapl_pairs(next_hapl_pairs, hapl_pair_dict, idx_pair, pair, t_prev_row, hapl_row):
	next = [-1, -1]
	lprob = [-1000.0, -1000.0]
	for lolz in range(2):
		prob, nn = find_next_with_prob(idx_pair[lolz], pair[lolz], t_prev_row, hapl_row)
		if nn == -1:
			maxn = -1
			maxlprob = -1000.0
			for ele in hapl_row.keys():
				if ele[2] == pair[lolz] and math.log10(hapl_row[ele]) - 8 > maxlprob:
					maxlprob = math.log10(hapl_row[ele]) - 8
					maxn = ele[1]
					break
			next[lolz] = maxn
			lprob[lolz] = maxlprob
		else:
			next[lolz] = nn
			lprob[lolz] = math.log10(prob)
	if next[0] == -1 or next[1] == -1:
		return
	# {(next_0, next_1) : [sum log probs, sym_0, sym_1, (cur_0, cur_1)]}
	if (next[0], next[1]) not in next_hapl_pairs:
		next_hapl_pairs[(next[0], next[1])] = \
			[hapl_pair_dict[idx_pair] + sum(lprob), pair[0], pair[1], idx_pair]
	elif hapl_pair_dict[idx_pair] + sum(lprob) > next_hapl_pairs[(next[0], next[1])][0]:
		next_hapl_pairs[(next[0], next[1])] = \
			[hapl_pair_dict[idx_pair] + sum(lprob), pair[0], pair[1], idx_pair]
	return

# computes the most likely haplotypes for each individual, conditional on 
# the HMM and genotype
def dipl_viterbi(pis, hapl_hmm, t_prev, psym, input_mat):
	input_t = reorient_genos(input_mat)
	glen = len(input_t)
	hlen = len(hapl_hmm)
	mp_hapls = []
	indices = [0, 1]
	# determine the start indices corresponding to the '1' and '0'
	plen = len(psym)
	for ii in range(plen):
		if psym[ii][1] == '0':
			indices[0] = ii
		else:
			indices[1] = ii
	# [[sum log probs, cur_0, cur_1], ...] <--- maybe not necessary
	# {(cur_0, cur_1) : sum_log_probs, (cur_0, cur_1) : idx_1,}
	# {(next_0, next_1) : [sum log probs, sym_0, sym_1, (cur_0, cur_1)]}
	# After constructing all of these, generate the next 
	# ok, I should maintain a series of backpointers
	# the first level should just contain the symbols associated
	# [{(0, 2) : ('0', '1')}, {(0, 2) : ('0', '1', (1, 2)), 
	# (2, 0) : ('1', '0', (2, 1))}, {(0, 0) : ('0', '1', (0, 2))}]
	# for each genotype, traverse the HMM to find the most likely pair of
	# haplotypes. merge them when they arrive at the same pair of vertices
	for ii in range(glen):
		# first, initialize the haplotype pair
		# also initialize the probability, a sum of logs
		# for each proposed path, include information in backpointers
		# then, for the next layer, simply calculate all the possible pairs of
		# next nodes. for each pair, find the maximum probability associated
		# with it. if we don't find any acceptable option, I guess we should
		# do the random jump to a compatible pair (refine later, perhaps)
		hapl_pair_dict = {}
		back_pair_ptrs = []
		if input_t[ii][0] == '2':
			idx = indices[1]
			back_pair_ptrs.append({(idx, idx) : ('1', '1')})
			hapl_pair_dict[(idx, idx)] = 2.0 * math.log10(pis[idx])
		elif input_t[ii][0] == '0':
			idx = indices[0]
			back_pair_ptrs.append({(idx, idx) : ('0', '0')})
			hapl_pair_dict[(idx, idx)] = 2.0 * math.log10(pis[idx])
		else:
			# oh, I can just arbitrarily calculate ('0', '1'): there's no point
			# in doing ('1', '0') here too
			ids = (indices[1], indices[0])
			back_pair_ptrs.append({ids : ('1', '0')})
			hapl_pair_dict[ids] = math.log10(pis[indices[1]]) + math.log10(pis[indices[0]])
		# now construct each layer
		for jj in range(hlen):
			# determine the best next_hapl_pairs
			next_hapl_pairs = {}
			t_prev_row = t_prev[jj]
			hapl_row = hapl_hmm[jj]
			for idx_pair in hapl_pair_dict:
				# Quick thought: we can penalize the transitions pretty easily!
				# Just add -8 when we find we're stuck and then pick a new node
				if input_t[ii][jj + 1] == '2':
					update_next_hapl_pairs(next_hapl_pairs, hapl_pair_dict, idx_pair, ['1', '1'], t_prev_row, hapl_row)
				elif input_t[ii][jj + 1] == '0':
					update_next_hapl_pairs(next_hapl_pairs, hapl_pair_dict, idx_pair, ['0', '0'], t_prev_row, hapl_row)
				else:
					update_next_hapl_pairs(next_hapl_pairs, hapl_pair_dict, idx_pair, ['1', '0'], t_prev_row, hapl_row)
					update_next_hapl_pairs(next_hapl_pairs, hapl_pair_dict, idx_pair, ['0', '1'], t_prev_row, hapl_row)
			# then update hapl_pair_dict and back_pair_ptrs
			hapl_pair_dict = {}
			tback = {}
			for pair in next_hapl_pairs:
				res = next_hapl_pairs[pair]
				hapl_pair_dict[pair] = res[0]
				tback[pair] = (res[1], res[2], res[3])
			back_pair_ptrs.append(tback)
		# then backtrack on back_pair_ptrs and then reverse the generated
		# haplotypes
		rev_hapls = [[], []]
		vertices = back_pair_ptrs[-1].keys()[0]
		for back_ptrs in reversed(back_pair_ptrs):
			res = back_ptrs[vertices]
			''' Leaving this here just in case I later discover a problem
			in backtracking
			try: 
				res = back_ptrs[vertices]
			except:
				print('Fix this stupid problem')
				res = back_ptrs.values()[0]'''
			rev_hapls[0].append(res[0])
			rev_hapls[1].append(res[1])
			if len(res) == 3:
				vertices = res[2]
		mp_hapls.append(rev_hapls[0][::-1])
		mp_hapls.append(rev_hapls[1][::-1])
	return mp_hapls

# Prints the haplotypes SNP by SNP in specified format
def print_haplotypes(h_t):
	with(open(outfile, "a+")) as f:
		rlen = len(h_t[0])
		row_cnt = len(h_t)
		for row in h_t:
			for ii in range(rlen - 1):
				f.write(row[ii] + ' ')
			f.write(row[rlen - 1] + '\n')
	return

def main(fname):
	print(fname)
	global outfile
	# Change this line when testing to avoid overwriting solution files
	outfile = fname[:-4] + "_sol.txt"
	print(outfile)
	input_mat = []
	with(open(fname, 'rb')) as f:
		reader = csv.reader(f, delimiter=' ')
		for row in reader:
			input_mat.append(row)
	# For each genotype, randomly generate n = dups compatible haplotype pairs
	hapls = rand_compat_hapls(input_mat, dups)
	local_graph = []
	hapl_hmm = []
	input_rev = input_mat[::-1]
	print('Constructed random haplotypes')
	iters = 6
	# It converges decently well with 4 iterations on large datasets, based on
	# analysis of example_data_3.txt performance
	if len(input_mat) > 50000:
		iters = 4
	for ii in range(iters):
		local_graph = localized_hapl_cluster(hapls)
		lens = []
		'''for ele in local_graph:
			lens.append(len(ele))
		print(lens)'''
		'''print('-------')
		print(len(local_graph))'''
		# In hindsight, this is just a haplotype HMM which is traversed by a pair
		pis, hapl_hmm, t_prev, psym = construct_dipl(local_graph)
		#if ii % 2 == 0:
		#	print_haplotypes(reorient_genos(dipl_viterbi(pis, hapl_hmm, t_prev, psym, input_mat)))
		# Reverse the direction with every iteration
		if ii % 2 == 0:
			hapls = generate_new_samples(hapl_hmm, pis, t_prev, psym, input_mat)
		else:
			hapls = generate_new_samples(hapl_hmm, pis, t_prev, psym, input_rev)
		hapls.reverse()
		print('Completed iteration {}'.format(ii))
	local_graph = localized_hapl_cluster(hapls)
	pis, hapl_hmm, t_prev, psym = construct_dipl(local_graph)
	print('Running Viterbi algorithm on final HMM')
	hapls = dipl_viterbi(pis, hapl_hmm, t_prev, psym, input_mat)
	hapls_t = reorient_genos(hapls)
	print_haplotypes(hapls_t)

if __name__=="__main__":
	if len(sys.argv) != 2:
		print("Incorrect number of arguments. Provide one filename.")
	else:
		main(sys.argv[1])
