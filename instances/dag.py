import random
import itertools
import os
import shutil

def matrix(dag):
	return "\n".join(" ".join(str(b) for b in a) for a in dag) + "\n"

def dagsave(dag, path):
	with open(path, "w") as f:
		f.write(matrix(dag))

def subgraph(dg, subset):
	n = len(subset)
	subdg = empty(n)
	
	for i in range(n):
		for j in range(n):
			subdg[i][j] = dg[subset[i]][subset[j]]
	
	return subdg

def empty(n):
	return [[0] * n for i in range(n)]

# Returns a random n x n bipartite DAG.
# The poset is constructed by splitting the elements into sets A and B and adding an edge (a,b)
# for all a in A, b in B with probability p.
def random_bipartite(n, p):
	dag = empty(2 * n)
	for i in range(n):
		for j in range(n):
			if random.random() < p:
				dag[i][n + j] = 1
	return dag


# Returns a random DAG.
# The DAG is constructed by choosing a uniform random order of the elements and
# adding each edge conforming to the order with probability p.
def random_density(n, p):
	dag = empty(n)
	order = list(range(n))
	random.shuffle(order)
	
	for i in range(n-1):
		for j in range(i+1, n):
			if random.random() < p:
				dag[order[i]][order[j]] = 1
	
	return dag


# Returns a random DAG.
# The DAG is constructed by choosing a uniform random order of the elements and
# adding each edge conforming to the order with probability p = k / (n-1), which
# ensures that the expected average degree of the DAG is k.
def random_average_degree(n, k):
	return random_density(n, k / (n-1))

def gen_source_network(net, k):
	def network_to_digraph(parents):
		dag = empty(n)
		for node in range(n):
			for parent in parents[node]:
				dag[node][parent] = 1
		return dag
	
	# read network
	n = len(net)
	
	if n < k:
		return None
		#raise Exception("Source network has less than %i nodes." % k)
	
	# switch to numerical nodes
	nodes = [nd[0] for nd in net]
	nodes = dict((nd,i) for (i,nd) in enumerate(nodes))
	parents = [[nodes[p] for p in ps] for nd,ps in net]
	
	dag = network_to_digraph(parents)
	
	subset = set()
	candidates = set([random.choice(range(n))])
	while len(subset) < k:
		# if there are candidates, add one of them at ranodm
		# otherwise add one of any nodes that have not yet been added
		if candidates:
			node = random.choice(list(candidates))
			candidates.remove(node)
		else:
			not_added = [nd for nd in range(n) if nd not in subset]
			node = random.choice(not_added)
		
		subset.add(node)
		candidates |= set(nd for nd in range(n) if nd not in subset and (dag[node][nd] or dag[nd][node]))
	
	return subgraph(dag, sorted(subset))


def save(name, poset):
	dagsave(poset, "%s.txt" % name)
