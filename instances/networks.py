
import random
import itertools



# NOTE: not a perfect parser but should work for any network from www.bnlearn.com/bnrepository
def read_net(path):
	with open(path, "r") as f:
		lines = f.read().split("\n")
	
	parents = {}
	
	for line in lines:
		line = line.strip()
		
		if not line.startswith("potential "):
			continue
		
		parts = line.split()
		assert(parts[1] == "(")
		assert(parts[-1] == ")")
		parts = parts[2:-1]
		node = parts[0]
		
		if len(parts) == 1:
			parents[node] = []
		else:
			assert(parts[1] == "|")
			parents[node] = parts[2:]
	
	def toposort(node):
		if node in order:
			return
		for parent in parents[node]:
			toposort(parent)
		order.append(node)
	
	order = []
	for node in sorted(parents):
		toposort(node)
	assert(len(order) == len(parents))
	
	return [[n, parents[n]] for n in order]







