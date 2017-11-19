#!/usr/bin/env python3

import random
import dag
import networks


sizes = [8, 10, 12, 16, 20, 24, 32, 40, 50, 64, 80, 100, 128, 160, 202, 256, 322, 406, 512]
half_sizes = [int(i/2) for i in sizes]

def generate_random_sparses():
	print("Generating random sparse posets...")
	for m in range(5):
		for k in [3,5]:
			for n in sizes:
				dag.save("avgdeg_%i_%03i_%i" % (k, n, m), dag.random_average_degree(n, k))

def generate_random_bipartites():
	print("Generating random bipartite posets...")
	for m in range(5):
		for p in [0.2, 0.5]:
			for n in half_sizes:
				dag.save("bipartite_%s_%03i_%i" % (str(p), 2 * n, m), dag.random_bipartite(n, p))
				dag.save("bipartite_%s_%03i_%i" % (str(p), 2 * n, m), dag.random_bipartite(n, p))

def generate_bayesian_network():
	print("Generating Bayesian network posets...")
	network_names = ["andes", "diabetes", "link", "pigs", "munin"]
	for network in network_names:
		path = "networks/%s.net" % network
		net = networks.read_net(path)
		
		for i in sizes:
			for m in range(5):
				dg = dag.gen_source_network(net, i)
				if dg is None:
					continue
				dag.save("bayesiannetwork_%s_%03i_%i" % (network, i, m), dg)


def generate_all():
	random.seed(0)
	generate_random_sparses()
	
	random.seed(0)
	generate_random_bipartites()
	
	random.seed(0)
	generate_bayesian_network()

generate_all()

