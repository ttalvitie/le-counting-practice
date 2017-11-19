#ifndef RELAXATION_H
#define RELAXATION_H

#define RELAXATION_ADAPTIVE 0
#define RELAXATION_MANUAL 1

#define RELAXATION_TRANSPOSE_NO 0
#define RELAXATION_TRANSPOSE_FULL 1
#define RELAXATION_TRANSPOSE_PARTS 2


struct RelaxationOptions
{
	int relaxation;
	int partition_max_size;
	
	int transpose;
	bool use_iterative_sampling;
	
	int n_partition;
	double partition_search_time;
	
	const char *visprefix;
	
	bool verbose;
	
	RelaxationOptions()
	{
		relaxation = RELAXATION_ADAPTIVE;
		partition_max_size = 40;
		
		transpose = RELAXATION_TRANSPOSE_FULL;
		use_iterative_sampling = true;
		
		n_partition = 0;
		partition_search_time = 0.0;
		
		visprefix = NULL;
		verbose = true;
	}
};



template <typename N>
struct Part
{
	unsigned n;
	int *elements;
	
	Part() {}
	
	RecursiveCounterAuto<N> *sampler;
	
	int *sample;
	int i;
	
	bool iterative;
	
	N n_extensions;
	
	// initializes sampling a subset of the given relaxation
	// returns true on success and false if memout was detected (all allocations made will be freed)
	bool init_sampling(digraph relaxation, RecursionOptions &opt, bool use_iterative)
	{
		// first count the linear extensions
		digraph poset(n, relaxation, elements);
		sampler = get_recursive_counter<N>(n, opt);
		
		try {
			n_extensions = sampler->count_linex(poset);
		} catch (std::bad_alloc) {
			sampler->deallocate();
			delete sampler;
			return false;
		}
		
		// decide whether to use iterative sampling
		iterative = use_iterative && !sampler->is_transposed();
		
		// initialize sampling
		if (!iterative) sample = new int[n];
		sampler->initialize_sampling();
		return true;
	}
	
	// uninitializes sampling, deallocating all data structures
	void uninit_sampling()
	{
		// uninitialize sampling
		sampler->uninitialize_sampling();
		
		// delete the sampling data structures
		if (!iterative) delete [] sample;
		sampler->deallocate();
		delete sampler;
	}
	
	// initializes drawing of a single sample, must be called once for every sample before invoking get_next_element()
	void init_draw_sample()
	{
		// if using iterative sampling, initialize the iterative sampler,
		// otherwise just draw the full sample and store it for sequential queries
		if (iterative) {
			sampler->init_draw_iterative_sample();
		} else {
			sampler->draw_sample(sample);
		}
		
		// the index of the next element in the sample
		i = 0;
	}
	
	// computes and returns the next element in the sample
	int get_next_element()
	{
		// update the index
		int index = i++;
		
		// if using iterative sampling, query the next element from the sampler,
		// otherwise return the i'th element in the precomputed sample
		if (iterative) {
			return sampler->draw_next_element();
		} else {
			return sample[index];
		}
	}
	
	std::string sampling_method()
	{
		if (iterative) return std::string("iterative");
		return std::string("complete");
	}
};


template <typename N>
struct RelaxationSampler
{
	RelaxationOptions options;
	
	unsigned n;
	digraph poset;
	digraph relaxation;
	
	int *sample;
	
	int *partition;
	
	N n_interleavings;
	N n_orders;
	
	Part<N> **parts;
	int n_parts;
	
	int **predecessors;
	int *n_predecessors;
	
	
	// *** partition modification functions
	
	// if p is not connected, it's modified to be a connected component and a new part is created for the remainder
	void split_into_connected_component(int p)
	{
		Part<N> *part = parts[p];
		
		// find a connected component
		int component[n];
		unsigned k = poset.find_connected_component(component, part->elements, part->n);
		
		// if the component contains all nodes
		if (k == part->n) return;
		
		bool subset[n];
		for (int i = 0; i < n; ++i) subset[i] = false;
		for (int i = 0; i < k; ++i) subset[component[i]] = true;
		
		Part<N> *new_part = new Part<N>();
		parts[n_parts++] = new_part;
		new_part->elements = part->elements + k;
		new_part->n = part->n - k;
		part->n = k;
		
		int i2 = k;
		for (int i1 = 0; i1 < k; ++i1) {
			if (subset[part->elements[i1]]) continue;
			while (!subset[part->elements[i2]]) ++i2;
			std::swap(part->elements[i1], part->elements[i2]);
		}
	}
	
	// splits an existing partition into connected components
	void split_partition_into_components()
	{
		// NOTE: n_parts may change during iteration and it works as intended
		for (int i = 0; i < n_parts; ++i) {
			split_into_connected_component(i);
		}
	}
	
	// splits part p in two, the first k elements remain in part p
	void split_part(int p, unsigned k)
	{
		int q = n_parts++;
		Part<N> *new_part = new Part<N>();
		parts[q] = new_part;
		
		new_part->elements = parts[p]->elements + k;
		new_part->n = parts[p]->n - k;
		parts[p]->n = k;
	}
	
	// splits part p into parts of given size (the last part can be smaller)
	// the number of resulting parts is returned and the part indices are stored in split_parts
	int split_part_max_size(int p, int size, int *split_parts)
	{
		int n_split_parts = 1;
		split_parts[0] = p;
		
		// while the remainder is larger
		while (parts[p]->n > size) {
			// split the part into a part of given size and the remainder
			split_part(p, size);
			p = n_parts - 1;
			split_parts[n_split_parts++] = p;
		}
		
		return n_split_parts;
	}
	
	
	// *** functions for optimizing a partition with element swaps
	
	// randomly redistributes elements between given parts
	void shuffle_partition(int *swap_parts, int n_swap_parts)
	{
		// get positions of all elements in the given parts
		int *elements[n];
		int n_swap_elements = 0;
		for (int i = 0; i < n_swap_parts; ++i) {
			Part<N> *part = parts[swap_parts[i]];
			for (int j = 0; j < part->n; ++j) {
				elements[n_swap_elements++] = &part->elements[j];
			}
		}
		
		// repermute the elements
		for (unsigned i = 0; i < n_swap_elements-1; ++i) {
			int j = i + rnd() * (n_swap_elements - i);
			std::swap(*elements[i], *elements[j]);
		}
	}
	
	// returns the number of comparable elements between parts i and j
	int comparable_elements(int i, int j)
	{
		int *elemA = parts[i]->elements;
		int *elemB = parts[j]->elements;
		unsigned nA = parts[i]->n;
		unsigned nB = parts[j]->n;
		int n_comparable = 0;
		
		for (unsigned a = 0; a < nA; ++a) {
			for (unsigned b = 0; b < nB; ++b) {
				if (poset.comparable(elemA[a], elemB[b])) ++n_comparable;
			}
		}
		
		return n_comparable;
	}
	
	// evaluates the partition of the given parts in terms of comparable pairs between parts
	int evaluate_partition(int *parts, int n_parts)
	{
		int n_comparable = 0;
		
		for (int i = 0; i < n_parts-1; ++i) {
			for (int j = i+1; j < n_parts; ++j) {
				n_comparable += comparable_elements(parts[i], parts[j]);
			}
		}
		
		return n_comparable;
	}
	
	// gets the pair of elements in parts pi and pj that results in the best improvement
	// in the score by swapping them
	int get_best_swap(int pi, int pj, int &best_a, int &best_b)
	{
		int *elemA = parts[pi]->elements;
		int *elemB = parts[pj]->elements;
		unsigned nA = parts[pi]->n;
		unsigned nB = parts[pj]->n;
		
		int swap_a_change[nA];
		int swap_b_change[nB];
		
		for (unsigned a = 0; a < nA; ++a) {
			swap_a_change[a] = 0;
			// if a in A is swapped to B, it will lose comparable pairs in B and gain comparable pairs in A
			for (unsigned i = 0; i < nB; ++i) if (poset.comparable(elemA[a], elemB[i])) --swap_a_change[a];
			for (unsigned i = 0; i < nA; ++i) if (poset.comparable(elemA[a], elemA[i])) ++swap_a_change[a];
		}
		
		for (unsigned b = 0; b < nB; ++b) {
			swap_b_change[b] = 0;
			// if b in B is swapped to A, it will lose comparable pairs in A and gain comparable pairs in B
			for (unsigned i = 0; i < nA; ++i) if (poset.comparable(elemB[b], elemA[i])) --swap_b_change[b];
			for (unsigned i = 0; i < nB; ++i) if (poset.comparable(elemB[b], elemB[i])) ++swap_b_change[b];
		}
		
		int best_change = 0;
		
		for (unsigned a = 0; a < nA; ++a) {
			for (unsigned b = 0; b < nB; ++b) {
				// the change in score after swapping a and b
				int swap_ab_change = swap_a_change[a] + swap_b_change[b];
				
				// if elemA[a] and elemB[b] are comparable, only the loss is counted (twice) above
				if (poset.comparable(elemA[a], elemB[b])) swap_ab_change += 2;
				
				if (swap_ab_change < best_change) {
					best_change = swap_ab_change;
					best_a = a;
					best_b = b;
				}
			}
		}
		
		return best_change;
	}
	
	// gets the best swap between all given parts (or 0 if no swap improves the score)
	int get_best_swap(int *swap_parts, int n_swap_parts, int &best_pa, int &best_pb, int &best_a, int &best_b)
	{
		int best_change = 0;
		for (int i = 0; i < n_swap_parts-1; ++i) {
			int pa = swap_parts[i];
			for (int j = i+1; j < n_swap_parts; ++j) {
				int pb = swap_parts[j];
				
				int a, b;
				int change = get_best_swap(pa, pb, a, b);
				if (change < best_change) {
					best_pa = pa;
					best_pb = pb;
					best_a = a;
					best_b = b;
					best_change = change;
				}
			}
		}
		
		return best_change;
	}
	
	// optimizes a k-partition greedily, by iteratively swapping elements between
	// parts in a way that results in the best local improvement in the score,
	// *not* guaranteed to find a global optimum
	int local_optimize_partition(int *swap_parts, int n_swap_parts)
	{
		int score = evaluate_partition(swap_parts, n_swap_parts);
		
		int pa, pb, a, b;
		int change;
		while ((change = get_best_swap(swap_parts, n_swap_parts, pa, pb, a, b)) < 0) {
			std::swap(parts[pa]->elements[a], parts[pb]->elements[b]);
			score += change;
		}
		
		return score;
	}
	
	void sample_k_partitions(int *swap_parts, int n_swap_parts)
	{
		N first_score = -1;
		N best_score = -1;
		int best[n];
		
		double start = seconds();
		
		if (options.verbose) {
			eprintf("Sampling partitions: ");
			eflush();
		}
		
		while (best_score == -1 || seconds() - start < options.partition_search_time) {
			shuffle_partition(swap_parts, n_swap_parts);
			local_optimize_partition(swap_parts, n_swap_parts);
			if (options.partition_search_time <= 0.0) return;
			
			if (!init_sampler()) continue;
			uninit_sampler();
			
			N score = n_orders;
			if (best_score == -1 || score < best_score) {
				if (options.verbose) eprintf("*");
				if (best_score == -1) first_score = score;
				best_score = score;
				std::copy(partition, partition + n, best);
			}
			eflush();
		}
		if (options.verbose) eprintf("\nImprovement over the first partition: %.2Lf\n", (long double)first_score / best_score);
		
		std::copy(best, best + n, partition);
	}
	
	
	
	// splits part p into parts of given size (the last part can be smaller)
	void choose_partition_relaxation_max(int p, int size)
	{
		int split_parts[n];
		int n_split_parts = split_part_max_size(p, size, split_parts);
		
		// optimize the partition between the split parts
		sample_k_partitions(split_parts, n_split_parts);
		
		if (options.verbose) eprintf("Comparable pairs: %i\n", evaluate_partition(split_parts, n_split_parts));
	}
	
	// chooses a partition relaxation on part p
	// (root function for choosing a partition-based relaxation)
	void choose_partition_relaxation(int p)
	{
		// if p has less at most 20 elements, don't partition it further
		if (parts[p]->n <= 20) return;
		
		choose_partition_relaxation_max(p, options.partition_max_size);
	}
	
	void remove_arcs_between(int i, int j)
	{
		int *elemA = parts[i]->elements;
		int *elemB = parts[j]->elements;
		unsigned nA = parts[i]->n;
		unsigned nB = parts[j]->n;
		
		for (unsigned a = 0; a < nA; ++a) {
			for (unsigned b = 0; b < nB; ++b) {
				relaxation.pair(elemA[a], elemB[b]) = false;
				relaxation.pair(elemB[b], elemA[a]) = false;
			}
		}
	}
	
	void remove_arcs_between_parts()
	{
		for (int i = 0; i < n_parts-1; ++i) {
			for (int j = i+1; j < n_parts; ++j) {
				remove_arcs_between(i, j);
			}
		}
	}
	
	void choose_partition_relaxation()
	{
		if (options.verbose) eprintf("Choosing a partition relaxation.\n");
		
		// choose a partition for the relaxation for each (current) part
		int n_parts_original = n_parts;
		for (int i = 0; i < n_parts_original; ++i) {
			choose_partition_relaxation(i);
		}
		
		// actually remove the arcs between parts and resplit just in case some parts are not connected
		remove_arcs_between_parts();
		split_partition_into_components();
	}
	
	void sort_parts_by_size() const
	{
		struct {
			bool operator() (Part<N> *p, Part<N> *q) {
				return p->n > q->n;
			}
		} descending_by_size;
		
		std::sort(parts, parts + n_parts, descending_by_size);
	}
	
	// the root function for choosing a relaxation via various strategies
	void choose_relaxation()
	{
		// start by making the relaxation a copy of the original poset
		relaxation.make_copy(poset);
		
		// detect connected components
		split_partition_into_components();
		
		// next choose a relaxation
		choose_partition_relaxation();
		
		sort_parts_by_size();
	}
	
	
	// initializes sampling for the currently chosen relaxation,
	// runs the exact counting of linear extensions for each part and builds up sampling data structures.
	// this returns true if sampling was initialized successfully and false if memout was detected,
	// in which case all work done will be uninitialized as if this function had not been called.
	bool init_sampler()
	{
		RecursionOptions opt;
		opt.sampling = true;
		opt.sort_topologically = false;
		opt.verbose = false;
		opt.transpose = options.transpose == RELAXATION_TRANSPOSE_PARTS ? TRANSPOSE_AUTO2 : TRANSPOSE_NO;
		
		// compute the number of interleavings between parts (n! / sum_i n_i!)
		Combinatorial<N> comb(n);
		comb.compute_factorials();
		n_interleavings = comb.factorial(n);
		for (int p = 0; p < n_parts; ++p) n_interleavings /= comb.factorial(parts[p]->n);
		
		// build a sampler for all parts and count the linear extensions
		n_orders = n_interleavings;
		for (int p = 0; p < n_parts; ++p) {
			if (options.verbose) eprintf("Building a sampler for part %i, size %i   \n", p, parts[p]->n);
			
			// build the sampler (if fails, uninitialize all previous parts)
			if (!parts[p]->init_sampling(relaxation, opt, options.use_iterative_sampling)) {
				if (options.verbose) eprintf("\nOut of memory, undoing sampling initialization.\n");
				for (int q = 0; q < p; ++q) parts[q]->uninit_sampling();
				return false;
			}
			
			// update the number of linear extensions
			n_orders *= parts[p]->n_extensions;
			
			if (options.verbose) eprintf("\x1b[A\r");
		}
		if (options.verbose) eprintf("\n");
		
		return true;
	}
	
	// uninitializes sampling
	void uninit_sampler()
	{
		// uninitialize the sampler for each part
		for (int p = 0; p < n_parts; ++p) {
			parts[p]->uninit_sampling();
		}
	}
	
	void preprocess_original_poset()
	{
		// the transitive closure will be necessary to maintain all information
		// of the order relation wihin arbitrary subsets of the poset
		poset.take_transitive_closure();
		
		// run a heuristic to determine whether it helps to tranpose the poset
		RecursionOptions opt;
		opt.verbose = false;
		RecursiveCounterAuto<N> *rca = get_recursive_counter<N>(n, opt);
		rca->set_poset(poset);
		if (options.transpose != RELAXATION_TRANSPOSE_NO && rca->transpose_heuristic2()) {
			if (options.verbose) eprintf("Transposing the original poset.\n")
			poset.invert();
		}
		delete rca;
		
		// compute predecessor lists
		digraph tr(poset);
		tr.remove_transitive_arcs();
		n_predecessors = new int[n];
		predecessors = new int*[n];
		for (unsigned i = 0; i < n; ++i) {
			predecessors[i] = new int[n];
			n_predecessors[i] = 0;
			for (unsigned j = 0; j < n; ++j) {
				if (tr.has(j, i)) predecessors[i][n_predecessors[i]++] = j;
			}
		}
	}
	
	RelaxationSampler(RelaxationOptions options_, const digraph &poset_) :
		options(options_),
		n(poset_.n),
		poset(poset_),
		relaxation(poset_.n),
		sample(new int[n]),
		partition(new int[n]),
		parts(new Part<N>*[n]),
		n_parts(1)
	{
		preprocess_original_poset();
		
		for (unsigned i = 0; i < n; ++i) partition[i] = i;
		parts[0] = new Part<N>();
		parts[0]->elements = partition;
		parts[0]->n = n;
		
		choose_relaxation();
		
		if (options.verbose) eprintf("==== Relaxation chosen ====\n");
		
		remove_arcs_between_parts();
		if (options.visprefix) visualize_partition();
	}
	
	~RelaxationSampler()
	{
		uninit_sampler();
		
		delete [] sample;
		delete [] partition;
		
		for (int p = 0; p < n_parts; ++p) {
			delete parts[p];
		}
		
		delete [] parts;
		
		for (unsigned i = 0; i < n; ++i) delete [] predecessors[i];
		delete [] n_predecessors;
		delete [] predecessors;
	}
	
	// attempts to draw a sample, returns true if the sample was accepted, false if rejected
	// on success, the sample is accessible via the 'sample' member
	bool draw_sample()
	{
		// initialize drawing a single sample for each part
		for (int p = 0; p < n_parts; ++p) {
			parts[p]->init_draw_sample();
		}
		
		bool selected[n];
		for (unsigned i = 0; i < n; ++i) selected[i] = false;
		
		// interleave the samples at random (uniformly)
		for (unsigned i = 0; i < n; ++i) {
			// we pick the next element in a part of the partition chosen at random
			// the probability of choosing each part is proportional to the number of its remaining elements
			int r = rnd() * (n - i);
			int cumul = 0;
			int u = -1;
			
			// loop through the parts until the cumulative probability exceeds r
			for (int p = 0; p < n_parts; ++p) {
				cumul += (parts[p]->n - parts[p]->i);
				if (cumul <= r) continue;
				
				u = parts[p]->elements[parts[p]->get_next_element()];
				break;
			}
			
			// add the element to the sample and mark it as selected
			sample[i] = u;
			selected[u] = true;
			
			// if a predecessor of u has not been selected yet, we can reject the sample
			for (unsigned j = 0; j < n_predecessors[u]; ++j) {
				int v = predecessors[u][j];
				if (poset.has(v, u) && !selected[v]) return false;
			}
		}
		
		// the sample respects the poset
		return true;
	}
	
	
	// display / visualization
	
	void visualize_partition()
	{
		if (options.verbose) eprintf("Visualizing... "); eflush();
		
		DOTOptions dot_opt;
		const char *node_colors[n];
		
		static const char *part_colors[] = {"#ff2222", "#6666ff", "#00bb00", "#ffaa00", "#22ffff", "#ff22ff", "#bb0000", "#00bbbb"};
		
		for (int p = 0; p < n_parts; ++p) {
			int *elements = parts[p]->elements;
			for (unsigned i = 0; i < parts[p]->n; ++i) {
				node_colors[elements[i]] = p < 8 ? part_colors[p] : "#000000";
			}
		}
		
		digraph reduction(relaxation);
		reduction.remove_transitive_arcs();
		
		dot_opt.node_colors = node_colors;
		dot_opt.edge_dir = "none";
		FILE *f = fopen((std::string(options.visprefix) + "partition.dot").c_str(), "w");
		reduction.print_dot(f, dot_opt);
		fclose(f);
		
		if (options.verbose) eprintf("Done.\n");
	}
	
	void display_relaxation()
	{
		eprintf("Part   Counting  Sampling     Size   Extensions\n");
		for (int p = 0; p < n_parts; ++p) {
			eprintf(" %2i:   %-6s    %-9s   %3i     ",
				p,
				parts[p]->sampler->short_description().c_str(),
				parts[p]->sampling_method().c_str(),
				parts[p]->n);
			std::cerr << parts[p]->sampler->get_extension_count() << std::endl;
		}
		
		eprintf("Interleavings: %Le\n", (long double)n_interleavings);
		eprintf("Relaxation extensions: %Le\n", (long double)n_orders);
	}
};


#endif
