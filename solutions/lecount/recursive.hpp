#ifndef RECURSIVE_H
#define RECURSIVE_H

#include <iostream>
#include <queue>
#include <unordered_map>
#include <cassert>
#include <cstring>
#include <cmath>
#include <cstdint>

#include <gmpxx.h>

#include "digraph.hpp"
#include "set.hpp"
#include "tools.hpp"
#include "cache.hpp"
#include "dot.hpp"



#define CCS_NONE 0
#define CCS_DFS 1
#define CCS_COVERS 2
#define CCS_AUTO 3

#define HUB_NONE 0
#define HUB_BEST 1
#define HUB_PERFECT 2

#define TRANSPOSE_NO 0
#define TRANSPOSE_YES 1
#define TRANSPOSE_AUTO 2
#define TRANSPOSE_AUTO2 3

#define SORT_NO 0
#define SORT_YES 1

#define LATTICE_NONE 0
#define LATTICE_IMAGES 1
#define LATTICE_NODES 2

#define NODE_SHAPE_DEFAULT 0
#define NODE_SHAPE_POINT 1



struct VisualizationOptions
{
	const char *prefix;
	bool subsets;
	int lattice;
	
	int edge_directions;
	int node_shape;
	int subset_box;
	
	VisualizationOptions()
	{
		prefix = NULL;
		subsets = false;
		lattice = LATTICE_NONE;
		edge_directions = true;
		node_shape = NODE_SHAPE_DEFAULT;
		subset_box = false;
	}
};


struct RecursionOptions
{
	int transpose;
	int minimum_transpose_n;
	int ccs;
	double ccs_auto_threshold;
	VisualizationOptions vis;
	// NOTE: disabled
// 	int print_subsets;
// 	int print_cache;
	int hub_splits;
	int sort_topologically;
	bool sampling;
	bool verbose;
	CacheOptions cache;
	
	RecursionOptions()
	{
		transpose = TRANSPOSE_AUTO2;
		minimum_transpose_n = 20;
		ccs = CCS_AUTO;
		ccs_auto_threshold = 3.3;
// 		print_subsets = 0;
// 		print_cache = 0;
		hub_splits = HUB_NONE;
		sort_topologically = SORT_YES;
		verbose = true;
		sampling = false;
	}
};



struct Statistics
{
	long long unsigned int recursive_calls;
	long long unsigned int component_splits;
	long long unsigned int cache_retrievals;
	long long unsigned int evaluated_subgraphs;
	long long unsigned int connected_subgraphs;
	long long unsigned int trivial_cases;
	long long unsigned int minimal_branchings;
	long long unsigned int hub_splits;
	long long unsigned int connectivity_checks;
	
	void init()
	{
		recursive_calls = 0;
		component_splits = 0;
		cache_retrievals = 0;
		evaluated_subgraphs = 0;
		connected_subgraphs = 0;
		trivial_cases = 0;
		minimal_branchings = 0;
		hub_splits = 0;
		connectivity_checks = 0;
	}
	
	void print()
	{
		long long unsigned n_connected = minimal_branchings
		                               + hub_splits;
		
		long long unsigned n_evaluations = trivial_cases
		                                 + component_splits
		                                 + connected_subgraphs;
		
		long long unsigned n_recursive_calls = cache_retrievals
		                                     + evaluated_subgraphs;
		
		assert(trivial_cases == 1);
		assert(n_connected == connected_subgraphs);
		assert(n_evaluations == evaluated_subgraphs);
		assert(n_recursive_calls == recursive_calls);
		
		eprintf("Recursive calls:         %12llu\n", recursive_calls);
		eprintf("- Cache retrievals:      %12llu\n", cache_retrievals);
		eprintf("- Evaluated subposets:   %12llu\n", evaluated_subgraphs);
		eprintf("  - Component splits:    %12llu\n", component_splits);
		eprintf("  - Minimal branchings:  %12llu\n", minimal_branchings);
		eprintf("  - Hub splits:          %12llu\n", hub_splits);
		eprintf("Connectivity checks:     %12llu\n", connectivity_checks);
	}
};



#define SUBPOSET_TYPE_TRIVIAL 0
#define SUBPOSET_TYPE_BRANCH 1
#define SUBPOSET_TYPE_DISCONNECTED 2

template <typename Set>
struct Subposet
{
	int type;
	int n_children;
	Set *children;
	
	~Subposet()
	{
		delete [] children;
	}
};



template <typename N, typename Set>
struct RecInfo
{
	Set *covers;
	unsigned cn;
	
	// search settings
	
	// use maximal downset intersections to determine connected components;
	// this is only guaranteed to work if all traversed subgraphs are upsets,
	// otherwise DFS must be used instead
	int ccs_method;
	
	// use the admissible partitions technique
	// this in particular searches through non-upsets and must use DFS
	int hub_splits;
};



// In some instances we need to convert the number of linear extensions
// to a floating point (unless it is already). By default we convert to
// long double, and the exceptions are specified below.
template <typename N>
struct fp_fallback {
	typedef long double n_orders;
};

// double can be simply kept as a double
template <>
struct fp_fallback<double> {
	typedef double n_orders;
};

// mpz_class is converted to double since conversion to long double is not supported
template <>
struct fp_fallback<mpz_class> {
	typedef double n_orders;
};





template <typename N>
struct RecursiveCounterAuto
{
	virtual ~RecursiveCounterAuto() {}
	
	virtual N count_linex(digraph &poset_) = 0;
	virtual void deallocate() = 0;
	virtual void draw_sample(int *sample) = 0;
	virtual void init_draw_iterative_sample() = 0;
	virtual int draw_next_element() = 0;
	virtual void display_statistics() = 0;
	virtual N get_extension_count() = 0;
	virtual Statistics &get_statistics() = 0;
	virtual std::string short_description() = 0;
	virtual bool is_transposed() = 0;
	virtual void initialize_sampling() = 0;
	virtual void uninitialize_sampling() = 0;
	virtual void set_poset(digraph &poset_) = 0;
	virtual bool transpose_heuristic() const = 0;
	virtual bool transpose_heuristic2() = 0;
};






template <typename N, typename Set>
struct RecursiveCounter : public RecursiveCounterAuto<N>
{
	digraph poset;
	unsigned n;
	
	unsigned *neighbors;
	unsigned *neighbors_n;
	
	Set *neighbor_sets;
	
	Set *predecessors;
	Set *successors;
	
	DOTDigraph *lattice;
	
	Combinatorial<N> comb;
	
	RecursionOptions options;
	
	RecInfo<N, Set> rec;
	
	// these are retained after the search
	Cache<Set, N> *cache;
	Cache<Set, Subposet<Set>*> *subposet_cache;
	Statistics stats;
	
	bool transposed;
	
	N n_extensions;
	
	
	
	void init(unsigned n)
	{
		neighbors = new unsigned[n * n];
		neighbors_n = new unsigned[n];
		neighbor_sets = new Set[n];
		predecessors = new Set[n];
		successors = new Set[n];
	}
	
	void uninit()
	{
		delete [] neighbors;
		delete [] neighbors_n;
		delete [] neighbor_sets;
		delete [] predecessors;
		delete [] successors;
	}
	
	void set_options(RecursionOptions &options_)
	{
		options = options_;
		
		if (options.sampling && options.hub_splits != HUB_NONE) {
			eprintf("  *** Sampling cannot be used together with --hub.\n");
			exit(1);
		}
		
	}
	
	// initializes the recursive linex counter for posets on n elements
	RecursiveCounter(unsigned n, RecursionOptions &options) : poset(n), n(n), comb(n)
	{
		init(n);
		comb.compute_binomials();
		set_options(options);
	}
	
	~RecursiveCounter()
	{
		uninit();
	}
	
	void print_dot(FILE *f, Set X, int mode, int vertex) const
	{
		DOTOptions opt;
		bool subset[n];
		const char *node_colors[n];
		const char *font_colors[n];
		opt.subset = subset;
		opt.node_colors = node_colors;
		opt.font_colors = font_colors;
		
		if (options.vis.node_shape == NODE_SHAPE_POINT) {
			opt.labels = false;
			opt.shape = "point";
		} else {
			opt.shape = "circle";
		}
		
		if (!options.vis.edge_directions) opt.edge_dir = "none";
		
		for (unsigned i = 0; i < n; ++i) {
			subset[i] = X[i];
			if (!X[i]) continue;
			
			node_colors[i] = "#0000ff";
			font_colors[i] = "#ffffff";
			
			if (mode == 2) {
				if (i == (unsigned)vertex) node_colors[i] = "#00dd00";
			} else if (mode == 3) {
				node_colors[i] = "#dddd00";
			}
		}
		
		poset.print_dot(f, opt);
	}
	
	void print_dot(const char *name, Set X, int mode, int vertex) const
	{
		FILE *f = fopen(name, "w");
		print_dot(f, X, mode, vertex);
		fclose(f);
	}
	
	
	void visualize_state(Set subset, int mode=0, int vertex=-1) const
	{
		if (options.vis.lattice) lattice_add_node(subset, mode);
		
		if (options.vis.subsets) {
			int pathlen = strlen(options.vis.prefix) + n + 32;
			char setname[pathlen];
			char filename[pathlen];
			for (unsigned i = 0; i < n; ++i) setname[i] = subset.has(i) ? '1' : '0';
			setname[n] = '\0';
			sprintf(filename, "%s%s.dot", options.vis.prefix, setname);
			
			print_dot(filename, subset, mode, vertex);
		}
	}
	
	void lattice_start()
	{
		lattice = new DOTDigraph();
	}
	
	void lattice_end()
	{
		FILE *f = fopen((std::string(options.vis.prefix) + "lattice.dot").c_str(), "w");
		lattice->print_to_file(f);
		fclose(f);
		
		delete lattice;
	}
	
	void lattice_add_node(Set set, int mode) const
	{
		static std::string mode_colors[4] = {"#0000ff", "#0000ff", "#00dd00", "#dddd00"};
		
		char setname[n+1];
		for (unsigned i = 0; i < n; ++i) setname[i] = set.has(i) ? '1' : '0';
		setname[n] = '\0';
		
		DOTNode *node = lattice->add_node(setname);
		node->add_attribute("label", "");
		
		if (options.vis.lattice == LATTICE_NODES) {
			node->add_attribute("shape", "circle");
		} else if (options.vis.subset_box) {
			node->add_attribute("shape", "rect");
		} else {
			node->add_attribute("shape", "none");
		}
		
		if (options.vis.lattice == LATTICE_IMAGES) {
			node->add_attribute("image", std::string(options.vis.prefix) + setname + ".png");
		} else if (options.vis.lattice == LATTICE_NODES) {
			node->add_attribute("style", "filled");
			node->add_attribute("fillcolor", mode_colors[mode]);
		}
	}
	
	void lattice_add_edge(Set set1, Set set2) const
	{
		char setname1[n+1];
		char setname2[n+1];
		for (unsigned i = 0; i < n; ++i) setname1[i] = set1.has(i) ? '1' : '0';
		setname1[n] = '\0';
		for (unsigned i = 0; i < n; ++i) setname2[i] = set2.has(i) ? '1' : '0';
		setname2[n] = '\0';
		
		DOTEdge *edge = lattice->add_edge(setname2, setname1);
		edge->add_attribute("dir", "back");
	}
	
	
	
	
	
	
	
	void init_subposet(Set subset, int type)
	{
		assert(subposet_cache->count(subset) == 0);
		Subposet<Set> *&subposet = (*subposet_cache)[subset];
		subposet = new Subposet<Set>();
		subposet->type = type;
		subposet->n_children = 0;
		subposet->children = NULL;
	}
	
	void add_subposet_child(Set subset, Set child)
	{
		Subposet<Set> *subposet = (*subposet_cache)[subset];
		
		int new_size = subposet->n_children + 1;
		Set *new_children = new Set[new_size];
		std::copy(subposet->children, subposet->children + subposet->n_children, new_children);
		delete [] subposet->children;
		subposet->children = new_children;
		subposet->n_children = new_size;
		
		subposet->children[new_size - 1] = child;
	}
	
	
	
	
	
	
	
	void child_DFS_nlist(Set &cmp, unsigned i) const
	{
		cmp.set(i);
		
		for (unsigned u = 0; u < n; u++) {
			if (cmp[u]) continue;
			if (poset.has(i, u)) child_DFS_nlist(cmp, u);
		}
	}
	
	Set find_descendants(unsigned i) const
	{
		Set cmp = Set::empty(n);
		child_DFS_nlist(cmp, i);
		return cmp;
	}
	
	void parent_DFS_nlist(Set &cmp, unsigned i) const
	{
		cmp.set(i);
		
		for (unsigned u = 0; u < n; u++) {
			if (cmp[u]) continue;
			if (poset.has(u, i)) parent_DFS_nlist(cmp, u);
		}
	}
	
	Set find_ancestors(unsigned i) const
	{
		Set cmp = Set::empty(n);
		parent_DFS_nlist(cmp, i);
		return cmp;
	}
	
	
	
	
	void compute_neighbor_lists()
	{
		for (unsigned i = 0; i < n; i++) {
			neighbors_n[i] = 0;
			neighbor_sets[i] = Set::empty(n);
		}
		
		for (unsigned i = 0; i < n; i++) {
			for (unsigned j = i+1; j < n; j++) {
				if (poset.has(i,j) || poset.has(j,i)) {
					neighbors[i*n+neighbors_n[i]] = j;
					neighbors[j*n+neighbors_n[j]] = i;
					neighbors_n[i]++;
					neighbors_n[j]++;
					neighbor_sets[i].set(j);
					neighbor_sets[j].set(i);
				}
			}
		}
	}
	
	void leaf_distance_visit(int *level, unsigned i) const
	{
		if (level[i] >= 0) return;
		
		int max_child_level = -1;
		for (unsigned j = 0; j < n; j++) {
			if (!poset.has(i, j)) continue;
			leaf_distance_visit(level, j);
			if (level[j] > max_child_level) max_child_level = level[j];
		}
		
		level[i] = max_child_level + 1;
	}
	
	void root_distance_visit(int *level, unsigned i) const
	{
		if (level[i] >= 0) return;
		
		int max_parent_level = -1;
		for (unsigned j = 0; j < n; j++) {
			if (!poset.has(j, i)) continue;
			root_distance_visit(level, j);
			if (level[j] > max_parent_level) max_parent_level = level[j];
		}
		
		level[i] = max_parent_level + 1;
	}
	
	// returns true iff the heuristic thinks the graph should be tranposed
	bool transpose_heuristic() const
	{
		if (n < 3) return false;
		
		int llevel[n];
		unsigned lcount[n];
		int rlevel[n];
		unsigned rcount[n];
		
		for (unsigned i = 0; i < n; i++) {
			llevel[i] = -1;
			lcount[i] = 0;
			rlevel[i] = -1;
			rcount[i] = 0;
		}
		
		for (unsigned i = 0; i < n; i++) {
			leaf_distance_visit(llevel, i);
			lcount[llevel[i]]++;
			root_distance_visit(rlevel, i);
			rcount[rlevel[i]]++;
		}
		
		unsigned k = n;
		while (lcount[k-1] == 0) k--;
		
		if (options.verbose) eprintf("    Minimal: %i   Maximal: %i\n", rcount[0], lcount[0]);
		bool t = rcount[0] >= lcount[0];
		
		return t;
	}
	
	// returns true iff element i has no predecessor
	bool is_minimal(unsigned i, Set X) const
	{
		for (unsigned j = 0; j < n; j++) {
			if (!X[j]) continue;
			if (poset.has(j, i)) return false;
		}
		
		return true;
	}
	
	Set get_minimal_elements(Set X) const
	{
		Set minimal = Set::empty(n);
		for (unsigned i = 0; i < n; i++) {
			if (!X[i]) continue;
			if (is_minimal(i, X)) minimal.set(i);
		}
		return minimal;
	}
	
	double estimate_hardness_connected(Set subset) const
	{
		Set minimals = get_minimal_elements(subset);
		return pow(2, minimals.cardinality(n)) + estimate_hardness(subset ^ minimals);
	}
	
	double estimate_hardness(Set subset) const
	{
		double estimate = 0.0;
		Set component;
		
		while (!subset.is_empty()) {
			find_connected_component(component, subset);
			estimate += estimate_hardness_connected(component);
			subset ^= component;
		}
		
		return estimate;
	}
	
	// returns true iff the heuristic thinks the graph should be tranposed
	bool transpose_heuristic2()
	{
		compute_neighbor_lists();
		double estimate_no = estimate_hardness(Set::complete(n));
		
		poset.invert();
		compute_neighbor_lists();
		double estimate_yes = estimate_hardness(Set::complete(n));
		
		if (options.verbose) eprintf("    No: %.2f   Yes: %.2f\n", estimate_no, estimate_yes);
		
		poset.invert();
		return estimate_no > estimate_yes;
	}
	
	
	
	// chooses the element with most comparable elements
	int find_central_vertex(Set subset) const
	{
		int d = -1;
		int d_comparable = -1;
		
		Set iter = subset;
		while (!iter.is_empty()) {
			int u = iter.next_element();
		
			Set downset = predecessors[u] & subset;
			Set upset = successors[u] & subset;
			
			Set comparable = downset | upset;
			int u_comparable = comparable.cardinality(n);
			
			if (u_comparable > d_comparable) {
				d = u;
				d_comparable = u_comparable;
			}
		}
		
		return d;
	}
	
	N count_linex_admissible_partitions(Set subset)
	{
		int d = find_central_vertex(subset);
		assert(d != -1);
		
		if (options.vis.prefix) visualize_state(subset, 2, d);
		
		Set A_base = predecessors[d] & subset;
		Set B_base = successors[d] & subset;
		
		Set C = subset ^ d ^ A_base ^ B_base;
		
		stats.hub_splits++;
		
		std::unordered_map<Set, unsigned char *> datas;
		std::queue<Set> ideals;
		
		Set empty = Set::empty(n);
		ideals.push(empty);
		datas[empty] = new unsigned char[n];
		
		for (unsigned v = 0; v < n; ++v) {
			if (!C.has(v)) continue;
			datas[empty][v] = 0;
			for (unsigned u = 0; u < n; ++u) {
				if (!C.has(u)) continue;
				if (poset.has(u, v)) ++datas[empty][v];
			}
		}
		
		N sum_all_partitions = 0;
		
		while (!ideals.empty()) {
			Set X = ideals.front(); ideals.pop();
			
			Set A = A_base | X;
			Set B = B_base | (C ^ X);
			
			if (options.vis.lattice) {
				lattice_add_edge(subset, A);
				lattice_add_edge(subset, B);
			}
			
			N part1 = count_linex_recursive(A);
			N part2 = count_linex_recursive(B);
			
			sum_all_partitions += part1 * part2;
			
			unsigned char *xdata = datas[X];
			
			for (unsigned u = 0; u < n; ++u) {
				if (!C.has(u)) continue;
				if (xdata[u] > 0) continue;
				Set Y = X | u;
				if (Y == X) continue;
				
				if (datas.count(Y) == 0) {
					unsigned char *ydata;
					ydata = new unsigned char[n];
					datas[Y] = ydata;
					ideals.push(Y);
					
					for (unsigned v = 0; v < n; ++v) {
						if (!C.has(v)) continue;
						if (!Y[v]) ydata[v] = xdata[v] - poset.has(u, v);
					}
				}
			}
			
			delete [] xdata;
			datas.erase(X);
		}
		
		return (*cache)[subset] = sum_all_partitions;
	}
	
	int find_hub_vertex(Set subset)
	{
		Set iter = subset;
		while (!iter.is_empty()) {
			int u = iter.next_element();
			
			Set A = predecessors[u] & subset;
			Set B = successors[u] & subset;
			
			if ((subset ^ u ^ A ^ B).is_empty()) return u;
		}
		
		return -1;
	}
	
	N count_linex_hub_split(Set subset)
	{
		int d = find_hub_vertex(subset);
		if (d == -1) return -1;

		Set A = predecessors[d] & subset;
		Set B = successors[d] & subset;
		
		if (options.vis.prefix) visualize_state(subset, 2, d);
		
		if (options.vis.lattice) {
			lattice_add_edge(subset, A);
			lattice_add_edge(subset, B);
		}
		
		stats.hub_splits++;
		
		N part1 = count_linex_recursive(A);
		N part2 = count_linex_recursive(B);
		
		return (*cache)[subset] = part1 * part2;
	}
	
	
	void find_connected_component(Set &component1, Set subset) const
	{
		Set fringe = Set::empty(n) | subset.first_element();
		component1 = Set::empty(n);
		
		for (;;) {
			int u = fringe.first_element();
			if (u == -1) break;
			
			component1.set(u);
			
			fringe = (fringe | (neighbor_sets[u] & subset)) - component1;
		}
	}
	
	N count_linex_connected_components_dfs(Set subset)
	{
		Set component1;
		find_connected_component(component1, subset);
		
		if (component1 == subset) return -1;
		
		if (options.vis.prefix) visualize_state(subset, 3);
		
		
		stats.component_splits++;
		
		Set component2 = subset ^ component1;
		
		if (options.vis.lattice) {
			lattice_add_edge(subset, component1);
			lattice_add_edge(subset, component2);
		}
		
		if (options.sampling) {
			init_subposet(subset, SUBPOSET_TYPE_DISCONNECTED);
			add_subposet_child(subset, component1);
			add_subposet_child(subset, component2);
		}
		
		unsigned n_elements1 = component1.cardinality(n);
		unsigned n_elements2 = component2.cardinality(n);
		
		N cmp1_ext_n = count_linex_recursive(component1, true);
		N cmp2_ext_n = count_linex_recursive(component2);
		
		N bin = comb.binomial(n_elements1 + n_elements2, n_elements2);
		
		return (*cache)[subset] = cmp1_ext_n * cmp2_ext_n * bin;
	}
	
	N count_linex_connected_components_covers(Set subset)
	{
		// project the covers to the subset (and remove empty ones)
		Set combined[rec.cn];
		unsigned k = 0;
		for (unsigned i = 0; i < rec.cn; i++) {
			combined[k] = rec.covers[i];
			combined[k] &= subset;
			if (!combined[k].is_empty()) ++k;
		}
		
		// combine covers that intersect with each other
		for (unsigned i = 0; i < k; i++) {
			bool changes;
			do {
				changes = false;
				for (unsigned j = i+1; j < k; j++) {
					Set cut = combined[i];
					cut &= combined[j];
					if (cut.is_empty()) continue;
					combined[i] |= combined[j];
					k--;
					combined[j] = combined[k];
					j--;
					changes = true;
				}
			} while (changes);
		}
		
		if (k == 1) return -1;
		
		if (options.vis.prefix) visualize_state(subset, 3);
		
		// if not connected, calculate components separately
		stats.component_splits++;
		
		
		N ext_n = 1;
		unsigned rem_n_elements = subset.cardinality(n);
		for (unsigned c = 0; c < k; c++) {
			Set &cmp = combined[c];
			
			unsigned cmp_n_elements = cmp.cardinality(n);
			N bin = comb.binomial(rem_n_elements, cmp_n_elements);
			
			N cmp_ext_n = count_linex_recursive(cmp, true);
			
			ext_n *= cmp_ext_n * bin;
			rem_n_elements -= cmp_n_elements;
			
			if (options.vis.lattice) lattice_add_edge(subset, cmp);
		}
		
		return (*cache)[subset] = ext_n;
	}
	
	
	// branches by removing an element from the graph
	N count_linex_recursive_branch(Set subset)
	{
		// if already cached (checked here to avoid modification to succ)
		if (cache->count(subset) != 0) {
			stats.recursive_calls++;
			stats.cache_retrievals++;
			return (*cache)[subset];
		}
		
		// count extensions of the subgraph
		return count_linex_recursive(subset);
	}
	
	// Counts linear extensions by removing minimal elements.
	// The resulting subproblems are computed by calling count_linex_recursive()
	N count_linex_branching(Set subset)
	{
		if (options.sampling) init_subposet(subset, SUBPOSET_TYPE_BRANCH);
		if (options.vis.prefix) visualize_state(subset, 1);
		
		// branch into subgraphs by removing minimal elements
		stats.minimal_branchings++;
		N sum = 0;
		
		// for all minimal elements v
		Set iter = subset;
		while (!iter.is_empty()) {
			int v = iter.next_element();
			
			// if v has predecessors in the subset, it's not minimal
			if (!(predecessors[v] & subset).is_empty()) continue;
			
			if (options.sampling) add_subposet_child(subset, subset ^ v);
			if (options.vis.lattice) lattice_add_edge(subset, subset ^ v);
			
			sum += count_linex_recursive_branch(subset ^ v);
		}
		
		return (*cache)[subset] = sum;
	}
	
	
	
	// stores and returns the number of linear extensions of a subgraph
	N count_linex_recursive(Set subset, bool known_to_be_connected=false)
	{
// 		if (options.print_subsets) subset.println(n);
		
		stats.recursive_calls++;
		
		// if already computed, return cached value
		if (cache->count(subset) != 0) {
			stats.cache_retrievals++;
			return (*cache)[subset];
		}
		
		stats.evaluated_subgraphs++;
		
		// If the subset is not already cached, we apply one of the
		// following rules to compute the number of posets:
		
		
		// if empty
		
		if (subset.is_empty()) {
			if (options.vis.prefix) visualize_state(subset, 1);
			if (options.sampling) init_subposet(subset, SUBPOSET_TYPE_TRIVIAL);
			stats.trivial_cases++;
			return (*cache)[subset] = 1;
		}
		
		
		// If there are multiple connected components, we count their extensions
		// separately. We can determine the components by looking at the intersections
		// of downsets of maximal elements. If partitioning is used, this doesn't work
		// and we use a depth-first search instead.
		
		if (rec.ccs_method > 0 && !known_to_be_connected) {
			N n_ext_cc;
			stats.connectivity_checks++;
			if (rec.ccs_method == CCS_COVERS) {
				n_ext_cc = count_linex_connected_components_covers(subset);
			} else {
				n_ext_cc = count_linex_connected_components_dfs(subset);
			}
			
			if (n_ext_cc != -1) return n_ext_cc;
		}
		
		stats.connected_subgraphs++;
		
		
		
		// If there are no connected components, we have two alternative methods,
		// admissible partitioning and branching by removing minimal elements.
		// First, we attempt partitioning, if enabled.
		
		if (rec.hub_splits != HUB_NONE) {
			N ext;
			if (rec.hub_splits == HUB_BEST) {
				ext = count_linex_admissible_partitions(subset);
			} else {
				ext = count_linex_hub_split(subset);
			}
			if (ext != -1) return ext;
		}
		
		
		// If partitioning is determined to be too inefficient, we branch by
		// removing minimal elements.
		
		return count_linex_branching(subset);
	}
	
	
	
	
	
	typename fp_fallback<N>::n_orders to_floating_point(long double d)
	{
		return d;
	}
	
	typename fp_fallback<N>::n_orders to_floating_point(double d)
	{
		return d;
	}
	
	typename fp_fallback<N>::n_orders to_floating_point(mpz_class mpz)
	{
		return mpz.get_d();
	}
	
	typename fp_fallback<N>::n_orders to_floating_point(uint64_t d)
	{
		return d;
	}
	
	int draw_sample(int *extension, Set subset)
	{
		Subposet<Set> *subposet = (*subposet_cache)[subset];
		
		if (subposet->type == SUBPOSET_TYPE_BRANCH) {
			
			typename fp_fallback<N>::n_orders r = to_floating_point((*cache)[subset]) * rnd();
			typename fp_fallback<N>::n_orders cumul_sum = 0;
			
			int c;
			for (c = 0; c < subposet->n_children; ++c) {
				cumul_sum += to_floating_point((*cache)[subposet->children[c]]);
				if (cumul_sum > r) break;
			}
			
			assert(c < subposet->n_children);
			
			extension[0] = (subset ^ subposet->children[c]).first_element();
			return draw_sample(extension + 1, subposet->children[c]) + 1;
			
		} else if (subposet->type == SUBPOSET_TYPE_DISCONNECTED) {
			
			assert(subposet->n_children == 2);
			
			// sample random extensions on the disjoint subposets
			int sub_extension1[n];
			int sub_extension2[n];
			int n1 = draw_sample(sub_extension1, subposet->children[0]);
			int n2 = draw_sample(sub_extension2, subposet->children[1]);
			int size = n1 + n2;
			
			// interleave the extensions randomly
			for (int i = 0, i1 = 0, i2 = 0; i < size; ++i) {
				if (rnd() < (double)(n1 - i1) / (size - i1 - i2)) {
					assert(i1 < n1);
					extension[i] = sub_extension1[i1++];
				} else {
					assert(i2 < n2);
					extension[i] = sub_extension2[i2++];
				}
			}
			
			return size;
		}
		
		assert(subposet->type == SUBPOSET_TYPE_TRIVIAL);
		
		return 0;
	}
	
	
	
	// an iterative sampler of linear extensions, computes and returns one element at a time in precedence order
	// this can save time in cases where some decision about the sample can be made based on the first few elements,
	// without having to compute the full linear extension
	struct IterativeSampler
	{
		RecursiveCounter<N, Set> *rc;
		
		Set subset;
		Subposet<Set> *subposet;
		int n;
		
		IterativeSampler *sampler1;
		IterativeSampler *sampler2;
		
		// checks if the subposet is connected
		// if not, the sampler is modified to sample indenpendetly from the components
		void check_connectivity()
		{
			subposet = (*rc->subposet_cache)[subset];
			if (subposet->type != SUBPOSET_TYPE_DISCONNECTED) return;
			
			assert(subposet->n_children == 2);
			
			// initialize subsamplers from the two components
			// NOTE: alternatively could reinitialize this parent as one of the new samplers
			sampler1 = rc->new_iterative_sampler();
			sampler2 = rc->new_iterative_sampler();
			sampler1->init(rc, subposet->children[0]);
			sampler2->init(rc, subposet->children[1]);
		}
		
		void init(RecursiveCounter<N, Set> *rc_, Set subset_)
		{
			rc = rc_;
			subset = subset_;
			
			n = subset.cardinality(rc->n);
			
			// check if the subposet is connected
			check_connectivity();
		}
		
		int next_element()
		{
			// if this subposet is disconnected
			if (subposet->type == SUBPOSET_TYPE_DISCONNECTED) {
				// choose one of the two components and draw one element from it (then decrement the size of this subposet)
				IterativeSampler *sampler = rnd() < (double)sampler1->n / n ? sampler1 : sampler2;
				--n;
				return sampler->next_element();
			}
			
			// if connected, choose one minimal element
			typename fp_fallback<N>::n_orders r = rc->to_floating_point((*rc->cache)[subset]) * rnd();
			typename fp_fallback<N>::n_orders cumul_sum = 0;
			
			int c;
			for (c = 0; c < subposet->n_children; ++c) {
				cumul_sum += rc->to_floating_point((*rc->cache)[subposet->children[c]]);
				if (cumul_sum > r) break;
			}
			
			int u = (subset ^ subposet->children[c]).first_element();
			
			// update the state of this sampler by removing the sampled element from the subset,
			// then check if this makes the subposet disconnected
			subset = subposet->children[c];
			--n;
			check_connectivity();
			
			return u;
		}
	};
	
	IterativeSampler *iterative_sampler_root;
	IterativeSampler *iterative_sampler_next;
	
	IterativeSampler *new_iterative_sampler()
	{
// 		assert(iterative_sampler_next - iterative_sampler_root < 2*n+1);
		return iterative_sampler_next++;
	}
	
	// initializes (iterative) sampling, must be called before any calls to init_draw_iterative_sample()
	void initialize_sampling()
	{
		if (transposed) return;
		
		// NOTE: probably overestimated, 2n+1 might be enough, or n+1 with implementation change
		iterative_sampler_root = new IterativeSampler[n * n + 1];
	}
	
	// uninitializes (iterative) sampling
	void uninitialize_sampling()
	{
		if (transposed) return;
		
		delete [] iterative_sampler_root;
	}
	
	// initializes iterative drawing of a sample, must be called once before each sample,
	// then followed by a number of calls to draw_iterative_sample()
	void init_draw_iterative_sample()
	{
		assert(!transposed);
		
		Set subset = Set::complete(n);
		iterative_sampler_next = iterative_sampler_root + 1;
		iterative_sampler_root->init(this, subset);
	}
	
	// draws the next element in a sample
	// this can be called up to n times after each call to init_draw_iterative_sample()
	int draw_next_element()
	{
		return iterative_sampler_root->next_element();
	}
	
	
	
	void draw_sample(int *sample)
	{
		Set subset = Set::complete(n);
		
		int extension[n];
		draw_sample(extension, subset);
		
		if (transposed) {
			for (unsigned i = 0; i < n; ++i) sample[i] = extension[n-i-1];
		} else {
			for (unsigned i = 0; i < n; ++i) sample[i] = extension[i];
		}
	}
	
	
	
	void print_cache(Cache<Set, N> *cache) const
	{
		for (unsigned i = 0; i < cache->used; ++i) {
			Item<Set, N> &item = cache->array[i];
			item.key.print(n);
			printf(" ");
			std::cout << item.data << std::endl;
		}
	}
	
	bool decide_transpose(RecursionOptions &options)
	{
		if (options.transpose == TRANSPOSE_NO) return false;
		if (options.transpose == TRANSPOSE_YES) return true;
		
		if (n < options.minimum_transpose_n) return false;
		
		bool heuristic;
		
		if (options.transpose == TRANSPOSE_AUTO) {
			if (options.verbose) eprintf("  Running simple transpose heuristic...\n");
			heuristic = transpose_heuristic();
		} else {
			if (options.verbose) eprintf("  Running advanced transpose heuristic...\n");
			heuristic = transpose_heuristic2();
		}
		
		if (heuristic) {
			if (options.verbose) eprintf("    Transpose predicted easier.\n");
			return true;
		}
		
		if (options.verbose) eprintf("    Transpose predicted harder.\n");
		return false;
	}
	
	void decide_ccs(double avg_deg)
	{
		// decide the CCS method
		if (options.sampling || options.hub_splits != HUB_NONE) {
			// sampling and hub splits require either DFS or none; AUTO selects DFS, anything else raises an error
			if (options.ccs == CCS_DFS || options.ccs == CCS_NONE) {
				rec.ccs_method = options.ccs;
			} else if (options.ccs == CCS_AUTO) {
				if (options.verbose) eprintf("  Autoselecting DFS.\n");
				rec.ccs_method = CCS_DFS;
			} else {
				eprintf("  *** Sampling and hub splits require either --ccs=dfs or --ccs=none.\n");
				exit(1);
			}
		} else {
			// otherwise auto uses the average degree heuristic to choose between DFS and covers
			if (options.ccs == CCS_AUTO) {
				if (options.verbose) eprintf("  Choosing CCS method by average degree: %f\n", avg_deg);
				
				if (avg_deg < options.ccs_auto_threshold) {
					rec.ccs_method = CCS_DFS;
				} else {
					rec.ccs_method = CCS_COVERS;
				}
			} else {
				rec.ccs_method = options.ccs;
			}
		}
		
		if (options.verbose) {
			eprintf("CCS method: ");
			if (rec.ccs_method == CCS_NONE) {
				eprintf("None\n");
			} else if (rec.ccs_method == CCS_COVERS) {
				eprintf("Maximal covers\n");
			} else if (rec.ccs_method == CCS_DFS) {
				eprintf("Depth-first search\n");
			}
		}
	}
	
	
	
	void set_poset(digraph &poset_)
	{
		assert(poset_.n == n);
		poset.make_copy(poset_);
	}
		
	N count_linex(digraph &poset_)
	{
		set_poset(poset_);
		stats.init();
		
		if (options.verbose) eprintf("Preprocessing...\n");
		
		double avg_deg = poset.average_degree();
		
		// transitive reduction
		if (options.verbose) eprintf("  Computing transitive reduction...\n");
		poset.remove_transitive_arcs();
		
		if (options.sort_topologically) {
			if (options.verbose) eprintf("  Sorting topologically...\n");
			poset.sort_topologically();
		}
		
		bool transpose = decide_transpose(options);
		
		if (transpose) {
			if (options.verbose) eprintf("  Transposing...\n");
			poset.invert();
			transposed = true;
		} else {
			transposed = false;
		}
		
		Set subset = Set::complete(n);
		
		compute_neighbor_lists();
		
		
		
		decide_ccs(avg_deg);
		
		if (options.hub_splits == HUB_BEST && n > 256) {
			eprintf("  *** --hub=best is usable only up 256 elements.\n");
			exit(1);
		}
		
		rec.hub_splits = options.hub_splits;
		
		Set covers[n];
		
		if (rec.ccs_method == CCS_COVERS) {
			if (options.verbose) eprintf("  Finding maximal elements...\n");
			
			unsigned maximals[n];
			
			rec.cn = poset.get_maximal_elements(maximals);
			
			if (options.verbose) eprintf("    Maximal elements: %i\n", rec.cn);
			if (options.verbose) eprintf("  Building covers...\n");
			
			rec.covers = covers;
			for (unsigned i = 0; i < rec.cn; i++) {
				rec.covers[i] = find_ancestors(maximals[i]);
			}
		}
		
		if (options.verbose) eprintf("  Building downsets and upsets...\n");
		
		for (unsigned u = 0; u < n; u++) {
			predecessors[u] = find_ancestors(u) ^ u;
			successors[u] = find_descendants(u) ^ u;
		}
		
		cache = new Cache<Set, N>(options.cache, NULL);
		
		if (options.sampling) {
			subposet_cache = new Cache<Set, Subposet<Set>*>(options.cache, NULL);
		}
		
		
		n_extensions = -1;
		
		if (options.verbose) eprintf("Solving...\n");
		if (options.vis.lattice) lattice_start();
		
		n_extensions = count_linex_recursive(subset);
		
		if (options.vis.lattice) lattice_end();
		if (options.verbose) eprintf("Solved. Deallocating...\n");
// 		if (options.print_cache) print_cache(cache);
		
		return n_extensions;
	}
	
	void display_statistics()
	{
		stats.print();
	}
	
	Statistics &get_statistics()
	{
		return stats;
	}
	
	void deallocate()
	{
		delete cache;
		
		if (options.sampling) {
			subposet_cache->delete_data();
			delete subposet_cache;
		}
	}
	
	N get_extension_count()
	{
		return n_extensions;
	}
	
	std::string short_description()
	{
		if (transposed) return std::string("DP (t)");
		return std::string("DP");
	}
	
	bool is_transposed()
	{
		return transposed;
	}
};




template <typename N>
RecursiveCounterAuto<N> *get_recursive_counter(unsigned n, RecursionOptions &options)
{
	// For small posets we use 32 or 64 bit integers to represent subsets of posets.
	// For larger posets we use the fixed_set class, precompiled for various sizes.
	// Simply add more lines here if you need fixed_set for larger posets.
	if      (n <= 32)  return new RecursiveCounter<N, set32>(n, options);
	else if (n <= 64)  return new RecursiveCounter<N, set64>(n, options);
	else if (n <= 128) return new RecursiveCounter<N, fixed_set<2>>(n, options);
	else if (n <= 256) return new RecursiveCounter<N, fixed_set<4>>(n, options);
	else if (n <= 512) return new RecursiveCounter<N, fixed_set<8>>(n, options);
	
	// To support posets up to size X = 64 * Y, add a line
//	else if (n <= X) return new RecursiveCounter<N, fixed_set<Y>>(n, options);
	
	// If the largest precompiled fixed_set is not large enough, we use dynamic_set,
	// which is in general much slower than a sufficiently large fixed_set would be.
	return new RecursiveCounter<N, dynamic_set>(n, options);
}



#endif
