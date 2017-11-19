#ifndef VEIE_H
#define VEIE_H

#include <cassert>

#include "digraph.hpp"
#include "tools.hpp"



struct VEIEOptions
{
	const char *elim_order_file;
	bool verbose;
	int n_random_orders;
	
	VEIEOptions()
	{
		elim_order_file = NULL;
		verbose = true;
		n_random_orders = 1000;
	}
};



struct OrderHeuristic
{
	digraph &poset;
	
	OrderHeuristic(digraph &poset) : poset(poset) {}
	virtual ~OrderHeuristic() {}
	
	virtual std::string name() const = 0;
	virtual void get_order(int *order) const = 0;
};

// use the default order of the elements
struct DefaultOrderHeuristic : OrderHeuristic
{
	DefaultOrderHeuristic(digraph &poset) : OrderHeuristic(poset) {}
	
	std::string  name() const
	{
		return std::string ("Default order heuristic");
	}
	
	void get_order(int *order) const
	{
		for (unsigned i = 0; i < poset.n; ++i) {
			order[i] = i;
		}
	}
};

// use the order induced by repeatedly removing elements with fewest neighors
struct LeastDegreeHeuristic : OrderHeuristic
{
	LeastDegreeHeuristic(digraph &poset) : OrderHeuristic(poset) {}
	
	std::string name() const
	{
		return std::string ("Least degree heuristic");
	}
	
	// degree of u in the elimination graph
	int degree(const digraph &eg, const bool *eliminated, unsigned u) const
	{
		int deg = 0;
		for (unsigned v = 0; v < eg.n; ++v) {
			if (!eliminated[v]) deg += eg.has(u, v);
		}
		return deg;
	}
	
	// choose next node to eliminate
	int get_next(const digraph &eg, const bool *eliminated) const
	{
		int v = -1;
		int vdeg = eg.n;
		
		for (unsigned u = 0; u < eg.n; ++u) {
			if (eliminated[u]) continue;
			int udeg = degree(eg, eliminated, u);
			if (udeg >= vdeg) continue;
			
			v = u;
			vdeg = udeg;
		}
		
		return v;
	}
	
	void get_order(int *order) const
	{
		// initialize an elimination graph (eg)
		unsigned n = poset.n;
		bool eliminated[n];
		digraph eg(n);
		
		for (unsigned u = 0; u < n; ++u) {
			eliminated[u] = false;
			for (unsigned v = 0; v < n; ++v) {
				eg.pair(u,v) = poset.has(u,v) || poset.has(v,u);
			}
		}
		
		// eliminate nodes one by one
		for (unsigned k = 0; k < n; ++k) {
			int u = order[k] = get_next(eg, eliminated);
			eliminated[u] = true;
			
			// connect all neighbors of the eliminated node to each other
			for (unsigned v = 0; v < n; ++v) {
				if (eliminated[v] || !eg.has(u, v)) continue;
				
				for (unsigned w = v+1; w < n; ++w) {
					if (eliminated[w] || !eg.has(u, w)) continue;
					eg.add(v, w);
					eg.add(w, v);
				}
			}
		}
	}
};

// use an order induced by a depth-first traversal
struct TreeHeuristic : OrderHeuristic
{
	TreeHeuristic(digraph &poset) : OrderHeuristic(poset) {}
	
	std::string name() const
	{
		return std::string("Tree heuristic");
	}
	
	void tree_sort_visit(std::vector<unsigned> &sorted, bool *visited, unsigned i) const
	{
		if (visited[i]) return;
		visited[i] = true;
		
		for (unsigned j = 0; j < poset.n; j++) {
			if (poset.has(j, i) || poset.has(i, j)) tree_sort_visit(sorted, visited, j);
		}
		
		sorted.push_back(i);
	}
	
	void get_order(int *order) const
	{
		std::vector<unsigned> sorted;
		bool visited[poset.n];
		
		for (unsigned i = 0; i < poset.n; i++) {
			visited[i] = false;
		}
		
		for (unsigned i = 0; i < poset.n; i++) {
			tree_sort_visit(sorted, visited, i);
		}
		
		for (unsigned i = 0; i < poset.n; i++) {
			order[i] = sorted[i];
		}
	}
};

// sample random orders and choose the best one
struct RandomHeuristic : OrderHeuristic
{
	int n_iterations;
	
	RandomHeuristic(digraph &poset, int n_iterations) : OrderHeuristic(poset), n_iterations(n_iterations) {}
	
	std::string name() const
	{
		return "Random heuristic (" + std::to_string(n_iterations) + ")";
	}
	
	void get_order(int *order) const
	{
		int sample[poset.n];
		int best_width = poset.n;
		
		for (unsigned i = 0; i < poset.n; ++i) {
			sample[i] = i;
		}
		
		for (int k = 0; k < n_iterations; ++k) {
			// permute the elements
			for (int i = 0; i < (int)poset.n-1; ++i) {
				int j = i + rand() % (poset.n - i);
				int swap = sample[i];
				sample[i] = sample[j];
				sample[j] = swap;
			}
			
			// check the induced width
			int width = poset.find_induced_width(sample);
			if (width < best_width) {
				best_width = width;
				for (unsigned i = 0; i < poset.n; ++i) {
					order[i] = sample[i];
				}
			}
		}
	}
};



template <typename N>
struct VEIE
{
	digraph poset;
	unsigned n;
	
	Combinatorial<N> comb;
	
	VEIEOptions options;
	
	VEIE(unsigned n, VEIEOptions &options) : poset(n), n(n), comb(n), options(options)
	{
		comb.compute_binomials();
	}
	
	~VEIE() {}
	
	
	struct Function
	{
		int domain_size;
		int range;
		int *domain;
		int n_values;
		N *values;
		
		int range_pow(int k)
		{
			int nn = 1;
			for (int i = 0; i < k; ++i) nn *= range;
			return nn;
		}
		
		Function(int domain_size, int range) : domain_size(domain_size), range(range)
		{
			domain = new int[domain_size];
			n_values = range_pow(domain_size);
			values = new N[n_values];
		}
		
		~Function()
		{
			delete [] domain;
			delete [] values;
		}
		
		int index(int *key)
		{
			int index = 0;
			int factor = 1;
			
			for (int i = 0; i < domain_size; ++i) {
				index += factor * key[i];
				factor *= range;
			}
			
			return index;
		}
		
		N &value(int *key)
		{
			return values[index(key)];
		}
		
		int map_index(int *key)
		{
			int index = 0;
			int factor = 1;
			
			for (int i = 0; i < domain_size; ++i) {
				index += factor * key[domain[i]];
				factor *= range;
			}
			
			return index;
		}
		
		N &map_value(int *map)
		{
			return values[map_index(map)];
		}
		
		bool domain_has(int k)
		{
			for (int i = 0; i < domain_size; ++i) {
				if (domain[i] == k) return true;
			}
			
			return false;
		}
		
		void print_domain()
		{
			if (domain_size == 0) {
				printf("Ø");
				return;
			}
			
			for (int i = 0; i < domain_size-1; ++i) printf("%i,", domain[i]);
			printf("%i", domain[domain_size-1]);
		}
		
		void print()
		{
			int nn = range_pow(domain_size);
			
			for (int i = 0; i < domain_size; ++i) printf("%i ", domain[i]);
			printf("\n");
			
			int *args = new int[domain_size+1];
			
			for (int i = 0; i < domain_size; ++i) args[i] = 0;
			
			for (int i = 0; i < nn; ++i) {
				for (int j = 0; j < domain_size; ++j) printf("%i ", args[j]);
				printf("  %i\n", value(args));
				
				int k = 0;
				while (++args[k] == range) args[k++] = 0;
			}
			
			delete [] args;
		}
	};
	
	
	Function *make_arc_function(int u, int v, int range)
	{
		Function *f = new Function(2, range);
		f->domain[0] = u;
		f->domain[1] = v;
		
		int args[2];
		int &i = args[0];
		int &j = args[1];
		
		for (i = 0; i < range; ++i) {
			for (j = 0; j < range; ++j) {
				f->value(args) = i < j ? 1 : 0;
			}
		}
		
		return f;
	}
	
	Function *eliminate(int var, Function **functions, int n_functions, int range)
	{
		int n_variables = n;
		
		bool in_domain[n_variables];
		for (int i = 0; i < n_variables; ++i) in_domain[i] = false;
		
		for (int i = 0; i < n_functions; ++i) {
			for (int j = 0; j < functions[i]->domain_size; ++j) {
				in_domain[functions[i]->domain[j]] = true;
			}
		}
		
		int domain_size = 0;
		for (int i = 0; i < n_variables; ++i) {
			if (in_domain[i]) ++domain_size;
		}
		
		int *domain = new int[domain_size];
		
		int j = 0;
		for (int i = 0; i < n_variables; ++i) {
			if (in_domain[i]) domain[j++] = i;
		}
		
		Function *f = new Function(domain_size - 1, range);
		
		j = 0;
		for (int i = 0; i < n_variables; ++i) {
			if (in_domain[i] && i != var) f->domain[j++] = i;
		}
		
		int *key = new int[n_variables];
		for (int i = 0; i < domain_size; ++i) key[domain[i]] = 0;
		
		delete [] domain;
		
		for (int i = 0; i < f->n_values; ++i) {
			N sum = 0;
			
			for (key[var] = 0; key[var] < range; ++key[var]) {
				N product = 1;
				for (int j = 0; j < n_functions; ++j) product *= functions[j]->map_value(key);
				sum += product;
			}
			
			f->map_value(key) = sum;
			
			int k = 0;
			while (k < f->domain_size && ++key[f->domain[k]] == range) {
				key[f->domain[k++]] = 0;
			}
		}
		
		delete [] key;
		
		return f;
	}
	
	
	N elimination(Function **functions, int n_functions, int *elim_order, int range)
	{
		int n_variables = n;
		
		Function **remainder = new Function*[n_functions];
		Function **bucket = new Function*[n_functions];
		
		for (int i = 0; i < n_functions; ++i) {
			remainder[i] = functions[i];
		}
		
		N empty = 1;
		
		for (int k = 0; k < n_variables; ++k) {
			if (options.verbose) eprintf(".");
			int var = elim_order[k];
			int bucket_size = 0;
			int remainder_size = 0;
			for (int i = 0; i < n_functions; ++i) {
				if (remainder[i]->domain_has(var)) {
					bucket[bucket_size++] = remainder[i];
				} else {
					remainder[remainder_size++] = remainder[i];
				}
			}
			
			if (bucket_size == 0) {
				empty *= range;
				continue;
			}
			
			Function *f = eliminate(var, bucket, bucket_size, range);
			
			for (int i = 0; i < bucket_size; ++i) {
				delete bucket[i];
			}
			
			remainder[remainder_size++] = f;
			n_functions = remainder_size;
		}
		
		for (int i = 0; i < n_functions; ++i) {
			assert(remainder[i]->domain_size == 0);
			empty *= remainder[i]->value(NULL);
		}
		
		for (int i = 0; i < n_functions; ++i) {
			delete remainder[i];
		}
		
		delete [] remainder;
		delete [] bucket;
		
		return empty;
	}
	
	void choose_elimination_order(int *elim_order)
	{
		if (options.elim_order_file != NULL) {
			if (options.verbose) eprintf("  Reading elimination order from file\n");
			FILE *f = fopen(options.elim_order_file, "r");
			assert(f);
			for (unsigned i = 0; i < n; ++i) assert(fscanf(f, "%i", &elim_order[i]) == 1);
			fclose(f);
			return;
		}
		
		if (options.verbose) eprintf("  Running elimination order heuristics:\n");
		
		std::vector<OrderHeuristic*> heuristics;
		heuristics.push_back(new DefaultOrderHeuristic(poset));
		heuristics.push_back(new LeastDegreeHeuristic(poset));
		heuristics.push_back(new TreeHeuristic(poset));
		heuristics.push_back(new RandomHeuristic(poset, options.n_random_orders));
		
		int order[n];
		int best_width = n;
		
		for (int i = 0; i < heuristics.size(); ++i) {
			heuristics[i]->get_order(order);
			int induced_width = poset.find_induced_width(order);
			
			if (options.verbose) eprintf("    %-30s: %i\n", heuristics[i]->name().c_str(), induced_width);
			
			if (induced_width < best_width) {
				best_width = induced_width;
				for (int j = 0; j < n; ++j) elim_order[j] = order[j];
			}
			
			delete heuristics[i];
		}
	}
	
	N count_extensions(digraph &poset_)
	{
		if (options.verbose) eprintf("Running VEIE...\n");
		
		poset.make_copy(poset_);
		
		// take the transitive reduction (cover graph)
		poset.remove_transitive_arcs();
		
		int elim_order[n];
		choose_elimination_order(elim_order);
		
		if (options.verbose) {
			eprintf("  Elimination order:");
			for (unsigned i = 0; i < n; ++i) {
				eprintf(" %i", elim_order[i]);
			}
			eprintf("\n");
			
			eprintf("  Induced width: %i\n", poset.find_induced_width(elim_order));
		}
		
		// perform the inclusion-exclusion steps
		N ext = 0;
		N sign = n % 2 == 0 ? 1 : -1;
		for (unsigned k = 0; k <= n; ++k) {
			if (options.verbose) eprintf("Round %i ", k);
			
			// create the initial set of functions, one for each arc in the cover graph
			Function *functions[n*n];
			int n_functions = 0;
			for (unsigned u = 0; u < n; ++u) {
				for (unsigned v = 0; v < n; ++v) {
					if (poset.has(u, v)) functions[n_functions++] = make_arc_function(u, v, k);
				}
			}
			
			N elim = elimination(functions, n_functions, elim_order, k);
			N binom = comb.binomial(n, k);
			
			ext += sign * binom * elim;
			sign = -sign;
			if (options.verbose) eprintf("\n");
		}
		
		return ext;
	}
};


#endif
