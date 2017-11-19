#ifndef DIGRAPH_H
#define DIGRAPH_H

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "tools.hpp"
#include "dot.hpp"


struct DOTOptions
{
	bool *subset;
	const char *shape;
	const char **node_colors;
	const char **node_edge_colors;
	bool labels;
	const char *font_size;
	const char **font_colors;
	const char *edge_dir;
	
	DOTOptions()
	{
		subset = NULL;
		shape = NULL;
		node_colors = NULL;
		node_edge_colors = NULL;
		labels = true;
		font_size = NULL;
		font_colors = NULL;
		edge_dir = "back";
	}
};



struct digraph
{
	unsigned n;
	bool *rel;
	
	void init(const bool *matrix)
	{
		for (unsigned i = 0; i < n; ++i) {
			for (unsigned j = 0; j < n; ++j) {
				pair(i, j) = matrix[i * n + j];
			}
		}
	}
	
	void make_empty()
	{
		for (unsigned i = 0; i < n; ++i) {
			for (unsigned j = 0; j < n; ++j) {
				pair(i, j) = false;
			}
		}
	}
	
	void make_copy(const digraph &dg)
	{
		init(dg.rel);
	}
	
	void make_subgraph(const digraph &dg, const int *elements)
	{
		for (unsigned u = 0; u < n; ++u) {
			for (unsigned v = 0; v < n; ++v) {
				pair(u, v) = dg(elements[u], elements[v]);
			}
		}
	}
	
	// constructs an uninitialized digraph of n elements
	digraph(unsigned n) : n(n), rel(new bool[n*n]) {}
	
	// constructs a digraph of n elements from given boolean matrix
	digraph(unsigned n, const bool *matrix) : n(n), rel(new bool[n*n])
	{
		init(matrix);
	}
	
	// constructs a copy of given digraph
	digraph(const digraph &dg) : n(dg.n), rel(new bool[n*n])
	{
		make_copy(dg);
	}
	
	// constructs an n-element subgrpah of given digraph
	// 'elements' is the list of elements in the given digraph to include in the subgraph
	digraph(unsigned n, const digraph &dg, const int *elements) : n(n), rel(new bool[n*n])
	{
		make_subgraph(dg, elements);
	}
	
	digraph(const char *fileName)
	{
		std::ifstream file(fileName);
		if (!file) {
			fprintf(stderr, "Error: Could not open file \"%s\".\n", fileName);
			exit(1);
		}
		std::string row;
		n = 0;
		
		// read the first row (and get the number of variables)
		std::vector<bool> firstRow(32);
		getline(file, row);
		std::istringstream rowStream(row);
		rowStream >> std::ws;
		while(!rowStream.eof()) {
			if (n >= firstRow.size())
				firstRow.resize(firstRow.size() * 2);
			int tmp;
			rowStream >> tmp;
			if (rowStream.fail()) {
				fprintf(stderr, "Error: Could not read the value on row 1 column %d.\n", n+1);
				exit(1);
			}
			firstRow[n] = (bool) tmp;
			rowStream >> std::ws;
			++n;
		}
		file >> std::ws;
		
		rel = new bool[n*n];
		
		// create the matrix and fill first row
		for (unsigned j = 0; j < n; ++j)
			if (firstRow[j]) {
				pair(0, j) = 1;
			} else {
				pair(0, j) = 0;
			}
		
		// load the data rest of the matrix
		for (unsigned i = 1; i < n; ++i) {
			getline(file, row);
			std::istringstream rowStream(row);
			for (unsigned j = 0; j < n; ++j) {
				int tmp;
				rowStream >> tmp;
				if (rowStream.fail()) {
					fprintf(stderr, "Error: Could not read the %dth value on row %d.\n", n+1, i+1);
					exit(1);
				}
				if (tmp) {
					pair(i, j) = 1;
				} else {
					pair(i, j) = 0;
				}
			}
			file >> std::ws;
		}
	}
	
	~digraph()
	{
		delete [] rel;
	}
	
	bool& pair(unsigned i, unsigned j)
	{
		return rel[i * n + j];
	}
	
	void flip(unsigned i, unsigned j)
	{
		pair(i, j) ^= 1;
	}
	
	void add(unsigned i, unsigned j)
	{
		pair(i, j) = 1;
	}
	
	void del(unsigned i, unsigned j)
	{
		pair(i, j) = 0;
	}
	
	bool has(unsigned i, unsigned j) const
	{
		return (*this)(i, j);
	}
	
	bool operator() (unsigned i, unsigned j) const
	{
		return rel[i * n + j];
	}
	
	bool comparable(unsigned i, unsigned j) const
	{
		return has(i, j) || has(j, i);
	}
	
	
	void print() const
	{
		for (unsigned j = 0; j < n; j++) {
			for (unsigned i = 0; i < n; i++) {
				printf("%s", has(i, j) ? "#" : "-");
				if (i < n - 1) printf(" ");
			}
			printf("\n");
		}
	}
	
	double average_degree() const
	{
		int edges = 0;
		for (unsigned i = 0; i < n; i++) {
			for (unsigned j = 0; j < n; j++) {
				if (has(i, j)) edges++;
			}
		}
		
		return 2.0 * edges / n;
	}
	
	// inverts the graph (in-place transpose)
	void invert()
	{
		for (unsigned i = 0; i < n-1; i++) {
			for (unsigned j = i+1; j < n; j++) {
				if (has(i, j) != has(j, i)) {
					flip(i, j);
					flip(j, i);
				}
			}
		}
	}
	
	void topological_sort_visit(std::vector<unsigned> &sorted, bool *visited, unsigned i) const
	{
		if (visited[i]) return;
		
		for (unsigned j = 0; j < n; j++) {
			if (has(j, i)) topological_sort_visit(sorted, visited, j);
		}
		
		visited[i] = true;
		sorted.push_back(i);
	}
	
	void find_topological_ordering(std::vector<unsigned> &sorted) const
	{
		bool visited[n];
		
		for (unsigned i = 0; i < n; i++) {
			visited[i] = false;
		}
		
		for (unsigned i = 0; i < n; i++) {
			if (!is_maximal(i)) continue;
			topological_sort_visit(sorted, visited, i);
		}
	}
	
	// sort the graph topologically in place
	void sort_topologically()
	{
		std::vector<unsigned> sorted;
		find_topological_ordering(sorted);
		
		bool new_rel[n * n];
		
		for (unsigned i = 0; i < n; i++) {
			for (unsigned j = 0; j < n; j++) {
				new_rel[i * n + j] = has(sorted[i], sorted[j]);
			}
		}
		
		for (unsigned i = 0; i < n; i++) {
			for (unsigned j = 0; j < n; j++) {
				pair(i,j) = new_rel[i * n + j];
			}
		}
	}
	
	// computes longest paths from node u to all other nodes recursively
	// visited keeps track of nodes already traversed to avoid recomputation
	// lpath is the array of found path lengths, -1 denoting no path
	void longest_paths(bool *visited, int *lpath, unsigned u) const
	{
		if (visited[u]) return;
		
		// for each successor v of u
		for (unsigned v = 0; v < n; v++) {
			if (!has(u, v)) continue;
			
			// first compute paths from v recursively
			longest_paths(visited, lpath, v);
			
			// then for all w that are reachable from v
			for (unsigned w = 0; w < n; w++) {
				if (lpath[v*n+w] == -1) continue;
				
				// update the distance from u to w
				int len = lpath[v*n+w] + 1;
				if (len > lpath[u*n+w]) lpath[u*n+w] = len;
			}
		}
		
		visited[u] = true;
	}
	
	void longest_paths(int *lpath) const
	{
		bool visited[n];
		
		for (unsigned u = 0; u < n; u++) {
			visited[u] = false;
			for (unsigned v = 0; v < n; v++) {
				lpath[u*n+v] = (u == v) ? 0 : -1;
			}
		}
		
		for (unsigned u = 0; u < n; u++) {
			longest_paths(visited, lpath, u);
		}
	}
	
	// as above but in-place
	void remove_transitive_arcs()
	{
		// find longest paths between all nodes
		int lpath[n*n];
		longest_paths(lpath);
		
		// include edges between nodes for which the longest path has length 1
		for (unsigned u = 0; u < n; u++) {
			for (unsigned v = 0; v < n; v++) {
				pair(u, v) = lpath[u*n+v] == 1;
			}
		}
	}
	
	void take_transitive_closure()
	{
		// find longest paths between all nodes
		int lpath[n*n];
		longest_paths(lpath);
		
		for (unsigned u = 0; u < n; u++) {
			for (unsigned v = 0; v < n; v++) {
				pair(u, v) = lpath[u*n+v] >= 1;
			}
		}
	}
	
	
	// returns true iff element u has no successor/child
	bool is_maximal(unsigned u) const
	{
		for (unsigned v = 0; v < n; v++) {
			if (has(u, v)) return false;
		}
		
		return true;
	}
	
	unsigned get_maximal_elements(unsigned *maximal) const
	{
		unsigned k = 0;
		for (unsigned u = 0; u < n; u++) {
			if (is_maximal(u)) {
				maximal[k] = u;
				k++;
			}
		}
		return k;
	}
	
	
	// returns the parents of i within set X
	template <typename Set>
	Set get_parents(unsigned i, Set X) const
	{
		Set ps = Set::empty(n);
		for (unsigned k = 0; k < n; k++) {
			if (!X[k]) continue;
			if (has(k, i)) ps.set(k);
		}
		return ps;
	}
	
	// returns the children of i within set X
	template <typename Set>
	Set get_children(unsigned i, Set X) const
	{
		Set cs = Set::empty(n);
		for (unsigned k = 0; k < n; k++) {
			if (!X[k]) continue;
			if (has(i, k)) cs.set(k);
		}
		return cs;
	}
	
	// returns the parents (predecessors) of u
	template <typename Set>
	Set get_parents(unsigned u) const
	{
		Set parents = Set::empty(n);
		for (unsigned v = 0; v < n; v++) {
			if (has(v, u)) parents.set(v);
		}
		return parents;
	}
	
	
	
	int find_induced_width(int *order) const
	{
		// undirected cover graph
		bool graph[n * n];
		for (unsigned i = 0; i < n; ++i) {
			for (unsigned j = i; j < n; ++j) {
				graph[j*n+i] = graph[i*n+j] = has(i,j) || has(j,i);
			}
		}
		
		int induced_width = 0;
		
		for (unsigned k = 0; k < n; ++k) {
			// vertex to be eliminated
			int e = order[k];
			
			int width = 0;
			for (unsigned i = 0; i < n; ++i) {
				if (graph[e*n+i]) ++width;
			}
			if (width > induced_width) induced_width = width;
			
			// search through all pairs of neighbors of e
			for (unsigned i = 0; i < n-1; ++i) {
				if (!graph[e*n+i]) continue;
				for (unsigned j = i+1; j < n; ++j) {
					if (i == j) continue;
					if (!graph[e*n+j]) continue;
					
					// add an edge between the neighbors
					graph[i*n+j] = graph[j*n+i] = 1;
				}
			}
			
			// remove e from the graph
			for (unsigned i = 0; i < n; ++i) {
				graph[e*n+i] = graph[i*n+e] = 0;
			}
		}
		
		return induced_width+1;
	}
	
	void print_dot(FILE *f, DOTOptions &options) const
	{
		DOTDigraph dot;
		
		// nodes
		for (unsigned i = 0; i < n; i++) {
			if (options.subset && !options.subset[i]) continue;
			
			DOTNode *node = dot.add_node(std::to_string(i));
			
			if (options.labels) {
				std::string label = std::to_string(i);
				node->add_attribute("label", label);
			} else {
				node->add_attribute("label", "");
			}
			
			if (options.shape) node->add_attribute("shape", options.shape);
			
			if (options.node_colors) {
				node->add_attribute("style", "filled");
				node->add_attribute("fillcolor", options.node_colors[i]);
			}
			
			if (options.node_edge_colors) {
				node->add_attribute("color", options.node_edge_colors[i]);
			}
			
			if (options.font_size) node->add_attribute("fontsize", options.font_size);
			if (options.font_colors) node->add_attribute("fontcolor", options.font_colors[i]);
		}
		
		// edges
		for (unsigned i = 0; i < n; i++) {
			for (unsigned j = 0; j < n;j++) {
				if (!has(i, j)) continue;
				if (options.subset && (!options.subset[i] || !options.subset[j])) continue;
				
				DOTEdge *edge = dot.add_edge(std::to_string(j), std::to_string(i));
				edge->add_attribute("dir", options.edge_dir);
			}
		}
		
		dot.print_to_file(f);
	}
	
	// finds a connected component in a subset of k elements, and writes the elements into the list component
	// returns the number of elements in the component
	unsigned find_connected_component(int *component, int *elements, unsigned k) const
	{
		bool added[k];
		for (unsigned i = 0; i < k; ++i) added[i] = false;
		
		int *head = component;
		int *tail = component + 1;
		*head = elements[0];
		added[0] = true;
		
		while (head < tail) {
			int u = *head;
			++head;
			
			for (unsigned i = 0; i < k; ++i) {
				if (added[i]) continue;
				int v = elements[i];
				if (has(u, v) || has(v, u)) {
					*tail = v;
					++tail;
					added[i] = true;
				}
			}
		}
		
		return tail - component;
	}
};



#endif
