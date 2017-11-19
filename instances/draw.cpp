#include <cassert>
#include "../solutions/lecount/digraph.hpp"

int main(int argc, char **argv)
{
	assert(argc == 3);
	
	digraph dg(argv[1]);
	dg.remove_transitive_arcs();
	
	DOTOptions options;
// 	options.shape = "circle";
	FILE *f = fopen(argv[2], "w");
	assert(f);
	dg.print_dot(f, options);
	fclose(f);
}
