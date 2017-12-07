#include <cstdio>
#include <cstring>

#include "lecount.hpp"


void print_usage(const char *cmd)
{
	printf("\nUsage: %s <input file> [OPTIONS ...]\n\n", cmd);
	printf("The main options are:\n");
	printf(" --algorithm=[dp|veie|armc]         the algorithm to use for counting\n");
	printf(" --numtype                          the datatype to use in counting\n");
	printf(" --visprefix                        the prefix (path) for visualization output\n");
	printf(" --seed                             RNG seed\n");
	printf(" --quiet                            don't output extra information\n");
	printf("\nEnter any option with no '=' or value to see a detailed description.\n");
	printf("Each algorithm accepts its own set of options:\n");
	printf("\nOptions for 'dp':\n");
	printf(" --transpose=[no|yes|auto|auto2]    whether to transpose the poset\n");
	printf(" --ccs=[none|dfs|covers|auto]       how to find connected components\n");
	printf(" --hub=[none|best|perfect]          whether to use hub splits\n");
	printf(" --sort=[yes|no]                    whether to sort the poset topologically\n");
	printf(" --visualize=[slLHpf]               visualization options\n");
	printf(" --sample                           number of linear extensions to sample after counting\n");
	printf(" --cache-params                     cache size parameters\n");
	printf("\nOptions for 'veie':\n");
	printf(" --elim-order-file                  file to read the elimination order from\n");
	printf(" --random-orders                    number of random orders to sample\n");
	printf("\nOptions for 'armc':\n");
	printf(" --partition-set-size               maximum size of sets in a partition relaxation\n");
	printf(" --partition-search-time            time to search for a partition relaxation\n");
	printf(" --iterative-sampling               whether to use iterative sampling\n");
	printf(" --relaxation-transpose             whether to transpose (parts of) the relaxation\n");
	printf(" --epsilon                          value of ε for an (ε,δ)-estimate\n");
	printf(" --delta                            value of δ for an (ε,δ)-estimate\n");
	printf("\n");
}

void parse_algorithm(const char *arg)
{
	if (!strcmp(arg, "dp")) {
		options.algorithm = ALGO_RECURSIVE;
	} else if (!strcmp(arg, "veie")) {
		options.algorithm = ALGO_BUCKET_ELIMINATION;
	} else if (!strcmp(arg, "armc")) {
		options.algorithm = ALGO_APPROXIMATE;
	} else {
		printf("The option 'algorithm' chooses the method for counting linear extensions.\n");
		printf("The alternatives are:\n");
		printf(" --algorithm=dp       dynamic programming (DP) over the upset lattice (default)\n");
		printf(" --algorithm=veie     variable elimination via inclusion-exclusion\n");
		printf(" --algorithm=armc     ARMC approximation scheme, uses DP as a subroutine\n");
		exit(1);
	}
}

void parse_hub(const char *arg)
{
	if (!strcmp(arg, "none")) {
		options.rec.hub_splits = HUB_NONE;
	} else if (!strcmp(arg, "best")) {
		options.rec.hub_splits = HUB_BEST;
	} else if (!strcmp(arg, "perfect")) {
		options.rec.hub_splits = HUB_PERFECT;
	} else {
		printf("The option 'hub' enables the 'admissible partitions' technique for counting\n");
		printf("linear extensions. Each step chooses a single hub element and splits the poset\n");
		printf("into its successors and predecessors in all possible ways.\n");
		printf("The options are:\n");
		printf(" --hub=none      no hub splits (default)\n");
		printf(" --hub=best      always choose the best hub, with most comparable elements\n");
		printf("                 (usable up to 256 elements)\n");
		printf(" --hub=perfect   do perfect splits only (all elements comparable with the hub)\n");
		exit(1);
	}
}

void parse_ccs(const char *arg)
{
	if (!strcmp(arg, "none")) {
		options.rec.ccs = CCS_NONE;
	} else if (!strcmp(arg, "dfs")) {
		options.rec.ccs = CCS_DFS;
	} else if (!strcmp(arg, "covers")) {
		options.rec.ccs = CCS_COVERS;
	} else if (!strcmp(arg, "auto")) {
		options.rec.ccs = CCS_AUTO;
	} else {
		printf("The option 'ccs' determines the method of detecting connected components.\n");
		printf("The options are:\n");
		printf(" --ccs=none      never split into connected components\n");
		printf(" --ccs=dfs       find components by depth-first search\n");
		printf(" --ccs=covers    find components by maximal cover intersections\n");
		printf(" --ccs=auto      use an average-degree heuristic to decide between dfs and covers (default)\n");
		exit(1);
	}
}

void parse_transpose(const char *arg)
{
	if (!strcmp(arg, "no")) {
		options.rec.transpose = TRANSPOSE_NO;
	} else if (!strcmp(arg, "yes")) {
		options.rec.transpose = TRANSPOSE_YES;
	} else if (!strcmp(arg, "auto")) {
		options.rec.transpose = TRANSPOSE_AUTO;
	} else if (!strcmp(arg, "auto2")) {
		options.rec.transpose = TRANSPOSE_AUTO2;
	} else {
		printf("The option 'transpose' determines whether to transpose the poset before counting.\n");
		printf("The options are:\n");
		printf(" --transpose=no     do not transpose\n");
		printf(" --transpose=yes    do tranpose\n");
		printf(" --transpose=auto   decide based on a simple heuristic\n");
		printf(" --transpose=auto2  decide based on a more advanced heuristic (default)\n");
		exit(1);
	}
}

void parse_numtype(const char *arg)
{
	if (!strcmp(arg, "auto")) {
		options.numtype = NUMTYPE_AUTO;
	} else if (!strcmp(arg, "auto-fp")) {
		options.numtype = NUMTYPE_AUTO_FLOATING_POINT;
	} else if (!strcmp(arg, "auto-int")) {
		options.numtype = NUMTYPE_AUTO_INTEGER;
	} else if (!strcmp(arg, "double")) {
		options.numtype = NUMTYPE_DOUBLE;
	} else if (!strcmp(arg, "long-double")) {
		options.numtype = NUMTYPE_LONG_DOUBLE;
	} else if (!strcmp(arg, "int64")) {
		options.numtype = NUMTYPE_UINT64;
	} else if (!strcmp(arg, "bigint")) {
		options.numtype = NUMTYPE_BIGINT;
	} else {
		printf("The option 'numtype' chooses the datatype for representing the number of linear extensions.\n");
		printf("By default a datatype is chosen that can contain the factorial of the number of elements.\n");
		printf("The ARMC scheme only works with floating points while VEIE only works with integers.\n");
		printf("Take care when choosing the type manually to ensure it is large/accurate enough.\n");
		printf("The alternatives are:\n");
		printf(" --numtype=auto         automatic choice based on algorithm and poset size (default)\n");
		printf(" --numtype=auto-fp      automatic choice of floating-point type based on poset size\n");
		printf(" --numtype=auto-int     automatic choice of integer type based on poset size\n");
		printf(" --numtype=double       double\n");
		printf(" --numtype=long-double  long double\n");
		printf(" --numtype=int64        64-bit integer\n");
		printf(" --numtype=bigint       GMP integer (unlimited size)\n");
		exit(1);
	}
}

void parse_sort(const char *arg)
{
	if (!strcmp(arg, "no")) {
		options.rec.sort_topologically = SORT_NO;
	} else if (!strcmp(arg, "yes")) {
		options.rec.sort_topologically = SORT_YES;
	} else {
		printf("The options for 'sort' are:\n");
		printf(" --sort=no    do not sort topologically\n");
		printf(" --sort=yes   sort topologically (default)\n");
		exit(1);
	}
}

void parse_elimination_order(const char *arg)
{
	if (*arg != '\0') {
		options.veie.elim_order_file = arg;
		return;
	}
	
	printf("The option 'elim-order-file' sets a file to read the elimination order from.\n");
	printf("The order is given as a list of elements separated by spaces. If no elimination\n");
	printf("order is specified in this way, it will be selected using a number of heuristics.\n");
	exit(1);
}

void parse_random_orders(const char *arg)
{
	if (sscanf(arg, "%i", &options.veie.n_random_orders) == 1) return;
	printf("The option 'random-orders' selects the number of random elimination orders\n");
	printf("to sample. The default is random-orders=1000.\n");
	exit(1);
}

void parse_visualize(const char *arg)
{
	if (*arg == '\0') {
		printf("The option 'visualize' sets flags that control visualization.\n");
		printf("It can only be used together with 'visprefix' option.\n");
		printf("The flags are:\n");
		printf("  s   create a .dot file for each subset of the poset evaluated\n");
		printf("  l   create a .dot file for the entire lattice of subsets evaluated\n");
		printf("  L   as above but use subset images as nodes of the lattice\n");
		printf("  H   draw the poset as a Hasse diagram (no edge directions)\n");
		printf("  p   draw each element as a point with no label\n");
		printf("  f   in the lattice, draw a frame around each subset\n");
		printf("\nExample usage: --visualize=sLpf\n");
		exit(1);
	}
	
	while (*arg != '\0') {
		int f = *arg;
		if (f == 's') {
			options.rec.vis.subsets = true;
		} else if (f == 'l') {
			options.rec.vis.lattice = LATTICE_NODES;
		} else if (f == 'L') {
			options.rec.vis.lattice = LATTICE_IMAGES;
			options.rec.vis.subsets = true;
		} else if (f == 'H') {
			options.rec.vis.edge_directions = false;
		} else if (f == 'p') {
			options.rec.vis.node_shape = NODE_SHAPE_POINT;
		} else if (f == 'f') {
			options.rec.vis.subset_box = true;
		} else {
			printf("Error: Unknown --visualize flag: %c\n", f);
			exit(1);
		}
		++arg;
	}
}

void parse_cache_params(const char *arg)
{
	unsigned n = sscanf(arg, "%d,%d,%lf,%d,%d",
		&options.rec.cache.initial_prime_index,
		&options.rec.cache.prime_index_increment,
		&options.rec.cache.max_load_factor,
		&options.rec.cache.initial_array_size,
		&options.rec.cache.array_resize_factor);
	
	if (n != 5) {
		printf("--cache-params requires 5 comma separated values:\n");
		printf(" - the initial cache size, given as the index of precomputed prime numbers\n");
		printf(" - the increment to the index on resize\n");
		printf(" - maximum load factor\n");
		printf(" - initial array size\n");
		printf(" - array resize factor\n");
		printf("\nThe default is --cache-params=4,3,0.5,1024,4\n");
		exit(1);
	}
	
	if (!(options.rec.cache.initial_prime_index >= 0 && options.rec.cache.initial_prime_index <= 25)) {
		printf("The option 'cache-init-prime-index' needs to be an integer 0..25.\n");
		exit(1);
	}
	
	if (!(options.rec.cache.prime_index_increment >= 1)) {
		printf("The option 'cache-prime-index-increment' needs to a positive integer.\n");
		exit(1);
	}
	
	if (!( options.rec.cache.max_load_factor > 0)) {
		printf("The option 'cache-max-load-factor' needs to be a positive number.\n");
		exit(1);
	}
	
	if (!(options.rec.cache.initial_array_size > 0)) {
		printf("The option 'cache-init-array-size' needs to be a positive integer.\n");
		exit(1);
	}
	
	if (!(options.rec.cache.array_resize_factor > 1)) {
		printf("The option 'cache-array-resize-factor' needs to be an integer larger than 1.\n");
		exit(1);
	}
	
	eprintf("Cache params: %d, %d, %lf, %d, %d\n",
		options.rec.cache.initial_prime_index,
		options.rec.cache.prime_index_increment,
		options.rec.cache.max_load_factor,
		options.rec.cache.initial_array_size,
		options.rec.cache.array_resize_factor);
}

void parse_rng_seed(const char *arg)
{
	if (sscanf(arg, "%i", &options.rng_seed) == 1) return;
	printf("The option 'seed' sets a positive seed for the random number generator.\n");
	printf("The default value --seed=0 always takes the seed from the clock.\n");
	exit(1);
}

void parse_sample_count(const char *arg)
{
	if (sscanf(arg, "%i", &options.sample_count) == 1) return;
	printf("The option 'sample' needs to be an integer.\n");
	exit(1);
}

void parse_partition_set_size(const char *arg)
{
	if (sscanf(arg, "%i", &options.approximation.relaxation.partition_max_size) == 1) {
		if (options.approximation.relaxation.partition_max_size >= 1) options.approximation.relaxation.relaxation = RELAXATION_MANUAL;
		return;
	}
	
	printf("The option 'partition-set-size' sets the maximum size of set in a partition relaxation.\n");
	printf("The default value, --partition-set-size=0, selects the size adaptively using the ARMC scheme.\n");
	exit(1);
}

void parse_partition_search_time(const char *arg)
{
	if (sscanf(arg, "%lf", &options.approximation.relaxation.partition_search_time) == 1) return;
	printf("The option 'partition-search-time' sets the number of seconds to spend searching for\n");
	printf("a good partition relaxation. The default value, partition-search-time=0, chooses the\n");
	printf("first relaxation found by greedy hillclimbing.\n");
	exit(1);
}

void parse_relaxation_transpose(const char *arg)
{
	if (!strcmp(arg, "no")) {
		options.approximation.relaxation.transpose = RELAXATION_TRANSPOSE_NO;
	} else if (!strcmp(arg, "full")) {
		options.approximation.relaxation.transpose = RELAXATION_TRANSPOSE_FULL;
	} else if (!strcmp(arg, "parts")) {
		options.approximation.relaxation.transpose = RELAXATION_TRANSPOSE_PARTS;
	} else {
		printf("Transposing a relaxation can speed up the DP phase significantly. In some cases \n");
		printf("it helps to transpose some connected components of the relaxation while leaving\n");
		printf("others unchanged; however, this option is incompatible with iterative sampling.\n");
		printf("The options are:\n");
		printf(" --relaxation-transpose=no     never transpose the relaxation\n");
		printf(" --relaxation-transpose=full   decide whether to transpose the entire relaxation (default)\n");
		printf(" --relaxation-transpose=parts  decide for each part independently whether to transpose them\n");
		exit(1);
	}
}

void parse_iterative_sampling(const char *arg)
{
	if (!strcmp(arg, "no")) {
		options.approximation.relaxation.use_iterative_sampling = false;
	} else if (!strcmp(arg, "yes")) {
		options.approximation.relaxation.use_iterative_sampling = true;
	} else {
		printf("Iterative sampling draws the elements of a linear order in lazy manner, one by one,\n");
		printf("which enables early rejection of orders that are not linear extensions.\n");
		printf("The options are:\n");
		printf(" --iterative-sampling=no    never use iterative sampling\n");
		printf(" --iterative-sampling=yes   use iterative sampling when possible (default)\n");
		exit(1);
	}
}

void parse_approximation_epsilon(const char *arg)
{
	if (sscanf(arg, "%lf", &options.approximation.epsilon) == 1 && options.approximation.epsilon > 0) return;
	printf("The option 'epsilon' needs to be a positive number.\n");
	exit(1);
}

void parse_approximation_delta(const char *arg)
{
	if (sscanf(arg, "%lf", &options.approximation.delta) == 1 && options.approximation.delta > 0) return;
	printf("The option 'delta' needs to be a positive number.\n");
	exit(1);
}



void parse_flag_argument(const char *flag, const char *arg)
{
	if (!strcmp(flag, "algorithm")) {
		parse_algorithm(arg);
	} else if (!strcmp(flag, "numtype")) {
		parse_numtype(arg);
	} else if (!strcmp(flag, "visprefix")) {
		options.visprefix = arg;
	} else if (!strcmp(flag, "seed")) {
		parse_rng_seed(arg);
	} else if (!strcmp(flag, "quiet")) {
		options.verbose = false;
	
	} else if (!strcmp(flag, "transpose")) {
		parse_transpose(arg);
	} else if (!strcmp(flag, "ccs")) {
		parse_ccs(arg);
	} else if (!strcmp(flag, "hub")) {
		parse_hub(arg);
	} else if (!strcmp(flag, "sort")) {
		parse_sort(arg);
	} else if (!strcmp(flag, "cache-params")) {
		parse_cache_params(arg);
	} else if (!strcmp(flag, "visualize")) {
		parse_visualize(arg);
	} else if (!strcmp(flag, "sample")) {
		parse_sample_count(arg);
	
	} else if (!strcmp(flag, "elim-order-file")) {
		parse_elimination_order(arg);
	} else if (!strcmp(flag, "random-orders")) {
		parse_random_orders(arg);
	
	} else if (!strcmp(flag, "partition-set-size")) {
		parse_partition_set_size(arg);
	} else if (!strcmp(flag, "partition-search-time")) {
		parse_partition_search_time(arg);
	} else if (!strcmp(flag, "relaxation-transpose")) {
		parse_relaxation_transpose(arg);
	} else if (!strcmp(flag, "iterative-sampling")) {
		parse_iterative_sampling(arg);
	} else if (!strcmp(flag, "epsilon")) {
		parse_approximation_epsilon(arg);
	} else if (!strcmp(flag, "delta")) {
		parse_approximation_delta(arg);
	
	} else {
		printf("Error: Unknown flag: %s\n", flag);
		exit(1);
	}
}

void parse_flag(const char *flag)
{
	const char *equals = strchr(flag, '=');
	
	if (equals == NULL) {
		parse_flag_argument(flag, "");
	} else {
		unsigned n = equals - flag;
		char name[n+1];
		strncpy(name, flag, n);
		name[n] = '\0';
		parse_flag_argument(name, equals+1);
	}
}

void parse_arguments(char **argv)
{
	const char *cmd = *argv++;
	
	while (char *arg = *argv) {
		if (arg[0] == '-' && arg[1] == '-') {
			parse_flag(arg+2);
		} else {
			if (options.dag_file != NULL) {
				printf("Error: Multiple input files were given.");
				exit(1);
			}
			options.dag_file = arg;
		}
		++argv;
	}
	
	if (options.dag_file == NULL) {
		print_usage(cmd);
		exit(1);
	}
	
	if (options.rec.vis.subsets || options.rec.vis.lattice) {
		if (options.visprefix == NULL) {
			printf("The --visprefix option must be set for visualization.\n");
			exit(1);
		}
	}
}

