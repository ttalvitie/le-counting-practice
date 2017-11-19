#include <cstdio>
#include <cstdint>

#include <gmpxx.h>

#include "lecount.hpp"
#include "set.hpp"
#include "recursive.hpp"
#include "veie.hpp"
#include "approximate.hpp"

LECountOptions options;


template <typename N>
void run_recursive(digraph &poset)
{
	unsigned n = poset.n;
	
	if (options.sample_count) options.rec.sampling = true;
	options.rec.vis.prefix = options.visprefix;
	options.rec.verbose = options.verbose;
	RecursiveCounterAuto<N> *rec = get_recursive_counter<N>(n, options.rec);
	
	N count = rec->count_linex(poset);
	if (options.rec.verbose) rec->display_statistics();
	
	int sample[n];
	for (unsigned k = 0; k < options.sample_count; ++k) {
		rec->draw_sample(sample);
		for (unsigned i = 0; i < n; ++i) {
			printf("%i ", sample[i]);
		}
		printf("\n");
	}
	
	rec->deallocate();
	delete rec;
	
	std::cout << count << std::endl;
}

template <typename N>
void run_armc(digraph &poset)
{
	options.approximation.relaxation.visprefix = options.visprefix;
	options.approximation.verbose = options.verbose;
	MonteCarloEstimator<N> mce(options.approximation, poset);
// 	double dp_time = seconds();
	N count = mce.approximate();
// 	double sample_time = seconds() - dp_time;
// 	std::cout << count << " " << dp_time << " " << sample_time << std::endl;
	std::cout << count << std::endl;
}

template <typename N>
void run_veie(digraph &poset)
{
	options.veie.verbose = options.verbose;
	VEIE<N> veie(poset.n, options.veie);
	N count = veie.count_extensions(poset);
	std::cout << count << std::endl;
}


void parse_arguments(char **argv);

int main(int, char **argv)
{
	// read command line arguments and the input DAG
	parse_arguments(argv);
	if (options.verbose) eprintf("DAG: %s\n", options.dag_file);
	digraph poset(options.dag_file);
	
	srand(options.rng_seed > 0 ? options.rng_seed : time(NULL));
	
	// if auto, choose auto type for the respective algorithm
	if (options.numtype == NUMTYPE_AUTO) {
		if (options.algorithm == ALGO_RECURSIVE || options.algorithm == ALGO_BUCKET_ELIMINATION) {
			options.numtype = NUMTYPE_AUTO_INTEGER;
		} else {
			options.numtype = NUMTYPE_AUTO_FLOATING_POINT;
		}
	}
	
	if (options.numtype == NUMTYPE_AUTO_INTEGER) {
		if (poset.n <= 20) {
			options.numtype = NUMTYPE_UINT64;
		} else {
			options.numtype = NUMTYPE_BIGINT;
		}
	} else if (options.numtype == NUMTYPE_AUTO_FLOATING_POINT) {
		if (poset.n <= 170) {
			options.numtype = NUMTYPE_DOUBLE;
// 		} else if (poset.n <= 1754) {
		} else {
			options.numtype = NUMTYPE_LONG_DOUBLE;
		}
	}
	
	if (options.algorithm == ALGO_RECURSIVE) {
		if      (options.numtype == NUMTYPE_DOUBLE)      run_recursive<double>(poset);
		else if (options.numtype == NUMTYPE_LONG_DOUBLE) run_recursive<long double>(poset);
		else if (options.numtype == NUMTYPE_UINT64)      run_recursive<uint64_t>(poset);
		else if (options.numtype == NUMTYPE_BIGINT)      run_recursive<mpz_class>(poset);
		else {
			eprintf("Invalid numtype: %i\n", options.numtype);
			return 1;
		}
	} else if (options.algorithm == ALGO_APPROXIMATE) {
		if      (options.numtype == NUMTYPE_DOUBLE)      run_armc<double>(poset);
		else if (options.numtype == NUMTYPE_LONG_DOUBLE) run_armc<long double>(poset);
		else {
			eprintf("Invalid numtype: %i\n", options.numtype);
			return 1;
		}
	} else if (options.algorithm == ALGO_BUCKET_ELIMINATION) {
		if      (options.numtype == NUMTYPE_UINT64) run_veie<uint64_t>(poset);
		else if (options.numtype == NUMTYPE_BIGINT) run_veie<mpz_class>(poset);
		else {
			eprintf("Invalid numtype: %i\n", options.numtype);
			return 1;
		}
	}
	
	return 0;
}

