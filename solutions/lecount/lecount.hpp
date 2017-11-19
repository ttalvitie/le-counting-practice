#ifndef LECOUNT_H
#define LECOUNT_H

#include "recursive.hpp"
#include "veie.hpp"
#include "approximate.hpp"

#define ALGO_RECURSIVE 0
#define ALGO_BUCKET_ELIMINATION 1
#define ALGO_APPROXIMATE 2

#define NUMTYPE_AUTO 0
#define NUMTYPE_AUTO_FLOATING_POINT 1
#define NUMTYPE_AUTO_INTEGER 2
#define NUMTYPE_DOUBLE 3
#define NUMTYPE_LONG_DOUBLE 4
#define NUMTYPE_UINT64 5
#define NUMTYPE_BIGINT 6


struct LECountOptions
{
	const char *dag_file;
	int algorithm;
	int numtype;
	RecursionOptions rec;
	VEIEOptions veie;
	ApproximationOptions approximation;
	int sample_count;
	const char *visprefix;
	unsigned int rng_seed;
	bool verbose;
	
	LECountOptions()
	{
		dag_file = NULL;
		algorithm = ALGO_RECURSIVE;
		numtype = NUMTYPE_AUTO;
		sample_count = 0;
		visprefix = NULL;
		rng_seed = 0;
		verbose = true;
	}
};

extern LECountOptions options;

#endif
