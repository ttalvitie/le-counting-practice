#include <cstdlib>
#include <ctime>


double rnd()
{
	return rand() / ((double)RAND_MAX + 1);
}

int max(int a, int b)
{
	return a >= b ? a : b;
}

double seconds()
{
	return clock() * 1e-6;
}

