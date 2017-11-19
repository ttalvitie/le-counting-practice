#ifndef TOOLS_H
#define TOOLS_H

#define eprintf(...) { fprintf(stderr, __VA_ARGS__); }
#define eflush() { fflush(stderr); }

double rnd();
int max(int a, int b);
double seconds();


// precomputed factorials and binomial coefficients
template <typename N>
struct Combinatorial
{
	unsigned n;
	N *factorials;
	N **binomials;
	
	Combinatorial(unsigned n) : n(n), factorials(NULL), binomials(NULL) {}
	
	~Combinatorial()
	{
		if (factorials) {
			delete [] factorials;
		}
		
		if (binomials) {
			for (unsigned i = 0; i <= n; ++i) {
				delete [] binomials[i];
			}
			delete [] binomials;
		}
	}
	
	void compute_factorials()
	{
		factorials = new N[n+1];
		
		factorials[0] = 1;
		for (unsigned i = 1; i <= n; ++i) {
			factorials[i] = factorials[i-1] * i;
		}
	}
	
	void compute_binomials()
	{
		binomials = new N*[n+1];
		
		for (unsigned i = 0; i <= n; ++i) {
			binomials[i] = new N[i/2 + 1];
			binomials[i][0] = 1;
			for (unsigned j = 1; j <= i/2; ++j) {
				binomials[i][j] = binomial(i-1, j-1) + binomial(i-1, j);
			}
		}
	}
	
	N binomial(unsigned m, unsigned k) const
	{
		return binomials[m][(k <= m - k) ? k : (m - k)];
	}
	
	N factorial(unsigned m) const
	{
		return factorials[m];
	}
};

#endif
