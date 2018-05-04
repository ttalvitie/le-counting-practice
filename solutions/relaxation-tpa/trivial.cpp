#include "partial_order.hpp"

thread_local mt19937 rng;

template <int W>
class TPA {
public:
	static double run(
		const PartialOrder<W>& po,
		double delta,
		double epsilon
	) {
		TPA<W> ctx(po);
		return lgamma(po.size() + 1) - ctx.tpa(delta, epsilon);
	}
	
private:
	static const int N = PartialOrder<W>::MaxVertCount;
	static const int MaxValue = (int)((uint32_t)-1 >> 1);
	static int randValue(mt19937& rng) {
		return (int)((uint32_t)rng() >> 1);
	}
	
	int n;
	int beta;
	
	int pool[2 * N * N];
	int* betaPredStart[N + 1];
	int* betaSuccStart[N + 1];
	
	int left[N];
	int right[N];
	
	TPA(const PartialOrder<W>& po) {
		n = po.size();
		
		BitSet<W> betaPred[N];
		BitSet<W> betaSucc[N];
		
		for(int v = 0; v < n; ++v) {
			BitSet<W> seen;
			po.topo(po.succ(v), [&](int x) {
				if(!seen[x]) {
					betaPred[x].add(v);
					betaSucc[v].add(x);
					seen |= po.succ(x);
				}
			});
		}
		
		int* pos = pool;
		betaPredStart[0] = pos;
		for(int v = 0; v < n; ++v) {
			betaPred[v].iterate([&](int x) {
				*pos++ = x;
			});
			betaPredStart[v + 1] = pos;
		}
		betaSuccStart[0] = pos;
		for(int v = 0; v < n; ++v) {
			betaSucc[v].iterate([&](int x) {
				*pos++ = x;
			});
			betaSuccStart[v + 1] = pos;
		}
	}
	
	void step(int v, int p, int* state) {
		if(p < state[v]) {
			for(int* x = betaPredStart[v]; x != betaPredStart[v + 1]; ++x) {
				if(state[*x] - p > beta) return;
			}
		} else {
			for(int* x = betaSuccStart[v]; x != betaSuccStart[v + 1]; ++x) {
				if(p - state[*x] > beta) return;
			}
		}
		state[v] = p;
	}
	
	void sample(long long steps = 1) {
		fill(left, left + n, 0);
		fill(right, right + n, MaxValue);
		
		uniform_int_distribution<int> vertDist(0, n - 1);
		
		mt19937 rngSaved = rng;
		for(long long stepIdx = 0; stepIdx < steps; ++stepIdx) {
			int v = vertDist(rng);
			int p = randValue(rng);
			
			step(v, p, left);
			step(v, p, right);
		}
		
		if(equal(left, left + n, right)) return;
		
		sample(2 * steps);
		
		for(long long stepIdx = 0; stepIdx < steps; ++stepIdx) {
			int v = vertDist(rngSaved);
			int p = randValue(rngSaved);
			
			step(v, p, left);
		}
	}
	
	void nextBeta() {
		sample();
		beta = 0;
		for(int v = 0; v < n; ++v) {
			for(int* x = betaSuccStart[v]; x != betaSuccStart[v + 1]; ++x) {
				beta = max(beta, left[v] - left[*x]);
			}
		}
	}
	
	double tpa(double delta, double epsilon) {
		auto run = [&](long long r) {
			long long k = 0;
			for(long long i = 0; i < r; ++i) {
				--k;
				beta = MaxValue;
				while(beta > 0) {
					nextBeta();
					++k;
				}
			}
			return (double)k / (double)r;
		};
		
		double A1 = run((long long)ceil(2.0 * log(2.0 / delta)));
		double s = log(1 + epsilon);
		return run((long long)ceil(2.0 * (A1 + sqrt(A1) + 2.0) * log(4.0 / delta) / (s * s * (1.0 - s))));
	}
};

struct Handler {
	template <int W>
	void operator()(PartialOrder<W> po) {
		double logResult = TPA<W>::run(po, 0.25, 1.0);
		long double result = exp((long double)logResult);
		
		cout << result << '\n';
	}
};

int main(int argc, char* argv[]) {
	if(argc != 2) {
		fail("Usage: <adjacency matrix filename>");
	}
	
	rng.seed(0);
	
	Handler handler;
	
	ifstream fp;
	fp.exceptions(fp.failbit | fp.badbit | fp.eofbit);
	fp.open(argv[1]);
	readPartialOrder(fp, handler);
	fp.close();
	
	return 0;
}
