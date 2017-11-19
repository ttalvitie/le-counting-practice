#include "sampling.hpp"

thread_local mt19937 rng;

const int preliminarySamples = 300;

template <int W>
class SortMC {
public:
	SortMC(PartialOrder<W> po, int exactLimit, LinextSamplingAlgorithm algo)
		: po(po),
		  count(1.0),
		  algo(algo),
		  exactLimit(exactLimit)
	{
		sort(po.vertMask(), -1);
		
		double k = factors.size();
		
		auto failureProb2 = [&](double samples, double minProb) {
			double probDiff = 0.5 - minProb;
			return 4.0 * (pow(1.0 + 1.0 / (minProb * samples), k) - 1.0) + 2.0 * k * exp(-2.0 * probDiff * probDiff * (double)preliminarySamples);
		};
		
		auto failureProb = [&](double samples) {
			double A = 0.0;
			double B = 0.5;
			double V = 1.0 / 0.0;
			while(B - A > 1e-4) {
				double M1 = (2.0 * A + B) / 3.0;
				double M2 = (A + 2.0 * B) / 3.0;
				double V1 = failureProb2(samples, M1);
				double V2 = failureProb2(samples, M2);
				V = min(V, V1);
				V = min(V, V2);
				if(V1 < V2) {
					B = M2;
				} else {
					A = M1;
				}
			}
			return V;
		};
		
		double A = 1.0;
		double B = 1.0;
		while(failureProb(B) > 0.25) {
			B *= 2;
		}
		
		while(B - A >= 1) {
			double M = 0.5 * (A + B);
			if(failureProb(M) > 0.25) {
				A = M;
			} else {
				B = M;
			}
		}
		
		int samples = (int)ceil(B);
		
		cerr << factors.size() << " factors, " << samples << " samples per factor.\n";
		
		for(const Factor& factor : factors) {
			double p = estimateOrderProb(factor.po, factor.verts, factor.a, factor.b, samples, algo);
			count /= p;
		}
	}
	
	double getCount() {
		return count;
	}
	
private:
	struct Factor {
		PartialOrder<W> po;
		BitSet<W> verts;
		int a;
		int b;
	};
	
	void sort(BitSet<W> verts, int pivot) {
		double logcoef = lgamma(verts.count() + 1);
		po.components(verts, [&](BitSet<W> comp) {
			sortComponent(comp, pivot);
			
			logcoef -= lgamma(comp.count() + 1);
		});
		count *= exp(logcoef);
	}
	
	bool trySplit(BitSet<W> verts, int pivot) {
		BitSet<W> A;
		BitSet<W> B = verts;
		int bestVal = -1;
		BitSet<W> bestA;
		BitSet<W> bestB;
		po.topo(verts, [&](int v) {
			A.add(v);
			B &= po.succ(v);
			if(B && verts == (A | B)) {
				int val = min(A.count(), B.count());
				if(val > bestVal) {
					bestVal = val;
					bestA = A;
					bestB = B;
				}
			}
		});
		if(bestVal != -1) {
			sort(bestA, pivot);
			sort(bestB, pivot);
			return true;
		}
		return false;
	}
	
	void sortComponent(BitSet<W> verts, int pivot) {
		if(verts.count() <= 1) {
			return;
		}
		if(trySplit(verts, pivot)) {
			return;
		}
		
		if(exactLimit) {
			long long exactExpansionsLeft = exactLimit;
			double res = po.tryLinextCount(verts, exactExpansionsLeft);
			if(exactExpansionsLeft) {
				count *= res;
				return;
			}
		}
		
		if(pivot == -1 || !verts[pivot]) {
			int bestVal = -1;
			verts.iterate([&](int v) {
				int pred = (po.pred(v) & verts).count();
				int succ = (po.succ(v) & verts).count();
				int val = min(pred, succ);
				if(val > bestVal) {
					bestVal = val;
					pivot = v;
				}
			});
		}
		
		BitSet<W> cand = verts;
		cand.del(pivot);
		cand &= ~po.pred(pivot);
		cand &= ~po.succ(pivot);
		int v = cand.random();
		order(verts, pivot, v);
		sortComponent(verts, pivot);
	}
	
	void order(BitSet<W> verts, int a, int b) {
		if(po(a, b) || po(b, a)) {
			return;
		}
		double p = estimateOrderProb(po, verts, a, b, preliminarySamples, algo);
		if(p < 0.5) {
			swap(a, b);
		}
		factors.push_back(Factor{po, verts, a, b});
		po.add(a, b);
	}
	
	PartialOrder<W> po;
	double count;
	LinextSamplingAlgorithm algo;
	vector<Factor> factors;
	long long exactLimit;
};

struct Handler {
	template <int W>
	void operator()(PartialOrder<W> po) {
		SortMC<W> smc(move(po), exactLimit, algo);
		cout << smc.getCount() << '\n';
	}
	
	long long exactLimit;
	LinextSamplingAlgorithm algo;
};

int main(int argc, char* argv[]) {
	Handler handler;
	handler.exactLimit = 0;
	
	if(argc < 3 || argc > 4) {
		fail("Usage: <adjacency matrix filename> <sampling algorithm> [<exact counter limit>]");
	}
	
	handler.algo = parseLinextSamplingAlgorithmName(argv[2]);
	if(argc > 3) {
		handler.exactLimit = fromString<long long>(argv[3]);
	}
	
	ifstream fp;
	fp.exceptions(fp.failbit | fp.badbit | fp.eofbit);
	fp.open(argv[1]);
	readPartialOrder(fp, handler);
	fp.close();
	
	return 0;
}
