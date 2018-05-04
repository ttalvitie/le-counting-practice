#include "sampling.hpp"

thread_local mt19937 rng;

template <typename NextBeta>
double tpa(NextBeta nextBeta, double startBeta = 1.0, double endBeta = 0.0, double delta = 0.25, double epsilon = 1.0) {
	if(startBeta < endBeta) {
		fail("TPA requires that startBeta is at least endBeta");
	}
	auto run = [&](long long r) {
		long long k = 0;
		for(long long i = 0; i < r; ++i) {
			--k;
			double beta = startBeta;
			while(beta > endBeta) {
				double newBeta = nextBeta(beta);
				if(newBeta > beta) {
					fail("Beta increases in TPA; this should not happen");
				}
				beta = newBeta;
				++k;
			}
		}
		return (double)k / (double)r;
	};
	
	double A1 = run((long long)ceil(2.0 * log(2.0 / delta)));
	double s = log(1 + epsilon);
	return run((long long)ceil(2.0 * (A1 + sqrt(A1) + 2.0) * log(4.0 / delta) / (s * s * (1.0 - s))));
}

template <int W>
struct NextBeta {
	int n;
	PartialOrder<W> po;
	int home[PartialOrder<W>::MaxVertCount];
	int invHome[PartialOrder<W>::MaxVertCount];
	
	int conf[PartialOrder<W>::MaxVertCount];
	
	NextBeta(const PartialOrder<W>& po)
		: n(po.size()),
		  po(po)
	{
		int i = 0;
		po.topo([&](int v) {
			home[i] = v;
			invHome[v] = i;
			++i;
		});
	}
	
	void tryFillMissing() {
		if(conf[n - 1] == n) {
			bool present[PartialOrder<W>::MaxVertCount + 1];
			fill(present, present + n + 1, false);
			for(int i = 0; i < n; ++i) {
				present[conf[i]] = true;
			}
			for(int i = 0; i < n; ++i) {
				if(!present[home[i]]) {
					conf[n - 1] = home[i];
					break;
				}
			}
		}
	}
	
	void step(int i, bool H, bool C, int ceilBeta) {
		if(H && (conf[i] == n || conf[i + 1] == n || !po(conf[i], conf[i + 1]))) {
			if(conf[i] == n || ((C && i + 1 - invHome[conf[i]] <= ceilBeta) || i + 1 - invHome[conf[i]] < ceilBeta)) {
				swap(conf[i], conf[i + 1]);
			}
		}
		tryFillMissing();
	}
	
	void sample(long long steps, int ceilBeta, double weightFactor) {
		fill(conf, conf + n, n);
		tryFillMissing();
		
		mt19937 rngSaved = rng;
		
		uniform_int_distribution<int> posDist(0, n - 2);
		uniform_real_distribution<double> unitDist(0.0, 1.0);
		
		for(long long stepIdx = 0; stepIdx < steps; ++stepIdx) {
			int i = posDist(rng);
			bool H = (bool)(rng() & 1);
			bool C = unitDist(rng) < weightFactor;
			step(i, H, C, ceilBeta);
		}
		
		bool done = true;
		for(int i = 0; i < n; ++i) {
			if(conf[i] == n) {
				done = false;
				break;
			}
		}
		
		if(done) {
			return;
		}
		
		sample(2 * steps, ceilBeta, weightFactor);
		
		for(long long stepIdx = 0; stepIdx < steps; ++stepIdx) {
			int i = posDist(rngSaved);
			bool H = (bool)(rngSaved() & 1);
			bool C = unitDist(rngSaved) < weightFactor;
			step(i, H, C, ceilBeta);
		}
	}
	
	double operator()(double beta) {
		int ceilBeta = (int)ceil(beta);
		double weightFactor = 1.0 + beta - (double)ceilBeta;
		
		sample(1, ceilBeta, weightFactor);
		
		int maxDiff = 0;
		int maxDiffCount = 0;
		for(int i = 0; i < n; ++i) {
			int diff = i - invHome[conf[i]];
			if(diff > maxDiff) {
				maxDiff = diff;
				maxDiffCount = 0;
			}
			if(diff == maxDiff) {
				++maxDiffCount;
			}
		}
		assert(maxDiff <= ceilBeta);
		
		if(maxDiff == 0) {
			return 0.0;
		}
		
		double value;
		if(maxDiff == ceilBeta) {
			value = uniform_real_distribution<double>(0.0, pow(weightFactor, (double)maxDiffCount))(rng);
		} else {
			value = uniform_real_distribution<double>(0.0, 1.0)(rng);
		}
		
		int newCeilBeta = maxDiff;
		double newWeightFactor = pow(value, 1.0 / (double)maxDiffCount);
		double newBeta = (double)newCeilBeta + newWeightFactor - 1.0;
		
		return newBeta;
	}
};

struct Handler {
template <int W>
void operator()(PartialOrder<W> po) {
	cout << exp(tpa(NextBeta<W>(po), (double)po.size(), 0.0)) << '\n';
}
};

int main(int argc, char* argv[]) {
	if(argc != 2) {
		fail("Usage: <adjacency matrix filename>");
	}
	
	ifstream fp;
	fp.exceptions(fp.failbit | fp.badbit | fp.eofbit);
	fp.open(argv[1]);
	readPartialOrder(fp, Handler());
	fp.close();
	
	return 0;
}
