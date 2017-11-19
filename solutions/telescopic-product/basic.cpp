#include "sampling.hpp"

thread_local mt19937 rng;

const int preliminarySamples = 300;

template <int W>
class SortMC {
public:
	SortMC(PartialOrder<W> po, LinextSamplingAlgorithm algo)
		: po(po),
		  count(1.0),
		  algo(algo)
	{
		vector<int> verts(po.size());
		for(int i = 0; i < po.size(); ++i) {
			verts[i] = i;
		}
		mergesort(verts);
		
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
			double p = estimateOrderProb(factor.po, factor.a, factor.b, samples, algo);
			count /= p;
		}
	}
	
	double getCount() {
		return count;
	}
	
private:
	struct Factor {
		PartialOrder<W> po;
		int a;
		int b;
	};
	
	void mergesort(vector<int>& verts) {
		if(verts.size() <= 1) {
			return;
		}
		
		vector<int> A(verts.begin(), verts.begin() + verts.size() / 2);
		vector<int> B(verts.begin() + verts.size() / 2, verts.end());
		mergesort(A);
		mergesort(B);
		
		int a = 0;
		int b = 0;
		
		verts.clear();
		while(a != (int)A.size() && b != (int)B.size()) {
			if(order(A[a], B[b])) {
				verts.push_back(A[a++]);
			} else {
				verts.push_back(B[b++]);
			}
		}
		while(a != (int)A.size()) {
			verts.push_back(A[a++]);
		}
		while(b != (int)B.size()) {
			verts.push_back(B[b++]);
		}
	}
	
	bool order(int a, int b) {
		if(po(a, b)) {
			return true;
		}
		if(po(b, a)) {
			return false;
		}
		double p = estimateOrderProb(po, a, b, preliminarySamples, algo);
		bool res = true;
		if(p < 0.5) {
			swap(a, b);
			res = false;
		}
		factors.push_back(Factor{po, a, b});
		po.add(a, b);
		return res;
	}
	
	PartialOrder<W> po;
	double count;
	LinextSamplingAlgorithm algo;
	vector<Factor> factors;
};

struct Handler {
	LinextSamplingAlgorithm algo;
	
	template <int W>
	void operator()(PartialOrder<W> po) {
		SortMC<W> smc(move(po), algo);
		cout << smc.getCount() << '\n';
	}
};

int main(int argc, char* argv[]) {
	Handler handler;
	
	if(argc != 3) {
		fail("Usage: <adjacency matrix filename> <sampling algorithm>");
	}
	
	handler.algo = parseLinextSamplingAlgorithmName(argv[2]);
	
	ifstream fp;
	fp.exceptions(fp.failbit | fp.badbit | fp.eofbit);
	fp.open(argv[1]);
	readPartialOrder(fp, handler);
	fp.close();
	
	return 0;
}
