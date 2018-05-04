#pragma once

#include "partial_order.hpp"

template <int W>
class SwapLinextSampler {
public:
	SwapLinextSampler(const PartialOrder<W>& po, BitSet<W> verts) : po(po) {
		verts &= po.vertMask();
		
		int order[PartialOrder<W>::MaxVertCount];
		int descendantCount[PartialOrder<W>::MaxVertCount];
		
		k = 0;
		verts.iterate([&](int v) {
			order[k] = v;
			descendantCount[v] = (po.succ(v) & verts).count() + 1;
			++k;
		});
		
		sort(order, order + k, [&](int a, int b) {
			return descendantCount[a] > descendantCount[b];
		});
		
		int p = 0;
		for(int i = 0; i < k - 1; ++i) {
			if(descendantCount[order[p]] >= k - i) {
				initialState[i] = order[p++];
			} else {
				initialState[i] = -1;
			}
		}
		initialState[k - 1] = order[p++];
		
		queueSize = 0;
		for(; p < k; ++p) {
			queue[queueSize++] = order[p];
		}
	}
	
	template <typename F>
	void sample(F f) {
		sample_(1);
		f((const int*)state);
	}
	
private:
	void sample_(int steps) {
		mt19937 rngSaved = rng;
		
		uniform_int_distribution<int> posDist(0, k - 2);
		
		copy(initialState, initialState + k, state);
		int queuePos = 0;
		
		for(int step = 0; step < steps; ++step) {
			bool d = (bool)(rng() & 1);
			int j = posDist(rng);
			
			if(d) {
				if(state[j] == -1 || state[j + 1] == -1 || !po(state[j], state[j + 1])) {
					swap(state[j], state[j + 1]);
					if(state[k - 1] == -1) {
						state[k - 1] = queue[queuePos++];
					}
				}
			}
		}
		
		if(queuePos == queueSize) {
			return;
		}
		
		sample_(2 * steps);
		
		for(int step = 0; step < steps; ++step) {
			bool d = (bool)(rngSaved() & 1);
			int j = posDist(rngSaved);
			
			if(d) {
				if(!po(state[j], state[j + 1])) {
					swap(state[j], state[j + 1]);
				}
			}
		}
	}
	
	const PartialOrder<W>& po;
	int k;
	int state[PartialOrder<W>::MaxVertCount];
	int initialState[PartialOrder<W>::MaxVertCount];
	int queue[PartialOrder<W>::MaxVertCount];
	int queueSize;
};

template <int W, typename F>
void sampleLinextSwap(const PartialOrder<W>& po, BitSet<W> verts, int sampleCount, F f) {
	SwapLinextSampler<W> sampler(po, verts);
	for(int i = 0; i < sampleCount; ++i) {
		sampler.sample(f);
	}
}

template <int W>
class GibbsLinextSampler {
public:
	GibbsLinextSampler(const PartialOrder<W>& po, BitSet<W> verts) {
		verts &= po.vertMask();
		
		n = po.size();
		k = 0;
		verts.iterate([&](int v) {
			order[k] = v;
			subVerts[k] = v;
			++k;
		});
		
		verts.iterate([&](int v) {
			BitSet<W> seen;
			po.topo(verts & po.succ(v), [&](int x) {
				if(!seen[x]) {
					pred[x].add(v);
					succ[v].add(x);
					seen |= po.succ(x);
				}
			});
		});
	}
	
	template <typename F>
	void sample(F f) {
		sample_(1);
		sort(order, order + k, [&](int a, int b) {
			return left[a] < left[b];
		});
		f((const int*)order);
	}
	
private:
	int n;
	int k;
	int order[PartialOrder<W>::MaxVertCount];
	int subVerts[PartialOrder<W>::MaxVertCount];
	BitSet<W> pred[PartialOrder<W>::MaxVertCount];
	BitSet<W> succ[PartialOrder<W>::MaxVertCount];
	uint32_t left[PartialOrder<W>::MaxVertCount];
	uint32_t right[PartialOrder<W>::MaxVertCount];
	
	void redistributeLeft_() {
		sort(order, order + k, [&](int a, int b) {
			return left[a] < left[b];
		});
		for(int i = 0; i < k; ++i) {
			right[i] = rng();
		}
		sort(right, right + k);
		for(int i = 0; i < k; ++i) {
			left[order[i]] = right[i];
		}
	}
	
	void sample_(int steps) {
		mt19937 rngSaved = rng;
		
		fill(left, left + n, 0);
		fill(right, right + n, (uint32_t)-1);
		
		auto runStep = [&](uint32_t* state, int j, uint32_t u) {
			uint32_t a = 0;
			uint32_t b = (uint32_t)-1;
			pred[j].iterate([&](int x) {
				if(state[x] > a) {
					a = state[x];
				}
			});
			succ[j].iterate([&](int x) {
				if(state[x] < b) {
					b = state[x];
				}
			});
			state[j] = a + (uint32_t)(((uint64_t)u * (uint64_t)(b - a)) >> 32);
		};
		
		for(int step = 0; step < steps; ++step) {
			int j = subVerts[uniform_int_distribution<int>(0, k - 1)(rng)];
			uint32_t u = rng();
			
			runStep(left, j, u);
			runStep(right, j, u);
		}
		
		sort(order, order + k, [&](int a, int b) {
			return left[a] < left[b];
		});
		
		bool ok = true;
		for(int i = 1; i < k; ++i) {
			if(right[order[i - 1]] >= left[order[i]]) {
				ok = false;
				break;
			}
		}
		if(ok) {
			redistributeLeft_();
			return;
		}
		
		sample_(2 * steps);
		
		for(int step = 0; step < steps; ++step) {
			int j = subVerts[uniform_int_distribution<int>(0, k - 1)(rngSaved)];
			uint32_t u = rngSaved();
			
			runStep(left, j, u);
		}
		redistributeLeft_();
	}
};

template <int W, typename F>
void sampleLinextGibbs(const PartialOrder<W>& po, BitSet<W> verts, int sampleCount, F f) {
	GibbsLinextSampler<W> sampler(po, verts);
	for(int i = 0; i < sampleCount; ++i) {
		sampler.sample(f);
	}
}

enum class LinextSamplingAlgorithm {
	SwapLinextSampler,
	GibbsLinextSampler
};
static const LinextSamplingAlgorithm DefaultLinextSamplingAlgorithm = LinextSamplingAlgorithm::GibbsLinextSampler;

inline LinextSamplingAlgorithm parseLinextSamplingAlgorithmName(const string& name) {
	if(name == "SwapLinextSampler") {
		return LinextSamplingAlgorithm::SwapLinextSampler;
	}
	if(name == "GibbsLinextSampler") {
		return LinextSamplingAlgorithm::GibbsLinextSampler;
	}
	fail("Unknown linext sampling algorithm name '", name, "'");
	return DefaultLinextSamplingAlgorithm;
}

template <int W, typename F>
void sampleLinext(const PartialOrder<W>& po, BitSet<W> verts, int sampleCount, LinextSamplingAlgorithm algo, F f) {
	if(algo == LinextSamplingAlgorithm::SwapLinextSampler) {
		sampleLinextSwap(po, verts, sampleCount, f);
	} else if(algo == LinextSamplingAlgorithm::GibbsLinextSampler) {
		sampleLinextGibbs(po, verts, sampleCount, f);
	} else {
		fail("sampleLinext: Unknown algorithm");
	}
}
template <int W, typename F>
void sampleLinext(const PartialOrder<W>& po, int sampleCount, LinextSamplingAlgorithm algo, F f) {
	sampleLinext(po, po.vertMask(), sampleCount, algo, f);
}
template <int W, typename F>
void sampleLinext(const PartialOrder<W>& po, BitSet<W> verts, int sampleCount, F f) {
	sampleLinext(po, verts, sampleCount, DefaultLinextSamplingAlgorithm, f);
}
template <int W, typename F>
void sampleLinext(const PartialOrder<W>& po, int sampleCount, F f) {
	sampleLinext(po, po.vertMask(), sampleCount, f);
}

template <int W>
double estimateOrderProb(const PartialOrder<W>& po, BitSet<W> verts, int a, int b, int sampleCount, LinextSamplingAlgorithm algo) {
	verts &= po.vertMask();
	if(!verts[a] || !verts[b]) {
		fail("estimateOrderProb: given vertices not in the set of vertices");
	}
	int s = verts.count();
	int abCount = 0;
	sampleLinext(po, verts, sampleCount, algo, [&](const int* order) {
		for(int i = 0; i < s; ++i) {
			if(order[i] == a) {
				++abCount;
				break;
			}
			if(order[i] == b) {
				break;
			}
		}
	});
	return (double)abCount / (double)sampleCount;
}
template <int W>
double estimateOrderProb(const PartialOrder<W>& po, int a, int b, int sampleCount, LinextSamplingAlgorithm algo) {
	return estimateOrderProb(po, po.vertMask(), a, b, sampleCount, algo);
}
template <int W>
double estimateOrderProb(const PartialOrder<W>& po, BitSet<W> verts, int a, int b, int sampleCount) {
	return estimateOrderProb(po, verts, a, b, sampleCount, DefaultLinextSamplingAlgorithm);
}
template <int W>
double estimateOrderProb(const PartialOrder<W>& po, int a, int b, int sampleCount) {
	return estimateOrderProb(po, po.vertMask(), a, b, sampleCount);
}
