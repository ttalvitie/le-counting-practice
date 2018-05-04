#pragma once

#include "partial_order.hpp"

double logNCR(double n, double k) {
	return lgamma(n + 1.0) - (lgamma(k + 1.0) + lgamma(n - k + 1.0));
}

typedef pair<double, vector<pair<int, int>>> SubRelaxationResult;

template <int W>
bool isTreeComponent(const PartialOrder<W>& po, BitSet<W> verts) {
	int count = 0;
	verts.iterate([&](int v) {
		BitSet<W> seen;
		po.topo(verts & po.succ(v), [&](int x) {
			if(!seen[x]) {
				++count;
				seen |= po.succ(x);
			}
		});
	});
	return count == verts.count() - 1;
}

template <int W>
pair<bool, SubRelaxationResult> tryComputeRelaxationByPerfectBipartiteDivision(
	const PartialOrder<W>& po,
	BitSet<W> verts
) {
	double best = -1.0;
	BitSet<W> bestLeft;
	BitSet<W> bestRight;
	
	BitSet<W> left;
	BitSet<W> right = verts;
	po.topo(verts, [&](int v) {
		left.add(v);
		right &= po.succ(v);
		if(right && (left | right) == verts) {
			double val = logNCR(verts.count(), left.count());
			if(val > best) {
				best = val;
				bestLeft = left;
				bestRight = right;
			}
		}
	});
	if(best == -1.0) {
		return make_pair(false, SubRelaxationResult());
	}
	
	left = bestLeft;
	right = bestRight;
	
	double leftLogCount, rightLogCount;
	vector<pair<int, int>> leftEdges, rightEdges;
	tie(leftLogCount, leftEdges) = computeRelaxation(po, left);
	tie(rightLogCount, rightEdges) = computeRelaxation(po, right);
	
	double logCount = leftLogCount + rightLogCount;
	vector<pair<int, int>> edges;
	edges.insert(edges.end(), leftEdges.begin(), leftEdges.end());
	edges.insert(edges.end(), rightEdges.begin(), rightEdges.end());
	
	left.iterate([&](int a) {
		right.iterate([&](int b) {
			edges.emplace_back(a, b);
		});
	});
	
	return make_pair(true, SubRelaxationResult(logCount, edges));
}

template <int W>
pair<bool, SubRelaxationResult> tryComputeExactLinextCount(
	const PartialOrder<W>& po,
	BitSet<W> verts
) {
	long long expansionsLeft = 10000;
	double linextCount = po.tryLinextCount(verts, expansionsLeft);
	if(!expansionsLeft || !isfinite(linextCount)) {
		return make_pair(false, SubRelaxationResult());
	}
	
	double logCount = log(linextCount);
	vector<pair<int, int>> edges;
	verts.iterate([&](int a) {
		(po.succ(a) & verts).iterate([&](int b) {
			edges.emplace_back(a, b);
		});
	});
	
	return make_pair(true, SubRelaxationResult(logCount, edges));
}

template <int W>
SubRelaxationResult computeRelaxationByImperfectBipartiteDivision(
	const PartialOrder<W>& po,
	BitSet<W> verts
) {
	double best = -1.0;
	BitSet<W> bestLeft;
	BitSet<W> bestRight;
	
	auto account = [&](BitSet<W> left, BitSet<W> right) {
		if(left && right) {
			double val = logNCR(left.count() + right.count(), left.count());
			if(val > best) {
				best = val;
				bestLeft = left;
				bestRight = right;
			}
			return val;
		} else {
			return -1.0;
		}
	};
	
	verts.iterate([&](int s) {
		BitSet<W> left;
		left.add(s);
		BitSet<W> right = po.succ(s) & verts;
		
		account(left, right);
		
		while(true) {
			double localBest = -1.0;
			int localBestAdd = -1;
			
			(verts & ~left).iterate([&](int v) {
				BitSet<W> newLeft = left;
				newLeft.add(v);
				BitSet<W> newRight = right & po.succ(v);
				double newVal = account(newLeft, newRight);
				if(newVal > localBest) {
					localBest = newVal;
					localBestAdd = v;
				}
			});
			
			if(localBestAdd == -1) {
				break;
			}
			
			left.add(localBestAdd);
			right &= po.succ(localBestAdd);
		}
	});
	
	assert(best != -1.0);
	
	BitSet<W> left = bestLeft;
	BitSet<W> right = bestRight;
	BitSet<W> other = verts & ~(left | right);
	
	double leftLogCount, rightLogCount, otherLogCount;
	vector<pair<int, int>> leftEdges, rightEdges, otherEdges;
	tie(leftLogCount, leftEdges) = computeRelaxation(po, left);
	tie(rightLogCount, rightEdges) = computeRelaxation(po, right);
	tie(otherLogCount, otherEdges) = computeRelaxation(po, other);
	
	double logCount = leftLogCount + rightLogCount + otherLogCount;
	logCount += logNCR(verts.count(), other.count());
	
	vector<pair<int, int>> edges;
	edges.insert(edges.end(), leftEdges.begin(), leftEdges.end());
	edges.insert(edges.end(), rightEdges.begin(), rightEdges.end());
	edges.insert(edges.end(), otherEdges.begin(), otherEdges.end());
	
	left.iterate([&](int a) {
		right.iterate([&](int b) {
			edges.emplace_back(a, b);
		});
	});
	
	return SubRelaxationResult(logCount, edges);
}

struct UnionFind {
	UnionFind(int n) : parent(n) {
		for(int i = 0; i < n; ++i) {
			parent[i] = i;
		}
	}
	
	int find(int x) {
		if(parent[x] != x) {
			parent[x] = find(parent[x]);
		}
		return parent[x];
	}
	
	void merge(int x, int y) {
		x = find(x);
		y = find(y);
		parent[x] = y;
	}
	
	vector<int> parent;
};

template <int W>
vector<long double> treeSpectrum(const vector<vector<pair<int, bool>>>& tree, int a, BitSet<W> nogo) {
	int b = -1;
	bool forward;
	for(auto edge : tree[a]) {
		if(!nogo[edge.first]) {
			b = edge.first;
			forward = edge.second;
			break;
		}
	}
	if(b == -1) {
		vector<long double> ret;
		ret.push_back(1.0l);
		return ret;
	}
	
	BitSet<W> newNogoA = nogo;
	newNogoA.add(b);
	vector<long double> A = treeSpectrum(tree, a, newNogoA);
	BitSet<W> newNogoB;
	newNogoB.add(a);
	vector<long double> B = treeSpectrum(tree, b, newNogoB);
	
	int As = A.size();
	int Bs = B.size();
	int Rs = As + Bs;
	
	vector<long double> R(Rs);
	
	if(forward) {
		for(int i = Bs - 2; i >= 0; --i) {
			B[i] += B[i + 1];
		}
		
		for(int r = 0; r < Rs; ++r) {
			long double val = 0.0;
			for(int i = max(0, r - Bs + 1); i < As && i < r + 1; ++i) {
				val += A[i] * exp((long double)(logNCR(r, i) + logNCR(Rs - r - 1, As - i - 1))) * B[r - i];
			}
			R[r] = val;
		}
	} else {
		for(int i = 1; i < Bs; ++i) {
			B[i] += B[i - 1];
		}
		
		for(int r = 0; r < Rs; ++r) {
			long double val = 0.0;
			for(int i = max(0, r - Bs); i < As && i < r; ++i) {
				val += A[i] * exp((long double)(logNCR(r, i) + logNCR(Rs - r - 1, As - i - 1))) * B[r - i - 1];
			}
			R[r] = val;
		}
	}
	
	return R;
}

template <int W>
long double treeLinextCount(const vector<vector<pair<int, bool>>>& tree, int s) {
	vector<long double> spectrum = treeSpectrum<W>(tree, s, BitSet<W>());
	long double ret = 0.0;
	for(long double x : spectrum) {
		ret += x;
	}
	return ret;
}

template <int W>
SubRelaxationResult computeTreeRelaxation(
	const PartialOrder<W>& po,
	BitSet<W> verts
) {
	auto oneTry = [&]() {
		vector<pair<int, int>> allEdges;
		vector<pair<int, int>> edges;
		
		verts.iterate([&](int v) {
			BitSet<W> seen;
			po.topo(verts & po.succ(v), [&](int x) {
				if(!seen[x]) {
					allEdges.emplace_back(v, x);
					seen |= po.succ(x);
				}
			});
		});
		
		shuffle(allEdges.begin(), allEdges.end(), rng);
		
		vector<vector<pair<int, bool>>> tree(po.size());
		
		UnionFind uf(po.size());
		for(auto edge : allEdges) {
			if(uf.find(edge.first) != uf.find(edge.second)) {
				edges.emplace_back(edge);
				tree[edge.first].emplace_back(edge.second, true);
				tree[edge.second].emplace_back(edge.first, false);
				uf.merge(edge.first, edge.second);
			}
		}
		
		double logCount = log(treeLinextCount<W>(tree, verts.min()));
		
		return SubRelaxationResult(logCount, edges);
	};
	
	SubRelaxationResult ret = oneTry();
	for(int i = 0; i < 5; ++i) {
		SubRelaxationResult ret2 = oneTry();
		if(ret2.first < ret.first) {
			swap(ret, ret2);
		}
	}
	return ret;
}

template <int W>
SubRelaxationResult computeRelaxationInComponent(
	const PartialOrder<W>& po,
	BitSet<W> verts
) {
	if(verts.count() <= 1) {
		return SubRelaxationResult(0.0, vector<pair<int, int>>());
	}
	
	{
		bool success;
		SubRelaxationResult result;
		
		tie(success, result) = tryComputeRelaxationByPerfectBipartiteDivision(po, verts);
		if(success) {
			return result;
		}
		
		if(isTreeComponent(po, verts)) {
			return computeTreeRelaxation(po, verts);
		}
		
		tie(success, result) = tryComputeExactLinextCount(po, verts);
		if(success) {
			return result;
		}
	}
	
	SubRelaxationResult a = computeRelaxationByImperfectBipartiteDivision(po, verts);
	SubRelaxationResult b = computeTreeRelaxation(po, verts);
	if(a.first < b.first) {
		return a;
	} else {
		return b;
	}
}

template <int W>
SubRelaxationResult computeRelaxation(
	const PartialOrder<W>& po,
	BitSet<W> verts
) {
	double logCount = 0.0;
	vector<pair<int, int>> edges;
	
	int vertsLeft = verts.count();
	po.components(verts, [&](BitSet<W> component) {
		double componentLogCount;
		vector<pair<int, int>> componentEdges;
		tie(componentLogCount, componentEdges) = computeRelaxationInComponent(po, component);
		
		logCount += componentLogCount;
		edges.insert(edges.end(), componentEdges.begin(), componentEdges.end());
		
		logCount += logNCR(vertsLeft, component.count());
		vertsLeft -= component.count();
	});
	
	return make_pair(logCount, edges);
}

template <int W>
pair<PartialOrder<W>, double> computeRelaxation(const PartialOrder<W>& po) {
	double logCount;
	vector<pair<int, int>> edges;
	tie(logCount, edges) = computeRelaxation(po, po.vertMask());
	
	PartialOrder<W> ret(po.size());
	for(auto p : edges) {
		ret.add(p.first, p.second);
	}
	
	return make_pair(ret, logCount);
}
