#pragma once

#include "bitset.hpp"

template <int W>
class PartialOrder {
public:
	static const int MaxVertCount = BitSet<W>::BitCount;
	
	PartialOrder() : size_(0), pred_(0), succ_(0) { }
	PartialOrder(int size) : size_(size), pred_(size), succ_(size) {
		if(size > MaxVertCount) {
			fail("PartialOrder too large for given BitSet size");
		}
		vertMask_ = BitSet<W>::ones(size);
	}
	
	int size() const {
		return size_;
	}
	BitSet<W> vertMask() const {
		return vertMask_;
	}
	
	BitSet<W> pred(int v) const {
		return pred_[v];
	}
	BitSet<W> succ(int v) const {
		return succ_[v];
	}
	
	bool operator()(int a, int b) const {
		return succ(a)[b];
	}
	
	void add(int a, int b) {
		if(a == b || (*this)(b, a)) {
			fail("PartialOrder: Cannot add edge");
		}
		if((*this)(a, b)) {
			return;
		}
		BitSet<W> A = pred(a);
		A.add(a);
		BitSet<W> B = succ(b);
		B.add(b);
		A.iterate([&](int v) {
			succ_[v] |= B;
		});
		B.iterate([&](int v) {
			pred_[v] |= A;
		});
	}
	
	template <typename F>
	void topo(BitSet<W> verts, F f) const {
		BitSet<W> todo = verts & vertMask_;
		while(todo) {
			topo_(f, todo, todo.min());
		}
	}
	template <typename F>
	void topo(F f) const {
		topo(vertMask_, f);
	}
	template <typename F>
	void reverseTopo(BitSet<W> verts, F f) const {
		BitSet<W> todo = verts & vertMask_;
		while(todo) {
			reverseTopo_(f, todo, todo.min());
		}
	}
	template <typename F>
	void reverseTopo(F f) const {
		reverseTopo(vertMask_, f);
	}
	
	template <typename F>
	void components(BitSet<W> verts, F f) const {
		BitSet<W> todo = verts & vertMask_;
		while(todo) {
			int v = todo.min();
			BitSet<W> seen;
			seen.add(v);
			BitSet<W> oldTodo = todo;
			while(true) {
				BitSet<W> d = seen & todo;
				if(!d) break;
				BitSet<W> dcopy = d;
				dcopy.iterate([&](int x) {
					todo.del(x);
					seen |= pred(x) | succ(x);
				});
			}
			f(todo ^ oldTodo);
		}
	}
	template <typename F>
	void components(F f) const {
		components(vertMask_, f);
	}
	
	BitSet<W> minimalVerts(BitSet<W> verts) const {
		verts &= vertMask_;
		BitSet<W> comp;
		verts.iterate([&](int v) {
			comp |= succ(v);
		});
		return verts & ~comp;
	}
	BitSet<W> minimalVerts() const {
		return minimalVerts(vertMask_);
	}
	BitSet<W> maximalVerts(BitSet<W> verts) const {
		verts &= vertMask_;
		BitSet<W> comp;
		verts.iterate([&](int v) {
			comp |= pred(v);
		});
		return verts & ~comp;
	}
	BitSet<W> maximalVerts() const {
		return maximalVerts(vertMask_);
	}
	
	void printDot(BitSet<W> verts, ostream& out = cout) const {
		out << "digraph {\n";
		verts.iterate([&](int v) {
			out << "  " << v + 1 << ";\n";
		});
		verts.iterate([&](int v) {
			BitSet<W> seen;
			topo(verts, [&](int x) {
				if(succ(v)[x] && !seen[x]) {
					out << "  " << v + 1 << " -> " << x + 1 << ";\n";
					seen |= succ(x);
				}
			});
		});
		out << "}\n";
	}
	void printDot(ostream& out = cout) const {
		printDot(vertMask_, out);
	}
	
	// Expand at most expansionsLeft states. If expansionsLeft is 0 in the end,
	// the computation did not succeed.
	double tryLinextCount(BitSet<W> verts, long long& expansionsLeft) const {
		verts &= vertMask_;
		unordered_map<BitSet<W>, double> mem;
		double result = 1.0;
		double logcoef = lgamma(verts.count() + 1);
		components(verts, [&](BitSet<W> comp) {
			mem.clear();
			
			if(minimalVerts(comp).count() < maximalVerts(comp).count()) {
				result *= linextCountRecursion_<&PartialOrder::minimalVerts>(mem, comp, expansionsLeft);
			} else {
				result *= linextCountRecursion_<&PartialOrder::maximalVerts>(mem, comp, expansionsLeft);
			}
			
			logcoef -= lgamma(comp.count() + 1);
		});
		return result * exp(logcoef);
	}
	
	double linextCount(BitSet<W> verts) const {
		long long expansionsLeft = -1;
		return tryLinextCount(verts, expansionsLeft);
	}
	double linextCount() const {
		return linextCount(vertMask_);
	}
	
private:
	template <typename F>
	void topo_(F& f, BitSet<W>& todo, int v) const {
		while(true) {
			BitSet<W> d = pred(v) & todo;
			if(!d) {
				break;
			}
			topo_(f, todo, d.min());
		}
		todo.del(v);
		f(v);
	}
	template <typename F>
	void reverseTopo_(F& f, BitSet<W>& todo, int v) const {
		while(true) {
			BitSet<W> d = succ(v) & todo;
			if(!d) {
				break;
			}
			reverseTopo_(f, todo, d.min());
		}
		todo.del(v);
		f(v);
	}
	
	template <BitSet<W> (PartialOrder::*RemovalFunc)(BitSet<W>) const>
	double linextCountRecursion_(
		unordered_map<BitSet<W>, double>& mem,
		BitSet<W> verts,
		long long& expansionsLeft
	) const {
		double result = 1.0;
		double logcoef = lgamma(verts.count() + 1);
		components(verts, [&](BitSet<W> comp) {
			auto it = mem.find(comp);
			if(it == mem.end()) {
				if(!expansionsLeft) {
					return;
				}
				--expansionsLeft;
				double val = 0.0;
				(this->*RemovalFunc)(comp).iterate([&](int v) {
					BitSet<W> subComp = comp;
					subComp.del(v);
					val += linextCountRecursion_<RemovalFunc>(mem, subComp, expansionsLeft);
				});
				it = mem.emplace(comp, val).first;
			}
			result *= it->second;
			
			logcoef -= lgamma(comp.count() + 1);
		});
		return result * exp(logcoef);
	}
	
	int size_;
	BitSet<W> vertMask_;
	vector<BitSet<W>> pred_;
	vector<BitSet<W>> succ_;
};

template <int W>
bool operator==(const PartialOrder<W>& a, const PartialOrder<W>& b) {
	if(a.size() != b.size()) {
		return false;
	}
	int n = a.size();
	for(int v = 0; v < n; ++v) {
		if(a.succ(v) != b.succ(v)) {
			return false;
		}
	}
	return true;
}
template <int W>
bool operator!=(const PartialOrder<W>& a, const PartialOrder<W>& b) {
	return !(a == b);
}

template <typename F>
void readPartialOrder(istream& in, F f) {
	string line;
	getline(in, line);
	stringstream liness(line);
	
	vector<pair<int, int>> edges;
	
	int n = 0;
	int x;
	while(liness >> x) {
		if(x) {
			edges.emplace_back(0, n);
		}
		++n;
	}
	
	for(int a = 1; a < n; ++a) {
		for(int b = 0; b < n; ++b) {
			if(!(in >> x)) fail("Invalid adjacency matrix");
			if(x) {
				edges.emplace_back(a, b);
			}
		}
	}
	
	if(n <= BitSet<1>::BitCount) {
		PartialOrder<1> po(n);
		for(pair<int, int> p : edges) {
			po.add(p.first, p.second);
		}
		f(move(po));
		return;
	}
	if(n <= BitSet<2>::BitCount) {
		PartialOrder<2> po(n);
		for(pair<int, int> p : edges) {
			po.add(p.first, p.second);
		}
		f(move(po));
		return;
	}
	if(n <= BitSet<3>::BitCount) {
		PartialOrder<3> po(n);
		for(pair<int, int> p : edges) {
			po.add(p.first, p.second);
		}
		f(move(po));
		return;
	}
	if(n <= BitSet<4>::BitCount) {
		PartialOrder<4> po(n);
		for(pair<int, int> p : edges) {
			po.add(p.first, p.second);
		}
		f(move(po));
		return;
	}
	if(n <= BitSet<5>::BitCount) {
		PartialOrder<5> po(n);
		for(pair<int, int> p : edges) {
			po.add(p.first, p.second);
		}
		f(move(po));
		return;
	}
	if(n <= BitSet<6>::BitCount) {
		PartialOrder<6> po(n);
		for(pair<int, int> p : edges) {
			po.add(p.first, p.second);
		}
		f(move(po));
		return;
	}
	if(n <= BitSet<7>::BitCount) {
		PartialOrder<7> po(n);
		for(pair<int, int> p : edges) {
			po.add(p.first, p.second);
		}
		f(move(po));
		return;
	}
	if(n <= BitSet<8>::BitCount) {
		PartialOrder<8> po(n);
		for(pair<int, int> p : edges) {
			po.add(p.first, p.second);
		}
		f(move(po));
		return;
	}
	fail("Partial order is too large");
}
template <typename F>
void readPartialOrder(F f) {
	readPartialOrder(cin, f);
}
