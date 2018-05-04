#pragma once

#include "common.hpp"

template <int W>
struct BitSet {
	static const int BitCount = W * 64;
	
	uint64_t words[W];
	
	BitSet<W>() : words{} { }
	
	static BitSet<W> ones() {
		BitSet<W> ret;
		for(int i = 0; i < W; ++i) {
			ret.words[i] = (uint64_t)-1;
		}
		return ret;
	}
	
	static BitSet<W> ones(int count) {
		if(count >= BitCount) {
			return ones();
		}
		BitSet<W> ret;
		int a = count >> 6;
		int b = count & 63;
		for(int i = 0; i < a; ++i) {
			ret.words[i] = (uint64_t)-1;
		}
		ret.words[a] = (bit64(b) - 1);
		return ret;
	}
	
	static BitSet<W> randomMask() {
		BitSet<W> ret;
		for(int i = 0; i < W; ++i) {
			ret.words[i] = uniform_int_distribution<uint64_t>()(rng);
		}
		return ret;
	}
	
	bool operator[](int idx) const {
		int a = idx >> 6;
		int b = idx & 63;
		return (bool)(words[a] & bit64(b));
	}
	
	void add(int idx) {
		int a = idx >> 6;
		int b = idx & 63;
		words[a] |= bit64(b);
	}
	
	void del(int idx) {
		int a = idx >> 6;
		int b = idx & 63;
		words[a] &= ~bit64(b);
	}
	
	void set(int idx, bool val) {
		int a = idx >> 6;
		int b = idx & 63;
		words[a] &= ~bit64(b);
		words[a] |= (uint64_t)val << b;
	}
	
	BitSet<W>& operator&=(BitSet<W> x) {
		for(int i = 0; i < W; ++i) {
			words[i] &= x.words[i];
		}
		return *this;
	}
	BitSet<W>& operator|=(BitSet<W> x) {
		for(int i = 0; i < W; ++i) {
			words[i] |= x.words[i];
		}
		return *this;
	}
	BitSet<W>& operator^=(BitSet<W> x) {
		for(int i = 0; i < W; ++i) {
			words[i] ^= x.words[i];
		}
		return *this;
	}
	
	BitSet<W> operator~() const {
		BitSet<W> ret;
		for(int i = 0; i < W; ++i) {
			ret.words[i] = ~words[i];
		}
		return ret;
	}
	
	template <typename F>
	void iterate(F f) const {
		for(int i = 0; i < W; ++i) {
			iterateOnes64(words[i], [&](int j) {
				f((i << 6) + j);
			});
		}
	}
	template <typename F>
	bool iterateAnd(F f) const {
		for(int i = 0; i < W; ++i) {
			if(!iterateOnes64And(words[i], [&](int j) {
				return (bool)f((i << 6) + j);
			})) {
				return false;
			}
		}
		return true;
	}
	
	int min() const {
		for(int i = 0; i < W; ++i) {
			if(words[i]) {
				return (i << 6) + ctz64(words[i]);
			}
		}
		return -1;
	}
	
	int max() const {
		for(int i = W - 1; i >= 0; --i) {
			if(words[i]) {
				return (i << 6) + 63 - clz64(words[i]);
			}
		}
		return -1;
	}
	
	int count() const {
		int ret = 0;
		for(int i = 0; i < W; ++i) {
			ret += popcount64(words[i]);
		}
		return ret;
	}
	
	explicit operator bool() const {
		for(int i = 0; i < W; ++i) {
			if(words[i]) {
				return true;
			}
		}
		return false;
	}
	
	int random() const {
		int t = uniform_int_distribution<int>(0, count() - 1)(rng);
		for(int i = 0; i < W; ++i) {
			int c = popcount64(words[i]);
			if(t < c) {
				return (i << 6) + randomMaskElement64(words[i]);
			}
			t -= c;
		}
		return -1;
	}
};

template <int W>
BitSet<W> operator&(BitSet<W> a, BitSet<W> b) {
	a &= b;
	return a;
}
template <int W>
BitSet<W> operator|(BitSet<W> a, BitSet<W> b) {
	a |= b;
	return a;
}
template <int W>
BitSet<W> operator^(BitSet<W> a, BitSet<W> b) {
	a ^= b;
	return a;
}

template <int W>
bool operator==(BitSet<W> a, BitSet<W> b) {
	for(int i = 0; i < W; ++i) {
		if(a.words[i] != b.words[i]) {
			return false;
		}
	}
	return true;
}
template <int W>
bool operator!=(BitSet<W> a, BitSet<W> b) {
	for(int i = 0; i < W; ++i) {
		if(a.words[i] != b.words[i]) {
			return true;
		}
	}
	return false;
}

template <int W>
bool operator<(BitSet<W> a, BitSet<W> b) {
	for(int i = W - 1; i >= 0; --i) {
		if(a.words[i] != b.words[i]) {
			return a.words[i] < b.words[i];
		}
	}
	return false;
}
template <int W>
bool operator<=(BitSet<W> a, BitSet<W> b) {
	for(int i = W - 1; i >= 0; --i) {
		if(a.words[i] != b.words[i]) {
			return a.words[i] < b.words[i];
		}
	}
	return true;
}
template <int W>
bool operator>(BitSet<W> a, BitSet<W> b) {
	for(int i = W - 1; i >= 0; --i) {
		if(a.words[i] != b.words[i]) {
			return a.words[i] > b.words[i];
		}
	}
	return false;
}
template <int W>
bool operator>=(BitSet<W> a, BitSet<W> b) {
	for(int i = W - 1; i >= 0; --i) {
		if(a.words[i] != b.words[i]) {
			return a.words[i] > b.words[i];
		}
	}
	return true;
}

namespace std {
template <int W>
struct hash<BitSet<W>> {
	size_t operator()(BitSet<W> x) const {
		size_t ret = 0;
		for(int i = 0; i < W; ++i) {
			ret = 31 * ret + (size_t)x.words[i];
		}
		return ret;
	}
};
}
