#pragma once

#include <unordered_map>
#include <functional>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <limits>
#include <random>
#include <vector>
#include <array>
#include <queue>

using namespace std;

extern thread_local mt19937 rng;

inline void stderrPrint() {
	cerr << "\n";
}
template <typename F, typename... T>
void stderrPrint(const F& f, const T&... p) {
	cerr << f;
	stderrPrint(p...);
}
template <typename... T>
void fail(const T&... p) {
	stderrPrint("FAIL: ", p...);
	abort();
}

template <typename T>
T fromString(const string& str) {
	T ret;
	stringstream ss(str);
	ss >> ret;
	if(ss.fail()) fail("fromString: Could not convert string '", str, "' to given type (typeid name '", typeid(T).name(), "').");
	return ret;
}

template <typename T>
string toString(const T& obj) {
	stringstream ss;
	ss << obj;
	return ss.str();
}

inline constexpr uint64_t ones64(int count) {
	return (count == 64)
		? (uint64_t)-1
		: (((uint64_t)1 << count) - (uint64_t)1);
}

inline constexpr uint64_t bit64(int pos) {
	return (uint64_t)1 << pos;
}

inline constexpr int popcount64(uint64_t x) {
	static_assert(sizeof(unsigned long long) == sizeof(uint64_t), "The case unsigned long long != uint64_t is not implemented");
	return __builtin_popcountll(x);
}

inline constexpr int ctz64(uint64_t x) {
	static_assert(sizeof(unsigned long long) == sizeof(uint64_t), "The case unsigned long long != uint64_t is not implemented");
	return __builtin_ctzll(x);
}

inline constexpr int clz64(uint64_t x) {
	static_assert(sizeof(unsigned long long) == sizeof(uint64_t), "The case unsigned long long != uint64_t is not implemented");
	return __builtin_clzll(x);
}

template <typename F>
void iterateOnes64(uint64_t mask, F f) {
	while(mask) {
		int i = ctz64(mask);
		f(i);
		mask ^= bit64(i);
	}
}
template <typename F>
bool iterateOnes64And(uint64_t mask, F f) {
	while(mask) {
		int i = ctz64(mask);
		if(!(bool)f(i)) {
			return false;
		}
		mask ^= bit64(i);
	}
	return true;
}

inline int randomMaskElement64(uint64_t mask) {
	int count = popcount64(mask);
	if(count == 0) {
		fail("Cannot compute random element from an empty bit mask");
	}
	int ret = 0;
	auto filter = [&](uint64_t bitMask, int outBit) {
		int count1 = popcount64(mask & bitMask);
		if(uniform_int_distribution<int>(0, count - 1)(rng) >= count1) {
			mask &= ~bitMask;
			count -= count1;
		} else {
			mask &= bitMask;
			ret |= outBit;
			count = count1;
		}
	};
	
	filter(0xFFFFFFFF00000000, 32);
	filter(0xFFFF0000FFFF0000, 16);
	filter(0xFF00FF00FF00FF00, 8);
	filter(0xF0F0F0F0F0F0F0F0, 4);
	filter(0xCCCCCCCCCCCCCCCC, 2);
	filter(0xAAAAAAAAAAAAAAAA, 1);
	
	return ret;
}
