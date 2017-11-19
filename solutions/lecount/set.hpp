#ifndef SET_H
#define SET_H

#include <cstdint>

#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>

struct dynamic_set
{
	boost::dynamic_bitset<> dbs;
	
	dynamic_set() {}
	dynamic_set(unsigned n) : dbs(n) {}
	dynamic_set(boost::dynamic_bitset<> dbs) : dbs(dbs) {}
	
	bool has(unsigned e) const
	{
		return dbs[e];
	}
	
	void set(unsigned e)
	{
		dbs[e] = true;
	}
	
	dynamic_set operator| (unsigned e)
	{
		boost::dynamic_bitset<> dbse = dbs;
		dbse[e] = true;
		return dbse;
	}
	
	dynamic_set operator^ (unsigned e)
	{
		boost::dynamic_bitset<> dbse = dbs;
		dbse.flip(e);
		return dbse;
	}
	
	void operator^= (unsigned e)
	{
		dbs[e] ^= true;
	}
	
	void flip(unsigned e)
	{
		dbs[e] ^= true;
	}
	
	bool operator== (const dynamic_set& other) const
	{
		return dbs == other.dbs;
	}
	
	bool operator!= (const dynamic_set& other) const
	{
		return dbs != other.dbs;
	}
	
	dynamic_set operator| (const dynamic_set& other) const
	{
		return dbs | other.dbs;
	}
	
	dynamic_set operator& (const dynamic_set& other) const
	{
		return dbs & other.dbs;
	}
	
	dynamic_set operator^ (const dynamic_set& other) const
	{
		return dbs ^ other.dbs;
	}
	
	dynamic_set operator- (const dynamic_set& other) const
	{
		return dbs - other.dbs;
	}
	
	void operator&= (const dynamic_set& other)
	{
		dbs &= other.dbs;
	}
	
	void operator|= (const dynamic_set& other)
	{
		dbs |= other.dbs;
	}
	
	void operator^= (const dynamic_set& other)
	{
		dbs ^= other.dbs;
	}
	
	dynamic_set operator~ () const
	{
		return ~dbs;
	}
	
	dynamic_set complement(unsigned) const
	{
		return ~*this;
	}
	
	unsigned count() const
	{
		return dbs.count();
	}
	
	unsigned count(unsigned) const
	{
		return count();
	}
	
	static dynamic_set empty(unsigned n)
	{
		return dynamic_set(n);
	}
	
	static dynamic_set complete(unsigned n)
	{
		boost::dynamic_bitset<> dbse(n);
		return dbse.set();
	}
	
	bool is_empty() const
	{
		return dbs.none();
	}
	
	bool operator[] (unsigned e) const
	{
		return has(e);
	}
	
	int first_element() const
	{
        return dbs.find_first();
	}
	
	int next_element()
	{
		int v = first_element();
		flip(v);
		return v;
	}
	
	unsigned cardinality(int n)
	{
		unsigned count = 0;
		for (int i = 0; i < n; i++) {
			if (has(i)) count++;
		}
		return count;
	}
	
	void print(int k)
	{
		for (int e = 0; e < k; e++) {
			printf("%i", has(e) ? 1 : 0);
		}
	}
	
	void println(int k)
	{
		print(k);
		putchar('\n');
	}
};

namespace std {
	template <>
	struct hash<dynamic_set> {
	public:
		size_t operator() (const dynamic_set& set) const {
			return boost::hash_value(set.dbs.m_bits);
		}
	};
}






template <typename T>
struct int_set
{
	T bits;
	
	int_set() {}
	int_set(T bits) : bits (bits) {}
	
	int_set operator& (int_set S)
	{
		return int_set(bits & S.bits);
	}
	
	int_set operator| (int_set S)
	{
		return int_set(bits | S.bits);
	}
	
	int_set operator^ (int_set S)
	{
		return int_set(bits ^ S.bits);
	}
	
	int_set operator| (unsigned e)
	{
		return int_set(bits | sing(e));
	}
	
	int_set operator^ (unsigned e)
	{
		return int_set(bits ^ sing(e));
	}
	
	int_set operator- (int_set S)
	{
		return int_set(bits & ~S.bits);
	}
	
	void operator&= (int_set S)
	{
		bits &= S.bits;
	}
	
	void operator|= (int_set S)
	{
		bits |= S.bits;
	}
	
	void operator^= (int_set S)
	{
		bits ^= S.bits;
	}
	
	void operator^= (unsigned e)
	{
		flip(e);
	}
	
	bool has(unsigned e) const
	{
		return bits & sing(e);
	}
	
	void set(unsigned e)
	{
		bits |= sing(e);
	}
	
	int_set complement(unsigned n)
	{
		return (~bits) & (((T)-1) >> (sizeof(T)*8 - n));
	}
	
	void flip(unsigned e)
	{
		bits ^= sing(e);
	}
	
	static T sing(unsigned e)
	{
		return (T)1 << e;
	}
	
	static int_set complete(unsigned n)
	{
		return n == 0 ? 0 : ((T)-1) >> (sizeof(T)*8 - n);
	}
	
	bool operator== (const int_set& other) const
	{
		return bits == other.bits;
	}
	
	bool operator!= (const int_set& other) const
	{
		return bits != other.bits;
	}
	
	static int_set empty(unsigned)
	{
		return int_set(0);
	}
	
	unsigned count(unsigned n) const
	{
		unsigned c = 0;
		for (unsigned k = 0; k < n; k++) {
			if (has(k)) c++;
		}
		return c;
	}
	
	bool is_empty() const
	{
		return bits == 0;
	}
	
	unsigned cardinality(int)
	{
		return __builtin_popcountll(bits);
	}
	
	bool operator[] (unsigned e) const
	{
		return has(e);
	}
	
	int first_element() const
	{
		return ffsll(bits) - 1;
	}
	
	int next_element()
	{
		int v = first_element();
		flip(v);
		return v;
	}
	
	unsigned get_elements(unsigned n, unsigned *elements)
	{
		unsigned k = 0;
		for (unsigned i = 0; i < n; ++i) {
			if (has(i)) elements[k++] = i;
		}
		
		return k;
	}
	
	void print(int k)
	{
		for (int e = 0; e < k; e++) {
			printf("%i", has(e) ? 1 : 0);
		}
	}
	
	void println(int k)
	{
		print(k);
		putchar('\n');
	}
};

namespace std {
	template <class T>
	struct hash<int_set<T>> {
	public:
		size_t operator() (const int_set<T>& set) const {
			std::hash<T> hash;
			return hash(set.bits);
		}
	};
}

typedef int_set<uint32_t> set32;
typedef int_set<uint64_t> set64;







// a set of size 64 * N
template <size_t N>
struct fixed_set
{
	long long unsigned int chunks[N];
	
	fixed_set()
	{
	}
	
	fixed_set(const fixed_set<N> &other)
	{
		for (unsigned i = 0; i < N; i++) {
			chunks[i] = other.chunks[i];
		}
	}
	
	static long long unsigned singleton(unsigned e)
	{
		return ((long long unsigned)1) << e;
	}
	
	
	void flip(unsigned e)
	{
		unsigned i = e >> 6;        // i = e / 64
		unsigned j = e - i * 64;    // j = e % 64
		chunks[i] ^= singleton(j);
	}
	
	bool has(unsigned e) const
	{
		unsigned i = e >> 6;        // i = e / 64
		unsigned j = e - i * 64;    // j = e % 64
		return chunks[i] & singleton(j);
	}
	
	void set(unsigned e)
	{
		unsigned i = e >> 6;        // i = e / 64
		unsigned j = e - i * 64;    // j = e % 64
		chunks[i] |= singleton(j);
	}
	
	fixed_set<N> operator^ (unsigned e) const
	{
		fixed_set<N> set(*this);
		set.flip(e);
		return set;
	}
	
	fixed_set<N> operator| (unsigned e) const
	{
		fixed_set<N> set(*this);
		set.set(e);
		return set;
	}
	
	bool operator== (const fixed_set<N> &other) const
	{
		for (unsigned c = 0; c < N; c++) {
			if (chunks[c] != other.chunks[c]) return false;
		}
		return true;
	}
	
	bool operator!= (const fixed_set<N> &other) const
	{
		for (unsigned c = 0; c < N; c++) {
			if (chunks[c] != other.chunks[c]) return true;
		}
		return false;
	}
	
	void operator&= (const fixed_set<N> &other)
	{
		for (unsigned c = 0; c < N; c++) {
			chunks[c] &= other.chunks[c];
		}
	}
	
	void operator|= (const fixed_set<N> &other)
	{
		for (unsigned c = 0; c < N; c++) {
			chunks[c] |= other.chunks[c];
		}
	}
	
	void operator^= (const fixed_set<N> &other)
	{
		for (unsigned c = 0; c < N; c++) {
			chunks[c] ^= other.chunks[c];
		}
	}
	
	fixed_set<N> operator^ (const fixed_set<N> &other) const
	{
		fixed_set<N> result(*this);
		for (unsigned c = 0; c < N; c++) {
			result.chunks[c] ^= other.chunks[c];
		}
		return result;
	}
	
	fixed_set<N> operator& (const fixed_set<N> &other) const
	{
		fixed_set<N> result(*this);
		for (unsigned c = 0; c < N; c++) {
			result.chunks[c] &= other.chunks[c];
		}
		return result;
	}
	
	fixed_set<N> operator| (const fixed_set<N> &other) const
	{
		fixed_set<N> result(*this);
		for (unsigned c = 0; c < N; c++) {
			result.chunks[c] |= other.chunks[c];
		}
		return result;
	}
	
	fixed_set<N> operator- (const fixed_set<N> &other) const
	{
		fixed_set<N> result(*this);
		for (unsigned c = 0; c < N; c++) {
			result.chunks[c] &= ~other.chunks[c];
		}
		return result;
	}
	
	void operator^= (unsigned e)
	{
		flip(e);
	}
	
	void operator|= (unsigned e)
	{
		set(e);
	}
	
	static fixed_set<N> empty()
	{
		fixed_set<N> set;
		for (unsigned i = 0; i < N; i++) {
			set.chunks[i] = 0;
		}
		return set;
	}
	
	static fixed_set<N> empty(unsigned)
	{
		return empty();
	}
	
	static fixed_set<N> complete(unsigned n)
	{
		fixed_set<N> set = fixed_set<N>::empty();
		for (unsigned e = 0; e < n; e++) {
			set.flip(e);
		}
		
		return set;
	}
	
	bool is_empty() const
	{
		for (unsigned c = 0; c < N; c++) {
			if (chunks[c] != 0) return false;
		}
		return true;
	}
	
	unsigned count(unsigned n) const
	{
		unsigned k = 0;
		for (unsigned e = 0; e < n; e++) {
			if (has(e)) k++;
		}
		return k;
	}
	
	unsigned cardinality(unsigned) const
	{
		unsigned ones = 0;
		for (unsigned c = 0; c < N; c++) {
			ones += __builtin_popcountll(chunks[c]);
		}
		
		return ones;
	}
	
	bool operator[] (unsigned e) const
	{
		return has(e);
	}
	
	int first_element() const
	{
		for (unsigned i = 0; i < N; i++) {
			for (unsigned j = 0; j < 64; ++j) {
				if (chunks[i]) return i * 64 + ffsll(chunks[i]) - 1;
			}
		}
		
		return -1;
	}
	
	int next_element()
	{
		int v = first_element();
		flip(v);
		return v;
	}
	
	unsigned get_elements(unsigned n, unsigned *elements)
	{
		unsigned k = 0;
		for (unsigned i = 0; i < n; ++i) {
			if (has(i)) elements[k++] = i;
		}
		
		return k;
	}
	
	
	void print(int k) const
	{
		for (int e = 0; e < k; e++) {
			printf("%i", has(e) ? 1 : 0);
		}
	}
	
	void println(int k) const
	{
		print(k);
		putchar('\n');
	}
};


namespace std {
	template <size_t N>
	struct hash<fixed_set<N>> {
	public:
		size_t operator() (const fixed_set<N>& set) const {
			return boost::hash_range(set.chunks, set.chunks + N);
		}
	};
}



#endif
