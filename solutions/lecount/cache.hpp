#ifndef CACHE_H
#define CACHE_H

static unsigned hash_primes[] = {
	53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317,
	196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843,
	50331653, 100663319, 201326611, 402653189, 805306457, 1610612741
};

struct CacheOptions
{
	int initial_prime_index;
	int prime_index_increment;
	double max_load_factor;
	
	int initial_array_size;
	int array_resize_factor;
	
	bool verbose;
	
	CacheOptions()
	{
		initial_prime_index = 4;
		prime_index_increment = 3;
		max_load_factor = 0.5;
		
		initial_array_size = 1024;
		array_resize_factor = 4;
		
		verbose = false;
	}
};

struct CacheAllocationCallback
{
	virtual void allocate(long long unsigned usage) = 0;
};

template <typename Key, typename Data>
struct Item
{
	Key key;
	Data data;
	int next;
	
	Item() {}
};

// A custom hash structure.
// The actual keys and data objects are stored in a contiguous array, while
// a separate hash table maps the hashes of keys to indices into the array.
// Collisions are handled with open addressing: Each object in the contiguous
// array has a pointer to the next object in the array with an identical hash.
// NOTE: the maximum number of elements stored is 2^31 - 1
template <typename Key, typename Data>
struct Cache
{
	// a hash table that maps keys to indices to a contiguous array
	int *table;
	
	// a contiguous array that contains the actual stored objects,
	// each object also contains next pointers to handle collisions
	Item<Key, Data> *array;
	
	// number of elements in table and array
	unsigned table_size;
	unsigned array_size;
	
	// index of the current hash prime in hash_primes
	unsigned hash_prime_index;
	
	// number of objects stored
	unsigned used;
	
	CacheOptions options;
	
	CacheAllocationCallback *allocation_callback;
	
	Cache(CacheOptions options, CacheAllocationCallback *callback=NULL) :
		table(NULL),
		array(NULL),
		options(options),
		allocation_callback(callback)
	{
		if (options.verbose) eprintf("Allocating cache (%i, %i, %f, %i, %i)...\n",
			options.initial_prime_index,
			options.prime_index_increment,
			options.max_load_factor,
			options.initial_array_size,
			options.array_resize_factor
		);
		
		used = 0;
		
		hash_prime_index = options.initial_prime_index;
		table_size = hash_primes[hash_prime_index];
		array_size = options.initial_array_size;
		
		table = new int[table_size];
		array = new Item<Key, Data>[array_size];
		
		for (unsigned i = 0; i < table_size; ++i) {
			table[i] = -1;
		}
	}
	
	~Cache()
	{
		delete [] table;
		delete [] array;
	}
	
	// invokes delete for all Data objects currently stored
	// (this may only be called if Data is a pointer type)
	void delete_data()
	{
		for (unsigned i = 0; i < used; ++i) {
			delete array[i].data;
		}
	}
	
	// resizes the hash table and rehashes all objects
	void rehash(unsigned new_size)
	{
		if (allocation_callback) {
			allocation_callback->allocate(memory_usage(new_size, array_size));
		}
		
		delete [] table;
		table = NULL;
		
		table_size = new_size;
		
		table = new int[table_size];
		
		for (unsigned i = 0; i < table_size; ++i) {
			table[i] = -1;
		}
		
		for (unsigned i = 0; i < used; ++i) {
			array[i].next = -1;
		}
		
		for (unsigned add = 0; add < used; ++add) {
			unsigned h1 = std::hash<Key>()(array[add].key) % table_size;
			
			int i = table[h1];
			int prev = -1;
			
			while (i != -1) {
				prev = i;
				i = array[i].next;
			}
			
			if (prev == -1) {
				table[h1] = add;
			} else {
				array[prev].next = add;
			}
		}
	}
	
	// allocates space for a new object in the array,
	// resizing the array and rehashing all objects if necessary
	void allocate_space()
	{
		// if the array is full, resize
		if (used == array_size) {
			int new_size = array_size * options.array_resize_factor;
			if (options.verbose) eprintf("  Reallocating to... %u\n", new_size);
			
			if (allocation_callback) {
				allocation_callback->allocate(memory_usage(table_size, array_size + new_size));
			}
			
			Item<Key, Data> *new_array = new Item<Key, Data>[new_size];
			std::copy(array, array + array_size, new_array);
			delete [] array;
			array = new_array;
			array_size = new_size;
		}
		
		// if the load factor is too large, rehash
		if (used >= table_size * options.max_load_factor) {
			// get the next hash prime
			hash_prime_index += options.prime_index_increment;
			
			if (hash_prime_index > 25) hash_prime_index = 25;
			
			unsigned size = hash_primes[hash_prime_index];
			
			if (options.verbose) eprintf("  Rehashing to %u...\n", size);
			rehash(size);
		}
	}
	
	int count(Key &key) const
	{
		std::size_t h1 = std::hash<Key>()(key) % table_size;
		
		int i = table[h1];
		
		while (i != -1) {
			if (array[i].key == key) return 1;
			i = array[i].next;
		}
		
		return 0;
	}
	
	// indexing of the structure
	Data& operator [] (Key &key)
	{
		allocate_space();
		
		// compute the hash
		std::size_t h1 = std::hash<Key>()(key) % table_size;
		
		// get the index of the first matching object in the array (-1 if the hash doesn't match any object)
		int i = table[h1];
		int prev = -1;
		
		// search through all objects that match the hash
		// (each matching object has a pointer (next) to the next object with the same hash)
		while (i != -1) {
			if (array[i].key == key) return array[i].data;
			prev = i;
			i = array[i].next;
		}
		
		// if didn't find the object, adding it now (used is the first free index)
		assert(used < 2147483648);
		int add = used++;
		array[add].key = key;
		array[add].next = -1;    // last object with this hash
		
		// if this is the first object with this hash, make the table point to this object,
		// otherwise make the previous object with the same hash point to this object
		if (prev == -1) {
			table[h1] = add;
		} else {
			array[prev].next = add;
		}
		
		// return the associated data
		return array[add].data;
	}
	
	// returns the (hypothetical) memory usage on given table_size and array_size
	long long unsigned memory_usage(long long unsigned table_size_, long long unsigned array_size_)
	{
		long long unsigned table_memory = table_size_ * sizeof(int);
		long long unsigned array_memory = array_size_ * sizeof(Item<Key, Data>);
		
		return table_memory + array_memory;
	}
	
	// returns the current memory usage based on current table_size and array_size
	long long unsigned memory_usage()
	{
		return memory_usage(table_size, array_size);
	}
};

#endif
