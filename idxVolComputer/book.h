// From the software distribution accompanying the textbook
// "A Practical Introduction to Data Structures and Algorithm Analysis,
// Third Edition (C++)" by Clifford A. Shaffer.
// Source code Copyright (C) 2007-2011 by Clifford A. Shaffer.

// A collection of various macros, constants, and small functions
// used for the software examples.
#ifndef BOOK_H
#define BOOK_H
// First, include all the standard headers that we need
#include <iostream>
#include <cstdlib>
#include <time.h>  // Used by timing functions
#include <vector>
#include <algorithm>
// Now all the standard names that we use
using std::cout;
using std::endl;
using std::string;
using std::ostream;
using std::vector;

const int defaultSize = 10; // Default size

// Return true iff "x" is even
inline bool EVEN(int x) { return (x % 2) == 0; }

// Return true iff "x" is odd
inline bool ODD(int x) { return (x % 2) != 0; }


// Swap two elements in a generic array
template<typename E>
inline void swap(E A[], int i, int j) {
	E temp = A[i];
	A[i] = A[j];
	A[j] = temp;
}
// Random number generator functions

inline void Randomize() // Seed the generator
{
	srand(1);
}

// Return a random value in range 0 to n-1
inline int Random(int n)
{
	return rand() % (n);
}


// Swap two integers
inline void swap(int& i, int& j) {
	int temp = i;
	i = j;
	j = temp;
}

// Swap two char*'s
inline void swap(char* i, char* j) {
	char* temp = i;
	i = j;
	j = temp;
}

// Your basic int type as an object.
class Int {
private:
	int val;
public:
	Int(int input = 0) { val = input; }
	// The following is for those times when we actually
	//   need to get a value, rather than compare objects.
	int key() const { return val; }
	// Overload = to support Int foo = 5 syntax
	Int operator= (int input) { val = input; return val; }
	// Let us print out Ints easily
	friend ostream& operator<<(ostream& s, const Int& i)
	{
		return s << i.key();
	}
	friend ostream& operator<<(ostream& s, const Int* i)
	{
		return s << i->key();
	}

};


static int s_discrimDim = 0;
template<typename E>
bool compareTD(const E& a1, const E& a2)
{
	return a1->coord()[s_discrimDim] < a2->coord()[s_discrimDim];
}

template<typename E>
E median(vector<E>& v, size_t& medianId, int discrim)
{
	size_t n = v.size() ;
	if (n == 1){
		medianId = 0;
		return v[medianId];
	}
		
	s_discrimDim = discrim;
	sort(v.begin(), v.end(), compareTD<E>);
	if (n % 2 == 0)
		medianId = n / 2 - 1;
	else
		medianId = n / 2;
	return v[medianId];
}

#endif