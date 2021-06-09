#ifndef TenColoredPointBitset_header
#define TenColoredPointBitset_header

#include "KPartiteKClique/kpkc.h"
using namespace std;

class Bitset2;

class Bitset2 : public Bitset {
    public:
        Bitset2(int n_vertices) : Bitset(n_vertices) {}
        void operator=(const Bitset2&);
        void clear();
        void union_assign(Bitset2& l, Bitset2& r);
        void flip_inplace();
};

#endif
