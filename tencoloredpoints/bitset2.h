#ifndef TenColoredPointBitset_header
#define TenColoredPointBitset_header

#include "KPartiteKClique/kpkc.h"
using namespace std;

class Bitset2;

class Bitset2 : public Bitset {
    public:
        Bitset2(int n_vertices);
        ~Bitset2();
        void init();
        void init(Bitset2& obj);
        void union_assign(Bitset2& l, Bitset2& r);
        void flip_inplace();
    private:
        uint64_t* mem;
};

#endif
