#ifndef TenColoredPointBitset_header
#define TenColoredPointBitset_header

#include "KPartiteKClique/kpkc.h"

namespace tencoloredpoints
{
class Bitset2;

class Bitset2 : public kpkc::Bitset {
    public:
        Bitset2(int n_vertices) : Bitset(n_vertices) {}
        void operator=(const Bitset2&);
        void clear();
        void union_assign(const Bitset2& l, const Bitset2& r);
        void flip_inplace();
};
}

#endif
