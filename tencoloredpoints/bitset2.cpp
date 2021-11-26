#include "bitset2.h"

#include <stdlib.h>
#include <cstring>
#include <cstdlib>

// Bitsets

namespace tencoloredpoints {

void Bitset2::clear(){
    for(int i=0; i<limbs; i++)
        data[i] = 0;
}

void Bitset2::operator=(const Bitset2& obj){
    for(int i=0; i<limbs; i++)
        data[i] = obj.data[i];
}

inline void Bitset2::flip_inplace(){
    for(int i=0; i<limbs; i++)
        data[i] = ~data[i];
}

inline void Bitset2::union_assign(const Bitset2& l, const Bitset2& r){
    for(int i = 0; i < limbs; i++)
        data[i] = l[i] | r[i];
}

}
