#include "bitset2.h"

#include <stdlib.h>
#include <cstring>
#include <cstdlib>

// Bitsets

void Bitset2::clear(){
    memset(data, 0, limbs*8);
}

void Bitset2::operator=(const Bitset2& obj){
    memcpy(data, obj.data, limbs*8);
}

inline void Bitset2::flip_inplace(){
    for(int i=0; i<limbs; i++)
        data[i] = ~data[i];
}

inline void Bitset2::union_assign(Bitset2& l, Bitset2& r){
    for(int i = 0; i < limbs; i++)
        data[i] = l[i] | r[i];
}
