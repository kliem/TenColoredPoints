#if __AVX2__
    #include <immintrin.h>
    #include <emmintrin.h>
#elif __SSE4_1__
    #include <emmintrin.h>
    #include <smmintrin.h>
#endif

#include "bitset2.h"

#include <stdlib.h>
#include <cstring>
#include <cstdlib>

// Bitset helpers.

#if __AVX2__
    const uintptr_t ALIGNMENT2 = 32;
#elif __SSE4_1__
    const uintptr_t ALIGNMENT2 = 16;
#else
    const uintptr_t ALIGNMENT2 = 8;
#endif

// Bitsets

Bitset2::Bitset2(int n_vertices) : Bitset(n_vertices){
    delete[] data;  // allocating and immediatly deallocating will be optimized out
    limbs = ((n_vertices-1)/(8*ALIGNMENT2) + 1)*(ALIGNMENT2/8);
    mem = new uint64_t[limbs + 4];
    uintptr_t x = (uintptr_t) mem;
    data = (uint64_t*) ((x + ALIGNMENT2 - 1) & ~(ALIGNMENT2 - 1));
}

Bitset2::~Bitset2(){
    delete[] mem;
    data = NULL;
}

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
    // Assumes all of same length.
    int i = 0;
#if __AVX2__
    for(i; i < limbs; i += 4){
        __m256i A = _mm256_load_si256((const __m256i*)&l[i]);
        __m256i B = _mm256_load_si256((const __m256i*)&r[i]);
        __m256i D = _mm256_or_si256(A, B);
        _mm256_store_si256((__m256i*)&data[i], D);
    }
#elif __SSE4_1__
    for(i; i < limbs; i += 2){
        __m128i A = _mm_loadu_si128((const __m128i*)&l[i]);
        __m128i B = _mm_loadu_si128((const __m128i*)&r[i]);
        __m256i D = _mm_or_si128(A, B);
        _mm_storeu_si128((__m128i*)&data[i], D);
    }
#else
    for(i; i < limbs; i++){
        data[i] = l[i] | r[i];
    }
#endif
}
