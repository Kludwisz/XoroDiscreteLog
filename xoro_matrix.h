#ifndef XORO_MATRIX_H
#define XORO_MATRIX_H

#include "xoroshiro.h"
#include "inttypes.h"

#define ONLY_TOP_BIT    (0x8000000000000000ULL)
#define FULL_64         (0xffffffffffffffffULL)

#ifdef __GNUC__
#include <x86intrin.h>
static inline int getParity(uint64_t val) {
    return __builtin_parityll(val); // aka: GCC, do your thing!
}
#else
#include <intrin.h>
static inline int getParity(uint64_t val) {
    return __popcnt64(val) & 1;
}
#endif

// structures for the fast xoroshiro transform
// matrix multiplication procedure

typedef struct FXTMatrix FXTMatrix;
struct FXTMatrix {
    uint64_t M[2][2][64];
};

// ------------------------------------------------------------------------------------

void copyFXTM(const FXTMatrix* from, FXTMatrix* to);
void transposeFXTM(const FXTMatrix *matrix, FXTMatrix *transposed);
void multiplyFXTM(const FXTMatrix *a, const FXTMatrix *bT, FXTMatrix *c);
void fastXoroMatrixPower(const FXTMatrix *matrix, uint64_t power, FXTMatrix *result);
void fastXoroMatrixPower128(const FXTMatrix *matrix, uint64_t powerUpper64, uint64_t powerLower64, FXTMatrix *result);
void advanceXoroshiroFXTM(Xoroshiro *state, const FXTMatrix *fxtm);
void fastAdvanceXoroshiroFXTM(Xoroshiro *state, const FXTMatrix* fxtmT);
void makeIdentityFXTM(FXTMatrix* fxtm);
void invertFXTM(const FXTMatrix *fxtm, FXTMatrix *inverted);

void printXoro(const Xoroshiro* state);

// -------------------------------------------------------------------------------------

extern const FXTMatrix XOROSHIRO_STANDARD_FXTM;
extern const FXTMatrix IDENTITY_FXTM;

#endif