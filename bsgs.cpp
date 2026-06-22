#include <cstdint>
#include <cstdio>
#include <unordered_map>
#include <cmath>

#include "bsgs.h"
#include "infint.h"

// uint64_t solve_discrete_log(FXTMatrix *A, FXTMatrix *M, FXTMatrix *Minv, uint64_t n)
// {
//     const uint64_t m = std::ceil(std::sqrt(n));
//
//     // Solving M^k = A for k.
//
//     FXTMatrix MT, Minv_m, Minv_mT;
//     transposeFXTM(M, &MT);
//     fastXoroMatrixPower(Minv, m, &Minv_m);
//     transposeFXTM(&Minv_m, &Minv_mT);
//
//     std::unordered_map<Xoroshiro, uint64_t, XoroshiroHash> baby_step;
//
//     //printf("baby step states:\n");
//     {
//         FXTMatrix mx[2];
//         int curr = 0;
//         makeIdentityFXTM(&(mx[0]));
//         for (uint64_t j = 0; j < m; j++) {
//             Xoroshiro x = {mx[curr].M[0][0][0], mx[curr].M[1][0][0]};
//             //printf("%016llx-%016llx\n", x.hi, x.lo);
//             baby_step[x] = j;
//             multiplyFXTM(&(mx[curr]), &MT, &(mx[curr ^ 1]));
//             curr ^= 1;
//         }
//     }
//
//     //printf("giant step states:\n");
//     {
//         FXTMatrix mx[2];
//         int curr = 0;
//         copyFXTM(A, &(mx[0]));
//         for (uint64_t i = 0; i < m; i++) {
//             Xoroshiro x = {mx[curr].M[0][0][0], mx[curr].M[1][0][0]};
//             //printf("%016llx-%016llx\n", x.hi, x.lo);
//             if (baby_step.find(x) != baby_step.end()) {
//                 return i * m + baby_step[x];
//             }
//             multiplyFXTM(&(mx[curr]), &Minv_mT, &(mx[curr ^ 1]));
//             curr ^= 1;
//         }
//     }
//
//     return 0;
// }

uint64_t solve_discrete_log_fast(FXTMatrix *A, FXTMatrix *M, FXTMatrix *Minv, uint64_t n)
{
    const uint64_t m = std::ceil(std::sqrt(n));

    // Solving M^k = A for k.

    FXTMatrix MT, Minv_m, Minv_mT;
    transposeFXTM(M, &MT);
    fastXoroMatrixPower(Minv, m, &Minv_m);
    transposeFXTM(&Minv_m, &Minv_mT);

    std::unordered_map<Xoroshiro, uint64_t, XoroshiroHash> baby_step;

    //printf("baby step states:\n");
    {
        Xoroshiro x = {12345, 54321};
        for (uint64_t j = 0; j < m; j++) {
            baby_step[x] = j;
            fastAdvanceXoroshiroFXTM(&x, &MT);
        }
    }

    //printf("giant step states:\n");
    {
        FXTMatrix AT;
        transposeFXTM(A, &AT);

        Xoroshiro x = {12345, 54321};
        fastAdvanceXoroshiroFXTM(&x, &AT);
        for (uint64_t i = 0; i < m; i++) {
            if (baby_step.find(x) != baby_step.end()) {
                return i * m + baby_step[x];
            }
            fastAdvanceXoroshiroFXTM(&x, &Minv_mT);
        }
    }

    return 0;
}

void project_into_subgroup(const FXTMatrix *A, const FXTMatrix *M, FXTMatrix *Ad, FXTMatrix *Md, FXTMatrix *Mdinv, uint64_t f) {
    InfInt two64 = InfInt(1LL << 32) * InfInt(1LL << 32);

    InfInt bigF = InfInt(f);
    InfInt bigN = two64 * two64 - InfInt(1);
    InfInt bigD = bigN / bigF;

    uint64_t lo = static_cast<uint64_t>((bigD % two64).toUnsignedLongLong());
    uint64_t hi = static_cast<uint64_t>((bigD / two64).toUnsignedLongLong());
    //printf("Projection: %016llx-%016llx\n", hi, lo);

    // project element
    fastXoroMatrixPower128(A, hi, lo, Ad); 
    // project state jump matrix & calculate its inverse
    fastXoroMatrixPower128(M, hi, lo, Md);
    invertFXTM(Md, Mdinv);
}
