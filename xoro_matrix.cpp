#include "xoro_matrix.h"
#include "utils.h"
//#include <intrin.h>

static const bool DEBUG_MODE = false;

// -----------------------------------------------------------

void copyFXTM(const FXTMatrix* from, FXTMatrix* to)
{
    for (int qi = 0; qi < 2; qi++)
    {
        for (int qj = 0; qj < 2; qj++)
        {
            for (int k = 0; k < 64; k++)
                (to->M)[qi][qj][k] = (from->M)[qi][qj][k];
        }
    }
}

void transposeFXTM(const FXTMatrix* matrix, FXTMatrix* transposed)
{
    for (int qi = 0; qi < 2; qi++) for (int qj = 0; qj < 2; qj++)
    {
        for (int j = 0; j < 64; j++)
        {
            const int sh = 63-j;
            uint64_t val = 0ULL;

            for (int i = 0; i < 64; i++)
            {
                val <<= 1;
                val |= ((matrix->M)[qi][qj][i] >> sh) & 1ULL;
            }
            (transposed->M)[qj][qi][j] = val;
        }
    }
}

// -----------------------------------------------------------


// bT is the FXTM "transposition" of b
void multiplyFXTM(const FXTMatrix* a, const FXTMatrix* bT, FXTMatrix* c)
{
    // calculate each bit in the result matrix separately
    for (int qi = 0; qi < 2; qi++) for (int qj = 0; qj < 2; qj++)
    {
        for (int i = 0; i < 64; i++)
        {
            uint64_t res = 0ULL;
            // these will be constant throughout the j-loop
            const uint64_t aLo = a->M[qi][0][i], aHi = a->M[qi][1][i];

            for (int j = 0; j < 64; j++)
            {
                res <<= 1;
                const uint64_t mulXor = (aLo & bT->M[qj][0][j]) ^ (aHi & bT->M[qj][1][j]);
                const int newBit = getParity(mulXor);
                res |= newBit; 
            }

            (c->M)[qi][qj][i] = res;
            // DEBUG("%llx\n", res);
        }
    }
}

void makeIdentityFXTM(FXTMatrix* fxtm) {
    copyFXTM(&IDENTITY_FXTM, fxtm);
}

void fastXoroMatrixPower(const FXTMatrix *matrix, uint64_t power, FXTMatrix *result)
{
    if (power == 0) {
        makeIdentityFXTM(result);
    }
    else {
        fastXoroMatrixPower128(matrix, 0ULL, power, result);
    }
}

// the two uint64_t fields give us the ability to advance by an arbitrary number of Xoroshiro states.
void fastXoroMatrixPower128(const FXTMatrix* matrix, uint64_t powerUpper64, uint64_t powerLower64,  FXTMatrix* result)
{
    FXTMatrix res1 = { 0 }, res2 = { 0 };
    FXTMatrix pow1 = { 0 }, pow2 = { 0 };
    FXTMatrix transposed = { 0 };
    FXTMatrix *temp;
    
    bool isResultZero = true;

    FXTMatrix *currentPower = &pow1, *nextPower = &pow2, *currentResult = &res1, *nextResult = &res2;
    copyFXTM(matrix, currentPower);

    uint64_t power = powerLower64;
    for (int bitID = 0; bitID < 128; bitID++)
    {
        if (power == 0ULL)
        {
            if (bitID > 63)
                break;  // already did upper powers
            else if (powerUpper64 == 0ULL)
                break;  // there is no upper power, just break
            // else continue until condition 3. reached, filling all the lower powers
        }

        bool transpositionDone = false;
        if (power & 1ULL)
        {
            DEBUG("if power & 1:  power = %llu\n", power);
            if (isResultZero)
            {
                copyFXTM(currentPower, nextResult);
                DEBUG("copied successfully\n");
                isResultZero = false;
            }
            else
            {
                transposeFXTM(currentPower, &transposed);
                transpositionDone = true;
                DEBUG("else: transposed successfully\n");
                multiplyFXTM(currentResult, &transposed, nextResult);
                DEBUG("else: performed fast mul successfully\n");
            }
                
            // swap next and current result, clearing space for the next operations
            temp = currentResult;
            currentResult = nextResult;
            nextResult = temp;
        }

        // create the transposition of currentPower
        DEBUG("calc next:  power = %llu\n", power)
        if (!transpositionDone)
            transposeFXTM(currentPower, &transposed);

        // perform a very fast, quasi-quadratic matrix multiplication
        multiplyFXTM(currentPower, &transposed, nextPower);
        power >>= 1;
        DEBUG("Multiplied, current bitID = %d\n", bitID)

        // swap next and current power, making space for the next operations
        temp = currentPower;
        currentPower = nextPower;
        nextPower = temp;

        if (bitID == 63)
        {
            // at key point, switch out bits to upper power
            DEBUG("Switch to upper power\n")
            power = powerUpper64;
        }
    }

    // copy result to output
    copyFXTM(currentResult, result);
}

// ---------------------------------------

void advanceXoroshiroFXTM(Xoroshiro *state, const FXTMatrix* fxtm)
{
    uint64_t oldState[2] = { state->lo, state->hi };
    uint64_t newState[2] = { 0ULL, 0ULL };

    for (int qi = 0; qi < 2; qi++) for (int qj = 0; qj < 2; qj++)
    {
        for (int i = 0; i < 64; i++)
        {
            const int sh = 63 - i;

            for (int j = 0; j < 64; j++)
            {
                const uint64_t bit1 = ((fxtm->M)[qi][qj][i] >> j) & 1ULL;
                const uint64_t bit2 = (oldState[qi] >> sh) & 1ULL;
                newState[qj] ^= (bit1 & bit2) << j;
            }
        }
    }

    state->lo = newState[0];
    state->hi = newState[1];
}

void fastAdvanceXoroshiroFXTM(Xoroshiro *state, const FXTMatrix* fxtmT)
{
    uint64_t oldState[2] = { state->lo, state->hi };
    uint64_t newState[2] = { 0ULL, 0ULL };

    for (int qi = 0; qi < 2; qi++) {
        for (int i = 0; i < 64; i++) {
            const uint64_t mulXor = (oldState[0] & fxtmT->M[qi][0][63-i]) ^ (oldState[1] & fxtmT->M[qi][1][63-i]);
            newState[qi] |= static_cast<uint64_t>(getParity(mulXor)) << i;
        }
    }

    state->lo = newState[0];
    state->hi = newState[1];
}

static void xorRows(FXTMatrix *left, FXTMatrix *right, int src_q, int src_i, int dst_q, int dst_i)
{
    for (int qj = 0; qj < 2; qj++) {
        left->M[dst_q][qj][dst_i] ^= left->M[src_q][qj][src_i];
        right->M[dst_q][qj][dst_i] ^= right->M[src_q][qj][src_i];
    }
}

void invertFXTM(const FXTMatrix *fxtm, FXTMatrix *inverted)
{
    FXTMatrix left, right;
    copyFXTM(fxtm, &left);
    makeIdentityFXTM(&right);

    for (int diag = 0; diag < 128; diag++) {
        int p_q = diag / 64;
        int p_i = diag % 64;
        int p_bit = 63 - p_i;

        int val = (left.M[p_q][p_q][p_i] >> p_bit) & 1;
        if (val == 0) {
            // Find another row below the row currently being fixed that has a 1 in the target col
            for (int r = diag + 1; r < 128; r++) {
                int r_q = r / 64;
                int r_i = r % 64;
                int r_val = (left.M[r_q][p_q][r_i] >> p_bit) & 1;
                if (r_val == 1) {
                    xorRows(&left, &right, r_q, r_i, p_q, p_i);
                    break;
                }
            }
        }

        for (int row = 0; row < 128; row++) {
            if (row == diag) continue;
            int r_q = row / 64;
            int r_i = row % 64;
            int r_val = (left.M[r_q][p_q][r_i] >> p_bit) & 1;
            if (r_val == 1) {
                xorRows(&left, &right, p_q, p_i, r_q, r_i);
            }
        }
    }

    copyFXTM(&right, inverted);
}

void printXoro(const Xoroshiro* state)
{
    printf("Current Xoroshiro state:   %016llx-%016llx\n", state->lo, state->hi);
}

/*
int main() {
    Xoroshiro a, b;
    a.lo = b.lo = 0;
    a.hi = b.hi = 1;

    FXTMatrix xInv, xInvInv;
    invertFXTM(&XOROSHIRO_STANDARD_FXTM, &xInv);
    invertFXTM(&xInv, &xInvInv);

    advanceXoroshiroFXTM(&a, &XOROSHIRO_STANDARD_FXTM);
    advanceXoroshiroFXTM(&b, &xInvInv);

    printXoro(&a);
    printXoro(&b);
}
*/

// -------------------------------------------------
// Base matrix defs

const FXTMatrix XOROSHIRO_STANDARD_FXTM = {{
    {
        {
        9223653511831486464ULL,
        4611826755915743232ULL,
        2305913377957871616ULL,
        1152956688978935808ULL,
        576478344489467904ULL,
        288239172244733952ULL,
        144119586122366976ULL,
        72059793061183488ULL,
        36029896530591744ULL,
        18014948265295872ULL,
        9007474132647936ULL,
        4503737066323968ULL,
        2251868533161984ULL,
        1125934266580992ULL,
        562967133290496ULL,
        281483566645248ULL,
        140741783322624ULL,
        70370891661312ULL,
        35185445830656ULL,
        17592722915328ULL,
        8796361457664ULL,
        9223376435035504640ULL,
        4611688217517752320ULL,
        2305844108758876160ULL,
        1152922054379438080ULL,
        576461027189719040ULL,
        288230513594859520ULL,
        144115256797429760ULL,
        72057628398714880ULL,
        36028814199357440ULL,
        18014407099678720ULL,
        9007203549839360ULL,
        4503601774919680ULL,
        2251800887459840ULL,
        1125900443729920ULL,
        562950221864960ULL,
        281475110932480ULL,
        140737555466240ULL,
        70368777733120ULL,
        35184388866560ULL,
        17592194433280ULL,
        8796097216640ULL,
        4398048608320ULL,
        2199024304160ULL,
        1099512152080ULL,
        549756076040ULL,
        274878038020ULL,
        137439019010ULL,
        68719509505ULL,
        9223372071214530560ULL,
        4611686035607265280ULL,
        2305843017803632640ULL,
        1152921508901816320ULL,
        576460754450908160ULL,
        288230377225454080ULL,
        144115188612727040ULL,
        72057594306363520ULL,
        36028797153181760ULL,
        18014398576590880ULL,
        9007199288295440ULL,
        4503599644147720ULL,
        2251799822073860ULL,
        1125899911036930ULL,
        562949955518465ULL
        },
        {
        134217728ULL,
        67108864ULL,
        33554432ULL,
        16777216ULL,
        8388608ULL,
        4194304ULL,
        2097152ULL,
        1048576ULL,
        524288ULL,
        262144ULL,
        131072ULL,
        65536ULL,
        32768ULL,
        16384ULL,
        8192ULL,
        4096ULL,
        2048ULL,
        1024ULL,
        512ULL,
        256ULL,
        128ULL,
        64ULL,
        32ULL,
        16ULL,
        8ULL,
        4ULL,
        2ULL,
        1ULL,
        9223372036854775808ULL,
        4611686018427387904ULL,
        2305843009213693952ULL,
        1152921504606846976ULL,
        576460752303423488ULL,
        288230376151711744ULL,
        144115188075855872ULL,
        72057594037927936ULL,
        36028797018963968ULL,
        18014398509481984ULL,
        9007199254740992ULL,
        4503599627370496ULL,
        2251799813685248ULL,
        1125899906842624ULL,
        562949953421312ULL,
        281474976710656ULL,
        140737488355328ULL,
        70368744177664ULL,
        35184372088832ULL,
        17592186044416ULL,
        8796093022208ULL,
        4398046511104ULL,
        2199023255552ULL,
        1099511627776ULL,
        549755813888ULL,
        274877906944ULL,
        137438953472ULL,
        68719476736ULL,
        34359738368ULL,
        17179869184ULL,
        8589934592ULL,
        4294967296ULL,
        2147483648ULL,
        1073741824ULL,
        536870912ULL,
        268435456ULL
        }
    },
    {
        {
        9223372036854775808ULL,
        4611686018427387904ULL,
        2305843009213693952ULL,
        1152921504606846976ULL,
        576460752303423488ULL,
        288230376151711744ULL,
        144115188075855872ULL,
        72057594037927936ULL,
        36028797018963968ULL,
        18014398509481984ULL,
        9007199254740992ULL,
        4503599627370496ULL,
        2251799813685248ULL,
        1125899906842624ULL,
        562949953421312ULL,
        281474976710656ULL,
        140737488355328ULL,
        70368744177664ULL,
        35184372088832ULL,
        17592186044416ULL,
        8796093022208ULL,
        9223376434901286912ULL,
        4611688217450643456ULL,
        2305844108725321728ULL,
        1152922054362660864ULL,
        576461027181330432ULL,
        288230513590665216ULL,
        144115256795332608ULL,
        72057628397666304ULL,
        36028814198833152ULL,
        18014407099416576ULL,
        9007203549708288ULL,
        4503601774854144ULL,
        2251800887427072ULL,
        1125900443713536ULL,
        562950221856768ULL,
        281475110928384ULL,
        140737555464192ULL,
        70368777732096ULL,
        35184388866048ULL,
        17592194433024ULL,
        8796097216512ULL,
        4398048608256ULL,
        2199024304128ULL,
        1099512152064ULL,
        549756076032ULL,
        274878038016ULL,
        137439019008ULL,
        68719509504ULL,
        34359754752ULL,
        17179877376ULL,
        8589938688ULL,
        4294969344ULL,
        2147484672ULL,
        1073742336ULL,
        536871168ULL,
        268435584ULL,
        134217792ULL,
        67108896ULL,
        33554448ULL,
        16777224ULL,
        8388612ULL,
        4194306ULL,
        2097153ULL
        },
        {
        134217728ULL,
        67108864ULL,
        33554432ULL,
        16777216ULL,
        8388608ULL,
        4194304ULL,
        2097152ULL,
        1048576ULL,
        524288ULL,
        262144ULL,
        131072ULL,
        65536ULL,
        32768ULL,
        16384ULL,
        8192ULL,
        4096ULL,
        2048ULL,
        1024ULL,
        512ULL,
        256ULL,
        128ULL,
        64ULL,
        32ULL,
        16ULL,
        8ULL,
        4ULL,
        2ULL,
        1ULL,
        9223372036854775808ULL,
        4611686018427387904ULL,
        2305843009213693952ULL,
        1152921504606846976ULL,
        576460752303423488ULL,
        288230376151711744ULL,
        144115188075855872ULL,
        72057594037927936ULL,
        36028797018963968ULL,
        18014398509481984ULL,
        9007199254740992ULL,
        4503599627370496ULL,
        2251799813685248ULL,
        1125899906842624ULL,
        562949953421312ULL,
        281474976710656ULL,
        140737488355328ULL,
        70368744177664ULL,
        35184372088832ULL,
        17592186044416ULL,
        8796093022208ULL,
        4398046511104ULL,
        2199023255552ULL,
        1099511627776ULL,
        549755813888ULL,
        274877906944ULL,
        137438953472ULL,
        68719476736ULL,
        34359738368ULL,
        17179869184ULL,
        8589934592ULL,
        4294967296ULL,
        2147483648ULL,
        1073741824ULL,
        536870912ULL,
        268435456ULL
        }
    }
}};

const FXTMatrix IDENTITY_FXTM = {{
    {
        {
        9223372036854775808ULL,
        4611686018427387904ULL,
        2305843009213693952ULL,
        1152921504606846976ULL,
        576460752303423488ULL,
        288230376151711744ULL,
        144115188075855872ULL,
        72057594037927936ULL,
        36028797018963968ULL,
        18014398509481984ULL,
        9007199254740992ULL,
        4503599627370496ULL,
        2251799813685248ULL,
        1125899906842624ULL,
        562949953421312ULL,
        281474976710656ULL,
        140737488355328ULL,
        70368744177664ULL,
        35184372088832ULL,
        17592186044416ULL,
        8796093022208ULL,
        4398046511104ULL,
        2199023255552ULL,
        1099511627776ULL,
        549755813888ULL,
        274877906944ULL,
        137438953472ULL,
        68719476736ULL,
        34359738368ULL,
        17179869184ULL,
        8589934592ULL,
        4294967296ULL,
        2147483648ULL,
        1073741824ULL,
        536870912ULL,
        268435456ULL,
        134217728ULL,
        67108864ULL,
        33554432ULL,
        16777216ULL,
        8388608ULL,
        4194304ULL,
        2097152ULL,
        1048576ULL,
        524288ULL,
        262144ULL,
        131072ULL,
        65536ULL,
        32768ULL,
        16384ULL,
        8192ULL,
        4096ULL,
        2048ULL,
        1024ULL,
        512ULL,
        256ULL,
        128ULL,
        64ULL,
        32ULL,
        16ULL,
        8ULL,
        4ULL,
        2ULL,
        1ULL
        },
        {
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL
        }
    },
    {
        {
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL,
        0ULL 
        },
        {
        9223372036854775808ULL,
        4611686018427387904ULL,
        2305843009213693952ULL,
        1152921504606846976ULL,
        576460752303423488ULL,
        288230376151711744ULL,
        144115188075855872ULL,
        72057594037927936ULL,
        36028797018963968ULL,
        18014398509481984ULL,
        9007199254740992ULL,
        4503599627370496ULL,
        2251799813685248ULL,
        1125899906842624ULL,
        562949953421312ULL,
        281474976710656ULL,
        140737488355328ULL,
        70368744177664ULL,
        35184372088832ULL,
        17592186044416ULL,
        8796093022208ULL,
        4398046511104ULL,
        2199023255552ULL,
        1099511627776ULL,
        549755813888ULL,
        274877906944ULL,
        137438953472ULL,
        68719476736ULL,
        34359738368ULL,
        17179869184ULL,
        8589934592ULL,
        4294967296ULL,
        2147483648ULL,
        1073741824ULL,
        536870912ULL,
        268435456ULL,
        134217728ULL,
        67108864ULL,
        33554432ULL,
        16777216ULL,
        8388608ULL,
        4194304ULL,
        2097152ULL,
        1048576ULL,
        524288ULL,
        262144ULL,
        131072ULL,
        65536ULL,
        32768ULL,
        16384ULL,
        8192ULL,
        4096ULL,
        2048ULL,
        1024ULL,
        512ULL,
        256ULL,
        128ULL,
        64ULL,
        32ULL,
        16ULL,
        8ULL,
        4ULL,
        2ULL,
        1ULL
        }
    }
}};
