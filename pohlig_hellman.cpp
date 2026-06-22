#include <cstdint>
#include <cstdio>
#include <iostream>

#include "pohlig_hellman.h"
#include "bsgs.h"
#include "xoro_matrix.h"

// factors of 2^128 - 1
static InfInt factors[9] = {
    InfInt(3),
    InfInt(5),
    InfInt(17),
    InfInt(257),
    InfInt(641),
    InfInt(65537),
    InfInt(274177),
    InfInt(6700417),
    InfInt(67280421310721),
};

// Adapted from the pseudocode in this Wikipedia article: https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
static InfInt extended_gcd_128(InfInt a, InfInt b, InfInt& s_out, InfInt& t_out) {
    InfInt old_r = a, r = b;
    InfInt old_s = InfInt(1), s = InfInt(0);
    InfInt old_t = InfInt(0), t = InfInt(1);
    
    while (r != 0) {
        InfInt quotient = old_r / r;
        InfInt new_r = old_r - quotient * r;
        InfInt new_s = old_s - quotient * s;
        InfInt new_t = old_t - quotient * t;
        old_r = r;
        old_s = s;
        old_t = t;
        r = new_r;
        s = new_s;
        t = new_t;
    }

    s_out = old_s;
    t_out = old_t;
    return old_r;
}

// Uses CRT to reconstruct the full discrete log value
static InfInt get_value_from_remainders(InfInt moduli[], InfInt remainders[], uint32_t num_eqns) {
    for (int32_t i = num_eqns - 2; i >= 0; i--) {
        InfInt m1, m2;
        InfInt n1 = moduli[i], n2 = moduli[i+1];
        InfInt a1 = remainders[i], a2 = remainders[i+1];
        (void)extended_gcd_128(n1, n2, m1, m2);
    
        InfInt mod_product = n1 * n2;
        InfInt val = (a1*m2*n2 + a2*m1*n1) % mod_product;
        if (val < 0) {
            val += mod_product;
        }

        remainders[i] = val;
        moduli[i] = n1 * n2;
    }
    return remainders[0];
}

void advance_infint(Xoroshiro *state, const InfInt& jump) {
    InfInt two64 = InfInt(1LL << 32) * InfInt(1LL << 32);
    uint64_t lo = (jump % two64).toUnsignedLongLong();
    uint64_t hi = (jump / two64).toUnsignedLongLong();

    FXTMatrix fxtm, fxtmT;
    fastXoroMatrixPower128(&XOROSHIRO_STANDARD_FXTM, hi, lo, &fxtm);
    transposeFXTM(&fxtm, &fxtmT);
    fastAdvanceXoroshiroFXTM(state, &fxtmT);
}

static InfInt xoroshiro_distance(Xoroshiro a, Xoroshiro b, const FXTMatrix *base_matrix, bool log_progress) {
    FXTMatrix Ks, Kt, Ksinv, A;

    /*
    To perform baby-step-giant-step, we recover the underlying matrix that advances the Xoroshiro state 
    from a to b: A=M^n by stacking multiple state transforms s*M^i --A--> t*M^i (for i=0..127) 
    into a single big operation S * A = T. Having that equation with known S and T, it's easy to see how 
    A = S^-1 * T. Once we have the matrix A, the problem becomes a standard discrete logarithm in a matrix group.
    */
    Xoroshiro curr_s = a;
    Xoroshiro curr_t = b;

    FXTMatrix MT;
    transposeFXTM(base_matrix, &MT);
    
    for (int row = 0; row < 128; row++) {
        int qj = row / 64;
        int j = row % 64;
        Ks.M[qj][0][j] = curr_s.lo;
        Ks.M[qj][1][j] = curr_s.hi;
        Kt.M[qj][0][j] = curr_t.lo;
        Kt.M[qj][1][j] = curr_t.hi;
        
        fastAdvanceXoroshiroFXTM(&curr_s, &MT);
        fastAdvanceXoroshiroFXTM(&curr_t, &MT);
    }

    // Reconstruct A = M^n
    invertFXTM(&Ks, &Ksinv);

    FXTMatrix KtT;
    transposeFXTM(&Kt, &KtT);
    multiplyFXTM(&Ksinv, &KtT, &A);

    //printf("Target A: %016llx-%016llx\n", XOROSHIRO_STANDARD_FXTM.M[0][0][0], XOROSHIRO_STANDARD_FXTM.M[1][0][0]);
    //printf("Calced A: %016llx-%016llx\n", A.M[0][0][0], A.M[1][0][0]);
    //return;

    InfInt remainders[9];
    InfInt work_factors[9];
    for (int i = 0; i < 9; i++) {
        work_factors[i] = factors[i];
    }

    for (int i = 0; i < 9; i++) {
        uint64_t n = work_factors[i].toUnsignedLongLong();
        FXTMatrix Ad, Md, Mdinv;
    
        project_into_subgroup(&A, base_matrix, &Ad, &Md, &Mdinv, n);
        uint64_t dl = solve_discrete_log_fast(&Ad, &Md, &Mdinv, n);

        remainders[i] = InfInt(dl);
        if (log_progress) {
            fprintf(stderr, "discrete log for subgroup N=%llu : %llu\n", n, dl);
        }
    }
    InfInt distance = get_value_from_remainders(work_factors, remainders, 9);
    return distance;
}

InfInt get_xoroshiro_state_distance(Xoroshiro a, Xoroshiro b, bool log_progress) {
    return xoroshiro_distance(a, b, &XOROSHIRO_STANDARD_FXTM, log_progress);
}

InfInt get_xoroshiro_jump_with_column(int matrix_column_index, Xoroshiro matrix_column, bool log_progress) {
    /*
    This uses pretty much the same discrete log process as xoroshiro state distance,
    as we can translate finding a set jump matrix column into the same kind of problem:
    
    C = M^d * U_n
    where U_n is a column vector with all zeros and a 1 at position n.
    Transposing both sides, we obtain

    C^T = (U_n)^T * (M^T)^d

    Now, both C^T and (U_n)^T are row vectors, so we essentially need the same calculations
    as we did for xoroshiro state distance. The only difference is that we're using the transposed
    Xoroshiro state advancement matrix instead of the regular one.
    */
    
    FXTMatrix xT;
    transposeFXTM(&XOROSHIRO_STANDARD_FXTM, &xT);

    bool is_low = matrix_column_index < 64;
    int bit_idx = matrix_column_index % 64;
    Xoroshiro a = {0, 0};
    if (is_low) {
        a.lo = 1ULL << bit_idx;
    }
    else {
        a.hi = 1ULL << bit_idx;
    }

    return xoroshiro_distance(matrix_column, a, &xT, log_progress);
}