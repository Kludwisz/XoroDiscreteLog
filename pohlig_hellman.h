#ifndef POHLIG_HELLMAN_H
#define POHLIG_HELLMAN_H

#include "infint.h"
#include "xoro_matrix.h"

/**
 * Returns d such that a * M^d = b, where M is the Xoroshiro128++ state advancement matrix.
 */
InfInt get_xoroshiro_state_distance(Xoroshiro a, Xoroshiro b, bool log_progress = false);

/**
 * Returns d such that the matrix M^d has the given matrix column at the given index.
 */
InfInt get_xoroshiro_jump_with_column(int matrix_column_index, Xoroshiro matrix_column, bool log_progress = false);

/**
 * Advances the Xoroshiro128++ state by 'jump' steps.
 */
void advance_infint(Xoroshiro *state, const InfInt& jump);

#endif