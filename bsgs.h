#ifndef BSGS_H
#define BSGS_H

#include "xoro_matrix.h"

uint64_t solve_discrete_log_fast(FXTMatrix *A, FXTMatrix *M, FXTMatrix *Minv, uint64_t n);
void project_into_subgroup(const FXTMatrix *A, const FXTMatrix *M, FXTMatrix *Ad, FXTMatrix *Md, FXTMatrix *Mdinv, uint64_t f);

#endif