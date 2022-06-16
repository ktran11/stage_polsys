/**
 * \file mat_pol.h
 * \brief shorts functions
 * \author Kevin Tran
 * \version 1.0
 * \date June 6 2022
 *
 * Add some functions availables on SageMath for the class nmod_poly_mat_t on flint
 * https://doc.sagemath.org/html/en/reference/matrices/sage/matrix/matrix_polynomial_dense.html
 *
 */

#ifndef M_BASIS_H
#define M_BASIS_H

#include "matpol.h"
#include "basis.h"
#include "perm_operator.h"
#include "nmod_list_poly_mat.h"

#include <stdint.h>
#include <nmod_mat.h>
#include <nmod_poly_mat.h>

void M_basis(nmod_poly_mat_t res, int64_t *res_shifts,
	     const nmod_poly_mat_t F, uint64_t sigma, const int64_t *shifts);
  
void M_basisII(nmod_poly_mat_t res, int64_t *res_shifts,
	       const nmod_poly_mat_t F, uint64_t sigma, const int64_t *shifts);

#endif /* M_BASIS_H */

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
