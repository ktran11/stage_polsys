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

#ifndef PERM_OPERATOR_H
#define PERM_OPERATOR_H

#include <stdint.h>
#include <nmod_mat.h>
#include <nmod_poly_mat.h>
#include <perm.h> 

/**
 * \brief Do the multiplication permutation matrix (P*A) for the function Basis
 *
 * \param A, the input and result matrix
 * \param P, the permutation 
 * \param n, length of P and row dimension of A
 */
void apply_perm_rows_to_matrix(nmod_mat_t mat, const slong * perm, slong rdim);

void apply_perm_cols_to_poly_matrix(nmod_poly_mat_t mat, const slong * perm, slong cdim);

void apply_perm_rows_to_poly_matrix(nmod_poly_mat_t mat, const slong *perm, slong rdim);

void apply_perm_to_vector(int64_t *res, const int64_t *initial_vect,
			  const slong *perm, slong length);

void sort_and_create_perm(slong *perm, const int64_t *vec, slong n);

#endif /* PERM_OPERATOR_H */

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
