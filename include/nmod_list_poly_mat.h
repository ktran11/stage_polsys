/**
 * \file mat_pol.h
 * \brief shorts functions
 * \author Kevin Tran
 * \version 1.0
 * \date June 6 2022
 *
 * Add some functions availables on SageMath for the class nmod_poly_mat_t on flint
 * https://doc.sagemath.org/html/en/reference/matrices/sage/matrix/matrix_polynomial_dense.html
 * All parameters are supposed init
 *
 */

#ifndef NMOD_LIST_POLY_MAT_H
#define NMOD_LIST_POLY_MAT_H

#include "matpol.h"

#include <stdint.h>
#include <nmod_mat.h>
#include <nmod_poly_mat.h>
#include <perm.h>
#include "perm_operator.h"

typedef struct
{
  slong degree;
  slong length;
  nmod_mat_struct *mat;
  slong rdim;
  slong cdim;
  mp_limb_t modulus;
} nmod_list_poly_mat_struct;

typedef nmod_list_poly_mat_struct nmod_list_poly_mat_t[1];

slong nmod_list_poly_mat_nrows(const nmod_list_poly_mat_t A);

slong nmod_list_poly_mat_ncols(const nmod_list_poly_mat_t A);

slong nmod_list_poly_mat_degree(const nmod_list_poly_mat_t A);

mp_limb_t nmod_list_poly_mat_modulus(const nmod_list_poly_mat_t A);

void nmod_list_poly_mat_print(const nmod_list_poly_mat_t A);

void nmod_list_poly_mat_get_coef(nmod_mat_t res, const nmod_list_poly_mat_t F, slong k);


void nmod_list_poly_mat_init(nmod_list_poly_mat_t res, slong degree, slong length,
			     slong rdim, slong cdim, mp_limb_t modulus);

void nmod_list_poly_mat_clear(nmod_list_poly_mat_t A);

void nmod_list_poly_mat_init_set(nmod_list_poly_mat_t res, const nmod_poly_mat_t F, slong length);

void nmod_list_poly_mat_set(nmod_list_poly_mat_t res, const nmod_poly_mat_t F);

void nmod_list_poly_mat_to_poly_mat(nmod_poly_mat_t res, const nmod_list_poly_mat_t F);

void nmod_list_poly_mat_naive_mul_coef(nmod_mat_t res, const nmod_list_poly_mat_t A,
				       const nmod_list_poly_mat_t B, slong k);

void structured_list_multiplication_blocks(nmod_list_poly_mat_t res, const nmod_mat_t A,
					    const slong *perm, slong rank, slong sigma);

void nmod_list_poly_mat_to_poly_mat(nmod_poly_mat_t res, const nmod_list_poly_mat_t F);


#endif /* NMOD_LIST_POLY_MAT_H */

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
