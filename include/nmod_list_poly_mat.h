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

#ifndef LIST_POLY_MAT_H
#define LIST_POLY_MAT_H

#include "matpol.h"

#include <stdint.h>
#include <nmod_mat.h>
#include <nmod_poly_mat.h>

typedef struct
{
  slong degree;
  nmod_mat_t *mat;
  slong rdim;
  slong cdim;
  mp_limb_t modulus;
} nmod_list_poly_mat_struct;

typedef nmod_list_poly_mat_struct nmod_list_poly_mat_t[1];

slong nmod_list_poly_mat_nrows(const nmod_poly_mat_t list_mat)
{
  return list_mat->rdim;
}

void nmod_list_poly_mat_init(nmod_list_poly_mat_t res);


#endif /* MATPOL_H */

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
