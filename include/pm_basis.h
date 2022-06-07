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

#ifndef PM_BASIS_H
#define PM_BASIS_H

#include <stdint.h>
#include <nmod_mat.h>
#include <nmod_poly_mat.h>
#include <perm.h> 

void Basis(nmod_poly_mat_t res, nmod_mat_t mat, int64_t *shift,
	   slong rdim, slong cdim, slong prime);
#endif /* PM_BASIS_H */

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
