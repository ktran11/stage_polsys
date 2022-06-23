/**
 * \file pm_basis.h
 * \brief pm basis functions
 * \author Kevin Tran
 * \version 1.0
 * \date June 6 2022
 *
 * Compute PM Basis
 *
 */

#ifndef PM_BASIS_H
#define PM_BASIS_H

#include "matpol.h"
#include "m_basis.h"
#include "perm_operator.h"
#include "nmod_list_poly_mat.h"

#include <stdint.h>
#include <nmod_mat.h>
#include <nmod_poly_mat.h>

#define INITIAL_CASE_BOUND 1

void PM_basis(nmod_poly_mat_t res, int64_t *res_shifts,
	      const nmod_poly_mat_t F, uint64_t sigma, const int64_t *shifts);

#endif /* PM_BASIS_H */

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
