/**
 * \file pm_basis.h
 * \brief pm basis functions
 * \author Kevin Tran
 * \version 1.0
 * \date June 6 2022
 *
 * Compute minimal polynomial of a matrix
 *
 */

#ifndef WIEDEMANN_ALGO_H
#define WIEDEMANN_ALGO_H

#include "matpol.h"
#include "pm_basis.h"

#include <stdint.h>
#include <nmod_mat.h>
#include <nmod_poly_mat.h>

void krylov_sequence(nmod_mat_t res, const nmod_mat_t A, const nmod_mat_t v,
		     slong n);

void Berlekamp_Massey(nmod_poly_t res, const nmod_poly_t P);

void Wiedemann(nmod_poly_t res, const nmod_mat_t A);

void Block_Wiedemann(nmod_poly_t res, const nmod_mat_t A);

#endif /* WIEDEMANN_ALGO_H */

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
