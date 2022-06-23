/**
 * \file m_basis.h
 * \brief m basis functions
 * \author Kevin Tran
 * \version 1.0
 * \date June 6 2022
 *
 * Compute M Basis
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


/** static void structured_multiplication_blocks(nmod_poly_mat_t res, const nmod_mat_t A,
 *                                               const slong *perm, slong rank)
 * 
 * This function compute the multiplication of specific polynomials matrix
 * It takes an intrger r, A a nmod_mat_t, res a nmod_poly_mat_t
 * and the permutation perm. 
 * It will compute the mutiplication of
 * M = perm^(-1) * [[x, 0], [A, 1]] * perm and res = [[R1],[R2]] 
 * Stocks in res the result
 *
 */
void structured_multiplication_blocks(nmod_poly_mat_t res, const nmod_mat_t A,
				      const slong *perm, slong rank);

/** void M_basis(nmod_poly_mat_t res, int64_t *res_shifts,
 *	         const nmod_poly_mat_t F, uint64_t sigma, const int64_t *shifts);
 * 
 * Compute the minimal approximant basis of F for the order sigma and the shifts shifts
 *
 * It used the structured multiplication blocks to compute x^{-k} P_{k-1} F mod x and to compute
 * P_{k} = M P_{k-1}
 * Use polynomial matrix multiplication
 *
 */
void M_basis(nmod_poly_mat_t res, int64_t *res_shifts,
	     const nmod_poly_mat_t F, uint64_t sigma, const int64_t *shifts);


/** void M_basisII(nmod_poly_mat_t res, int64_t *res_shifts,
 *	         const nmod_poly_mat_t F, uint64_t sigma, const int64_t *shifts);
 * 
 * Compute the minimal approximant basis of F for the order sigma and the shifts shifts
 *
 * F \in K[x]^{nxm} <-> F_prime K^{mxn}[x] will be FIX and
 * compute iteratively P_{k} \in K[x]^{mxm}, then will transform to P_prime_{k} \in K^{mxm}[x] 
 * and compute x^{-k} P_{k-1} F mod x with a naive polynomial multiplication
 *
 * Use naive polynomial multiplication
 */
void M_basisII(nmod_poly_mat_t res, int64_t *res_shifts,
	       const nmod_poly_mat_t F, uint64_t sigma, const int64_t *shifts);

/** void M_basisIII(nmod_poly_mat_t res, int64_t *res_shifts,
 *	         const nmod_poly_mat_t F, uint64_t sigma, const int64_t *shifts);
 * 
 * Compute the minimal approximant basis of F for the order sigma and the shifts shifts
 *
 * F \in K[x]^{nxm} <-> F_prime K^{mxn}[x]  and
 * Compute P_{k-1} F iteratively with a the list_structured_multiplication_blocks 
 * (in nmod_list_poly_mat.h)
 * Compute P_{k-1} \in K[x]^{mxm} with the structured_multiplication_blocks
 *
 * Use structured_multiplication_blocks and list_structured_multiplication_blocks 
 */
void M_basisIII(nmod_poly_mat_t res, int64_t *res_shifts,
		const nmod_poly_mat_t F, uint64_t sigma, const int64_t *shifts);

/** void M_basisIV(nmod_poly_mat_t res, int64_t *res_shifts,
 *	         const nmod_poly_mat_t F, uint64_t sigma, const int64_t *shifts);
 * 
 * Compute the minimal approximant basis of F for the order sigma and the shifts shifts
 *
 * F \in K[x]^{nxm} <-> F_prime K^{mxn}[x]  and
 * Compute P_{k-1} F iteratively with a the list_structured_multiplication_blocks 
 * (in nmod_list_poly_mat.h)
 * Compute P_{k-1} \in K^{mxm}[x] with the list_structured_multiplication_blocks
 * Transform P_{sigma - 1} as a poly_mat_t 
 *
 * Use list_structured_multiplication_blocks 
 */
void M_basisIV(nmod_poly_mat_t res, int64_t *res_shifts,
	       const nmod_poly_mat_t F, uint64_t sigma, const int64_t *shifts);


/** void M_basisIII(nmod_poly_mat_t res, int64_t *res_shifts,
 *	         const nmod_poly_mat_t F, uint64_t sigma, const int64_t *shifts);
 * 
 * Compute the minimal approximant basis of F for the order sigma and the shifts shifts
 *
 * F \in K[x]^{nxm} <-> F_prime K^{mxn}[x] FIX and
 * Compute P_{k-1} \in K[x]^{mxm} with the list structured_multiplication_blocks
 * Compute x^{-k} P_{k-1} F mod x iteratively with a naive polynomial multiplication  
 *
 * Use naive polynomial multiplication and list_structured_multiplication_blocks 
 */
void M_basisV(nmod_poly_mat_t res, int64_t *res_shifts,
	      const nmod_poly_mat_t F, uint64_t sigma, const int64_t *shifts);
#endif /* M_BASIS_H */

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
