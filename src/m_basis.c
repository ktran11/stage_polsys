#include "m_basis.h"

/** static void structured_multiplication_blocks(nmod_poly_mat_t res, const nmod_mat_t A,
 *                                               const slong *perm, slong rank,
 *                                               slong rdim, slong cdim, slong prime)
 * 
 * This function compute the multiplication of specific polynomials matrix
 * It takes an intrger r, A a nmod_mat_t, res a nmod_poly_mat_t
 * and the permutation perm. 
 * It will compute the mutiplication of
 * M = perm^(-1) * [[x, 0], [A, 1]] * perm and res = [[R1],[R2]] 
 * Stocks in res the result
 *
 */
static void structured_multiplication_blocks(nmod_poly_mat_t res, const nmod_mat_t A,
					     const slong *perm, slong rank,
					     slong rdim, slong cdim,slong prime)
{
  nmod_poly_mat_t R1, R2, R1_cp, R2_cp, A_poly;
  nmod_poly_t X;

  slong *inv_perm = _perm_init(rdim);
  
  // Multiplication for matrix top 

  apply_perm_rows_to_poly_matrix(res, perm, rdim);

  nmod_poly_mat_window_init(R1, res, 0, 0, rank, cdim);
  nmod_poly_mat_init_set(R1_cp, R1);

  nmod_poly_init(X, prime);
  nmod_poly_set_coeff_ui(X, 1, 1);

  nmod_poly_mat_shift(R1, 1);
  // Multiplication for matrix bottom
  // (nmod_mat_t) A to (nmod_poly_mat_t) A_poly
  nmod_poly_set_coeff_ui(X, 1, 0);

  slong alloc;
  nmod_poly_mat_init(A_poly, rdim - rank, rank, prime);
  for (slong i = 0; i < rdim - rank; i++)
    for (slong j = 0; j < rank; j++)
      {
	alloc = (slong) nmod_mat_get_entry(A, i, j);
	nmod_poly_set_coeff_ui(X, 0, alloc);
	nmod_poly_set(nmod_poly_mat_entry(A_poly, i, j), X); 
      }
  
  nmod_poly_clear(X);
  
  nmod_poly_mat_window_init(R2, res, rank, 0, rdim, cdim);
  nmod_poly_mat_init_set(R2_cp, R2);

  nmod_poly_mat_mul(R2, A_poly, R1_cp);
  nmod_poly_mat_clear(A_poly);
  nmod_poly_mat_clear(R1_cp);
  
  nmod_poly_mat_add(R2, R2, R2_cp);

  nmod_poly_mat_clear(R2_cp);
  nmod_poly_mat_window_clear(R1);
  nmod_poly_mat_window_clear(R2);

  _perm_inv(inv_perm, perm, rdim);
  apply_perm_rows_to_poly_matrix(res, inv_perm, rdim);
  _perm_clear(inv_perm);
}

void M_basis(nmod_poly_mat_t res, int64_t *res_shifts,
	     const nmod_poly_mat_t F, uint64_t sigma, const int64_t *shifts)

{
  slong rdim = F->r, cdim = F->c;
  mp_limb_t prime = F->modulus;
  nmod_poly_mat_t F_prime;
  nmod_mat_t A_k, constant_mat; 
  slong rank, *perm;
  perm = _perm_init(rdim);

  nmod_mat_init(constant_mat, rdim, cdim, prime);

  nmod_poly_mat_init_set(F_prime, F);
  coefficient_matrix(constant_mat, F_prime, 0); //Compute F mod x

  rank = Basis_for_M_basis(A_k, res_shifts, perm, constant_mat, shifts);

  nmod_poly_mat_one(res);
  structured_multiplication_blocks(res, A_k, perm, rank, rdim, rdim, prime); //Compute P0

  for (uint64_t k = 1; k < sigma; k++)
    {
      // doing the operation x^(-k) F_prime mod x     
      structured_multiplication_blocks(F_prime, A_k, perm, rank, rdim, cdim, prime);
      coefficient_matrix(constant_mat, F_prime, k);

      rank = Basis_for_M_basis(A_k, res_shifts, perm, constant_mat, res_shifts);
 
      structured_multiplication_blocks(res, A_k, perm, rank, rdim, rdim, prime);

    }

  nmod_poly_mat_clear(F_prime);
  nmod_mat_clear(constant_mat);
  nmod_mat_clear(A_k);
  _perm_clear(perm);
  
}

void M_basisII(nmod_poly_mat_t res, int64_t *res_shifts,
	       const nmod_poly_mat_t F, uint64_t sigma, const int64_t *shifts)
{
  slong rdim = F->r, cdim = F->c;
  mp_limb_t prime = F->modulus;
  nmod_list_poly_mat_t F_prime, res_list_repr;
  nmod_mat_t A_k, constant_mat; 
  slong rank, *perm;

  perm = _perm_init(rdim);

  nmod_mat_init(constant_mat, rdim, cdim, prime);

  nmod_list_poly_mat_init_set(F_prime, F, 0);

  nmod_list_poly_mat_get_coef(constant_mat, F_prime, 0);

  rank = Basis_for_M_basis(A_k, res_shifts, perm, constant_mat, shifts);

  nmod_poly_mat_one(res);
  structured_multiplication_blocks(res, A_k, perm, rank, rdim, rdim, prime);
  //Compute P0 = inv_perm * [[x,0],[A0,1]]* perm  * res

  for (uint64_t k = 1; k < sigma; k++)
    {
      nmod_list_poly_mat_init_set(res_list_repr, res, 0);

      nmod_list_poly_mat_naive_mul_coef(constant_mat, res_list_repr, F_prime, k);

      rank = Basis_for_M_basis(A_k, res_shifts, perm, constant_mat, res_shifts);

      structured_multiplication_blocks(res, A_k, perm, rank, rdim, rdim, prime);
      nmod_list_poly_mat_clear(res_list_repr);
    
    }
  
  nmod_list_poly_mat_clear(F_prime);
  nmod_mat_clear(constant_mat);
  nmod_mat_clear(A_k);
  _perm_clear(perm);  
}

void M_basisIII(nmod_poly_mat_t res, int64_t *res_shifts,
	       const nmod_poly_mat_t F, uint64_t sigma, const int64_t *shifts)
{
  slong rdim = F->r, cdim = F->c;
  mp_limb_t prime = F->modulus;
  nmod_list_poly_mat_t F_prime;
  nmod_mat_t A_k, constant_mat; 
  slong rank, *perm;

  perm = _perm_init(rdim);
  nmod_mat_init(constant_mat, rdim, cdim, prime);

  nmod_list_poly_mat_init_set(F_prime, F, sigma);

  nmod_list_poly_mat_get_coef(constant_mat, F_prime, 0);

  rank = Basis_for_M_basis(A_k, res_shifts, perm, constant_mat, shifts);

  nmod_poly_mat_one(res);
  structured_multiplication_blocks(res, A_k, perm, rank, rdim, rdim, prime);
  //Compute P0 = inv_perm * [[x,0],[A0,1]]* perm  * res

  for (uint64_t k = 1; k < sigma; k++)
    {
      structured_list_multiplication_blocks(F_prime, A_k, perm, rank, sigma);
      nmod_list_poly_mat_get_coef(constant_mat, F_prime, k);

      rank = Basis_for_M_basis(A_k, res_shifts, perm, constant_mat, res_shifts);

      structured_multiplication_blocks(res, A_k, perm, rank, rdim, rdim, prime);
    }
  nmod_list_poly_mat_clear(F_prime);
  nmod_mat_clear(constant_mat);
  nmod_mat_clear(A_k);
  _perm_clear(perm);  
}

void M_basisIV(nmod_poly_mat_t res, int64_t *res_shifts,
	       const nmod_poly_mat_t F, uint64_t sigma, const int64_t *shifts)
{
  slong rdim = F->r, cdim = F->c;
  mp_limb_t prime = F->modulus;
  nmod_list_poly_mat_t F_prime, res_prime;
  nmod_mat_t A_k, constant_mat; 
  slong rank, *perm;

  perm = _perm_init(rdim);
  nmod_mat_init(constant_mat, rdim, cdim, prime);

  nmod_list_poly_mat_init_set(F_prime, F, sigma);

  nmod_list_poly_mat_get_coef(constant_mat, F_prime, 0);

  rank = Basis_for_M_basis(A_k, res_shifts, perm, constant_mat, shifts);

  nmod_poly_mat_one(res);
  structured_multiplication_blocks(res, A_k, perm, rank, rdim, rdim, prime);
  nmod_list_poly_mat_init_set(res_prime, res, sigma);
  //Compute P0 = inv_perm * [[x,0],[A0,1]]* perm  * res

  for (uint64_t k = 1; k < sigma; k++)
    {
      structured_list_multiplication_blocks(F_prime, A_k, perm, rank, sigma);
      nmod_list_poly_mat_get_coef(constant_mat, F_prime, k);

      rank = Basis_for_M_basis(A_k, res_shifts, perm, constant_mat, res_shifts);

      structured_list_multiplication_blocks(res_prime, A_k, perm, rank, sigma);
    }

  nmod_list_poly_mat_to_poly_mat(res, res_prime);
  nmod_list_poly_mat_clear(F_prime);
  nmod_list_poly_mat_clear(res_prime);
  nmod_mat_clear(constant_mat);
  nmod_mat_clear(A_k);
  _perm_clear(perm);  
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
