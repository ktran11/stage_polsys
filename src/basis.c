#include "basis.h"

/** solve XR = B
 *  where R is tril
 */
static void nmod_mat_solve_tril_right(nmod_mat_t X, const nmod_mat_t R,
				      const nmod_mat_t B, slong prime)
{
  nmod_mat_t R_inv; 
  slong rank = nmod_mat_nrows(R);
  
  nmod_mat_init(R_inv, rank, rank, prime);
  nmod_mat_inv(R_inv, R); 
  nmod_mat_mul(X, B, R_inv);
  
  nmod_mat_neg(X, X);
  
  nmod_mat_clear(R_inv);
}

void Basis(nmod_poly_mat_t res, int64_t *res_shifts,
	   const nmod_mat_t mat, const int64_t *shifts,
	   slong rdim, slong cdim, slong prime)
{
  slong rank, *P, *P_inv, *comp, *comp_inv, *perm = _perm_init(rdim);
  nmod_mat_t G, Lr, prod_G_Lr_inv, mat_cp;
  nmod_poly_t One, constant;

  sort_and_create_perm(perm, shifts, rdim);

  /* LU operation on pi*mat and extraction  L = [ [Lr, 0], [G, I] ] */

  nmod_mat_init_set(mat_cp, mat);

  apply_perm_rows_to_matrix(mat_cp, perm, rdim);
  
  P = _perm_init(rdim);
  rank = nmod_mat_lu(P, mat_cp, 0);

  nmod_mat_init(Lr, rank, rank, prime);
  nmod_mat_one(Lr);
  for(slong i = 1; i < rank; i++)
    for(slong j = 0; j < i; j++)
      nmod_mat_set_entry(Lr, i, j, nmod_mat_entry(mat_cp, i, j));

  nmod_mat_window_init(G, mat_cp, rank, 0, rdim, cdim);

  /* Computation of -G*Lr^(-1) */
  nmod_mat_init(prod_G_Lr_inv, rdim - rank, rank, prime);

  nmod_mat_solve_tril_right(prod_G_Lr_inv, Lr, G, prime);
  
  nmod_mat_clear(Lr); 
  nmod_mat_window_clear(G);
  
  /* Computation of the block matrix  [ [xIr, 0], [-G*Lr^(-1), I_{rdim - rank}] ] */

  nmod_poly_mat_zero(res);
  
  nmod_poly_init(One, prime);
  nmod_poly_set_coeff_ui(One, 0, 1); // 1
  
  for (slong i = rank; i < rdim; i++)
    nmod_poly_set(nmod_poly_mat_entry(res, i, i), One);

  nmod_poly_shift_left(One, One, 1); // x
  for (slong i = 0; i < rank; i++)
    nmod_poly_set(nmod_poly_mat_entry(res, i, i), One);

  nmod_poly_clear(One);

  slong alloc;
  nmod_poly_init(constant, prime);
  for (slong i = 0; i < rdim - rank; i++)
    for (slong j = 0; j < rank; j++)
      {
	alloc = (slong) nmod_mat_get_entry(prod_G_Lr_inv, i, j);
	nmod_poly_set_coeff_ui(constant, 0, alloc);
	nmod_poly_set(nmod_poly_mat_entry(res, i + rank, j), constant); 
      }

  nmod_poly_clear(constant);
  nmod_mat_clear(prod_G_Lr_inv);

  /* Multiply by the permutations */
  P_inv = _perm_init(rdim);
  comp = _perm_init(rdim);
  comp_inv = _perm_init(rdim);
  
  _perm_inv(P_inv, P, rdim);
  
  _perm_compose(comp, P_inv, perm, rdim);

  _perm_inv(comp_inv, comp, rdim);
  
  _perm_clear(perm);
  _perm_clear(P);
  _perm_clear(P_inv);
 
  
  apply_perm_cols_to_poly_matrix(res, comp, rdim);  
  apply_perm_rows_to_poly_matrix(res, comp_inv, rdim);
  
  /* Compute shift */
  slong temp[rdim];
    
  apply_perm_to_vector(temp, shifts, comp, rdim);
  _perm_clear(comp);

  for (slong i = 0; i < rank; i++)
    temp[i] += 1;

  apply_perm_to_vector(res_shifts, temp, comp_inv, rdim);
 
  _perm_clear(comp_inv);
}


slong Basis_for_M_basis(nmod_mat_t res, int64_t *res_shifts, slong *res_perm,
		       const nmod_mat_t mat, const int64_t *shifts,
			slong rdim, slong cdim, slong prime)
{
  slong rank = 0, *P, *P_inv, *comp_inv, *perm = _perm_init(rdim);
  nmod_mat_t G, Lr, Lr_inv, mat_cp;

  sort_and_create_perm(perm, shifts, rdim);

  nmod_mat_init_set(mat_cp, mat);

  /* LU operation on pi*mat and extraction  L = [ [Lr, 0], [G, I] ] */
  apply_perm_rows_to_matrix(mat_cp, perm, rdim);

  P = _perm_init(rdim);
  rank = nmod_mat_lu(P, mat_cp, 0);

  nmod_mat_init(Lr, rank, rank, prime);
  nmod_mat_one(Lr);
  for(slong i = 1; i < rank; i++)
    for(slong j = 0; j < i; j++)
      nmod_mat_set_entry(Lr, i, j, nmod_mat_entry(mat_cp, i, j));

  nmod_mat_window_init(G, mat_cp, rank, 0, rdim, cdim);

  /* Computation of -G*Lr^(-1) */
  nmod_mat_init(Lr_inv, rank, rank, prime);
  nmod_mat_inv(Lr_inv, Lr); // inv can be improved with op in tril

  nmod_mat_clear(Lr);

  nmod_mat_init(res, rdim - rank, rank, prime); 
  nmod_mat_mul(res, G, Lr_inv);
  // can be improved with choice of algorithm (default = Strassen) and multithread
  
  nmod_mat_neg(res, res);

  nmod_mat_clear(Lr_inv);
  nmod_mat_clear(G);

  /* Multiply by the permutations */

  P_inv = _perm_init(rdim);  
  _perm_inv(P_inv, P, rdim);
  
  _perm_compose(res_perm, P_inv, perm, rdim);
  
  _perm_clear(perm);
  _perm_clear(P);
  _perm_clear(P_inv);
  /* Compute shift */
   slong temp[rdim];
    
  apply_perm_to_vector(temp, shifts, res_perm, rdim);

  for (slong i = 0; i < rank; i++)
    temp[i] += 1;
  
  comp_inv = _perm_init(rdim) ;

  _perm_inv(comp_inv, res_perm, rdim);
  apply_perm_to_vector(res_shifts, temp, comp_inv, rdim);
 
  _perm_clear(comp_inv);

  return rank;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
