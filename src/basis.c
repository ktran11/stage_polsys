#include "basis.h"

void Basis(nmod_poly_mat_t res, int64_t *res_shifts,
	   const nmod_mat_t mat, const int64_t *shifts)
{
  slong rdim = mat->r;
  mp_limb_t prime = mat->mod.n;
  
  slong  i, j, rank_mat, rank_kernel, *P, *P_inv, *comp, *comp_inv, *perm = _perm_init(rdim);
  nmod_mat_t K, mat_cp;
  nmod_poly_t One, constant;

  sort_and_create_perm(perm, shifts, rdim);

  /* LU operation on pi*mat and extraction  L = [ [Lr, 0], [G, I] ] */

  nmod_mat_init_set(mat_cp, mat);

  apply_perm_rows_to_matrix(mat_cp, perm, rdim);

  P = _perm_init(rdim);

  rank_kernel = nmod_mat_left_nullspace_compact(K, P, mat_cp);
  rank_mat = rdim - rank_kernel;
  nmod_mat_clear(mat_cp);
  
  P_inv = _perm_init(rdim);
  for (i = 0; i < rank_kernel; i++)
    P_inv[rank_mat + i] = P[i];
  for (i = 0; i < rank_mat; i++)
    P_inv[i] = P[rank_kernel + i];

  _perm_inv(P_inv, P_inv, rdim);
  
  /* Computation of the block matrix  [ [xIr, 0], [K, I_{rdim - rank}] ] */

  nmod_poly_mat_zero(res);
  
  nmod_poly_init(One, prime);
  nmod_poly_set_coeff_ui(One, 0, 1); // 1
  
  for (i = rank_mat; i < rdim; i++)
    nmod_poly_set(nmod_poly_mat_entry(res, i, i), One);

  nmod_poly_shift_left(One, One, 1); // x
  for (i = 0; i < rank_mat; i++)
    nmod_poly_set(nmod_poly_mat_entry(res, i, i), One);

  nmod_poly_clear(One);

  slong alloc;
  nmod_poly_init(constant, prime);
  for (i = 0; i < rank_kernel; i++)
    for (j = 0; j < rank_mat; j++)
      {
	alloc = (slong) nmod_mat_get_entry(K, i, j);
	nmod_poly_set_coeff_ui(constant, 0, alloc);
	nmod_poly_set(nmod_poly_mat_entry(res, i + rank_mat, j), constant); 
      }

  nmod_mat_clear(K);
  nmod_poly_clear(constant);
 
  /* Multiply by the permutations */
  comp = _perm_init(rdim);
  comp_inv = _perm_init(rdim);
    
  _perm_compose(comp, P_inv, perm, rdim);

  _perm_inv(comp_inv, comp, rdim);
  
  _perm_clear(perm);
  _perm_clear(P);
  _perm_clear(P_inv);
  printf("\nBasis comp permutation\n");
  _perm_print(comp, rdim);

  apply_perm_cols_to_poly_matrix(res, comp, rdim);  
  apply_perm_rows_to_poly_matrix(res, comp_inv, rdim);
  
  /* Compute the new shift */
  slong temp[rdim];
    
  apply_perm_to_vector(temp, shifts, comp, rdim);
  _perm_clear(comp);

  for (i = 0; i < rank_mat; i++)
    temp[i] += 1;

  apply_perm_to_vector(res_shifts, temp, comp_inv, rdim);
  _perm_clear(comp_inv);
}

slong Basis_for_M_basis(nmod_mat_t res, int64_t *res_shifts, slong *res_perm,
		       const nmod_mat_t mat, const int64_t *shifts)

{
  slong rdim = mat->r;  
  slong  i, rank_mat, rank_kernel;
  slong *P, *P_inv, *comp_inv, *perm = _perm_init(rdim);
  nmod_mat_t K, mat_cp;

  sort_and_create_perm(perm, shifts, rdim);

  /* compute left kernel of perm*mat with rank profile */

  nmod_mat_init_set(mat_cp, mat);
  
  apply_perm_rows_to_matrix(mat_cp, perm, rdim);

  P = _perm_init(rdim);
  
  rank_kernel = nmod_mat_left_nullspace_compact(K, P, mat_cp);
  rank_mat = rdim - rank_kernel;

  nmod_mat_clear(mat_cp);
  
  P_inv = _perm_init(rdim);
  for (i = 0; i < rank_kernel; i++)
    P_inv[rank_mat + i] = P[i];
  for (i = 0; i < rank_mat; i++)
    P_inv[i] = P[rank_kernel + i];

  _perm_inv(P_inv, P_inv, rdim);

  nmod_mat_init_set(res, K);
  nmod_mat_clear(K);

  /* Multiply by the permutations */
  comp_inv = _perm_init(rdim);
    
  _perm_compose(res_perm, P_inv, perm, rdim);
  
  _perm_clear(perm);
  _perm_clear(P);
  _perm_clear(P_inv);
 
  /* Compute the new shift */
    slong temp[rdim];
    
  apply_perm_to_vector(temp, shifts, res_perm, rdim);

  for (slong i = 0; i < rank_mat; i++)
    temp[i] += 1;
  
  comp_inv = _perm_init(rdim) ;

  _perm_inv(comp_inv, res_perm, rdim);
  apply_perm_to_vector(res_shifts, temp, comp_inv, rdim);
 
  _perm_clear(comp_inv);
  return rank_mat;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
