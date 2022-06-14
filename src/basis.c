#include "basis.h"
#include "matpol.h"
#include <flint.h>


/* Type only for the function Basis */
typedef struct
{
  int64_t value;
  slong ord;
} int_tuple;

static void print_vect_uint64(const uint64_t *vect, slong len)
{
    if (!vect)
        return;
    for (slong i = 0; i < len; i++)
        printf("%lu ", *(vect+i));
    printf("\n");
    return;
}
////////////////// Permutation from sorting vector //////////////////////////

/* Parameter for quicksort */
static int compare(const void *a, const void *b)
{
  int_tuple a_prime = * (const int_tuple *) a;
  int_tuple b_prime = * (const int_tuple *) b;
  return a_prime.value - b_prime.value;
}

/**
 * \brief Will only be used in Basis.
 *
 *  Creates a permutation from the sorting of a shift 
 *
 * \param perm, the result his length must be n
 * \param vec, the shift we want to sort increasly
 * \param n, length of perm and vec
 */
static void sort_and_create_perm(slong *perm, const int64_t *vec, slong n)
{
  int_tuple temp[n];
  for (slong i = 0; i < n; i++)
    {
      temp[i].value = vec[i];
      temp[i].ord = i;
    }

  qsort(temp, n, sizeof(int_tuple), compare);
  for (slong i = 0; i < n; i++)
    perm[temp[i].ord] = i;
    
}

/**
 * \brief A function fund in the doc of nmod_mat_t
 *
 */
static void nmod_poly_mat_swap_cols(nmod_poly_mat_t mat, slong * perm,
			     slong r, slong s)
{
  if (r != s && !nmod_poly_mat_is_empty(mat))
    {
      slong t;
      if (perm)
        {
	  t = perm[s];
	  perm[s] = perm[r];
	  perm[r] = t;
        }

      for (t = 0; t < mat->r; t++)
	nmod_poly_swap(&mat->rows[t][r],&mat->rows[t][s]);
    }
}

static void nmod_poly_mat_swap_rows(nmod_poly_mat_t mat, slong * perm, slong r, slong s)
{
  if (r != s && !nmod_poly_mat_is_empty(mat))
    {
      nmod_poly_struct *u;
      slong t;

      if (perm)
        {
	  t = perm[s];
	  perm[s] = perm[r];
	  perm[r] = t;
        }

      u = mat->rows[s];
      mat->rows[s] = mat->rows[r];
      mat->rows[r] = u;
    }
}


void apply_perm_rows_to_matrix(nmod_mat_t mat, const slong * perm, slong rdim)
{
  slong *perm_copy = _perm_init(rdim); 
  _perm_set(perm_copy, perm, rdim);
  for (slong i = 0; i < rdim; i++)
      for (slong j = 0; j < rdim; j++)
	{
	  if (i == perm_copy[i])
	    break;
	  nmod_mat_swap_rows(mat, perm_copy, i, perm_copy[i]);
	}
	
  _perm_clear(perm_copy);
}

void apply_perm_cols_to_poly_matrix(nmod_poly_mat_t mat, const slong * perm, slong cdim)
{
  slong *perm_copy = _perm_init(cdim);
  _perm_inv(perm_copy, perm, cdim);
  for (slong i = 0; i < cdim; i++)
    {
      for (slong j = 0; j < cdim; j++)
      {
	if (i == perm_copy[i])
	  break;
	nmod_poly_mat_swap_cols(mat, perm_copy, i, perm_copy[i]);
      }
    }
  _perm_clear(perm_copy);
}

void apply_perm_rows_to_poly_matrix(nmod_poly_mat_t mat, const slong *perm, slong rdim)
{
  slong *perm_copy = _perm_init(rdim); 
  _perm_set(perm_copy, perm, rdim);

  for (slong i = 0; i < rdim; i++)
    for (slong j = 0; j < rdim; j++)
      {
	if (i == perm_copy[i])
	  break;
	nmod_poly_mat_swap_rows(mat, perm_copy, i, perm_copy[i]);
      }
  _perm_clear(perm_copy);    
}



void apply_perm_to_vector(int64_t *res, const int64_t *initial_vect,
			  const slong *perm, slong length)
{
  for (slong i = 0; i < length; i++)
    res[perm[i]] = initial_vect[i];
}

void Basis(nmod_poly_mat_t res, int64_t *res_shifts,
	   const nmod_mat_t mat, const int64_t *shifts,
	   slong rdim, slong cdim, slong prime)
{
  slong rank, *P, *P_inv, *comp, *comp_inv, *perm = _perm_init(rdim);
  nmod_mat_t G, Lr, Lr_inv, prod_G_Lr_inv, mat_cp;
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
  nmod_mat_init(Lr_inv, rank, rank, prime);
  nmod_mat_inv(Lr_inv, Lr); // inv can be improved with op in tril

  nmod_mat_clear(Lr);

  nmod_mat_init(prod_G_Lr_inv, rdim - rank, rank, prime);
  nmod_mat_mul(prod_G_Lr_inv, G, Lr_inv);
  // can be improved with choice of algorithm (default = Strassen) and multithread
  
  nmod_mat_neg(prod_G_Lr_inv, prod_G_Lr_inv);
  
  nmod_mat_clear(Lr_inv);
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
