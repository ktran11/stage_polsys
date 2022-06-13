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


//////////// Permutation on columns of Poly Matrix or rows Matrix or Vector //////////////////////

/**
 * \brief A function fund in the doc of nmod_mat_t
 *
 */
void nmod_poly_mat_swap_cols(nmod_poly_mat_t mat, slong * perm,
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

void nmod_poly_mat_swap_rows(nmod_poly_mat_t mat, slong * perm, slong r, slong s)
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


/**
 * \brief Do the multiplication permutation matrix (P*A) for the function Basis
 *
 * \param A, the input and result matrix
 * \param P, the permutation 
 * \param n, length of P and row dimension of A
 */
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
  _perm_set(perm_copy, perm, cdim);  
  for (slong i = 0; i < cdim; i++)
    for (slong j = 0; j < cdim; j++)
      {
	if (i == perm_copy[i])
	  break;
	nmod_poly_mat_swap_cols(mat, perm_copy, i, perm_copy[i]);
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

//////////////////////////////////////////////////////////////////////////////////////////////

/**
 * \fn void Basis(nmod_poly_mat_t res, int64_t *res_shift,
	   nmod_mat_t mat, int64_t *shift,
	   slong rdim, slong cdim, slong prime);
 * \brief set on res the minimal approximant basis of mat for order 1
 and set res_shift the shift row degree of res
 * 
 *
 * \param All parameters are supposed init
 * \param res a polynomial matrix of dimensions (rdim, rdim)
 * \param res_shift a vector of length rdim
 * \param mat a polynomial matrix (!) mat is modified (!)
 * \param shift a vector of length rdim
 * \param rdim the row's number of mat
 * \param cdim the column's number of mat
 * \param prime the modulus
 *
 */
void Basis(nmod_poly_mat_t res, int64_t *res_shifts,
	   const nmod_mat_t mat, const int64_t *shifts,
	   slong rdim, slong cdim, slong prime);

slong Basis_for_M_basis(nmod_mat_t res, int64_t *res_shifts, slong *res_perm,
		       const nmod_mat_t mat, const int64_t *shifts,
			slong rdim, slong cdim, slong prime);

#endif /* PM_BASIS_H */

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
