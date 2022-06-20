#include "nmod_list_poly_mat.h"


slong nmod_list_poly_mat_nrows(const nmod_list_poly_mat_t A)
{
  return A->rdim;
}

slong nmod_list_poly_mat_ncols(const nmod_list_poly_mat_t A)
{
  return A->cdim;
}

slong nmod_list_poly_mat_degree(const nmod_list_poly_mat_t A)
{
  return A->degree;
}


mp_limb_t nmod_list_poly_mat_modulus(const nmod_list_poly_mat_t A)
{
  return A->modulus;
}

void nmod_list_poly_mat_init(nmod_list_poly_mat_t res, slong degree, slong length,
			     slong rdim, slong cdim, mp_limb_t modulus)
{

 slong i;

    if (rdim > 0)
      {
	res->mat = (nmod_mat_struct *) flint_malloc(length * sizeof(nmod_mat_struct));
	for (i = 0; i < length; i++)
	  nmod_mat_init(res->mat + i, rdim, cdim, modulus); 
      }
    else
      res->mat = NULL;
    

    res->length = length;
    res->degree = degree;
    res->modulus = modulus;
    res->rdim = rdim;
    res->cdim = cdim;
}

void nmod_list_poly_mat_clear(nmod_list_poly_mat_t A)
{
  if (A == NULL)
    return;
  if (A->mat != NULL)
    {
      for (slong i = 0; i < A->length; i++)
	nmod_mat_clear(A->mat + i);

      flint_free(A->mat);
    }
}

void nmod_list_poly_mat_print(const nmod_list_poly_mat_t A)
{
  for (slong i = 0; i <= A->degree; i++)
    nmod_mat_print_pretty(A->mat + i);
}

void nmod_list_poly_mat_init_set(nmod_list_poly_mat_t res, const nmod_poly_mat_t F, slong length)
{
  if (!F)
    {
      printf("Error your polynomial matrix is NULL");
      return;
    }
  
  slong degree, rdim, cdim;
  mp_limb_t modulus;
  degree = nmod_poly_mat_degree(F);
  rdim = nmod_poly_mat_nrows(F);
  cdim = nmod_poly_mat_ncols(F);
  modulus = nmod_poly_mat_modulus(F);

  if (length < degree)
    nmod_list_poly_mat_init(res, degree, degree + 1, rdim, cdim, modulus);
  else
    nmod_list_poly_mat_init(res, degree, length, rdim, cdim, modulus);

  nmod_mat_t mat;
  nmod_mat_init(mat, rdim, cdim, modulus);
  for (slong i = 0; i <= degree; i++)
    {
      coefficient_matrix(mat, F, i);
      nmod_mat_set(res->mat + i, mat);
    }
    
  
  nmod_mat_clear(mat);
}

void nmod_list_poly_mat_to_poly_mat(nmod_poly_mat_t res, const nmod_list_poly_mat_t F)
{
  if (!F)
    {
      printf("Error your list polynomial matrix is NULL");
      return;
    }
  
  slong degree, rdim, cdim;
  mp_limb_t modulus;
  degree = F->degree;
  rdim = F->rdim;
  cdim = F->cdim;
  modulus = F->modulus;

  nmod_poly_mat_t mat;
  
  nmod_poly_mat_clear(mat);
}

void nmod_list_poly_mat_get_coef(nmod_mat_t res, const nmod_list_poly_mat_t F, slong k)
{
  nmod_mat_set(res, F->mat + k);
}

void nmod_list_poly_mat_naive_mul_coef(nmod_mat_t res, const nmod_list_poly_mat_t A,
				       const nmod_list_poly_mat_t B, slong k)

{
  nmod_mat_t temp;
  
  slong A_rdim, A_cdim, B_rdim, B_cdim;
  mp_limb_t A_modulus, B_modulus;

  A_rdim = nmod_list_poly_mat_nrows(A);
  A_cdim = nmod_list_poly_mat_ncols(A);
  A_modulus = nmod_list_poly_mat_modulus(A);
  
  B_rdim = nmod_list_poly_mat_nrows(B);
  B_cdim = nmod_list_poly_mat_ncols(B);
  B_modulus = nmod_list_poly_mat_modulus(B);
  if (A_modulus != B_modulus)
    {
      printf("\nERROR! Wrong modulus: nmod_list_poly_mat_naive_mul_coef\n");
      return;
    }
  if (A_cdim != B_rdim)
    {
      printf("\nERROR! Wrong shape: nmod_list_poly_mat_naive_mul_coef\n");
      printf("shape A = (%ld, %ld), shape B = (%ld, %ld)\n", A_rdim, A_cdim, B_rdim, B_cdim);
      return;
    }
 
  nmod_mat_zero(res);

  nmod_mat_init(temp, A_rdim, B_cdim, A_modulus);
  for (slong i = 0; i <= k; i++)
    {
      nmod_mat_mul(temp, A->mat + i, B->mat + (k - i));
      nmod_mat_add(res, res, temp);
    }	
}


 void structured_list_multiplication_blocks(nmod_list_poly_mat_t res, const nmod_mat_t A,
					    const slong *perm, slong rank, slong sigma)
{
  slong  rdim = res->rdim, cdim = res->cdim;
  mp_limb_t modulus = res->modulus;
  slong i;
  nmod_mat_t zero, previous_R1, R1, R2, R1_cp, R2_cp;
  nmod_mat_struct *r_i;
  slong *inv_perm = _perm_init(rdim);

  nmod_mat_init(R1_cp, rank, cdim, modulus);
  nmod_mat_init(previous_R1, rank, cdim, modulus);
  nmod_mat_init(R2_cp, rdim - rank, cdim, modulus);
  
  nmod_mat_init(zero, rank, cdim, modulus);  
  nmod_mat_zero(zero);  
  
  for (i = 0; i < sigma; i++)
    {
      r_i = res->mat + i;
      apply_perm_rows_to_matrix(r_i, perm, rdim);
      nmod_mat_window_init(R1, r_i, 0, 0, rank, cdim);
      nmod_mat_init_set(R1_cp, R1);
      
      if (i != 0)
	nmod_mat_set(R1, previous_R1); //it will set the top block of r_{i-1} on r_{i}
      else
	{
	  nmod_mat_set(R1, zero);
	  nmod_mat_clear(zero);
	}
      
      nmod_mat_set(previous_R1, R1_cp);

      nmod_mat_window_init(R2, r_i, rank, 0, rdim, cdim);
      nmod_mat_set(R2_cp, R2);

      nmod_mat_mul(R2, A, R1_cp);
      nmod_mat_add(R2, R2, R2_cp);
    }
  
  nmod_mat_clear(R1_cp);
  nmod_mat_clear(R2_cp);
  nmod_mat_clear(previous_R1);

  nmod_mat_window_clear(R1);
  nmod_mat_window_clear(R2);
  

  _perm_inv(inv_perm, perm, rdim);
  for (i = 0; i <= res->degree; i++)
    {
      r_i = res->mat + i;
      apply_perm_rows_to_matrix(r_i, inv_perm, rdim);
    }

  _perm_clear(inv_perm);
}
