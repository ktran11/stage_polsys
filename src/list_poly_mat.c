#include "nmod_list_poly_mat.h"

void nmod_list_poly_mat_init(nmod_list_poly_mat_t res, slong degree,
			     slong rdim, slong cdim, mp_limb_t modulus)
{

 slong i;

    if (rdim > 0)
      {
	res->list_mat = (nmod_mat_struct *) flint_malloc((degree + 1) * sizeof(nmod_mat_struct));
	for (i = 0; i <= degree; i++)
	  nmod_mat_init(res->list_mat[i], rdim, cdim, modulus); 
      }
    else
      res->list_mat = NULL;
    

    res->degree = degree;
    res->modulus = modulus;
    res->rdim = rdim;
    res->cdim = cdim;
}

void nmod_list_poly_mat_set(nmod_list_poly_mat_t res, const nmod_poly_mat_t F)
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
  
  nmod_list_poly_mat_init(res, degree, rdim, cdim, modulus);

  nmod_mat_t mat;
  nmod_mat_init(mat, rdim, cdim, modulus);
  for (slong i = 0; i <= degree; i++)
    {
      coefficient_matrix(mat, F, i);
      res->list_mat[i] = mat;
    }
  nmod_mat_clear(mat);
}

void nmod_list_poly_mat_set_constant(nmod_list_poly_mat_t res, const nmod_mat_t F)
{
  if (!F)
    {
      printf("ERROR! Your matrix is NULL: nmod_list_poly_mat_set_constant");
      return;
    }
  slong rdim, cdim;
  mp_limb_t modulus;
  rdim = nmod_mat_nrows(F);
  cdim = nmod_mat_ncols(F);
  modulus = nmod_mat_modulus(F);
  
  nmod_list_poly_mat_init(res, 0, rdim, cdim, modulus);

  res->list_mat[0] = mat;
}

void nmod_list_poly_mat_naive_mul_coef(nmod_mat_t res, const nmod_list_poly_mat_t A,
				       const nmod_list_poly_mat_t B, slong k)

{
  nmod_mat_t res, temp;
  
  slong A_degree, A_rdim, A_cdim, B_degree, B_rdim, B_cdim;
  mp_limb_t A_modulus, B_modulus;

  A_degree = nmod_poly_mat_degree(A);
  A_rdim = nmod_poly_mat_nrows(A);
  A_cdim = nmod_poly_mat_ncols(A);
  A_modulus = nmod_poly_mat_modulus(A);
  
  B_degree = nmod_poly_mat_degree(B);
  B_rdim = nmod_poly_mat_nrows(B);
  B_cdim = nmod_poly_mat_ncols(B);
  B_modulus = nmod_poly_mat_modulus(B);
  if (A_modulus != B_modulus)
    {
      printf("ERROR! Wrong modulus: nmod_list_poly_mat_naive_mul_coef");
      return;
    }
  if (A_cdim != B_rdim)
    {
      printf("ERROR! Wrong shape: nmod_list_poly_mat_naive_mul_coef");
      return;
    }
 
  nmod_mat_zero(res);

  nmod_mat_init(temp, A_rdim, B_cdim, A_modulus);
  for (slong i = 0; i < A_degree; i++)
    for (slong j = 0; j < B_degree; j++)
      if (i + j = k)
	{
	  nmod_mat_mul(temp, A->list_mat[i], B->list_mat[j]);
	  nmod_mat_add(res, res, temp);
	}
}
