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

void nmod_list_poly_mat_init(nmod_list_poly_mat_t res, slong degree,
			     slong rdim, slong cdim, mp_limb_t modulus)
{

 slong i;

    if (rdim > 0)
      {
	res->mat = (nmod_mat_struct *) flint_malloc((degree + 1) * sizeof(nmod_mat_struct));
	for (i = 0; i <= degree; i++)
	  nmod_mat_init(&res->mat[i], rdim, cdim, modulus); 
      }
    else
      res->mat = NULL;
    

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
      for (slong i = 0; i <= A->degree; i++)
	nmod_mat_clear(A->mat + i);

      flint_free(A->mat);
    }
}

void nmod_list_poly_mat_print(const nmod_list_poly_mat_t A)
{
  for (slong i = 0; i <= A->degree; i++)
    {
      nmod_mat_print_pretty(A->mat + i);
    }
}

void nmod_list_poly_mat_init_set(nmod_list_poly_mat_t res, const nmod_poly_mat_t F)
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
      nmod_mat_set(&res->mat[i], mat);
    }
  nmod_mat_clear(mat);
}

void nmod_list_poly_mat_get_coef(nmod_mat_t res, const nmod_list_poly_mat_t F, slong k)
{
  nmod_mat_set(res, &F->mat[k]);
}

void nmod_list_poly_mat_set_constant(nmod_list_poly_mat_t res, const nmod_mat_t F,
				     mp_limb_t modulus)
{
  if (!F)
    {
      printf("ERROR! Your matrix is NULL: nmod_list_poly_mat_set_constant");
      return;
    }
  slong rdim, cdim;
  rdim = nmod_mat_nrows(F);
  cdim = nmod_mat_ncols(F);
  
  nmod_list_poly_mat_init(res, 0, rdim, cdim, modulus);

  nmod_mat_set(&res->mat[0], F);
}

void nmod_list_poly_mat_naive_mul_coef(nmod_mat_t res, const nmod_list_poly_mat_t A,
				       const nmod_list_poly_mat_t B, slong k)

{
  nmod_mat_t temp;
  
  slong A_degree, A_rdim, A_cdim, B_degree, B_rdim, B_cdim;
  mp_limb_t A_modulus, B_modulus;

  A_degree = nmod_list_poly_mat_degree(A);
  A_rdim = nmod_list_poly_mat_nrows(A);
  A_cdim = nmod_list_poly_mat_ncols(A);
  A_modulus = nmod_list_poly_mat_modulus(A);
  
  B_degree = nmod_list_poly_mat_degree(B);
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
  for (slong i = 0; i <= A_degree; i++)
    for (slong j = 0; j <= B_degree; j++)
      if (i + j == k)
	{
	  nmod_mat_mul(temp, A->mat + i, B->mat + j);
	  nmod_mat_add(res, res, temp);
	}
}
