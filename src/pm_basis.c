#include "pm_basis.h"

void middle_product(nmod_poly_mat_t res, const nmod_poly_mat_t A,
		    const nmod_poly_mat_t B, slong d, slong h)
{
  slong rdim = A->r, cdim = B->c;
  mp_limb_t prime = A->modulus;
  nmod_poly_mat_t C;
  nmod_poly_struct *P; 
  nmod_poly_t R; 
  
  nmod_poly_mat_init(C, rdim, cdim, prime);

  nmod_poly_mat_mul(C, A, B);

  nmod_poly_init(R, prime);
  for (slong k = 0; k < rdim; k++)
    for (slong l = 0; l < cdim; l++)
      {
	P = nmod_poly_mat_entry(C, k, l);
	for (slong i = d - 1; i < h - d; i++)
	  nmod_poly_set_coeff_ui(R, i - d + 1,
				 nmod_poly_get_coeff_ui(P, i));
	nmod_poly_set(P, R); 
      }
}


void PM_basis(nmod_poly_mat_t res, int64_t *res_shifts,
	      const nmod_poly_mat_t F, uint64_t sigma, const int64_t *shifts)
{
  slong rdim = F->r, cdim = F->c;
  mp_limb_t prime = F->modulus;
  nmod_poly_mat_t Pl, Ph, F_prime;
  slong ul[rdim];
  
  if (sigma == INITIAL_CASE_BOUND)
    M_basis(res, res_shifts, F, sigma, shifts);

  nmod_poly_mat_init(Pl, rdim, rdim, prime);
  PM_basis(Pl, ul, F, sigma / 2, shifts);

  nmod_poly_mat_init(F_prime, rdim, cdim, prime);
  middle_product(F_prime, Pl, F, sigma / 2 + 1, sigma);

  nmod_poly_mat_init(Ph, rdim, rdim, prime);
  if (sigma % 2 == 0)
    PM_basis(Ph, res_shifts, F_prime, sigma / 2, shifts);
  else
    PM_basis(Ph, res_shifts, F_prime, sigma / 2 + 1, shifts);

  nmod_poly_mat_mul(res, Ph, Pl);

  nmod_poly_mat_clear(Ph);
  nmod_poly_mat_clear(Pl);
  nmod_poly_mat_clear(F_prime);   
}
