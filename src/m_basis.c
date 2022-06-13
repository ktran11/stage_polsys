#include "m_basis.h"

#include <profiler.h>

#include "matpol.h"
#include "pm_basis.h"

/**
 * This function compute the multiplication of specific polynomials matrix
 * It takes an intrger r, A a nmod_mat_t, res a nmod_poly_mat_t
 * and the permutation perm. 
 * It will compute the mutiplication of
 * M = perm^(-1) * [[x, 0], [A, 1]] * perm and res = [[R1],[R2]] 
 */
static void structured_multiplication_blocks(nmod_poly_mat_t res, const nmod_mat_t A,
					     const slong *perm, slong rank,
					     slong rdim, slong cdim, slong prime)
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

  nmod_poly_mat_scalar_mul_nmod_poly(R1, R1_cp, X);

  // Multiplication for matrix bottom
  // (nmod_mat_t) A to (nmod_poly_mat_t) A_poly
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
	     const nmod_poly_mat_t F, const uint64_t sigma, int64_t *shifts,
	     slong rdim, slong cdim, slong prime)

{
  nmod_poly_mat_t F_prime;
  nmod_mat_t A_k, constant_mat; 
  slong rank, *perm;
  
  perm = _perm_init(rdim);

  nmod_mat_init(constant_mat, rdim, cdim, prime);

  nmod_poly_mat_init_set(F_prime, F);
  coefficient_matrix(constant_mat, F_prime, 0); //Compute F mod x 

  rank = Basis_for_M_basis(A_k, res_shifts, perm, constant_mat, shifts, rdim, cdim, prime);

  nmod_poly_mat_one(res);
  structured_multiplication_blocks(res, A_k, perm, rank, rdim, cdim, prime); //Compute P0
  
  for (uint64_t k = 1; k < sigma; k++)
    {
      structured_multiplication_blocks(F_prime, A_k, perm, rank, rdim, cdim, prime);
      coefficient_matrix(constant_mat, F_prime, k); // doing the operation x^(-k) F_prime mod x     
      rank = Basis_for_M_basis(A_k, res_shifts, perm, constant_mat, res_shifts, rdim, cdim, prime);
      structured_multiplication_blocks(res, A_k, perm, rank, rdim, cdim, prime);
    }

  nmod_poly_mat_clear(F_prime);
  nmod_mat_clear(constant_mat);
  nmod_mat_clear(A_k);
  _perm_clear(perm);
  
}

int main(void)
{
  nmod_poly_mat_t mat, copy_mat, first_res, second_res;
  nmod_mat_t A;
  slong rdim = 100, cdim = 50, prime = 101, rank = 10, len = 100, *shifts,
    first_shifts[rdim], second_shifts[rdim];  
  timeit_t t0;
  
  nmod_poly_mat_init(mat, rdim, cdim, prime);
  //nmod_mat_init(A, rdim - rank, rank, prime);
  nmod_mat_init(A, rdim, rdim, prime);
  
  flint_rand_t state;
  flint_randinit(state);

  nmod_poly_mat_randtest(mat, state, len);
  nmod_poly_mat_init_set(copy_mat, mat);

  //nmod_mat_randtest(A, state);
  
  shifts= _perm_init(rdim);
  _perm_randtest(shifts, rdim, state);


  M_basis(first_res, first_shifts, mat, 0, shifts, rdim, cdim, prime);

  char *x = malloc(150);

  nmod_poly_mat_print(first_res, x);

  coefficient_matrix(A, mat, 0);

  Basis(second_res, second_shifts, A, shifts, rdim, cdim, prime);
  
  nmod_poly_mat_print(second_res,x);
  /**
  timeit_start(t0);

  structured_multiplication_blocks(mat, A, perm, rank, rdim, cdim, prime);

  timeit_stop(t0);

  flint_printf("cpu = %wd ms  wall = %wd ms\n", t0->cpu, t0->wall);

  nmod_poly_mat_init(G, rdim, rdim, prime);
  nmod_poly_mat_zero(G);

  nmod_poly_t One;
  
  nmod_poly_init(One, prime);
  nmod_poly_set_coeff_ui(One, 0, 1); // 1
  
  for (slong i = rank; i < rdim; i++)
    nmod_poly_set(nmod_poly_mat_entry(G, i, i), One);

  nmod_poly_shift_left(One, One, 1); // x
  for (slong i = 0; i < rank; i++)
    nmod_poly_set(nmod_poly_mat_entry(G, i, i), One);

  nmod_poly_clear(One);

  nmod_poly_t constant;
  slong alloc;
  nmod_poly_init(constant, prime);
  for (slong i = 0; i < rdim - rank; i++)
    for (slong j = 0; j < rank; j++)
      {
	alloc = (slong) nmod_mat_get_entry(A, i, j);
	nmod_poly_set_coeff_ui(constant, 0, alloc);
	nmod_poly_set(nmod_poly_mat_entry(G, i + rank, j), constant); 
      }

  timeit_stop(t0);
  
  apply_perm_rows_to_poly_matrix(copy_mat, perm, rdim);
  nmod_poly_mat_mul(copy_mat, G, copy_mat);

  flint_printf("cpu = %wd ms  wall = %wd ms\n", t0->cpu, t0->wall);

  printf("\n mat prod is equal %d\n", nmod_poly_mat_equal(mat, copy_mat));;
  nmod_poly_mat_clear(mat);
  nmod_poly_mat_clear(copy_mat);

  nmod_mat_clear(A);
  _perm_clear(perm);
    */
  return 0;
}

