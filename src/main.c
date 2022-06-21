#include "main.h"
#include "nmod_list_poly_mat.h"
#include "m_basis.h"
#include "matpol.h"
#include "basis.h"
#include "perm_operator.h"
#include <time.h>
#include <profiler.h>

#define GREAT_PRIME_30_BITS 536870923
#define GREAT_PRIME_60_BITS 576460752303423619


/** Test
 * Verify if all versions of M_Basis give the same results, test time
 */
int test_mbasis(void)
{
  nmod_poly_mat_t mat, first_res, second_res, third_res, fourth_res, fifth_res, sixth_res;
  nmod_mat_t A;
  slong rdim, cdim, prime = 3, sigma;
  
  printf("Saisir rdim: ");
  scanf("%ld", &rdim);
  printf("Saisir cdim: ");
  scanf("%ld", &cdim);
  printf("Saisir une valeur > 0 de sigma: ");
  scanf("%ld", &sigma);
  slong len = sigma;
  slong shifts[rdim];
  slong  first_shifts[rdim], second_shifts[rdim], third_shifts[rdim], fourth_shifts[rdim],
    fifth_shifts[rdim], sixth_shifts[rdim];

  
  // Creates random poly mat
  nmod_poly_mat_init(mat, rdim, cdim, prime);
  
  flint_rand_t state;
  
  flint_randinit(state);
  
  nmod_poly_mat_randtest(mat, state, len);
  printf("\nStudied matrix \n");
  printf("m = %ld\nn = %ld \nF = GF(%ld)\nsigma = %ld\n", rdim, cdim, prime, sigma);
  printf("A = Ms(");
  nmod_poly_mat_print_sage(mat);
  printf(")\n\n");
  // Creates random shifts
  //_perm_randtest(shifts, rdim, state);
  for (slong i = 0; i < rdim; i++)
    shifts[i] = 0;
  printf("shifts = ");
  int64_print_sage(shifts, rdim);

  printf("\n");
  timeit_t t0;
  
  nmod_poly_mat_init(first_res, rdim, rdim, prime);

  timeit_start(t0);
  M_basis(first_res, first_shifts, mat, sigma, shifts);
  timeit_stop(t0);
  flint_printf("M_basis: cpu = %wd ms  wall = %wd ms\n", t0->cpu, t0->wall);

  /** print
  printf("MBasis matrix representation\n");
  nmod_poly_mat_print_sage(first_res);
  printf("shifts result");
  int64_print_sage(first_shifts, rdim);
  */
  if (sigma < 2)
    {
      nmod_mat_init(A, rdim, cdim, prime);
      coefficient_matrix(A, mat, 0);
      nmod_poly_mat_init(second_res, rdim, rdim, prime);
      Basis(second_res, second_shifts, A, shifts);
      printf("\nBasis()\n");
      nmod_poly_mat_print_sage(second_res);
      printf("res_shifts basis\n");
      int64_print_sage(second_shifts, rdim);
    }

  nmod_poly_mat_init(third_res, rdim, rdim, prime);

  timeit_start(t0);
  M_basisII(third_res, third_shifts, mat, sigma, shifts);
  timeit_stop(t0);
  flint_printf("M_basisII cpu = %wd ms  wall = %wd ms\n", t0->cpu, t0->wall);

  /** print
  printf("MBasis II  polynomial representation for F_prime and convert P_k into list poly mat\n");
  nmod_poly_mat_print_sage(third_res);
  printf("shifts result");
  int64_print_sage(third_shifts, rdim);
  */
  
  nmod_poly_mat_init(fourth_res, rdim, rdim, prime);

  timeit_start(t0);
  M_basisIII(fourth_res, fourth_shifts, mat, sigma, shifts);
  timeit_stop(t0);
  flint_printf("M_basisIII cpu = %wd ms  wall = %wd ms\n", t0->cpu, t0->wall);

  /**
  printf("MBasis polynomial representation for F_prime recalculated for each loop and P_k poly mat\n");
  nmod_poly_mat_print_sage(fourth_res);
  printf("shifts result");
  int64_print_sage(fourth_shifts, rdim);
  */
  nmod_poly_mat_init(fifth_res, rdim, rdim, prime);

  timeit_start(t0);
  M_basisIV(fifth_res, fifth_shifts, mat, sigma, shifts);
  timeit_stop(t0);
  flint_printf("M_basisV cpu = %wd ms  wall = %wd ms\n", t0->cpu, t0->wall);

  /**
  printf("MBasis polynomial representation for F_prime and P_k recalculated for each loop \n");
  nmod_poly_mat_print_sage(fifth_res);
  printf("shifts result");
  int64_print_sage(fifth_shifts, rdim);
  */
  
  nmod_poly_mat_init(sixth_res, rdim, rdim, prime);

  timeit_start(t0);
  M_basisV(sixth_res, sixth_shifts, mat, sigma, shifts);
  timeit_stop(t0);
  flint_printf("M_basisV cpu = %wd ms  wall = %wd ms\n", t0->cpu, t0->wall);

  /**
  printf("MBasis polynomial representation for F_prime and P_k recalculated for each loop \n");
  nmod_poly_mat_print_sage(sixth_res);
  printf("shifts result");
  int64_print_sage(sixth_shifts, rdim);
  */
  if (!nmod_poly_mat_equal(first_res, third_res))
    printf("ERROR not equals m_basis I and II\n");

  if (!nmod_poly_mat_equal(fourth_res, third_res))
    printf("ERROR not equals m_basis II and III\n");

  if (!nmod_poly_mat_equal(fourth_res, first_res))
    printf("ERROR not equals m_basis I and III\n");

  if (!nmod_poly_mat_equal(fifth_res, first_res))
    printf("ERROR not equals m_basis I and IV\n");
  
  if (!nmod_poly_mat_equal(sixth_res, first_res))
    printf("ERROR not equals m_basis I and V\n");

  if (!nmod_poly_mat_equal(sixth_res, fifth_res))
   printf("ERROR not equals m_basis IV and V\n");
 
  return 0;
}

/** Test 
 * Verify if nmod_poly_mat_mul + coefficient_matrix and 
 * nmod_list_poly_mat_naive_mul_coef 
 * return the same result.
 * Verify if nmod_poly_mat -> nmod_list_poly_mat -> nmod_poly_mat
 * unchanged the initial matrix
 */
int test_nmod_list_poly_mat(void)
{
  nmod_mat_t prod, prod_2;
  nmod_poly_mat_t mat, mat_2, res, res_2, poly_prod;
  nmod_list_poly_mat_t mat_repr, mat_repr_2; 
  slong rdim = 6, cdim = 2, prime = 7,
    len = 10, len_2 = 20, k = 28;  
  char *x = malloc(1);

  /** init **/
  nmod_poly_mat_init(mat, rdim, cdim, prime);
  nmod_poly_mat_init(mat_2, cdim, rdim, prime);

  nmod_poly_mat_init(res, rdim, cdim, prime);
  nmod_poly_mat_init(res_2, cdim, rdim, prime);
  
  nmod_poly_mat_init(poly_prod, rdim, rdim, prime);

  nmod_mat_init(prod, rdim, rdim, prime);
  nmod_mat_init(prod_2, rdim, rdim, prime);
  
  flint_rand_t state;
  flint_randinit(state);
  srand(time(NULL));
  flint_randseed(state, rand(), rand());

  /** test poly_mat <-> list_poly_mat **/
  nmod_poly_mat_randtest(mat, state, len);
  nmod_poly_mat_randtest(mat_2, state, len_2);

  nmod_list_poly_mat_init_set(mat_repr, mat, 10);
  nmod_list_poly_mat_init_set(mat_repr_2, mat_2, 9);

  nmod_list_poly_mat_to_poly_mat(res, mat_repr);
  nmod_list_poly_mat_to_poly_mat(res_2, mat_repr_2);
  
  if (!nmod_poly_mat_equal(mat, res))
    printf("error mat not equals");
  if (!nmod_poly_mat_equal(mat_2, res_2))
    printf("\nerror mat_2 not equals\n");

  /** test mul: poly_mat <-> list_poly_mat **/
  nmod_poly_mat_mul(poly_prod, res, res_2);
  coefficient_matrix(prod, poly_prod, k);

  nmod_poly_mat_print(poly_prod, x);

  nmod_list_poly_mat_naive_mul_coef(prod_2, mat_repr, mat_repr_2, k);


  //print
  nmod_mat_print_pretty(prod);
  nmod_mat_print_pretty(prod_2);
  
  if (!nmod_mat_equal(prod_2, prod))
    printf("error mul not equals");

  /** clear **/
  nmod_list_poly_mat_clear(mat_repr);
  nmod_poly_mat_clear(mat);
  nmod_poly_mat_clear(res);
  
  nmod_list_poly_mat_clear(mat_repr_2);
  nmod_poly_mat_clear(mat_2);
  nmod_poly_mat_clear(res_2);

  nmod_poly_mat_clear(poly_prod);
  nmod_mat_clear(prod);
  nmod_mat_clear(prod_2);
  
  flint_randclear(state);
  free(x);
  
  return 0;
}

/** Test 
 * Verify if structured_multiplication blocks and 
 * list_structured_multiplication_blocks and basis + nmod_poly_mat_mul
 * return the same result
 */
int test_structured_mul_blocks(void)
{
  nmod_poly_mat_t mat, res, res_2, resid, resid_2, resid_3;
  nmod_list_poly_mat_t mat_repr;
  nmod_mat_t A, constant_mat;
  slong rank, rdim = 8, cdim = 4, prime = 101, len = 10, sigma = 30,
    shifts[rdim], res_shifts[rdim],
    *perm = _perm_init(rdim);

  /** init **/
  nmod_poly_mat_init(mat, rdim, cdim, prime);
  
  nmod_poly_mat_init(resid, rdim, cdim, prime);
  nmod_poly_mat_init(resid_2, rdim, cdim, prime);
  nmod_poly_mat_init(resid_3, rdim, cdim, prime);

  nmod_poly_mat_init(res, rdim, rdim, prime);
  nmod_poly_mat_init(res_2, rdim, rdim, prime);

  flint_rand_t state;
  flint_randinit(state);
  srand(time(NULL));
  flint_randseed(state, rand(), rand());
  
  nmod_poly_mat_randtest(mat, state, len);
  nmod_poly_mat_set(resid, mat);
  nmod_poly_mat_set(resid_2, mat);
  nmod_list_poly_mat_init_set(mat_repr, mat, sigma);

  
  for (slong i = 0; i < rdim; i++)
    shifts[i] =  rand() % 10 - 5;
  
  nmod_mat_init(constant_mat, rdim, cdim, prime);
 
  coefficient_matrix(constant_mat, mat, 0);
  rank = Basis_for_M_basis(A, res_shifts, perm, constant_mat, shifts);
  Basis(res_2, res_shifts, constant_mat, shifts);
  nmod_poly_mat_mul(resid_2, res_2, resid_2);

  nmod_poly_mat_one(res);
  structured_multiplication_blocks(res, A, perm, rank);

  structured_multiplication_blocks(resid, A, perm, rank);
   
  structured_list_multiplication_blocks(mat_repr, A, perm, rank, sigma);

  nmod_list_poly_mat_to_poly_mat(resid_3, mat_repr); 
  
  // print
  char *x = malloc(1);
  nmod_poly_mat_print(res, x); 
  nmod_poly_mat_print(res_2, x); 
  nmod_poly_mat_print(resid, x); 
  nmod_poly_mat_print(resid_2, x); 
  nmod_poly_mat_print(resid_3, x);
  free(x);

  // test
  if (!nmod_poly_mat_equal(res, res_2))
    printf("Basis and Basis for M Basis not same result");

  if (!nmod_poly_mat_equal(resid, resid_2))
    printf("Basis and Basis for M Basis not same result");

  if (!nmod_poly_mat_equal(resid, resid_3))
    printf("Basis and list_struct_mul not same result");
  
  /** clear **/
  nmod_mat_clear(A);
  nmod_mat_clear(constant_mat);

  nmod_poly_mat_clear(mat);
  nmod_poly_mat_clear(res);
  nmod_poly_mat_clear(resid);
  nmod_poly_mat_clear(res_2);
  nmod_poly_mat_clear(resid_2);
  nmod_poly_mat_clear(resid_3);
  
  nmod_list_poly_mat_clear(mat_repr);

  _perm_clear(perm);
  flint_randclear(state);
  return 1;
}


/** Test
 * verify if there are leaks in basis.c functions
 **/
int test_basis(void)
{
  nmod_mat_t mat;
  slong rdim = 9, cdim = 3, prime = 7;

  nmod_mat_init(mat, rdim, cdim, prime);

   
  flint_rand_t state;
  flint_randinit(state);
  srand(time(NULL));
  flint_randseed(state, rand(), rand());
  
  nmod_mat_randtest(mat, state);

  
  nmod_poly_mat_t res; 
  nmod_poly_mat_init(res, rdim, rdim, prime);
  nmod_poly_mat_zero(res);

  int64_t res_shift[rdim];

  int64_t shift[rdim];
  
  for (slong i = 0; i < rdim; i++)
    shift[i] =  rand() % 10 - 5;
  
  Basis(res, res_shift, mat, shift);


  printf("Matrix\n");
  nmod_mat_print_sage(mat);
  printf("Result Basis\n");
  nmod_poly_mat_print_sage(res);
  nmod_poly_mat_clear(res);
  flint_randclear(state);

  nmod_mat_t res2;
  int64_t res_shift2[rdim];
  slong res_perm2[rdim];
  slong rank = Basis_for_M_basis(res2, res_shift2, res_perm2, mat, shift);
    
  nmod_mat_clear(mat);
  
  printf("\nMatrix A from basis = [[x,0],[A,1]]\n");
  nmod_mat_print_pretty(res2);
  
  printf("\nshifts = ");
  for(slong i = 0; i < rdim; i++)
    printf("%ld ", res_shift2[i]);

  printf("\npermutation = ");
  for(slong i = 0; i < rdim; i++)
    printf("%ld ", res_perm2[i]);
  

  printf("\nrank = %ld\n", rank);
  
  nmod_mat_clear(res2);
  return 0;
}


/** Test
 * many function of matpol, verify if there are leaks 
 **/
int test_matpol(void)
{
    nmod_poly_mat_t A;
    nmod_poly_mat_init(A,3,3,4);

    slong cdim = A->c, rdim = A->r;
    printf("number of columns: %lu, and rows: %lu\n", cdim, rdim);

    flint_rand_t seed;
    flint_randinit(seed);
    srand(time(NULL));
    flint_randseed(seed, rand(), rand());
    
   
    nmod_poly_mat_randtest(A, seed, 5);

    printf("A\n");
    nmod_poly_mat_print_sage(A);

    
    
    nmod_mat_t B;
    nmod_mat_init(B, rdim, cdim, 4);
    coefficient_matrix(B, A, 2);
    printf("\ncoefficient matrix for degree 2 of A\n");
    nmod_mat_print_sage(B);

    int64_t shifts[cdim];

    shifts[0] = 1; shifts[1] = 2; shifts[2] = 3;
    printf("shifts: ");
    int64_print_sage(shifts, cdim);

    int64_t cols_deg[cdim];
    column_degrees(cols_deg, A, shifts);
    int64_t rows_deg[rdim];
    row_degrees(rows_deg, A, shifts);

    printf("\ncolumn degree of A with shift \n");
    int64_print_sage(cols_deg, rdim);
    printf("\nrow degree of A with shift \n");
    int64_print_sage(rows_deg, cdim);
    
    slong deg_A = nmod_poly_mat_degree(A);
    printf("\nA's degree: %ld\n", deg_A);

    
    matrix_wise row_wise = 0;
    int64_t *mat_deg = malloc(sizeof(int64_t) * rdim * cdim);
 
    degree_matrix(mat_deg, A, shifts, row_wise);
    printf("\n");
    int64_mat_print(mat_deg, rdim, cdim);
    
    int64_t lead_pos[cdim];
    leading_positions(lead_pos, A, shifts, row_wise);
    printf("leading position: ");
    int64_print_sage(lead_pos, cdim);

    leading_matrix(B, A, shifts, row_wise);
    printf("leading matrix for shift shifts of A\n");
    nmod_mat_print_sage(B);

    nmod_poly_mat_clear(A);
    nmod_mat_clear(B);
    flint_randclear(seed);
    free(mat_deg);
    
    return 1;
}


int main(void)
{
  //test_nmod_list_poly_mat();
  test_mbasis();
  //test_matpol();
  //test_basis();
  //test_structured_mul_blocks();
  return EXIT_SUCCESS;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
