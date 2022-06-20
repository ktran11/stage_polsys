#include "main.h"
#include "nmod_list_poly_mat.h"
#include "m_basis.h"
#include "matpol.h"
#include "basis.h"
#include "perm_operator.h"
#include <time.h>
#include <profiler.h>

/**
int test_mul_special_matrix(void)
{
  //nmod_mat_init(A, rdim - rank, rank, prime);
  //nmod_mat_randtest(A, state);

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
    
}
*/

int test_mbasis(void)
{
  nmod_poly_mat_t mat, first_res, second_res, third_res, fourth_res;
  nmod_mat_t A;
  slong rdim, cdim, prime, sigma;
  
  printf("Saisir rdim: ");
  scanf("%ld", &rdim);
  printf("Saisir cdim: ");
  scanf("%ld", &cdim);
  printf("Saisir prime: ");
  scanf("%ld", &prime);
  printf("Saisir une valeur > 0 de sigma: ");
  scanf("%ld", &sigma);
  slong len = sigma;
  slong shifts[rdim];
  slong  first_shifts[rdim], second_shifts[rdim], third_shifts[rdim], fourth_shifts[rdim];

  
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
  _perm_randtest(shifts, rdim, state);
  printf("shifts = ");
  int64_print_sage(shifts, rdim);


  // timeit_t t0;
  
  nmod_poly_mat_init(first_res, rdim, rdim, prime);

  // timeit_start(t0);
  M_basis(first_res, first_shifts, mat, sigma, shifts);
  //timeit_stop(t0);
  //flint_printf("cpu = %wd ms  wall = %wd ms\n", t0->cpu, t0->wall);

  printf("MBasis matrix representation\n");
  nmod_poly_mat_print_sage(first_res);
  printf("shifts result");
  int64_print_sage(first_shifts, rdim);

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

  // timeit_start(t0);
  M_basisII(third_res, third_shifts, mat, sigma, shifts);
  //timeit_stop(t0);
  //flint_printf("cpu = %wd ms  wall = %wd ms\n", t0->cpu, t0->wall);

  printf("MBasis II  polynomial representation for F_prime and convert P_k into list poly mat\n");
  nmod_poly_mat_print_sage(third_res);
  printf("shifts result");
  int64_print_sage(third_shifts, rdim);

  
  nmod_poly_mat_init(fourth_res, rdim, rdim, prime);

  // timeit_start(t0);
  M_basisIII(fourth_res, fourth_shifts, mat, sigma, shifts);
  //timeit_stop(t0);
  //flint_printf("cpu = %wd ms  wall = %wd ms\n", t0->cpu, t0->wall);

  printf("MBasis polynomial representation for F_prime recalculated for each loop and P_k poly mat\n");
  nmod_poly_mat_print_sage(fourth_res);
  printf("shifts result");
  int64_print_sage(fourth_shifts, rdim);
  
  if (!nmod_poly_mat_equal(first_res, third_res))
    printf("ERROR not equals m_basis I and II");

  if (!nmod_poly_mat_equal(fourth_res, third_res))
    printf("ERROR not equals m_basis II and III");

  if (!nmod_poly_mat_equal(fourth_res, first_res))
    printf("ERROR not equals m_basis I and III");

  return 0;
}

int test_nmod_list_poly_mat(void)
{
  nmod_poly_mat_t mat;
  nmod_list_poly_mat_t mat_repr; 
  slong rdim = 6, cdim = 2, prime = 7, len = 10;  
  char *x = malloc(1);

  // Creates random poly mat
  nmod_poly_mat_init(mat, rdim, cdim, prime);
  
  flint_rand_t state;
  flint_randinit(state);

  nmod_poly_mat_randtest(mat, state, len);
  printf("\nStudied matrix\n");
  nmod_poly_mat_print(mat, x);

  nmod_list_poly_mat_init_set(mat_repr, mat, 0);
  nmod_list_poly_mat_print(mat_repr);

  nmod_list_poly_mat_clear(mat_repr);
  nmod_poly_mat_clear(mat);

  return 0;
}


int test_left_kernel(void)
{
 
  nmod_mat_t mat;
  slong *perm, rdim = 16, cdim = 8, prime = 3;
  nmod_poly_mat_t res;
  int64_t res_shifts[rdim], shifts[rdim];
  
  nmod_mat_init(mat, rdim, cdim, prime);

  flint_rand_t state;
  flint_randinit(state);
  srand(time(NULL));
  flint_randseed(state, rand(), rand());

  nmod_mat_randtest(mat, state);

  for (slong i = 0; i < rdim; i++)
    shifts[i] =  rand() % 10 - 5;

  printf("Matrix\n");
  nmod_mat_print_sage(mat);

  printf("Shifts");
  int64_print_sage(shifts, rdim);

  printf("\n");
  perm = _perm_init(rdim);

  nmod_poly_mat_init(res, rdim, rdim, prime);

  printf("Basis\n");
  Basis(res, res_shifts, mat, shifts);

  
  //nmod_poly_mat_print_pretty(res, rdim, rdim);
  //printf("New shifts ");
  //int64_print(res_shifts, rdim);


  nmod_mat_t res_II;
  int64_t res_II_shifts[rdim];
  slong *res_perm_II;
  res_perm_II = _perm_init(rdim);

  Basis_for_M_basis(res_II, res_II_shifts, res_perm_II, mat, shifts);
  printf("\n Result Basis for M Basis \n");
  printf("Mat\n");
  nmod_mat_print_sage(res_II);
  printf("\nshifts");
  int64_print_sage(res_II_shifts, rdim);
  printf("\n perm ");
  _perm_print(res_perm_II, rdim);


  nmod_mat_clear(res_II);
  _perm_clear(res_perm_II);  
  nmod_poly_mat_clear(res);
  _perm_clear(perm);
  nmod_mat_clear(mat);
  flint_randclear(state);
  return 0;
  
}



int test_basis(void)
{
  nmod_mat_t mat, copy_mat;
  slong rdim = 9, cdim = 3, prime = 7;

  nmod_mat_init(mat, rdim, cdim, prime);

   
  flint_rand_t state;
  flint_randinit(state);

  nmod_mat_randtest(mat, state);
  nmod_mat_init_set(copy_mat, mat);
  
  nmod_poly_mat_t res; 
  nmod_poly_mat_init(res, rdim, rdim, prime);
  nmod_poly_mat_zero(res);

  int64_t res_shift[rdim];

  int64_t shift[rdim];

  for (int64_t i = 0; i < rdim; i++)
    shift[i] = i;
  
  Basis(res, res_shift, mat, shift);
  
  nmod_poly_mat_clear(res);

  nmod_mat_print_pretty(copy_mat);
  
  nmod_mat_t res2;
  int64_t res_shift2[rdim];
  slong res_perm2[rdim];
  slong rank = Basis_for_M_basis(res2, res_shift2, res_perm2, copy_mat, shift);

  nmod_mat_print_pretty(res2);

  nmod_mat_clear(res2);

  nmod_mat_randtest(mat,state);

  printf("\nnew_mat\n");
  nmod_mat_print_pretty(mat);
  rank = Basis_for_M_basis(res2, res_shift2, res_perm2, copy_mat, shift);

  printf("new_res2\n");
  nmod_mat_print_pretty(res2);
  
  printf("\nshifts = ");
  for(slong i = 0; i < rdim; i++)
    printf("%ld ", res_shift2[i]);

  printf("\npermutation = ");
  for(slong i = 0; i < rdim; i++)
    printf("%ld ", res_perm2[i]);
  

  printf("\nrank = %ld\n", rank);
  nmod_mat_clear(res2);
  
  /**
  rank = nmod_mat_lu(P, mat, 0);
  nmod_mat_print_pretty(mat);

  nmod_mat_init(L, rank, rank, prime);
  nmod_mat_one(L);
  for(slong i = 1; i < rank; i++)
    for(slong j = 0; j < i; j++)
      nmod_mat_set_entry(L, i, j, nmod_mat_entry(mat, i, j));
  
  nmod_mat_print_pretty(L);

  printf("rank mat = %ld", rank);

  nmod_mat_window_init(G, mat, rank, 0, rdim, cdim);

  nmod_mat_print_pretty(G);
  */

  /** petit test sort and create perm
  slong *p = _perm_init(rdim);
  _perm_randtest(p, rdim, state);
  printf("\nperm \n");
  _perm_print(p, rdim);

  nmod_mat_print_pretty(mat);
  apply_perm_rows_to_matrix(mat, p, rdim);
  nmod_mat_print_pretty(mat);

  printf("perm after\n");
  _perm_print(p, rdim);

  printf("\n");
  int64_t s[] = { 4, 1, -2, 3, 5, 1, 2, 5 };
  sort_and_create_perm(p, s, rdim);
  
  _perm_print(p, rdim);
  _perm_clear(p);
  */
  nmod_mat_clear(mat);
  flint_randclear(state);
  return 0;
}


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
  //test_left_kernel();
  test_mbasis();
  //test_matpol();
  return EXIT_SUCCESS;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
