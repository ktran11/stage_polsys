#include "main.h"
#include "m_basis.h"
#include "matpol.h"
#include "basis.h"

#include <profiler.h>

static void print_vect_uint64(const uint64_t *vect, slong len)
{
    if (!vect)
        return;
    for (slong i = 0; i < len; i++)
        printf("%lu ", *(vect+i));
    printf("\n");
    return;
}

static void print_mat_uint64(const uint64_t *mat, slong rows, slong cols)
{
    if (!mat)
        return;
    for(slong i = 0; i < rows; i++)
    {
        for(slong j = 0; j < cols; j++)
            printf("%ld ", *(mat+(i*rows)+j)) ;
        printf("\n" );
    }
    return;
}

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
  nmod_poly_mat_t mat, copy_mat, first_res, second_res;
  nmod_mat_t A;
  slong rdim = 6, cdim = 3, prime = 7, rank = 3, len = 10, *shifts,
    first_shifts[rdim], second_shifts[rdim];  
  char *x = malloc(150);

  
  nmod_poly_mat_init(mat, rdim, cdim, prime);
  
  flint_rand_t state;
  flint_randinit(state);

  nmod_poly_mat_randtest(mat, state, len);
  nmod_poly_mat_init_set(copy_mat, mat);

  
  shifts= _perm_init(rdim);
  _perm_randtest(shifts, rdim, state);

  nmod_poly_mat_init(first_res, rdim, rdim, prime);
  printf("shifts: \n");
  print_vect_uint64(shifts, rdim);
  
  M_basis(first_res, first_shifts, mat, 1, shifts, rdim, cdim, prime);

  
  printf("M_basis(sigma = 1)\n");
  nmod_poly_mat_print(first_res, x);
  printf("res_shifts\n");
  print_vect_uint64(first_shifts, rdim);
  
  nmod_mat_init(A, rdim, cdim, prime);
  
  coefficient_matrix(A, mat, 0);
  nmod_poly_mat_init(second_res, rdim, rdim, prime);

  Basis(second_res, second_shifts, A, shifts, rdim, cdim, prime);

  printf("\nBasis()\n");
  nmod_poly_mat_print(second_res,x);
  print_vect_uint64(second_shifts, rdim);
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
  
  Basis(res, res_shift, mat, shift, rdim, cdim, prime);
  
  nmod_poly_mat_clear(res);

  nmod_mat_print_pretty(copy_mat);
  
  nmod_mat_t res2;
  int64_t res_shift2[rdim];
  slong res_perm2[rdim];
  slong rank = Basis_for_M_basis(res2, res_shift2, res_perm2, copy_mat, shift, rdim, cdim, prime);

  nmod_mat_print_pretty(res2);

  nmod_mat_clear(res2);

  nmod_mat_randtest(mat,state);

  printf("\nnew_mat\n");
  nmod_mat_print_pretty(mat);
  rank = Basis_for_M_basis(res2, res_shift2, res_perm2, copy_mat, shift, rdim, cdim, prime);

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

    slong cols = nmod_poly_mat_ncols(A), rows = nmod_poly_mat_nrows(A);
    printf("number of columns: %lu, and rows: %lu\n", cols,rows);

    flint_rand_t seed;
    flint_randinit(seed);
    nmod_poly_mat_randtest(A, seed, 5);
    char *x = malloc(150);

    printf("A");
    nmod_poly_mat_print(A, x);

    nmod_mat_t B;
    nmod_mat_init(B, rows, cols, 4);
    coefficient_matrix(B, A, 2);
    printf("coefficient matrix for degree 2 of A");
    nmod_mat_print_pretty(B);

    uint64_t shifts[cols];

    shifts[0] = 1; shifts[1] = 2; shifts[2] = 3;
    printf("shifts: ");
    print_vect_uint64(shifts, cols);

    uint64_t cols_deg[cols];
    column_degrees(cols_deg, A, shifts);
    uint64_t rows_deg[rows];
    row_degrees(rows_deg, A, shifts);

    printf("column degree of A with shift \n");
    print_vect_uint64(cols_deg, rows);
    printf("row degree of A with shift \n");
    print_vect_uint64(rows_deg, cols);

    slong deg_A = nmod_poly_mat_degree(A);
    printf("A's degree: %ld\n", deg_A);

    matrix_wise row_wise = 0;
    uint64_t *mat_deg = malloc(sizeof(uint64_t) * cols * rows);
    degree_matrix(mat_deg, A, shifts, row_wise);
    printf("\n");
    print_mat_uint64(mat_deg, rows, cols);

    uint64_t lead_pos[cols];
    leading_positions(lead_pos, A, shifts, row_wise);
    printf("leading position: ");
    print_vect_uint64(lead_pos, cols);

    leading_matrix(B, A, shifts, row_wise);
    printf("leading matrix for shift shifts of A");
    nmod_mat_print_pretty(B);

    return EXIT_SUCCESS;
}


int main(int argc, char *argv)
{
  test_mbasis();
  return EXIT_SUCCESS;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
