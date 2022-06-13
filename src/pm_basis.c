#include "pm_basis.h"
#include "matpol.h"
#include <flint.h>


/* Type only for the function Basis */
typedef struct
{
  int64_t value;
  slong ord;
} int_tuple;

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
    {
      perm[temp[i].ord] = i;
    }
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
  apply_perm_to_vector(res_shifts, shifts, comp, rdim);

  _perm_clear(comp);
  
  for (slong i = 0; i < rank; i++)
    res_shifts[i] += 1;

  apply_perm_to_vector(res_shifts, res_shifts, comp_inv, rdim);

  _perm_clear(comp_inv);
}


slong Basis_for_M_basis(nmod_mat_t res, int64_t *res_shifts, slong *res_perm,
		       const nmod_mat_t mat, const int64_t *shifts,
			slong rdim, slong cdim, slong prime)
{
  if (res != NULL)
    nmod_mat_clear(res);

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
  apply_perm_to_vector(res_shifts, shifts, res_perm, rdim);
  
  for (slong i = 0; i < rank; i++)
    res_shifts[i] += 1;

  comp_inv = _perm_init(rdim);
  _perm_inv(comp_inv, res_perm, rdim);
  apply_perm_to_vector(res_shifts, res_shifts, comp_inv, rdim);

  _perm_clear(comp_inv);
  return rank;
}



int main(void)
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
  
  printf("Matrix M \n");
  nmod_mat_print_pretty(mat);
  
  Basis(res, res_shift, mat, shift, rdim, cdim, prime);

  char *x = malloc(150);
  printf("\nBasis result\n");
  nmod_poly_mat_print(res, x);

  printf("shift res\n");
  for (slong i = 0; i <rdim; i++)
    printf("%ld ", res_shift[i]);
  
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
