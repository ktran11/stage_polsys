#include "pm_basis.h"
#include "matpol.h"
#include <flint.h>


/* Type only for the function Basis */
typedef struct
{
  int64_t value;
  slong ord;
} int_tuple;

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


/**
 * \brief A function fund in the doc of nmod_mat_t
 *
 */
static void nmod_poly_mat_swap_cols(nmod_poly_mat_t mat, slong * perm, slong r, slong s)
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
        {
            mp_limb_t c = mat->rows[t][r];
            mat->rows[t][r] = mat->rows[t][s];
            mat->rows[t][s] = c;
        }
    }
}

/**
 * \brief Do the multiplication permutation matrix (P*A) for the function Basis
 *
 * \param A, the input and result matrix
 * \param P, the permutation 
 * \param n, length of P and row dimension of A
*/
static void apply_perm_to_matrix(nmod_poly_mat_t A, slong *P, slong n)
{
  for (slong i = 0; i < n; i++)
    nmod_mat_swap_cols(A, P, 0, P[0]); 
}

void Basis(nmod_poly_mat_t res, nmod_mat_t mat, int64_t *shift,
	   slong rdim, slong cdim, slong prime)
{

  slong rank, *P, *P_inv, *perm_inv, *comp, *perm = _perm_init(rdim);
  nmod_mat_t G, Lr, Lr_inv, prod_G_Lr_inv;
  nmod_poly_t One, constant;

  sort_and_create_perm(perm, shift, rdim);

  /* LU operation on pi*mat and extraction  L = [ [Lr, 0], [G, I] ] */
  apply_perm_to_matrix(mat, perm, rdim);
  rank = nmod_mat_lu(P, mat, 0);
  
  nmod_mat_init(Lr, rank, rank, prime);
  nmod_mat_one(Lr);
  for(slong i = 1; i < rank; i++)
    for(slong j = 0; j < i; j++)
      nmod_mat_set_entry(Lr, i, j, nmod_mat_entry(mat, i, j));

  nmod_mat_window_init(G, mat, rank, 0, rdim, cdim);
  

  /* Computation of -G*Lr^(-1) */
  nmod_mat_init(Lr_inv, rank, rank, prime);
  nmod_mat_inv(Lr_inv, Lr); // inv can be improved with op in tril

  nmod_mat_clear(Lr);

  nmod_mat_init(prod_G_Lr_inv, rdim - rank, rank, prime);
  nmod_mat_mul(prod_G_Lr_inv, G, Lr_inv); // can be improved with choice of algorithm (default = Strassen) and multithread
  
  nmod_mat_neg(prod_G_L_inv, G_Lr_inv);
  
  nmod_mat_clear(Lr_inv);
  nmod_mat_clear(G);

  /* Computation of the block matrix  [ [xIr, 0], [-G*Lr^(-1), I_{rdim - rank}] ] */
  // We suppose res set at [0]_{rdim x rdim}

  nmod_poly_init(One, prime);
  nmod_poly_set_coeff_ui(One, 0, 1); // 1
  
  for (slong i = rank; i < rdim; i++)
    nmod_poly_mat_entry(res, i, i) = One;

  nmod_poly_shift_left(One, One, 1); // x
  for (slong i = 0; i < rank; i++)
    nmod_poly_mat_entry(res, i, i) = One;

  nmod_poly_clear(One);

  for (slong i = 0; i < rdim - rank; i++)
    for (slong j = 0; j < rank; j++)
      {
	nmod_poly_init2(constant, ); // TO DO //
	nmod_poly_mat_entry(res, i + rank, j) = constant; 
      }

  nmod_mat_clear(G_Lr_inv);
  /* Multiply by the permutations */

  P_inv = _perm_init(rdim);
  perm_inv = _perm_init(rdim);
  comp =  _perm_init(rdim);
  
  _perm_inv(perm_inv, perm, rdim);
  _perm_inv(P_inv, P, rdim);
  
  _perm_compose(comp, P_inv, perm, rdim);

  apply_perm_to_matrix(res, comp, rdim);

  /* Compute shift */

  
    
  _perm_clear(perm);
  _perm_clear(P);
  _perm_clear(perm_inv);
  _perm_clear(P_inv);
}   

int main(void)
{
  nmod_mat_t mat, G, L;
  slong rdim = 8, cdim = 4, prime = 7, rank;
  slong *P = _perm_init(rdim);

  
  flint_rand_t state;
  flint_randinit(state);

  nmod_mat_init(mat, rdim, cdim, prime);
  nmod_mat_randtest(mat, state);
  nmod_mat_print_pretty(mat);

  /** LU test 
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
  _perm_print(p, rdim);

  apply_perm_to_matrix(mat, p, rdim);
  nmod_mat_print_pretty(mat);

  int64_t s[] = { 4, 1, -2, 3 };
  sort_and_create_perm(p, s, rdim);
  
  _perm_print(p, rdim);
  _perm_clear(p);  
  */
  
  nmod_mat_clear(mat);
  flint_randclear(state);
  return EXIT_SUCCESS;
}
