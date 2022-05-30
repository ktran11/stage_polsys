#include <fmpz.h>
#include <matpol.h>


void coefficient_matrix(nmod_mat_t res, const nmod_poly_mat_t mat, int degree)
{
  slong rows = nmod_mat_nrows(res), cols = nmod_mat_ncols(res);
  nmod_poly_t *P;
  for(slong i = 0; i < rows; i++)
    for(slong j = 0; j < cols; j++)
      {
	P = nmod_poly_mat_entry(mat, i, j);
	nmod_mat_entry(res, i, j) = nmod_poly_get_coeff_ui(P, degree);
      }
}


fmpz *column_degrees(nmod_poly_mat_t mat, fmpz *shifts)
{
  slong cols = nmod_poly_mat_ncols(mat), rows = nmod_poly_mat_nrows(mat);
  fmpz *cols_deg = _fmpz_vec_init(cols);
  nmod_poly_t *P;
  fmpz max;
  slong d;
  for (slong i = 0; i < cols; i++)
    {
      max = 0;
      for (slong j = 0; j < rows; j++)
	{
	  P = nmod_poly_mat_entry(mat, j, i);
	  d = nmod_poly_degree(P) + shifts[j];
	  if (max < d)
	    max = d;
	}
      cols_deg[i] = max;
    }
  return cols_deg;
}

fmpz *row_degrees(nmod_poly_mat_t mat, fmpz *shifts)
{
  slong cols = nmod_poly_mat_ncols(mat), rows = nmod_poly_mat_nrows(mat);
  fmpz *rows_deg = _fmpz_vec_init(cols);
  nmod_poly_t *P;
  fmpz max;
  slong d;
  for (slong i = 0; i < rows; i++)
    {
      max = 0;
      for (slong j = 0; j < cols; j++)
	{
	  P = nmod_poly_mat_entry(mat, i, j);
	  d = nmod_poly_degree(P) + shifts[j];
	  if (max < d)
	    max = d;
	}
      rows_deg[i] = max;
    }
  return rows_deg;
}

fmpz nmod_mat_pol_degree(nmod_poly_mat_t mat)
{
  slong cols = nmod_poly_mat_ncols(mat), rows = nmod_poly_mat_nrows(mat);
  nmod_poly_t *P;
  fmpz d, deg = 0;
  for(slong i = 0; i < rows; i++)
    for(slong j = 0; j < cols; j++)
      {
	P = nmod_poly_mat_entry(mat, i, j);
	d = nmod_poly_degree(P);
	if (deg < d)
	  deg = d;
	
      }
  return deg;
}

fmpz_mat_t *degree_matrix(nmod_poly_mat_t mat, fmpz *shifts, bool row_wise)
{
  slong cols = nmod_poly_mat_ncols(mat), rows = nmod_poly_mat_nrows(mat);
  fmpz_mat_t D;
  fmpz_mat_init(D, rows, cols);
  nmod_poly_t *P;
  fmpz d;

  if (row_wise)
    {
      for(slong i = 0; i < rows; i++)
	for(slong j = 0; j < cols; j++)
	  {
	    P = nmod_poly_mat_entry(mat, i, j);
	    d = nmod_poly_degree(P) + shifts[j];
	    D[i,j] = d
	  }
      return D;
    }
  
  for(slong i = 0; i < cols; i++)
    for(slong j = 0; j < rows; j++)
      {
	P = nmod_poly_mat_entry(mat, j, i);
	d = nmod_poly_degree(P) + shifts[j];
	D[i,j] = d;
      }
  return D;
}

bool is_constant(nmod_poly_mat_t mat)
{
  slong cols = nmod_poly_mat_ncols(mat), rows = nmod_poly_mat_nrows(mat);
  nmod_poly_t *P;
  fmpz d;
  for(slong i = 0; i < rows; i++)
    for(slong j = 0; j < cols; j++)
      {
	P = nmod_poly_mat_entry(mat, i, j);
	d = nmod_poly_degree(P);
	if (d > 0)
	  return false;
      }
  return true;
 
}

nmod_mat_t *leading_matrix(nmod_poly_mat_t mat, fmpz *shifts, bool row_wise)
{
  
}

fmpz *leading_positions(nmod_poly_mat_t mat, fmpz *shifts, bool row_wise)
{
  slong cols = nmod_poly_mat_ncols(mat), rows = nmod_poly_mat_nrows(mat);
  nmod_poly_t *P;
  fmpz max;
  slong d;
  fmpz ind; 
  if (row_wise)
    {
      fmpz *ind_lead = _fmpz_vec_init(rows);

      for (slong i = 0; i < rows; i++)
	{
	max = 0;
	for (slong j = 0; j < cols; j++)
	  {
	    P = nmod_poly_mat_entry(mat, i, j);
	    d = nmod_poly_degree(P) + shifts[i];
	    if (max < d)
	      {
		max = d;
		ind = j;
	      }
	  }
	ind_lead[i] = ind;
      }
    return ind_lead;
      
    }
  {
    fmpz *ind_lead = _fmpz_vec_init(cols);
    for (slong i = 0; i < cols; i++)
      {
	max = 0;
	for (slong j = 0; j < rows; j++)
	  {
	    P = nmod_poly_mat_entry(mat, j, i);
	    d = nmod_poly_degree(P) + shifts[i];
	    if (max < d)
	      {
		max = d;
		ind = j;
	      }
	  }
	ind_lead[i] = ind;
      }
    return ind_lead;
  }
}


nmod_poly_mat_t *reverse(nmod_poly_mat_t mat, bool row_wise);

nmod_poly_mat_t *shift(nmod_poly_mat_t mat, bool row_wise);

nmod_poly_mat_t *truncate(nmod_poly_mat_t mat, int degree, bool row_wise);

int main(void)
{
  nmod_poly_mat_t A;
  nmod_poly_mat_init(A,3,3,4);
  
  slong cols = nmod_poly_mat_ncols(A), rows = nmod_poly_mat_nrows(A);
  printf("%lu, %lu\n", cols,rows);

  flint_rand_t seed;
  flint_randinit(seed);
  nmod_poly_mat_randtest(A, seed, 5);
  char *x = malloc(150);
  nmod_poly_mat_print(A,x);

  nmod_mat_t B;
  nmod_mat_init(B, 3, 3, 4);
  coefficient_matrix(B, A, 2);

  nmod_mat_print_pretty(B);
  
  fmpz *shifts = _fmpz_vec_init(3);

  shifts[0] = 1; shifts[1] = 2; shifts[2] = 3;
  
  fmpz *cols_deg = column_degrees(A, shifts);
  fmpz *rows_deg = row_degrees(A, shifts);
  _fmpz_vec_print(cols_deg, 3);
  printf("\n");
  _fmpz_vec_print(rows_deg, 3);
  printf("\n");

  fmpz deg_A = nmod_mat_pol_degree(A);
  printf("%ld\n", deg_A);

  fmpz *lead_pos = leading_positions(A, shifts, true);
  _fmpz_vec_print(lead_pos, 3);
  printf("\n");
  return EXIT_SUCCESS;
}
