#include <matpol.h>

void int64_print(const int64_t *shifts, slong length)
{
  printf("[");
  for (slong i = 0; i < length - 1; i++)
    printf("%ld,", shifts[i]);
  printf("%ld]\n", shifts[length - 1]);
}

void nmod_poly_mat_print_pretty(const nmod_poly_mat_t mat,
				slong rdim, slong cdim)
{
  nmod_poly_t P;
  slong length;
  mp_limb_t mod = nmod_poly_mat_modulus(mat);

  nmod_poly_init(P, mod);
  printf("[");
  for (slong i = 0; i < rdim; i++)
    {
      printf("[");
      for (slong j = 0; j < cdim; j++)
      {
	nmod_poly_set(P, nmod_poly_mat_entry(mat, i, j));
	length = nmod_poly_length(P);
	if (length == 0)
	  {
	    if (j != cdim - 1)
	      printf("0, ");
	    else
	      printf("0");
	  }
	else
	  {
	    for (slong k = 0; k < length; k++)
	      {
		if (k != length - 1) 
		  {
		    if (k == 0)
		      printf("%ld +", nmod_poly_get_coeff_ui(P, k));
		    else
		      printf("%ld*x^%ld + ", nmod_poly_get_coeff_ui(P, k), k);
		  }
		else
		  {
		    if (j != cdim - 1)
		      {
		
			if (k == 0)
			  printf("%ld,", nmod_poly_get_coeff_ui(P, k));
			else
			  printf("%ld*x^%ld,", nmod_poly_get_coeff_ui(P, k), k);	

		      }
		    else
		      {
			
			if (k == 0)
			  printf("%ld", nmod_poly_get_coeff_ui(P, k));
			else
			  printf("%ld*x^%ld", nmod_poly_get_coeff_ui(P, k), k);
		      }
		  }
	      }
	  }
      }
      if (i != rdim -1)
	printf("],\n");
      else
	printf("]");
    }
  printf("]\n");
}

void nmod_mat_print(const nmod_mat_t mat,
		    slong rdim, slong cdim)
{
  printf("[");
  for (slong i = 0; i < rdim; i++)
    {
      printf("[");
      for (slong j = 0; j < cdim; j++)
      {
	if (j != cdim - 1)
	  printf("%ld, ",  nmod_mat_get_entry(mat, i, j));
	else
	  printf("%ld",  nmod_mat_get_entry(mat, i, j));

      }
      if (i != rdim -1)
	printf("],\n");
      else
	printf("]");
    }
  printf("]\n");
}


void coefficient_matrix(nmod_mat_t res, const nmod_poly_mat_t mat, int degree)
{
    slong rows = nmod_mat_nrows(res), cols = nmod_mat_ncols(res);
    nmod_poly_struct * P;
    for(slong i = 0; i < rows; i++)
        for(slong j = 0; j < cols; j++)
        {
            P = nmod_poly_mat_entry(mat, i, j);
            nmod_mat_entry(res, i, j) = nmod_poly_get_coeff_ui(P, degree);
        }
}

void column_degrees(uint64_t *res, const nmod_poly_mat_t mat, uint64_t *shifts)
{
    // test len(res) == cols and len(shifts) == rows
    slong cols = nmod_poly_mat_ncols(mat), rows = nmod_poly_mat_nrows(mat);
    nmod_poly_struct *P;
    slong max, d;
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
        res[i] = (uint64_t) max;
    }
}

void row_degrees(uint64_t *res, const nmod_poly_mat_t mat, uint64_t *shifts)
{
    // test len(res) == rows and len(shifts) == cols 
    slong cols = nmod_poly_mat_ncols(mat), rows = nmod_poly_mat_nrows(mat);
    nmod_poly_struct *P;
    slong max, d;
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
        res[i] = (uint64_t) max;
    }
}

slong nmod_poly_mat_degree(const nmod_poly_mat_t mat)
{
    slong cols = nmod_poly_mat_ncols(mat), rows = nmod_poly_mat_nrows(mat);
    nmod_poly_struct *P;
    slong d, deg = 0;
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


void degree_matrix(uint64_t *res, const nmod_poly_mat_t mat, uint64_t *shifts,
                   matrix_wise row_wise)
{
    //test len(res) == rows and len(res[0]) == cols
    slong cols = nmod_poly_mat_ncols(mat), rows = nmod_poly_mat_nrows(mat);
    nmod_poly_struct *P;
    slong d;

    if (row_wise)
    {
        for(slong i = 0; i < rows; i++)
        {
            for(slong j = 0; j < cols; j++)
            {

                P = nmod_poly_mat_entry(mat, i, j);
                d = nmod_poly_degree(P) + shifts[j];
                *(res + (i * rows) + j) = d;
            }
        }
        return;
    }
    for(slong i = 0; i < cols; i++)
        for(slong j = 0; j < rows; j++)
        {
            P = nmod_poly_mat_entry(mat, j, i);
            d = nmod_poly_degree(P) + shifts[j];
            *(res + i  + (j * cols)) = d;
        }
}


void leading_matrix(nmod_mat_t res, const nmod_poly_mat_t mat,
                    uint64_t *shifts, matrix_wise row_wise)
{
    slong cols = nmod_poly_mat_ncols(mat), rows = nmod_poly_mat_nrows(mat);
    nmod_poly_t P;
    if (row_wise)
    {
        uint64_t rdeg[rows];
        row_degrees(rdeg, mat, shifts);
        for(slong i = 0; i < rows; i++)
            for(slong j = 0; j < cols; j++)
            {
	      nmod_poly_set(P, nmod_poly_mat_entry(mat, i, j));
	      nmod_mat_set_entry(res, i, j, nmod_poly_get_coeff_ui(P, rdeg[i] - shifts[j]));
            }
        return;
    }
    {
        uint64_t cdeg[cols];
        column_degrees(cdeg, mat, shifts);
        for(slong i = 0; i < cols; i++)
            for(slong j = 0; j < rows; j++)
            {
	      nmod_poly_set(P, nmod_poly_mat_entry(mat, j, i));
	      nmod_mat_set_entry(res, j, i, nmod_poly_get_coeff_ui(P, cdeg[i] - shifts[j]));
            }
    }
}

void leading_positions(uint64_t *res, const nmod_poly_mat_t mat,
                       uint64_t *shifts, matrix_wise row_wise)
{
    // test len(res) == rows if row_wise == ROW_WISE or len(res) == COLUMN_WISE otherwise
    slong cols = nmod_poly_mat_ncols(mat), rows = nmod_poly_mat_nrows(mat);
    nmod_poly_struct *P;
    slong max;
    slong d;
    uint64_t ind; 
    if (row_wise)
    {
        for (slong i = 0; i < rows; i++)
        {
            max = 0;
            for (slong j = 0; j < cols; j++)
            {
                P = nmod_poly_mat_entry(mat, i, j);
                d = nmod_poly_degree(P) + shifts[j];
                if (max < d)
                {
                    max = d;
                    ind = j;
                }
            }
            res[i] = ind;
        }
        return;
    }
    {
        for (slong i = 0; i < cols; i++)
        {
            max = 0;
            for (slong j = 0; j < rows; j++)
            {
                P = nmod_poly_mat_entry(mat, j, i);
                d = nmod_poly_degree(P) + shifts[j];
                if (max < d)
                {
                    max = d;
                    ind = j;
                }
            }
            res[i] = ind;
        }
    }
}

int is_hermite(const nmod_poly_mat_t mat, matrix_wise row_wise)
{
    slong cols = nmod_poly_mat_ncols(mat), rows = nmod_poly_mat_nrows(mat);
    slong deg_mat = nmod_poly_mat_degree(mat);

    if (row_wise)
    {
        uint64_t shifts[cols];
        for (slong i = 0; i < cols; i++)
            shifts[i] = (cols - i) * (deg_mat + 1);
        return is_popov(mat, shifts, row_wise, 0);
    }

    uint64_t shifts[rows];
    for (slong i = 0; i < rows; i++)
        shifts[i] = i * (deg_mat + 1);
    return is_popov(mat, shifts, row_wise, 0);

}

int is_popov(const nmod_poly_mat_t mat, uint64_t *shifts, matrix_wise row_wise, int ordered)
{
    if (!is_weak_popov(mat, shifts, row_wise, ordered))
        return 0;
    slong cols = nmod_poly_mat_ncols(mat), rows = nmod_poly_mat_nrows(mat), pivot_deg, d;
    nmod_poly_struct *P, *pivot;

    if (row_wise)
    {
        uint64_t lead_pos[rows];
        leading_positions(lead_pos, mat, shifts, row_wise);

        for(slong i = 0; i < rows; i++)
        {
            pivot = nmod_poly_mat_entry(mat, i, lead_pos[i]);
	    pivot_deg = nmod_poly_degree(pivot);
            if (nmod_poly_get_coeff_ui(pivot, nmod_poly_degree(pivot)) != 1)
                return 0;
            for(slong j = 0; j < rows ; j++)
            {
                P = nmod_poly_mat_entry(mat, j, lead_pos[i]);
                d = nmod_poly_degree(P);
                if (d >= pivot_deg)
                    return 0;
            }
        }
        return 1;
    }

    uint64_t lead_pos[cols];
    leading_positions(lead_pos, mat, shifts, row_wise);

    for(slong i = 0; i < cols; i++)
    {
        pivot = nmod_poly_mat_entry(mat, lead_pos[i], i);
	pivot_deg = nmod_poly_degree(pivot);
	if (nmod_poly_get_coeff_ui(pivot, nmod_poly_degree(pivot)) != 1)
            return 0;
        for(slong j = 0; j < cols ; j++)
        {
            P = nmod_poly_mat_entry(mat, lead_pos[i], j);
            d = nmod_poly_degree(P);
            if (d >= pivot_deg)
                return 0;
        }
    }
    return 1;

}

int is_reduced(const nmod_poly_mat_t mat, uint64_t *shifts, matrix_wise row_wise)
{
    slong cols = nmod_poly_mat_ncols(mat), rows = nmod_poly_mat_nrows(mat);
    nmod_mat_t B;
    nmod_mat_init(B, rows, cols, nmod_poly_mat_modulus(mat));
    leading_matrix(B, mat, shifts, row_wise);
    slong rank_lead = nmod_mat_rank(B);
    return (int) (rows == rank_lead);
}

static int intComparator ( const void * first, const void * second ) {
    int firstInt = * (const int *) first;
    int secondInt = * (const int *) second;
    return firstInt - secondInt;
}


int is_weak_popov(const nmod_poly_mat_t mat, uint64_t *shifts, matrix_wise row_wise, int ordered)
{
  if (!is_reduced(mat, shifts, row_wise))
    return 0;
  slong cols = nmod_poly_mat_ncols(mat), rows = nmod_poly_mat_nrows(mat);
  
  if (row_wise)
    {
      uint64_t lead_pos[rows];
      leading_positions(lead_pos, mat, shifts, row_wise);
      
      if (!ordered)
	qsort(lead_pos, rows, sizeof(uint64_t), intComparator);
      
      for (slong i = 0; i < rows - 1; i++)
        {
	  if (lead_pos[i] > lead_pos[i+1])
	    return 0;
        }
      return 1;
    }

  uint64_t lead_pos[cols];
  leading_positions(lead_pos, mat, shifts, row_wise);

  if (!ordered)
    {
      qsort(lead_pos, cols, sizeof(uint64_t), intComparator);
	
      for (slong i = 0; i < cols; i++)
	{
	  if (lead_pos[i] > lead_pos[i+1])
	    return 0;
	}
    }
  return 1;
}

int is_zero_mod_xk(const nmod_poly_mat_t mat, int64_t k)
{
  nmod_poly_t P;
  nmod_poly_init(P, mat->modulus);

  for(slong i = 0; i < mat->r; i++)
    for(slong j = 0; j < mat->c; j++)
      {
	nmod_poly_set(P, nmod_poly_mat_entry(mat, i, j));
	nmod_poly_shift_right(P, P, k+1);
	if (!nmod_poly_is_zero(P))
	  {
	    nmod_poly_clear(P);
	    return 0;
	  }
      }
  nmod_poly_clear(P);
  return 1;
}

int is_minimal_approximant_basis(const nmod_poly_mat_t base,
				 const nmod_mat_t mat, int64_t order,
				 int64_t *shifts)
{
  slong rdim = mat->r, cdim = mat->c;
  mp_limb_t prime = mat->mod.n;
  nmod_poly_t constant;
  nmod_poly_mat_t mat_poly, res_mul;

  if (base->c != base->r)
    {
      printf("not basis: wrong shape");
      return 0;
    }
  
  nmod_poly_mat_init(mat_poly, rdim, cdim, prime);
  slong alloc;
  nmod_poly_init(constant, prime); 
  for (slong i = 0; i < rdim; i++)
    for (slong j = 0; j < cdim; j++)
      {
	alloc = (slong) nmod_mat_get_entry(mat, i, j);
	nmod_poly_set_coeff_ui(constant, 0, alloc);
	nmod_poly_set(nmod_poly_mat_entry(mat_poly, i, j), constant); 
      }
  nmod_poly_mat_init(res_mul, rdim, cdim, prime);  
  nmod_poly_mat_mul(res_mul, base, mat_poly);
  
  if (! is_zero_mod_xk(res_mul, order))
    {
      printf("not zero");
      return 0;
    }
  uint64_t lead_pos[rdim];
  leading_positions(lead_pos, base, shifts, ROW_WISE);
  printf("\nleading positions\n");
  for (slong i = 0; i < rdim; i++)
    printf("%lu ", lead_pos[i]); 
  return 1;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
