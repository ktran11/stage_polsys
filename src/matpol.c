#include <matpol.h>

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
    nmod_poly_struct *P;
    if (row_wise)
    {
        uint64_t rdeg[rows];
        row_degrees(rdeg, mat, shifts);
        for(slong i = 0; i < rows; i++)
            for(slong j = 0; j < cols; j++)
            {
                P = nmod_poly_mat_entry(mat, i, j);
                nmod_mat_entry(res, i, j) = nmod_poly_get_coeff_ui(P, rdeg[i] - shifts[j]);
            }
        return;
    }
    {
        uint64_t cdeg[cols];
        column_degrees(cdeg, mat, shifts);
        for(slong i = 0; i < cols; i++)
            for(slong j = 0; j < rows; j++)
            {
                P = nmod_poly_mat_entry(mat, j, i);
                nmod_mat_entry(res, j, i) = nmod_poly_get_coeff_ui(P, cdeg[i] - shifts[j]);
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
        qsort(lead_pos, cols, sizeof(uint64_t), intComparator);

    for (slong i = 0; i < cols; i++)
    {
        if (lead_pos[i] > lead_pos[i+1])
            return 0;
    }
    return 1;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
