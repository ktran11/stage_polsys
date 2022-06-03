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
static void reverse(int arr[], int n)
{
    for (int low = 0, high = n - 1; low < high; low++, high--)
    {
        int temp = arr[low];
        arr[low] = arr[high];
        arr[high] = temp;
    }
}
**/

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

int main(void)
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

    qsort(lead_pos, rows, sizeof(uint64_t), intComparator);
    printf("leading position (sorted): ");
    print_vect_uint64(lead_pos, cols);

    leading_matrix(B, A, shifts, row_wise);
    printf("leading matrix for shift shifts of A");
    nmod_mat_print_pretty(B);

    return EXIT_SUCCESS;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
