#ifndef MATPOL_H
#define MATPOL_H

#include <stdbool.h>


nmod_mat_t coefficient_matrix(nmod_poly_mat_t mat, int degree, bool row_wise);

fpmz *column_degrees(nmod_poly_mat_t mat, fpmz *shifts);

fpmz *row_degrees(nmod_poly_mat_t mat, fpmz *shifts);

int degree(nmod_poly_mat_t mat);

fmpz_mat_t degree_matrix(nmod_poly_mat_t mat, fpmz *shifts, bool row_wise);

nmod_poly_mat_t is_constant(nmod_poly_mat_t mat);

nmod_poly_mat_t leading_matrix(nmod_poly_mat_t mat, fpmz *shifts, bool row_wise);

fpmz *leading_positions(nmod_poly_mat_t mat, fpmz *shifts, bool row_wise, bool return_degree);

nmod_poly_mat_t reverse(nmod_poly_mat_t mat, bool row_wise);

nmod_poly_mat_t shift(nmod_poly_mat_t mat, bool row_wise);

nmod_poly_mat_t truncate(nmod_poly_mat_t mat, int degree, bool row_wise);


#endif /* MATPOL_H */
