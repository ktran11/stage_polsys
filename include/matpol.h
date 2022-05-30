#ifndef MATPOL_H
#define MATPOL_H

#include <stdbool.h>
#include <nmod_mat.h>
#include <nmod_poly_mat.h>
#include <fmpz.h>
#include <fmpz_mat.h>
#include <flint.h>


void coefficient_matrix(nmod_mat_t res, const nmod_poly_mat_t mat, int degree);


fmpz *column_degrees(nmod_poly_mat_t mat, fmpz *shifts);

fmpz *row_degrees(nmod_poly_mat_t mat, fmpz *shifts);

int degree(nmod_poly_mat_t mat);

fmpz_mat_t *degree_matrix(nmod_poly_mat_t mat, fmpz *shifts, bool row_wise);

bool is_constant(nmod_poly_mat_t mat);

nmod_poly_mat_t *leading_matrix(nmod_poly_mat_t mat, fmpz *shifts, bool row_wise);

fmpz *leading_positions(nmod_poly_mat_t mat, fmpz *shifts, bool row_wise);

nmod_poly_mat_t *reverse(nmod_poly_mat_t mat, bool row_wise);

nmod_poly_mat_t *shift(nmod_poly_mat_t mat, bool row_wise);

nmod_poly_mat_t *truncate(nmod_poly_mat_t mat, int degree, bool row_wise);


#endif /* MATPOL_H */
