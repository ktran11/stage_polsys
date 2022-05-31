#ifndef MATPOL_H
#define MATPOL_H

#include <stdint.h>
#include <nmod_mat.h>
#include <nmod_poly_mat.h>


typedef enum
  {
    COLUMN_WISE = 0,
    ROW_WISE = 1

  } matrix_wise; 



void coefficient_matrix(nmod_mat_t res, const nmod_poly_mat_t mat, int degree);

void column_degrees(uint64_t *res, const nmod_poly_mat_t mat, uint64_t *shifts);

void row_degrees(uint64_t *res, const nmod_poly_mat_t mat, uint64_t *shifts);

slong nmod_poly_mat_degree(const nmod_poly_mat_t mat);

void degree_matrix(uint64_t *res, const nmod_poly_mat_t mat, uint64_t *shifts,
			  matrix_wise row_wise);

int is_hermite(const nmod_poly_mat_t mat, matrix_wise row_wise);

int is_popov(const nmod_poly_mat_t mat, uint64_t *shifts, matrix_wise row_wise, int ordered);

int is_reduced(const nmod_poly_mat_t mat, uint64_t *shifts, matrix_wise row_wise);

int is_weak_popov(const nmod_poly_mat_t mat, uint64_t *shifts, matrix_wise row_wise, int ordered);

void leading_matrix(nmod_mat_t res, const nmod_poly_mat_t mat, uint64_t *shifts,
				matrix_wise row_wise);

void leading_positions(uint64_t *res, const nmod_poly_mat_t mat, uint64_t *shifts,
			matrix_wise row_wise);

#endif /* MATPOL_H */
