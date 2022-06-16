/**
 * \file mat_pol.h
 * \brief shorts functions
 * \author Kevin Tran
 * \version 1.0
 * \date June 6 2022
 *
 * Add some functions availables on SageMath for the class nmod_poly_mat_t on flint
 * https://doc.sagemath.org/html/en/reference/matrices/sage/matrix/matrix_polynomial_dense.html
 * All parameters are supposed init
 *
 */

#ifndef MATPOL_H
#define MATPOL_H

#include <stdint.h>
#include <nmod_mat.h>
#include <nmod_poly_mat.h>


/**
 * \enum matrix_wise
 * \brief way to see a polynomial matrix 
 *
 * Permits the user to use function in different perspectives 
 *
 */
typedef enum
{
    COLUMN_WISE = 0,
    ROW_WISE = 1

} matrix_wise;

void int64_print(const int64_t *shifts, slong length);

int is_minimal_approximant_basis(const nmod_poly_mat_t base,
				 const nmod_mat_t mat, int64_t order,
				 int64_t *shifts);
void nmod_mat_print(const nmod_mat_t mat,
		    slong rdim, slong cdim);

void nmod_poly_mat_print_pretty(const nmod_poly_mat_t mat,
				slong rdim, slong cdim);

/**
 * \fn void coefficient_matrix(nmod_mat_t res, const nmod_poly_mat_t mat, int degree)
 * \brief set on res the coefficient matrix of degree degree of mat
 * 
 * if mat = sum_{i = 0}^{k} m_i X^i, then res = m_{degree}
 *
 * \param res a matrix on the same ring as mat and has the same dimensions (empty)
 * \param mat a polynomial matrix 
 * \param degree the degree the user need to study 
 *
 */
void coefficient_matrix(nmod_mat_t res, const nmod_poly_mat_t mat, int degree);

/**
 * \fn void column_degrees(uint64_t *res, const nmod_poly_mat_t mat, uint64_t *shifts)
 * \brief set on res the columns degree of mat with the shifts shifts
 * 
 * if mat = (m_{i,j}) then res = (max_i(m_{i,j}.degree + shifts[i]))_j
 *
 * \param res a pointer of an array of length the number of columns of mat (empty)
 * \param mat a polynomial matrix 
 * \param shifts a pointer of an array of lenght the number of rows of mat   
 *
 */
void column_degrees(uint64_t *res, const nmod_poly_mat_t mat, uint64_t *shifts);

/**
 * \fn void row_degrees(uint64_t *res, const nmod_poly_mat_t mat, uint64_t *shifts)
 * \brief Stock on res the rows degree of mat with the shifts shifts
 * 
 * if mat = (m_{i,j}) then res = (max_j(m_{i,j}.degree + shifts[j]))_i
 *
 * \param res a pointer of an array of length the number of columns of mat (empty)
 * \param mat a polynomial matrix 
 * \param shifts a pointer of an array of lenght the number of columns of mat   
 *
 */
void row_degrees(uint64_t *res, const nmod_poly_mat_t mat, uint64_t *shifts);

/**
 * \fn slong nmod_poly_mat_degree(const nmod_poly_mat_t mat)
 * \brief Gives the maximal polynomial degree for all the entry of mat
 * 
 * if mat = (m_{i,j}) then res = max(mat.column_degrees()) = max(mat.row_degrees())
 *
 * \param mat a polynomial matrix
 * \return the degree of the matrix mat
 */
slong nmod_poly_mat_degree(const nmod_poly_mat_t mat);


/**
 * \fn void degree_matrix(uint64_t *res, const nmod_poly_mat_t mat, uint64_t *shifts,
                   matrix_wise row_wise)
 * \brief Stock on res the degree matrix of mat with the shifts shifts see depending of the matrix_wise
 * 
 * if mat = (m_{i,j}) then res = ( degree(m_{i,j}) + shifts[i or j] ) 
 *
 *
 * \param res a pointer on the first value of the degree_matrix which being stock line by line (empty)
 * \param mat a polynomial matrix 
 * \param shifts a pointer of an array of lenght the number of rows/columns of mat depending of row_wise   
 * \param row_wise gives two options for the shifts
 */
void degree_matrix(uint64_t *res, const nmod_poly_mat_t mat, uint64_t *shifts,
                   matrix_wise row_wise);

int is_hermite(const nmod_poly_mat_t mat, matrix_wise row_wise);

/**
 * \fn int is_popov(const nmod_poly_mat_t mat, uint64_t *shifts, matrix_wise row_wise, int ordered);
 * \brief Verify if mat is  popov
 * 
 * Verify if mat is weak popov and the degree of each polynoms 
 * are strictly inferior the pivot of the columns or rows
 * 
 * \param mat a polynomial matrix 
 * \param shifts a pointer of an array of lenght the number of rows/columns of mat depending of row_wise   
 * \param row_wise gives two options for the shifts
 * \param ordered, 1 if the leading positions must be increasing, 0 else
 * \return 1 if mat is popov, 0 else
 */
int is_popov(const nmod_poly_mat_t mat, uint64_t *shifts, matrix_wise row_wise, int ordered);

/**
 * \fn  int is_reduced(const nmod_poly_mat_t mat, uint64_t *shifts, matrix_wise row_wise)
 * \brief Verify if mat is reduced
 * 
 * Verify if the leading matrix is full rank (e.g for any unimodular polynomial matrix U, rdeg_s(mat) <= rdeg_s(U*mat) if row_wise == ROW_WISE)   
 * 
 * \param mat a polynomial matrix 
 * \param shifts a pointer of an array of lenght the number of rows/columns of mat depending of row_wise   
 * \param row_wise gives two options for the shifts
 * \return 1 if mat is reduced, 0 else
 */
int is_reduced(const nmod_poly_mat_t mat, uint64_t *shifts, matrix_wise row_wise);

/**
 * \fn int is_weak_popov(const nmod_poly_mat_t mat, uint64_t *shifts, matrix_wise row_wise, int ordered)
 * \brief Verify if mat is weak popov
 * 
 * Verify if the leading positions of mat are distinct pairwise   
 * 
 * \param mat a polynomial matrix 
 * \param shifts a pointer of an array of lenght the number of rows/columns of mat depending of row_wise   
 * \param row_wise gives two options for the shifts
 * \param ordered, 1 if the leading positions must be increasing, 0 else
 * \return 1 if mat is weak popov, 0 else
 */
int is_weak_popov(const nmod_poly_mat_t mat, uint64_t *shifts, matrix_wise row_wise, int ordered);

/**
 * \fn void leading_matrix(nmod_mat_t res, const nmod_poly_mat_t mat, uint64_t *shifts,
                    matrix_wise row_wise)
 * \brief Stock on res the leading matrix of mat with the shifts shifts see depending of the matrix_wise
 * 
 * Will use the row or column degrees of mat to compute the leading matrix of mat,
 * For the i-st line, it will be stock the coefficient for the degree row/column_degrees[i] for the line.
 * 
 * \param res a nmod_mat_t where will be stock the result of dimensions mat (empty)
 * \param mat a polynomial matrix 
 * \param shifts a pointer of an array of lenght the number of rows/columns of mat depending of row_wise   
 * \param row_wise gives two options for the shifts
 */
void leading_matrix(nmod_mat_t res, const nmod_poly_mat_t mat, uint64_t *shifts,
                    matrix_wise row_wise);

/**
 * \fn void leading_positions(uint64_t *res, const nmod_poly_mat_t mat, uint64_t *shifts,
                       matrix_wise row_wise)
 * \brief Stock on res the leading positions of mat with the shifts shifts see depending of the matrix_wise
 * 
 * If we see mat with a ROW_WISE, the leading positions is the first column where the row degrees is reach starting by the right for each rows of mat.
 *
 * \param res a pointer with length the number of rows/columns of mat depending on row_wise (empty)
 * \param mat a polynomial matrix 
 * \param shifts a pointer of an array of lenght the number of columns/rows of mat depending of row_wise   
 * \param row_wise gives two options for the shifts
 */
void leading_positions(uint64_t *res, const nmod_poly_mat_t mat, uint64_t *shifts,
                       matrix_wise row_wise);

#endif /* MATPOL_H */

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
