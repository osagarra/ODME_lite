#ifndef ULA_GEN_FUNCS_H
#define ULA_GEN_FUNCS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

/****** File handling ********/
FILE* open_file(char *mode,char *nom);

int* safe_int_realloc(int* vect, int len, int new_len, int zero_val);
void* xrealloc(void* vector, size_t midaenbytes);
/****** Number handling ********/
int mineq_int(int a, int b);
int maxeq_int(int a, int b);
/****** Vector handling ********/
double* cast_vec_double(int len);
int* cast_vec_int(int len);
void reverse_vec_int(int* inver_a,int j);
void reverse_vec_double(double* inver_a,int j);
int count_appearances_int(int *vec, int val, int len);
int count_appearances_double(double *vec, double val, double tol, int len);
int find_value_int(int * vec, int val, int len);
double sum_vec_double(double* vec, int len);
int sum_vec_int(int* vec, int len);
int max_value_int(int* vect, int len);
int min_value_int(int* vect, int len);
double min_value_double(double* vect, int len);
double max_value_double(double* vect, int len);
int* vec_double_to_int(double* vect, int len);
double * vec_int_to_double(int* vect, int len);
void scale_vec_int(int* vect, int factor, int len);
void scale_vec_double(double* vect, double factor, int len);
/****** Matrix handling ********/
void free_mat_int(int **mat_d, int rows);
void free_mat_double(double **mat_d, int rows);
double** cast_mat_double(int rows, int cols);
int** cast_mat_int(int rows, int cols);
double * flatten_matrix_triangular_double(double** w, int rows, int * length);
int * flatten_sparse_matrix_int(int** w, int rows, int cols, int zeros, int * length);
double * flatten_sparse_matrix_double(double** w, int rows, int cols, int zeros, int * length);
int * flatten_matrix_int(int** w, int rows, int cols, int diag, int * length);
double * flatten_matrix_double(double** w, int rows, int cols, int diag, int * length);
double matrix_max_value_double(double** vect, int row, int col);
double matrix_min_value_double(double** vect, int row, int col);
int matrix_max_value_int(int** vect, int row, int col);
int matrix_min_value_int(int** vect, int row, int col);
void average_matrix(double** matrix, int rows, int cols, int reps);
void scale_matrix(double** matrix, int rows, int cols, double reps);
int sum_matrix_int(int** matrix, int rows, int cols);
double sum_matrix_double(double** matrix, int rows, int cols);

#endif
