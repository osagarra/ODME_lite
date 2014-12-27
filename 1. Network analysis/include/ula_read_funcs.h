#ifndef ULA_READ_FUNCS_H
#define ULA_READ_FUNCS_H
	
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_fit.h>
#include "ula_gen_funcs.h"

/******* Read funcs ********/
double** read_distances(char *input_name,int num_nodes, int header, int log_opt);
int**    read_node_list_int(char *input_name,int num_nodes, int header);
int* read_node_list_int_undir(char *input_name,int num_nodes, int header);
double** read_node_list_double(char *input_name,int num_nodes, int header);
double* read_node_list_double_undir(char *input_name,int num_nodes, int header);
int**    read_net_list(char *input_name, int num_nodes, int header);
double** read_net_list_double(char *input_name, int num_nodes, int header);
//void write_indices(double i1, double i11, double i2, double i22, double i3, double i33, double gamma,int reps, char *output);
/******* Print funcs ******/
//void     print_w_2_file_double(char *input_name, int N_nodes, double **w, double **w2, double eps);
//void     print_w_2_file_int(char *input_name, int N_nodes, int **w);
/******* Index funcs ******/
/*
double* z_score_theor(int **w_real, double **x, int N_nodes, double eps, double wsum);
double* z_score_theor_rad(int **w_real, double **x, int* s, int N_nodes, double eps, double wsum);
double* z_score_trips(int **w_real, double **w_sim, double **w2_sim,  int N_nodes, int *missing, double eps);
double* rel_error_trips(int **w_real, double **w_sim, int N_nodes, int E);
double   sorensen_index(int **w_real, int **w_sim, int N_nodes);
double   nmae_index(int **w_real, int **w_sim, int N_nodes);
double   nrmse_index(int **w_real, int **w_sim, int N_nodes);
double   pnzt_index(int **w_real, int **w_sim, int N_nodes);
*/
/******* Net funcs from adj ******/
/*
int ** compute_s(int ** w, int N_nodes);
int ** compute_k(int ** w, int N_nodes);
double * compute_w_double(int** w, int N_nodes, int* num_edges, int zeros);
int * compute_w_int(int** w, int N_nodes, int* num_edges, int zeros);
double * compute_d_edges(int** w, double ** d, int N_nodes, int* num_edges);
double * compute_d_trips(int** w, double ** d, int N_nodes, int* num_trips);
gsl_histogram * compute_pij(int ** w, double ** d, int bins, int N_nodes, double max_d_edge);
int compute_sij(int  **s, double ** dist, int origin, int dest, int N_nodes, gsl_rng * randgsl, int perturb);
double *compute_sij_edges(int ** w, double **d, int N_nodes, gsl_rng * randgsl);
double *compute_s_out_edges(int ** w, int N_nodes );
double *compute_s_in_edges(int ** w, int N_nodes);
double *compute_s_out_in_edges(int ** w, int N_nodes);
double** s_nn(int** w , int N_nodes);
*/
/******* Accumulators ******/
/*
void acc_update_all(gsl_histogram** accs, int** w_false, double ** d, int N_nodes, int* E2, int* E1, int* T1);
void acc2d_update_all(gsl_histogram2d ** accs, int ** w, double ** dist, int N_nodes, gsl_rng * randgsl);
*/
/******* Tests ******/
/*
void test_k_s(int **w, int N_nodes);
void test_w(int **w, int N_nodes);
void test_d(int **w, double ** d, int N_nodes, int bins);
void test_pij(int ** w_real, double ** d, int N_nodes, int bins);
*/

#endif
