#ifndef OTHERS_NULL_MODELS_H
#define OTHERS_NULL_MODELS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>

#include "ula_gen_funcs.h"
#include "ula_read_funcs.h"
#include "ula_graph_structs.h"
#include "ula_w_graph_funcs.h"

double compute_gamma_frenchies(double surface);
//double compute_p_french(int *s_in, double **d, double gamma, int dest, int origin, int* ind_in , int * ind_out, int len, int self_opt, int verbose);
void* compute_p_french_multinomial(double* p, int* s_in, double **d, double gamma, int real_origin, int* ind_in ,int len, int N_nodes, int self_opt, int verbose);
double compute_p_french_bernouilli(int* s_in, double **d, double gamma, int origin, int dest, int* ind_out, int* ind_in ,int len, int self_opt);
//double compute_p_french_ini(int* s_in, double **d, double gamma, int dest, int origin, int len, int self_opt, int verbose);
//void check_norm_french(int *s_in, double **d, double gamma, int origin, int* ind_in, int *ind_out, int len, int self_opt, int verbose);
//W_GRAPH* w_graph_seq_gravity_directed(int N_nodes,int **s, double **d, double gamma, gsl_rng* randgsl, int Max_fails, int self_opt, int verbose, int ps_opt, double**ps);
W_GRAPH* w_graph_seq_gravity_directed(int N_nodes,int **s, double **d, double gamma, gsl_rng* randgsl,int self_opt, int verbose);

/***********/
int destination_rad_model(int** s, double** d, int origin, int N_nodes, gsl_rng * randgsl);
int compute_sij_rad(int  **s, double **dist, int origin, int dest, int N_nodes);
W_GRAPH*  w_graph_radiation_model_stochastic_directed(int N_nodes,int **s, double **d, gsl_rng* randgsl, int self_opt, int verbose);
W_GRAPH* w_graph_radiation_model_multinomial_directed(int N_nodes,int **s, double **d, gsl_rng* randgsl, int self_opt, int verbose, int ps_opt, double**ps);
/***********/
W_GRAPH* gravity_poisson_directed_graph2(double**x,  int N_nodes , double ** d, double gamma, gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* gravity_poisson_undirected_graph2(double*x,  int N_nodes , double **d, double gamma, gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* w_graph_seq_gravity_multinomial_directed(int N_nodes,int **s, double **d, double gamma, gsl_rng* randgsl, int self_opt, int verbose);
W_GRAPH* w_graph_seq_gravity_bernouilli_directed(int N_nodes,int **s, double **d, double gamma, gsl_rng* randgsl, int self_opt, int verbose);
/***** Prob Choice ******/
double * prob_mult_Cs_undir(int N_nodes, double* x, double gamma, double** d);
double * prob_mult_Cs_dir(int N_nodes, double **x, double gamma, double** d);
double * prob_mult_C_undir(int N_nodes, double gamma, double** d);
double * prob_mult_C_dir(int N_nodes, double gamma, double** d);

#endif
