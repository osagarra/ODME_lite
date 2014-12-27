#ifndef ULA_NULL_MODELS_H
#define ULA_NULL_MODELS_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_roots.h>

#include "ula_gen_funcs.h"
#include "ula_read_funcs.h"
#include "ula_graph_structs.h"
#include "ula_w_graph_funcs.h"


/***********/
W_GRAPH* uncorrelated_multinomial_directed_graph(double* ps, double**x, int N_nodes , int T, gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* uncorrelated_multinomial_undirected_graph(double* ps, double*x, int N_nodes , int T, gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* uncorrelated_computational_directed_graph(int** s, int N_nodes , gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* uncorrelated_computational_undirected_graph(int* s, int N_nodes , gsl_rng* randgsl, int max_trials, int verbose, int self_opt);
W_GRAPH* uncorrelated_poisson_multinomial_directed_graph(double* ps, double**x, int N_nodes , int T, gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* uncorrelated_poisson_multinomial_undirected_graph(double* ps, double*x,  int N_nodes , int T, gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* uncorrelated_poisson_undirected_graph2(double*x,  int N_nodes , gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* uncorrelated_poisson_directed_graph2(double**x,  int N_nodes , gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* fixedEs_poisson_undirected_graph(double*x,  double lam, int N_nodes , gsl_rng* randgsl, int verbose, int self_opt, int max_reps);
W_GRAPH* fixedEs_poisson_directed_graph(double**x,  double lam, int N_nodes , gsl_rng* randgsl, int verbose, int self_opt, int max_reps);
W_GRAPH* custompij_poisson_directed_graph(double**pij,  int N_nodes , gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* custompij_poisson_undirected_graph(double**pij,  int N_nodes , gsl_rng* randgsl, int verbose, int self_opt);
/*********Distros of p's*******/
double * prob_mult_s_undir(double* x, int N_nodes, int self_opt);
double * prob_mult_s_dir(double** x, int N_nodes, int self_opt);
/********* Aux funcs *******/
int compute_T(int N_nodes, double av_k, int* x);

#endif
