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
#include <gsl/gsl_sf_lambert.h>

#include "ula_gen_funcs.h"
#include "ula_read_funcs.h"
#include "ula_graph_structs.h"
#include "ula_w_graph_funcs.h"



/*********** LInear constraints arbitrary p **************/
W_GRAPH* multinomial_directed_graph(double* ps, int N_nodes , int T, gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* multinomial_undirected_graph(double* ps, int N_nodes , int T, gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* poisson_multinomial_directed_graph(double* ps, int N_nodes , int T, gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* poisson_multinomial_undirected_graph(double* ps, int N_nodes , int T, gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* custompij_poisson_directed_graph(double**pij,  int N_nodes , gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* custompij_poisson_undirected_graph(double**pij,  int N_nodes , gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* custompij_geometric_directed_graph(double**pij,  int N_nodes , gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* custompij_geometric_undirected_graph(double**pij,  int N_nodes , gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* custompij_ZIP_undirected_graph(double**pij_b, double**pij,  int N_nodes , gsl_rng* randgsl, int verbose, int self_opt, int max_reps);
W_GRAPH* custompij_ZIP_directed_graph(double**pij_b, double**pij, int N_nodes , gsl_rng* randgsl, int verbose, int self_opt, int max_reps);
W_GRAPH* custompij_ZIG_undirected_graph(double**pij_b, double**pij,  int N_nodes , gsl_rng* randgsl, int verbose, int self_opt, int max_reps);
W_GRAPH* custompij_ZIG_directed_graph(double**pij_b, double**pij, int N_nodes , gsl_rng* randgsl, int verbose, int self_opt, int max_reps);
/***************** Micro-canonical configuration model ****************************/
W_GRAPH* fixeds_computational_directed_graph(int** s, int N_nodes , gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* fixeds_computational_undirected_graph(int* s, int N_nodes , gsl_rng* randgsl, int max_trials, int verbose, int self_opt);
/*********** Linear constraints **************/
//** Distinguishable weights **//
W_GRAPH* fixeds_poisson_undirected_graph2(double*x,  int N_nodes , gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* fixeds_poisson_directed_graph2(double**x,  int N_nodes , gsl_rng* randgsl, int verbose, int self_opt);
//** Undistinguishable weights **//
W_GRAPH* fixeds_geometric_undirected_graph2(double*x,  int N_nodes , gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* fixeds_geometric_directed_graph2(double**x,  int N_nodes , gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* fixeds_negbinomial_undirected_graph2(double*x,  int N_nodes , int layers, gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* fixeds_negbinomial_directed_graph2(double**x,  int N_nodes , int layers, gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* fixeds_binomial_undirected_graph2(double*x,  int N_nodes , int layers, gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* fixeds_binomial_directed_graph2(double**x,  int N_nodes , int layers, gsl_rng* randgsl, int verbose, int self_opt);

/*********** Binary constraints **************/
W_GRAPH* fixedk_poisson_undirected_graph(double*x,  double mu, int N_nodes , gsl_rng* randgsl, int verbose, int self_opt, int max_reps);
W_GRAPH* fixedk_poisson_directed_graph(double**x, double mu,  int N_nodes , gsl_rng* randgsl, int verbose, int self_opt, int max_reps);
W_GRAPH* fixedk_geometric_undirected_graph(double*x,  double mu, int N_nodes , gsl_rng* randgsl, int verbose, int self_opt, int max_reps);
W_GRAPH* fixedk_geometric_directed_graph(double**x, double mu,  int N_nodes , gsl_rng* randgsl, int verbose, int self_opt, int max_reps);
W_GRAPH* fixedk_negbinom_undirected_graph(double*x,  double mu,  int N_nodes ,int layers, gsl_rng* randgsl, int verbose, int self_opt, int max_reps);
W_GRAPH* fixedk_negbinom_directed_graph(double**x, double mu, int N_nodes , int layers, gsl_rng* randgsl, int verbose, int self_opt, int max_reps);
W_GRAPH* fixedk_binom_directed_graph(double**x, double mu,  int N_nodes , int layers, gsl_rng* randgsl, int verbose, int self_opt, int max_reps);
W_GRAPH* fixedk_binom_undirected_graph(double*x,  double mu, int N_nodes , int layers, gsl_rng* randgsl, int verbose, int self_opt, int max_reps);
W_GRAPH* fixedk_bernouilli_directed_graph(double**x, int N_nodes , gsl_rng* randgsl, int verbose, int self_opt);
W_GRAPH* fixedk_bernouilli_undirected_graph(double*x, int N_nodes , gsl_rng* randgsl, int verbose, int self_opt);

/*********** Linear&Binary constraints **************/
W_GRAPH* fixedEs_poisson_undirected_graph(double*x,  double lam, int N_nodes , gsl_rng* randgsl, int verbose, int self_opt, int max_reps);
W_GRAPH* fixedEs_poisson_directed_graph(double**x,  double lam, int N_nodes , gsl_rng* randgsl, int verbose, int self_opt, int max_reps);

W_GRAPH* fixedks_poisson_undirected_graph(double**x, int N_nodes , gsl_rng* randgsl, int verbose, int self_opt, int max_reps);
W_GRAPH* fixedks_poisson_directed_graph(double**x, int N_nodes , gsl_rng* randgsl, int verbose, int self_opt, int max_reps);
W_GRAPH* fixedks_binomial_undirected_graph(double**x, int N_nodes , int layers, gsl_rng* randgsl, int verbose, int self_opt, int max_reps);
W_GRAPH* fixedks_binomial_directed_graph(double**x, int N_nodes , int layers, gsl_rng* randgsl, int verbose, int self_opt, int max_reps);
/*********Distros of p's*******/
double * prob_mult_s_undir(double* x, int N_nodes, int self_opt);
double * prob_mult_s_dir(double** x, int N_nodes, int self_opt);
/********* Aux funcs *******/
int compute_T(int N_nodes, double av_k, int* x);

#endif
