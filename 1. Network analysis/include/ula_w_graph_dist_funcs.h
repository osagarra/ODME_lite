#ifndef ULA_W_GRAPH_DIST_FUNCS_H
#define ULA_W_GRAPH_DIST_FUNCS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_fit.h>
#include "ula_gen_funcs.h"
#include "ula_stat_funcs.h"
#include "ula_graph_structs.h"
#include "ula_w_graph_funcs.h"


/******* Net funcs from W_GRAPH ******/
double** w_graph_dist_compute_k_analitic_from_xygamma_directed(double** x2, int N_nodes, int self_opt,double gamma, double** dist);
double * w_graph_dist_compute_k_analitic_from_xgamma_undirected(double* x, int N_nodes, int self_opt, double gamma, double** dist);
int w_graph_dist_compute_sij(int  **s, double ** dist, int origin, int dest, int N_nodes, gsl_rng * randgsl, int perturb);
double * w_graph_dist_compute_d_edges(W_GRAPH* WG, double ** d, int N_nodes, int* num_edges);
double * w_graph_dist_compute_d_trips(W_GRAPH* WG, double ** d, int N_nodes, int* num_edges);
gsl_histogram * w_graph_dist_compute_pij(W_GRAPH* WG, double ** d, int bins, int N_nodes, double max_d_edge);
/******* Net funcs from adj ******/
double *w_graph_dist_compute_sij_edges(W_GRAPH* WG, double **d, int N_nodes, gsl_rng * randgsl);
double *w_graph_dist_compute_s_out_edges(W_GRAPH* WG, int N_nodes);
double *w_graph_dist_compute_s_in_edges(W_GRAPH* WG, int N_nodes);
int w_graph_dist_compute_sij(int  **s, double ** dist, int origin, int dest, int N_nodes, gsl_rng * randgsl, int perturb);
/******* Stats ******/
void w_graph_dist_all_stats(W_GRAPH* WG, int N_nodes, int run, double bin_exp, double av_k, int opt_dir, int self_opt, double** dist, gsl_rng * randgsl, double dmax);
/******* Ensemble Stats ******/
gsl_histogram ** w_graph_dist_all_stats_ensemble_allocate(int dir, double d_max, int w_max);
void w_graph_dist_all_stats_ensemble_update(gsl_histogram** acc, W_GRAPH* WG, int N_nodes, int dir,double** dist);
void w_graph_dist_all_stats_ensemble_print(gsl_histogram** acc, int len, int reps, int N_nodes, double av_k, int opt_dir);

#endif
