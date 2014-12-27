#ifndef ULA_W_GRAPH_FUNCS_H
#define ULA_W_GRAPH_FUNCS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_fit.h>
#include "ula_gen_funcs.h"
#include "ula_stat_funcs.h"
#include "ula_graph_structs.h"


/******* Read funcs ********/
W_GRAPH* w_graph_dist_read_edge_list(char *input_name, int num_nodes, int opt_dir, int header);
/******* Allocation ********/
W_GRAPH* w_graph_alloc(int N_nodes);
void w_graph_free_nodes(W_GRAPH* WG, int N_nodes);

/******* Add edge ********/
void w_graph_add_multi_link(W_GRAPH * WG, int N_nodes, int origin, int dest, int weight);
void w_graph_add_multi_link_undirected(W_GRAPH * WG, int N_nodes, int origin, int dest, int weight);

/******* Printing ********/
void w_graph_print_adj_list(W_GRAPH* WG, int N_nodes, char* output);

/******* S & k's ********/
int ** w_graph_compute_s(W_GRAPH* WG, int N_nodes);
int ** w_graph_compute_k(W_GRAPH* WG, int N_nodes);
double ** w_graph_compute_k_analitic(W_GRAPH* WG, int N_nodes, int self_opt);
double ** w_graph_compute_k_analitic_from_s_directed(int** s, int N_nodes, int self_opt);
double ** w_graph_compute_k_analitic_from_s_undirected(int* s, int N_nodes, int self_opt);

/******* Snn & knn's ********/
double ** w_graph_compute_s_nn(W_GRAPH* WG, int N_nodes, int weight, int opt_dir);
/******* Y2 ********/
double ** w_graph_compute_Y2(W_GRAPH * WG, int N_nodes, int opt_dir);
//double ** w_graph_compute_xy(W_GRAPH * WG, int N_nodes);

/******* Clsutering ********/
double ** w_graph_compute_clust(W_GRAPH * WG, int N_nodes);
/******* w's ********/
int * w_graph_compute_w(W_GRAPH* WG, int N_nodes, int* aux, int zeros);
double** w_graph_compute_p_w_analitic_from_s_undirected(int maxt, double binn, int* s, int N_nodes, int self_opt, int* len);
double** w_graph_compute_p_w_analitic_from_s_directed(int maxt, double binn, int** s, int N_nodes, int self_opt, int* len);
double * w_graph_compute_w_ss(W_GRAPH* WG, int N_nodes, int weight);
/******* Distances *****/
//double ** w_graph_compute_dists(W_GRAPH* WG, int N_nodes);
/*******  Entropies *****/
double w_graph_entropy(W_GRAPH* WG, int N_nodes);
void w_graph_print_entropy(double* seq,  int len,char* output);
double w_graph_loglikelyhood(W_GRAPH* WG,int N_nodes, double** wij);
/******* All_stats *****/
void w_graph_node_stats_list(W_GRAPH* WG, int N_nodes, int run, double av_k, int opt_dir, int opt_clust, int self_opt);
int w_graph_total_weight( W_GRAPH* WG, int N_nodes);
int w_graph_total_edges( W_GRAPH* WG, int N_nodes);
void w_graph_all_stats(W_GRAPH* WG, int N_nodes, int run, double bin_exp,double av_k, int opt_dir, int self_opt, int w_anal);

/******* Ensemble stats *****/
void w_graph_node_stats_ensemble(W_GRAPH* WG, int N_nodes, double** container, double** container2, int** WG_nonzero, double* T_container, int opt_dir, int opt_clust );
void w_graph_node_stats_ensemble_print(int reps, int N_nodes, double* Tcont, double** cont, double ** cont2, int** WG_nonzero, double av_k, double bin_exp, int len_acc, int opt_dir);

gsl_histogram ** w_graph_all_stats_ensemble_allocate(int dir, int s_min, int s_max, int k_min, int k_max, int w_max);
void w_graph_all_stats_ensemble_update(gsl_histogram** acc, W_GRAPH* WG, int N_nodes, int dir);
void w_graph_all_stats_ensemble_print(gsl_histogram** acc, int len, int reps, int N_nodes, double av_k, int opt_dir);
/****************************************************************************
 * Graph algebra *
 ****************************************************************************/
void w_graph_sum_graphs(W_GRAPH* WG1, W_GRAPH* WG2, int N_nodes);
void w_graph_normalize(W_GRAPH* WG, int reps, int N_nodes);

#endif
