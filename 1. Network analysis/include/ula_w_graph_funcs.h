#ifndef ULA_W_GRAPH_FUNCS_H
#define ULA_W_GRAPH_FUNCS_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include "ula_gen_funcs.h"
#include "ula_stat_funcs.h"
#include "ula_graph_structs.h"


/******* Read funcs ********/
W_GRAPH* w_graph_read_edge_list(char *input_name, int num_nodes, int opt_dir, int header, int verbose);
/******* Allocation ********/
W_GRAPH* w_graph_alloc(int N_nodes);
void w_graph_free_nodes(W_GRAPH* WG, int N_nodes);

/******* Add edge ********/
void w_graph_add_multi_link(W_GRAPH * WG, int N_nodes, int origin, int dest, int weight);
void w_graph_add_multi_link_undirected(W_GRAPH * WG, int N_nodes, int origin, int dest, int weight);

/******* Printing ********/
void w_graph_print_adj_list(W_GRAPH* WG, int N_nodes, char* output);

/******* node features ********/
int ** w_graph_compute_s(W_GRAPH* WG, int N_nodes);
int ** w_graph_compute_k(W_GRAPH* WG, int N_nodes);
int w_graph_nonempty_nodes(W_GRAPH* WG);
double ** w_graph_compute_k_analitic(W_GRAPH* WG, int N_nodes, int self_opt);
double w_graph_compute_E_binary_undirected(double* x, int N_nodes, int self_opt);
double * w_graph_compute_k_binary_undirected(double* x, int N_nodes, int self_opt);
double w_graph_compute_E_binary_directed(double** x, int N_nodes, int self_opt);
double ** w_graph_compute_k_binary_directed(double** x, int N_nodes, int self_opt);
double ** w_graph_compute_k_analitic_from_s_directed(int** s, int N_nodes, int self_opt);
double ** w_graph_compute_k_analitic_from_s_undirected(int* s, int N_nodes, int self_opt);
double w_graph_compute_T_analytic_from_xy(double** x2, int N_nodes, int layers, double gamma, int opt_dir, int opt_self, int opt_indist, int opt_agg, int fixk);
/******* Snn & knn's ********/
double ** w_graph_compute_s_nn(W_GRAPH* WG, int N_nodes, int weight, int opt_dir);
/******* Y2 ********/
double ** w_graph_compute_Y2(W_GRAPH * WG, int N_nodes, int opt_dir);
//double ** w_graph_compute_xy(W_GRAPH * WG, int N_nodes);

/******* Clsutering ********/
double ** w_graph_compute_clust(W_GRAPH * WG, int N_nodes);
/******* w's ********/
int * w_graph_compute_w(W_GRAPH* WG, int N_nodes, int* aux, int zeros);
int * w_graph_compute_p(W_GRAPH* WG, int N_nodes, int* aux);
double** w_graph_compute_p_w_analitic_from_s_undirected(int maxt, double binn, int* s, int N_nodes, int self_opt, int* len);
double** w_graph_compute_p_w_analitic_from_s_directed(int maxt, double binn, int** s, int N_nodes, int self_opt, int* len);
double * w_graph_compute_wp_ss(W_GRAPH* WG, int N_nodes, int weight);
double * w_graph_compute_w_ss(W_GRAPH* WG, int N_nodes, int weight);
double * w_graph_compute_w_s(W_GRAPH* WG, int N_nodes, int weight, int direction);
double * w_graph_compute_wp_s(W_GRAPH* WG, int N_nodes, int weight, int direction);
/******* Distances *****/
//double ** w_graph_compute_dists(W_GRAPH* WG, int N_nodes);
/*******  Entropies *****/
double w_graph_ML_ensemble_entropy_multinomial(W_GRAPH* WG, int N_nodes, int opt_dir);
double w_graph_surprise_poisson_pij(W_GRAPH* WG,int N_nodes, double** pij, int opt_self, int opt_dir);
double w_graph_surprise_poisson(W_GRAPH* WG, double** x,int N_nodes, int opt_self, int opt_dir);
double w_graph_surprise_geometric(W_GRAPH* WG, double**x, int N_nodes, int opt_self, int opt_dir);
double w_graph_surprise_binomial(W_GRAPH* WG, double**x, int N_nodes, int layers, int opt_self, int opt_dir);
double w_graph_surprise_negbinomial(W_GRAPH* WG, double**x, int N_nodes, int layers, int opt_self, int opt_dir);
double w_graph_surprise_bernouilli(W_GRAPH* WG, double** x, int N_nodes, int opt_self, int opt_dir);

double w_graph_surprise_ZIP(W_GRAPH* WG, double** x, int N_nodes, double gamma, int opt_self, int opt_dir);
double w_graph_surprise_ZIP2(W_GRAPH* WG, double** x, int N_nodes, int opt_self, int opt_dir);
double w_graph_surprise_ZIB2(W_GRAPH* WG, double** x, int N_nodes, int layers, int opt_self, int opt_dir);

void w_graph_print_entropy(double* seq,  int len,char* output);

/*******  Loglikelyhoods *****/
double w_graph_loglikelyhood_poisson(W_GRAPH* WG,int N_nodes,double** wij, int opt_self, int opt_dir);
double w_graph_loglikelyhood_poisson_xy(W_GRAPH* WG, double** x,int N_nodes, int opt_self, int opt_dir);
double w_graph_loglikelyhood_geometric_xy(W_GRAPH* WG, double**x, int N_nodes, int opt_self, int opt_dir);
double w_graph_loglikelyhood_binomial_xy(W_GRAPH* WG, double**x, int N_nodes, int layers, int opt_self, int opt_dir);
double w_graph_loglikelyhood_negbinomial_xy(W_GRAPH* WG, double**x, int N_nodes, int layers, int opt_self, int opt_dir);
double w_graph_loglikelyhood_bernouilli_xy(W_GRAPH* WG, double** x, int N_nodes, int opt_self, int opt_dir);
double w_graph_loglikelyhood_ZIP_xy(W_GRAPH* WG, double** x, int N_nodes, double gamma, int opt_self, int opt_dir);
double w_graph_loglikelyhood_ZIP2_xy(W_GRAPH* WG, double** x, int N_nodes, int opt_self, int opt_dir);
double w_graph_loglikelyhood_ZIB2_xy(W_GRAPH* WG, double** x, int N_nodes, int layers, int opt_self, int opt_dir);

/******* All_stats *****/
void w_graph_node_stats_list(W_GRAPH* WG, int N_nodes, int run, int opt_dir, int opt_clust, int self_opt, int verbose);
int w_graph_total_weight( W_GRAPH* WG, int N_nodes);
int w_graph_total_edges( W_GRAPH* WG, int N_nodes);
int w_graph_total_edgepairs( W_GRAPH* WG, int N_nodes);
void w_graph_all_stats(W_GRAPH* WG, int N_nodes, int run, double bin_exp,int opt_dir, int self_opt, int w_anal, int verbose);

/******* Ensemble stats *****/
void w_graph_node_stats_ensemble(W_GRAPH* WG, int N_nodes, double** container, double** container2, int** WG_nonzero, double* T_container, int opt_dir, int opt_clust );
void w_graph_node_stats_ensemble_print(int reps, int N_nodes, double* Tcont, double** cont, double ** cont2, int** WG_nonzero, double bin_exp, int len_acc, int opt_dir);

gsl_histogram ** w_graph_all_stats_ensemble_allocate(int dir, int s_min, int s_max, int k_min, int k_max, int w_max);
void w_graph_all_stats_ensemble_update(gsl_histogram** acc, W_GRAPH* WG, int N_nodes, int dir);
void w_graph_all_stats_ensemble_print(gsl_histogram** acc, int len, int reps, int N_nodes, int opt_dir);
/******* Other funcs *********/
double w_graph_compute_rho(double E_av, double T, int opt_indist);
/******* Transformations ********/
double** w_graph_to_adj_matrix(W_GRAPH* WG, int N_nodes);
/****************************************************************************
 * Graph algebra *
 ****************************************************************************/
void w_graph_sum_graphs(W_GRAPH* WG1, W_GRAPH* WG2, int N_nodes);
void w_graph_normalize(W_GRAPH* WG, int reps, int N_nodes);


/****************************************************************************
 * Graph compare *
 ****************************************************************************/
double w_graph_compute_sorensen(W_GRAPH* WG, W_GRAPH* WGoriginal, int N_nodes);
double w_graph_compute_sorensen_av(W_GRAPH* WGoriginal, double** pij, int N_nodes, double T);

/****************************************************************************
 * Graph Filtering *
 ****************************************************************************/
W_GRAPH* w_graph_filter_xij(W_GRAPH* WG, double* x, double* y, int N_nodes, double gamma, int mode, int M, int verbose);
double find_tmintmax_xy(int* tt,double xy, double gamma, int mode, int M);
#endif
