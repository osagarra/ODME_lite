#ifndef ULA_STAT_FUNCS_H
#define ULA_STAT_FUNCS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_randist.h>
#include "ula_gen_funcs.h"

/****** Basci stats ********/
double mean_vec_int(int * vec, int len);
double mean_vec_double(double * vec, int len);
double var_vec_int(int * vec, int len);
double var_vec_double(double * vec, int len);
double** y_of_x(double* vectx, double* vecty, double* xrange, int len, int len_ranges);
/****** Histograms and Cumulatives ********/
int * histogram_int(int* vect, int minx, int maxx, int len );
double * normalize_hist_int(int* hist, int len);
/****** GSL ********/
/* 1 d*/
gsl_histogram * histogram_double_log(double* vect, double minx, double maxx, double expo, int len);
gsl_histogram * histogram_double(double* vect, double minx, double maxx, int nbins, int len);
void normalize_gsl_hist(gsl_histogram * hist);
void normalize_gsl_log_hist(gsl_histogram * hist);
/* 2 d*/
gsl_histogram2d * histogram_2d_double_log(double* vectx, double* vecty, double minx, double maxx, double miny, double maxy, double expx, double expy, int lenx, int leny);
gsl_histogram2d * histogram_2d_double_linlog(double* vectx, double* vecty, double minx, double maxx, double miny, double maxy, int xbins, double expy, int lenx, int leny);
gsl_histogram2d * histogram_2d_double(double* vectx, double* vecty, double minx, double maxx, double miny, double maxy, int xbins, int ybins, int lenx, int leny);
void histogram_2d_mean(gsl_histogram2d * hist, int axis, double * mean_v, double * std_v, double* xrange1, double* xrange2, int to_int);
/****** Bins ********/
double * log_bins_double(double minx, double maxx, double expo, int * numbins);
/************ Accumulators *************/
/************ 1d *************/
gsl_histogram * set_acc_int(int minv, int maxv);
gsl_histogram * set_acc_double(double minv, double maxv, int bins);
void update_int_acc(int * vect, int len, gsl_histogram * acc_k, gsl_histogram * acc_k2, int norm);
void update_double_acc(double * k, int len, gsl_histogram * acc_k, gsl_histogram * acc_k2, int norm);

void print_acc(char *input_name, gsl_histogram * acc1, gsl_histogram * acc2);
/************ 2d *************/
gsl_histogram2d * set_acc2d_linlog(double minx, double maxx, double miny, double maxy, int binsx, double exp_y);
gsl_histogram2d * set_acc2d_log(double minx, double maxx, double miny, double maxy, double exp_x, double exp_y);
void acc2d_compute_std(gsl_histogram2d ** accs,int len, double norm);
void acc2d_normalize(gsl_histogram2d ** accs,int len, double norm);
void acc2d_allocate_all(gsl_histogram2d ** accs, double ** d, int ** s, int N_nodes, double exp_s, double exp_t, int binsx);
void acc2d_free_all(gsl_histogram2d ** accs, int len);

void print_acc2d(char *input_name, gsl_histogram2d * acc1, gsl_histogram2d * acc2);
/****** Accumulator of accumulators ********/
void acc_compute_std(gsl_histogram ** accs,int len, double norm);
void acc_normalize(gsl_histogram ** accs,int len, double norm);
void acc_allocate_all(gsl_histogram ** accs, double ** d, int ** s, int N_nodes);
void acc_free_all(gsl_histogram ** accs, int len);

/****** Print Hists ********/
void print_hist_double(char *input_name, int len, double *av_hist, double *av_hist2);
void print_hist_int(char *input_name, int len, int *av_hist);
void print_hist2d_mean(char *input_name, double * h_mean, double * h_std, double * xrange, int len);
/************ CHECK FUNCS *************/
void test_hist_double(unsigned int samples);
void test_hist_int(unsigned int samples,double prob, int trials);
void test_hist_log(unsigned int samples, double expo);

#endif
