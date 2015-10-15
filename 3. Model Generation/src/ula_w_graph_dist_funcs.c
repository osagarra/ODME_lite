/************************************************************
 *
 *                    W_Graph_dist Library
 *
 *		Functions useful when dealing with spatial directed weighted networks
 *		which are spatial or have an associated cost matrix.
 *
 *************************************************************/



#include "ula_w_graph_dist_funcs.h"

/****************************************************************************
 * Net functions *
 ****************************************************************************/

double** w_graph_dist_compute_k_analitic_from_xygamma_directed(double** x2, int N_nodes, int self_opt,double gamma, double** dist){
    double** k=cast_mat_double(2,N_nodes);
	double d;
    int i,j;
    for(i=0;i<N_nodes;i++)
    {
        k[0][i]=(double)N_nodes;
        k[1][i]=(double)N_nodes;
        if(self_opt<=0)
        {
            k[0][i]-=1;
            k[1][i]-=1;
        }
        for(j=0;j<N_nodes;j++)
        {
            if(j!=i)
            {
				d=dist[maxeq_int(i,j)][mineq_int(i,j)];
                k[0][i]-=exp(-(double)x2[0][i]*((double)x2[1][j])*exp(-gamma*d));
                k[1][i]-=exp(-(double)x2[1][i]*((double)x2[0][j])*exp(-gamma*d));
            }else{
				if(self_opt>0)
				{
		    		k[0][i]-=exp(-(double)x2[0][i]*((double)x2[1][j]));
		    		k[1][i]-=exp(-(double)x2[1][i]*((double)x2[0][j]));
				}
	    	}
        }
    }
    return k;
}

double * w_graph_dist_compute_k_analitic_from_xgamma_undirected(double* x, int N_nodes, int self_opt, double gamma, double** dist){
    double* k=cast_vec_double(N_nodes);
	double d;
    int i,j;
    for(i=0;i<N_nodes;i++)
    {
        k[i]=(double)N_nodes;
        if(self_opt<=0) k[i]-=1;
        for(j=0;j<N_nodes;j++)
        {
            if(j!=i)
            {
				d=dist[maxeq_int(i,j)][mineq_int(i,j)];
                k[i]-=exp(-(double)x[i]*((double)x[j]*exp(-gamma*d)));
            }else{
				if(self_opt>0)
				{
		    	k[i]-=exp(-(double)x[i]*((double)x[i]));
				}
	    	}
        }
    }
    return k;
}



/****************************************************************************
 * Distance handling * // 
 ****************************************************************************/


/****************************************************************************
 * Net functions (from adj matrix) * // --> Edge averages <w_ij> (from integer weights)


 ****************************************************************************/
double *w_graph_dist_compute_sij_edges(W_GRAPH* WG, double **d, int N_nodes, gsl_rng * randgsl){
	int E1,i,j,aux,dest;
	E1 = w_graph_total_edges(WG, N_nodes);
	double* d_list=cast_vec_double(E1);
	int** s=w_graph_compute_s(WG,N_nodes);
	aux=0;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<WG->node[i].kout;j++)
        {
			dest = WG->node[i].out[j];
            d_list[aux]=(double)w_graph_dist_compute_sij(s,d,i,dest,N_nodes,randgsl,1); // s_ij origin
            aux++;
            if(WG->node[i].out[j]==i) // count twice
            {
                d_list[aux]=0;
                aux++;
            }
        }
    }
	assert(aux==E1);
	return d_list;
}

double *w_graph_dist_compute_s_out_edges(W_GRAPH* WG, int N_nodes){
	int E1,i,j,aux,dest;
	E1 = w_graph_total_edges(WG, N_nodes);
	double* d_list=cast_vec_double(E1);
	aux=0;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<WG->node[i].kout;j++)
        {
			dest = WG->node[i].out[j];
			d_list[aux]=(double)WG->node[dest].sout;
            aux++;
            if(WG->node[i].out[j]==i) // count twice
            {
				d_list[aux]=(double)WG->node[i].sout;
                aux++;
            }
        }
    }
	assert(aux==E1);
	return d_list;
}

double *w_graph_dist_compute_s_in_edges(W_GRAPH* WG, int N_nodes){
	int E1,i,j,aux,dest;
	E1 = w_graph_total_edges(WG, N_nodes);
	double* d_list=cast_vec_double(E1);
	aux=0;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<WG->node[i].kout;j++)
        {
			dest = WG->node[i].out[j];
			d_list[aux]=(double)WG->node[dest].sin;
            aux++;
            if(WG->node[i].out[j]==i) // count twice
            {
				d_list[aux]=(double)WG->node[i].sin;
                aux++;
            }
        }
    }
	assert(aux==E1);
	return d_list;
}


/****************************************************************************
 * Net functions (from adj matrix) * // --> Distance edge-trips
 ****************************************************************************/
/*********************************************/
int w_graph_dist_compute_sij(int  **s, double ** dist, int origin, int dest, int N_nodes, gsl_rng * randgsl, int perturb){
	// perturb si perturb >0
	int i;
	int sij=0;
	double dd=dist[maxeq_int(origin,dest)][mineq_int(origin,dest)];
	double distt;
	for(i=0;i<N_nodes;i++)
	{
		distt=dist[maxeq_int(i,origin)][mineq_int(i,origin)];
		if(perturb>0) distt+=(gsl_rng_uniform(randgsl)-0.5);
		if((i!=origin)&&(i!=dest)&&(distt<dd)) // if distance minor or equal,
		{
			/*if(fabs(distt-dd)<eps)//if equal undraw
			{
				r=gsl_rng_uniform(randgsl);
				if(r>0.5)
				{
					sij+=s[i][1];
				}
			}else{
				sij+=s[i][1];
			}*/
			sij+=s[1][i];
		}
	}
	return sij;
}


/********************** Edges and trips distances ****************************************/


double * w_graph_dist_compute_d_edges(W_GRAPH* WG, double ** d, int N_nodes, int* num_edges){
	int i,j,aux,dest;
    (*num_edges)=w_graph_total_edges(WG, N_nodes);
	double * d_edges=cast_vec_double(*num_edges);
	aux=0;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<WG->node[i].kout;j++)
        {
			dest = WG->node[i].out[j];
			d_edges[aux]=d[maxeq_int(i,dest)][mineq_int(i,dest)];;
            aux++;
            if(WG->node[i].out[j]==i) // count twice
            {
                d_edges[aux]=0;
                aux++;
            }
        }
    }
    assert(aux==*num_edges);
	//d_edges=realloc(d_edges,aux*sizeof(double)); // realloc memory
	//(*num_edges)=aux;
	return d_edges;
}
double * w_graph_dist_compute_d_trips(W_GRAPH* WG, double ** d, int N_nodes, int* num_edges){
	int i,j,k,aux,dest;
	double dd;
    (*num_edges)=w_graph_total_weight(WG, N_nodes);
	double * d_edges=cast_vec_double(*num_edges);
	aux=0;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<WG->node[i].kout;j++)
        {
			dest = WG->node[i].out[j];
			dd = d[maxeq_int(i,dest)][mineq_int(i,dest)];
            for(k=0;k<WG->node[i].w_out[j];k++)
            {
				d_edges[aux]=dd;
            	aux++;
            }
        }
    }
    assert(aux==*num_edges);
	//d_edges=realloc(d_edges,aux*sizeof(double)); // realloc memory
	//(*num_edges)=aux;
	return d_edges;
}

/********************** p_ij ****************************************/
gsl_histogram * w_graph_dist_compute_pij(W_GRAPH* WG, double ** d, int bins, int N_nodes, double max_d_edge, int opt_self, int opt_dir){
	int E,D;
    double d_min;
	double * d_edges=w_graph_dist_compute_d_edges(WG, d, N_nodes, &E);
	double * dists= flatten_matrix_triangular_double(d, N_nodes,&D);
	//double max_d_edge= ceil(max_value_double(dists, D));
	//double max_d= max_value_double(dists, D);
    if(opt_self>0)
    {
        d_min = 1;
    }else{
        d_min = 0;
    }
	gsl_histogram * d_edge_hist=histogram_double(d_edges, d_min, max_d_edge , bins, E);
	gsl_histogram * d_hist=histogram_double(dists, d_min, max_d_edge , bins, D);
	gsl_histogram_div(d_edge_hist, d_hist);
    if(opt_dir>0)
    {
        gsl_histogram_scale(d_edge_hist,0.5); // multiply by 0.5
    }
	return d_edge_hist;
}


/*************************************************************************/
/********************** Net stats ****************************************/
/*************************************************************************/
void w_graph_dist_all_stats(W_GRAPH* WG, int N_nodes, int run, double bin_exp, int opt_dir, int self_opt, double** dist,gsl_rng * randgsl, double dmax, int verbose){
    if(verbose>0) printf("... Computing graph edge stats...\n");
    //w(sin sout), w(kin,kout), w
    int **k=w_graph_compute_k(WG, N_nodes);
    int **s=w_graph_compute_s(WG, N_nodes);
    int E,E2,E3,q,L;
    int *w=w_graph_compute_w(WG, N_nodes, &E, -1);
    int *w_zeros=w_graph_compute_w(WG, N_nodes, &L, 1);
    int *p_zeros=w_graph_compute_p(WG, N_nodes, &L);

    double *wss=w_graph_compute_wp_ss(WG, N_nodes, 1);
    double *wkk=w_graph_compute_wp_ss(WG, N_nodes, -1);
    double *wss_zeros=w_graph_compute_w_ss(WG, N_nodes, 1);

	
	double* d_edges= w_graph_dist_compute_d_edges(WG, dist, N_nodes, &E2);
	double* d_trips= w_graph_dist_compute_d_trips(WG, dist, N_nodes, &E3);
	double* s_in_d = w_graph_dist_compute_s_in_edges(WG, N_nodes);
	double* s_out_d = w_graph_dist_compute_s_in_edges(WG, N_nodes);
	//printf("d_sij\n"); fflush(stdout); --> extremely slow!
	//double* s_ij_d = w_graph_dist_compute_sij_edges(node, dist, N_nodes, randgsl);  
    char cadena[100];

	//printf("diff: %d\n", (int)(E-E2));fflush(stdout);
	//assert(E==E2); // security
	//int TT = sum_vec_int(s[0],N_nodes);
    //assert(TT == E3); // security debugging
    gsl_histogram* h1;
    
    double* sout;
    double* xranges;
    int xbins;
    double** yy;
    
    /// w histogram /////
    sout=vec_int_to_double(w,E);
    q=max_value_int(w,E);
    if(verbose>0) printf("\t Max weight: %d | av existing weight: %.3f+-%.3f \n",q,mean_vec_int(w, E),sqrt(var_vec_int(w, E)));
    h1=histogram_double(sout,0,q,q,E);
    free(sout);
    //sprintf(cadena,"run_%dN%d_w.hist",run,N_nodes);
	if(opt_dir>0)
	{
	    sprintf(cadena,"N%d_w.hist",N_nodes);		
	}else{
	    sprintf(cadena,"N%d_undir_w.hist",N_nodes);		
	}
    print_acc(cadena, h1, h1);
    gsl_histogram_free(h1);

    /// edge distance histogram /////
    //double q2=max_value_double(d_edges,E);
    //h1=histogram_double(d_edges,0,q2,200,E2);
	h1=histogram_double(d_edges,0,dmax,100,E2);
	//normalize_gsl_hist(h1);	// we don't normalize it so we compare trips to raw edges
    //sprintf(cadena,"run_%dN%d_w.hist",run,N_nodes);
	if(opt_dir>0)
	{
	    sprintf(cadena,"N%d_d_edges.hist",N_nodes);
	}else{
	    sprintf(cadena,"N%d_undir_d_edges.hist",N_nodes);
	}
    print_acc(cadena, h1, h1);
    gsl_histogram_free(h1);

    /// trip distance histogram /////
    //h1=histogram_double(d_trips,0,q2,200,E3);
    h1=histogram_double(d_trips,0,dmax,100,E3);
	//normalize_gsl_hist(h1);	// we don't normalize it so we compare trips to raw edges
    //sprintf(cadena,"run_%dN%d_w.hist",run,N_nodes);
	if(opt_dir>0)
	{
	    sprintf(cadena,"N%d_d_trips.hist",N_nodes);
	}else{
	    sprintf(cadena,"N%d_undir_d_trips.hist",N_nodes);
	}
    print_acc(cadena, h1, h1);
    gsl_histogram_free(h1);


    /// connection probability as func of distance /////
	h1 = w_graph_dist_compute_pij(WG, dist, 100, N_nodes, dmax, self_opt, opt_dir);
	if(opt_dir>0)
	{
		sprintf(cadena,"N%d_pij.hist",N_nodes);
	}else{
		sprintf(cadena,"N%d_undir_pij.hist",N_nodes);
	}
	print_acc(cadena, h1, h1);
	gsl_histogram_free(h1);


    /// exsiting weight as func of ss /////
    sout=vec_int_to_double(w,E);
    xranges=log_bins_double(0, max_value_double(wss,E) , bin_exp, &xbins);
    yy=y_of_x(wss, sout, xranges,  E,  xbins);
	if(opt_dir>0)
	{
	    sprintf(cadena,"N%d_wp_s_oi.hist",N_nodes);
	}else{
	    sprintf(cadena,"N%d_undir_w_s_oi.hist",N_nodes);
	}
    print_hist2d_mean(cadena, yy[1], yy[2], yy[0], xbins-1);
    free(sout);
    free(wss);
    free(xranges);
    free_mat_double(yy,4);

    /// exsiting weight as func of kk /////
    sout=vec_int_to_double(w,E);
    xranges=log_bins_double(0, max_value_double(wkk,E) , bin_exp, &xbins);
    yy=y_of_x(wkk, sout, xranges,  E,  xbins);
	if(opt_dir>0)
	{
	    sprintf(cadena,"N%d_wp_k_oi.hist",N_nodes);
	}else{
	    sprintf(cadena,"N%d_undir_w_k_oi.hist",N_nodes);
	}
    print_hist2d_mean(cadena, yy[1], yy[2], yy[0], xbins-1);
    free(sout);
    free(wkk);
    free(xranges);
    free_mat_double(yy,4);

    /// average weight as func of ss /////
    sout=vec_int_to_double(w_zeros,L);
    xranges=log_bins_double(0, max_value_double(wss_zeros,L) , bin_exp, &xbins);
    yy=y_of_x(wss_zeros, sout, xranges,  L,  xbins);
	if(opt_dir>0)
	{
	    sprintf(cadena,"N%d_w_s_oi.hist",N_nodes);
	}else{
	    sprintf(cadena,"N%d_undir_w_s_oi.hist",N_nodes);
	}
    print_hist2d_mean(cadena, yy[1], yy[2], yy[0], xbins);
    free(sout);
    free(xranges);
    free_mat_double(yy,4);

    /// conn prob as func of ss /////
    sout=vec_int_to_double(p_zeros,L);
    xranges=log_bins_double(0, max_value_double(wss_zeros,L) , bin_exp, &xbins);
    yy=y_of_x(wss_zeros, sout, xranges,  L,  xbins);
	if(opt_dir>0)
	{
	    sprintf(cadena,"N%d_p_s_oi.hist",N_nodes);
	}else{
	    sprintf(cadena,"N%d_undir_p_s_oi.hist",N_nodes);
	}
    print_hist2d_mean(cadena, yy[1], yy[2], yy[0], xbins);
    free(sout);
    free(wss_zeros);
    free(xranges);
    free_mat_double(yy,4);  



    /// exsiting weight as func of s_in /////
    sout=vec_int_to_double(w,E);
    xranges=log_bins_double(0, max_value_double(s_in_d,E) , bin_exp, &xbins);
    yy=y_of_x(s_in_d, sout, xranges,  E,  xbins);
	if(opt_dir>0)
	{
	    sprintf(cadena,"N%d_w_s_i.hist",N_nodes);
	}else{
	    sprintf(cadena,"N%d_undir_w_s_i.hist",N_nodes);
	}
    print_hist2d_mean(cadena, yy[1], yy[2], yy[0], xbins-1);
    free(sout);
    free(s_in_d);
    free(xranges);
    free_mat_double(yy,4);

    /// exsiting weight as func of s_out /////
    sout=vec_int_to_double(w,E);
    xranges=log_bins_double(0, max_value_double(s_out_d,E) , bin_exp, &xbins);
    yy=y_of_x(s_out_d, sout, xranges,  E,  xbins);
	if(opt_dir>0)
	{
	    sprintf(cadena,"N%d_w_s_o.hist",N_nodes);
	}else{
	    sprintf(cadena,"N%d_undir_w_s_o.hist",N_nodes);
	}
    print_hist2d_mean(cadena, yy[1], yy[2], yy[0], xbins-1);
    free(sout);
    free(s_out_d);
    free(xranges);
    free_mat_double(yy,4);

	/*
    sout=vec_int_to_double(w,E);
    xranges=log_bins_double(0, max_value_double(s_ij_d,E) , 1.05, &xbins);
    yy=y_of_x(s_ij_d, sout, xranges,  E,  xbins);
    sprintf(cadena,"N%d_w_s_ij.hist",N_nodes);
    print_hist2d_mean(cadena, yy[1], yy[2], yy[0], xbins-1);
    free(sout);
    free(s_ij_d);
    free(xranges);
    free_mat_double(yy,4);
	*/

    /// exsiting weight as func of d /////
    sout=vec_int_to_double(w,E);
    xranges=log_bins_double(0, dmax , 1.05, &xbins);
    yy=y_of_x(d_edges, sout, xranges,  E,  xbins);
	if(opt_dir>0)
	{
	    sprintf(cadena,"N%d_w_dij.hist",N_nodes);
	}else{
	    sprintf(cadena,"N%d_undir_w_dij.hist",N_nodes);
	}
    print_hist2d_mean(cadena, yy[1], yy[2], yy[0], xbins-1);
    free(sout);
    //free(d_edges);
    free(xranges);
    free_mat_double(yy,4);

//// free all
    free_mat_int(s,2);
    free_mat_int(k,2);
    free(w);
    free(w_zeros);
    free(p_zeros);
	free(d_trips);
	free(d_edges);
    return;
}

/****************************************************************************
 * Ensembles averaging stats *
 ****************************************************************************/
gsl_histogram ** w_graph_dist_all_stats_ensemble_allocate(int dir, double d_max, int w_max){
	int len_acc=6; // P(w),p(delta_r),p(delta_r)_Edges,
	int bins;
	if(w_max>1e5)
	{
		bins = 100000;
	}else{
		bins=w_max;
	}
	//P(w)
	gsl_histogram ** acc = (gsl_histogram**)malloc(sizeof(gsl_histogram*)*len_acc);
	acc[0] = set_acc_double(0,w_max,bins);
	acc[1] = set_acc_double(0,w_max,bins);
	bins = (int)(d_max/500.); //minimal grid
	if(bins<5 || bins>200)
	{
		bins=50; // final tunning
	}
	//P(delta_r)
	acc[2] = set_acc_double(0,d_max,bins);
	acc[3] = set_acc_double(0,d_max,bins);
	//P(delta_r_edges)
	acc[4] = set_acc_double(0,d_max,bins);
	acc[5] = set_acc_double(0,d_max,bins);
		
	return acc;
}

void w_graph_dist_all_stats_ensemble_update(gsl_histogram** acc, W_GRAPH* WG, int N_nodes, int dir,double** dist){
    int E,E2,E1;
	int *w=w_graph_compute_w(WG, N_nodes, &E, -1);
	//int T = sum_vec_int(w,E);
	double* d_edges= w_graph_dist_compute_d_edges(WG, dist, N_nodes, &E1);
	assert(E1==E);
	double* d_trips= w_graph_dist_compute_d_trips(WG, dist, N_nodes, &E2);
	//printf("Lengths: E:%d T:%d E2:%d wei:%d\n",E,T,E1,E2);fflush(stdout);
	//assert(E2 == T);
	update_int_acc(w, E, acc[0], acc[1], 1);
	update_double_acc(d_trips, E2, acc[2], acc[3], 1);
	update_double_acc(d_edges, E1, acc[4], acc[5], 1);
	free(w);
	free(d_edges);
	free(d_trips);
	return;
}


void w_graph_dist_all_stats_ensemble_print(gsl_histogram** acc, int len, int reps, int N_nodes, int opt_dir){
	char	cadena[100];
	//Normalize and std
	acc_normalize(acc,len, 1./(double)reps);
	acc_compute_std(acc,len, 1.);
	// Print all
	/*
	if(dir==1)
	{
		sprintf(cadena,"N%dexpo%.2f_ens_r%d_sIN.hist",N_nodes,expo-1,reps);
		print_acc(cadena, acc[0], acc[1]);
		sprintf(cadena,"N%dexpo%.2f_ens_r%d_sOUT.hist",N_nodes,expo-1,reps);
		print_acc(cadena, acc[2], acc[3]);
		sprintf(cadena,"N%dexpo%.2f_ens_r%d_kIN.hist",N_nodes,expo-1,reps);
		print_acc(cadena, acc[4], acc[5]);
		sprintf(cadena,"N%dexpo%.2f_ens_r%d_kOUT.hist",N_nodes,expo-1,reps);
		print_acc(cadena, acc[6], acc[7]);
		sprintf(cadena,"N%dexpo%.2f_ens_r%d_w.hist",N_nodes,expo-1,reps);
		print_acc(cadena, acc[8], acc[9]);
	}else{
		sprintf(cadena,"N%dexpo%.2f_ens_r%d_sOUT.hist",N_nodes,expo-1,reps);
		print_acc(cadena, acc[0], acc[1]);
		sprintf(cadena,"N%dexpo%.2f_ens_r%d_kOUT.hist",N_nodes,expo-1,reps);
		print_acc(cadena, acc[2], acc[3]);
		sprintf(cadena,"N%dexpo%.2f_ens_r%d_w.hist",N_nodes,expo-1,reps);
		print_acc(cadena, acc[4], acc[5]);
	}
	*/
	if(opt_dir>0)
	{
		sprintf(cadena,"N%d_ens_r%d_w.hist",N_nodes,reps);		
	}else{
		sprintf(cadena,"N%d_undir_ens_r%d_w.hist",N_nodes,reps);		
	}

	print_acc(cadena, acc[0], acc[1]);
	if(opt_dir>0)
	{
		sprintf(cadena,"N%d_ens_r%d_d_edges.hist",N_nodes,reps);
	}else{
		sprintf(cadena,"N%d_undir_ens_r%d_d_edges.hist",N_nodes,reps);
	}
	print_acc(cadena, acc[4], acc[5]);
	if(opt_dir>0)
	{
		sprintf(cadena,"N%d_ens_r%d_d_trips.hist",N_nodes,reps);
	}else{
		sprintf(cadena,"N%d_undir_ens_r%d_d_trips.hist",N_nodes,reps);
	}
	print_acc(cadena, acc[2], acc[3]);

	//Free all
	acc_free_all(acc, len);
	return;
}

/**************************/
/******* Net entropy ******/
/**************************/

double w_graph_entropy_poisson_dist(W_GRAPH* WG, double** x,int N_nodes, double** dist, double gamma, int opt_self, int opt_dir){
	double S,p,mu;
	int i,j,t;
	int out_d,in_d;
    double *y;
    double *yy;
	S = 0;
    if(opt_dir>0)
    {
        y = x[0];
        yy = x[1];
    }else{
        y = x[0];
        yy = x[0];
    }
	for(i=0;i<N_nodes;i++) // all edges
	{
		for(j=0;j<N_nodes;j++)
		{
            if((opt_self>0)||(i!=j))
            {
                t=find_value_int(WG->node[i].out, j, WG->node[i].kout); // check if link exists
                out_d=maxeq_int(i,j);
				in_d=mineq_int(i,j);
                if(t>=0)
                {
					t = WG->node[i].w_out[t];
				}else{
					t = 0;
				}
                mu = y[i]*yy[j]*exp(-gamma*dist[out_d][in_d]);
                if(mu>0)
				{
					p = gsl_ran_poisson_pdf (t, mu);
					//printf("Mu:%f p:%f t:%d",mu,p,t);fflush(stdout);
				}else{
					p = 0;
				}
				S+= p*log(p);
			}
        }
	}
	if(opt_dir<=0) S = S/2.;
	return -S;    
}
