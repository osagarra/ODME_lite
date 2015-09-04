/************************************************************
 *
 *                    w_Graph Library
 *
 *		Functions useful when dealing with directed weighted networks in edge list format
 *
 *
 *************************************************************/



#include "ula_w_graph_funcs.h"

/****************************************************************************
 * I-O *
 ****************************************************************************/
W_GRAPH* w_graph_read_edge_list(char *input_name, int N_nodes, int opt_dir, int header){
	int i,j,tij;
	int n,E;
	n=E=0;
	printf("reading edge list file and converting to Weighted adjacency matrix...\n");
    // alloc graph
	W_GRAPH* WG = w_graph_alloc(N_nodes);
	FILE* input=open_file("r",input_name);
	i=j=tij=0;
	int k=0;
	int g;
	int opt_self;
	opt_self = 0;
	char dummy[100];
	while(k<header)
	{
		 //printf("I am here\n"); fflush(stdout);
		 fgets(dummy, 100, input);
		 //printf("%s",dummy);
		 k++;
	}
	while (!feof(input))
	{		///we start reading
		g = fscanf(input, "%d %d %d\n", &i, &j, &tij);
		if(g!=3)
		{
			printf("Error reading edge list\n");
			abort();
		}else{
			//printf("%d %d %d %d %d\n",k,g,i,j,tij);
			//printf("%d %d %d %d\n",k,i,j,tij);fflush(stdout);
			if(tij>0)
			{
				if(opt_dir>0)
				{
					w_graph_add_multi_link(WG, N_nodes, i, j, tij);
				}else{
					w_graph_add_multi_link_undirected(WG, N_nodes, i, j, tij);
				}
				if(i==j) opt_self = 1; // self_loop flag
			}
			n+=tij;
			E++;
		}
  	}
	printf(" -- Total num of trips %i \t Total number of edges:%d. Selfloops? %d\n",n,E,opt_self);
	fclose(input);
	WG->opt_dir = opt_dir;
	WG->opt_self = 0;
	return WG;
}


/****************************************************************************
 * Allocation *
 ****************************************************************************/
W_GRAPH* w_graph_alloc(int N_nodes){
    int i;
    W_GRAPH* WG = (W_GRAPH*)malloc(sizeof(W_GRAPH));
    WG->node = (NODE*)malloc(sizeof(NODE)*N_nodes);
    WG->N_nodes = N_nodes;
    for(i=0;i<N_nodes;i++)
    {
        WG->node[i].idnum=i;
        //node[i].x=x[0][i];
        //node[i].y=x[1][i];
        WG->node[i].kin=0;
        WG->node[i].kout=0;
        WG->node[i].kout=0;
        WG->node[i].sin=0;
        WG->node[i].sout=0;
        WG->node[i].mem_in=1;
        WG->node[i].mem_out=1;
        WG->node[i].out=calloc(1,sizeof(int));
        WG->node[i].out[0]=-1;
        WG->node[i].in=calloc(1,sizeof(int));
        WG->node[i].in[0]=-1;
        WG->node[i].w_in=calloc(1,sizeof(int));
        WG->node[i].w_out=calloc(1,sizeof(int));
    }
    // not setting these values, so error jumps if user does not set them
    //WG->opt_dir = 1; // directed by defatult (values set for backwards compatibliity)
    //WG->opt_self = 1; // self loops accepted by default 
    return WG;
}

void w_graph_free_nodes(W_GRAPH* WG, int N_nodes){ // warning, does not free edges!
    int i;
    for(i=0;i<N_nodes;i++)
    {
        free(WG->node[i].out);
        free(WG->node[i].in);
		free(WG->node[i].w_out);
		free(WG->node[i].w_in);
    }
    free(WG->node);
    return;
}

/****************************************************************************
 * Add edges *
 ****************************************************************************/
/// warning: Self-loops accepted ///
void w_graph_add_multi_link(W_GRAPH * WG, int N_nodes, int origin, int dest, int weight){
	int neigh,dummy;
	if(weight > 0)
	{
		WG->node[origin].sout +=weight;
		WG->node[dest].sin	 +=weight;
		neigh = find_value_int(WG->node[origin].out, dest, WG->node[origin].kout);
		if(neigh<0)// new connection
		{
			WG->node[origin].kout++;
			if(WG->node[origin].kout > WG->node[origin].mem_out)
			{
				dummy = WG->node[origin].mem_out;
				WG->node[origin].mem_out = dummy*2;
				WG->node[origin].out		= safe_int_realloc(WG->node[origin].out, dummy, WG->node[origin].mem_out, -1);
				//safe_int_realloc(WG->node[origin].out,dummy,node[origin].mem_out);
				WG->node[origin].w_out	= safe_int_realloc(WG->node[origin].w_out, dummy, WG->node[origin].mem_out, -1);
				//safe_int_realloc(WG->node[origin].w_out,dummy,node[origin].mem_out);
			}
			WG->node[origin].out[WG->node[origin].kout-1] 	= dest;
			WG->node[origin].w_out[WG->node[origin].kout-1] 	= weight;
		}else{ // existing connection
			WG->node[origin].w_out[neigh]					+=	weight;
		}
		neigh = find_value_int(WG->node[dest].in, origin, WG->node[dest].kin);
		if(neigh<0)// new connection
		{
			WG->node[dest].kin++;
			if(WG->node[dest].kin > WG->node[dest].mem_in)
			{
				dummy 					= WG->node[dest].mem_in;
				WG->node[dest].mem_in 	= 2*dummy;
				WG->node[dest].in 		= safe_int_realloc(WG->node[dest].in, dummy, WG->node[dest].mem_in, -1);
				//safe_int_realloc(WG->node[dest].in,dummy,node[dest].mem_in);
				WG->node[dest].w_in 		= safe_int_realloc(WG->node[dest].w_in, dummy, WG->node[dest].mem_in, -1);
				//safe_int_realloc(WG->node[dest].w_in,dummy,node[dest].mem_in);
			}
			WG->node[dest].in[WG->node[dest].kin-1]	= origin;
			WG->node[dest].w_in[WG->node[dest].kin-1]	= weight;
		}else{ // existing connection
			WG->node[dest].w_in[neigh] 				+= weight;
		}
	}
    return;
}


void w_graph_add_multi_link_undirected(W_GRAPH * WG, int N_nodes, int origin, int dest, int weight){
    if(weight > 0)
    {
		// keep sin-sout for compatibility issues and add twice if self_loop
		// update s
		WG->node[origin].sout	+= weight;
		WG->node[origin].sin		+= weight;
		if(dest!=origin)
		{
			WG->node[dest].sout 	+= weight;
			WG->node[dest].sin 	+= weight;
		}
		// check new connection
		int neigh,dummy;
		neigh = find_value_int(WG->node[origin].out, dest, WG->node[origin].kout);
		if(neigh<0)// new connection
		{
		// update k
			WG->node[origin].kout++;
			WG->node[origin].kin++;
		// if not enough space, allocate more
			if(WG->node[origin].kout>WG->node[origin].mem_out)
			{
				dummy 					= WG->node[origin].mem_out;
				WG->node[origin].mem_out	= 2*dummy;
				WG->node[origin].out		= safe_int_realloc(WG->node[origin].out, dummy, WG->node[origin].mem_out, -1);
				//safe_int_realloc(WG->node[origin].out,dummy,node[origin].mem_out);
				//printf("lalal am here , node: %d, dest: %d, k:%d, mem:%d\n", origin, dest, node[origin].kout, node[origin].mem_out); 				fflush(stdout);
				//safe_int_realloc(WG->node[origin].w_out,dummy,node[origin].mem_out);
				WG->node[origin].w_out	= safe_int_realloc(WG->node[origin].w_out, dummy, WG->node[origin].mem_out, -1);
			}
			// store neighbor
			WG->node[origin].out[WG->node[origin].kout-1] = dest;
			WG->node[origin].w_out[WG->node[origin].kout-1]		= weight;
			if(dest!=origin)
			{
				WG->node[dest].kout++;
				WG->node[dest].kin++;
				if(WG->node[dest].kout > WG->node[dest].mem_out)
				{
					dummy 					= WG->node[dest].mem_out;
					WG->node[dest].mem_out	= 2*dummy;
					WG->node[dest].out		= safe_int_realloc(WG->node[dest].out, dummy, WG->node[dest].mem_out,-1);
					//safe_int_realloc(WG->node[dest].out,dummy,node[dest].mem_out);
					WG->node[dest].w_out		= safe_int_realloc(WG->node[dest].w_out, dummy, WG->node[dest].mem_out,-1);
					//safe_int_realloc(WG->node[dest].w_out,dummy,node[dest].mem_out);
				}
				WG->node[dest].out[WG->node[dest].kout-1]   = origin;
				WG->node[dest].w_out[WG->node[dest].kout-1] = weight;
			}
		}else{ // existing connection
			WG->node[origin].w_out[neigh]   += weight;
			if(origin != dest)
			{
				neigh = find_value_int(WG->node[dest].out, origin, WG->node[dest].kout);
				assert(neigh>=0);
				WG->node[dest].w_out[neigh] += weight;
			}
		}
	}
    return;
}

/****************************************************************************
 * Printing *
 ****************************************************************************/
void w_graph_print_adj_list(W_GRAPH* WG, int N_nodes, char* output){
    int i,j;
    FILE* fil=open_file("w", output);
    fprintf(fil,"## Adjancenncy list (directed weighted format): Node id_source Node_id target weight \n");fflush(stdout);
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<WG->node[i].kout;j++)
        {
            fprintf(fil,"%d %d %d\n",WG->node[i].idnum, WG->node[WG->node[i].out[j]].idnum, WG->node[i].w_out[j]);
        }
    }
    fclose(fil);    
    return;
}
/****************************************************************************
 * S and K's *
 ****************************************************************************/

int ** w_graph_compute_s(W_GRAPH* WG, int N_nodes){
    int** s=cast_mat_int(2,N_nodes);
    int i;
    for(i=0;i<N_nodes;i++)
    {
        s[0][i] = WG->node[i].sout;
        s[1][i] = WG->node[i].sin;
    }
    return s;
}


int ** w_graph_compute_k(W_GRAPH* WG, int N_nodes){
    int** k=cast_mat_int(2,N_nodes);
    int i;
    for(i=0;i<N_nodes;i++)
    {
        k[0][i]=WG->node[i].kout;
        k[1][i]=WG->node[i].kin;
    }
    return k;
}


double ** w_graph_compute_k_analitic(W_GRAPH* WG, int N_nodes, int self_opt){
    double** k=cast_mat_double(2,N_nodes);
    int **s=w_graph_compute_s(WG, N_nodes);
    long int T=sum_vec_int(s[0],N_nodes);
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
                k[0][i]-=exp(-(double)WG->node[i].sout*((double)WG->node[j].sin/(long double)T));
                k[1][i]-=exp(-(double)WG->node[i].sin*((double)WG->node[j].sout/(long double)T));
            }else{
            	if(self_opt>0)
				{
	                k[0][i]-=exp(-(double)WG->node[i].sout*((double)WG->node[j].sin/(long double)T));
	                k[1][i]-=exp(-(double)WG->node[i].sin*((double)WG->node[i].sout/(long double)T));					
				}
            }
        }
    }
    free_mat_int(s,2);
    return k;
}


double * w_graph_compute_k_binary_undirected(double* x, int N_nodes, int self_opt){
    int i,j;
    double* k = cast_vec_double(N_nodes);
    for(i=0;i<N_nodes;i++)
    {
		k[i] = 0;
		for(j=0;j<N_nodes;j++)
		{
			if(i!=j || self_opt>0)
			{
				k[i]+= x[j] / (1.+x[i]*x[j]);
			}
		}
		k[i] = k[i]*x[i];
	}
	return k;
}

double w_graph_compute_E_binary_undirected(double* x, int N_nodes, int self_opt){
	double* kk = w_graph_compute_k_binary_undirected(x, N_nodes, self_opt);
	double E = sum_vec_double(kk,N_nodes);
	free(kk);
	return E;
}

double ** w_graph_compute_k_binary_directed(double** x, int N_nodes, int self_opt){
    int i,j;
    double** k = cast_mat_double(2,N_nodes);
    for(i=0;i<N_nodes;i++)
    {
		k[0][i] = 0;
		k[1][i] = 0;
		for(j=0;j<N_nodes;j++)
		{
			if(i!=j || self_opt>0)
			{
				k[0][i]+= x[1][j] / (1.+x[0][i]*x[1][j]);
				k[1][i]+= x[0][j] / (1.+x[1][i]*x[0][j]);
			}
		}
		k[0][i] = k[0][i]*x[0][i];
		k[1][i] = k[1][i]*x[1][i];
	}
	return k;
}

double w_graph_compute_E_binary_directed(double** x, int N_nodes, int self_opt){
	double** kk = w_graph_compute_k_binary_directed(x, N_nodes, self_opt);
	double E = sum_vec_double(kk[0],N_nodes);
	//printf("E:%f\n",E);fflush(stdout);
	free_mat_double(kk,2);
	return E;
}



double ** w_graph_compute_k_analitic_from_s_directed(int** s, int N_nodes, int self_opt){
    double** k=cast_mat_double(4,N_nodes); // 2 firts out-in, then sigma.
    long int T=(long int)sum_vec_int(s[0],N_nodes);
    int i,j;
    for(i=0;i<N_nodes;i++)
    {
        k[0][i]=(double)N_nodes;
        k[1][i]=(double)N_nodes;
        k[2][i]=(double)N_nodes;
        k[3][i]=(double)N_nodes;
        if(self_opt<=0)
        {
            k[0][i]-=1;
            k[1][i]-=1;
            k[2][i]-=1;
            k[3][i]-=1;
        }
        for(j=0;j<N_nodes;j++)
        {
            if(j!=i)
            {
                k[0][i]-=exp(-(double)s[0][i]*((double)s[1][j]/(long double)T));
                k[1][i]-=exp(-(double)s[1][i]*((double)s[0][j]/(long double)T));
				k[2][i]-=exp(-2*(double)s[0][i]*((double)s[1][j]/(long double)T));
				k[3][i]-=exp(-2*(double)s[1][i]*((double)s[0][j]/(long double)T));
            }else{
				if(self_opt>0)
				{
		    		k[0][i]-=exp(-(double)s[0][i]*((double)s[1][i]/(long double)T));
		    		k[1][i]-=exp(-(double)s[1][i]*((double)s[0][i]/(long double)T));
		    		k[2][i]-=exp(-2*(double)s[0][i]*((double)s[1][i]/(long double)T));
		    		k[3][i]-=exp(-2*(double)s[1][i]*((double)s[0][i]/(long double)T));

				}
	    	}
        }
		k[2][i] = k[2][i]-k[0][i]; // sigma = N- k - exp 2!
		k[3][i] = k[3][i]-k[1][i];
    }
    return k;
}

double ** w_graph_compute_k_analitic_from_s_undirected(int* s, int N_nodes, int self_opt){
    double** k=cast_mat_double(2,N_nodes);
    long int T=(long int)sum_vec_int(s,N_nodes);
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
                k[0][i]-=exp(-(double)s[i]*((double)s[j]/(long double)T));
                k[1][i]-=exp(-2*(double)s[i]*((double)s[j]/(long double)T));

            }else{
				if(self_opt>0)
				{
		    		k[0][i]-=exp(-(double)s[i]*((double)s[j]/(long double)T));
		    		k[1][i]-=exp(-2*(double)s[i]*((double)s[j]/(long double)T));
				}
	    	}
        }
		k[1][i] = k[1][i]-k[0][i];
    }
    return k;
}


double w_graph_compute_T_analytic_from_xy(double** x2, int N_nodes, int layers, double gamma, int opt_dir, int opt_self, int opt_indist, int opt_agg, int fixk){
    int i,j,ind2,ind3,ind4;
    double t,T;
    // if fixk=0 consider the case where s is fixed, fixk=1 case with fixed number of edges, fixk=2 s and k are fixed
    T=0;
    if(opt_dir>0)
    {
        ind2=1;
        ind3=2;
        ind4=3;
    }else{
        ind2=0;
        ind3=1;
        ind4=1;
    }
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<N_nodes;j++)
        {
            if((opt_indist<=0)&&(opt_agg<=0))
            { // ME
                if(fixk<=0)
                {
                    t = x2[0][i]*x2[ind2][j];                
                }else if(fixk==1){
                    t = exp(x2[0][i]*x2[ind2][j])*x2[0][i]*x2[ind2][j]*gamma/(1.+gamma*(exp(x2[0][i]*x2[ind2][j])-1));
                }else{
                    t = exp(x2[0][i]*x2[ind2][j])*x2[ind3][i]*x2[ind4][j]/(1.+x2[ind3][i]*x2[ind4][j]*(exp(x2[0][i]*x2[ind2][j])-1));
                }
            }else{
                if(opt_indist>0)
                { // W
                    if(fixk<=0)
                    {
                        t = layers*x2[0][i]*x2[ind2][j]/(1.-x2[0][i]*x2[ind2][j]);
                    }else if(fixk==1){
                        t = layers/pow(1.-x2[0][i]*x2[ind2][j],layers)*x2[0][i]*x2[ind2][j]*gamma/(1.+gamma*(1./pow(1.-x2[0][i]*x2[ind2][j],layers)-1));
                    }else{
                        t = layers/pow(1.-x2[0][i]*x2[ind2][j],layers)*x2[ind3][i]*x2[ind4][j]/(1.+x2[ind3][i]*x2[ind4][j]*(1./pow(1.-x2[0][i]*x2[ind2][j],layers)-1));
                    }                    
                }else{ // Binary
                    if(fixk<=0)
                    {
                        t = layers*x2[0][i]*x2[ind2][j]/(1.+x2[0][i]*x2[ind2][j]);                    
                    }else if(fixk==1){
                        t = layers*pow(1.+x2[0][i]*x2[ind2][j],layers)*x2[0][i]*x2[ind2][j]*gamma/(1.+gamma*(pow(1.+x2[0][i]*x2[ind2][j],layers)-1));
                    }else{
                        t = layers*pow(1.+x2[0][i]*x2[ind2][j],layers)*x2[ind3][i]*x2[ind4][j]/(1.+x2[ind3][i]*x2[ind4][j]*(pow(1.+x2[0][i]*x2[ind2][j],layers)-1));
                    }                    
                }
            }
            T+=t;
        }    
    }
    return T;
}



/****************************************************************************
 * Snn and Knn's *
 ****************************************************************************/
double ** w_graph_compute_s_nn(W_GRAPH* WG, int N_nodes, int weight, int opt_dir){
    double ** s_nn=cast_mat_double(4, N_nodes);
    int i,j;
    //FILE* prova=open_file("a","check.dat");
    for(i=0;i<N_nodes;i++)
    {
        s_nn[0][i]=s_nn[1][i]=s_nn[2][i]=s_nn[3][i]=0;
        if(WG->node[i].kout>0)
	{
	for(j=0;j<WG->node[i].kout;j++)
        {
            if(weight>0)
            {
                s_nn[1][i]+=(double)(WG->node[i].w_out[j])*(double)WG->node[WG->node[i].out[j]].sin; // average out-neigh, in strength
                s_nn[0][i]+=(double)(WG->node[i].w_out[j])*(double)WG->node[WG->node[i].out[j]].sout; // average out-neigh, out
            }else if(weight==0){
                s_nn[1][i]+=(double)(WG->node[i].w_out[j])*(double)WG->node[WG->node[i].out[j]].kin; // weighted degree
                s_nn[0][i]+=(double)(WG->node[i].w_out[j])*(double)WG->node[WG->node[i].out[j]].kout; // weighted degree
            }else{
                s_nn[1][i]+=(double)WG->node[WG->node[i].out[j]].kin; // average out-neigh, in strength
                s_nn[0][i]+=(double)WG->node[WG->node[i].out[j]].kout; // average out-neigh, out                                
            }
        }
    	}
	if((WG->node[i].kin>0) && (opt_dir>0))
	{
        for(j=0;j<WG->node[i].kin;j++)
        {        
            if(weight>0)
            {
                s_nn[3][i]+=WG->node[i].w_in[j]*(double)WG->node[WG->node[i].in[j]].sin; // average in-neigh, in
                s_nn[2][i]+=WG->node[i].w_in[j]*(double)WG->node[WG->node[i].in[j]].sout; // average in-neigh, out
            }else if(weight==0){
                s_nn[3][i]+=WG->node[i].w_in[j]*(double)WG->node[WG->node[i].in[j]].kin; // average in-neigh, in
                s_nn[2][i]+=WG->node[i].w_in[j]*(double)WG->node[WG->node[i].in[j]].kout; // average in-neigh, out
            }else{
                s_nn[3][i]+=(double)WG->node[WG->node[i].in[j]].kin; // average out-neigh, in strength
                s_nn[2][i]+=(double)WG->node[WG->node[i].in[j]].kout; // average out-neigh, out                
            }

        }
	}
        if(weight>=0)
        {
	    if(WG->node[i].kout>0)
	    {
            s_nn[0][i]/=(double)WG->node[i].sout;
            s_nn[1][i]/=(double)WG->node[i].sout;
	    }
	    if((WG->node[i].kin>0) && (opt_dir>0))
	    {
	    	s_nn[2][i]/=(double)WG->node[i].sin;
            s_nn[3][i]/=(double)WG->node[i].sin;
	    }
        //fprintf(prova,"%d %d %lf\n",i,node[i].sout,s_nn[0][i]);
        }else{
	    if(WG->node[i].kout>0)
	    {
            s_nn[0][i]/=(double)WG->node[i].kout;
            s_nn[1][i]/=(double)WG->node[i].kout;
	    }
	    if((WG->node[i].kin>0) && (opt_dir>0))
	    {
            s_nn[2][i]/=(double)WG->node[i].kin;
            s_nn[3][i]/=(double)WG->node[i].kin;                
	    }
        }
    }
    //fclose(prova);
    return s_nn;
}

/****************************************************************************
 * Clustering *
 ****************************************************************************/

double ** w_graph_compute_clust(W_GRAPH * WG, int N_nodes){ // 2 cols: unweighted, weighted
    //printf("#### Computing c ####\n"); fflush(stdout);
    int i,j,k,l;
    int c_glo, c_glo_den,e;
    double ee;
    c_glo=c_glo_den=0;
    double** c=cast_mat_double(2,N_nodes);
    ///for per cada node i
    for(i=0;i<N_nodes;++i)
    {
	if(WG->node[i].kout>0)
	{
    //printf("start, %d, %d \n",i,node[i].kin);fflush(stdout);
    	e=ee=0;
	//printf("node=%i\n",i);
	///for per cada vei j
	for(j=0;j<WG->node[i].kout;++j)
	{
	    //printf("vei=%i amb degree=%i\n",node[i].out[j],node[WG->node[i].out[j]].k);
	    ///for per els altres nodes veins de i. comencem a j+1 per no repetir
	    for(k=j+1;k<WG->node[i].kout;++k)
	    {
		//printf("segon_vei=%i\n",node[i].out[k]);
		///un for per veure els veins del node j i veure si esta conectat amb el node k
		for (l=0;l<WG->node[WG->node[i].out[j]].kout;++l)
		{
		    //printf("l=%i\n",l);
		    if(WG->node[WG->node[i].out[j]].out[l]==WG->node[i].out[k]) 
		    {
			e=e+1;
			ee+=(WG->node[i].w_out[j]+WG->node[i].w_out[k]);
		    } 
		    //printf("   node=%i vei=%i veivei=%i veicomparar=%i  e=%i\n",i,node[i].out[j],node[WG->node[i].out[j]].out[l],node[i].out[k],e);
		}
	    }
	}                        
	if(e!=0)
        {
	    c_glo=c_glo+e; // triangles
	    c_glo_den=c_glo_den+WG->node[i].kout*(WG->node[i].kout-1)/2; // pairs of neighbours
	    //cc=cc+(double)e*2./node[i].k/(WG->node[i].k-1); // l
	    c[0][i]= (double)e*2./WG->node[i].kout/(WG->node[i].kout-1); // clustering
	    c[1][i]= (double)ee/(2.*WG->node[i].sout)/(WG->node[i].kout-1); // weighted clust
	}
    	}
        //printf("stop, %d\n",i);fflush(stdout);
    }
    ///and we take the mean
    //cc=cc/(n-number_k[1]);
    //printf("global clustering= %f\n",(double)c_glo/(double)c_glo_den);
    return c;
}



/****************************************************************************
 * Weight funcs *
 ****************************************************************************/

int w_graph_total_weight( W_GRAPH* WG, int N_nodes){
    int i;
    int T;
    T=0;
    for(i=0;i<N_nodes;i++)
    {
        T+=(int)WG->node[i].sout;
    }
    WG->T = T;
    return T;
}

int w_graph_total_edges( W_GRAPH* WG, int N_nodes){
    // counts self loops properly (bit slower... but ok)
    int i,j;
    //int T,aux2;
    int aux2;
    aux2=0;
    for(i=0;i<N_nodes;i++)
    {
        //T+=WG->node[i].kout;
        for(j=0;j<WG->node[i].kout;j++)
        {
			aux2++;
			if(WG->node[i].out[j]==i) // if selfs count twice
			{
				aux2++;
			}
		}
    }
    //return T;
    WG->E = aux2;
    return aux2;
}

int w_graph_total_edgepairs( W_GRAPH* WG, int N_nodes){
    // counts total number of edge pairs
    int L;
    if (WG->opt_self>0)
    {
		L = N_nodes*N_nodes+N_nodes; // we coutn self-loops twice to adapt to undirected

	}else{
		L = N_nodes*(N_nodes-1);
	}
	WG->L = L;
    return L;
}

int * w_graph_compute_p(W_GRAPH* WG, int N_nodes, int* aux){
    int E;
    int* w;
    int i,j,aux2,mem,t,dest;
    aux2=0;
    
	E = w_graph_total_edgepairs( WG, N_nodes);
	w = cast_vec_int(E);
	mem=E;
	for(i=0;i<N_nodes;i++) // all edges
	{
		for(j=0;j<N_nodes;j++)
		{
			if((WG->opt_self>0)||(i!=j))
			{
				dest=find_value_int(WG->node[i].out, j, WG->node[i].kout);
				if(dest>=0)
				{
					t = 1;
				}else{
					t = 0;
				}
				w[aux2] = (int)t;
				aux2++;
				if(j==i) // count twice
				{
					w[aux2]=(int)t;
					aux2++;
				}
				//printf("Edgepair: %d. I:%d J:%d | Connection: %d Dest:%d \n",aux2,i,j,t,dest);fflush(stdout);
			}
		}
	}
	assert(aux2==E);
	*aux=aux2;
    return w;
}


int * w_graph_compute_w(W_GRAPH* WG, int N_nodes, int* aux, int zeros){
    //zeros not implemented
    int E;
    int* w;
    int i,j,aux2,mem,t,dest;
    aux2=0;

	int L2;
    L2 =  w_graph_total_edgepairs(WG, N_nodes);
    
    if(zeros>0) // add also zeros
    {
		E = w_graph_total_edgepairs( WG, N_nodes);
		w = cast_vec_int(E);
		mem=E;
		for(i=0;i<N_nodes;i++) // all edges
		{
			for(j=0;j<N_nodes;j++)
			{
				if((WG->opt_self>0)||(i!=j))
				{
					dest=find_value_int(WG->node[i].out, j, WG->node[i].kout);
					if(dest>=0)
					{
						t = WG->node[i].w_out[dest];
					}else{
						t = 0;
					}
					w[aux2] = (int)t;
					aux2++;
					if(j==i) // count twice
					{
						w[aux2]=(int)t;
						aux2++;
					}
					//printf("Edgepair: %d. I:%d J:%d | Connection: %d Dest:%d \n",aux2,i,j,t,dest);fflush(stdout);
				}
			}
		}
		assert(aux2==L2);
    }else{
		E = w_graph_total_edges( WG, N_nodes);
		w = cast_vec_int(E);
		mem=E;
        for(i=0;i<N_nodes;i++)
        {
            for(j=0;j<WG->node[i].kout;j++)
            {
                if(aux2+1>mem)
                {
					w=safe_int_realloc(w, mem, 2*mem, 0);
                    mem=2*mem;
                }
                w[aux2]=(int)WG->node[i].w_out[j];
                aux2++;
                if(WG->node[i].out[j]==i)
                {
                    if(aux2+1>mem)
                    {
						w=safe_int_realloc(w, mem, 2*mem, 0);
                        mem=2*mem;
                    }
                    w[aux2]=(int)WG->node[i].w_out[j];
                    aux2++;
                }
                
            }
        }
		w=safe_int_realloc(w, mem, aux2, 0);// fix final size
    }
    //printf("E: %d, aux2:%d delta:%d\n", E, aux2, E-aux2); fflush(stdout);
    *aux=aux2;
    return w;
}

double** w_graph_compute_p_w_analitic_from_s_undirected(int maxt, double binn, int* s, int N_nodes, int self_opt, int* len){
    // computes already normalized p
    printf("Computing w, tmax:%d\n",maxt); fflush(stdout);
    int i,j,t,aux;
    int T=sum_vec_int(s,N_nodes);
    double norm,mu;
    double** p=cast_mat_double(2,N_nodes);
    mu=0;
    norm=0;
    aux=0;
    t=1;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<i;j++)
        {
            norm+=1.-(double)exp(-(double)s[i]*(double)s[j]/(double)T);
        }
        if(self_opt>0)
        {
            norm+=1.-(double)exp(-(double)s[i]*(double)s[i]/(double)T);
        }
    }
    while(t<maxt+1)
    {
        p[0][aux]=t;
        p[1][aux]=0;
        //fact=factorial(t);
        for(i=0;i<N_nodes;i++)
        {
            for(j=0;j<i;j++)
            {
                mu = (double)s[i]*(double)s[j]/(double)T;
                p[1][aux]+=gsl_ran_poisson_pdf (t, mu);
                if(t==0)norm+=(double)exp(-mu);
            }
            if(self_opt>0)
            {
                mu = (double)s[i]*(double)s[i]/(double)T;
                if(t==0)norm+=(double)exp(-mu);
                p[1][aux]+=gsl_ran_poisson_pdf (t, mu);
            }
        }
        p[1][aux]=p[1][aux]/norm;
        //printf("done %d\n",t); fflush(stdout);
        if(t<10){
            t++;
        }else{
            t=t*binn;
        }
        aux++;
    }
    (*len)=aux;
    printf("done: %d bins\n",aux); fflush(stdout);
    p[0] = safe_double_realloc(p[0],N_nodes,aux,0);
    p[1] = safe_double_realloc(p[1],N_nodes,aux,0);
    return p;
}

double** w_graph_compute_p_w_analitic_from_s_directed(int maxt, double binn, int** s, int N_nodes, int self_opt, int* len){
    // computes already normalized p
    //printf("Computeing w, tmax:%d\n",maxt); fflush(stdout);
    printf("Computing w, tmax:%d\n",maxt); fflush(stdout);
    int i,j,t,aux;
    int T=sum_vec_int(s[0],N_nodes);
    double norm,mu;
    double** p=cast_mat_double(2,N_nodes);
    mu=0;
    norm=(double)N_nodes;
    aux=0;
    t=1;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<N_nodes;j++)
        {
            mu = (double)s[0][i]*(double)s[1][j]/(double)T;
            if(i!=j)norm+=1.-(double)exp(-mu);
        }
        if(self_opt>0)
        {
            mu = (double)s[0][i]*(double)s[1][i]/(double)T;
           norm+=1.-(double)exp(-mu);
        }
    }
    while(t<maxt+1)
    {
        p[0][aux]=t;
        p[1][aux]=0;
        for(i=0;i<N_nodes;i++)
        {
            for(j=0;j<N_nodes;j++)
            {
                if(j!=i)
                {
                    mu = (double)s[0][i]*(double)s[1][j]/(double)T;
                    p[1][aux]+=gsl_ran_poisson_pdf (t, mu);
                }
            }
            if(self_opt>0)
            {
                mu = (double)s[0][i]*(double)s[1][j]/(double)T;
                p[1][aux]+=gsl_ran_poisson_pdf (t, mu);
            }
        }
        p[1][aux]=p[1][aux]/norm;
        //printf("done %d\n",t); fflush(stdout);
        if(t<10){
            t++;
        }else{
            t=t*binn;
        }
        aux++;
    }
    (*len)=aux;
    //printf("done\n"); fflush(stdout);
    p[0] = safe_double_realloc(p[0],N_nodes,aux,0);
    p[1] = safe_double_realloc(p[1],N_nodes,aux,0);
    return p;
}

double * w_graph_compute_wp_ss(W_GRAPH* WG, int N_nodes, int weight){
	// computes existing weight average //
    int E;
    E=w_graph_total_edges(WG,N_nodes);
    //int* ww=cast_vec_int(E);
    double* ww=cast_vec_double(E);
    int i,j,aux;
    aux=0;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<WG->node[i].kout;j++)
        {
            if (weight>0)
            {
                ww[aux]=((double)WG->node[i].sout)*((double)WG->node[WG->node[i].out[j]].sin);
            }else{
                ww[aux]=((double)WG->node[i].kout)*((double)WG->node[WG->node[i].out[j]].kin);
            }
            aux++;
            if(WG->node[i].out[j]==i) // if self loop, count twice
            {
                if (weight>0)
                {
                    ww[aux]=((double)WG->node[i].sout)*((double)WG->node[WG->node[i].out[j]].sin);
                }else{
                    ww[aux]=((double)WG->node[i].kout)*((double)WG->node[WG->node[i].out[j]].kin);
                }
                aux++;
            }
        }
    }
    assert(aux==E);
    return ww;
}

double * w_graph_compute_w_ss(W_GRAPH* WG, int N_nodes, int weight){
	// computes total weight average //
	int opt_self = WG->opt_self;
    int E;
    E=w_graph_total_edgepairs(WG,N_nodes);
    double* ww=cast_vec_double(E);
    int i,j,aux;
    aux=0;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<N_nodes;j++)
        {
			if((opt_self>0)||(i!=j))
			{
				if (weight>0)
				{
					ww[aux]=((double)WG->node[i].sout)*((double)WG->node[j].sin);
				}else{
					ww[aux]=((double)WG->node[i].kout)*((double)WG->node[j].kin);
				}
				aux++;
				if(i==j)
				{
					if (weight>0)
					{
						ww[aux]=((double)WG->node[i].sout)*((double)WG->node[j].sin);
					}else{
						ww[aux]=((double)WG->node[i].kout)*((double)WG->node[j].kin);
					}
					aux++;
				}
            }
        }
    }
    assert(aux==E);
    return ww;
}
//double ** w_graph_compute_xy(w_graph * node, int N_nodes){
    
    //int i,j;
    //double** y2=cast_mat_double(3,N_nodes);
    //for(i=0;i<N_nodes;i++)
    //{
	//if(WG->node[i].kout>0)
	//{
        //y2[0][i]=0;
	    //for(j=0;j<node[i].kout;j++)
	    //{
	    	//y2[0][i]+=((double)node[i].w_out[j])*((double)node[i].w_out[j]); // x
	    //}
	//y2[1][i]=node[i].sout*node[i].sout; // y
	//y2[2][i]=y2[1][i]*y2[0][i]; // xy
    	//}
    //}
    //return y2;
//}


double ** w_graph_compute_Y2(W_GRAPH * WG, int N_nodes, int opt_dir){
    
    int i,j;
    double** y2=cast_mat_double(2,N_nodes);
    for(i=0;i<N_nodes;i++)
    {
	y2[0][i]=0;
	if(WG->node[i].kout>0)
	{
	    for(j=0;j<WG->node[i].kout;j++)
	    {
	    	y2[0][i]+=((double)WG->node[i].w_out[j])*((double)WG->node[i].w_out[j]);
/*            if((WG->node[i].w_out[j]<1) )
            {
                printf("%f\n",(double)WG->node[i].w_out[j]);
            }
 */
	    }
	    y2[0][i]/=(((double)WG->node[i].sout)*((double)WG->node[i].sout));
    }
	y2[1][i]=0;
	if((WG->node[i].kin>0) && (opt_dir>0))
	{        
	    for(j=0;j<WG->node[i].kin;j++)
	    {
            	y2[1][i]+=((double)WG->node[i].w_in[j])*((double)WG->node[i].w_in[j]);            
	    }
	    y2[1][i]/=(((double)WG->node[i].sin)*((double)WG->node[i].sin));
    }
/*    if((y2[0][i] >0)|( y2[1][i] >0))
       {
           printf("out y2:%f \t s: %d \t k: %d\n in y2:%f \t s: %d \t k: %d\n\n",y2[0][i],WG.WG->node[i].sout,WG->node[i].kout,y2[1][i],WG->node[i].sin,WG->node[i].kin); fflush(stdout);
       }
*/
    }
    return y2;
}

/****************************************************************************
 * Entropy *
 ****************************************************************************/
double w_graph_entropy_multinomial(W_GRAPH* WG, int N_nodes, int opt_dir){ // not entirely correct (esactly)... but...
	double S,p;
	int i,j,t;
	int T = w_graph_total_weight(WG,N_nodes);
	S = 0;
	for(i=0;i<N_nodes;i++) // all edges
	{
		for(j=0;j<WG->node[i].kout;j++)
		{
			t = WG->node[i].w_out[j];
			p = (double)t/(double)T;
			S+= p*log(p);
		}		
	}
	if(opt_dir<=0) S=S/2.;
	return -S;
}


double w_graph_entropy_poisson(W_GRAPH* WG, double** x,int N_nodes, int opt_self, int opt_dir){
	double S,p,mu;
	int i,j,t;
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
                if(t>=0)
                {
					t = WG->node[i].w_out[t];
				}else{
					t = 0;
				}
                mu = y[i]*yy[j];
                if(mu>0)
				{
					p = gsl_ran_poisson_pdf (t, mu);
					//printf("Mu:%f p:%f t:%d",mu,p,t);fflush(stdout);
				}else{
					p = 1./0.; // infinity! (not compatible)
					if(t==0) p = 0;
				}
				if(p>0) S+= p*log(p);
			}
        }
	}
	if(opt_dir<=0) S = S/2.;
	return -S;    
}


double w_graph_entropy_geometric(W_GRAPH* WG, double**x, int N_nodes, int opt_self, int opt_dir){
	double S,p,mu;
	int i,j,t;
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
                if(t>=0)
                {
					t = WG->node[i].w_out[t];
				}else{
					t = 0;
				}
                mu = y[i]*yy[j]; // for gsl this is 1-p
                if(mu>0)
				{
					p = gsl_ran_geometric_pdf (t, 1.-mu)*mu; // returns p (1-p)^(k-1) [so must multiply by (1-p)]
					//printf("Mu:%f p:%f t:%d",mu,p,t);fflush(stdout);
				}else{
					p = 1./0.; // infinity! (not compatible)
					if(t==0) p = 0;
				}
				if(p>0) S+= p*log(p);
			}
        }
	}
	if(opt_dir<=0) S = S/2.;
	return -S;    
}

double w_graph_entropy_binomial(W_GRAPH* WG, double**x, int N_nodes, int layers, int opt_self, int opt_dir){
	double S,p,mu;
	int i,j,t;
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
                if(t>=0)
                {
					t = WG->node[i].w_out[t];
				}else{
					t = 0;
				}
                mu = y[i]*yy[j]/(1+y[i]*yy[j]);
                if(mu>0)
				{
					gsl_ran_binomial_pdf(t,mu,layers);
					//printf("Mu:%f p:%f t:%d",mu,p,t);fflush(stdout);
				}else{
					p = 1./0.; // infinity! (not compatible)
					if(t==0) p = 0;
				}
				if(p>0) S+= p*log(p);
			}
        }
	}
	if(opt_dir<=0) S = S/2.;
	return -S;    
}

double w_graph_entropy_negbinomial(W_GRAPH* WG, double**x, int N_nodes, int layers, int opt_self, int opt_dir){
	double S,p,mu;
	int i,j,t;
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
                if(t>=0)
                {
					t = WG->node[i].w_out[t];
				}else{
					t = 0;
				}
                mu = y[i]*yy[j]/(1.+y[i]*yy[j]);
                if(mu>0)
				{
					gsl_ran_negative_binomial_pdf(t,1.-mu,layers); // for gsl p is inverted
					//printf("Mu:%f p:%f t:%d",mu,p,t);fflush(stdout);
				}else{
					p = 1./0.; // infinity! (not compatible)
					if(t==0) p = 0;
				}
				if(p>0) S+= p*log(p);
			}
        }
	}
	if(opt_dir<=0) S = S/2.;
	return -S;    
}
double w_graph_entropy_bernouilli(W_GRAPH* WG, double** x, int N_nodes, int opt_self, int opt_dir){
	double S,p,mu;
	int i,j,t;
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
                if(t>=0)
                {
					t = WG->node[i].w_out[t];
				}else{
					t = 0;
				}
                mu = y[i]*yy[j]/(1+y[i]*yy[j]);
                if(mu>0)
				{
					p = mu;
					//printf("Mu:%f p:%f t:%d",mu,p,t);fflush(stdout);
				}else{
					p = 1./0.; // infinity! (not compatible)
					if(t==0) p = 0;
				}
				if(p>0) S+= p*log(p);
			}
        }
	}
	if(opt_dir<=0) S = S/2.;
	return -S;    
}


double w_graph_entropy_ZIP(W_GRAPH* WG, double** x, int N_nodes, double gamma, int opt_self, int opt_dir){
	double S,p;
	int i,j,t,dest;
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
				t=find_value_int(WG->node[i].out, j, WG->node[i].kout);
                if(t>-1) // connection   
                {
					dest = WG->node[i].out[t];
					t = WG->node[i].w_out[dest];
					p = gamma*exp(y[i]*yy[dest])/(gamma*(exp(y[i]*yy[dest])-1)+1) * gsl_ran_poisson_pdf (t, y[i]*yy[dest]);//
				}else{   
                    p = 1 - gamma*(exp(y[i]*yy[j])-1)/(gamma*(exp(y[i]*yy[j])-1)+1);
				}
            }		
		    S+= p*log(p);
        }
	}
	if(opt_dir<0) S =S/2.;
	return -S;    
}



double w_graph_entropy_ZIP2(W_GRAPH* WG, double** x,int N_nodes,  int opt_self, int opt_dir){
	double S,p;
	int i,j,t,dest;
    double *y;
    double *yy;
    double * z;
    double * zz;
	S = 0;
    if(opt_dir>0)
    {
        y = x[0];
        yy = x[1];
        z = x[2];
        zz = x[3];
    }else{
        y = x[0];
        yy = x[0];
        z = x[1];
        zz = x[1];
    }
	for(i=0;i<N_nodes;i++) // all edges
	{
 		for(j=0;j<N_nodes;j++)
		{
            if((opt_self>0)||(i!=j))
            {
				t=find_value_int(WG->node[i].out, j, WG->node[i].kout);
                if(t>-1) // connection   
                {
					dest = WG->node[i].out[t];
					t = WG->node[i].w_out[dest];
					p = z[i]*zz[dest]*exp(y[i]*yy[dest])/(z[i]*zz[dest]*(exp(y[i]*yy[dest])-1)+1) * gsl_ran_poisson_pdf (t, y[i]*yy[dest]);//
				}else{   
                    p = 1 - z[i]*zz[j]*(exp(y[i]*yy[j])-1)/(z[i]*zz[j]*(exp(y[i]*yy[j])-1)+1);
				}
            }		
		    S+= p*log(p);
        }
	}
	if(opt_dir<0) S =S/2.;
	return -S;    
}


double w_graph_entropy_ZIB2(W_GRAPH* WG, double** x, int N_nodes, int layers, int opt_self, int opt_dir){
	double S,p;
	int i,j,t,dest;
    double *y;
    double *yy;
    double * z;
    double * zz;
	S = 0;
    if(opt_dir>0)
    {
        y = x[0];
        yy = x[1];
        z = x[2];
        zz = x[3];
    }else{
        y = x[0];
        yy = x[0];
        z = x[1];
        zz = x[1];
    }
	for(i=0;i<N_nodes;i++) // all edges
	{
 		for(j=0;j<N_nodes;j++)
		{
            if((opt_self>0)||(i!=j))
            {
				t=find_value_int(WG->node[i].out, j, WG->node[i].kout);
                if(t>-1) // connection   
                {
					dest = WG->node[i].out[t];
					t = WG->node[i].w_out[dest];
					p = z[i]*zz[dest]*pow(1+y[i]*yy[dest],layers)/(z[i]*zz[dest]*(pow(1+y[i]*yy[dest],layers)-1)+1) * gsl_ran_binomial_pdf (t, y[i]*yy[dest]/(1+y[i]*yy[dest]),layers);//
				}else{   
                    p = 1 - z[i]*zz[j]*(pow(1+y[i]*yy[j],layers)-1)/(z[i]*zz[j]*(pow(1+y[i]*yy[j],layers)-1)+1);
				}
            }		
		    S+= p*log(p);
        }
	}
	if(opt_dir<0) S =S/2.;
	return -S;  
}



void w_graph_print_entropy(double* seq,int len,char* output){
	int bins;
	double mins,maxs,eps;
	bins = (int)len/20.;
	if(bins<5)
	{
		bins=5;
	}else if(bins>50){
		bins=50;
	}
	mins = min_value_double(seq,len);
	maxs = max_value_double(seq,len);
	eps = (maxs-mins)/1000.;
    mins = mins - mins*eps;
    maxs = maxs+maxs*eps;
    if(mins<0) mins = 0;
    if(maxs<0) maxs = 0;
    if(mins>=maxs) maxs = mins+1;
    gsl_histogram* h1=histogram_double(seq,mins,maxs,bins,len);
	print_acc(output,h1,h1);
	return;
}


/****************************************************************************
 *LIkelyhood funcs *
 ****************************************************************************/
double w_graph_loglikelyhood_poisson(W_GRAPH* WG,int N_nodes,double** wij, int opt_self, int opt_dir){
	double L,p,mu;
	int i,j,t;
	L = 0;
	for(i=0;i<N_nodes;i++) // all edges
	{
		for(j=0;j<N_nodes;j++)
		{
            if((opt_self>0)||(i!=j))
            {
                t=find_value_int(WG->node[i].out, j, WG->node[i].kout);
                if(t>=0)
                {
					t = WG->node[i].w_out[t];
				}else{
					t = 0;
				}
                mu = wij[i][j];
                if(mu>0)
				{
					p = gsl_ran_poisson_pdf (t, mu);
					//printf("Mu:%f p:%f t:%d",mu,p,t);fflush(stdout);
				}else{
					p = 1./0.; // infinity! (not compatible)
					if(t==0) p = 0;
				}				
				if(p>0) L+= log(p);
			}
		}
	}
	if(opt_dir<=0) L = L/2.;
	return L;    
}


double w_graph_loglikelyhood_multinomial(W_GRAPH* WG, int N_nodes, int opt_dir){ // not entirely correct (esactly)... but...
	double S,p;
	int i,j,t;
	int T = w_graph_total_weight(WG,N_nodes);
	S = 0;
	for(i=0;i<N_nodes;i++) // all edges
	{
		for(j=0;j<WG->node[i].kout;j++)
		{
			t = WG->node[i].w_out[j];
			p = (double)t/(double)T;
			S+= log(p);
		}		
	}
	if(opt_dir<=0) S=S/2.;
	return -S;
}


double w_graph_loglikelyhood_poisson_xy(W_GRAPH* WG, double** x,int N_nodes, int opt_self, int opt_dir){
	double S,p,mu;
	int i,j,t;
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
                if(t>=0)
                {
					t = WG->node[i].w_out[t];
				}else{
					t = 0;
				}
                mu = y[i]*yy[j];
                if(mu>0)
				{
					p = gsl_ran_poisson_pdf (t, mu);
					//printf("Mu:%f p:%f t:%d",mu,p,t);fflush(stdout);
				}else{
					p = 1./0.; // infinity! (not compatible)
					if(t==0) p = 0;
				}
				if(p>0)S+= log(p);
			}
        }
	}
	if(opt_dir<=0) S = S/2.;
	return -S;    
}


double w_graph_loglikelyhood_geometric_xy(W_GRAPH* WG, double**x, int N_nodes, int opt_self, int opt_dir){
	double S,p,mu;
	int i,j,t;
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
                if(t>=0)
                {
					t = WG->node[i].w_out[t];
				}else{
					t = 0;
				}
                mu = y[i]*yy[j]; // for gsl this is 1-p
                if(mu>0)
				{
					p = gsl_ran_geometric_pdf (t+1, 1.-mu); // returns p (1-p)^(k-1) [so must multiply by (1-p)]
					//printf("Mu:%f p:%f t:%d",mu,p,t);fflush(stdout);
				}else{
					p = 1./0.; // infinity! (not compatible)
					if(t==0) p = 0; // prob is 1, but it will not count anyway
				}
				if(p>0) S+= log(p);
			}
        }
	}
	if(opt_dir<=0) S = S/2.;
	return -S;    
}

double w_graph_loglikelyhood_binomial_xy(W_GRAPH* WG, double**x, int N_nodes, int layers, int opt_self, int opt_dir){
	double S,p,mu;
	int i,j,t;
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
                if(t>=0)
                {
					t = WG->node[i].w_out[t];
				}else{
					t = 0;
				}
                mu = y[i]*yy[j]/(1+y[i]*yy[j]);
                if(mu>0)
				{
					p = gsl_ran_binomial_pdf(t,mu,layers);
					//printf("Mu:%f p:%f t:%d",mu,p,t);fflush(stdout);
				}else{
					p = 1./0.; // infinity! (not compatible)
					if(t==0) p = 0;
				}
				if(p>0) S+= log(p);
			}
        }
	}
	if(opt_dir<=0) S = S/2.;
	return -S;    
}

double w_graph_loglikelyhood_negbinomial_xy(W_GRAPH* WG, double**x, int N_nodes, int layers, int opt_self, int opt_dir){
	double S,p,mu;
	int i,j,t;
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
                if(t>=0)
                {
					t = WG->node[i].w_out[t];
				}else{
					t = 0;
				}
                mu = y[i]*yy[j]/(1.+y[i]*yy[j]);
                if(mu>0)
				{
					p = gsl_ran_negative_binomial_pdf(t,1.-mu,layers); // for gsl p is inverted
					//printf("Mu:%f p:%f t:%d",mu,p,t);fflush(stdout);
				}else{
					p = 1./0.; // infinity! (not compatible)
					if(t==0) p = 0;
				}
				if(p>0) S+= log(p);
			}
        }
	}
	if(opt_dir<=0) S = S/2.;
	return -S;    
}
double w_graph_loglikelyhood_bernouilli_xy(W_GRAPH* WG, double** x, int N_nodes, int opt_self, int opt_dir){
	double S,p,mu;
	int i,j,t;
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
                if(t>=0)
                {
					t = WG->node[i].w_out[t];
				}else{
					t = 0;
				}
                mu = y[i]*yy[j]/(1+y[i]*yy[j]);
                if(mu>0)
				{
					p = mu;
					//printf("Mu:%f p:%f t:%d",mu,p,t);fflush(stdout);
				}else{
					p = 1./0.; // infinity! (not compatible)
					if(t==0) p = 0;
				}
				if(p>0) S+= log(p);
			}
        }
	}
	if(opt_dir<=0) S = S/2.;
	return -S;    
}


double w_graph_loglikelyhood_ZIP_xy(W_GRAPH* WG, double** x, int N_nodes, double gamma, int opt_self, int opt_dir){
	double S,p;
	int i,j,t,dest;
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
				t=find_value_int(WG->node[i].out, j, WG->node[i].kout);
                if(t>-1) // connection   
                {
					dest = WG->node[i].out[t];
					t = WG->node[i].w_out[dest];
					p = gamma*exp(y[i]*yy[dest])/(gamma*(exp(y[i]*yy[dest])-1)+1) * gsl_ran_poisson_pdf (t, y[i]*yy[dest]);//
				}else{   
                    p = 1 - gamma*(exp(y[i]*yy[j])-1)/(gamma*(exp(y[i]*yy[j])-1)+1);
				}
            }		
		    if(p>0) S+= log(p);
        }
	}
	if(opt_dir<0) S =S/2.;
	return -S;    
}



double w_graph_loglikelyhood_ZIP2_xy(W_GRAPH* WG, double** x,int N_nodes,  int opt_self, int opt_dir){
	double S,p;
	int i,j,t,dest;
    double *y;
    double *yy;
    double * z;
    double * zz;
	S = 0;
    if(opt_dir>0)
    {
        y = x[0];
        yy = x[1];
        z = x[2];
        zz = x[3];
    }else{
        y = x[0];
        yy = x[0];
        z = x[1];
        zz = x[1];
    }
	for(i=0;i<N_nodes;i++) // all edges
	{
 		for(j=0;j<N_nodes;j++)
		{
            if((opt_self>0)||(i!=j))
            {
				t=find_value_int(WG->node[i].out, j, WG->node[i].kout);
                if(t>-1) // connection   
                {
					dest = WG->node[i].out[t];
					t = WG->node[i].w_out[dest];
					p = z[i]*zz[dest]*exp(y[i]*yy[dest])/(z[i]*zz[dest]*(exp(y[i]*yy[dest])-1)+1) * gsl_ran_poisson_pdf (t, y[i]*yy[dest]);//
				}else{   
                    p = 1 - z[i]*zz[j]*(exp(y[i]*yy[j])-1)/(z[i]*zz[j]*(exp(y[i]*yy[j])-1)+1);
				}
            }		
		    if(p>0) S+= log(p);
        }
	}
	if(opt_dir<0) S =S/2.;
	return -S;    
}


<<<<<<< HEAD
double w_graph_loglikelyhood_ZIB2_xy(W_GRAPH* WG, double** x, int N_nodes, int layers, int opt_self, int opt_dir){
	double S,p;
	int i,j,t,dest;
    double *y;
    double *yy;
    double * z;
    double * zz;
	S = 0;
    if(opt_dir>0)
    {
        y = x[0];
        yy = x[1];
        z = x[2];
        zz = x[3];
    }else{
        y = x[0];
        yy = x[0];
        z = x[1];
        zz = x[1];
    }
	for(i=0;i<N_nodes;i++) // all edges
	{
 		for(j=0;j<N_nodes;j++)
		{
            if((opt_self>0)||(i!=j))
            {
				t=find_value_int(WG->node[i].out, j, WG->node[i].kout);
                if(t>-1) // connection   
                {
					dest = WG->node[i].out[t];
					t = WG->node[i].w_out[dest];
					p = z[i]*zz[dest]*pow(1+y[i]*yy[dest],layers)/(z[i]*zz[dest]*(pow(1+y[i]*yy[dest],layers)-1)+1) * gsl_ran_binomial_pdf (t, y[i]*yy[dest]/(1+y[i]*yy[dest]),layers);//
				}else{   
                    p = 1 - z[i]*zz[j]*(pow(1+y[i]*yy[j],layers)-1)/(z[i]*zz[j]*(pow(1+y[i]*yy[j],layers)-1)+1);
				}
            }		
		    if(p>0) S+= log(p);
        }
	}
	if(opt_dir<0) S =S/2.;
	return -S;  
=======
void w_graph_print_entropy(double* seq,int len,char* output){
	int bins;
	double mins,maxs,eps;
	bins = (int)len/20.;
	if(bins<5)
	{
		bins=5;
	}else if(bins>50){
		bins=50;
	}
	mins = min_value_double(seq,len);
	maxs = max_value_double(seq,len);
	eps = (maxs-mins)/1000.;
    mins = mins - mins*eps;
    maxs = maxs+maxs*eps;
    if(mins<0) mins = 0;
    if(maxs<0) maxs = 0;
    if(mins>=maxs) maxs = mins+1;
    gsl_histogram* h1=histogram_double(seq,mins,maxs,bins,len);
	print_acc(output,h1,h1);
	return;
}

double w_graph_loglikelyhood_poisson(W_GRAPH* WG,int N_nodes,double** wij){
	double L,p,mu;
	int i,j,t;
	L = 0;
	for(i=0;i<N_nodes;i++) // all edges
	{
		for(j=0;j<N_nodes;j++)
		{
            if((opt_self>0)||(i!=j))
            {
                t=find_value_int(WG->node[i].out, j, WG->node[i].kout);
                if(t>=0)
                {
					t = WG->node[i].w_out[t];
				}else{
					t = 0
				}
                mu = wij[i][j];
				if(mu>0)
				{
					p = gsl_ran_poisson_pdf (t, mu);
					//printf("Mu:%f p:%f t:%d",mu,p,t);fflush(stdout);
					L+= log(p);
				}
			}
		}
	}
	return L;    
>>>>>>> e4d232dadd44cc2975bd6a834f9b613f25a7c94b
}








/****************************************************************************
 * Indices *
 ****************************************************************************/
double w_graph_compute_sorensen(W_GRAPH* WG, W_GRAPH* WGoriginal, int N_nodes){
	int i,j,ncc;
	double soren=0;
	int norm = w_graph_total_weight(WG,N_nodes);
	norm += w_graph_total_weight(WGoriginal,N_nodes);
	int* ts = cast_vec_int(N_nodes); // dummy vector
	for(i=0;i<N_nodes;i++)
	{
		ncc = 0;
		for(j=0;j<WGoriginal->node[i].kout;j++)
		{
			ts[WGoriginal->node[i].out[j]] =  WGoriginal->node[i].w_out[j];
		}
		for(j=0;j<WG->node[i].kout;j++)
		{
			ncc+= mineq_int(WG->node[i].w_out[j],ts[WG->node[i].out[j]]);
		}
		for(j=0;j<N_nodes;j++)
		{
			ts[j]=0;
		}
		soren +=(double)ncc;
	}
	return 2*soren/(double)norm;
}
double w_graph_compute_sorensen_av(W_GRAPH* WGoriginal, double** pij, int N_nodes, double T){
	int i,j,ncc;
	double soren=0;
	double norm = sum_matrix_double(pij, N_nodes, N_nodes);
	double norm2 = 0;
	double t;
	scale_matrix(pij, N_nodes, N_nodes, T/norm); // scale and normalize	
	for(i=0;i<N_nodes;i++)
	{
		ncc = 0;
		for(j=0;j<WGoriginal->node[i].kout;j++)
		{
			t = fmin((double)WGoriginal->node[i].w_out[j],pij[i][WGoriginal->node[i].out[j]]);
			norm2 += pij[i][WGoriginal->node[i].out[j]];
			norm2 += (double)WGoriginal->node[i].w_out[j];
			soren +=t;
		}
	}
	return 2*soren/norm2;
}
/****************************************************************************
 * aLL STATS *
 ****************************************************************************/
void w_graph_node_stats_list(W_GRAPH* WG, int N_nodes, int run, double av_k, int opt_dir, int opt_clust, int self_opt){
    //p(k),p(s),p(w),w(sin sout), s_nn, k_nn, s(k)
    int **k=w_graph_compute_k(WG, N_nodes);
    int **s=w_graph_compute_s(WG, N_nodes);
    //int E;
    //long int T=sum_vec_int(s[0],N_nodes);
    //int *wkk2=w_graph_compute_w_ss(WG->node, N_nodes, 0);
    double ** ss_n=w_graph_compute_s_nn(WG, N_nodes,1, opt_dir);// 4 rows
    double ** kk_n=w_graph_compute_s_nn(WG, N_nodes,-1, opt_dir);// 4 rows
    double ** kkw_n=w_graph_compute_s_nn(WG, N_nodes,0, opt_dir);// 4 rows
    double **y2 = w_graph_compute_Y2(WG, N_nodes, opt_dir); // 2 rows
    char cadena[100];
    int i; 
    // Get analytical predictions //
    // k(s) //
	if(opt_dir>0)
	{
		sprintf(cadena,"N%davs%8.5fnode_list.list",N_nodes,av_k);
	}else{
		sprintf(cadena,"N%davs%8.5f_undir_node_list.list",N_nodes,av_k);
	}
    
    FILE* fil=open_file("w", cadena);
    if(opt_dir==1)
    {
    	double ** k_anal=w_graph_compute_k_analitic(WG, N_nodes, self_opt);
		// Note for indices on s_nn:
		// [0] -> \sum t_{ij} s_out_j / s_out_i
		// [1] -> \sum t_{ij} s_in_j / s_out_i
		// [2] -> \sum t_{ji} s_out_j / s_in_i
		// [3] -> \sum t_{ji} s_in_j / s_in_i
    	fprintf(fil,"# Node_num\tk\tk_anal\ts\tY2\tk_nn\tk^w_nn\ts^w_nn (out) then (in) # \n");
    	for(i=0;i<N_nodes;i++)
    	{
        	fprintf(fil,"%d %d %.3f %d %f %f %f %f %d %.3f %d %f %f %f %f\n",i,k[0][i],k_anal[0][i],s[0][i],
			y2[0][i],kk_n[1][i],kkw_n[1][i],ss_n[1][i],k[1][i],k_anal[1][i],s[1][i],y2[1][i],kk_n[2][i],kkw_n[2][i],ss_n[2][i]);
    	}
	free_mat_double(k_anal,2);
    }else{
	double ** k_anal=w_graph_compute_k_analitic(WG, N_nodes, self_opt);
	if(opt_clust==1)
	{
	    double ** c= w_graph_compute_clust(WG, N_nodes); // 2 rows
	    fprintf(fil,"# Node_num\tk\tk_anal\ts\tY2\tk_nn\tk^w_nn\ts^w_nn\tc\tc^w# \n");
	    for(i=0;i<N_nodes;i++)
	    {
        	fprintf(fil,"%d %d %.3f %d %f %f %f %f %f %f\n",i,k[0][i],k_anal[0][i],s[0][i],y2[0][i],kk_n[0][i],
			kkw_n[0][i],ss_n[0][i],c[0][i],c[1][i]);
	    }
	    free_mat_double(c,2);
	}else{
	    fprintf(fil,"# Node_num\tk\tk_anal\ts\tY2\tk_nn\tk^w_nn\ts^w_nn\n");
	    for(i=0;i<N_nodes;i++)
	    {
        	fprintf(fil,"%d %d %.3f %d %f %f %f %f\n",i,k[0][i],k_anal[0][i],s[0][i],y2[0][i],kk_n[0][i],
			kkw_n[0][i],ss_n[0][i]);
	    }
	free_mat_double(k_anal,2);
	}
    }
    fclose(fil);
    free_mat_int(s,2);
    free_mat_int(k,2);
    free_mat_double(ss_n,4);
    free_mat_double(kk_n,4);
    free_mat_double(kkw_n,4);
    free_mat_double(y2,2);
    return;
}


void w_graph_all_stats(W_GRAPH* WG, int N_nodes, int run, double bin_exp, double av_k, int opt_dir, int self_opt, int w_anal){
    //w(sin sout), w(kin,kout), w
    int **k=w_graph_compute_k(WG, N_nodes);
    int **s=w_graph_compute_s(WG, N_nodes);
    int E,L;


    int *w_zeros=w_graph_compute_w(WG, N_nodes, &L, 1);
    int *p_zeros=w_graph_compute_p(WG, N_nodes, &L);
    int *w=w_graph_compute_w(WG, N_nodes, &E, -1);

	double *wss=w_graph_compute_wp_ss(WG, N_nodes, 1);
    double *wkk=w_graph_compute_wp_ss(WG, N_nodes, -1);
    double *wss_zeros=w_graph_compute_w_ss(WG, N_nodes, 1);

    char cadena[100];
    
    gsl_histogram* h1;
    
    double* sout;
    double* xranges;
    int xbins;
    double** yy;
    
    /// w histogram /////
    sout=vec_int_to_double(w,E);
    int q=max_value_int(w,E);
    h1=histogram_double(sout,0,q,q,E);
    free(sout);
    //sprintf(cadena,"run_%dN%d_w.hist",run,N_nodes);
	if(opt_dir>0)
	{
	    sprintf(cadena,"N%davs%8.5f_w.hist",N_nodes,av_k);
	}else{
	    sprintf(cadena,"N%davs%8.5f_undir_w.hist",N_nodes,av_k);
	}

    print_acc(cadena, h1, h1);
    gsl_histogram_free(h1);

    /// analitical w /////
    if(w_anal>0)
    {
        double** pp;
        int lenn;
        if(opt_dir>0)
        {
            pp = w_graph_compute_p_w_analitic_from_s_directed(10*q,1.5,s,N_nodes, self_opt, &lenn);
			sprintf(cadena,"N%davs%8.5f_w_anal.hist",N_nodes,av_k);
        }else{
            pp = w_graph_compute_p_w_analitic_from_s_undirected(10*q,1.5,s[0],N_nodes, self_opt, &lenn);
			sprintf(cadena,"N%davs%8.5f_undir_w_anal.hist",N_nodes,av_k);
        }
        FILE* fil=open_file("w", cadena);
        fprintf(fil,"# t p(t) # \n");
        int i;
        for(i=0;i<lenn;i++)
        {
            fprintf(fil,"%.8f %.8f\n",pp[0][i],pp[1][i]);
        }
        fclose(fil);
        //free(pp[0]);
        //free(pp[1]);
        //free(pp);
    }

    /// exsiting weight as func of ss /////
    sout=vec_int_to_double(w,E);
    xranges=log_bins_double(0, max_value_double(wss,E) , bin_exp, &xbins);
    yy=y_of_x(wss, sout, xranges,  E,  xbins);
	if(opt_dir>0)
	{
	    sprintf(cadena,"N%davs%8.5f_wp_s_oi.hist",N_nodes,av_k);
	}else{
	    sprintf(cadena,"N%davs%8.5f_undir_wp_s_oi.hist",N_nodes,av_k);
	}
    print_hist2d_mean(cadena, yy[1], yy[2], yy[0], xbins);
    free(sout);
    free(wss);
    free(xranges);
    free_mat_double(yy,4);

    /// existing weight as func of kk /////
    sout=vec_int_to_double(w,E);
    xranges=log_bins_double(0, max_value_double(wkk,E) , bin_exp, &xbins);
    yy=y_of_x(wkk, sout, xranges,  E,  xbins);
	if(opt_dir>0)
	{
	    sprintf(cadena,"N%davs%8.5f_wp_k_oi.hist",N_nodes,av_k);
	}else{
	    sprintf(cadena,"N%davs%8.5f_undir_wp_k_oi.hist",N_nodes,av_k);
	}

    print_hist2d_mean(cadena, yy[1], yy[2], yy[0], xbins);
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
	    sprintf(cadena,"N%davs%8.5f_w_s_oi.hist",N_nodes,av_k);
	}else{
	    sprintf(cadena,"N%davs%8.5f_undir_w_s_oi.hist",N_nodes,av_k);
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
	    sprintf(cadena,"N%davs%8.5f_p_s_oi.hist",N_nodes,av_k);
	}else{
	    sprintf(cadena,"N%davs%8.5f_undir_p_s_oi.hist",N_nodes,av_k);
	}
    print_hist2d_mean(cadena, yy[1], yy[2], yy[0], xbins);
    free(sout);
    free(wss_zeros);
    free(xranges);
    free_mat_double(yy,4);    

//// free all
    free_mat_int(s,2);
    free_mat_int(k,2);
    free(w);
    free(w_zeros);
    free(p_zeros);
    return;
}


/****************************************************************************
 * Ensembles averaging stats *
 ****************************************************************************/

// Nodes
void w_graph_node_stats_ensemble(W_GRAPH* WG, int N_nodes, double** container, double ** container2, int** node_nonzero ,double* T_container, int opt_dir, int opt_clust ){
    //p(k),p(s),p(w),w(sin sout), s_nn, k_nn, s(k)
    int **k=w_graph_compute_k(WG, N_nodes);
    int **s=w_graph_compute_s(WG, N_nodes);
    //int E;
    int T=sum_vec_int(s[0],N_nodes);
    //int *wkk2=w_graph_compute_w_ss(WG->node, N_nodes, 0);
    double ** ss_n=w_graph_compute_s_nn(WG, N_nodes,1, opt_dir);// 4 rows
    double ** kk_n=w_graph_compute_s_nn(WG, N_nodes,-1, opt_dir);// 4 rows
    double ** kkw_n=w_graph_compute_s_nn(WG, N_nodes,0, opt_dir);// 4 rows
    double **y2 = w_graph_compute_Y2(WG, N_nodes, opt_dir); // 2 rows
    //double **xy = w_graph_compute_xy(WG->node, N_nodes); // 2 rows
    //char cadena[100];
		// Note for indices on s_nn:
		// [0] -> \sum t_{ij} s_out_j / s_out_i
		// [1] -> \sum t_{ij} s_in_j / s_out_i
		// [2] -> \sum t_{ji} s_out_j / s_in_i
		// [3] -> \sum t_{ji} s_in_j / s_in_i

    int i;
    if(opt_dir==1)
    { 
    	for(i=0;i<N_nodes;i++)
    	{
	    if(s[0][i] > 0)
	    {
		//assert(s[0][i]>=k[0][i]);
		node_nonzero[i][0]+=1;
		container[i][0]+=(double)k[0][i];
		container2[i][0]+=(double)k[0][i]*k[0][i];
		container[i][1]+=(double)s[0][i];
		container2[i][1]+=(double)s[0][i]*s[0][i];
		container[i][2]+=(double)y2[0][i];
		container2[i][2]+=(double)y2[0][i]*y2[0][i];
		container[i][3]+=(double)kk_n[1][i];
		container2[i][3]+=(double)kk_n[1][i]*kk_n[1][i];
		container[i][4]+=(double)kkw_n[1][i];
		container2[i][4]+=(double)kkw_n[1][i]*kkw_n[1][i];
		container[i][5]+=(double)ss_n[1][i];
		container2[i][5]+=(double)ss_n[1][i]*ss_n[1][i];

	    }
	    if(s[1][i] > 0)
	    {
		//assert(s[1][i]>=k[1][i]);
		node_nonzero[i][1]+=1;
		container[i][6]+=(double)k[1][i];
		container2[i][6]+=(double)k[1][i]*k[1][i];
		container[i][7]+=(double)s[1][i];
		container2[i][7]+=(double)s[1][i]*s[1][i];
		container[i][8]+=(double)y2[1][i];
		container2[i][8]+=(double)y2[1][i]*y2[1][i];
		container[i][9]+=(double)kk_n[2][i];
		container2[i][9]+=(double)kk_n[2][i]*kk_n[2][i];
		container[i][10]+=(double)kkw_n[2][i];
		container2[i][10]+=(double)kkw_n[2][i]*kkw_n[2][i];
		container[i][11]+=(double)ss_n[2][i];
		container2[i][11]+=(double)ss_n[2][i]*ss_n[2][i];
	    }
    	}
    }else{
		if(opt_clust==1)
		{
	    	double** c=w_graph_compute_clust(WG, N_nodes);
	    	for(i=0;i<N_nodes;i++)
	    	{
	    		if(s[0][i] > 0)
	    		{		
					container[i][6]+=c[0][i];
					container2[i][6]+=c[0][i]*c[0][i];
					container[i][7]+=c[1][i];
					container2[i][7]+=c[1][i]*c[1][i];
	    		}
	    	}
	    	free_mat_double(c,2);
		}
		for(i=0;i<N_nodes;i++)
		{
			if(s[0][i] > 0)
	    	{
				//assert(s[0][i]>=k[0][i]);		
				node_nonzero[i][0]+=1;
				container[i][0]+=(double)k[0][i];
				container2[i][0]+=(double)k[0][i]*k[0][i];
				container[i][1]+=(double)s[0][i];
				container2[i][1]+=(double)s[0][i]*s[0][i];
				container[i][2]+=(double)y2[0][i];
				container2[i][2]+=(double)y2[0][i]*y2[0][i];
	    		container[i][3]+=(double)kk_n[0][i];
				container2[i][3]+=(double)kk_n[0][i]*kk_n[0][i];
				container[i][4]+=(double)kkw_n[0][i];
				container2[i][4]+=(double)kkw_n[0][i]*kkw_n[0][i];
				container[i][5]+=(double)ss_n[0][i];
				container2[i][5]+=(double)ss_n[0][i]*ss_n[0][i];
/*
		container[i][4]+=(double)xy[0][i]; // x
		container2[i][4]+=(double)xy[0][i]*(double)xy[0][i];
		container[i][5]+=(double)xy[1][i]; // y
		container2[i][5]+=(double)xy[1][i]*(double)xy[1][i];
		container[i][6]+=(double)xy[2][i]; //xy
		container2[i][6]+=(double)xy[2][i]*(double)xy[2][i];
*/
	    }
	}
	}
    T_container[0]+=T;
    T_container[1]+=T*T;
    free_mat_int(k,2);
    free_mat_int(s,2);
    free_mat_double(ss_n,4);
    free_mat_double(kk_n,4);
    free_mat_double(kkw_n,4);
    free_mat_double(y2,2);
    return;
}


void w_graph_node_stats_ensemble_print(int reps, int N_nodes, double* Tcont, double** cont, double ** cont2, int** node_nonzero, double av_k, double bin_exp, int len_acc, int opt_dir){
    char cadena[100];
    int i,j;
    //printf("i am printing\n");fflush(stdout);
    scale_vec_double(Tcont,1./(double)reps,2);
    //average_matrix(cont, N_nodes, len_acc, reps);
    //average_matrix(cont2, N_nodes, len_acc, reps);
    ///// Only over existing but not for degrees or strengths!!!! ////
	//printf("alohaaa %d\n\n",len_acc/2);fflush(stdout);
    for(i=0;i<N_nodes;i++)
    {
	if (opt_dir==1)
	{
	    for(j=0;j<2;j++) // s and k
	    {
		cont[i][j]=cont[i][j]/(double)reps;
		cont2[i][j]=cont2[i][j]/(double)reps;
	    }	    
	    for(j=len_acc/2;j<len_acc/2+2;j++) // s and k
	    {
		cont[i][j]=cont[i][j]/(double)reps;
		cont2[i][j]=cont2[i][j]/(double)reps;
	    }
	    
	    for(j=2;j<len_acc/2;j++)
	    {
			if(node_nonzero[i][0]>0)
			{
				cont[i][j]=cont[i][j]/(double)node_nonzero[i][0];
				cont2[i][j]=cont2[i][j]/(double)node_nonzero[i][0];
			}else{
				cont[i][j]=0;
				cont2[i][j]=0;
			}
	    }	    
	    for(j=len_acc/2+2;j<len_acc;j++)
	    {
			if(node_nonzero[i][1]>0)
			{
				cont[i][j]=cont[i][j]/(double)node_nonzero[i][1];
				cont2[i][j]=cont2[i][j]/(double)node_nonzero[i][1];
			}else{
				cont[i][j]=0;
				cont2[i][j]=0;
			}
	    }	
	}else{
	    for(j=0;j<2;j++)
	    {
		cont[i][j]=cont[i][j]/(double)reps;
		cont2[i][j]=cont2[i][j]/(double)reps;
	    }	    
	    for(j=2;j<len_acc;j++)
	    {
			if(node_nonzero[i][0]>0)
			{
				cont[i][j]=cont[i][j]/(double)node_nonzero[i][0];
				cont2[i][j]=cont2[i][j]/(double)node_nonzero[i][0];
			}else{
				cont[i][j]=0;
				cont2[i][j]=0;
			}
	    }
	}
    }
    if(opt_dir==1)
    {
		sprintf(cadena,"N%davs%.5f_ens_r%dnode_list.list",N_nodes,av_k,reps);
	}else{
		sprintf(cadena,"N%davs%.5f_undir_ens_r%d_node_list.list",N_nodes,av_k,reps);		
	}
    FILE* fil=open_file("w", cadena);
    fprintf(fil,"# <T>=%f+-%f # \n",Tcont[0] ,sqrt(Tcont[1]-Tcont[0]*Tcont[0]));
    if(opt_dir==1)
    {
    	fprintf(fil,"# Node_num\tk\tsset\tY2\tk_nn\tk^w_nn\ts^w_nn (out) then (in) # \n");
    }else{
	fprintf(fil,"# Node_num\tk\tsset\tY2\tk_nn\tk^w_nn\ts^w_nn (optionally \tc\tc^w) # \n");
    }
    //int i,j;
    //double** node_atts=cast_mat_double(len_acc,N_nodes);
    for(i=0;i<N_nodes;i++)
    {
        fprintf(fil,"%d",i);
		//assert(cont[i][1]<=cont[i][2]);
		//assert(cont[i][7]<=cont[i][9]);
    	for(j=0;j<len_acc;j++)        
    	{
	    fprintf(fil," %f %f",cont[i][j],sqrt(cont2[i][j]-cont[i][j]*cont[i][j]));
	    //node_atts[j][i] = cont[i][j];
	}
	    fprintf(fil,"\n");
    }
    fclose(fil);
    /*
    gsl_histogram * h1;
    // str in
    h1=histogram_double_log(WG->node_atts[2],0,max_value_double(WG->node_atts[2],N_nodes),bin_exp,N_nodes);
    sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%dnode_sIN.hist",N_nodes,av_k,expo-1,reps);
    print_acc(cadena, h1, h1);
    gsl_histogram_free(h1);
    // degrees in
    h1=histogram_double_log(WG->node_atts[0],0,max_value_double(WG->node_atts[0],N_nodes),bin_exp,N_nodes);
    sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%dnode_kIN.hist",N_nodes,av_k,expo-1,reps);
    print_acc(cadena, h1, h1);
    gsl_histogram_free(h1);

    if(opt_dir==1)
    {
    // degrees out
    h1=histogram_double_log(WG->node_atts[7],0,max_value_double(WG->node_atts[7],N_nodes),bin_exp,N_nodes);	
    sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%dnode_kOUT.hist",N_nodes,av_k,expo-1,reps);
    print_acc(cadena, h1, h1);
    gsl_histogram_free(h1);
    // str out
    h1=histogram_double_log(WG->node_atts[9],0,max_value_double(WG->node_atts[9],N_nodes),bin_exp,N_nodes);
    sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%dnode_sOUT.hist",N_nodes,av_k,expo-1,reps);
    print_acc(cadena, h1, h1);
    gsl_histogram_free(h1);
    }
    */
    return;
}
// P(k),P(s),P(w)

gsl_histogram ** w_graph_all_stats_ensemble_allocate(int dir, int s_min, int s_max, int k_min, int k_max, int w_max){
	int len_acc=2;
	int bins;
	if(w_max>1e7)
	{
		bins = 10000000;
	}else{
		bins=w_max;
	}	/*
	if(dir==1) 
	{
		len_acc=10;
	}else{
		len_acc=6;
	}
	*/
	gsl_histogram ** acc = (gsl_histogram**)malloc(sizeof(gsl_histogram*)*len_acc);
	acc[0] = set_acc_double(0,w_max,bins);
	acc[1] = set_acc_double(0,w_max,bins);
	/*
	if(dir==1)
	{
	
		acc[0] = set_acc_int(s_min,s_max);
		acc[1] = set_acc_int(s_min,s_max);
		acc[2] = set_acc_int(s_min,s_max);
		acc[3] = set_acc_int(s_min,s_max);

		acc[4] = set_acc_int(k_min,k_max);
		acc[5] = set_acc_int(k_min,k_max);
		acc[6] = set_acc_int(k_min,k_max);
		acc[7] = set_acc_int(k_min,k_max);

		acc[8] = set_acc_int(0,w_max);
		acc[9] = set_acc_int(0,w_max);

	}else{
		acc[0] = set_acc_int(s_min,s_max);
		acc[1] = set_acc_int(s_min,s_max);
		acc[2] = set_acc_int(k_min,k_max);
		acc[3] = set_acc_int(k_min,k_max);
		acc[4] = set_acc_int(0,w_max);
		acc[5] = set_acc_int(0,w_max);
	}
	*/
	return acc;
}

void w_graph_all_stats_ensemble_update(gsl_histogram** acc, W_GRAPH* WG, int N_nodes, int dir){
    	int E;
	//int **k=w_graph_compute_k(WG->node, N_nodes);
	//int **s=w_graph_compute_s(WG->node, N_nodes);
	int *w=w_graph_compute_w(WG, N_nodes, &E, -1);
	//long int T=sum_vec_int(s[0],N_nodes);
	/*
	if(dir==1)
	{
		update_int_acc(s[0], N_nodes, acc[0], acc[1], -1);
		update_int_acc(s[1], N_nodes, acc[2], acc[3], -1);
		update_int_acc(k[0], N_nodes, acc[4], acc[5], -1);
		update_int_acc(k[1], N_nodes, acc[6], acc[7], -1);
		update_int_acc(w, N_nodes, acc[8], acc[9], -1);
	}else{
		update_int_acc(s[0], N_nodes, acc[0], acc[1], -1);
		update_int_acc(k[0], N_nodes, acc[2], acc[3], -1);
		update_int_acc(w, N_nodes, acc[4], acc[5], -1);
	}
	*/
	update_int_acc(w, E, acc[0], acc[1], 1);
	free(w);
	return;
}
void w_graph_all_stats_ensemble_print(gsl_histogram** acc, int len, int reps, int N_nodes, double av_k, int opt_dir){
	char	cadena[100];
	//Normalize and std
	acc_normalize(acc,len, 1./(double)reps);
	acc_compute_std(acc,len, 1.);
	// Print all
	/*
	if(dir==1)
	{
		sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%d_sIN.hist",N_nodes,av_k,expo-1,reps);
		print_acc(cadena, acc[0], acc[1]);
		sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%d_sOUT.hist",N_nodes,av_k,expo-1,reps);
		print_acc(cadena, acc[2], acc[3]);
		sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%d_kIN.hist",N_nodes,av_k,expo-1,reps);
		print_acc(cadena, acc[4], acc[5]);
		sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%d_kOUT.hist",N_nodes,av_k,expo-1,reps);
		print_acc(cadena, acc[6], acc[7]);
		sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%d_w.hist",N_nodes,av_k,expo-1,reps);
		print_acc(cadena, acc[8], acc[9]);
	}else{
		sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%d_sOUT.hist",N_nodes,av_k,expo-1,reps);
		print_acc(cadena, acc[0], acc[1]);
		sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%d_kOUT.hist",N_nodes,av_k,expo-1,reps);
		print_acc(cadena, acc[2], acc[3]);
		sprintf(cadena,"N%davs%8.5fexpo%.2f_ens_r%d_w.hist",N_nodes,av_k,expo-1,reps);
		print_acc(cadena, acc[4], acc[5]);
	}
	*/
    if(opt_dir==1)
    {
		sprintf(cadena,"N%davs%8.5f_ens_r%d_w.hist",N_nodes,av_k,reps);
	}else{
		sprintf(cadena,"N%davs%8.5f_undir_ens_r%d_w.hist",N_nodes,av_k,reps);
	}
	print_acc(cadena, acc[0], acc[1]);
	//Free all
	acc_free_all(acc, len);
	return;
}

/****************************************************************************
 * Other funcs *
 ****************************************************************************/

double w_graph_compute_rho(double E_av, double T, int opt_indist){
    double rho;
	if(opt_indist<=0) // ME
	{
		double x = - (double)T/E_av * exp(-(double)T/E_av);
		//printf("x:%f tplus:%f\n",gsl_sf_lambert_W0 (x),(double)T/E_av);fflush(stdout);
		rho = gsl_sf_lambert_W0 (x) + (double)T/E_av;
	}else{ //W,AW
		rho = 1. - E_av / (double)T; // 1-p = E/T
	}
	return rho;
}


/****************************************************************************
 * Transformations *
 ****************************************************************************/
double** w_graph_to_adj_matrix(W_GRAPH* WG, int N_nodes){
	double** wij = cast_mat_double(N_nodes,N_nodes);
	int i,j;
	for(i=0;i<N_nodes;i++)
	{
		for(j=0;j<WG->node[i].kout;j++)
			wij[i][WG->node[i].out[j]] = (double)WG->node[i].w_out[j];
	}
	return wij;
}



/****************************************************************************
 * Graph Filtering *
 ****************************************************************************/
W_GRAPH* w_graph_filter_xij(W_GRAPH* WG, double* x, double* y, int N_nodes, double gamma, int mode, int M){
    // filters graph if t in in confidence bounds according to C.I gamma //
	int i,j;
    int id_node;
    double mu;
    int t;
    int tot = 0;
    int etot = 0;
    int* l;
    int theta;
    assert(gamma<=1);
    assert(gamma>=0);
    W_GRAPH* WGf = w_graph_alloc(N_nodes);
    if((mode<0)||(mode>2))
    {
        printf("Invalid mode. Must be 0 (ME), 1(B) or 2(W). Aborting...\n");
        abort();
    }
    double T = (double)w_graph_total_weight(WG, N_nodes);
    double E = (double)w_graph_total_edges(WG, N_nodes);
	for(i=0;i<N_nodes;i++)
	{
		for(j=0;j<WG->node[i].kout;j++)
        {
            id_node = WG->node[i].out[j];
            t = WG->node[i].w_out[j];
            mu = x[i]*y[id_node];
            // if theta==1 keep, else, not keep
            if(mu==0)
            {
                theta = 1; // surely not in interval, keep
            }else{
                l = find_tmintmax_xy(mu,gamma,mode,M);                
                if((l[0]<0)&&(l[1]<0))
                {
                    theta = 0; // surely in interval, do not keep
                }else{
                    if((t>=l[0])&&(t<=l[1])) // in interval (at most with gamma% chance)
                    {
                        // |_1 gamma |_2 --> prob outside <1-gamma
                        theta=0; // do not keep
                    }else{
                        theta=1; // keep
                    }
                }
            }
            if(theta==1) // if not in interval, keep
            {
                //if(t>50)printf("t:%d mu:%f tmin:%d tmax:%d\n",t,mu,l[0],l[1]);
                w_graph_add_multi_link(WGf, N_nodes, i, id_node, t);
                tot+=t;
                etot+=1;
            }
        }
    }
    printf("\tTotal number of units after filtering: events: %d (f:%.5f) | edges:%d (f:%.5f) \n",tot,(double)tot/T,etot,(double)etot/E);
    return WGf;
}



int* find_tmintmax_xy(double xy, double gamma, int mode, int M)
{
    assert(gamma<=1);
    assert(gamma>=0);
    // finds T such that p(x<=T) = gamma;
    int* tt = cast_vec_int(2);
    double mu,p,P0;
    int l1,l2,k,aux;
    int t;
    //int c;
    // mode selection //
    if (mode==0)
    {
        mu = xy;
        t = (int)round(mu);
        tt[0] = t;
        tt[1] = t;        
        if(mu<=t) // forward
        {
            k=1;
            l1 = t;
            l2 = t-1;
        }else{
            k=-1;
            l1 = t+1;
            l2 = t;
        }
        P0 = 0;
        //c = 0;
        //if(mu>100)printf("c :%d t:%d xy:%f k:%d aux:%d l1:%d l2:%d P0:%f Gamma:%f\n",c,t,xy,k,aux,l1,l2,P0,gamma);
        while(P0<gamma)
        {
            if(k>0) // forward
            {
                aux = l1;
                l1= l1 + 1;
                k=-k;
            }else{ // backwards
                aux = l2;
                l2= l2 - 1;
                if(l2<0)
                {
                    k=1; // if non negative, keep positive
                }else{
                    k=-k;
                }
            }
            P0 += gsl_ran_poisson_pdf(aux,mu);                
            //c++;
            if(fabs(aux-mu)>10*sqrt(mu)) // 10 std away!
            {
                l1 = -1;
                l2 = -1;
                break; //break if too large
            }
            //if(mu>100)printf("c :%d t:%d aux:%d xy:%f k:%d P0:%f Gamma:%f\n",c,t,aux,xy,k,P0,gamma);
        }
    }else if(mode==1){
        mu = M*(xy/(1+xy));        
        p = xy/(1+xy);
        t = (int)round(mu);
        tt[0] = t;
        tt[1] = t;        
        if(mu<=t) // forward
        {
            k=1;
            l1 = t;
            l2 = t-1;
        }else{
            k=-1;
            l1 = t+1;
            l2 = t;
        }
        P0 = 0;
        while(P0<gamma)
        {
            if(k>0) // forward
            {
                aux = l1;
                l1+=1;
                k=-k;
            }else{ // backwards
                aux = l2;
                l2-=1;
                if(l2<0)
                {
                    k=1; // if non negative, keep positive
                }else{
                    k=-k;
                }
            }
            P0 += gsl_ran_binomial_pdf(aux,p,M);
            if(fabs(aux-mu)>10*sqrt(mu/(1+xy))) // 10 std away!
            {
                l1 = -1;
                l2 = -1;
                break; //break if too large
            }
        }
    }else{
        mu = M*xy/(1.-xy);        
        p = (1.-xy);
        t = (int)round(mu-xy/p);
        tt[0] = t;
        tt[1] = t;        
        if(mu<=t) // forward
        {
            k=1;
            l1 = t;
            l2 = t-1;
        }else{
            k=-1;
            l1 = t+1;
            l2 = t;
        }
        P0 = 0;
        while(P0<gamma)
        {
            if(k>0) // forward
            {
                aux = l1;
                l1+=1;
                k=-k;
            }else{ // backwards
                aux = l2;
                l2-=1;
                if(l2<0)
                {
                    k=1; // if non negative, keep positive
                }else{
                    k=-k;
                }
            }
            P0 += gsl_ran_negative_binomial_pdf(aux,p,M); 
            // --> think about this! p(k) = {\Gamma(n + k) \over \Gamma(k+1) \Gamma(n) } p^n (1-p)^k
            if(fabs(aux-mu)>10*sqrt(mu/p)) // 10 std away!
            {
                l1 = -1;
                l2 = -1;
                break; //break if too large
            }
        }
    }
	//printf("mu:%f t:%d xy:%f k:%d aux:%d l1:%d l2:%d P0:%f Gamma:%f\n",mu,t,xy,k,aux,l1,l2,P0,gamma);
    //printf("FInal values: %d %d %f\n",tt[0]-l2,tt[1]+l1,mu);
    tt[0] = l2;
    if(l2<0) tt[0]=0;
    tt[1] = l1;
    return tt;
}
