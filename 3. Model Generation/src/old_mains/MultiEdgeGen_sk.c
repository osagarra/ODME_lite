/********************************************************************************
 * 
 *			GENERATE Maximum entropy ensembles for networks
 *
 *	This program generates networks (and statistics) with pre-defined node strenghts and degrees using a
 *  macro-canonical approach.
 *
 *   More details in
 * [1] O. Sagarra, C. J. Pérez Vicente, and A. Díaz-Guilera, "Statistical mechanics of multiedge networks" Phys. Rev. E 88, 062806 (2013)
 * 
 * Author: Oleguer Sagarra, 2013.
 * 
 * 
 *	Output:
 *		
 *		*.tr: A file stored in list format of flows between given nodes
 * 			node i node j t_ij
 *		*.hist: Diverse statistics on the network (indicative names on the files)
 *		*.list: Node attributes average with average deviation (see files).
 * 
 * ******************************************************************************/



#include "main.h"


int main(int argc, char *argv[]){

	printf(\
	"################################\n"
	"########## Maximum entropy multi-edge networks: Fixed degrees, fixed strengths ###########\n"
	"####################################\n");
 /***********************************************************************
	 we read the parameters and initialize the random generator	 N_nodes, Reps, Seed, s_sequence file, d_file),
 ************************************************************************/
	int N_nodes		= -1;
	int T	 		= -1;
	char* file_s;
	int seed		= 1;      			
    int opt_dir		= -1;			
    double bin_exp 	= -1;				
	int verbose		= 0;				
	int opt_clust 	= 0;			
	int self_opt	= 1;
	int header		= 1;
    int print_tr	=0;
    int reps		=100;
	int max_reps	= 100;
	
	int ch;
	        while ((ch = getopt(argc, argv, "N:s:d:f:x:v:c:l:h:T:r:p:e:")) != -1) {
	             switch (ch) {
	             case 's': /* seed */
	                     seed=atoi(optarg);
	                     break;
			     case 'N': /* N_nodes */
			             N_nodes=atoi(optarg);
						 break;
			     case 'p': /* print opt */
			             print_tr=atoi(optarg);
			             break;
				 case 'T': /* TOtal trips */
				 		 T=atoi(optarg);
			             break;
				case 'r': /* reps */
			             reps=atoi(optarg);
					     break;
	            case 'd': /* dir_opt */
	                     opt_dir=atoi(optarg);
	                     break;
		 	    case 'e': /* max_reps */
		 	             max_reps=atoi(optarg);
		 	             break;
	             case 'f': /* file adj */
	                     file_s=optarg;
	                     break;
	             case 'x': /* exp */
	                     bin_exp=atoi(optarg);
	                     break;
			     case 'v': /* verbose */
			             verbose=atoi(optarg);
			             break;
			     case 'c': /* clust */
			             opt_clust=atoi(optarg);
			             break;
			     case 'l': /* self loop */
			             self_opt=atoi(optarg);
			             break;
				case 'h': /* self loop */
					     header=atoi(optarg);
					     break;
	             default:
				 {
	                     printf("Unknown flag %c\n", ch);
	                     exit(EXIT_FAILURE);
	             }
	}
	}
	if((N_nodes<0)||(opt_dir<0)|| (T<0))
	{
 		fprintf(stderr,	"\nCorrect usage is: ./simus -args \n\nWhere:\n\n"
 				" *  Compulosry items:\n"
 				" *		-N N_nodes. Number of nodes (int)\n"
 				" *		-d dir_opt. Undirected (0) or Directed (1)\n"
 			    " *		-f file_xyzw Path to file with with lagrange multipliers id x y z w\n"
				" *		-T Number of trips (int)  \n"
// 			    " *		-a file_dist Path to distance list in format node_i node_j d_ij\n"
 				" *  Optional items: \n"
				" *		-e Max_reps: Maximum repetitions to generate zero inflated poisson [default=100] \n"
 				" *		-s seed.initial seed for random generator (int) (default=1)\n"
 				" *		-x Exponent for log-binning (-1 for no log binning) (Default=-1)\n"
 				" *		-v Verbose (1 for on, 0 for off) (Default 0)\n"
 				" *		-c Clustering option (1 for yes) (warning: Depending on av_s makes simulations orders of magnitude slower) (Default=0)\n"
 				" *		-l Self-loop option (>0 for accepting them) (Default =1) \n"
				" *		-h Number of header lines in input files (default=1)\n"
				" *		-r Number of reps for averaging (default=100)\n"
//				" *		-m Maximum distance for binning (default= 20000) [in meters]\n"
				" *		-p Printing option of sample adjacency list (default=0)\n\n"
 				"Please, read the DOCS/README file for more info!\n");
 		return 0;
 	}

	/****** Check all in params are good ******/
	if(bin_exp<=1) bin_exp=1.05;
    if((opt_dir!=1)&&(opt_dir!=0))
    {
		printf("Select directed or undirected!, aborting...\n");
		abort();
	}
	if((print_tr != 0) && (print_tr != 1))
	{
		printf("Select if you want to print the adj list, aborting ....\n");
		abort();
	}
	if(opt_clust==1)
	{
		if(opt_dir==1)
		{
			printf("Ignoring clustering option, only defined for undirected networks ....\n");
			opt_clust=0;
		}else{
			printf("CLustering option selected! This may cause low performance for high <s>! ....\n");
		}
	}
  /*** Set rand generator (uses GLS THAU) ***/ 
	//seed = seed + time(NULL); // Change for trully random generation
	gsl_rng * randgsl = gsl_rng_alloc(gsl_rng_taus);	/// we initialize the random generator of the gsl
	gsl_rng_set(randgsl,seed);

 /*** Print some info ***/

	printf("SEED=%i\n",seed);
	printf("Average over %d reps\n", reps);
	//printf("Iterations=%i\n",Reps);


	double X;
	//long double Tfake;
	double ** x2;
	double ** x;
//	double** dist;
	double av_k;
 /***********************************************************************
 	Allocating memory + reading distro
 ************************************************************************/ 	
	//dist = read_distances(file_d,N_nodes,header,opt_log);
	if(opt_dir==1)
	{
		x2 = read_node_list_xatts_double(file_s, N_nodes, 4, header);
		X=sum_vec_double(x2[0],N_nodes); 
	}else{
		x = read_node_list_xatts_double(file_s, N_nodes, 2, header);
		X=sum_vec_double(x[0],N_nodes); 
	}
	av_k=(double)T/(double)N_nodes;
	//Tfake = (long int)T;
	if(verbose==1)
	{
		if(opt_dir==1)
		{
			printf("-  x_max/X : (out) =%lf (in) %lf\n",max_value_double(x2[0],N_nodes)/X, max_value_double(x2[1],N_nodes)/X );
		}else{
			printf("-  x_max/X = %lf\n",max_value_double(x[0],N_nodes)/X );
		}
		if(self_opt>0)
		{
			printf("Warning! Self-loops accepted! \n");
			fflush(stdout);
		}
	}
	char cadena[100];
/***********************************************************************
	Preparing ensemble reps!
 ************************************************************************/ 	
	int r, len_acc_nodes;
	W_GRAPH* WG;
	double ** node_cont; //node info container for averaging
	double ** node_cont2; // node info std for averaging
	int ** node_nonzero; // node count non-zero strength/degree
	double * Tcont; // Total T container (av and std)
	int w_max;
	double xmax,ymax,zmax,wmax;
	//int s_min;
	if(opt_dir==1)
	{
		len_acc_nodes=14; // in-out variables
		xmax = max_value_double(x2[0],N_nodes);
		ymax = max_value_double(x2[1],N_nodes);
		zmax = max_value_double(x2[2],N_nodes);
		wmax = max_value_double(x2[3],N_nodes);
		w_max=(int)(50*xmax*ymax*exp(xmax*ymax)*zmax*wmax/(1+zmax*wmax*(exp(xmax*ymax)-1)));
	}else{
		xmax = max_value_double(x[0],N_nodes);
		zmax = max_value_double(x[1],N_nodes);
		w_max=(int)(50*xmax*xmax*exp(xmax*xmax)*zmax*zmax/(1+zmax*zmax*(exp(xmax*xmax)-1)));
		if(opt_clust==1)
		{
			len_acc_nodes=9;
		}else{
			len_acc_nodes=7;
		}
	}
	w_max=w_max*2;
	if(w_max<1)
	{
		w_max=100;
	}
	if(verbose==1)
	{
		printf("-- wmax for histogram: %d\n",w_max);
	}
	node_cont=cast_mat_double(N_nodes,len_acc_nodes);
	node_cont2=cast_mat_double(N_nodes,len_acc_nodes);
	node_nonzero=cast_mat_int(N_nodes,2); // to count the number of times a node is non-connected (in-out)
	
	Tcont=cast_vec_double(2); // mean and std deviation
	//double* entropy_seq = cast_vec_double(reps); // entropy histogram
	//gsl_histogram ** acc_ensemble= w_graph_dist_all_stats_ensemble_allocate(opt_dir, dmax, w_max);
	gsl_histogram ** acc_ensemble= w_graph_all_stats_ensemble_allocate(opt_dir, 1, w_max*w_max/T, 1, w_max*w_max/T, w_max);
		
	/***********************************************************************
		Fill analytic calculation for degree (and optionally clustering)
	 ************************************************************************/ 	
	int i;
	if(opt_dir==1)
	{
		//double** k=w_graph_compute_k_analitic_from_s_directed(xx2,N_nodes, self_opt);
		//double** k=w_graph_dist_compute_k_analitic_from_xygamma_directed(x2,N_nodes, self_opt,gamma,dist);
		double ** k = cast_mat_double(2,N_nodes);
		for(i=0;i<N_nodes;i++)
		{
			node_cont[i][1]=k[0][i]*reps;
			node_cont2[i][1]=k[0][i]*k[0][i]*reps;
			node_cont[i][7]=k[1][i]*reps;
			node_cont2[i][7]=k[1][i]*k[1][i]*reps;
		}
		free_mat_double(k,2);
	}else{
		//double * k=w_graph_compute_k_analitic_from_s_undirected(xx,N_nodes, self_opt);
		//double * k=w_graph_dist_compute_k_analitic_from_xgamma_undirected(x,N_nodes, self_opt,gamma,dist);
		double * k = cast_vec_double(N_nodes);
		// needs new function (compute_k_analitic_from_hidden) <- gamma
		for(i=0;i<N_nodes;i++)
		{
			node_cont[i][1]=k[i]*reps;
			node_cont2[i][1]=k[i]*k[i]*reps;
		}
		free(k);
	}

	
/***********************************************************************
	Start of ensemble reps!
 ************************************************************************/

	for(r=0;r<reps;r++)
	{
		if(verbose==1) printf("============### ZIP model ####===========\n"); fflush(stdout);
		if(opt_dir==1)
		{
			WG = fixedks_poisson_directed_graph(x2, N_nodes, randgsl, verbose, self_opt, max_reps);
		}else{
			WG = fixedks_poisson_undirected_graph(x, N_nodes, randgsl, verbose, self_opt, max_reps);
		}
		if(r==0) // if first rep, store all stats
		{
			/// extended stats with distance ////
			if (verbose>0) printf("... Graph stats ... \n");
			//w_graph_dist_all_stats(node, N_nodes, 0,  bin_exp, av_k, opt_dir,self_opt,dist,randgsl,dmax);
			w_graph_all_stats(WG, N_nodes, r, bin_exp, av_k, opt_dir, self_opt, -1);
			if (verbose>0) printf("... Node stats ... \n");
			//w_graph_node_stats_list(node,N_nodes,0, av_k, opt_dir, opt_clust, self_opt);
			w_graph_node_stats_list(WG,N_nodes,0, av_k, opt_dir, opt_clust, self_opt);
            if(print_tr==1)
            {
                if(verbose>0)
                {
                    printf("Printing adj matrix\n");
                    fflush(stdout);
                }
                sprintf(cadena,"N%d_cust.tr",N_nodes);
                w_graph_print_adj_list(WG, N_nodes, cadena);
            }
		}
		if (verbose>0) printf("... Ensemble stats ... \n");
		w_graph_node_stats_ensemble(WG,N_nodes,node_cont,node_cont2,node_nonzero, Tcont,opt_dir,opt_clust);
		w_graph_all_stats_ensemble_update(acc_ensemble, WG, N_nodes, opt_dir);
		//w_graph_dist_all_stats_ensemble_update(acc_ensemble,node, N_nodes, opt_dir,dist);
		//entropy_seq[r]= w_graph_entropy(node,N_nodes);
		w_graph_free_nodes(WG, N_nodes);
	}
/***********************************************************************
	 Print output
************************************************************************/ 	
	int len=6;
	/*
	if(opt_dir==1)
	{
		len=10;
	}else{
		len=6;
	}
	*/
	/// extended stats with distance ////
	//fflush(stdout);
	if (verbose>0) printf("... Printing ensemble net stats ... \n");
	w_graph_all_stats_ensemble_print(acc_ensemble, 2, reps, N_nodes, av_k, opt_dir);
	//w_graph_dist_all_stats_ensemble_print(acc_ensemble, len, reps, N_nodes,av_k,opt_dir);
	if (verbose>0) printf("... Printing ensemble node stats ... \n");	
	w_graph_node_stats_ensemble_print(reps, N_nodes, Tcont, node_cont, node_cont2, node_nonzero, av_k, bin_exp,len_acc_nodes, opt_dir);
/*
	if(reps>10)
	{
		sprintf(cadena,"N%davs%8.5fentropies.hist",N_nodes,av_k);
		w_graph_print_entropy(entropy_seq,reps,cadena);
	}else{
		printf("If you want entropy histogram, do at least 10 reps!\n");
	}
	free(entropy_seq);
*/	gsl_rng_free (randgsl);
	return 0;
}
