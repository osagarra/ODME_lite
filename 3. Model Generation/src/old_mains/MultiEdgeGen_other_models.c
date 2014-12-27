/********************************************************************************
 * 
 *			Network randomizer
 *
 * Randomizes given network using radiation or sequential gravity model
 * Seq:
 * 	[1] Lenormand M, Huet S, Gargiulo F, Deffuant G (2012) A Universal Model of
 *    Commuting Networks. PLoS ONE 7(10): e45985. doi:10.1371/journal.pone.0045985
 * Rad:
 *  [2] F. Simini, M. C. González, A. Maritan, and A.-L. Barabási, Nature 484, 96 (2012).
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
	"########## Multi-edge network Randomizer: Non entropy models ###########\n"
	"####################################\n");
 
 /***********************************************************************
	 we read the parameters and initialize the random generator	 N_nodes, Reps, Seed, s_sequence file, d_file),
 ************************************************************************/
	int N_nodes		= -1;
	char* file_s;
	char* file_d;
	int meth		= -1;
	int seed		= -32;      			
    int opt_dir		= -1;			
    double bin_exp 	= -1;				
	int verbose		= 0;				
	int opt_clust 	= 0;			
	int self_opt	= 1;
	int header		= 1;
	double dmax 	= 20000;
    int print_tr	=0;
	double gamma	=-1;
    int reps		=100;
	int opt_log		= 0; // logarithmic option

	int ch;
	        while ((ch = getopt(argc, argv, "N:s:d:f:a:x:v:c:l:h:m:g:r:p:e:L:")) != -1) {
	             switch (ch) {
	             case 's': /* seed */
	                     seed=atoi(optarg);
	                     break;
			     case 'e': /* method */
			             meth=atoi(optarg);
			             break;
			     case 'N': /* N_nodes */
			             N_nodes=atoi(optarg);
						 break;
			     case 'p': /* print opt */
			             print_tr=atoi(optarg);
			             break;
				case 'r': /* reps */
			             reps=atoi(optarg);
					     break;
				case 'g': /* reps */
			             gamma=atof(optarg);
					     break;		 
	            case 'd': /* dir_opt */
	                     opt_dir=atoi(optarg);
	                     break;
	             case 'f': /* file adj */
	                     file_s=optarg;
	                     break;
			     case 'a': /* file Dist */
			             file_d=optarg;
			             break;
	             case 'x': /* exp */
	                     bin_exp=atof(optarg);
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
				case 'm': /* max_Distance */
						 dmax = atof(optarg);
						 break;
		 		case 'L': /* Log opt */
		 				opt_log=atoi(optarg);
		 				break;

	             default:
				 {
	                     printf("Unknown flag %c\n", ch);
	                     exit(EXIT_FAILURE);
	             }
	}
	}
	if(opt_log>0)
	{
		dmax = log(dmax);
	}
	if((N_nodes<0)||(opt_dir<0)|| (meth<0))
	{
 		fprintf(stderr,	"\nCorrect usage is: ./simus -args \n\nWhere:\n\n"
 				" *  Compulosry items:\n"
 				" *	       -N N_nodes. Number of nodes (int)\n"
 				" *        -d dir_opt. Undirected (0) or Directed (1)\n"
 			    " *        -f file_adj Path to file with node strength, formaT: node_id node_sout node_sin\n"
 			    " *        -a file_dist Path to distance list in format node_i node_j d_ij (for sequential case needs to be in meters!) \n"
				" *        -e Method: 0 for Radiation model multinomial, 1 for Radiation model stochastic, 2 for seq. gravity model multinomial, 3 for seq. gravity  model bernouilli  \n"
 				" *  Optional items: \n"
 				" *        -s seed.initial seed for random generator (int) (default= Use random seed from processor clock)\n"
 				" *        -x Exponent for log-binning (-1 for no log binning) (Default=-1)\n"
 				" *        -v Verbose (1 for on, 0 for off) (Default 0)\n"
 				" *        -c Clustering option (1 for yes) (warning: Depending on av_s makes simulations orders of magnitude slower) (Default=0)\n"
 				" *        -l Self-loop option (>0 for accepting them) (Default =1) \n"
				" *        -h Number of header lines in input files (default=1)\n"
				" *        -m Maximum distance for binning (default= 20000) [in meters]\n"
				" *        -g Gamma value (for the seq. gravity model only, if not given automatically approximated assuming cost function is given in meters, ignored in the radiation case)\n"
				" *        -r Number of reps for averaging (default=100)\n"
				" *        -p Printing option of sample adjacency list (default=0)\n"
				" *        -h Number of lines to skip while reading the files (header, default=1)\n"
				" *        -L Log-dist option (to compute the logarithm of the cost matrix) [default=0]\n\n"
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
	if (opt_dir!=1)
	{
		printf("Sorry, not implemented for undirected, aborting...\n");
		abort();
	}
	if((meth>3) || (meth<0))
	{
		printf("Select apropiate method rad (0 or 1) or seq. gravity (2)), aborting ....\n");
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
	if(seed==-32)
	{
		printf("Default SEED selected, using random seed from clock\n");
		seed = seed + time(NULL); // Change for trully random generation
	}
	gsl_rng * randgsl = gsl_rng_alloc(gsl_rng_taus);	/// we initialize the random generator of the gsl
	gsl_rng_set(randgsl,seed);
	unsigned long int maxi =  gsl_rng_max (randgsl);
	unsigned long int mini = gsl_rng_min (randgsl);
 /*** Print some info ***/

	printf("SEED=%i\n",seed);
	printf("Average over %d reps\n", reps);
	printf("Max value for RAND generator (taus) : %ld, Min: %ld\n",maxi,mini);


	int T;
	//long double Tfake;
	int** xx2;
	int* xx;
	double av_k;
	double** dist;
	double** ps=cast_mat_double(N_nodes,N_nodes);
	int ps_opt;
 /***********************************************************************
 	Allocating memory + reading distro + fixing gamma
 ************************************************************************/ 	
	dist = read_distances(file_d,N_nodes,header,opt_log);

	if(meth>1)
	{
		printf("Seq. Gravity model selected \n");
		if(gamma<0)
		{
			double* surf = nearest_neighbor_dist(dist,N_nodes); // average
			double av_surf = sum_vec_double(surf,N_nodes)*1./(double)N_nodes; // average inter-nodal nearest distance
			printf("\t Gamma not specified, using default from paper:\n \t Surface calculated using as proxy for the average surface the squared " 
			"average distance between nearest neighbors: %f\n",av_surf); // considering disc
			gamma = compute_gamma_frenchies((av_surf*1e-3/2.)*(av_surf*1e-3/2.)*3.141593); // surface in meters needs to be transformed in km^2
			free(surf);
		}
		printf("\t Selected inverse gamma: %f\n",1./gamma);
	}else{
		printf("Ignoring gamma value, Radiation model selected\n");
	}
	//printf("Iterations=%i\n",Reps);
	if(opt_dir==1)
	{
		xx2 = read_node_list_int(file_s, N_nodes, header); // strenght sequence (ints)
		T=sum_vec_int(xx2[0],N_nodes); // T is \sum_i s_i (for all cases)
	}else{
		xx = read_node_list_int_undir(file_s, N_nodes, header); // strenght sequence (ints)
		T=sum_vec_int(xx,N_nodes); // T is \sum_i s_i (for all cases)
		if(T%2!=0)
		{
			if(verbose>0)
			{
                printf("Warning: T is not even, aborting\n");
                //abort();
			}
		}
	}
	av_k=(double)T/(double)N_nodes;
	//Tfake = (long int)T;
	if(verbose==1)
	{
		if(opt_dir==1)
		{
			printf("- T to be sorted: %d -> <s>=%f || s_max/T : (out) =%f (in) %f\n",T,(double)T/(double)N_nodes, (double)max_value_int(xx2[0],N_nodes)/(double)T, (double)max_value_int(xx2[1],N_nodes)/(double)T );
		}else{
			printf("- T to be sorted: %d -> <s>=%f || s_max/2T = %f\n",T,(double)T/(double)N_nodes, (double)max_value_int(xx,N_nodes)/(double)(T) );
		}
		if(self_opt>0)
		{
			printf("Warning! Self-loops accepted! \n");
			fflush(stdout);
		}
	}
/***********************************************************************
	Preparing ensemble reps!
 ************************************************************************/ 	
	int r, len_acc_nodes;
	W_GRAPH* WG;
	double ** node_cont; //node info container for averaging
	double ** node_cont2; // node info std for averaging
	int ** node_nonzero; // node count non-zero strength/degree
	double * Tcont; // Total T container (av and std)
	double w_max;

	//int s_min;
	if(opt_dir==1)
	{
		len_acc_nodes=14; // in-out variables
		w_max=50*(double)max_value_int(xx2[0],N_nodes)*(double)max_value_int(xx2[1],N_nodes)/T;
	}else{
		w_max=50*(double)max_value_int(xx,N_nodes)*(double)max_value_int(xx,N_nodes)/T;
		if(opt_clust==1)
		{
			len_acc_nodes=9;
		}else{
			len_acc_nodes=7;
		}
	}
	w_max=w_max*2;
	printf("-- wmax for histogram: %d\n",(int)w_max);
	if(w_max<1)
	{
		w_max=10;
	}
	node_cont=cast_mat_double(N_nodes,len_acc_nodes);
	node_cont2=cast_mat_double(N_nodes,len_acc_nodes);
	node_nonzero=cast_mat_int(N_nodes,2); // to count the number of times a node is non-connected (in-out)
	Tcont=cast_vec_double(2); // mean and std deviation
	gsl_histogram ** acc_ensemble= w_graph_dist_all_stats_ensemble_allocate(opt_dir, 40000, (int)w_max);
		
	/***********************************************************************
		Fill analytic calculation for degree (and optionally clustering)
	 ************************************************************************/ 	
	int i;
	if(opt_dir==1)
	{
		for(i=0;i<N_nodes;i++)
		{
			node_cont[i][1]=0;
			node_cont2[i][1]=0;
			node_cont[i][7]=0;
			node_cont2[i][7]=0;
		}
	}else{
		// needs new function (compute_k_analitic_from_hidden) <- gamma
		for(i=0;i<N_nodes;i++)
		{
			node_cont[i][1]=0;
			node_cont2[i][1]=0;
		}
	}

	
/***********************************************************************
	Start of ensemble reps!
 ************************************************************************/
	// w_graph* w_graph_seq_gravity_directed(int N_nodes,int **s, double **d, double gamma, gsl_rng* randgsl, int Max_fails, int self_opt, int verbose){
	for(r=0;r<reps;r++)
	{
		if(r==0)
		{
			ps_opt = -1;
		}else{
			ps_opt =1;
		}
		if(meth==2)
		{
			if(verbose==1) printf("============### Seq model multinomial ####===========\n"); fflush(stdout);
			if(opt_dir==1)
			{
				WG = w_graph_seq_gravity_multinomial_directed(N_nodes, xx2, dist, gamma, randgsl, self_opt, verbose);
			}else{
				printf("Not defined for undirected case... \n"); 
			}
		}else if(meth==3){
			if(verbose==1) printf("============### Seq model bernouilli ####===========\n"); fflush(stdout);
			if(opt_dir==1)
			{
				WG = w_graph_seq_gravity_bernouilli_directed(N_nodes, xx2, dist, gamma, randgsl, self_opt, verbose);
			}else{
				printf("Not defined for undirected case... \n"); 
			}
		}else if (meth==1){
			if(verbose==1)printf("============### Radiation Stochastic model ####===========\n"); fflush(stdout);
			if(opt_dir==1)
			{
				WG =   w_graph_radiation_model_stochastic_directed(N_nodes,xx2, dist, randgsl,self_opt,verbose);
			}else{
				printf("Not implemented for undirected case... \n");
			}
		}else{
			if(verbose==1)printf("============### Radiation Multinomial model ####===========\n"); fflush(stdout);
			if(opt_dir==1)
			{
				WG =   w_graph_radiation_model_multinomial_directed(N_nodes,xx2, dist, randgsl,self_opt,verbose,ps_opt,ps);
			}else{
				printf("Not implemented for undirected case... \n");
			}			
		}
		if(r==0) // if first rep, store all stats
		{
			/// extended stats with distance ////
			if (verbose>0) printf("... Graph stats ... \n");
			w_graph_dist_all_stats(WG, N_nodes, 0,  bin_exp, av_k, opt_dir,self_opt,dist,randgsl,dmax);
			if (verbose>0) printf("... Node stats ... \n");
			w_graph_node_stats_list(WG,N_nodes,0, av_k, opt_dir, opt_clust, self_opt);
            char cadena[100];
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
			/// extended stats with distance ////
		w_graph_node_stats_ensemble(WG,N_nodes,node_cont,node_cont2,node_nonzero, Tcont,opt_dir,opt_clust);
		w_graph_dist_all_stats_ensemble_update(acc_ensemble,WG, N_nodes, opt_dir,dist);
		w_graph_free_nodes(WG, N_nodes);
	}
	free_mat_double(ps,N_nodes);
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
	w_graph_dist_all_stats_ensemble_print(acc_ensemble, len, reps, N_nodes,av_k,opt_dir);
	w_graph_node_stats_ensemble_print(reps, N_nodes, Tcont, node_cont, node_cont2, node_nonzero, av_k, bin_exp,len_acc_nodes, opt_dir);
	gsl_rng_free (randgsl);
	return 0;
}
