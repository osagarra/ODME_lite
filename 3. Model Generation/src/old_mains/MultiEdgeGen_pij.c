/********************************************************************************
 * 
 *	GENERATE Maximum entropy ensembles for networks
 *
 *	This program generates networks using a grand canonical/canonical approach on a network of N nodes, with a given
 *  pij matrix of size NxN.
 *   pij: t_ij/T [needs to be normalized]
 *  
 *  More details in
 * [1] O. Sagarra, C. J. Pérez Vicente, and A. Díaz-Guilera, "Statistical mechanics of multiedge networks" Phys. Rev. E 88, 062806 (2013)
 * 
 * 
 * 
 * Author: Oleguer Sagarra, 2014.
 * 
 *  Inputs:
 * 		*.dists: A distance file on format node_ori [int] node_dest [int] dist [float]
 *		*.pij : A pij file on format: node_id [int] node_id_dest [int] pij [float]
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
	"########## Maximum entropy multi-edge networks: Fixed collection p_ij ###########\n"
	"####################################\n");
 /***********************************************************************
	 we read the parameters and initialize the random generator	 N_nodes, Reps, Seed, s_sequence file, d_file),
 ************************************************************************/
	int N_nodes		= -1;
	int T	 		= -1;
	char* file_s;
	char* file_d;
	char* file_ori;
    int meth		= 2;
	int seed		= 1;      			
    int opt_dir		= -1;			
    double bin_exp 	= -1;				
	int verbose		= 0;				
	int opt_clust 	= 0;			
	int self_opt	= 1;
	int header		= 1;
	double dmax 		= 20000;
    int print_tr	=0;
    int reps		=100;
	int opt_log		= 0; // logarithmic option
	int opt_soren	=-1;
	
	int ch;
	        while ((ch = getopt(argc, argv, "N:s:d:f:a:x:v:c:l:h:m:T:e:r:p:L:O:Q:")) != -1) {
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
	   			case 'e': /* method */
	   				 	 meth=atoi(optarg);
	   			         break;
				case 'r': /* reps */
			             reps=atoi(optarg);
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
				case 'O': /* Sorensen calculation */
						opt_soren = atoi(optarg);
						break;
				case 'Q': /* Path to original file */
						file_ori = optarg;
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
	if((N_nodes<0)||(opt_dir<0)|| (T<0))
	{
 		fprintf(stderr,	"\nCorrect usage is: ./simus -args \n\nWhere:\n\n"
 				" *  Compulosry items:\n"
 				" *		-N N_nodes. Number of nodes (int)\n"
 				" *		-d dir_opt. Undirected (0) or Directed (1)\n"
 			    " *		-f file_pij Path to file with with Pij matrix in format: node_id_ori node_id_dest p_ij\n"
 			    " *		-a file_dist Path to distance list in format node_i node_j d_ij\n"
				" *		-T Number of trips (int) \n"
 				" *  Optional items: \n"
				" *		-e Ensemble method (0 canonical, 1 Grand Canonical Multinomial, 2 Grand Canonical POisson) (default=2)\n"				
 				" *		-s seed.initial seed for random generator (int) (default=1)\n"
 				" *		-x Exponent for log-binning (-1 for no log binning) (Default=-1)\n"
 				" *		-v Verbose (1 for on, 0 for off) (Default 0)\n"
 				" *		-c Clustering option (1 for yes) (warning: Depending on av_s makes simulations orders of magnitude slower) (Default=0)\n"
 				" *		-l Self-loop option (>0 for accepting them) (Default =1) \n"
				" *		-h Number of header lines in input files (default=1)\n"
				" *		-m Maximum distance for binning (default= 20000) [in meters]\n"
				" *		-r Number of reps for averaging (default=100)\n"
				" *		-p Printing option of sample adjacency list (default=0)\n"
				" *		-L Log-dist option (to compute the logarithm of the cost matrix) [default=0]\n"
				" *		-Q Set >0 if you want to calculate average indices comparing to real data [default=-1]\n"
				" *		-o If -Q>1, indicate path to original data in adj list format node_ori node_dest t_ij [all integers]\n\n"
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
	if((meth>2) || (meth<0))
	{
		printf("Select apropiate method (canonical, grandcanonical...), aborting ....\n");
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
	if(seed==1)
	{
		seed = seed + time(NULL); // Change for trully random generation
		if(verbose>1) printf("Default seed given, using processor time...\n");
	}
	gsl_rng * randgsl = gsl_rng_alloc(gsl_rng_taus);	/// we initialize the random generator of the gsl
	gsl_rng_set(randgsl,seed);
 /*** Print some info ***/

	printf("SEED=%i\n",seed);
	printf("Average over %d reps\n", reps);
	printf("Total number of trips:%d\n", T);
	//printf("Iterations=%i\n",Reps);

	//double ** x2 = cast_mat_double(2,2); // dummy in this case
	//double * x = cast_vec_double(1); // dummy in this case
	double** dist;
	double** pij;
	double * pij_flat;
	double av_k;
 /***********************************************************************
 	Allocating memory + reading distro
 ************************************************************************/ 	
	pij = read_net_list_double(file_s, N_nodes, header);
	double norm = sum_matrix_double(pij, N_nodes, N_nodes);
	scale_matrix(pij, N_nodes, N_nodes, 1./norm);
	dist = read_distances(file_d,N_nodes,header,opt_log);
	av_k=(double)T/(double)N_nodes;
	char cadena[100];
/***********************************************************************
	Preparing ensemble reps!
 ************************************************************************/ 	
	int r, len_acc_nodes;
	int L=0; // dummy length
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

	}else{
		if(opt_clust==1)
		{
			len_acc_nodes=9;
		}else{
			len_acc_nodes=7;
		}
	}
	w_max= (T*matrix_max_value_double(pij, N_nodes,N_nodes));
	w_max=w_max*100;
	if(w_max<1)
	{
		w_max=100;
	}
	if(verbose==1)
	{
		printf("-- wmax for histogram: %d\n",(int)w_max);
	}
	node_cont=cast_mat_double(N_nodes,len_acc_nodes);
	node_cont2=cast_mat_double(N_nodes,len_acc_nodes);
	node_nonzero=cast_mat_int(N_nodes,2); // to count the number of times a node is non-connected (in-out)
	
	Tcont=cast_vec_double(2); // mean and std deviation
	double* entropy_seq = cast_vec_double(reps); // entropy histogram
	gsl_histogram ** acc_ensemble= w_graph_dist_all_stats_ensemble_allocate(opt_dir, dmax, (int)w_max);
		
	/***********************************************************************
		Fill analytic calculation for degree (and optionally clustering)
	 ************************************************************************/ 	
	int i;
	if(opt_dir==1)
	{
		//double** k=w_graph_compute_k_analitic_from_s_directed(xx2,N_nodes, self_opt);
		//double** k=w_graph_dist_compute_k_analitic_from_xygamma_directed(x2,N_nodes, self_opt,gamma,dist);
		// not defined for this case ... //
		for(i=0;i<N_nodes;i++)
		{
			node_cont[i][1]=0;
			node_cont2[i][1]=0;
			node_cont[i][8]=0;
			node_cont2[i][8]=0;
		}
		//free_mat_double(k,2);
	}else{
		//double * k=w_graph_compute_k_analitic_from_s_undirected(xx,N_nodes, self_opt);
		//double * k=w_graph_dist_compute_k_analitic_from_xgamma_undirected(x,N_nodes, self_opt,gamma,dist);
		// needs new function (compute_k_analitic_from_hidden) <- gamma
		// not defined for this case ... //
		for(i=0;i<N_nodes;i++)
		{
			node_cont[i][1]=0;
			node_cont2[i][1]=0;
		}
		//free(k);
	}
	// scale pij to actual values if poisson!
	if(meth==2)
	{
		scale_matrix(pij,N_nodes,N_nodes,(double)T);
	}else{ // flatten matrix (normalized)
		pij_flat=flatten_matrix_double(pij, N_nodes, N_nodes, self_opt, &L);
		free_mat_double(pij,N_nodes);
	// if canonical, flatten pij	
	}
	
/***********************************************************************
	Start of ensemble reps!
 ************************************************************************/
	double soren = 0;
	double soren2 = 0;
	double sorfake;
	W_GRAPH* ori;
	if(opt_soren>0)
	{
		ori = w_graph_dist_read_edge_list(file_ori, N_nodes, opt_dir,header);
	}
	for(r=0;r<reps;r++)
	{
		if(meth==2)
		{
			if(verbose==1) printf("============### Poisson model ####===========\n"); fflush(stdout);
			if(opt_dir==1)
			{
				WG = custompij_poisson_directed_graph(pij, N_nodes, randgsl, verbose, self_opt);
			}else{
				WG = custompij_poisson_undirected_graph(pij, N_nodes, randgsl, verbose, self_opt);
			}
		}else{
			if(meth==0)
			{
				if(verbose==1)printf("============### Multinomial model ####===========\n"); fflush(stdout);
				if(opt_dir==1)
				{
					WG = multinomial_directed_graph(pij_flat, N_nodes, T, randgsl,verbose, self_opt);
				}else{
					WG = multinomial_undirected_graph(pij_flat, N_nodes, T, randgsl,verbose, self_opt);
				}
			}else{
				if(verbose==1)printf("============### Poisson Multinomial model ####===========\n"); fflush(stdout);
				if(opt_dir==1)
				{
					WG = poisson_multinomial_directed_graph(pij_flat, N_nodes, T, randgsl, verbose, self_opt);
				}else{
					WG = poisson_multinomial_undirected_graph(pij_flat, N_nodes, T, randgsl, verbose, self_opt);
				}
			}
		}
		if(r==0) // if first rep, store all stats
		{
			/// extended stats with distance ////
			if (verbose>0) printf("... Graph stats ... \n");
			w_graph_dist_all_stats(WG, N_nodes, 0,  bin_exp, av_k, opt_dir,self_opt,dist,randgsl,dmax);
			if (verbose>0) printf("... Node stats ... \n");
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
			/// extended stats with distance ////
		w_graph_node_stats_ensemble(WG,N_nodes,node_cont,node_cont2,node_nonzero, Tcont,opt_dir,opt_clust);
		w_graph_dist_all_stats_ensemble_update(acc_ensemble,WG, N_nodes, opt_dir,dist);
		entropy_seq[r]= w_graph_entropy(WG,N_nodes);
		if(opt_soren>0)
		{
			sorfake = w_graph_compute_sorensen(WG, ori, N_nodes);
			soren += sorfake;
			soren2 += sorfake*sorfake;
		}
		w_graph_free_nodes(WG, N_nodes);
	}
	//printf("i am here!\n");fflush(stdout);
	if(meth==2)
	{
		free(pij);		
	}else{
		free(pij_flat);
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
	if (verbose>0) printf("... Printing ensemble stats ... \n");
	w_graph_dist_all_stats_ensemble_print(acc_ensemble, len, reps, N_nodes,av_k,opt_dir);
	w_graph_node_stats_ensemble_print(reps, N_nodes, Tcont, node_cont, node_cont2, node_nonzero, av_k, bin_exp,len_acc_nodes, opt_dir);
	if(reps>10)
	{
		sprintf(cadena,"N%davs%8.5fentropies.hist",N_nodes,av_k);
		w_graph_print_entropy(entropy_seq,reps,cadena);
	}else{
		printf("If you want entropy histogram, do at least 10 reps!\n");
	}
	if (opt_soren>0)
	{
		soren=soren/(double)reps;
		printf("... Average Sorensen:%f+-%f ... \n",soren,sqrt(soren2/(double)reps - (soren*soren)));
		w_graph_free(ori,N_nodes);
	}
	gsl_rng_free (randgsl);
	free(entropy_seq);
	return 0;
}
