/********************************************************************************
 * 
 * Network analyzer
 * 
 * Author: Oleguer Sagarra, 2014.
 * 
 * 
 *	Output:
 *		
 *		*.hist: Diverse statistics on the network (indicative names on the files)
 *		*.list: Node attributes average with average deviation (see files).
 * 
 * ******************************************************************************/



#include <stdio.h>
#include <math.h>
#include "main.h"


int main(int argc, char *argv[]){
	
	printf(\
	"################################\n"
	"########## Multi-Edge Analyzer ###########\n"
	"####################################\n");
 /***********************************************************************
	 we read the parameters and initialize the random generator	 N_nodes, Reps, Seed, s_sequence file, d_file),
 ************************************************************************/
	int  	N_nodes=-1	;
	char* file_s;
	char* file_d;
	int opt_dist = 0;
	int  	seed=1;      			
    int opt_dir=-1;			
    double bin_exp = -1;				
	int verbose	= 0;				
	int opt_clust = 0;			
	int self_opt = 1;
	int header=1;
	double dmax = 20000;
	int opt_log		= 0; // logarithmic option
    char* xypath; // path to lagrange multipliers
    int opt_filter = 0; // filter option
    int cases = 0 ; // default ME
    double ci = 0.95; // confidence level
    int M=1; // layers
	
	int ch;
	        while ((ch = getopt(argc, argv, "N:s:d:f:a:x:v:c:l:h:m:z:L:F:C:M:Z:X:")) != -1) {
	             switch (ch) {
				 case 'z': /*distance analysis */
				 		 opt_dist=atoi(optarg);
						 break;
	             case 's': /* seed */
	                     seed=atoi(optarg);
	                     break;
			     case 'N': /* N_nodes */
			             N_nodes=atoi(optarg);
			             break;
	             case 'd': /* dir_opt */
	                     opt_dir=atoi(optarg);
	                     break;
	             case 'f': /* file adj */
	                     file_s=optarg;
	                     break;
			     case 'a': /* file Dist */
				 		 //opt_dist=1;
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
                case 'F': /* option for filtering*/
                         opt_filter = atoi(optarg);
                         break;
                case 'C': /* case*/
                         cases = atoi(optarg);
                         break;
                case 'M': /* layers*/
                         M = atoi(optarg);
                         break;
                case 'Z': /*confidence interval*/
                         ci = atof(optarg);
                         break;
                case 'X': /*lagrange mult path */
                         xypath =optarg;
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

	if((N_nodes<0)||(opt_dir<0))
	{
 		fprintf(stderr,	"\nCorrect usage is: ./MultiEdgeAnalyzer -args \n\nWhere:\n\n"
 				" *  Compulosry items:\n"
 				" *		-N N_nodes. Number of nodes (int)\n"
 				" *		-d dir_opt. Undirected (0) or Directed (1)\n"
 			    " *		-f file_adj Path to file with adj format node_i node_j t_ij\n"
 				" *  Optional items: \n"				
				" *		-z Distance_option: Include distances in analysis? [default=0, 1 for yes]\n"
	 			" *		-a file_dist Path to distance list in format node_i node_j d_ij [default, no distance]\n"
 				" *		-s seed.initial seed for random generator (int) (default=1)\n"
 				" *		-x Exponent for log-binning (-1 for no log binning) (Default=-1)\n"
 				" *		-v Verbose (1 for on, 0 for off) (Default 0)\n"
 				" *		-c Clustering option (1 for yes) (warning: Depending on av_s makes simulations orders of magnitude slower) (Default=0)\n"
 				" *		-l Self-loop option (>0 for accepting them) (Default =1) \n"
				" *		-h Number of header lines in input files (default=1)\n"
				" *		-m Maximum distance for binning (default= 20000) [in meters]\n"
				" *		-L Log-dist option (to compute the logarithm of the cost matrix) [default=0]\n"
                " *     -F Implement Graph-filtering (applies Filter to graph according to null model with fixed strengths) [default=0]\n"
                " *         -C case Null model type (ME=0, B=1, W=2) [default=0]\n"
                " *         -M layers Number of layers [default=1] \n"
                " *         -Z alpha Confidence level [default=0.95] \n"
                " *         -X Path to file containing lagrange multiplier values xy\n"
 				"Please, read the DOCS/README file for more info!\n");
 		return 0;
 	}
	/****** Check all in params are good ******/
	if(opt_dist>0)
	{
		printf("Distance given, spatial network analysis...\n");
	}else
	{
		printf("No distance given, regular network analysis...\n");		
	}
	if(bin_exp<=1) bin_exp=1.05;
    if((opt_dir!=1)&&(opt_dir!=0))
    {
		printf("Select directed or undirected!, aborting...\n");
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
    if(opt_filter>0)
    {
        printf("Applying filter in analysis...\n");
        if(cases==0)
        {
            printf("\tME Case selected\n");
        }else if(cases==1)
        {
            printf("\tB Case selected\n");
        }else if(cases==2){
            printf("\tW Case selected\n");
        }else{
            printf("\tCase must be ME (0), B(1) or W(2). Aborting... \n");
            abort();
        }
        if(M<0)
        {
            printf("\tLayers must be positive! Aborting...\n");
            abort();
        }else{
            printf("\tSelected %d layers\n",M);
        }
        //if((ci<=0) ||(ci>=1))
        if((ci<=0) ||(ci>1))
        {
            printf("\tC.I. cannot be larger or equal than 1 or smaller than 0! Aborting...\n");
            abort();
        }else{
            printf("\tSelected C.I of %.5f\n",ci);
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
	//printf("SEED=%i\n",seed);
	//printf("Iterations=%i\n",Reps);


	double** dist;
	double av_k;
	int* xx;
	int** xx2;
	int T;
    char cadena[100]; // to print files
 /***********************************************************************
 	Allocating memory + reading distro
 ************************************************************************/ 	
	if(opt_dist>0)
	{
		dist = read_distances(file_d,N_nodes,header,opt_log, verbose);
	}
	W_GRAPH* WG;
    W_GRAPH* WGf;
    double** x2 = cast_mat_double(2,N_nodes);
    if(opt_filter>0)// if necessary, apply filter
    {
        WGf = w_graph_read_edge_list(file_s, N_nodes, opt_dir,header, verbose);
        // load xy's
        if(cases==0)
        {
            //xx2 = w_graph_compute_s(WGf, N_nodes);
            //x2[0] = vec_int_to_double(xx2[0], N_nodes);
            //x2[1] = vec_int_to_double(xx2[1], N_nodes);
            //T = w_graph_total_weight(WGf, N_nodes);
            //scale_matrix(x2, 2, N_nodes, 1./sqrt((double)T));
            x2 = read_node_list_xatts_double(xypath, N_nodes, 2, header, verbose);
        }else{
            x2 = read_node_list_xatts_double(xypath, N_nodes, 2, header, verbose);
        }
        printf("Filtering graph (may take a while)...\n");
        if(opt_dir==0)
        {
            WG = w_graph_filter_xij(WGf, x2[0], x2[0], N_nodes, ci, cases, M, verbose);
        }else{
            WG = w_graph_filter_xij(WGf, x2[0], x2[1], N_nodes, ci, cases, M, verbose);
        }
        free_mat_double(x2,2);
		T = w_graph_total_weight(WG, N_nodes);
        //abort(); // only for special case!!!
		if(T>0)
		{
			printf("\tPrinting filtered adj matrix\n");
			sprintf(cadena,"N%d_filteredCI%.3f.tr",N_nodes,ci);
			w_graph_print_adj_list(WG, N_nodes, cadena);
		}else{
			printf("\t Filter returned an empty net!\n");
			abort();
		}
        //abort(); // only for special case!!!
    }else{
        WG = w_graph_read_edge_list(file_s, N_nodes, opt_dir,header, verbose);
    } 
	xx2 = w_graph_compute_s(WG, N_nodes);
	T = w_graph_total_weight(WG, N_nodes);
	av_k=(double)T/(double)N_nodes;

	//Tfake = (long int)T;
	if(verbose==1)
	{
		if(opt_dir==1)
		{
			printf("-  x_max/X : (out) =%lf (in) %lf\n",max_value_int(xx2[0],N_nodes)/(double)T, max_value_int(xx2[1],N_nodes)/(double)T );
		}else{
			xx = xx2[0];
			printf("-  x_max/X = %lf\n",max_value_int(xx,N_nodes)/(double)T );
		}
	}
/***********************************************************************
	Analizing network
 ************************************************************************/ 	
	int w_max;
	//int s_min;
	int E;
	int * w = w_graph_compute_w(WG, N_nodes, &E, -1);
	w_max= max_value_int(w,E);
	if (verbose>0) printf("-- wmax for histogram: %d\n",w_max);
	/// extended stats with distance ////
	if (verbose>0) printf("... Graph stats ... \n");
	if (opt_dist>0)
	{
		w_graph_dist_all_stats(WG, N_nodes, 0,  bin_exp, opt_dir, self_opt,dist,randgsl,dmax, verbose);
	}else{
		w_graph_all_stats(WG, N_nodes, 0, bin_exp, opt_dir, self_opt, 0, verbose);
	}
	if (verbose>0) printf("... Node stats ... \n");
	w_graph_node_stats_list(WG,N_nodes,0, opt_dir, opt_clust, self_opt, verbose);
	if (verbose>0) printf("# ML ensemble Multi-Edge Net entropy per event: %f\n",w_graph_ML_ensemble_entropy_multinomial(WG,N_nodes,opt_dir));
	w_graph_free_nodes(WG, N_nodes);
	free(WG);
}
