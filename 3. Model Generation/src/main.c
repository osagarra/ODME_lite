/********************************************************************************
 * 
 *  GENERATE instances for various Network models
 *
 *  This program generates networks according to different models and cases. Cases are described in the readme file.
 *  
 * 
 * Author: Oleguer Sagarra, 2014.
 * 
 *  Output:
 *      
 *      *.tr: A file stored in list format of flows between given nodes
 *          node i node j t_ij
 *      *.hist: Diverse statistics on the network (indicative names on the files)
 *      *.list: Node attributes average with average deviation (see files).
 * 
 * ******************************************************************************/



#include "main.h"


int main(int argc, char *argv[]){

    printf(\
    "################################\n"
    "########## General Network Model Generator ###########\n"
    "####################################\n");
 /***********************************************************************
     we read the parameters and initialize the random generator  N_nodes, Reps, Seed, s_sequence file, d_file),
 ************************************************************************/
    // params
    int N_nodes         = -1;
    double T            = -1; // total number of events (sometimes ignored) 
    double gamma        = -1; // gamma (lambda for case of fixed k)
    int cases           = -1; // case
    // files
    char* file_s;   // input lagrange multipliers file
    char* file_d; // distance file (if computing distances)
    char* file_ori; // original network (to compute loglikelyhood, entropy, sorensen
    // options
    int opt_dir         = -1;           
    int opt_dist        = 0; // distance option
    int opt_log         = 0; // logarithmic option on distance if on
    int opt_indist      = 0; // indist. option (only for fixed k,s,sk)
    int opt_agg         = 0; // agg. option (only for fixed k,s,sk)
    int opt_clust       = 0; // compute clustering?         
    int opt_self        = 1; // self loops?
    int opt_soren       = 0; // sorensen compute?
    int opt_entropy     = 0; // entropy compute?
    int opt_print_tr    = 0; // print sample of adj matrix?
    int opt_verbose     = 0;                
    // more opts
    int layers          = 1; // number of layers
    int seed            = 1;                
    int meth            = -1;
    double bin_exp      = -1;               
    int header          = 1;
    double dmax         = 20000;
    int reps            = 100;
    int max_reps        = 100;

    int ch;
            //while ((ch = getoapt(argc, argv, "N:T:g:C:f:a:Q:d:z:L:D:c:l:S:p:v:s:e:x:h:m:r:R:")) != -1) {
            while ((ch = getopt(argc, argv, "N:T:g:C:f:a:Q:d:z:L:I:A:c:l:S:p:v:M:s:e:x:h:m:r:R:E:")) != -1) {
                 switch (ch) {
                 case 'N': /* N_nodes */
                         N_nodes=atoi(optarg);
                         break;
                 case 'T': /* TOtal trips */
                         T=atoi(optarg);
                         break;
                 case 'g': /* gamma */
                         gamma=atof(optarg);
                         break;     
                 case 'C': /* case */
                         cases=atoi(optarg);
                         break; 
                 case 'f': /* file adj */
                         file_s=optarg;
                         break;
                 case 'a': /* file Dist */
                         file_d=optarg;
                         break;
                 case 'Q': /* Path to original file */
                         file_ori = optarg;
                         break;                                                                      
                 case 'd': /* dir_opt */
                         opt_dir=atoi(optarg);
                         break;                      
                 case 'z': /*distance analysis */
                         opt_dist=atoi(optarg);
                         break;
                 case 'L': /* Log opt */
                         opt_log=atoi(optarg);
                         break;
                 case 'I': /* Distinguishable option */
                         opt_indist=atoi(optarg);
                         break;
                 case 'A': /* Aggregate option */
                         opt_agg=atoi(optarg);
                         break;
                 case 'c': /* clust */
                         opt_clust=atoi(optarg);
                         break;                      
                 case 'l': /* self loop */
                         opt_self=atoi(optarg);
                         break;                      
                 case 'S': /* Sorensen calculation */
                        opt_soren = atoi(optarg);
                        break;                                                                                               
                 case 'E': /* Entropy calculation */
                        opt_entropy = atoi(optarg);
                        break;                                                                                               
                 case 'p': /* print opt */
                         opt_print_tr=atoi(optarg);
                         break;
                 case 'v': /* opt_verbose */
                         opt_verbose=atoi(optarg);
                         break;
                 case 'M': /* layers */
                         layers = atoi(optarg);
                         break;
                 case 's': /* seed */
                         seed=atoi(optarg);
                         break;
                 case 'e': /* method */
                         meth=atoi(optarg);
                         break;
                 case 'x': /* exp */
                         bin_exp=atof(optarg);
                         break;
                case 'h': /* self loop */
                         header=atoi(optarg);
                         break;
                case 'm': /* max_Distance */
                         dmax = atof(optarg);
                         break;
                 case 'r': /* reps */
                         reps=atoi(optarg);
                         break;
                case 'R': /* max_reps */
                         max_reps=atoi(optarg);
                         break;
                 default:
                 {
                         printf("Unknown flag %c\n", ch);
                         exit(EXIT_FAILURE);
                 }
        }
    }
    if((N_nodes<0)||(opt_dir<0) || (cases<0))
    {
        fprintf(stderr, "\nCorrect usage is: ./GenNetGen -args \n\nWhere:\n\n"
                " *  Compulsory items:\n"
                " *     -N N_nodes. Number of nodes (int)\n"
                " *     -d dir_opt. Undirected (0) or Directed (1)\n"
                " *     -f file_input Path to file with appropiate Lagrange Multipliers in each case or srtength sequence (see README)\n"
                " *     -C Case: Case to apply to the model, see README for more info and description of each case (int)\n"
                " *  Optional (in some cases might be compulsory) items: \n"
                " *     -e Ensemble method to use in each case [only available in some cases, see README] \n"
                " *     -g Gamma Value: Deferrence value for distances or number of edges (according to each case) \n"
                " *     -T Number of trips (if needed, see DOCS) (int) \n"
                " *     -z distance option, If >0, compute cost statistics with given cost matrix (int, default = 0) \n"
                " *         -a Path to cost matrix list in format node_i node_j d_ij [must be triangular matrix symmetric] \n"
                " *         -L Log-dist option (to compute the logarithm of the cost matrix) [default=0]\n"
                " *         -m Maximum distance for binning (default= 20000) [in meters]\n"
                " *     -S Sorensen option: Compute sorensen with given original network as benchmark (int, default=0) \n"
                " *         (warning: This option slows down the simulation)\n"
                " *         -Q Path to original network to compare to, if wanted in edge list format: node_origin node_dest t_ij (int int int)\n"
                " *     -E Entropy option: Compute surprise histogram for the various repetitions of the networks (int, default=0), only for more than 10 reps \n"
                " *         (warning: This option slows down the simulation)\n"
                " *     -I Indistiguishable option (only for case of fixed degrees or fixed strengths see README) (>0 for indistinguishable weights) (Default =0) \n"
                " *     -A Aggregated binary option (only for case of fixed degrees, strenghts or both, see README) (>0 for aggregated binary weights) (Default =0) \n"
                " *         -M Number of layers from aggregation, if -I or -A >0 (Default=1) \n"
                " *     -l Self-loop option (>0 for accepting them) (Default =1) \n"
                " *     -c Clustering option (1 for yes) (warning: Depending on av_s makes simulations orders of magnitude slower) (Default=0)\n"
                " *     -v Verbose (1 for on, 0 for off) (Default 0)\n"
                " *     -r Number of reps for averaging (default=100)\n"
                " *     -p Printing option of sample adjacency list (default=0, >0 for yes)\n"
                " *     -s seed.initial seed for random generator (int) (default=1)\n"
                " *     -x Exponent for log-binning (-1 for no log binning) (Default=-1)\n"
                " *     -h Number of header lines in input files (default=1)\n"
                " *     -R Number of maximum trials reps for binary linking (default=100) [only for binary cases, see README]\n"
                "Please, read the README file for more info or the DOCS folder for deatils on simulations!\n");
        return 0;
    }

    /****** OPtion checking (only basics) ******/
    printf("## Checking options...\n");
    if((meth<0)&&(cases==0))
    {
        printf("\tSelect apropriate method flag -e, aborting...\n");
        abort();
    }
    if(gamma<0)
    {
        if((cases==2)||(cases==3))
        {
            printf("\tSelect appropriate gamma value, aborting...\n");
            abort();
        }
    }
    if((opt_dir!=1)&&(opt_dir!=0))
    {
        printf("\tSelect directed or undirected!, aborting...\n");
        abort();
    }
    if(opt_clust==1)
    {
        if(opt_dir>0)
        {
            printf("\tIgnoring clustering option, only defined for undirected networks ....\n");
            opt_clust=0;
        }else{
            printf("\tClustering option selected! This may cause low performance for high <s>! ....\n");
        }
    }   
    if(opt_log>0)
    {
        dmax = log(dmax);
    }
    if((layers>1)&&(opt_agg<=0)&&(opt_indist<=0))
    {
        printf("\t Ignoring number of layers for ME case, they are all equal! ....\n");
        layers=1;
    }
    if((opt_indist>0)&&(opt_agg>0))
    {
        printf("\t Agg option refers to binary aggregated formats and is incompatible with the Indisinguishable option... Chose one or the other! \n");
        abort();
    }
    if(bin_exp<=1) bin_exp=1.05;

    /*** Set rand generator (uses GLS THAU) ***/ 
    if(seed==1)
    {
        seed = seed + time(NULL); // Change for trully random generation
        if(opt_verbose>1) printf("Initial seed not given, using processor time...\n");
    }
    gsl_rng * randgsl = gsl_rng_alloc(gsl_rng_taus);    /// we initialize the random generator of the gsl
    gsl_rng_set(randgsl,seed);


/****** Common variables definition ******/
    char cadena[100]; // to print files
    double** dist; // distance matrix
    double** ps; // probabilities (pij fixed case)
    double * pij_flat; // probs (canonical case)
    double** pij;
    double av_k; // average node strength (for naming purposes)
    double norm; // normalization
    int ps_opt; // radiation model option (else ignored)
    // for entropy models 
    double ** x2;
    // for others
    int** xx2;
    int* xx;
    double rho,E_av;    // for fixed degree
    double xmax,ymax,zmax,wmax; // max values for lagrange multipliers
    double w_max; // max value to bin histogram
    double LogL; // loglikelyhood
    LogL = -1;
    /****** Distance reading ******/
    if(opt_dist>0) dist = read_distances(file_d,N_nodes,header,opt_log, opt_verbose);
/*    int i,j;
    if(fabs(gamma)<1e-12)
    {
        dist = cast_mat_double(N_nodes,N_nodes);
        if(opt_log>0)
        {
            for(i=0;i<N_nodes;i++)
            {
                for(j=0;j<N_nodes;j++)
                {
                    dist[i][j]=1;
                }
            }
        }else{
            for(i=0;i<N_nodes;i++)
            {
                for(j=0;j<N_nodes;j++)
                {
                    dist[i][j]=0;
                }
            }
            
        }
    }
*/

 /***********************************************************************
    Cases selection + distance reading + allocation
 ************************************************************************/ 


    printf("## Case selection\n");
    if(cases>5)
    {
        printf("\tSelect appropriate case (0 to 5). Aborting...\n");
        abort();
    }else{
        if(cases==0) // other models case
        {
            ps=cast_mat_double(N_nodes,N_nodes);
            printf("\t0:: Selected Case with non entropy models\n");
            if(opt_entropy>0) printf("\t[Entropy/surprise calculation is not implemented]\n");
            if(opt_dist<=0)
            {
                printf("\t\t Must provide cost matrix! Aborting...\n");
                abort();
            }           
            if(meth>1)
            {
                printf("\t\tSeq. Gravity model selected \n");
                if(gamma<0)
                {
                    double* surf = nearest_neighbor_dist(dist,N_nodes); // average
                    double av_surf = sum_vec_double(surf,N_nodes)*1./(double)N_nodes; // average inter-nodal nearest distance
                    printf("\t\t Gamma not specified, using default from paper:\n \t Surface calculated using as proxy for the average surface the squared " 
                    "average distance between nearest neighbors: %f\n",av_surf); // considering disc
                    gamma = compute_gamma_frenchies((av_surf*1e-3/2.)*(av_surf*1e-3/2.)*3.141593); // surface in meters needs to be transformed in km^2
                    free(surf);
                }
                printf("\t\t Selected inverse gamma: %f\n",1./gamma);
            }else{
                printf("Ignoring gamma value, Radiation model selected\n");
            }
            if(opt_dir>0)
            {
                xx2 = read_node_list_int(file_s, N_nodes, header, opt_verbose); // strenght sequence (ints)
                T= sum_vec_int(xx2[0],N_nodes); // T is \sum_i s_i (for all cases)
                xmax = (double)max_value_int(xx2[0],N_nodes);
                ymax = (double)max_value_int(xx2[1],N_nodes);
                w_max=10.*xmax*ymax/T;
            }else{
                xx = read_node_list_int_undir(file_s, N_nodes, header, opt_verbose); // strenght sequence (ints)
                T= sum_vec_int(xx,N_nodes); // T is \sum_i s_i (for all cases)
                xmax = (double)max_value_int(xx,N_nodes);
                w_max= 10*xmax*xmax/T;
            }       
        }
        else if(cases==1) // fixed pij case
        {
            printf("\t1:: Selected Case with provided custom pij matrix\n");
            pij = read_net_list_double(file_s, N_nodes, header, opt_verbose);
            norm = sum_matrix_double(pij, N_nodes, N_nodes);
            if((T<(double)N_nodes/100)||(T<10))
            {
                printf("\t\t Specify enough events with flag T. Current:%d! Aborting...\n",(int)T);
                abort();
            }
            scale_matrix(pij, N_nodes, N_nodes, 1./norm);
            w_max= (T*matrix_max_value_double(pij, N_nodes,N_nodes));
        }else{ // all other entropy cases
            // number of attributes to read //
            if(cases==2)
            {
                if((opt_indist<=0)&&(opt_agg<=0)&&(fabs(gamma)>1e-12))
                {
                    printf("\t2:: Selected Case with fixed strength sequence and average cost\n");          
                }else{
                    printf("\t2:: Selected Case with fixed strength sequence (gamma=0)\n");
                }
                if(opt_dist<=0) // not dist provided
                {
                    if((opt_indist<=0)&&(opt_agg<=0)&&(fabs(gamma)>1e-12)) // ME and gamma provided
                    {
                        printf("\t\t Must provide cost matrix! Aborting...\n");
                        abort();
                    }else{ // all other cases
                        dist = cast_mat_double(N_nodes,N_nodes); 
                    }
                }
                // read inputs //
                if(opt_dir>0)
                {
                    x2 = read_node_list_xatts_double(file_s, N_nodes, 2, header, opt_verbose);
                    if(opt_self>0)
                    {
						xmax = max_value_double(x2[0],N_nodes);
						ymax = max_value_double(x2[1],N_nodes);
					}else{
						xmax = max_value_double_k(x2[0],N_nodes,2);
						ymax = max_value_double_k(x2[1],N_nodes,2);
					}
                    if((opt_indist<=0)&&(opt_agg<=0))
                    { // ME case
                        w_max= xmax*ymax;

                    }else{
                        if(opt_indist>0)
                        { // weighted
							//printf("xm :%f ym %f layers %d\n",xmax,ymax,layers);
							assert(xmax*ymax<1);
							w_max=layers*xmax*ymax/(1.-xmax*ymax);                           
                        }else{ // binary
                            w_max=layers*xmax*ymax/(1.+xmax*ymax);
                        }
                    }
                }else{
                    x2 = read_node_list_xatts_double(file_s, N_nodes, 1, header, opt_verbose);
					xmax = max_value_double(x2[0],N_nodes);
                    if(opt_self>0)
                    {
						ymax = xmax;
					}else{
						ymax = max_value_double_k(x2[0],N_nodes,2);
					}
                    if((opt_indist<=0)&&(opt_agg<=0))
                    { // ME case
                        w_max= xmax*ymax;
                    }else{
                        if(opt_indist>0)
                        { // weighted
							assert(xmax*ymax<1);
                            w_max=layers*xmax*ymax/(1.-xmax*ymax); 
                        }else{ // binary
                            w_max=layers*ymax*xmax/(1.+ymax*xmax);
                        }
                    }
                }
                T = w_graph_compute_T_analytic_from_xy(x2, N_nodes, layers, 0, opt_dir, opt_self, opt_indist, opt_agg, 0);
            }else if(cases==3){
                meth = 2;
                printf("\t3:: Selected Case with fixed strength sequence and total number of binary edges\n");
                if(opt_dir>0)
                {
                    x2 = read_node_list_xatts_double(file_s, N_nodes, 2, header, opt_verbose);
                    xmax = max_value_double(x2[0],N_nodes);
                    ymax = max_value_double(x2[1],N_nodes);
                    w_max=xmax*ymax*exp(xmax*ymax)*gamma/(1+gamma*(exp(xmax*ymax)-1));
                }else{
                    x2 = read_node_list_xatts_double(file_s, N_nodes, 1, header, opt_verbose);
                    xmax = max_value_double(x2[0],N_nodes);
                    w_max=xmax*xmax*exp(xmax*xmax)*gamma/(1+gamma*(exp(xmax*xmax)-1));
                }
                if((gamma<0)&&(opt_indist<=0)&&(opt_agg<=0)) // if ME and not gamma
                {
                    printf("\t\t Specify appropriate (positive) gamma value! Aborting...\n");
                    abort();
                }
                T = w_graph_compute_T_analytic_from_xy(x2, N_nodes, layers, gamma, opt_dir, opt_self, opt_indist, opt_agg, 1);
            }else if(cases==4){
                meth = 2;
                printf("\t4:: Selected Case with fixed strength sequence and degree sequence\n");                       
                // number of attributes to read: Exception (2 sets of lagrange multipliers) //
                if(opt_dir>0)
                {
                    x2 = read_node_list_xatts_double(file_s, N_nodes, 4, header, opt_verbose);
                    xmax = max_value_double(x2[0],N_nodes);
                    ymax = max_value_double(x2[1],N_nodes);
                    zmax = max_value_double(x2[2],N_nodes);
                    wmax = max_value_double(x2[3],N_nodes);
                    w_max=xmax*ymax*exp(xmax*ymax)*zmax*wmax/(1+zmax*wmax*(exp(xmax*ymax)-1));
                }else{
                    x2 = read_node_list_xatts_double(file_s, N_nodes, 2, header, opt_verbose);
                    xmax = max_value_double(x2[0],N_nodes);
                    zmax = max_value_double(x2[1],N_nodes);
                    w_max=xmax*xmax*exp(xmax*xmax)*zmax*zmax/(1+zmax*zmax*(exp(xmax*xmax)-1));
                }
                T = w_graph_compute_T_analytic_from_xy(x2, N_nodes, layers, 0, opt_dir, opt_self, opt_indist, opt_agg, 2);
            }else{
                meth = 2;
                printf("\t2:: Selected Case with fixed degree sequence\n");                     
                printf("\t\t Warning: The surprise calculation will only include the terms corresponding to the binary structure \n");                       
                if(opt_dir>0)
                {
                    x2 = read_node_list_xatts_double(file_s, N_nodes, 2, header, opt_verbose); // degree lagrange sequence (double)
                    xmax = max_value_double(x2[0],N_nodes);
                    ymax = max_value_double(x2[1],N_nodes);
                    E_av = w_graph_compute_E_binary_directed(x2, N_nodes, opt_self);    // compute expected degree, stop if E>T:
                }else{
                    x2 = read_node_list_xatts_double(file_s, N_nodes, 1, header, opt_verbose); // sdegree lagrange sequence sequence (double)
                    xmax = max_value_double(x2[0],N_nodes);
                    E_av = w_graph_compute_E_binary_undirected(x2[0], N_nodes, opt_self); // compute expected degree, stop if E>T:
                }
                if(T<1)
                {
                    printf("\t\t Specify enough events with flag T! Aborting...\n");
                    abort();
                }
                if(layers==1) // only simplex
                {
                    if(opt_agg>1)
                    {
                        rho = 1.; // if binary
                    }else{
                        rho = w_graph_compute_rho(E_av,(int)T,opt_indist); // compute lagrange multiplier (analytical)
                    }
                }else{
                    if(gamma>0)
                    {
                        rho = gamma;
                    }else{
                        printf("\t\t If agg or indist and layers>1, need to provide lagrange mult. for number of edges. Aborting...\n");
                        abort();
                    }
                }
                printf("\t\t Total expected edges: %f, Total sampling:%f, average edge weight: %f\n",E_av,T,T/E_av);
                w_max = T/E_av ;
            }
        }
    }
    /****** Info printing ******/
    printf("## Information about simulation\n");
    if(opt_dir==0)
    {
        printf("\t Undirected simulation\n");
    }else{
        printf("\t Directed simulation\n");
    }
    if(opt_self>0)
    {
        printf("\t Self-loops accepted... \n");
    }else{
        printf("\t Self-loops prohibited... \n");       
    }
    if(opt_agg>0)
    {
        printf("\t Generating case Aggregated Binary (AG) with %d layers\n",layers);                
    }else{
            if(opt_indist>0)
            {
                printf("\t Generating case Aggregated Weighted (AW) with %d layers\n",layers);              
            }else{
                printf("\t Generating case Multi-Edge (ME)\n");             
            }
    }
    if(opt_dist>0)
    {
        printf("\t Analizing distance attributes \n");
    }
    if(opt_soren>0)
    {
        printf("\t Computing sorensen index \n");   
    }   
    if((opt_entropy>0)&&(reps>10))
    {
        printf("\t Computing surprise histogram \n");   
    }else{
		opt_entropy=0;   
	}
    if(opt_print_tr>0)
    {
        printf("\t Printing sample of adj. matrix for first iteration\n");              
    }   
    if((cases<2)||(cases==5))
    {
        printf("\t Total trips to sample :%f\n",T);
    }
    printf("\t SEED=%i\n",seed);
    printf("\t Average over %d reps\n", reps);

    av_k=T/(double)N_nodes; // used in some cases, in others not


/***********************************************************************
    Preparing ensemble reps!
 ************************************************************************/  
    w_max=w_max*10; // x10 to be sure
    if(w_max<1)
    {
        w_max=10;
    }
    if(gsl_finite(w_max)!=1)
    {
        w_max=1e7; //just in case		
	}
    printf("\t Wmax for histogram: %d\n",(int)w_max);
    // length of containers //
    int len_acc_nodes;
    if(opt_dir>0)
    {
        len_acc_nodes=12; // in-out variables

    }else{
        if(opt_clust==1)
        {
            len_acc_nodes=8;
        }else{
            len_acc_nodes=6;
        }
    }
    // setting up node containers //
    gsl_histogram ** acc_ensemble;
    if(opt_dist>0)
    {
        acc_ensemble = w_graph_dist_all_stats_ensemble_allocate(opt_dir, dmax, (int)w_max);
    }else{
        acc_ensemble= w_graph_all_stats_ensemble_allocate(opt_dir, 1, (int)w_max ,1, N_nodes ,(int)w_max);        
    }
    double* entropy_seq;
    if(opt_entropy>0) entropy_seq = cast_vec_double(reps); // entropy histogram   
    double ** node_cont; //node info container for averaging
    double ** node_cont2; // node info std for averaging
    int ** node_nonzero; // node count non-zero strength/degree
    double * Tcont; // Total T container (av and std)
    node_cont=cast_mat_double(N_nodes,len_acc_nodes);
    node_cont2=cast_mat_double(N_nodes,len_acc_nodes);
    node_nonzero=cast_mat_int(N_nodes,2); // to count the number of times a node is non-connected (in-out)  
    Tcont=cast_vec_double(2); // mean and std deviation
    int L=0; // dummy length
    
    if(cases==1)
    {
        // scale pij to actual values if poisson!
        if(meth==2)
        {
            scale_matrix(pij,N_nodes,N_nodes,(double)T);
        }else{ // flatten matrix (normalized)
            pij_flat=flatten_matrix_double(pij, N_nodes, N_nodes, opt_self, &L);
            free_mat_double(pij,N_nodes);
            // if canonical, flatten pij    
        }
    }
    if(cases==2)
    {
        /** Probabilities choice (could be extended to other types) **/
        if(opt_dir>0)
        {
            ps = cast_mat_double(N_nodes,1);
            if (meth!=2) 
            {
                ps[0] = prob_mult_Cs_dir(N_nodes, x2, gamma, dist);
                //ps = prob_mult_s_dir(x2, N_nodes, opt_self); // needs gamma input!
            }
        }else{
            if (meth!=2)
            {
                //ps = prob_mult_s_undir(x, N_nodes, opt_self); // needs gamma input!
                ps[0] =  prob_mult_Cs_undir(N_nodes, x2[0], gamma, dist);
            }
        }
    }
/*********** Read original file if given ******************************/    
    double soren = 0;
    double soren2 = 0;
    double sorfake;
    W_GRAPH* ori;
    W_GRAPH* WG;
    if(opt_soren>0)
    {
        printf("\t Loading original file to compute Sorensen\n");
        ori = w_graph_read_edge_list(file_ori, N_nodes, opt_dir,header, opt_verbose);
    }
/***********************************************************************
    Start of ensemble reps!
 ************************************************************************/

    int r;

    for(r=0;r<reps;r++)
    {
        if(cases==0) // other models
        {
            if(r==0)
            {
                ps_opt = -1;
            }else{
                ps_opt =1;
            }
            if(meth==2)
            {
                if(opt_verbose==1) printf("============### Seq model multinomial ####===========\n"); fflush(stdout);
                if(opt_dir>0)
                {
                    WG = w_graph_seq_gravity_multinomial_directed(N_nodes, xx2, dist, gamma, randgsl, opt_self, opt_verbose);
                }else{
                    printf("Not defined for undirected case... Aborting\n"); 
                    abort();
                }
            }else if(meth==3){
                if(opt_verbose==1) printf("============### Seq model bernouilli ####===========\n"); fflush(stdout);
                if(opt_dir>0)
                {
                    WG = w_graph_seq_gravity_bernouilli_directed(N_nodes, xx2, dist, gamma, randgsl, opt_self, opt_verbose);
                }else{
                    printf("Not defined for undirected case... \n"); 
                    abort();
                }
            }else if (meth==1){
                if(opt_verbose==1)printf("============### Radiation Stochastic model ####===========\n"); fflush(stdout);
                if(opt_dir>0)
                {
                    WG =   w_graph_radiation_model_stochastic_directed(N_nodes,xx2, dist, randgsl,opt_self,opt_verbose);
                }else{
                    printf("Not implemented for undirected case... \n");
                    abort();
                }
            }else{
                if(opt_verbose==1)printf("============### Radiation Multinomial model ####===========\n"); fflush(stdout);
                if(opt_dir>0)
                {
                    WG =   w_graph_radiation_model_multinomial_directed(N_nodes,xx2, dist, randgsl,opt_self,opt_verbose,ps_opt,ps);
                }else{
                    printf("Not implemented for undirected case... \n");
                    abort();
                }           
            }
            //if(opt_entropy) entropy_seq[r]= w_graph_surprise_multinomial(WG,N_nodes,opt_dir); // assuming multinomial case
        }else if(cases==1){ // custom pij
            if(meth==2)
            {
                if(opt_verbose==1) printf("============### Poisson model ####===========\n"); fflush(stdout);
                if(opt_dir>0)
                {
                    WG = custompij_poisson_directed_graph(pij, N_nodes, randgsl, opt_verbose, opt_self);
                }else{
                    WG = custompij_poisson_undirected_graph(pij, N_nodes, randgsl, opt_verbose, opt_self);
                }
            }else{
                if(meth==0)
                {
                    if(opt_verbose==1)printf("============### Multinomial model ####===========\n"); fflush(stdout);
                    if(opt_dir>0)
                    {
                        WG = multinomial_directed_graph(pij_flat, N_nodes, T, randgsl,opt_verbose, opt_self);
                    }else{
                        WG = multinomial_undirected_graph(pij_flat, N_nodes, T, randgsl,opt_verbose, opt_self);
                    }
                }else{
                    if(opt_verbose==1)printf("============### Poisson Multinomial model ####===========\n"); fflush(stdout);
                    if(opt_dir>0)
                    {
                        WG = poisson_multinomial_directed_graph(pij_flat, N_nodes, T, randgsl, opt_verbose, opt_self);
                    }else{
                        WG = poisson_multinomial_undirected_graph(pij_flat, N_nodes, T, randgsl, opt_verbose, opt_self);
                    }
                }
            if(opt_entropy) entropy_seq[r]= w_graph_surprise_poisson_pij(WG,N_nodes,pij,opt_self, opt_dir); // assuming poisson case, multinomial not implemented
            }
        }else if(cases==2){ // fixed s or fixed s and C
            if(meth==2)
            {
                if((opt_indist<=0)&&(opt_agg<=0)) // ME case
                {
                    if(opt_verbose==1) printf("============### Poisson model ####===========\n"); fflush(stdout);
                    if(fabs(gamma)>1e-12)
                    {
						if(opt_dir>0)
						{
							WG = gravity_poisson_directed_graph2(x2, N_nodes, dist, gamma, randgsl, opt_verbose, opt_self);
						}else{
							WG = gravity_poisson_undirected_graph2(x2[0], N_nodes, dist, gamma, randgsl, opt_verbose, opt_self);
						}
						if(opt_entropy) entropy_seq[r]= w_graph_surprise_poisson_dist(WG,x2,N_nodes,dist,gamma,opt_self,opt_dir); // assuming poisson case
					}else{
						if(opt_dir>0)
						{
							WG = fixeds_poisson_directed_graph2(x2, N_nodes , randgsl, opt_verbose, opt_self);
						}else{
							WG = fixeds_poisson_undirected_graph2(x2[0], N_nodes, randgsl, opt_verbose, opt_self);
						}
						if(opt_entropy) entropy_seq[r]= w_graph_surprise_poisson(WG,x2,N_nodes,opt_self,opt_dir); // assuming poisson case
					}
                }else{
                    if(layers==1) // single layer
                    {
                        if(opt_indist>0) // weighted case (geometric)
                        {
                            if(opt_verbose==1) printf("============### Geometric model ####===========\n"); fflush(stdout);
                            if(opt_dir>0)
                            {
                                //WG = fixeds_geometric_directed_graph2(x2, N_nodes , randgsl, opt_verbose, opt_self);
                                WG = fixeds_negbinomial_directed_graph2(x2, N_nodes, 1,  randgsl, opt_verbose, opt_self);
                            }else{
                                //WG = fixeds_geometric_undirected_graph2(x2[0],  N_nodes , randgsl, opt_verbose, opt_self);
                                WG = fixeds_negbinomial_undirected_graph2(x2[0], N_nodes, layers , randgsl, opt_verbose, opt_self);
                            }
                            if(opt_entropy) entropy_seq[r]= w_graph_surprise_geometric(WG,x2,N_nodes,opt_self,opt_dir); // assuming poisson case
                        }else{ // binary
                            printf("Cannot fix the strength sequence of a binary network, use case 5 instead. Aborting... \n");
                            abort();
                        }
                    }else{ // multilayer
                        if(opt_indist>0) // neg.binomial
                        {
                            if(opt_verbose==1) printf("============### Negative Binomial model ####===========\n"); fflush(stdout);
                            if(opt_dir>0)
                            {
                                WG = fixeds_negbinomial_directed_graph2(x2, N_nodes, layers,  randgsl, opt_verbose, opt_self);
                            }else{
                                WG = fixeds_negbinomial_undirected_graph2(x2[0], N_nodes, layers , randgsl, opt_verbose, opt_self);
                            }
                            if(opt_entropy) entropy_seq[r]= w_graph_surprise_negbinomial(WG,x2,N_nodes,layers,opt_self,opt_dir); // assuming poisson case
                        }else{ // binomial
                            if(opt_verbose==1) printf("============### Binomial model ####===========\n"); fflush(stdout);
                            if(opt_dir>0)
                            {
                                WG = fixeds_binomial_directed_graph2(x2, N_nodes, layers, randgsl, opt_verbose, opt_self);
                            }else{
                                WG = fixeds_binomial_undirected_graph2(x2[0], N_nodes, layers, randgsl, opt_verbose, opt_self);
                            }
                            if(opt_entropy) entropy_seq[r]= w_graph_surprise_binomial(WG,x2,N_nodes,layers,opt_self,opt_dir); // assuming poisson case
                        }
                    }
                }
            }else{
                if((opt_indist<=0)&&(opt_agg<=0))
                {
                    if(meth==0)
                    {
                        if(opt_verbose==1)printf("============### Multinomial model ####===========\n"); fflush(stdout);
                        if(opt_dir>0)
                        {
                            WG = multinomial_directed_graph(ps[0], N_nodes, T, randgsl,opt_verbose, opt_self);
                        }else{
                            WG = multinomial_undirected_graph(ps[0], N_nodes, T, randgsl,opt_verbose, opt_self);
                        }
                    }else{
                        if(opt_verbose==1)printf("============### Poisson Multinomial model ####===========\n"); fflush(stdout);
                        if(opt_dir>0)
                        {
                            WG = poisson_multinomial_directed_graph(ps[0], N_nodes, T, randgsl, opt_verbose, opt_self);
                        }else{
                            WG = poisson_multinomial_undirected_graph(ps[0], N_nodes, T, randgsl, opt_verbose, opt_self);
                        }
                    }
                    //if(opt_entropy) entropy_seq[r]= w_graph_surprise_multinomial(WG,N_nodes,opt_dir); // assuming poisson case
                }else{
                    printf("If indist or agg set to True, only available method is Grand-Canonical. Aborting... \n");
                    abort();
                }
            }
        }else if(cases==3){ // fixed s and E
            if(opt_verbose==1) printf("============### ZIP model ####===========\n"); fflush(stdout);
            if(opt_dir>0)
            {
                WG = fixedEs_poisson_directed_graph(x2, N_nodes, gamma, randgsl, opt_verbose, opt_self, max_reps);
            }else{
                WG = fixedEs_poisson_undirected_graph(x2[0], N_nodes, gamma, randgsl, opt_verbose, opt_self, max_reps);
            }
            if(opt_entropy) entropy_seq[r] = w_graph_surprise_ZIP(WG, x2, N_nodes, gamma, opt_self, opt_dir);
        }else if(cases==4){ // fixed s and k
            if(opt_agg<=0)
            {
                if(opt_verbose==1) printf("============### ZIP model ####===========\n"); fflush(stdout);
                if(opt_dir>0)
                {
                    WG = fixedks_poisson_directed_graph(x2, N_nodes, randgsl, opt_verbose, opt_self, max_reps);
                }else{
                    WG = fixedks_poisson_undirected_graph(x2, N_nodes, randgsl, opt_verbose, opt_self, max_reps);
                }
                if(opt_entropy) entropy_seq[r] = w_graph_surprise_ZIP2(WG, x2, N_nodes, opt_self, opt_dir);
            }else{
                if(layers>1)
                {
                    if(opt_verbose==1) printf("============### ZIB model ####===========\n"); fflush(stdout);
                    if(opt_dir>0)
                    {
                        WG = fixedks_binomial_directed_graph(x2, N_nodes, layers, randgsl, opt_verbose, opt_self, max_reps);
                    }else{
                        WG = fixedks_binomial_undirected_graph(x2, N_nodes, layers, randgsl, opt_verbose, opt_self, max_reps);
                    }
                    if(opt_entropy) entropy_seq[r] = w_graph_surprise_ZIB2(WG, x2, N_nodes, layers, opt_self, opt_dir);                  
                }else{
                    printf("Cannot fix the strength and degree sequence of a single binary network, use case 5 instead. Aborting... \n");
                    abort();
                }
            }
        }else{ // fixed k
            if((opt_indist<=0)&&(opt_agg<=0)) // ME case
            {
                if(opt_verbose==1) printf("============### ZIP Uniform model ####===========\n"); fflush(stdout);
                if(opt_dir>0)
                {
                    WG = fixedk_poisson_directed_graph(x2, rho, N_nodes, randgsl, opt_verbose, opt_self, max_reps);                    
                }else{
                    WG = fixedk_poisson_undirected_graph(x2[0], rho, N_nodes, randgsl, opt_verbose, opt_self, max_reps);                
                }
            }else{ // others
                if(layers==1) // simplex
                {
                    if(opt_indist>0) // weighted
                    {
                        if(opt_verbose==1) printf("============### ZIG Uniform model ####===========\n"); fflush(stdout);
                        if(opt_dir>0)
                        {
                            WG = fixedk_geometric_directed_graph(x2, rho, N_nodes, randgsl, opt_verbose, opt_self, max_reps);

                        }else{
                            WG = fixedk_geometric_undirected_graph(x2[0], rho, N_nodes, randgsl, opt_verbose, opt_self, max_reps);                
                        }
                    }
                    if(opt_agg>0) // binary bernouilli
                    {
                        if(opt_verbose==1) printf("============### Bernouilli model ####===========\n"); fflush(stdout);
                        if(opt_dir>0)
                        {
                            WG = fixedk_bernouilli_directed_graph(x2, N_nodes, randgsl, opt_verbose, opt_self);                    
                        }else{
                            WG = fixedk_bernouilli_undirected_graph(x2[0], N_nodes, randgsl, opt_verbose, opt_self);                
                        }                       
                    }
                }else{
                    if(opt_indist>0)
                    {
                        if(opt_verbose==1) printf("============### ZINB Uniform model ####===========\n"); fflush(stdout);
                        if(opt_dir>0)
                        {
                            WG = fixedk_negbinom_directed_graph(x2, rho, N_nodes, layers, randgsl, opt_verbose, opt_self, max_reps);

                        }else{
                            WG = fixedk_negbinom_undirected_graph(x2[0], rho, N_nodes, layers, randgsl, opt_verbose, opt_self, max_reps);                
                        }                       
                    }
                    if(opt_agg>0)
                    {
                        if(opt_verbose==1) printf("============### ZIG uniform model ####===========\n"); fflush(stdout);
                        if(opt_dir>0)
                        {
                            WG = fixedk_binom_directed_graph(x2, rho, N_nodes, layers, randgsl, opt_verbose, opt_self, max_reps);                    
                        }else{
                            WG = fixedk_binom_undirected_graph(x2[0], rho, N_nodes, layers, randgsl, opt_verbose, opt_self, max_reps);                
                        }                           
                    }                   
                }
            }
            if(opt_entropy) entropy_seq[r] = w_graph_surprise_bernouilli(WG, x2, N_nodes, opt_self, opt_dir); // not counting additional constant terms                  
        }

        if(r==0) // if first rep, store all stats
        {
            /// extended stats with distance ////
            if (opt_verbose>0) printf("... Graph stats ... \n");
            if(opt_dist>0)
            {
                w_graph_dist_all_stats(WG, N_nodes, 0,  bin_exp, opt_dir,opt_self,dist,randgsl,dmax,opt_verbose);
            }else{
				w_graph_all_stats(WG, N_nodes, r, bin_exp, opt_dir, opt_self, -1, opt_verbose);
            }
            if (opt_verbose>0) printf("... Node stats ... \n");
            w_graph_node_stats_list(WG,N_nodes,0, opt_dir, opt_clust, opt_self, opt_verbose);
            if(opt_print_tr==1)
            {
                if(opt_verbose>0)printf("Printing adj matrix\n");
                sprintf(cadena,"N%d_cust.tr",N_nodes);
                w_graph_print_adj_list(WG, N_nodes, cadena);
            }
        }
        if (opt_verbose>0) printf("... Ensemble stats ... \n");
        /// extended stats with distance ////
        w_graph_node_stats_ensemble(WG,N_nodes,node_cont,node_cont2,node_nonzero, Tcont,opt_dir,opt_clust);
        if(opt_dist>0)
        {
            w_graph_dist_all_stats_ensemble_update(acc_ensemble,WG, N_nodes, opt_dir,dist);
        }else{
            w_graph_all_stats_ensemble_update(acc_ensemble, WG, N_nodes, opt_dir);
        }
        if(opt_soren>0)
        {
            sorfake = w_graph_compute_sorensen(WG, ori, N_nodes);
            soren += sorfake;
            soren2 += sorfake*sorfake;
        }
        w_graph_free_nodes(WG, N_nodes);
    }
    if(opt_soren>0)
    {
		sorfake = -1;
        // if cases= 1,2
        if((cases==1)||(cases==2))
        {
			if(cases==1) // case 1 all poisson from pij
			{
				if(meth==2)
				{
					LogL = w_graph_loglikelyhood_poisson(ori, N_nodes, pij, opt_self, opt_dir);
					sorfake = w_graph_compute_sorensen_av(ori, pij, N_nodes, T);
				}
			}else{ // case 2
				if((opt_indist<=0)&&(opt_agg<=0)) // ME case
				{
					LogL = w_graph_loglikelyhood_poisson_xy(ori, x2, N_nodes, opt_self, opt_dir);
				}else{
					if(opt_indist>0)
					{
						if(layers>1)
						{
							LogL = w_graph_loglikelyhood_negbinomial_xy(ori, x2, N_nodes, layers, opt_self, opt_dir);							
						}else{
							LogL = w_graph_loglikelyhood_geometric_xy(ori, x2, N_nodes, opt_self, opt_dir);														
						}
					}else{
							LogL = w_graph_loglikelyhood_binomial_xy(ori, x2, N_nodes, layers, opt_self, opt_dir);	
					}
				}
			}
            
            printf("\t Minus Loglikelyhood of original file, given the model:%.1f \n",LogL); // for the moment only implemented here              
        }
    }
    if(cases==0)free_mat_double(ps,N_nodes);        
    if(cases==1)
    {
        if(meth==2)
        {
            free(pij);      
        }else{
            free(pij_flat);
        }
    }
    if(cases==2)
    {
        if(meth<2) free(ps[0]);
    }

/***********************************************************************
     Print output
************************************************************************/   
    int len=6;
    /// extended stats with distance ////
    //fflush(stdout);
    if (opt_verbose>0) printf("... Printing ensemble stats ... \n");
    if(opt_dist>0)
    {
        w_graph_dist_all_stats_ensemble_print(acc_ensemble, len, reps, N_nodes,opt_dir);
    }else{
        w_graph_all_stats_ensemble_print(acc_ensemble, 2, reps, N_nodes, opt_dir);
    }
    w_graph_node_stats_ensemble_print(reps, N_nodes, Tcont, node_cont, node_cont2, node_nonzero, bin_exp,len_acc_nodes, opt_dir);
    if(opt_entropy)
    {
        sprintf(cadena,"N%davs%8.5fentropies.hist",N_nodes,av_k);
        w_graph_print_entropy(entropy_seq,reps,cadena);
    }
    if (opt_soren>0)
    {
        soren=soren/(double)reps;
        printf("... Average Sorensen:%f+-%f ...\n",soren,sqrt(soren2/(double)reps - (soren*soren)));
        if(sorfake>0) printf("... Sorensen of average :%f\n",sorfake);
        w_graph_free_nodes(ori,N_nodes);
    }
    gsl_rng_free (randgsl);
    free(entropy_seq);
    return 0;
}

