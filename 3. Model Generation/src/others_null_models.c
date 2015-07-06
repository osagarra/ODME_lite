/************************************************************
 *
 *                    Transport Models
 *
 *	Null models for weighted transport nets
 *
 *************************************************************/

#include "others_null_models.h"

/*********************************************/
/************** Lenormand et altr. (PLOS ONE) ********************/
/*
*    Lenormand M, Huet S, Gargiulo F, Deffuant G (2012) A Universal Model of
*    Commuting Networks. PLoS ONE 7(10): e45985. doi:10.1371/journal.pone.0045985
*/
/*********************************************/

double compute_gamma_frenchies(double surface){
	//double alpha=0.02731 ; //
	double alpha=3.15e-4;// (surface in km)
    double nu=0.1777;
	double res;
	res=alpha*pow(surface,-nu);
	return res; // in m^-1
}
/**********************************/
double compute_p_french_bernouilli(int* s_in, double **d, double gamma, int origin, int dest, int* ind_out, int* ind_in ,int len, int self_opt){
	int i,in_d,out_d,real_dest,real_origin,aux_dest;
	double norm,pp;
	// get origin and all to zero //
	real_origin = ind_out[origin];
	real_dest   = ind_in [dest];
	out_d=maxeq_int(real_dest,real_origin);
	in_d=mineq_int(real_dest,real_origin);
	if((real_dest!=real_origin) || (self_opt>0)) // avoid self-loops (or not)	
	{
		pp=(double)s_in[dest]*exp(-d[out_d][in_d]*gamma);
	}else{
		pp=0;
	}
	norm=0;

	for(i=0;i<len;i++) // for avaiable ins, compute norm
	{
		aux_dest=ind_in[i];
		if((aux_dest!=real_origin) || (self_opt>0)) // avoid self-loops (or not)
		{
			out_d=maxeq_int(aux_dest,real_origin);
			in_d=mineq_int(aux_dest,real_origin);
			norm +=(double)s_in[i]*exp(-d[out_d][in_d]*gamma);
		}
	}
	//printf("Norm: %f| ori: %d %d\n",norm,real_origin,origin);				
	pp = pp/norm;
	return pp;
}

void* compute_p_french_multinomial(double* p, int* s_in, double **d, double gamma, int real_origin, int* ind_in ,int len, int N_nodes, int self_opt, int verbose){
	int i,in_d,out_d,real_dest;
	//double norm;
	// get origin and all to zero //
	//norm=0;
	for(i=0;i<N_nodes;i++)
	{
		p[i]=0;
	}
	for(i=0;i<len;i++)
	{
		real_dest=ind_in[i];
		if((real_dest!=real_origin) || (self_opt>0)) // avoid self-loops (or not)
		{
			out_d=maxeq_int(real_dest,real_origin);
			in_d=mineq_int(real_dest,real_origin);
			p[i]=(double)s_in[i]*exp(-d[out_d][in_d]*gamma);
			//norm += p[i];
			//printf("p:%f, dgamma:%f, s:%d\n",p[i],d[out_d][in_d]*gamma,s_in[i]);fflush(stdout);
		}else{
			p[i]=0;
		}
	}
	return p;
}

W_GRAPH* w_graph_seq_gravity_bernouilli_directed(int N_nodes,int **s, double **d, double gamma, gsl_rng* randgsl, int self_opt, int verbose){
    if(verbose>0) printf("-- Warning: This takes a while, if you are in a hurry, use the multinomial version (asymptotically equivalent)...-- \n");
    // alloc w_graph
    W_GRAPH* WG  = w_graph_alloc(N_nodes);
    WG->opt_self = self_opt;
    WG->opt_dir = 1;
    int i;
    int out_av,in_av;

    // allocs (aux_vectors)
	int* s_in=cast_vec_int(N_nodes);
	int* s_out=cast_vec_int(N_nodes);
	int* in_inds=cast_vec_int(N_nodes);
	int* out_inds=cast_vec_int(N_nodes);
    
    // aux variables
	int tot_trips2,tot_trips3,aux1,aux2;
    aux1=in_av=0;
	aux2=out_av=0;
	tot_trips2=tot_trips3=0;
	
    // check all good
    for(i=0;i<N_nodes;i++)
	{
		if(s[1][i]!=0)
		{
			s_in[aux1]=s[1][i];
			in_inds[aux1]=i;
			aux1++;
			tot_trips2+=s[1][i];
		}else{
			in_inds[N_nodes-in_av-1]=i;
			in_av++;
		}
		if(s[0][i]!=0)
		{
			s_out[aux2]=s[0][i];
			out_inds[aux2]=i;
			aux2++;
			tot_trips3+=s[0][i];
		}else{
			out_av++;
			out_inds[N_nodes-out_av-1]=i;
		}
	}
	//tot_trips2=sum_vec_int(s_outf,N_nodes);
	//tot_trips3=sum_vec_int(s_inf,N_nodes);
	for(i=0;i<N_nodes;i++)
	{
		aux1=count_appearances_int(in_inds,i,N_nodes);
		aux2=count_appearances_int(out_inds,i,N_nodes);
		if((aux1!=1) || (aux2!=1)) // check 1 to 1 mapping
		{
			if(verbose>0)printf("Incorrect mapping for i: %i %i %i\n",i,aux1,aux2);
			abort();
		}
	}
	if(verbose>0)
    {
        if(tot_trips2-tot_trips3 != 0)printf("Warning, total outs (%i) and total ins (%i) different!\n",tot_trips2,tot_trips3);
        printf("Available in %i, Available out %i\n",N_nodes-in_av,N_nodes-out_av);
    }

    //go! //
    int origin,dest,flag,acc,rej;
	double r,p2;
	//double* p = cast_vec_double(N_nodes);
	//double nn;
	in_av=N_nodes-in_av;
	out_av=N_nodes-out_av;

    flag=-1;
    tot_trips2=0;
    do{// While available outs
		rej = 0;
		acc = -1;
        do
        {// while not accepted
			origin  = (int)gsl_rng_uniform_int(randgsl,(unsigned int)out_av); //out from out_av
			aux1=out_inds[origin];
			aux2=aux1;
			//compute_p_french_multinomial(p,s_in, d, gamma, aux1, in_inds, in_av, N_nodes, self_opt, verbose);
			//nn = sum_vec_double(p,in_av);
			// select random origin node
			if(self_opt>0)
			{
				dest  = (int)gsl_rng_uniform_int(randgsl,(unsigned int)in_av); //out from out_av
				aux2=in_inds[dest];
			}else{
				do
				{
					dest  = (int)gsl_rng_uniform_int(randgsl,(unsigned int)in_av); //out from out_av
					aux2=in_inds[dest];
				}while(aux1==aux2);
			}
			if(in_av>1) // if available ins
			{
				p2 = compute_p_french_bernouilli(s_in, d, gamma, origin, dest, out_inds, in_inds ,in_av, self_opt);				
				r=gsl_rng_uniform(randgsl); // rand trial
				//if(r<=p[dest]/nn)
				if(r<=p2)
				{
					//assert((p2-p[dest]/nn)<1e-8);
					acc=1;
					//if(rej>1000)printf("Rejs: %d %f\n",rej,p2);
					//if(rej>1000)printf("Rejs: %d %f\n",rej,p[dest]/nn);
				}else{
					rej++;
				}
			}else{
				acc=1; // always accept if only one left
			}
		}while(acc==-1);
		// add multi-link //
		//printf(" ### %d %d\n",aux1,aux2);
		w_graph_add_multi_link(WG, N_nodes, aux1, aux2, 1);
		tot_trips2++;
		s_in[dest]--;
		s_out[origin]--;
		//printf("Accepted trip! From %i to %i, s_out %i s_in %i To go %i\n",aux1,aux2,s_out[origin],s_in[dest],tot_trips3-tot_trips2);fflush(stdout);
		if(s_in[dest]<1)
		{
			s_in[dest] = s_in[in_av-1];
			s_in[in_av-1] = 0;
			in_inds[dest] = in_inds[in_av-1];
			in_av--;
			if(verbose>0)
			{
				printf("%d Exhausted. Available in %d, Available out %d  2go: %d\n",aux2,in_av,out_av,tot_trips3-tot_trips2);fflush(stdout);
			}
		}
		if(s_out[origin]<1)
		{
			s_out[origin] = s_out[out_av-1];
			s_out[out_av-1] = 0;
			out_inds[origin] = out_inds[out_av-1];
			out_av--;
			if(verbose>0) printf("%d Exhausted. Available in %i, Available out %i 2go: %d\n",aux1,in_av,out_av,tot_trips3-tot_trips2);fflush(stdout);
		}
		if(out_av<=0)
		{
			if(in_av<=0)
			{
				flag=0;
			}else{
				flag=1;
				printf("Aborted: There are still available out-trips but no available in trips!\n");
				abort();
			}
		}
		if(in_av<=0)
		{
			if(out_av<=0)
			{
				flag=0;
			}else{
				flag=1;
				printf("Aborted: There are still available in-trips but no available out trips!\n");
				abort();
			}
		}
		if((out_av==1) && (in_av==1) && (out_inds[0]==in_inds[0]))
		{
			if(self_opt<=0)
			{
				flag=1;// avoid last self loop
				if(verbose>0) printf("Last trips are self-loops, ignoring last %i trips!\n",tot_trips3-tot_trips2);
			}
		}
		if(tot_trips3-tot_trips2<0)
		{
			flag=3;
			if(verbose>0) printf("Aborted: Something is wrong... too many trips!\n");
			abort();
		}
		if((tot_trips2%50000==0)&&(verbose>0))
		{
			printf("Done %d trips, %f frac\n",tot_trips2,tot_trips2/(double)tot_trips3);
		}
	}while(flag<0);
	free(s_in);
	free(s_out);
	free(in_inds);
	free(out_inds);
    if(verbose>0) printf("Total trips for this round: %i, Suposed to be: %i\n",tot_trips2,tot_trips3);
	return WG;
}

// ula_multinomial_implementation //
W_GRAPH* w_graph_seq_gravity_multinomial_directed(int N_nodes,int **s, double **d, double gamma, gsl_rng* randgsl, int self_opt, int verbose){
    // alloc w_graph
    W_GRAPH* WG  = w_graph_alloc(N_nodes);
    WG->opt_self = self_opt;
    WG->opt_dir = 1;

    int i;
    int out_av,in_av;

    // allocs (aux_vectors)
	int* s_in=cast_vec_int(N_nodes);
	int* s_out=cast_vec_int(N_nodes);
	int* in_inds=cast_vec_int(N_nodes);
	int* out_inds=cast_vec_int(N_nodes);
    
    // aux variables
	int tot_trips2,tot_trips3,aux1,aux2;

    aux1=in_av=0;
	aux2=out_av=0;
	tot_trips2=tot_trips3=0;
	
    // check all good
    for(i=0;i<N_nodes;i++)
	{
		if(s[1][i]!=0)
		{
			s_in[aux1]=s[1][i];
			in_inds[aux1]=i;
			aux1++;
			tot_trips2+=s[1][i];
		}else{
			in_inds[N_nodes-in_av-1]=i;
			in_av++;
		}
		if(s[0][i]!=0)
		{
			s_out[aux2]=s[0][i];
			out_inds[aux2]=i;
			aux2++;
			tot_trips3+=s[0][i];
		}else{
			out_inds[N_nodes-out_av-1]=i;
			out_av++;
		}
	}
	//tot_trips2=sum_vec_int(s_outf,N_nodes);
	//tot_trips3=sum_vec_int(s_inf,N_nodes);
	for(i=0;i<N_nodes;i++)
	{
		aux1=count_appearances_int(in_inds,i,N_nodes);
		aux2=count_appearances_int(out_inds,i,N_nodes);
		if((aux1!=1) || (aux2!=1)) // check 1 to 1 mapping
		{
			if(verbose>0)printf("Incorrect mapping for i: %i s_out: %d sin: %d %i %i\n",i,s[1][i],s[0][i],aux1,aux2);
			abort();
		}
	}
	if(verbose>0)
    {
        if(tot_trips2-tot_trips3 != 0)printf("Warning, total outs (%i) and total ins (%i) different!\n",tot_trips2,tot_trips3);
        printf("Available in %i, Available out %i\n",N_nodes-in_av,N_nodes-out_av);
    }

    //go! //
    int origin,dest,flag;

    // multinomial probabilities and events //
	double* p = cast_vec_double(N_nodes);
	//unsigned int* dummy_dest = cast_vec_uint(N_nodes);
	//double* ddd = cast_vec_double(N_nodes);
	//double rr,dum,dum2;
	double norm,sum_p;
	unsigned int nnn;
	in_av=N_nodes-in_av;
	out_av=N_nodes-out_av;

    flag=-1;
    tot_trips2=0;
    do{// While available outs
        aux1=aux2=-1;
        // select random origin node
        //origin = (int)round((out_av-1)*gsl_rng_uniform(randgsl)); // maxime's method
        origin  = gsl_rng_uniform_int(randgsl,(unsigned int)out_av); //out from out_av 
        aux1=out_inds[origin];
        // set to zero multinomial //
        //for(i=0;i<N_nodes;i++)
        //{
		//	dummy_dest[i]=0;
		//}
        // compute probabilities && allocate using multinomial
        if(in_av>1)
		{
			compute_p_french_multinomial(p,s_in, d, gamma, aux1, in_inds, in_av, N_nodes, self_opt, verbose);
			dest = -1;
			///////// maxime's way //////////
			/*dest =in_av-1;
			rr = gsl_rng_uniform(randgsl);
			dum = 0;
			dum2 = sum_vec_double(p,in_av);
			ddd[0]=0;
			for(i=1;i<in_av;i++)
			{
				dum = dum+p[i-1]/dum2;
				ddd[i]=dum;
			}
			for(i=0;i<(in_av-1);i++)
			{			
				if( (rr>=ddd[i]) && (rr<ddd[i+1]))
				{
					dest = i;
				}
			}
			//printf("### %d, %f %f\n",dest,ddd[in_av-1],ddd[dest]);
			/////////////////////////////
			*/
			
			///// my way 1///////////
			norm = sum_vec_double(p,in_av);
			sum_p = 0;
			for (i = 0; i < in_av; i++)
			{
				nnn =0;
				if (p[i] > 0.0)
				{
					nnn = gsl_ran_bernoulli (randgsl, p[i] / (norm-sum_p));
					sum_p += p[i];
				}
				if (nnn==1) 
				{
					dest = i;
					break;
				}
			}
			///// my way 2/////////// --> equivalent to my way 1
			//gsl_ran_multinomial (randgsl, N_nodes, 1, p, dummy_dest);
			//dest = find_value_uint(dummy_dest, 1, N_nodes);
			////////////////////////
			if(dest<0) //if not working, chose at random
			{
				if(verbose>0)
				{
					printf("# Something is not right, choosing at random... check your gamma \n");
				}
				//dest = (int)gsl_rng_uniform_int(randgsl,(unsigned int)in_av); //out from out_av				
				assert(dest>=0);
			}
			aux2 = in_inds[dest];
			//printf("%d %d %f\n",aux2,dest,p[dest]);

		}
		else{
			aux2 = in_inds[0]; // only one left
			dest = 0;
			if(self_opt<=0) assert(aux2!=aux1);
        }
        if((self_opt>0) || (aux1!=aux2)) // if different or accepting self-loops
        {      
			// add multi-link //
			//printf(" ### %d %d\n",aux1,aux2);
			w_graph_add_multi_link(WG, N_nodes, aux1, aux2, 1);
			tot_trips2++;
			s_in[dest]--;
			s_out[origin]--;
			//printf("Accepted trip! From %i to %i, s_out %i s_in %i To go %i\n",aux1,aux2,s_out[origin],s_in[dest],tot_trips3-tot_trips2);fflush(stdout);
			if(s_in[dest]<1)
			{
				s_in[dest] = s_in[in_av-1];
				s_in[in_av-1] = 0;
				in_inds[dest] = in_inds[in_av-1];
				in_av--;
				if(verbose>0)
				{
					printf("%d Exhausted. Available in %d, Available out %d  2go: %d\n",aux2,in_av,out_av,tot_trips3-tot_trips2);fflush(stdout);
				}
			}
			if(s_out[origin]<1)
			{
				s_out[origin] = s_out[out_av-1];
				s_out[out_av-1] = 0;
				out_inds[origin] = out_inds[out_av-1];
				out_av--;
				if(verbose>0) printf("%d Exhausted. Available in %i, Available out %i 2go: %d\n",aux1,in_av,out_av,tot_trips3-tot_trips2);fflush(stdout);
			}
			if(out_av<=0)
			{
				if(in_av<=0)
				{
					flag=0;
				}else{
					flag=1;
					printf("Aborted: There are still available out-trips but no available in trips!\n");
					abort();
				}
			}
			if(in_av<=0)
			{
				if(out_av<=0)
				{
					flag=0;
				}else{
					flag=1;
					printf("Aborted: There are still available in-trips but no available out trips!\n");
					abort();
				}
			}
			if((out_av==1) && (in_av==1) && (out_inds[0]==in_inds[0]))
			{
				if(self_opt<=0)
				{
					flag=1;// avoid last self loop
					if(verbose>0) printf("Last trips are self-loops, ignoring last %i trips!\n",tot_trips3-tot_trips2);
				}
			}
			if(tot_trips3-tot_trips2<0)
			{
				flag=3;
				if(verbose>0) printf("Aborted: Something is wrong... too many trips!\n");
				abort();
			}
		}
	}while(flag<0);
	//free(dummy_dest);
	free(p);
	free(s_in);
	free(s_out);
	free(in_inds);
	free(out_inds);
    if(verbose>0) printf("Total trips for this round: %i, Suposed to be: %i\n",tot_trips2,tot_trips3);
	return WG;
}
/*********************************************/
/************** Radiation Model ********************
 Simini et altr. Nature
 [1] F. Simini, M. C. González, A. Maritan, and A.-L. Barabási, Nature 484, 96 (2012).
*********************************************/
int destination_rad_model(int** s, double** d, int origin, int N_nodes, gsl_rng * randgsl)
{
	int i,j,dest;
	double max_val,r,node_val,dist2,dist=1e9,eps=1e-10;
	max_val=-1;
	//int finite=((double)s[0][origin]/(1.-((double)s[1][origin]/(double)T)));
	for(j=0;j<s[0][origin];j++)
	{//generate randoms for original node
		r=gsl_rng_uniform(randgsl);
		if(r>max_val)max_val=r;
	}
	node_val=max_val;
	dest=origin;
	for(i=0;i<N_nodes;i++)
	{//for each node (others)
		if(i!=origin)
		{
			max_val=-1;
			for(j=0;j<s[1][i];j++)
			{//generate s_in randoms
				r=gsl_rng_uniform(randgsl);
				if(r>max_val)max_val=r;
			}
			if(max_val>=node_val)
			{//keep the ones that exceed threshold
				dist2=d[maxeq_int(origin,i)][mineq_int(origin,i)]; //store dist
				if((dist2+eps)<=dist) // if closer
				{
					if(fabs(dist-dist2)<eps)//if equal undraw
					{
						r=gsl_rng_uniform(randgsl);
						if(r>0.5)
						{
							dist=dist2;
							dest=i;
						}
					}else{
						dist=dist2;
						dest=i;
					}
				}
			}
		}
	}
	return dest;
}
/*********************************************/
W_GRAPH*  w_graph_radiation_model_stochastic_directed(int N_nodes,int **s, double **d, gsl_rng* randgsl, int self_opt, int verbose)
{
	if(verbose>0) printf("-- Warning: This takes a while, if you are in a hurry, use the multinomial version(asymptotically equivalent)...-- \n");
    int i;
    // alloc graph
    W_GRAPH* WG = w_graph_alloc(N_nodes);
    WG->opt_self = self_opt;
    WG->opt_dir = 1;


	int tot_trips2,tot_trips3;
	int given_trips,dest;
	tot_trips2=tot_trips3=0;
	for(i=0;i<N_nodes;i++)
	{
		tot_trips2+=s[0][i];
		tot_trips3+=s[1][i];
	}
	if(tot_trips2-tot_trips3 != 0){
        if(verbose>0)printf("Warning, total outs (%i) and total ins (%i) different!\n",tot_trips2,tot_trips3);
    }
    tot_trips2=0;
    for(i=0;i<N_nodes;i++) // for each node (they are indep)
    {
        given_trips=s[0][i];
        //printf("Doing %i %i\n",i,given_trips);fflush(stdout);
        while(given_trips>0) // if it has indeed out-going trips!
        {
            dest = 	destination_rad_model(s, d, i, N_nodes,randgsl);
            given_trips--; //regardless of success
            if((dest!=i) || self_opt>0) //if success
            {
                tot_trips2++;
                w_graph_add_multi_link(WG, N_nodes, i, dest, 1);
            }
        }
        if((i%10==0) && (verbose>0))printf("---Performed %i nodes out of %i, %i trips to go---\n",i,N_nodes,tot_trips3-tot_trips2);
    }
    if(verbose>0) printf("Total trips for this round: %i, Suposed to be: %i\n",tot_trips2,tot_trips3);
	return WG;
}
/*********************************************/
W_GRAPH* w_graph_radiation_model_multinomial_directed(int N_nodes,int **s, double **d, gsl_rng* randgsl, int self_opt, int verbose, int ps_opt, double**ps)
{
	int i,j; // aux indices
    // alloc graph
    W_GRAPH* WG = w_graph_alloc(N_nodes);
    WG->opt_self = self_opt;
    WG->opt_dir = 1;

    int tot_trips2,tot_trips3,sij;
	tot_trips2=tot_trips3=0;
    /// prepare probs ///
	if(ps_opt<1)
	{
		for(i=0;i<N_nodes;i++)
		{
			tot_trips2+=s[0][i];
			tot_trips3+=s[1][i];
		}
		if((tot_trips2-tot_trips3 != 0) && (verbose>0)){
	        printf("Warning, total outs (%i) and total ins (%i) different!\n",tot_trips2,tot_trips3);
	    }
	
		for(i=0;i<N_nodes;i++) // each node
		{
			for(j=0;j<N_nodes;j++) // each available node
			{
				if(((s[0][i]>0) && (s[1][j]>0) ) && (i!=j)) // if sout_i and sin_j != 0
				{
					sij=compute_sij_rad(s,d,i,j,N_nodes);
					ps[i][j]=((double)s[1][i]*(double)s[1][j])/(((double)s[1][i]+(double)sij)*((double)s[1][i]+(double)sij+(double)s[1][j])); // probabilities
					//ps[i][j]=((double)s[0][i]*(double)s[1][j])/(((double)s[0][i]+(double)sij)*((double)s[0][i]+(double)sij+(double)s[1][j])); // probabilities
					ps[i][j]=ps[i][j]/(1.-((double)s[1][i]/(double)tot_trips2));
					assert(ps[i][j]>=0);
				}
				else
				{
					ps[i][j]=0;
				}
			}
	        if(self_opt>0)
	        {
				if(((s[0][i]>0) && (s[1][j]>0) )) // if sout_i and sin_j != 0
	            {
					sij=0;
					ps[i][i]=(double)(s[0][i]*s[1][i])/((double)(s[0][i]+sij)*(double)(s[0][i]+sij+s[1][i])); // probabilities
					//ps[i][i]=ps[i][i]*1./(1.-((double)s[1][i]/(double)tot_trips2));
				}else{
					ps[i][i] = 0;
				}
	        }else{
	            ps[i][i]=0.;
	        }
			//printf("%d\n",i);fflush(stdout);
			if((i%100==0) && (verbose>0))printf("---Performed %i nodes out of %i ---\n",i,N_nodes);
		}
	}
    /// on to simulate ///
    tot_trips2=0;
	unsigned int * ps_false = cast_vec_uint(N_nodes);
    for(i=0;i<N_nodes;i++)
    {
		if(s[0][i]>0)
		{
			gsl_ran_multinomial(randgsl, N_nodes, s[0][i], ps[i], ps_false); // run multinomial
			assert(s[0][i]==sum_vec_uint(ps_false,N_nodes));
			for(j=0;j<N_nodes;j++)
			{
				if (ps_false[j]>0)
				{
						w_graph_add_multi_link(WG, N_nodes, i, j, ps_false[j]);
						tot_trips2+=ps_false[j];
						ps_false[j]=0;
				}
			}
		}
    }
    // security checks //
    int ** s_aux = w_graph_compute_s(WG, N_nodes);
    int s1 = sum_vec_int(s_aux[0],N_nodes);
    int s2 = sum_vec_int(s[0],N_nodes);
    assert(s1==s2);
    //printf("ss : %d %d %d\n",s1,s2,tot_trips2);
    free_mat_int(s_aux,2);
    free(ps_false);
    if(verbose>0) printf("Total trips for this round: %i, Suposed to be: %i\n",tot_trips2,sum_vec_int(s[0],N_nodes));
	return WG;
}
/*********************************************/
int compute_sij_rad(int  **s, double ** dist, int origin, int dest, int N_nodes){
	int i;
	int sij=0;
	double dd=dist[maxeq_int(origin,dest)][mineq_int(origin,dest)];
	for(i=0;i<N_nodes;i++)
	{
		if((i!=origin)&&(i!=dest)&&(dist[maxeq_int(i,origin)][mineq_int(i,origin)]<dd))
		{
			//printf("d_ori_dest: %f d_actual:%f ori: %d dest: %d other:%d\n",dd,dist[maxeq_int(i,origin)][mineq_int(i,origin)],origin,dest,i);
			sij+=s[1][i];
		}
	}
/*
    if(sij==0)
    {
		printf("%d %d %d %f\n",sij,origin,dest,dd);
	}
*/  
	return sij;
}



/*********************************************/
/*********** Wilson Models  ********************
* Multinomial:
    [1] A. Wilson, Geogr. Anal. 42, 364 (2010).
* Poisson model:
    [2] Sagarra O., Pérez-Vicente C. and Diaz-Guilera A., PRE 88-6 1013.
*********************************************/
/*********************************************/
W_GRAPH* gravity_poisson_undirected_graph2(double*x,  int N_nodes , double **d, double gamma, gsl_rng* randgsl, int verbose, int self_opt){
    //assert(N_nodes%2==0); // assert N even
    unsigned int ps_false,aux2;
    double mu;
    int i,j,out_d, in_d;
    W_GRAPH* WG = w_graph_alloc(N_nodes);
    WG->opt_self = self_opt;
    WG->opt_dir = 0;

    aux2=0;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<i;j++)
        {
            out_d=maxeq_int(i,j);
            in_d=mineq_int(i,j);
            mu=x[j]*x[i]*exp(-gamma*d[out_d][in_d]);
            ps_false= gsl_ran_poisson(randgsl, mu);
            aux2+=ps_false;
            if(ps_false>0)
            {
                w_graph_add_multi_link(WG, N_nodes, i, j, ps_false);
                w_graph_add_multi_link(WG, N_nodes, j, i, ps_false);
            }
        }
        if(self_opt>0)
        {
            mu=x[i]*x[i];
            ps_false= gsl_ran_poisson(randgsl, mu);
            aux2+=ps_false;
            if(ps_false>0)
            {
                w_graph_add_multi_link(WG, N_nodes, i, i, ps_false);
            }
        }
    }
    double sout_mean,sout_std;
    double kout_mean,kout_std;
    int nulls_out;
    sout_mean=kout_mean=kout_std=sout_std=0;
    nulls_out=0;
    for(i=0;i<N_nodes;i++)
    {
        kout_mean+=(double)WG->node[i].kout;
        kout_std+=(double)WG->node[i].kout*WG->node[i].kout;
        sout_mean+=(double)WG->node[i].sout;
        sout_std+=(double)WG->node[i].sout*WG->node[i].sout;
        if(WG->node[i].kout==0)
        {
            nulls_out++;
        }
    }
    sout_mean/=(double)(N_nodes-nulls_out);
    sout_std=sqrt(sout_std/(double)(N_nodes-nulls_out)-sout_mean*sout_mean);
    kout_mean/=(double)(N_nodes-nulls_out);
    kout_std=sqrt(kout_std/(double)(N_nodes-nulls_out)-kout_mean*kout_mean);
	if(verbose==1)
	{
		printf("# I generated an uncorrelated weighted net\n");
		printf("# <s>_out=%.3lf+-%.3lf, <k>_out=%.3lf+-%.3lf\n",sout_mean, sout_std, kout_mean, kout_std);
		printf("# Number of Empty nodes: (out) %d\n",nulls_out);
		printf("# Total T:%d \n",aux2);
		printf("# Total Entropy S: %f\n",w_graph_entropy_multinomial(WG,N_nodes,self_opt));
	}
    return WG;
    
}

W_GRAPH* gravity_poisson_directed_graph2(double**x,  int N_nodes , double ** d, double gamma, gsl_rng* randgsl, int verbose, int self_opt){
    //assert(N_nodes%2==0); // assert N even
    //unsigned int* ps_false=(unsigned int*)malloc(sizeof(unsigned int)*N_nodes*(N_nodes-1)/2);
    unsigned int ps_false,aux2;
    double mu;
    int i,j,out_d, in_d;
    //gsl_ran_multinomial(randgsl, N_nodes*(N_nodes-1)/2,  T, ps, ps_false);
    //free(ps);
    W_GRAPH* WG = w_graph_alloc(N_nodes);
    WG->opt_self = self_opt;
    WG->opt_dir = 1;

    aux2=0;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<N_nodes;j++)
        {
            out_d=maxeq_int(i,j);
            in_d=mineq_int(i,j);
            mu=x[1][j]*x[0][i]*exp(-gamma*d[out_d][in_d]);
            ps_false= gsl_ran_poisson(randgsl, mu);
            aux2+=ps_false;
            if((i!=j) || (self_opt>0))
            {
                if(ps_false>0)
                {
                    w_graph_add_multi_link(WG, N_nodes, i, j, ps_false);
                }
            }
        }
    }
    double sout_mean,sout_std,sin_mean,sin_std;
    double kout_mean,kout_std,kin_mean,kin_std;
    int nulls_out,nulls_in;
    sout_mean=kout_mean=kout_std=sout_std=sin_mean=kin_mean=sin_std=kin_std=0;
    nulls_out=nulls_in=0;
    for(i=0;i<N_nodes;i++)
    {
        kout_mean+=(double)WG->node[i].kout;
        kin_mean+=(double)WG->node[i].kin;
        kout_std+=(double)WG->node[i].kout*WG->node[i].kout;
        kin_std+=(double)WG->node[i].kin*WG->node[i].kin;
        sout_mean+=(double)WG->node[i].sout;
        sin_mean+=(double)WG->node[i].sin;
        sout_std+=(double)WG->node[i].sout*WG->node[i].sout;
        sin_std+=(double)WG->node[i].sin*WG->node[i].sin;
        if(WG->node[i].kout==0)
        {
            nulls_out++;
        }
        if(WG->node[i].kin==0)
        {
            nulls_in++;
        }
    }
    sout_mean/=(double)(N_nodes-nulls_out);
    sin_mean/=(double)(N_nodes-nulls_in);
    sout_std=sqrt(sout_std/(double)(N_nodes-nulls_out)-sout_mean*sout_mean);
    sin_std=sqrt(sin_std/(double)(N_nodes-nulls_in)-sin_mean*sin_mean);
    kout_mean/=(double)(N_nodes-nulls_out);
    kin_mean/=(double)(N_nodes-nulls_in);
    kout_std=sqrt(kout_std/(double)(N_nodes-nulls_out)-kout_mean*kout_mean);
    kin_std=sqrt(kin_std/(double)(N_nodes-nulls_in)-kin_mean*kin_mean);
	if(verbose==1)
	{
		printf("# I generated an uncorrelated weighted net\n");
		printf("# (out) <s>=%.3lf+-%.3lf, <k>=%.3lf+-%.3lf\n",sout_mean, sout_std, kout_mean, kout_std);
		printf("# (in)  <s>=%.3lf+-%.3lf, <k>=%.3lf+-%.3lf\n",sin_mean, sin_std, kin_mean, kin_std);
		printf("# Number of Empty nodes: (out) %d \t (in) %d \n",nulls_out,nulls_in);
		printf("# Total T:%d \n",aux2);
		printf("# Total Entropy S: %f\n",w_graph_entropy_multinomial(WG,N_nodes,self_opt));

	}
    return WG;
    
}

/*********************************************/
/************ Probabilities choice ***********/
/*********************************************/
double * prob_mult_C_undir(int N_nodes, double gamma, double** d){
    double* ps=(double*)malloc(sizeof(double)*N_nodes*(N_nodes-1)/2);
    int i,j,aux,a,b;
    aux=0;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<i;j++)
        {   
            a=maxeq_int(i,j);
            b=mineq_int(i,j);    
            ps[aux]=exp(-gamma*d[a][b]);
            aux++;
        }
    }
    //assert(fabs((double)T-check)<1e-10);
    assert(aux==N_nodes*(N_nodes-1)/2);    
    return ps;
}
double * prob_mult_C_dir(int N_nodes, double gamma, double** d){
    double* ps=(double*)malloc(sizeof(double)*N_nodes*(N_nodes-1));
    int i,j,aux,a,b;
    aux=0;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<N_nodes;j++)
        {
            if(j!=i)
            {
                a=maxeq_int(i,j);
                b=mineq_int(i,j);    
                ps[aux]=exp(-gamma*d[a][b]);
                aux++;
            }
        }
    }
    //assert(fabs((double)T-check)<1e-10);
    assert(aux==N_nodes*(N_nodes-1));    
    return ps;
}

double * prob_mult_Cs_undir(int N_nodes, double* x, double gamma, double** d){
    double* ps=(double*)malloc(sizeof(double)*N_nodes*(N_nodes-1)/2);
    int i,j,aux,a,b;
    aux=0;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<i;j++)
        {   
            a=maxeq_int(i,j);
            b=mineq_int(i,j);    
            ps[aux]=x[i]*x[j]*exp(-gamma*d[a][b]);
            aux++;
        }
    }
    //assert(fabs((double)T-check)<1e-10);
    assert(aux==N_nodes*(N_nodes-1)/2);    
    return ps;
}
double * prob_mult_Cs_dir(int N_nodes, double **x, double gamma, double** d){
    double* ps=(double*)malloc(sizeof(double)*N_nodes*(N_nodes-1));
    int i,j,aux,a,b;
    aux=0;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<N_nodes;j++)
        {
            if(j!=i)
            {
                a=maxeq_int(i,j);
                b=mineq_int(i,j);    
                ps[aux]=x[0][i]*x[1][j]*exp(-gamma*d[a][b]);
                aux++;
            }
        }
    }
    //assert(fabs((double)T-check)<1e-10);
    assert(aux==N_nodes*(N_nodes-1));    
    return ps;
}
