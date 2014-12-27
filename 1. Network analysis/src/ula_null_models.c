/************************************************************
 *
 *                    Null Models for multi-edge 
 *                   undirected &directed networks
 *
 *
 *
 *************************************************************/

#include "ula_null_models.h"

/*********************************************/
/************** Generation of probabilities  ********************/
/*********************************************/

double * prob_mult_s_undir(double* x, int N_nodes, int self_opt){
    double* ps;
    if(self_opt>0)
    {
        ps=(double*)malloc(sizeof(double)*(int)(N_nodes*(N_nodes+1)/2.));
    }else{
        ps=(double*)malloc(sizeof(double)*(int)(N_nodes*(N_nodes-1)/2.));        
    }
    int i,j,aux;
    aux=0;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<i;j++)
        {
            ps[aux]=x[j]*x[i];
            aux++;
        }
        //printf("%d %d\n",i,aux);fflush(stdout);
        if (self_opt>0)
        {
            ps[aux]=x[i]*x[i];
            aux++;
        }
    }
    if(self_opt>0)
    {
        assert(aux==(int)(N_nodes*(N_nodes+1)/2.));
    }else{
        assert(aux==(int)(N_nodes*(N_nodes-1)/2.));
    }
    return ps;
}
double * prob_mult_s_dir(double** x, int N_nodes, int self_opt){
    double* ps;
    if(self_opt>0)
    {
        ps=(double*)malloc(sizeof(double)*N_nodes*N_nodes);
    }else{
        ps=(double*)malloc(sizeof(double)*N_nodes*(N_nodes-1));        
    }
    int i,j,aux;
    aux=0;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<N_nodes;j++)
        {
            if(j!=i)
            {
                ps[aux]=x[1][j]*x[0][i];
                aux++;
            }else{
                if(self_opt>0)
                {
                    ps[aux]=x[1][j]*x[0][i];
                    aux++;
                }
            }
        }
    }
    //assert(fabs((double)T-check)<1e-10);
    if(self_opt>0)
    {
        assert(aux==N_nodes*N_nodes);
    }else{
        assert(aux==N_nodes*(N_nodes-1));        
    }
    return ps;
}

int compute_T(int N_nodes, double av_k, int* x){
    int X=sum_vec_int(x,N_nodes);
    long long int X2=(long long int)X*(long long int)X;
    double mu=((double)N_nodes*av_k)/(2.*X2);
    long int T= (int)((mu+1)*log(mu+1)*0.5*X2);
    printf("T to be sorted: %ld -> <s>=%lf <k>=%lf mu=%lf X=%d X*X=%lld \n",T,(2.*T)/N_nodes,av_k,mu,X,X2);fflush(stdout);
    return T;
}


    
/*********************************************/
/************** Uncorrelated graphs  ********************/
/*********************************************/
W_GRAPH* uncorrelated_multinomial_directed_graph(double* ps, double**x, int N_nodes , int T, gsl_rng* randgsl, int verbose, int self_opt){
    unsigned int* ps_false;
    unsigned int dim;
    if(self_opt>0)
    {
        dim=N_nodes*N_nodes;
    }else{
        dim=N_nodes*(N_nodes-1);
    }
    int i,j,aux;
    ps_false=(unsigned int*)malloc(sizeof(unsigned int)*dim);
    gsl_ran_multinomial(randgsl, dim,  T, ps, ps_false);
    W_GRAPH* WG = w_graph_alloc(N_nodes);

	aux=0;
    for(i=0;i<N_nodes;i++)
    {    
        for(j=0;j<N_nodes;j++)
        {
            if(ps_false[aux]>0)
            {
                if((i!=j) || self_opt>0)
                {
                    w_graph_add_multi_link(WG, N_nodes, i, j, ps_false[aux]);
                }
            }
			aux++;
        }
    }
    free(ps_false);
    double sin_mean,sout_mean,sin_std,sout_std;
    double kin_mean,kout_mean,kin_std,kout_std;
    int nulls_in,nulls_out;
    sin_mean=sout_mean=kin_mean=kout_mean=kin_std=kout_std=sin_std=sout_std=0;
    nulls_in=nulls_out=0;
    for(i=0;i<N_nodes;i++)
    {
        kin_mean+=(double)WG->node[i].kin;
        kin_std+=(double)WG->node[i].kin*WG->node[i].kin;
        kout_mean+=(double)WG->node[i].kout;
        kout_std+=(double)WG->node[i].kout*WG->node[i].kout;
        sin_mean+=(double)WG->node[i].sin;
        sin_std+=(double)WG->node[i].sin*WG->node[i].sin;
        sout_mean+=(double)WG->node[i].sout;
        sout_std+=(double)WG->node[i].sout*WG->node[i].sout;
        if(WG->node[i].kin==0)
        {
            nulls_in++;
        }
        if(WG->node[i].kout==0)
        {
            nulls_out++;
        }
    }
	assert(sin_mean==sout_mean);
	assert(kin_mean==kout_mean);
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
	}
    return WG;
    
}


W_GRAPH* uncorrelated_poisson_multinomial_directed_graph(double* ps, double**x, int N_nodes , int T, gsl_rng* randgsl, int verbose, int self_opt){
    unsigned int T_prime = gsl_ran_poisson(randgsl, T);
    if(verbose==1) printf("## Poisson sorting of T: %d\n",T_prime);
    fflush(stdout);
    return uncorrelated_multinomial_directed_graph(ps, x, N_nodes, (int)T_prime, randgsl, verbose, self_opt);
    }

W_GRAPH* uncorrelated_multinomial_undirected_graph(double* ps, double*x,  int N_nodes , int T, gsl_rng* randgsl, int verbose, int self_opt){
    unsigned int* ps_false;
    unsigned int dim;
    if(self_opt>0)
    {
        dim=N_nodes*(N_nodes+1)/2;
    }else{
        dim=N_nodes*(N_nodes-1)/2;
    }
    ps_false=(unsigned int*)malloc(sizeof(unsigned int)*dim);
    int i,j,aux;
    T= T/2;
    if(verbose>0) printf("## Multinomial sorting of T: %d\n",T);fflush(stdout);
    gsl_ran_multinomial(randgsl, dim,  T, ps, ps_false);
    //free(ps);
    W_GRAPH* WG = w_graph_alloc(N_nodes);
    aux=0;

    for(i=0;i<N_nodes;i++)
    {    
        for(j=0;j<i;j++)
        {
            if(ps_false[aux]>0)
            {
                //w_graph_add_multi_link(node, N_nodes, i, j, ps_false[aux]);
                //w_graph_add_multi_link(node, N_nodes, j, i, ps_false[aux]);
                w_graph_add_multi_link_undirected(WG, N_nodes, j, i, ps_false[aux]);
            }
            aux++;
        }
        if(self_opt>0) // if accepting self-loops
        {
            if(ps_false[aux]>0)
            {
                w_graph_add_multi_link_undirected(WG, N_nodes, i, i, ps_false[aux]);
                //w_graph_add_multi_link(node, N_nodes, i, i, ps_false[aux]);
                //w_graph_add_multi_link(node, N_nodes, i, i, ps_false[aux]);
            }
            aux++;
        }
    }
    free(ps_false);
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
		printf("# Total Entropy S: %f\n",w_graph_entropy(WG,N_nodes));
    }
	return WG;
    
}

W_GRAPH* uncorrelated_poisson_multinomial_undirected_graph(double* ps, double*x,  int N_nodes , int T, gsl_rng* randgsl, int verbose, int self_opt){
    unsigned int T_prime = gsl_ran_poisson(randgsl, T);
    if(verbose==1) printf("## Poisson sorting of T: %u\n",T_prime);
    fflush(stdout);
    return uncorrelated_multinomial_undirected_graph(ps, x, N_nodes, (int)T_prime, randgsl, verbose, self_opt);
    }


W_GRAPH* uncorrelated_poisson_undirected_graph2(double*x,  int N_nodes , gsl_rng* randgsl, int verbose, int self_opt){
    //assert(N_nodes%2==0); // assert N even
    unsigned int ps_false,aux2;
    double mu;
    int i,j;
    int flag;
    flag=0;
    double T=sum_vec_double(x,N_nodes);
    W_GRAPH* WG = w_graph_alloc(N_nodes);
    aux2=0;
    for(i=0;i<N_nodes;i++)
    {    
        for(j=0;j<i;j++)
        {
            mu=x[j]*(x[i]/(double)T);
            ps_false= gsl_ran_poisson(randgsl, mu);
            aux2+=ps_false;
            if(ps_false>0)
            {
                w_graph_add_multi_link_undirected(WG, N_nodes, i, j, ps_false);
                //w_graph_add_multi_link(node, N_nodes, j, i, ps_false);
                //w_graph_add_multi_link(node, N_nodes, i, j, ps_false);
            }
        }
        if(self_opt>0)
        {
            mu=x[i]*(x[i]/(double)T);
            ps_false= gsl_ran_poisson(randgsl, mu);
            aux2+=ps_false;
            if(ps_false>0)
            {
                flag+=ps_false;
                w_graph_add_multi_link_undirected(WG, N_nodes, i, i, ps_false);
                //w_graph_add_multi_link(node, N_nodes, i, i, ps_false);
                //w_graph_add_multi_link(node, N_nodes, i, i, ps_false);
            }
        }
    }
    if(verbose==1) printf("Warning: There are %d self-loops\n",flag);fflush(stdout);
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
		printf("# Total T:%d (original : %f) \n",aux2,T/2);
		printf("# Total Entropy S: %f\n",w_graph_entropy(WG,N_nodes));
	}
    return WG;
    
}

W_GRAPH* uncorrelated_poisson_directed_graph2(double**x,  int N_nodes , gsl_rng* randgsl, int verbose, int self_opt){
    //assert(N_nodes%2==0); // assert N even
    //unsigned int* ps_false=(unsigned int*)malloc(sizeof(unsigned int)*N_nodes*(N_nodes-1)/2);
    unsigned int ps_false,aux2;
    double mu;
    int i,j;
    double T=sum_vec_double(x[0],N_nodes);
    //gsl_ran_multinomial(randgsl, N_nodes*(N_nodes-1)/2,  T, ps, ps_false);
    //free(ps);
    W_GRAPH* WG = w_graph_alloc(N_nodes);
    aux2=0;
    for(i=0;i<N_nodes;i++)
    {    
        for(j=0;j<N_nodes;j++)
        {
            mu=x[1][j]*(x[0][i]/(double)T);
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
	assert(sin_mean==sout_mean);
	assert(kin_mean==kout_mean);
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
		printf("# Total T:%d (original: %f)\n",aux2,T/2);
	}
    return WG;
    
}





W_GRAPH* uncorrelated_computational_directed_graph(int** s, int N_nodes , gsl_rng* randgsl, int verbose, int self_opt){
    int T=sum_vec_int(s[0], N_nodes);
    assert(sum_vec_int(s[1],N_nodes) == T);
    unsigned int* stubs_out=(unsigned int*)malloc(sizeof(unsigned int)*T);
    unsigned int* stubs_in=(unsigned int*)malloc(sizeof(unsigned int)*T);
    int i,j,aux1,aux2;
    aux1=aux2=0;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<s[0][i];j++)
        {
            stubs_out[aux1]=i;
            aux1++;
        }
        for(j=0;j<s[1][i];j++)
        {
            stubs_in[aux2]=i;
            aux2++;
        }
    }
    W_GRAPH* WG = w_graph_alloc(N_nodes);
    int togo=T;
    int dest,origin;
    int e,f;
    int trials,flag;
    int fake;
    flag=0;
    for(i=0;i<T;i++)
    {
        e  = gsl_rng_uniform_int(randgsl,togo); //stub from remaining trips
        origin = stubs_out[e];
        dest=origin;
        trials=0;
        if(self_opt<1){
            do
            {
                f=gsl_rng_uniform_int(randgsl,togo);
                dest  = stubs_in[f]; //stub from remaining trips
                trials++;
            }while((origin==dest) && (trials<=10*togo)); // we exclude the same node
            if(trials>togo)
            {
                flag++;
            }
        }else{
            f=gsl_rng_uniform_int(randgsl,togo);
            dest  = stubs_in[f]; //stub from remaining trips
        }
        fake=stubs_out[togo-1];
        stubs_out[togo-1]=origin;
        stubs_out[e]=fake;
        fake=stubs_in[togo-1];
        stubs_in[togo-1]=dest;
        stubs_in[f]=fake;
        togo--;
        //printf("got here %d!\n",T-togo);fflush(stdout);
        if(origin!=dest)
        {
            w_graph_add_multi_link(WG, N_nodes, origin, dest, 1) ;// add edge
        }else{
            if(self_opt>0)
            {
                w_graph_add_multi_link(WG, N_nodes, origin, dest, 1) ;// add edge
                flag+=1;
            }
        }
    }
    free(stubs_in);
    free(stubs_out);
    if(verbose==1) printf("Warning: There are %d self-loops\n",flag);fflush(stdout);
    double sin_mean,sout_mean,sin_std,sout_std;
    double kin_mean,kout_mean,kin_std,kout_std;
    int nulls_in,nulls_out;
    sin_mean=sout_mean=kin_mean=kout_mean=kin_std=kout_std=sin_std=sout_std=0;
    nulls_in=nulls_out=0;
	//printf("Start read: out: %d in : %d\n",s[0][0],s[1][0]);fflush(stdout);
	//printf("Start read: out: %d in : %d\n",node[0].sout,node[0].sin);fflush(stdout);
    for(i=0;i<N_nodes;i++)
    {
        kin_mean+=(double)WG->node[i].kin;
        kin_std+=(double)WG->node[i].kin*WG->node[i].kin;
        kout_mean+=(double)WG->node[i].kout;
        kout_std+=(double)WG->node[i].kout*WG->node[i].kout;
        sin_mean+=(double)WG->node[i].sin;
        sin_std+=(double)WG->node[i].sin*WG->node[i].sin;
        sout_mean+=(double)WG->node[i].sout;
        sout_std+=(double)WG->node[i].sout*WG->node[i].sout;
        if(WG->node[i].sin==0)
        {
            nulls_in++;
        }
        if(WG->node[i].sout==0)
        {
            nulls_out++;
        }
    }
	assert(sin_mean==sout_mean);
	assert(kin_mean==kout_mean);
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
    }
	return WG;
}



W_GRAPH* uncorrelated_computational_undirected_graph(int* s, int N_nodes , gsl_rng* randgsl, int max_trials, int verbose, int self_opt){
    int T=sum_vec_int(s, N_nodes);
    unsigned int* stubs_out=(unsigned int*)malloc(sizeof(unsigned int)*T);
    int i,j,aux1;
    //,aux2;
    aux1=0;
    for(i=0;i<N_nodes;i++)
    {
        for(j=0;j<s[i];j++)
        {
            stubs_out[aux1]=i;
            aux1++;
        }
    }
    W_GRAPH* WG = w_graph_alloc(N_nodes);
    int togo=T;
    int dest,origin;
    int e,f;
    int trials,flag;
    //int fake;
    flag=0;
    aux1=0;
    while(togo>0)
    {
        e  = gsl_rng_uniform_int(randgsl,togo); //stub from remaining trips
        origin = stubs_out[e];
        dest=origin;
        trials=0;
        if(self_opt<=0)
        {
            do
            {
                f=gsl_rng_uniform_int(randgsl,togo);
                dest  = stubs_out[f]; //stub from remaining trips
                trials++;
            }while(((dest==origin) && (trials<=10*togo)));
            //printf("to go: %d\n",togo);fflush(stdout); 
            if(trials>togo)
            {
                flag++;
            }
        }else{
            f=e;
            do
            {
                f=gsl_rng_uniform_int(randgsl,togo);
            }while(e==f);
            dest = stubs_out[f]; //stub from remaining trips
        }
        /*else if(trials>max_trials){
            printf("Failed rewiring steps, aborting!"); fflush(stdout);
            abort();
        }*/
        // substitution
        // cases! self loop or not
        if((dest!=origin) || (self_opt>0))
        {
            ///we copy the last two graph to the position where the graph just connected were
            
            if(e==togo-2 )
            {///en el cas que p1 sigui sum_k-2 primer copiem el vector[sum_k-2] a vector[p2] per borrarlo
                stubs_out[f]=stubs_out[togo-1];
            }else if(f==togo-1){
                stubs_out[e]=stubs_out[togo-2];
            }else if(f==togo-2){
                stubs_out[e]=stubs_out[togo-1];
            }else if(e==togo-1){
                stubs_out[f]=stubs_out[togo-2];
            }else{
                stubs_out[e]=stubs_out[togo-1];
                stubs_out[f]=stubs_out[togo-2];
            }
            togo-=2;
            w_graph_add_multi_link_undirected(WG, N_nodes, origin, dest, 1) ;// add edge
            //w_graph_add_multi_link(node, N_nodes, dest, origin, 1) ;// add edge
            //w_graph_add_multi_link(node, N_nodes, origin, dest, 1) ;// add edge
            if((self_opt>0)&&(origin==dest))
            {
                flag+=2;
                w_graph_add_multi_link_undirected(WG, N_nodes, origin, origin, 1) ;// add edge
            }
            //aux1++;
        }
    }
    free(stubs_out);
    if(verbose==1) printf("Warning: There are %d self-loops\n",flag);fflush(stdout);
    double sout_mean,sout_std;
    double kout_mean,kout_std;
    int nulls_out;
    sout_mean=kout_mean=kout_std=sout_std=0;
    nulls_out=0;
    for(i=0;i<N_nodes;i++)
    {
        kout_mean+=WG->node[i].kout;
        kout_std+=WG->node[i].kout*WG->node[i].kout;
        sout_mean+=WG->node[i].sout;
        sout_std+=WG->node[i].sout*WG->node[i].sout;
        if(WG->node[i].sout==0)
        {
            nulls_out++;
        }
    }
    //printf("T - aux: %d, aux-sout_mean:%f <s>=%f\n",T-aux1,T-sout_mean,(double)aux1/N_nodes); fflush(stdout);
    sout_mean/=(double)(N_nodes-nulls_out);
    sout_std=sqrt(sout_std/(double)(N_nodes-nulls_out)-sout_mean*sout_mean);
    kout_mean/=(double)(N_nodes-nulls_out);
    kout_std=sqrt(kout_std/(double)(N_nodes-nulls_out)-kout_mean*kout_mean);
    if(verbose==1)
	{
        printf("# I generated an uncorrelated weighted net\n");
		printf("# <s>_out=%.3lf+-%.3lf, <k>_out=%.3lf+-%.3lf\n",sout_mean, sout_std, kout_mean, kout_std);
		printf("# Number of Empty nodes: (out) %d\n",nulls_out);
		printf("# Total Entropy S: %f\n",w_graph_entropy(WG,N_nodes));
	}
    return WG;
}


//// custom p_ij poisson //
W_GRAPH* custompij_poisson_undirected_graph(double**pij,  int N_nodes , gsl_rng* randgsl, int verbose, int self_opt){
	// assuming pij is NOT normalized //
    //assert(N_nodes%2==0); // assert N even
    unsigned int ps_false,aux2;
    double mu;
    int i,j;
    int flag;
    flag=0;
    double T=sum_matrix_double(pij,N_nodes,N_nodes);
    W_GRAPH* WG = w_graph_alloc(N_nodes);
    aux2=0;
    for(i=0;i<N_nodes;i++)
    {    
        for(j=0;j<i;j++)
        {
            mu=pij[i][j];
            ps_false= gsl_ran_poisson(randgsl, mu);
            aux2+=ps_false;
            if(ps_false>0)
            {
                w_graph_add_multi_link_undirected(WG, N_nodes, i, j, ps_false);
                //w_graph_add_multi_link(node, N_nodes, j, i, ps_false);
                //w_graph_add_multi_link(node, N_nodes, i, j, ps_false);
            }
        }
        if(self_opt>0)
        {
            mu=(double)pij[i][j];
            ps_false= gsl_ran_poisson(randgsl, mu);
            aux2+=ps_false;
            if(ps_false>0)
            {
                flag+=ps_false;
                w_graph_add_multi_link_undirected(WG, N_nodes, i, i, ps_false);
                //w_graph_add_multi_link(node, N_nodes, i, i, ps_false);
                //w_graph_add_multi_link(node, N_nodes, i, i, ps_false);
            }
        }
    }
    if(verbose==1) printf("Warning: There are %d self-loops\n",flag);fflush(stdout);
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
		printf("# Total T:%d (original : %f) \n",aux2,T/2);
		printf("# Total Entropy S: %f\n",w_graph_entropy(WG,N_nodes));
	}
    return WG;
    
}
W_GRAPH* custompij_poisson_directed_graph(double**pij,  int N_nodes , gsl_rng* randgsl, int verbose, int self_opt){
	// assuming pij is NOT normalized //
    //assert(N_nodes%2==0); // assert N even
    //unsigned int* ps_false=(unsigned int*)malloc(sizeof(unsigned int)*N_nodes*(N_nodes-1)/2);
    unsigned int ps_false,aux2;
    double mu;
    int i,j;
    double T=sum_matrix_double(pij,N_nodes,N_nodes);
    //gsl_ran_multinomial(randgsl, N_nodes*(N_nodes-1)/2,  T, ps, ps_false);
    //free(ps);
    W_GRAPH* WG = w_graph_alloc(N_nodes);
    aux2=0;
    for(i=0;i<N_nodes;i++)
    {    
        for(j=0;j<N_nodes;j++)
        {
            mu=pij[i][j];
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
	assert(sin_mean==sout_mean);
	assert(kin_mean==kout_mean);
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
		printf("# Total T:%d (original: %f)\n",aux2,T);
	}
    return WG;
    
}



//// fixedEs poisson //
W_GRAPH* fixedEs_poisson_undirected_graph(double*x,  double lam, int N_nodes , gsl_rng* randgsl, int verbose, int self_opt, int max_reps){
    //assert(N_nodes%2==0); // assert N even
    unsigned int ps_false,aux2;
    double mu,p,pij;
    int i,j,reps;
    int flag;
    flag=0;
    double T=sum_vec_double(x,N_nodes);
    W_GRAPH* WG = w_graph_alloc(N_nodes);
    aux2=0;
    for(i=0;i<N_nodes;i++)
    {    
        for(j=0;j<i;j++)
        {
			// bernouilli trial first //
			pij = lam * (exp(x[i]*x[j])-1.) / (1.+lam * (exp(x[i]*x[j])-1.));
			p = gsl_rng_uniform(randgsl);
			if(p>pij)
			{
	            mu= x[i]*x[j]*lam*exp(x[i]*x[j]) / (1.+lam * (exp(x[i]*x[j])-1.));
				reps=0;
				ps_false=0;
				do{
	            	ps_false= gsl_ran_poisson(randgsl, mu);
					reps++;
				}while(ps_false<=0  || reps<max_reps);
	            aux2+=ps_false;
	            w_graph_add_multi_link_undirected(WG, N_nodes, i, j, ps_false);
	        }				
		}

        if(self_opt>0)
        {
			// bernouilli trial first //
			pij = lam * (exp(x[i]*x[i])-1.) / (1.+lam * (exp(x[i]*x[i])-1.));
			p = gsl_rng_uniform(randgsl);
			if(p>pij)
			{
				mu= x[i]*x[i]*lam*exp(x[i]*x[i]) / (1.+lam * (exp(x[i]*x[i])-1.));
				reps=0;
				ps_false=0;
				do{
	            	ps_false= gsl_ran_poisson(randgsl, mu);
					reps++;
				}while(ps_false<=0 || reps<max_reps);
	            aux2+=ps_false;
	            w_graph_add_multi_link_undirected(WG, N_nodes, i, j, ps_false);
	        }
		}
    }
    if(verbose==1) printf("Warning: There are %d self-loops\n",flag);fflush(stdout);
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
		printf("# I generated a zip weighted net\n");
		printf("# <s>_out=%.3lf+-%.3lf, <k>_out=%.3lf+-%.3lf\n",sout_mean, sout_std, kout_mean, kout_std);
		printf("# Number of Empty nodes: (out) %d\n",nulls_out);
		printf("# Total T:%d (original : %f) \n",aux2,T/2);
		printf("# Total Entropy S: %f\n",w_graph_entropy(WG,N_nodes));
	}
    return WG;
    
}
W_GRAPH* fixedEs_poisson_directed_graph(double**x, double lam,  int N_nodes , gsl_rng* randgsl, int verbose, int self_opt, int max_reps){
    //assert(N_nodes%2==0); // assert N even
    //unsigned int* ps_false=(unsigned int*)malloc(sizeof(unsigned int)*N_nodes*(N_nodes-1)/2);
    unsigned int ps_false,aux2,reps;
    double mu,p,pij;
    int i,j;
    double T=sum_vec_double(x[0],N_nodes);
    //gsl_ran_multinomial(randgsl, N_nodes*(N_nodes-1)/2,  T, ps, ps_false);
    //free(ps);
    W_GRAPH* WG = w_graph_alloc(N_nodes);
    aux2=0;
    for(i=0;i<N_nodes;i++)
    {    
        for(j=0;j<N_nodes;j++)
        {
			// bernouilli trial first //
			pij = lam * (exp(x[0][i]*x[1][j])-1.) / (1.+lam * (exp(x[0][i]*x[1][j])-1.));
			p = gsl_rng_uniform(randgsl);
			if(p>pij)
			{
	            mu= x[0][i]*x[1][j]*lam*exp(x[0][i]*x[1][j]) / (1.+lam * (exp(x[0][i]*x[1][j])-1.));
				reps=0;
				ps_false=0;
				do{
	            	ps_false= gsl_ran_poisson(randgsl, mu);
					reps++;
				}while(ps_false<=0 || reps<max_reps);
	            aux2+=ps_false;
	        }				
            if((i!=j) || (self_opt>0))
            {
                w_graph_add_multi_link(WG, N_nodes, i, j, ps_false);
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
	assert(sin_mean==sout_mean);
	assert(kin_mean==kout_mean);
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
		printf("# I generated a zip weighted net\n");
		printf("# (out) <s>=%.3lf+-%.3lf, <k>=%.3lf+-%.3lf\n",sout_mean, sout_std, kout_mean, kout_std);
		printf("# (in)  <s>=%.3lf+-%.3lf, <k>=%.3lf+-%.3lf\n",sin_mean, sin_std, kin_mean, kin_std);
		printf("# Number of Empty nodes: (out) %d \t (in) %d \n",nulls_out,nulls_in);
		printf("# Total T:%d (original: %f)\n",aux2,T/2);
	}
    return WG;
    
}




