/************************************************************
 *
 *                    Read Library
 *
 *		Functions useful when ...Reading node or edge attributes from files
 *
 *
 *************************************************************/



#include "ula_read_funcs.h"

/****************************************************************************
 *  ...Reading functions *
 ****************************************************************************/

double** read_distances(char *input_name, int num_nodes, int header, int opt_log){
	printf("...Reading distance file...\n");
	FILE* input=open_file("r",input_name);
	int i,j,di,dj;
	double dij;
	double** d=(double**)malloc(sizeof(double*)*num_nodes);
	for(i=0;i<num_nodes;++i)d[i]=(double*)malloc(sizeof(double)*(i+1)); // triangular matrix!
	int n=0;
	int k=0;
	char dummy[100];	
	while(k<header)
	{
		 fgets(dummy, 100, input);
		 k++;
	}
	if(opt_log>0)
	{
		while (!feof(input))
		{///we start ...Reading
			if(fscanf(input, "%d %d %lf\n", &i, &j, &dij)!=3)
			{
				printf("Problem in ...Reading\n");
				abort();
			}else{
				di = maxeq_int(i,j); // for distance (it is triangular)
				dj = mineq_int(i,j); // idem
				//printf("%i %i %lf\n",i,j,dij);fflush(stdout);
	  			d[di][dj]=log(dij);
	  			n++;
			}
		}		
	}else{
		while (!feof(input))
		{///we start ...Reading
			if(fscanf(input, "%d %d %lf\n", &i, &j, &dij)!=3)
			{
				printf("Problem in ...Reading\n");
				abort();
			}else{
				di = maxeq_int(i,j); // for distance (it is triangular)
				dj = mineq_int(i,j); // idem
				//printf("%i %i %lf\n",i,j,dij);fflush(stdout);
	  			d[di][dj]=dij;
	  			n++;
			}
		}
		
	}
	/*
	if(n!=num_nodes*num_nodes && n!=num_nodes*(num_nodes+1)/2 && n!=num_nodes*(num_nodes-1)/2){
		int dec_nodes=(int)(-1+sqrt(4*n+1))/2;
		printf("Check your distance list, its different (size %d, makes %d nodes) than the declared number of nodes (%d)!\n",n,dec_nodes,num_nodes);
		//abort();
		}
	*/
	//printf("...Read %i distances...\n",n);
	fclose(input);
	return d;
}




int** read_net_list(char *input_name, int num_nodes, int header){
	printf("...Reading net (NxN) list integer file...\n");
	FILE* input=open_file("r",input_name);
	int** d=cast_mat_int(num_nodes,num_nodes);
	int i,j,dij,k;
	int n=0;
	char dummy[100];
	k=0;	
	while(k<header)
	{
		 fgets(dummy, 100, input);
		 k++;
	}
	while (!feof(input))
	{		///we start ...Reading	
		if(fscanf(input, "%d %d %d\n", &i, &j, &dij)!=3)
		{
			printf("error while reading. Aborting.\n");
			abort();
		}else{
	  		d[i][j]=dij;
	  		n+=dij;
	  		//printf("%i",n);fflush(stdout);			
		}
		k++;
  	}
	printf("Total sum of attribute %i\n",n);
	fclose(input);
	return d;
}

double** read_net_list_double(char *input_name, int num_nodes, int header){
	printf("...Reading net (NxN) list file with float arguments ...\n");
	FILE* input=open_file("r",input_name);
	double** d=cast_mat_double(num_nodes,num_nodes);
	int i,j,k;
	double dij;
	double n=0;
	char dummy[100];
	k=0;	
	while(k<header)
	{
		 fgets(dummy, 100, input);
		 k++;
	}
	
	while (!feof(input))
	{		///we start ...Reading	
		if(fscanf(input, "%d %d %lf\n", &i, &j, &dij)!=3)
		{
			printf("error while reading. Aborting.\n");
			abort();
		}else{
			//printf("%lf\n",dij);fflush(stdout);
	  		d[i][j]=dij;
	  		n+=dij;
	  		//printf("%i",n);fflush(stdout);			
		}
  	}
	printf("Total sum of attribute %lf\n",n);
	fclose(input);
	return d;
}

/* Bad one */
/*
double*** read_edge_list_double(char *input_name, int num_nodes, int header){ // check
	printf("...Reading edge list file with float arguments ...\n");
	FILE* input=open_file("r",input_name);
	double** d=cast_mat_double(num_nodes,num_nodes);
	double** d2=cast_mat_double(num_nodes,num_nodes);
	int i,j,k;
	double dij,dij2;
	double n=0;
	char dummy[100];
	k=0;	
	while(k<header)
	{
		 fgets(dummy, 100, input);
		 k++;
	}
	
	while (!feof(input))
	{		///we start ...Reading	
		if(fscanf(input, "%d %d %lf %lf\n", &i, &j, &dij, &dij2)!=4)
		{
			printf("error while reading. Aborting.\n");
			abort();
		}else{
	  		d[i][j]=dij;
	  		d2[i][j]=dij2;
	  		n+=dij;
	  		//printf("%i",n);fflush(stdout);			
		}
  	}
	printf("Total average num of trips %lf\n",n);
	fclose(input);
	double *** dtot=(double***)malloc(sizeof(double**)*2);
	dtot[0]=d;
	dtot[1]=d2;
	return dtot;
}
*/

/************************************/

int** read_node_list_int(char *input_name,int num_nodes,int header){
	printf("...Reading node directed attribute file...\n");
	FILE* input=open_file("r",input_name);
	int** s=(int**)malloc(sizeof(int*)*2);
	s[0]=(int*)malloc(sizeof(int)*num_nodes);
	s[1]=(int*)malloc(sizeof(int)*num_nodes);
	int n,out,in,k;
	int m=0;
	int att_in,att_out;
	int flag=0;
	att_in=att_out=0;
	k=0;
	char dummy[100];	
	while(k<header)
	{
		 fgets(dummy, 100, input);
		 k++;
	}
	
	while (!feof(input))
	{		///we start ...Reading	
		if(fscanf(input, "%d %d %d\n", &n, &out, &in)!=3)
		{
			printf("error while reading. Aborting.\n");
			abort();
		}else{
			//assert(fscanf(input, "%d %d %d\n", &n, &out, &in)==3);
  			if(m!=n)flag=1;
  			s[0][m]=out;
  			s[1][m]=in;
  			m++;
  			att_in+=in;
  			att_out+=out;
		}
 	}
	if(num_nodes!=m){
		printf("Check your node list, its different (size %d) than the declared number of nodes (%d)!\n",m,num_nodes);
		abort();
	}else if(flag==1){printf("Node numbers are not congruent with given list!");
	}
	printf("Total number of nodes %i \t Total in_att %i \t Total out_att %i Difference %i\n",m,att_in,att_out,abs(att_in-att_out));
	fclose(input);
	return s;
}

int* read_node_list_int_undir(char *input_name,int num_nodes,int header){
	printf("...Reading node integer undirected attribute file...\n");
	FILE* input=open_file("r",input_name);
	int* s=(int*)malloc(num_nodes*sizeof(int));
	int n,out,k;
	int m=0;
	int att_out;
	int flag=0;
	att_out=0;
	k=0;
	char dummy[100];	
	while(k<header)
	{
		 fgets(dummy, 100, input);
		 k++;
	}
	while (!feof(input))
	{		///we start ...Reading	
		if(fscanf(input, "%d %d\n", &n, &out)!=2)
		{
			abort();
		}else{
			//assert(fscanf(input, "%d %d %d\n", &n, &out, &in)==3);
	  		if(m!=n)flag=1;
	  		s[m]=out;
	  		m++;
	  		att_out+=out;			
		}
 	}
	if(num_nodes!=m){
		printf("Check your node list, its different (size %d) than the declared number of nodes (%d)!\n",m,num_nodes);
		abort();
	}else if(flag==1){printf("Node numbers are not congruent with given list!");
	}
	printf("Total number of nodes %i \t Total att %i \n",m,att_out);
	fclose(input);
	return s;
}


double* read_node_list_double_undir(char *input_name,int num_nodes, int header){
	printf("...Reading node undirected float attribute file...\n");
	FILE* input=open_file("r",input_name);
	double* s=(double*)malloc(num_nodes*sizeof(double));
	int n,k;
	int m=0;
	double att_out,out;
	int flag=0;
	att_out=0;
	k=0;
	char dummy[100];	
	while(k<header)
	{
		 fgets(dummy, 100, input);
		 k++;
	}
	while (!feof(input))
	{		///we start ...Reading	
		if(fscanf(input, "%d %lf\n", &n, &out)!=2)
		{
			abort();
		}else{
			//assert(fscanf(input, "%d %d %d\n", &n, &out, &in)==3);
	  		if(m!=n)flag=1;
	  		s[m]=out;
	  		m++;
	  		att_out+=out;			
		}
 	}
	if(num_nodes!=m){
		printf("Check your node list, its different (size %d) than the declared number of nodes (%d)!\n",m,num_nodes);
		abort();
	}else if(flag==1){printf("Node numbers are not congruent with given list!");
	}
	printf("Total number of nodes %i \t Total att %lf \n",m,att_out);
	fclose(input);
	return s;
}



/************************************/

double** read_node_list_double(char *input_name,int num_nodes, int header){
	printf("...Reading node float directed attribute file...\n");
	FILE* input=open_file("r",input_name);
	double ** s=(double**)malloc(sizeof(double*)*2);
	s[0]=(double*)malloc(sizeof(double)*num_nodes);
	s[1]=(double*)malloc(sizeof(double)*num_nodes);
	int n,k;
	int m=0;
	double att_in,att_out,in,out;
	int flag=0;
	att_in=att_out=0;
	k=0;
	char dummy[100];	
	while(k<header)
	{
		 fgets(dummy, 100, input);
		 k++;
	}
	while (!feof(input))
	{		///we start ...Reading	
		if(fscanf(input, "%d %lf %lf\n", &n, &out, &in)!=3)
		{
			abort();
		}else{
			//assert(fscanf(input, "%d %lf %lf\n", &n, &out, &in)==3);
	  		if(m!=n)flag=1;
	  		s[0][m]=out;
	  		s[1][m]=in;
	  		m++;
	  		att_in+=in;
	  		att_out+=out;			
		}
 	}
	if(num_nodes!=m){
		printf("Check your node list, its different (size %d) than the declared number of nodes (%d)!\n",m,num_nodes);
		abort();
	}else if(flag==1){printf("Node numbers are not congruent with given list!");
	}
	printf("Total number of nodes %i \t Total in_att %lf \t Total out_att %lf Difference %lf\n",m,att_in,att_out,fabs(att_in-att_out));
	fclose(input);
	return s;
}

double** read_node_list_xatts_double(char *input_name,int num_nodes, int num_atts, int header){
	printf("...Reading node float attribute file with %d atts...\n",num_atts);
	FILE* input=open_file("r",input_name);
	double ** s= cast_mat_double(num_atts,num_nodes);
	double aux;
	//double * dum= cast_vec_double(num_atts-1);
	int n,k,kk;
	int m=0;
	int flag=0;
	k=0;
	char dummy[100];	
	while(k<header)
	{
		 fgets(dummy, sizeof(dummy), input);
		 k++;
	}
	while (!feof(input))
	{		///we start ...Reading	
		fgets(dummy, sizeof(dummy), input); // read to dummy
		//printf("%s",dummy);
		sscanf(dummy,"%d%[^\n]",&n,dummy); // read node
		//printf("%d # %s\n",n,dummy);
		for(kk=0;kk<num_atts;kk++)
		{
			sscanf(dummy,"%lf%[^\n]",&aux,dummy);
			//printf("%lf ",aux);
			s[kk][m] = aux;
		}
		//printf("\n");
		//assert(fscanf(input, "%d %lf %lf\n", &n, &out, &in)==3);
	  	if(m!=n)flag=1;
	  	m++;
	  	if(m>=num_nodes) break;
 	}
	if(flag==1){
		printf("Node numbers are not congruent with given list!");
	}
	fclose(input);
	return s;
}

int** read_node_list_xatts_int(char *input_name,int num_nodes, int num_atts, int header){
	printf("...Reading node int attribute file with %d atts...\n",num_atts);
	FILE* input=open_file("r",input_name);
	int ** s= cast_mat_int(num_nodes, num_atts);
	//double * dum= cast_vec_double(num_atts-1);
	int n,k,kk;
	int m=0;
	int flag=0;
	k=0;
	char dummy[100];	
	while(k<header)
	{
		 fgets(dummy, sizeof(dummy), input);
		 k++;
	}
	while (!feof(input))
	{		///we start ...Reading	
		fgets(dummy, sizeof(dummy), input); // read to dummy
		sscanf(dummy,"%d ",&n); // read node
		for(kk=0;kk<num_atts;kk++)
		{
			sscanf(dummy,"%d ",&s[kk][m]);
		}
		//assert(fscanf(input, "%d %lf %lf\n", &n, &out, &in)==3);
	  	if(m!=n)flag=1;
	  	m++;
 	}
	if(num_nodes!=m){
		printf("Check your node list, its different (size %d) than the declared number of nodes (%d)!\n",m,num_nodes);
		abort();
	}else if(flag==1){
		printf("Node numbers are not congruent with given list!");
	}
	fclose(input);
	return s;
}


/****************************************************************************
 *  Distance functions *
 ****************************************************************************/

double* nearest_neighbor_dist(double ** dist, int num_nodes)
{
	int i,j;
	int di,dj;
	double* near=cast_vec_double(num_nodes);
	double mind,av;
	av = 0;
	for(i=0;i<num_nodes;i++)
	{
		mind =  1e9;
		for(j=0;j<num_nodes;j++)
		{
			if(i!=j)
			{
				di = maxeq_int(i,j); // for distance (it is triangular)
				dj = mineq_int(i,j); // idem
				mind = fmin(mind,dist[di][dj]);
			}
		}
		near[i] = mind;
		av+=mind;
		//printf("near: %f\n",near[i]);
	}
	//printf("Average: %f\n",av/(double)num_nodes);
	return near;
}

int* nearest_neighbor(double ** dist, int num_nodes)
{
	int i,j,k;
	int di,dj;
	double mind;
	int* near=cast_vec_int(num_nodes);
	for(i=0;i<num_nodes;i++)
	{
		mind =  1e9;
		k=-1;
		for(j=0;j<num_nodes;j++)
		{
			if(i!=j)
			{
				di = maxeq_int(i,j); // for distance (it is triangular)
				dj = mineq_int(i,j); // idem
				if(mind<fmin(mind,dist[di][dj]))
				{
					k=j;
				}
				mind = fmin(mind,dist[di][dj]);
			}
		}
		near[i] = k;
	}
	return near;
}
