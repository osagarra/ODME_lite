/************************************************************
 *
 *                    General functions library
 *
 *		General purpouse programming functions
 *
 *
 *************************************************************/

#include "ula_gen_funcs.h"


/********************************************************************************
 *  FILE HANDLING
********************************************************************************/


FILE* open_file(char *mode,char *nom){
	
	FILE *fitxer;
	fitxer=fopen(nom,mode);
		if(fitxer==NULL)
	   {
	   printf( "cannot open file %s.\n",nom );
	   abort();
	   }
	
	return fitxer;

}

/********************************************************************************
 ********************************************************************************/
/********************************************************************************
 *  REALLOC HANDLING
 ********************************************************************************/

///un realloc segur que també copia
void* xrealloc(void* vector, size_t midaenbytes){
    void *tmpvec=NULL;
    tmpvec=realloc(vector,midaenbytes);
    if (tmpvec==NULL)
    {
        fprintf(stderr,"# Error al fer realloc. Avortant\n");
        free(vector);
        abort();
    }
    return tmpvec;
}

int* safe_int_realloc(int* vect, int len, int new_len, int zero_val){
	assert(new_len>len);
	int i;
	int * new_vect = cast_vec_int(new_len);
	//vect = xrealloc(vect,sizeof(int)*new_len);
	for(i=0;i<len;i++)
	{
		new_vect[i]=vect[i];
	}
	for(i=len;i<new_len;i++)
	{
		new_vect[i]=zero_val;
	}	
	free(vect);
	return new_vect;
}


/********************************************************************************
 *  MATRIX HANDLING
 ********************************************************************************/
int** cast_mat_int(int rows, int cols){
	int** imat;
	int i,j;
	imat=(int**)malloc(sizeof(int*)*rows);
	for(i=0;i<rows;i++)imat[i]=(int*)malloc(sizeof(int)*cols);
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
		imat[i][j]=0;
		}
	}
	return imat;
}

double** cast_mat_double(int rows, int cols){
	double** dmat;
	int i,j;
	//printf("rows:%d,cols:%d\n",rows, cols); fflush(stdout);
	dmat=(double**)malloc(sizeof(double*)*rows);
	for(i=0;i<rows;i++)
	{
		dmat[i]=(double*)malloc(sizeof(double)*cols);
		for(j=0;j<cols;j++)
		{
		dmat[i][j]=0.;
		//printf("%d,%d\n",i,j);fflush(stdout);
		}
	}
	return dmat;
}

/********************************************************************************
 ********************************************************************************/

void free_mat_double(double **mat_d, int rows){
	int i;
	for(i=0;i<rows;i++)
	{
		free(mat_d[i]);
	}
	free(mat_d);
	return;
}

void free_mat_int(int **mat_d, int rows){
	int i;
	for(i=0;i<rows;i++)
	{
		free(mat_d[i]);
	}
	free(mat_d);
	return;
}
/********************************************************************************
 ********************************************************************************/
double * flatten_matrix_triangular_double(double** w, int rows, int * length){
	int i,j,aux;
	double * flat_w=cast_vec_double(rows*rows);
	aux=0;
	for(i=0;i<rows;i++)
	{
		for(j=0;j<i;j++)
		{
			flat_w[aux]=w[i][j];
			aux++;
		}
	}
	*length=aux;
	flat_w=realloc(flat_w,aux*sizeof(double)); // realloc memory
	return flat_w;
}


int * flatten_sparse_matrix_int(int** w, int rows, int cols, int zeros, int * length){
	int i,j,aux;
	int * flat_w=cast_vec_int(rows*cols);
	aux=0;
	if(zeros>0)
	{
		for(i=0;i<rows;i++)
		{
			for(j=0;j<cols;j++)
			{
				flat_w[aux]=w[i][j];
				aux++;
			}
		}
	}else{
		for(i=0;i<rows;i++)
		{
			for(j=0;j<cols;j++)
			{
				if(w[i][j]>0)
				{
					flat_w[aux]=w[i][j];
					aux++;
				}
			}
		}
		flat_w=realloc(flat_w,aux*sizeof(int)); // realloc memory
	}
	(*length)=aux;
return flat_w;
}

double * flatten_sparse_matrix_double(double** w, int rows, int cols, int zeros, int * length){
	int i,j,aux;
	double * flat_w=cast_vec_double(rows*cols);
	aux=0;
	if(zeros>0)
	{
		for(i=0;i<rows;i++)
		{
			for(j=0;j<cols;j++)
			{
				flat_w[aux]=w[i][j];
				aux++;
			}
		}
	}else{
		for(i=0;i<rows;i++)
		{
			for(j=0;j<cols;j++)
			{
				if(w[i][j]>0)
				{
					flat_w[aux]=w[i][j];
					aux++;
				}
			}
		}
		flat_w=realloc(flat_w,aux*sizeof(double)); // realloc memory
	}
	(*length)=aux;
return flat_w;
}


int * flatten_matrix_int(int** w, int rows, int cols, int diag, int * length){
	int i,j,aux;
	int * flat_w=cast_vec_int(rows*cols);
	aux=0;
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			if(diag>0 || i!=j)
			{
				flat_w[aux]=w[i][j];
				aux++;
			}
		}
	}
	flat_w=realloc(flat_w,aux*sizeof(int)); // realloc memory
	(*length)=aux;
return flat_w;
}

double * flatten_matrix_double(double** w, int rows, int cols, int diag, int * length){
	int i,j,aux,tot;
	double * flat_w;
	if(diag>0)
	{
		tot = rows*cols;
	}else{
		tot = rows*cols-maxeq_int(rows,cols);
	}
	flat_w = cast_vec_double(tot);
	aux=0;
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			if(diag>0 || i!=j)
			{
				flat_w[aux]=w[i][j];
				//printf("%.15lf\n",w[i][j]); fflush(stdout);
				aux++;
			}
		}
	}
	//flat_w=realloc(flat_w,aux*sizeof(double)); // realloc memory
	assert(aux==tot);
	(*length)=aux;
return flat_w;
}


/********************************************************************************
 ********************************************************************************/
double matrix_max_value_double(double** vect, int row, int col){
	int i,j;
	double max_v=0;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			max_v=fmax(vect[i][j],max_v);
		}
	}
	return max_v;
}
double matrix_min_value_double(double** vect, int row, int col){
	int i,j;
	double max_v=0;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			max_v=fmin(vect[i][j],max_v);
		}
	}
	return max_v;
}

int matrix_max_value_int(int** vect, int row, int col){
	int i,j;
	int max_v=0;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			max_v=maxeq_int(vect[i][j],max_v);
		}
	}
	return max_v;
}
int matrix_min_value_int(int** vect, int row, int col){
	int i,j;
	int max_v=0;
	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
		{
			max_v=mineq_int(vect[i][j],max_v);
		}
	}
	return max_v;
}
/********************************************************************************
 ********************************************************************************/

// this is legacy code //
void average_matrix(double** matrix, int rows, int cols, int reps){
	int i,j;
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			matrix[i][j]=matrix[i][j]/reps;
		}
	}
	return;
}
// up to here //
void scale_matrix(double** matrix, int rows, int cols, double reps){
	int i,j;
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			matrix[i][j]=matrix[i][j]*reps;
		}
	}
	return;
}

/********************************************************************************
 ********************************************************************************/
int sum_matrix_int(int** matrix, int rows, int cols){
	int sum=0;
	int i,j;
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			sum+=matrix[i][j];
		}
	}
	return sum;
}

double sum_matrix_double(double** matrix, int rows, int cols){
	double sum=0;
	int i,j;
	for(i=0;i<rows;i++)
	{
		for(j=0;j<cols;j++)
		{
			sum+=matrix[i][j];
		}
	}
	return sum;
}


/********************************************************************************
 *  VECTOR HANDLING
 ********************************************************************************/

int* cast_vec_int(int len){
	int* ivec;
	int i;
	ivec=(int*)malloc(sizeof(int)*len);
	for(i=0;i<len;i++)ivec[i]=0;
	return ivec;
}

double* cast_vec_double(int len){
	double* dvec;
	int i;
	dvec=(double*)malloc(sizeof(double)*len);
	for(i=0;i<len;i++)dvec[i]=0;
	return dvec;
}


/********************************************************************************
 ********************************************************************************/

double sum_vec_double(double* vec, int len){
	double su=0;
	int i;
	for(i=0;i<len;i++)
	{
		su+=vec[i];
	}
	return su;
}

int sum_vec_int(int* vec, int len){
	int su=0;
	int i;
	for(i=0;i<len;i++)
	{
		su+=vec[i];
	}
	return su;
}

/********************************************************************************
 ********************************************************************************/
int count_appearances_int(int *vec, int val, int len){
	int apps=0;
	int i;
	for(i=0;i<len;i++)
	{
		if(vec[i]==val)apps+=1;
	}
	return apps;
}


int count_appearances_double(double *vec, double val, double tol, int len){
	int apps=0;
	int i;
	for(i=0;i<len;i++)
	{
		if(fabs(vec[i]-val)<tol)apps+=1;
	}
	return apps;
}

int find_value_int(int * vec, int val, int len){
	int pos=-1;
	int i;
	for(i=0;i<len;i++)
	{
		if(vec[i]==val)
		{
			pos=i;
		}
	}
	return pos;
}
/********************************************************************************
 ********************************************************************************/

void reverse_vec_int(int* inver_a,int j)
{
   int i,temp;
   j--;
   for(i=0;i<(j/2);i++)
   {
      temp=inver_a[i];
      inver_a[i]=inver_a[j];
      inver_a[j]=temp;
      j--;
   }
}


void reverse_vec_double(double* inver_a,int j)
{
   int i;
   double temp;
   j--;
   for(i=0;i<(j/2);i++)
   {
      temp=inver_a[i];
      inver_a[i]=inver_a[j];
      inver_a[j]=temp;
      j--;
   }
}

/********************************************************************************
 ********************************************************************************/
double max_value_double(double* vect, int len){
	int i;
	double max_v=vect[0];
	for(i=0;i<len;i++)
	{
		max_v=fmax(vect[i],max_v);
	}
	return max_v;
}
double min_value_double(double* vect, int len){
	int i;
	double max_v=vect[0];
	for(i=0;i<len;i++)
	{
		max_v=fmin(vect[i],max_v);
	}
	return max_v;
}
int max_value_int(int* vect, int len){
	int i,max_v;
	max_v=vect[0];
	for(i=0;i<len;i++)
	{
		max_v=maxeq_int(vect[i],max_v);
	}
	return max_v;
}
int min_value_int(int* vect, int len){
	int i,max_v;
    max_v=vect[0];
	for(i=0;i<len;i++)
	{
		max_v=mineq_int(vect[i],max_v);
	}
	return max_v;
}

/********************************************************************************
 ********************************************************************************/
void scale_vec_int(int* vect, int factor, int len){
    int i;
    for(i=0;i<len;i++)
    {
        vect[i]=(int) (vect[i]*factor);
    }
    return;
}

void scale_vec_double(double* vect, double factor, int len){
    int i;
    for(i=0;i<len;i++)
    {
        vect[i]= vect[i]*factor;
    }
    return;
}


/********************************************************************************
 ********************************************************************************/
double* vec_int_to_double(int* vect, int len){
	 int i;
	 double* newvec=cast_vec_double(len);
	 for(i=0;i<len;i++)
	 {
		newvec[i]=(double)vect[i];
	 }
	 return newvec;
}

int* vec_double_to_int(double* vect, int len){
	 int i;
	 int * newvec=cast_vec_int(len);
	 for(i=0;i<len;i++)
	 {
		newvec[i]=(int)vect[i];
	 }
	 return newvec;
}


/********************************************************************************
 *  Single number HANDLING
 ********************************************************************************/
int maxeq_int(int a, int b){
	int c;
	if(a<b)
	{
		c=b;
	}else{
		c=a;
	}
	return c;
}
		
/********************************************************************************
 ********************************************************************************/
int mineq_int(int a, int b){
	int c;
	if(a<b)
	{
		c=a;
	}else{
		c=b;
	}
	return c;
}
