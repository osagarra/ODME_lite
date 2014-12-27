/************************************************************
 *
 *                    Stat Library
 *
 *		Statistical functions (histograms and so on)
 *
 *
 *************************************************************/
 
 #include "ula_stat_funcs.h"


/*****************************************************/
/************ Basic stats on vectors *************/ 
/*****************************************************/
double mean_vec_int(int * vec, int len){
    int i,check;
    int acc=0;
    check=0;
    double mean;
    for(i=0;i<len;i++)
    {
	acc+=vec[i];
	check++;
    }
    assert(check==len);
    mean=(double)acc/(double)len;
    return mean;
}
double mean_vec_double(double * vec, int len){
    int i,check;
    double acc=0;
    check=0;
    double mean;
    for(i=0;i<len;i++)
    {
	acc+=vec[i];
	check++;
    }
    assert(check==len);
    mean=(double)acc/(double)len;
    return mean;
}
double var_vec_int(int * vec, int len){
    int i,check;
    int acc,acc2;
    double mean,var;
    acc=acc2=0;
    check=var=0;
    for(i=0;i<len;i++)
    {
	acc+=vec[i];
	acc2+=vec[i]*vec[i];
	check++;
    }
    assert(check==len);
    mean=(double)acc/(double)len;
    var=(double)acc2/(double)len-mean*mean;
    return var;
}   
double var_vec_double(double * vec, int len){
    int i,check;
    double acc,acc2;
    double mean,var;
    acc=acc2=0;
    check=var=0;
    for(i=0;i<len;i++)
    {
	acc+=vec[i];
	acc2+=vec[i]*vec[i];
	check++;
    }
    assert(check==len);
    mean=acc/(double)len;
    var=acc2/(double)len-mean*mean;
    return var;
}



/*****************************************************/
/************ Histograms and cumulatives *************/ 
/*****************************************************/
int * histogram_int(int* vect, int minx, int maxx, int len ){
    // only for positive values!
    assert((maxx-minx)>0);assert(maxx>0);assert(minx>=0);
    int i,max_len=0;
    int * hist=cast_vec_int(maxx);
    for(i=0;i<len;i++)
    {
	assert(vect[i]>=0);
	if((vect[i]<=maxx) || (vect[i]>=minx))
	    {
		hist[vect[i]]++;
		max_len=maxeq_int(vect[i],max_len);
	    }
    }
    hist=realloc(hist,sizeof(int)*max_len);
    return hist;
}
/********************/
double * normalize_hist_int(int* hist, int len){
	int i;
	double* norm_hist=cast_vec_double(len);
	int norm=sum_vec_int(hist,len);
	double check=0;
	for(i=0;i<len;i++)
	{
		norm_hist[i]=(double)hist[i]/(double)norm;
		check+=norm_hist[i];
	}
	printf(" Normalizing hist ... norm=%.4f \n",check);
	return norm_hist;
}

void normalize_gsl_hist(gsl_histogram * hist){
    double N= 1./gsl_histogram_sum(hist);
    assert(N!=0);
    gsl_histogram_scale(hist, N);
    return;
}

void normalize_gsl_log_hist(gsl_histogram * hist){
	// not tested yet //
	int n= gsl_histogram_bins (hist);
	int c,i;
	double dx1,dx2;
	double norm = 0.;
	for(i=0;i<n;i++)
	{
		c=gsl_histogram_get (hist, i);	
		gsl_histogram_get_range (hist, i, &dx1, &dx2);
		norm += (double)c*(dx2-dx1);
	}
    gsl_histogram_scale(hist, 1./norm);
    return;	
}

/************ Using GSL **********/
/****** 1 D ******/
gsl_histogram * histogram_double(double* vect, double minx, double maxx, int nbins, int len){
    gsl_histogram * hist = gsl_histogram_alloc (nbins);
    gsl_histogram_set_ranges_uniform(hist, minx, maxx);
    int i;
    for(i=0;i<len;i++)
    {
        gsl_histogram_increment (hist, vect[i]);
    }
    return hist;
}
// This histogram does not work properly (needs additional normalization by bin-width)
gsl_histogram * histogram_double_log(double* vect, double minx, double maxx, double expo, int len){ // not working correctly!!!!!
    int i=0;
    int numbins=0;
    double * ranges=log_bins_double(minx, maxx, expo, &numbins);
    //double* aux=cast_vec_double(numbins-1);
    gsl_histogram * hist = gsl_histogram_alloc (numbins-1);
    /*double x=0;
    while(x<maxx)
    {
        printf("%lf %lf\n",x, ranges[i]);fflush(stdout);
        x+=ranges[i];
        i++;
    }*/
    /*for(i=0;i<numbins;i++)
    {
        printf("%d %lf\n",i,ranges[i]);fflush(stdout);
    }
    */
    gsl_histogram_set_ranges(hist, ranges, numbins);
    for(i=0;i<len;i++)
    {
        gsl_histogram_increment ( hist, vect[i]);        
    }
    // Dirty hack for log bins --> normalize bin //
    /*for(i=0;i<numbins-1;i++)
    {
        aux[i] = gsl_histogram_get (hist, i);
        aux[i] = aux[i]/(ranges[i+1]-ranges[i]);
    }
    
    gsl_histogram_reset (hist);
    for(i=0;i<numbins-1;i++)
    {
        gsl_histogram_accumulate (hist, aux[i]);
    }    
    free(aux);
    */
    free(ranges);
    return hist;
}
/****** 2 D ******/
double** y_of_x(double* vectx, double* vecty, double* xrange, int len, int len_ranges){
    double** yy=cast_mat_double(4,len_ranges);
    int i,j;
    for(i=0;i<len;i++)
    {
	for(j=0;j<len_ranges;j++)
	{
	    if(vectx[i]-xrange[j]<0)
	    {
		//printf("##### %lf %lf\n",vectx[i],xrange[j]);fflush(stdout);
	    	yy[1][j-1]+=vecty[i];
	    	yy[2][j-1]+=vecty[i]*vecty[i];
	    	yy[3][j-1]+=1;
		break;
	    }
	}
    }
    for(i=0;i<len_ranges-1;i++)
    {
    	yy[0][i]=xrange[i];
	//yy[0][i]=(yy[0][i+1]-yy[0][i])/2. + yy[0][i];
	yy[1][i]/=yy[3][i];
	yy[2][i]=sqrt(yy[2][i]/yy[3][i]-(yy[1][i]*yy[1][i]));
	//printf("%lf %lf %lf %lf\n",yy[0][i],yy[1][i],yy[2][i],yy[3][i]);fflush(stdout);
    }
    return yy;
}


gsl_histogram2d * histogram_2d_double(double* vectx, double* vecty, double minx, double maxx, double miny, double maxy, int xbins, int ybins, int lenx, int leny){
    gsl_histogram2d * hist = gsl_histogram2d_alloc (xbins,ybins);
    gsl_histogram2d_set_ranges_uniform(hist, minx, maxx, miny, maxy);
    int i,j;
    for(i=0;i<lenx;i++)
    {
    	for(j=0;j<leny;j++)
    	{
        gsl_histogram2d_increment (hist, vectx[i],vecty[j]);
    	}
    }
    return hist;
}

gsl_histogram2d * histogram_2d_double_log(double* vectx, double* vecty, double minx, double maxx, double miny, double maxy, double expx, double expy, int lenx, int leny){
    int j;
    int xbins=0;
    int ybins=0;
    double * xranges=log_bins_double(minx, maxx, expx, &xbins);
    double * yranges=log_bins_double(miny, maxy, expy, &ybins);
    gsl_histogram2d * hist = gsl_histogram2d_alloc (xbins-1,ybins-1);
    gsl_histogram2d_set_ranges( hist, xranges, xbins, yranges, ybins);
    free(xranges);
    free(yranges);
    for(j=0;j<leny;j++)
    {
    	gsl_histogram2d_increment (hist, vectx[j],vecty[j]);
	//printf("%lf,%lf\n",vectx[j],vecty[j]);fflush(stdout);
    }
    //printf("<x> %lf, <y> %lf, bin_max_x: %lf bin_max_y: %lf\n",gsl_histogram2d_xmean (hist),gsl_histogram2d_ymean (hist),gsl_histogram2d_xmax (hist), gsl_histogram2d_ymax (hist));fflush(stdout);
    return hist;
}

gsl_histogram2d * histogram_2d_double_linlog(double* vectx, double* vecty, double minx, double maxx, double miny, double maxy, int xbins, double expy, int lenx, int leny){
    int i,j;
    int ybins=0;
    double x,deltax;
    double * yranges=log_bins_double(miny, maxy, expy, &ybins);
    double* xranges=cast_vec_double(xbins+1);
    x=minx-1e-15;
    deltax=(maxx-minx+1e-15)/(double)(xbins);
    for(i=0;i<xbins+1;i++)
    {
	xranges[i]=x;
	x+=deltax;
    }
    //printf("%lf, max %lf, min %lf x %lf\n",xranges[xbins],maxx,minx,x);fflush(stdout);
    gsl_histogram2d * hist = gsl_histogram2d_alloc (xbins,ybins-1);
    gsl_histogram2d_set_ranges( hist, xranges, xbins+1, yranges, ybins);
    free(xranges);
    free(yranges);
    assert(lenx==leny);
    for(j=0;j<leny;j++)
    {
    	gsl_histogram2d_increment (hist, vectx[j],vecty[j]);
	//printf("%lf,%lf\n",vectx[j],vecty[j]);fflush(stdout);
    }
    //printf("<x> %lf, <y> %lf, bin_max_x: %lf bin_max_y: %lf\n",gsl_histogram2d_xmean (hist),gsl_histogram2d_ymean (hist),gsl_histogram2d_xmax (hist), gsl_histogram2d_ymax (hist));fflush(stdout);
    return hist;
}

void histogram_2d_mean(gsl_histogram2d * hist, int axis, double * mean_v, double * std_v, double* xrange1, double* xrange2, int to_int){// axis=0 for <> over column, 1 for rows
    int lenx=gsl_histogram2d_nx(hist);
    int leny=gsl_histogram2d_ny(hist);
    double xmin,xmax,deltax,norm;
    int i,j;
    double mean,std;
    assert(axis==1 || axis==0);
    if(axis==0) // col average
    {
    	for(i=0;i<lenx;i++)
    	{
	    mean=std=norm=0;
	    for(j=0;j<leny;j++)
	    {
		gsl_histogram2d_get_yrange (hist, j, &xmin, &xmax);
		if(to_int>0)
		{
		    deltax=floor((xmax-xmin)/2.+xmin);
		}else{
		    deltax=(xmax-xmin)/2.+xmin;
		}
		//deltax=xmax;
		norm+=gsl_histogram2d_get(hist, i, j);
		mean+=gsl_histogram2d_get(hist, i, j)*deltax;
		//printf(" i:%d j:%d x: %lf m:%lf\n",i,j,mean,gsl_histogram2d_get(hist, i, j));fflush(stdout);
	    	std+=gsl_histogram2d_get(hist, i, j)*deltax*deltax;
	    	gsl_histogram2d_get_xrange (hist, i, &xrange1[i], &xrange2[i]);
	    }
	    mean/=(double)norm;
	    std=sqrt(std/(double)norm-mean*mean);
	    mean_v[i]=mean;
	    std_v[i]=std;
    	}
    }else{
    	for(i=0;i<leny;i++)
    	{
	    mean=std=norm=0;
	    for(j=0;j<lenx;j++)
	    {
		gsl_histogram2d_get_yrange (hist, j, &xmin, &xmax);
		if(to_int>0)
		{
		    deltax=floor(xmax);
		}else{
		    deltax=(xmax-xmin)/2.+xmin;
		}
		//deltax=xmax-xmin;
		//norm+=deltax;
		norm+=gsl_histogram2d_get(hist, i, j);		
	    	mean+=gsl_histogram2d_get(hist, j, i)*deltax;
	    	std+=gsl_histogram2d_get(hist, j, i)*gsl_histogram2d_get(hist, j, i)*deltax*deltax;
		gsl_histogram2d_get_yrange (hist, i, &xrange1[i], &xrange2[i]);
	    }
	    mean/=(double)norm;
	    std=sqrt(std/(double)norm-mean*mean);
	    mean_v[i]=mean;
	    std_v[i]=std;
    	}
    }
    return;
}
/*****************************************************/
/************ Bins *************/
/*****************************************************/
double * log_bins_double(double minx, double maxx, double expo, int * numbins){
    assert(minx<maxx);
    int i,check;
    int max_ind=1000;
    double pos;
    double * bins= cast_vec_double(max_ind);
    double deltax=1.;
    pos=minx;
    check=i=0;
    if(minx==0)
    {
	bins[0]=0;
    	pos=1.0000001;
    	i=1;
	check++;
    }
    do
    {
	bins[i]=pos;
	deltax*=expo;
	pos+=(deltax*expo);
	i++;
	check++;
	if(check>max_ind)
	{
	    bins=realloc(bins,sizeof(double)*check*2); // doblem espai 
	}
    	//printf("%lf,%lf %d\n",pos,deltax,check);fflush(stdout);
	//printf("%lf\n",bins [i-1]);fflush(stdout);
	}while(pos<maxx);
    //bins[i-1]=maxx+maxx*0.1*(maxx-bins[i-2]);
    bins[i-1]=maxx+1;
    //printf("%lf %lf\n",bins[i-1],maxx);fflush(stdout);
    bins=realloc(bins,sizeof(double)*check);
    //printf("# Done %d bins\n",check);fflush(stdout);
    //printf("# Bin 0: %lf\n",bins[0]);fflush(stdout);
    //printf("# Bin 1: %lf\n",bins[1]);fflush(stdout);
    *numbins=check;
    return bins;
}
/********************/
/*****************************************************/
/************ Accumulators *************/
/*****************************************************/
/************** 1 D **************************/
// Sets and updates accumulators to aggregate during simulation --> valid for all
gsl_histogram * set_acc_int(int minv, int maxv){
	double* dummy=cast_vec_double(2);
	dummy[0]=dummy[1]=0;
	if(minv>0) minv-=1;
	gsl_histogram *  acc_k=histogram_double(dummy, (double)minv, (double)(maxv), (double)maxv+1, 2);
	gsl_histogram_reset (acc_k);
	return acc_k;
}
gsl_histogram * set_acc_double(double minv, double maxv, int bins){
	double* dummy=cast_vec_double(2);
	dummy[0]=dummy[1]=0;
	gsl_histogram *  acc_k=histogram_double(dummy, minv-0.000001, maxv, bins , 2);
	gsl_histogram_reset (acc_k);
	return acc_k;
}
// valid for w,k,s
void update_int_acc(int * vect, int len, gsl_histogram * acc_k, gsl_histogram * acc_k2, int norm){
	double * k=vec_int_to_double(vect,len);
	update_double_acc(k, len, acc_k, acc_k2, norm);
	free(k);
	return ;
}

void update_double_acc(double * k, int len, gsl_histogram * acc_k, gsl_histogram * acc_k2, int norm){
	gsl_histogram *  k_hist= histogram_double(k,gsl_histogram_min(acc_k), gsl_histogram_max(acc_k),(int)gsl_histogram_bins(acc_k) ,len);
	if(norm>0)
	{
	    normalize_gsl_hist(k_hist);
	}
	assert(gsl_histogram_equal_bins_p (acc_k, k_hist) ==1 );
	gsl_histogram_add (acc_k, k_hist);
	gsl_histogram_mul(k_hist,k_hist);
	gsl_histogram_add (acc_k2, k_hist);
	gsl_histogram_free(k_hist);
	return ;
}


/*****************************************************/
void acc_compute_std(gsl_histogram ** accs,int len, double norm){
    gsl_histogram * cl;
    int i;
    for(i=0;i<len;i+=2)
    {
    	cl= gsl_histogram_clone (accs[i]);
    	gsl_histogram_mul(cl,cl);
    	gsl_histogram_scale(cl, norm);
    	gsl_histogram_sub(accs[i+1],cl);
    	gsl_histogram_free(cl);
    }
    return;
}

void acc_normalize(gsl_histogram ** accs,int len, double norm){
    int i;
    for(i=0;i<len;i++)
    {
	//normalize_gsl_hist(accs[i]);
	gsl_histogram_scale (accs[i], norm);
    }
    return;
}
/*****************************************************/
void acc_allocate_all(gsl_histogram ** accs, double ** d, int ** s, int N_nodes){
    int maxs=maxeq_int(max_value_int(s[0],N_nodes),max_value_int(s[1],N_nodes));
    int max_length_d;
    double * d_flat=flatten_matrix_triangular_double(d, N_nodes, &max_length_d);
    double maxd=ceil(max_value_double(d_flat,max_length_d));
    free(d_flat);
    // histogram allocation
    int i;
    for(i=0;i<6;i++) // distances
    {
    	accs[i]=set_acc_double( 0,  maxd, (int)(N_nodes/10));
    }
    for(i=6;i<12;i++) // w,kout,kin
    {
    	accs[i]=set_acc_int( 0,  maxs); // way over
    }
    for(i=12;i<16;i++) // sout,sin
    {
	accs[i]=set_acc_int( 0, maxs+1000);
    }
    return;
}

void acc_free_all(gsl_histogram ** accs, int len){
    int i;
    for(i=0;i<len;i++)
    {
	gsl_histogram_free(accs[i]);
    }
    return;
}
/************** 2 D **************************/

gsl_histogram2d * set_acc2d_linlog(double minx, double maxx, double miny, double maxy, int binsx, double exp_y){
    double* dummy1=cast_vec_double(2);
    double* dummy2=cast_vec_double(2);
    dummy1[0]=dummy1[1]=0.5;
    dummy2[0]=dummy2[1]=0.5;
    gsl_histogram2d* w_r=histogram_2d_double_linlog(dummy1, dummy2, minx , maxx, miny, maxy, binsx, exp_y, 2, 2);
    gsl_histogram2d_reset(w_r);
    free(dummy1);
    free(dummy2);
    return w_r;
}

gsl_histogram2d * set_acc2d_log(double minx, double maxx, double miny, double maxy, double exp_x, double exp_y){
    double* dummy1=cast_vec_double(2);
    double* dummy2=cast_vec_double(2);
    dummy1[0]=dummy1[1]=0.5;
    dummy2[0]=dummy2[1]=0.5;
    gsl_histogram2d* w_r=histogram_2d_double_log(dummy1,dummy2,minx,maxx,miny,maxy,exp_x,exp_y, 2,2);
    gsl_histogram2d_reset(w_r);
    free(dummy1);
    free(dummy2);
    return w_r;
}

/*****************************************************/

void acc2d_compute_std(gsl_histogram2d ** accs,int len, double norm){
    gsl_histogram2d * cl;
    int i;
    for(i=0;i<len;i+=2)
    {
    	cl= gsl_histogram2d_clone (accs[i]);
    	gsl_histogram2d_mul(cl,cl);
    	gsl_histogram2d_scale(cl, norm);
    	gsl_histogram2d_sub(accs[i+1],cl);
    	gsl_histogram2d_free(cl);
    }
    return;
}

void acc2d_normalize(gsl_histogram2d ** accs,int len, double norm){
    int i;
    for(i=0;i<len;i++)
    {
	//normalize_gsl_hist(accs[i]);
	gsl_histogram2d_scale (accs[i], norm);
    }
    return;
}




/*****************************************************/
void acc2d_allocate_all(gsl_histogram2d ** accs, double ** d, int ** s, int N_nodes, double exp_s, double exp_t, int binsx){
    int maxs=maxeq_int(max_value_int(s[0],N_nodes),max_value_int(s[1],N_nodes));
    int t_trips=sum_vec_int(s[0],N_nodes);
    int max_length_d;
    double * d_flat=flatten_matrix_triangular_double(d, N_nodes, &max_length_d);
    double maxd=ceil(max_value_double(d_flat,max_length_d));
    free(d_flat);
    // histogram allocation
    accs[0]=set_acc2d_linlog( 0,  maxd, 0, maxs, binsx, exp_t); // w_r
    accs[1]=set_acc2d_linlog( 0,  maxd, 0, maxs, binsx, exp_t); // w_r
    accs[2]=set_acc2d_log( 0,  maxs, 0, maxs, exp_s, exp_t); // w_s_out
    accs[3]=set_acc2d_log( 0,  maxs, 0, maxs, exp_s, exp_t); // w_s_out
    accs[4]=set_acc2d_log( 0,  maxs, 0, maxs, exp_s, exp_t); // w_s_in
    accs[5]=set_acc2d_log( 0,  maxs, 0, maxs, exp_s, exp_t); // w_s_in
    accs[6]=set_acc2d_log( 0,  maxs*maxs, 0, maxs, exp_s, exp_t); // w_s_in_out
    accs[7]=set_acc2d_log( 0,  maxs*maxs, 0, maxs, exp_s, exp_t); // w_s_in_out
    accs[8]=set_acc2d_log( 0,  t_trips, 0, maxs, exp_s, exp_t); // w_s_ij
    accs[9]=set_acc2d_log( 0,  t_trips, 0, maxs, exp_s, exp_t); // w_s_ij
    return;
}

void acc2d_free_all(gsl_histogram2d ** accs, int len){
    int i;
    for(i=0;i<len;i++)
    {
	gsl_histogram2d_free(accs[i]);
    }
    return;
}

/*****************************************************/
/************ Print Hists *************/
/*****************************************************/
void print_hist_double(char *input_name, int len, double *av_hist, double *av_hist2){
    // includes first and second moment over repetitions of integer or real values (already normalized)
	int m;
	double norm;
	FILE* input=open_file("w",input_name);
	//printf("... Printing Hist list...\n");
	norm=0;
	for(m=0;m<len;m++)
	{
		fprintf(input, "%d %.4lf %.4lf\n",m,av_hist[m],av_hist2[m]);
		norm+=av_hist[m];
	}
	fflush(input);
	//printf("Norm: %.8lf\n",norm);
	fclose(input);
	return;
}
/************************************/
void print_hist_int(char *input_name, int len, int *av_hist){
	// Prints normalized hist ! //
	int m,norm;
	double check;
	FILE* input=open_file("w",input_name);
	//printf("... Printing Hist list...\n");
	norm=check=0;
	for(m=0;m<len;m++)
	{
	    //printf("%d %d\n",m,av_hist[m]);fflush(stdout);
	    norm+=av_hist[m];
	}
	fprintf(input,"### Bin, Num_events/Cumulative #### ");
	fprintf(input,"### Norm Factor: %d\n",norm);
	for(m=0;m<len;m++)
	{
	    fprintf(input, "%d %.4lf\n",m,(double)av_hist[m]/(double)norm);
	    check+=(double)av_hist[m]/(double)norm;
	}	
	fflush(input);
	//printf("Norm: %.8lf Factor:%d\n",check,norm);
	fclose(input);
	return;
}
/*****************************************************/
void print_hist2d_mean(char *input_name, double * h_mean, double * h_std, double * xrange, int len){
    FILE* input=open_file("w",input_name);
    fprintf(input,"# X value, <Y>(x), sigma #\n");
    int i;
    for(i=0;i<len;i++)
    {
	fprintf(input,"%lf %lf %lf\n",xrange[i]+(xrange[i+1]-xrange[i])/2.,h_mean[i],h_std[i]);
    }
    fclose(input);
    return;
} 

/*****************************************************/
void print_acc(char *input_name, gsl_histogram * acc1, gsl_histogram * acc2){
	//acc2 corresponds to variances!
	FILE* input=open_file("w",input_name);
	//printf("... Printing Hist list...\n");
	fprintf(input,"# Mean:%.4lf Std:%.4lf Norm:%.4lf Max:%.4lf Min:%.4lf Max_bin:%.4lf Min_bin:%.4lf\n",(double)gsl_histogram_mean(acc1),(double)gsl_histogram_sigma(acc1),(double)gsl_histogram_sum(acc1),(double)gsl_histogram_max_val(acc1),(double)gsl_histogram_min_val(acc1),(double) gsl_histogram_max_bin(acc1), (double)gsl_histogram_min_bin(acc1));
	fprintf(input,"#Bin_id Bin_min Bin_max Bin_val Bin_std CCDF\n");
	int m;
	int len=gsl_histogram_bins (acc1);
	double bmax,bmin,bin_val;
	double CDF=0;
	double norm=0;
	for(m=0;m<len;m++)
	{
	    gsl_histogram_get_range (acc1, m, &bmin, &bmax);
	    bin_val=gsl_histogram_get (acc1, m);
	    norm+=bin_val*(bmax-bmin);
	}
	for(m=0;m<len;m++)
	{
	    gsl_histogram_get_range (acc1, m, &bmin, &bmax);
	    bin_val=gsl_histogram_get (acc1, m);
	    if(bin_val>1e-15)
	    {
	    	fprintf(input, "%d %lf %lf %lf %lf %lf\n",m, bmin, bmax, bin_val ,sqrt( gsl_histogram_get(acc2,m)+1e-15), 1.-CDF/norm);
	    }
	    CDF+=bin_val*(bmax-bmin);
	}
	fclose(input);
	return;
}

void print_acc2d(char *input_name, gsl_histogram2d * acc1, gsl_histogram2d * acc2){
    //acc2 corresponds to variances!
    FILE* input=open_file("w",input_name);
    //printf("... Printing Hist list...\n");
    fprintf(input,"#x Mean:%.4lf x Std:%.4lf x Norm:%.4lf x Max:%.4lf x Min:%.4lf\n",(double)gsl_histogram2d_xmean(acc1),(double)gsl_histogram2d_xsigma(acc1),(double)gsl_histogram2d_sum(acc1),(double)gsl_histogram2d_max_val(acc1),(double)gsl_histogram2d_min_val(acc1));
    fprintf(input,"#y Mean:%.4lf y Std:%.4lf y \n",(double)gsl_histogram2d_ymean(acc1),(double)gsl_histogram2d_ysigma(acc1));
    fprintf(input,"#Bin_idx Bin_idy Bin_minx Bin_maxx Bin_miny Bin_maxy Bin_val Bin_std\n");
    int m,n;
    int lenx=gsl_histogram2d_nx(acc1);
    int leny=gsl_histogram2d_ny(acc1);
    double bmaxx,bminx,bin_val;
    double bmaxy,bminy;
    //double CDF=0;
    //double norm=0;
    for(m=0;m<lenx;m++)
    {
    	for(n=0;n<leny;n++)
    	{
	    gsl_histogram2d_get_xrange (acc1, m, &bminx, &bmaxx);
	    gsl_histogram2d_get_yrange (acc1, n, &bminy, &bmaxy);
	    bin_val=gsl_histogram2d_get (acc1, m, n);
	    //CDF+=bin_val*(bmax-bmin);
	    if(bin_val>1e-15)
	    {
	    	fprintf(input, "%d %d %lf %lf %lf %lf %lf %lf\n",m, n,  bminx, bmaxx, bminy, bmaxy, bin_val ,sqrt(gsl_histogram2d_get(acc2,m,n)+1e-15)); //1.-CDF/norm
	    }
    	}
	fprintf(input,"\n");
    }
    fclose(input);
    return;
}
