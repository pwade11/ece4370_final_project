#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <string.h>

#define BILLION 1000000000L
//#define CHUNK_SIZE 32768
//#define CHUNK_SIZE 8192
#define CHUNK_SIZE 16384



//gcc -O3 -o project_sim_2d_c project_sim_2d.c

int main(void){
	

	int SecondWG = 0;


	struct timespec start_gen, start_processing, end_processing, end;
	
	clock_gettime(CLOCK_MONOTONIC, &start_gen);

	long double lambd = 1.5*1E-6;
	long double k0 = 2*M_PI/lambd;

	long double complex nCladding = 1;
	long double complex nCore = 1.1;

	long double n_bar = (nCladding + nCore)/2;

	long double inputWGLength = 100*1E-6;
	long double inputWGWidth = 1.5*1E-6;
	long double widthDomain = 15*1E-6;
	long double lengthDomain = inputWGLength;

	long double coupGap = 0.5*1E-6;
	long double sig = 0.5*1.2*1E-6;


	//long double del_z = 0.25*1E-6;
	long double del_z = 0.01*1E-6;

	int Nzpts = round(lengthDomain/del_z);

	long double del_x = 0.01*1E-6;
	//long double del_x_2 = del_x*del_x;
	int N = round(widthDomain / del_x);	
	


	long double alpha = 0.5;

	//all the 1d arrays
	long double complex * u = (long double complex *) malloc (N * sizeof(long double complex));
	long double complex * b = (long double complex *) malloc (N * sizeof(long double complex));
	long double complex * gamma = (long double complex *) malloc (N * sizeof(long double complex));
	long double complex * r = (long double complex *) malloc (N * sizeof(long double complex));
	long double * x = (long double *) malloc (N * sizeof(long double));
	long double complex * kappa = (long double complex *) malloc (N * sizeof(long double complex));

	int min_core_ind = N;
	int max_core_ind = 0;
	int min_core_ind_2 = N;
	int max_core_ind_2 = 0;
	//linspace as well as a lot of allocation

	long double widthAbsEdge = 3*1E-6;
	long double kappa_max = -0.2;
		
	
	long double step = widthDomain/((long double)(N-1));
	for(int i = 0; i <N; i++){
		x[i] = -1*widthDomain/2 + ((long double)i*step);
		u[i] = expl(-1*(x[i]/sig)*(x[i]/sig));
		
		if((x[i] <= inputWGWidth/2) && (x[i] >= -1*inputWGWidth/2)){
			//THIS ONLY DOES THE SINGLE WAVEGUIDE FOR NOW
			if(i < min_core_ind){
				min_core_ind = i;
			}
			if(i > max_core_ind){
				max_core_ind = i;
			}
		}
		if((x[i] <= (3*inputWGWidth/2 + coupGap)) && (x[i]>=(inputWGWidth/2 + coupGap))){
			if(i < min_core_ind_2){
				min_core_ind_2 = i;
			}
			if(i > max_core_ind_2){
				max_core_ind_2 = i;
			}
		}
		if(x[i] < -1*widthDomain/2 + widthAbsEdge){
			kappa[i] = ((widthDomain/2 - widthAbsEdge + x[i])/widthAbsEdge);
			kappa[i] = kappa_max*kappa[i]*kappa[i]; //square self and mult kappa max
		}
		else if(x[i] > widthDomain/2 - widthAbsEdge){
			kappa[i] = (x[i] - widthDomain/2 + widthAbsEdge)/widthAbsEdge;
			kappa[i] = kappa_max*kappa[i]*kappa[i];
		}
		else{
			kappa[i] = 0;
		}

		
	}

	//printf("%d, %d, \n", max_core_ind, min_core_ind);

	//all the 2d arrays
	long double complex ** allFeild = (long double complex **) malloc((Nzpts + 1)*sizeof(long double complex *));
	long double complex ** n = (long double complex **) malloc((Nzpts + 1)*sizeof(long double complex *));
	

	for (int m = 0; m < Nzpts+1; m++){
		allFeild[m] = (long double complex *) malloc(N * sizeof(long double complex));
		n[m] = (long double complex *) malloc(N * sizeof(long double complex));
		for (int j =0; j < N; j++){
			if(j <= max_core_ind && j >= min_core_ind){
				n[m][j] = nCore;
			}
			else{
				if(SecondWG == 1 && j <= max_core_ind_2 && j >= min_core_ind_2){
					n[m][j] = nCore;
				}
				else{
					n[m][j] = nCladding;
				}
			}
			if(m ==0){
				allFeild[m][j] = u[j];
			}
			else{
				allFeild[m][j] = 0;
			}



		}
	}


	

	clock_gettime(CLOCK_MONOTONIC, &start_processing);

	//memoisation with the k0*k0*(n[m+1][0]*n[m+1][0] - kappa[0]*kappa[0] + 2*I*n[m+1][0]*kappa[0] - n_bar*n_bar) term
	long double complex * next_nmj_sqr_term = (long double complex *) malloc (N * sizeof(long double complex));



	long double complex a = -alpha/(del_x*del_x);


	next_nmj_sqr_term[0] = k0*k0*(n[0][0]*n[0][0] - kappa[0]*kappa[0] + 2*I*n[0][0]*kappa[0] - n_bar*n_bar);
	for (int j = 1; j < N-1; j++){
		next_nmj_sqr_term[j] = k0*k0*(n[0][j]*n[0][j] - kappa[j]*kappa[j] + 2*I*n[0][j]*kappa[j] - n_bar*n_bar);
	}
	next_nmj_sqr_term[N-1] = k0*k0*(n[0][N-1]*n[0][N-1] - kappa[N-1]*kappa[N-1] + 2*I*n[0][N-1]*kappa[N-1] - n_bar*n_bar);

	long double complex shared_r_b_term = 2*I*k0*n_bar/del_z;
	long double complex later_part_of_r = -1*(2*(1-alpha)/(del_x*del_x)) + shared_r_b_term;
	for (int m = 0; m < Nzpts; m ++){


		r[0] = ((1-alpha)/(del_x*del_x)) * (0+ u[1]) + ((1 - alpha)*next_nmj_sqr_term[0] + later_part_of_r)*u[0];
		next_nmj_sqr_term[0] = k0*k0*(n[m+1][0]*n[m+1][0] - kappa[0]*kappa[0] + 2*I*n[m+1][0]*kappa[0] - n_bar*n_bar);
		b[0] = (2*alpha/(del_x*del_x)) - alpha*next_nmj_sqr_term[0] + shared_r_b_term;

		

		for (int j = 1; j < N-1; j++){
			
			r[j] = ((1-alpha)/(del_x*del_x))*(u[j-1] + u[j+1]) + ((1 - alpha)*next_nmj_sqr_term[j] + later_part_of_r)*u[j];
			next_nmj_sqr_term[j] = k0*k0*(n[m+1][j]*n[m+1][j] - kappa[j]*kappa[j] + 2*I*n[m+1][j]*kappa[j] - n_bar*n_bar);
			b[j] = (2*alpha/(del_x*del_x)) - alpha*next_nmj_sqr_term[j] + shared_r_b_term;
			
		}

		r[N-1] = ((1-alpha)/(del_x*del_x))*(u[N-2] + 0) + ((1 - alpha)*next_nmj_sqr_term[N-1] + later_part_of_r)*u[N-1];
		next_nmj_sqr_term[N-1] = k0*k0*(n[m+1][N-1]*n[m+1][N-1] - kappa[N-1]*kappa[N-1] + 2*I*n[m+1][N-1]*kappa[N-1] - n_bar*n_bar);
		b[N-1] = (2*alpha/(del_x*del_x)) - alpha*next_nmj_sqr_term[N-1] + shared_r_b_term;
		
		/*
		long double complex beta = b[0];
		u[0] = r[0] / beta;
		for (int j =1; j<N; j++){
			gamma[j] = a/beta;
			beta = b[j] - a*gamma[j];
			u[j] = (r[j] - a*u[j-1])/beta;
		}

		for (int j =0; j <N-1; j++){
			int ktemp = N-j-2;
			u[ktemp] = u[ktemp] - gamma[ktemp+1]*u[ktemp+1];
			
		}
		*/
		//try in form the first paper gave
		//this way is substancially faster
		u[0] = r[0] / b[0]; //equivalent to d1 <-- d1/b1 in literature I think
		gamma[0] = a/b[0]; //eqivalent to their c1 part
		long double complex temp_val; //r in the paper 
		for (int j =1; j<N; j++){
			temp_val = 1/(b[j] - a*gamma[j-1]);
			u[j] = temp_val*(r[j] - a*u[j-1]);
			gamma[j] = temp_val*a;

		}

		for (int j =0; j <N-1; j++){
			int ktemp = N-j-2;
			u[ktemp] = u[ktemp] - gamma[ktemp]*u[ktemp+1];
			//allFeild[m+1][ktemp] = u[ktemp];
			//printf("%d\n",ktemp);
			
		}
		for (int j =0; j < N; j++){
			allFeild[m+1][j] = u[j];
		}
	}

	free(next_nmj_sqr_term);

	clock_gettime(CLOCK_MONOTONIC, &end_processing);
	
	

	//output
	/*
	FILE * fp;
	fp = fopen("project_sim_output.csv","w+");


	for (int m = 0; m < Nzpts + 1; m ++){
		for (int j =0; j < N; j++){
			//fprintf(fp, "%Lf",cabsl(allFeild[m][j])*cabsl(allFeild[m][j]));
			fprintf(fp, "%.*Le",DECIMAL_DIG, cabsl(allFeild[m][j])*cabsl(allFeild[m][j]));
			if(j != N-1){
				fprintf(fp,	",");
			}
		}
		if(m != Nzpts){
			fprintf(fp,	"\n");
		}
	}
	fclose(fp);
	*/

	//changing writes based on this:
	//https://stackoverflow.com/questions/41210227/fastest-way-to-write-integer-to-file-in-c
	
	FILE * fp;
	fp = fopen("project_sim_output.csv","w");

	char file_buffer[CHUNK_SIZE + 64] ;
	int buffer_count = 0 ;
	//int num_write = 0;

	for (int m = 0; m < Nzpts + 1; m ++){
		for (int j =0; j < N; j++){
			//buffer_count += sprintf( &file_buffer[buffer_count], "%.*Le", DECIMAL_DIG, (long double) (allFeild[m][j]*conjl(allFeild[m][j])));
			buffer_count += sprintf( &file_buffer[buffer_count], "%0.14Lf", (long double) (allFeild[m][j]*conjl(allFeild[m][j])));
			//the commented out one above specifies printing at full precision, but it its not necessary and this is much faster
			//printf("%d\n",buffer_count);
			

			if(j != N-1){
				//fprintf(fp,	",");
				buffer_count += sprintf( &file_buffer[buffer_count], ",");
			}

			// if the chunk is big enough, write it.
			if( buffer_count >= CHUNK_SIZE )
			{
				fwrite( file_buffer, CHUNK_SIZE, 1, fp ) ;
				buffer_count -= CHUNK_SIZE ;
				memcpy( file_buffer, &file_buffer[CHUNK_SIZE], buffer_count ) ;
				//num_write ++;
			}
		}
		if(m != Nzpts){
			buffer_count += sprintf( &file_buffer[buffer_count], "\n");
		}
		if( buffer_count >= CHUNK_SIZE )
		{
			fwrite( file_buffer, CHUNK_SIZE, 1, fp ) ;
			buffer_count -= CHUNK_SIZE ;
			memcpy( file_buffer, &file_buffer[CHUNK_SIZE], buffer_count ) ;
			//num_write ++;
		}
	}

	// Write remainder
	if( buffer_count > 0 )
	{
		fwrite( file_buffer, 1, buffer_count, fp ) ; 
		//num_write++;   
	}
	
	fclose(fp);
	






	clock_gettime(CLOCK_MONOTONIC, &end);

	
	double time_from_init = BILLION *(end.tv_sec - start_gen.tv_sec) +(end.tv_nsec - start_gen.tv_nsec);
	time_from_init = time_from_init / BILLION;

	printf("Elapsed time from init: %lf seconds\n", time_from_init);

	double time_from_proc = BILLION *(end.tv_sec - start_processing.tv_sec) +(end.tv_nsec - start_processing.tv_nsec);
	time_from_proc = time_from_proc / BILLION;

	printf("Elapsed time from processing: %lf seconds\n", time_from_proc);
	


	double time_proc = BILLION *(end_processing.tv_sec - start_processing.tv_sec) +(end_processing.tv_nsec - start_processing.tv_nsec);
	time_proc = time_proc / BILLION;
	printf("Elapsed time processing: %lf seconds\n", time_proc);
	//printf("%lf\n", time_proc);

	double time_saving = BILLION *(end.tv_sec - end_processing.tv_sec) +(end.tv_nsec - end_processing.tv_nsec);
	time_saving = time_saving / BILLION;
	printf("Elapsed time saving: %lf seconds\n", time_saving);



	//printf("num write = %d\n",num_write);


	for (int m = 0; m < Nzpts + 1; m ++){
		free(allFeild[m]);
		free(n[m]);
	}

	free(allFeild);
	free(n);
	free(r);
	free(x);
	free(kappa);
	free(b);
	free(gamma);
	free(u);


	//printf("%lu\n", sizeof(long double complex));
	//printf("%d\n", N);
	//printf("%d\n", Nzpts+1);


}








