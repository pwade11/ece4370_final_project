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

void generate_n_grid(long double complex nCore, long double complex ** n, int N, int Nzpts, 
	long double * x, long double * z, long double w_B, long double w, long double s, long double R, long double angle, 
	long double taper_start, long double taper_end, long double turn_start, long double tip_width){

	//assumes that n is aready created to be all nCladding before this function, since that can be done during initialization
	//this likely could be substancially optimized, becasue doing a lot more comparisons than required, but just starting this way
	//for now, since dont think will be a substancial time contributer

	int bus_start_ind = N;
	int bus_end_ind = 0;
	for(int j = 0; j<N; j++){
		if((x[j] >= -1*w_B/2) && j <bus_start_ind)
			bus_start_ind = j;
		if((x[j] <= w_B/2) && j > bus_end_ind)
			bus_end_ind = j;
	}

	long double z0 = turn_start + (R+w)*sinl(angle);
	long double x_start_angle = w_B/2 + s + R + w - cosl(angle)*(R+w);
	for (int m = 0; m < Nzpts+1; m++){

		for (int j =bus_start_ind; j <= bus_end_ind; j++){
			n[m][j] = nCore;
			//the bus waveguide
		}

		//know that the other waveguide stuff will start after the bus, and nothing above taper start
		
		if(z[m] >= taper_start){
			for(int j = bus_end_ind; j < N; j++){
				if(z[m] <= taper_end){
					//taper part
					if((x[j] >= w_B/2 + s) && (x[j] <= w_B/2 + s + tip_width))
						n[m][j] = nCore;
					if((x[j] >= w_B/2 + s + tip_width) && (x[j] <= w_B/2 + s + tip_width + (z[m] - taper_start)*(w - tip_width)/(taper_end - taper_start)))
						n[m][j] = nCore;
				}
				else if(z[m] <= turn_start){
					//straight part
					if((x[j] >= w_B/2 + s) && (x[j] <= w_B/2 + s + w))
						n[m][j] = nCore;
				}
				else if(z[m] <= turn_start + (R + w)*sinl(angle)){
					//turn part
					if((x[j] >= w_B/2 + s + R + w - sqrtl((R + w)*(R + w) - (z[m] - turn_start)*(z[m] - turn_start))) && (x[j] <= w_B/2 + s + R + w - sqrtl(R*R - (z[m] - turn_start)*(z[m] - turn_start))))
						n[m][j] = nCore;
				}
				else{
					//diagonal part
					if((x[j] >= x_start_angle + tanl(angle)*(z[m] - z0)) && (x[j] <= x_start_angle + w + tanl(angle)*(z[m] - z0)))
						n[m][j] = nCore;
				}
			}
		}
				
	}
	return;

}

void save_data(int argc, char **argv, long double complex ** allFeild, int Nzpts, int N){
	//changing writes based on this:
	//https://stackoverflow.com/questions/41210227/fastest-way-to-write-integer-to-file-in-c
	
	FILE * fp;
	if(argc == 3)
		fp = fopen(argv[2],"w");
	else
		fp = fopen("project_sim_output.csv","w");
	if(fp == NULL){
		printf("Could not open the output file\nAborting\n");
		exit(0);
	}

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
}


int main(int argc, char **argv){
	
	if((argc != 1) && (argc != 2) && (argc != 3)){
		printf("Usage project_sim_2d_c <input file (optional)> <output file(optional)>\n Aborting...\n");
		exit(0);
	}

	int save_n = 0; //if you want to dump the n data, only use this to check the file, will not optimize that saving 

	struct timespec start_gen, start_processing, end_processing, end;
	
	clock_gettime(CLOCK_MONOTONIC, &start_gen);
	


	
	
	//All of the parameters that you may want to adjust, should adjust them in the file

	long double complex nCladding;	//cladding index
	long double complex nCore;		//core (effective) index
	long double lambd;				//wavelength (um)
	long double widthDomain;		//width of total thing simulating (um)
	long double lengthDomain; 		//length of total thing simulating (um)
	long double w_B;				//width of the bus waveguide (um)
	long double s;					//separation between the two waveguides (um)
	long double w;					//width of the output waveguide (um)
	long double taper_start;		//top of the taper for the output waveguide (um)
	long double taper_end;			//bottom point of the taper for the output waveguide (um)
	long double turn_start;			//startpoint for where the output waveguide turns
	long double R;					//radius of the curve (um)
	long double angle;				//angle the output waveguide curves to (degrees)
	long double tip_width;			//the width of the tip of the taper of the output waveguide (um)
	long double sig;				//width of the initial amplitude distribution (um)
	long double del_x;				//discritization in x (um)
	long double del_z;				//discritization in z (um)
	long double alpha;				//mix between the two methods of simulation

	if(argc == 1){
		nCladding = 1.4446; //silicon dioxide at 1500 nm https://refractiveindex.info/?shelf=main&book=SiO2&page=Malitson
		nCore = 1.5927508350083501; //n effective from the python code, based on n = 1.747 for aluminum oxide at 1500 nm https://refractiveindex.info/?shelf=main&book=Al2O3&page=Malitson-o
		lambd = 1.5*1E-6;
		widthDomain = 15*1E-6;
		lengthDomain = 100*1E-6; //inputWGLength;
		w_B = 0.22*1E-6;
		s = 0.5*1E-6;
		w = 0.3*1E-6;
		taper_start = lengthDomain/4;
		taper_end = lengthDomain/2;
		turn_start = 3*lengthDomain/4;
		R = 100*1E-6;
		angle = 10*M_PI/180;
		tip_width = 0.15*1E-6;
		sig = 0.4*1.2*1E-6;
		del_x = 0.01*1E-6;
		del_z = 0.1*1E-6;
		alpha = 0.5;	
	}
	else{
		//reading file based on https://stackoverflow.com/questions/3501338/c-read-file-line-by-line
		FILE * fp_in;
		char * line = NULL;
		size_t len = 0;
		ssize_t read;
		fp_in = fopen(argv[1], "r");
		if(fp_in == NULL){
			printf("Could not read the input file\nAborting\n");
			exit(0);
		}

		int vals_ind = 0;
		while((read = getline(&line, &len, fp_in)) != -1){
			//read gives the length of the line
			//line is the data 
			char * tmp_char_ptr = strtok(line, " "); //split the line based on a space
			tmp_char_ptr = strtok(NULL, " ");	//for some reason need to do this to get second element in the string
			long double tmp_val = atof(tmp_char_ptr); //turn the rest of the line into a value
			
			//assign all the values
			if(vals_ind == 0)
				nCladding = tmp_val;
			else if(vals_ind == 1)
				nCore = tmp_val;
			else if(vals_ind == 2)
				lambd = tmp_val*1E-6;
			else if(vals_ind == 3)
				widthDomain = tmp_val*1E-6;
			else if(vals_ind == 4)
				lengthDomain = tmp_val*1E-6;
			else if(vals_ind == 5)
				w_B = tmp_val*1E-6;
			else if(vals_ind == 6)
				s = tmp_val*1E-6;
			else if(vals_ind == 7)
				w = tmp_val*1E-6;
			else if(vals_ind == 8)
				taper_start = tmp_val*1E-6;
			else if(vals_ind == 9)
				taper_end = tmp_val*1E-6;
			else if(vals_ind == 10)
				turn_start = tmp_val*1E-6;
			else if(vals_ind == 11)
				R = tmp_val*1E-6;
			else if(vals_ind == 12)
				angle = tmp_val * (long double)M_PI/((long double)180);
			else if(vals_ind == 13)
				tip_width = tmp_val*1E-6;
			else if(vals_ind == 14)
				sig = tmp_val*1E-6;
			else if(vals_ind == 15)
				del_x = tmp_val*1E-6;
			else if(vals_ind == 16)
				del_z = tmp_val*1E-6;
			else if(vals_ind == 17)
				alpha = tmp_val;
			else{
				printf("Extra line in source file\n Data = %s\n Aborting \n", line);
				exit(0);
			}
			vals_ind++;
		}
		fclose(fp_in);
		if(line)
			free(line);
	}


	//long double del_z = 0.25*1E-6;
	long double k0 = 2*M_PI/lambd;
	long double n_bar = (nCladding + nCore)/2;

	int Nzpts = round(lengthDomain/del_z);

	
	//long double del_x_2 = del_x*del_x;
	int N = round(widthDomain / del_x);	

	printf("\nN=%d\nNzpts=%d\nnCladding=%0.9Lf\nnCore=%0.9Lf\n", N, Nzpts, creall(nCladding), creall(nCore));
	printf("lambd(um)=%0.9Lf\nwidthDomain(um)=%0.9Lf\nlengthDomain(um)=%0.9Lf\n", lambd*1E6,widthDomain*1E6,lengthDomain*1E6);
	printf("w_B(um)=%0.9Lf\ns(um)=%0.9Lf\nw(um)=%0.9Lf\ntaper_start(um)=%0.9Lf\n",w_B*1E6,s*1E6,w*1E6,taper_start*1E6);
	printf("taper_end(um)=%0.9Lf\nturn_start(um)=%0.9Lf\nR(um)=%0.9Lf\nangle(rad)=%0.9Lf\n",taper_end*1E6,turn_start*1E6,R*1E6,angle);
	printf("tip_width(um)=%0.9Lf\nsig(um)=%0.9Lf\ndel_x(um)=%0.9Lf\ndel_z(um)=%0.9Lf\nalpha=%0.9Lf\n\n",tip_width*1E6,sig*1E6,del_x*1E6,del_z*1E6,alpha);


	//all the 1d arrays
	long double complex * u = (long double complex *) malloc (N * sizeof(long double complex));
	long double complex * b = (long double complex *) malloc (N * sizeof(long double complex));
	long double complex * gamma = (long double complex *) malloc (N * sizeof(long double complex));
	long double complex * r = (long double complex *) malloc (N * sizeof(long double complex));
	long double * x = (long double *) malloc (N * sizeof(long double));
	long double complex * kappa = (long double complex *) malloc (N * sizeof(long double complex));
	//linspace as well as a lot of allocation

	long double widthAbsEdge = 3*1E-6;
	long double kappa_max = -0.2;
		

	//all the 2d arrays
	long double complex ** allFeild = (long double complex **) malloc((Nzpts + 1)*sizeof(long double complex *));
	long double complex ** n = (long double complex **) malloc((Nzpts + 1)*sizeof(long double complex *));
	long double * z = (long double *) malloc ((Nzpts + 1) * sizeof(long double));
	long double z_step = del_z*((long double)Nzpts)/((long double)(Nzpts+1));

	for (int m = 0; m < Nzpts+1; m++){
		allFeild[m] = (long double complex *) malloc(N * sizeof(long double complex));
		n[m] = (long double complex *) malloc(N * sizeof(long double complex));
		z[m] = (long double) m*z_step;
		for (int j =0; j < N; j++){
			n[m][j] = nCladding;
		}
	}
	
	long double x_step = widthDomain/((long double)(N-1));
	for(int i = 0; i <N; i++){
		x[i] = -1*widthDomain/2 + ((long double)i*x_step);
		u[i] = expl(-1*(x[i]/sig)*(x[i]/sig));
		allFeild[0][i] = u[i];
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


	generate_n_grid(nCore, n, N, Nzpts, x, z, w_B, w, s, R, angle, taper_start, taper_end, turn_start, tip_width);

	if(save_n){
		//output the N data to check that it is right
		FILE * fp_n;
		fp_n = fopen("project_sim_n.csv","w+");

		
		for (int m = 0; m < Nzpts + 1; m ++){
			for (int j =0; j < N; j++){
				//fprintf(fp, "%Lf",cabsl(allFeild[m][j])*cabsl(allFeild[m][j]));
				fprintf(fp_n, "%0.7Lf", creall(n[m][j]));
				if(j != N-1){
					fprintf(fp_n,	",");
				}
			}
			if(m != Nzpts){
				fprintf(fp_n,	"\n");
			}
		}
		fclose(fp_n);
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
	
	save_data(argc, argv, allFeild, Nzpts, N);

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
	free(z);
	free(kappa);
	free(b);
	free(gamma);
	free(u);


	//printf("%lu\n", sizeof(long double complex));
	//printf("%d\n", N);
	//printf("%d\n", Nzpts+1);


}




