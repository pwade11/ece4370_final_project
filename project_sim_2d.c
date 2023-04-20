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
#define SWEEP_MODE 1 //toggle it into a mode where only prints out the in and out powers, and takes the s parameter as input


//gcc -O3 -o project_sim_2d_c project_sim_2d.c

void generate_n_grid(double complex nCore, double complex ** n, int N, int Nzpts, 
	double * x, double * z, double w_B, double w, double s, double R, double angle, 
	double taper_start, double taper_end, double turn_start, double tip_width){

	//assumes that n is aready created to be all nCladding before this function, since that can be done during initialization
	//this likely could be substancially optimized, becasue doing a lot more comparisons than required, but just starting this way
	//for now, since dont think will be a substancial time contributer

	int bus_start_ind = N;
	int bus_end_ind = 0;
	int coupler_start_ind = N; 
	int coupler_end_ind = 0;
	for(int j = 0; j<N; j++){
		if((x[j] >= -1*w_B/2) && j <bus_start_ind)
			bus_start_ind = j;
		if((x[j] <= w_B/2) && j > bus_end_ind)
			bus_end_ind = j;
		if((x[j] >= w_B/2 + s) && j < coupler_start_ind)
			coupler_start_ind = j;
		if((x[j] <= w_B/2 + s + w) && j > coupler_end_ind)
			coupler_end_ind = j;
	}

	double z0 = turn_start + (R+w)*sin(angle);
	double x_start_angle = w_B/2 + s + R + w - cos(angle)*(R+w);
	for (int m = 0; m < Nzpts+1; m++){

		for (int j =bus_start_ind; j <= bus_end_ind; j++){
			n[m][j] = nCore;
			//the bus waveguide
		}

		//know that the other waveguide stuff will start after the bus, and nothing above taper start
		
		if((z[m] >= taper_start) && (z[m] <= turn_start)){
			for(int j = coupler_start_ind; j <= coupler_end_ind; j++){
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
			}
		}
		else if(z[m] >= turn_start){
			for(int j = coupler_start_ind; j < N; j++){	
				if(z[m] <= turn_start + (R + w)*sin(angle)){
					//turn part
					if((x[j] >= w_B/2 + s + R + w - sqrt((R + w)*(R + w) - (z[m] - turn_start)*(z[m] - turn_start))) && (x[j] <= w_B/2 + s + R + w - sqrt(R*R - (z[m] - turn_start)*(z[m] - turn_start))))
						n[m][j] = nCore;
					if(x[j] > w_B/2 + s + R + w - sqrt(R*R - (z[m] - turn_start)*(z[m] - turn_start))){
						j = N; //stop it early if can
					}
				}
				else{
					//diagonal part
					if((x[j] >= x_start_angle + tan(angle)*(z[m] - z0)) && (x[j] <= x_start_angle + w + tan(angle)*(z[m] - z0)))
						n[m][j] = nCore;
					if(x[j] > x_start_angle + w + tan(angle)*(z[m] - z0)){
						j = N; //stop it early if can
					}
				}
			}
		}
				
	}
	return;
}

void save_data(int argc, char **argv, double complex ** allFeild, int Nzpts, int N){
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
			//buffer_count += sprintf( &file_buffer[buffer_count], "%.*Le", DECIMAL_DIG, (double) (allFeild[m][j]*conj(allFeild[m][j])));
			buffer_count += sprintf( &file_buffer[buffer_count], "%0.14f", (double) (allFeild[m][j]*conj(allFeild[m][j])));
			//the commented out one above specifies printing at full precision, but it its not necessary and this is much faster
			
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

void solve_tridiag(int N, int m, double complex * u, double complex * gamma, double complex a, double complex * b, double complex * r){
	/*
	double complex beta = b[0];
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
	//try in form these papers give: 
	//https://upcommons.upc.edu/bitstream/handle/2117/360430/TFM_SHARDOOL_KULKARNI.pdf?sequence=1&isAllowed=y
	//https://www.sciencedirect.com/science/article/pii/S001046552030357X?via%3Dihub
	//this way is substancially faster
	u[0] = r[0] / b[0]; //equivalent to d1 <-- d1/b1 in literature I think
	gamma[0] = a/b[0]; //eqivalent to their c1 part
	double complex temp_val; //r in the paper 
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
}

void read_input_file(int argc, char **argv, double complex * nCladding, double complex *nCore, 
	double * lambd, double * widthDomain, double * lengthDomain, double * w_B, double * s, 
	double * w, double * taper_start, double * taper_end, double * turn_start, double * R, 
	double * angle, double * tip_width, double * sig, double * del_x, double * del_z, 
	double * alpha){
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
		double tmp_val = atof(tmp_char_ptr); //turn the rest of the line into a value
		
		//assign all the values
		if(vals_ind == 0)
			*nCladding = tmp_val;
		else if(vals_ind == 1)
			*nCore = tmp_val;
		else if(vals_ind == 2)
			*lambd = tmp_val*1E-6;
		else if(vals_ind == 3)
			*widthDomain = tmp_val*1E-6;
		else if(vals_ind == 4)
			*lengthDomain = tmp_val*1E-6;
		else if(vals_ind == 5)
			*w_B = tmp_val*1E-6;
		else if(vals_ind == 6)
			*s = tmp_val*1E-6;
		else if(vals_ind == 7)
			*w = tmp_val*1E-6;
		else if(vals_ind == 8)
			*taper_start = tmp_val*1E-6;
		else if(vals_ind == 9)
			*taper_end = tmp_val*1E-6;
		else if(vals_ind == 10)
			*turn_start = tmp_val*1E-6;
		else if(vals_ind == 11)
			*R = tmp_val*1E-6;
		else if(vals_ind == 12)
			*angle = tmp_val * (double)M_PI/((double)180);
		else if(vals_ind == 13)
			*tip_width = tmp_val*1E-6;
		else if(vals_ind == 14)
			*sig = tmp_val*1E-6;
		else if(vals_ind == 15)
			*del_x = tmp_val*1E-6;
		else if(vals_ind == 16)
			*del_z = tmp_val*1E-6;
		else if(vals_ind == 17)
			*alpha = tmp_val;
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

void save_n_data(int Nzpts, int N, double complex ** n){
	FILE * fp_n;
	fp_n = fopen("project_sim_n.csv","w+");

	
	for (int m = 0; m < Nzpts + 1; m ++){
		for (int j =0; j < N; j++){
			//fprintf(fp, "%Lf",cabsl(allFeild[m][j])*cabsl(allFeild[m][j]));
			fprintf(fp_n, "%0.7f", creal(n[m][j]));
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

void init_x_u_kappa(int N, double widthDomain, double sig, double * x, double complex * u, 
	double complex ** allFeild, double widthAbsEdge, double kappa_max, double complex * kappa){
	double x_step = widthDomain/((double)(N-1));
	for(int i = 0; i <N; i++){
		x[i] = -1*widthDomain/2 + ((double)i*x_step);
		u[i] = exp(-1*(x[i]/sig)*(x[i]/sig));
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
}

void perform_beam_prop(int N, int Nzpts, double alpha, double del_x, double k0, double complex ** n, 
	double complex * kappa, double n_bar, double del_z, double complex * u, double complex * b,
	double complex * r, double complex * gamma, double complex ** allFeild){

	//memoisation with the k0*k0*(n[m+1][0]*n[m+1][0] - kappa[0]*kappa[0] + 2*I*n[m+1][0]*kappa[0] - n_bar*n_bar) term
	double complex * next_nmj_sqr_term = (double complex *) malloc (N * sizeof(double complex));


	double complex a = -alpha/(del_x*del_x);


	next_nmj_sqr_term[0] = k0*k0*(n[0][0]*n[0][0] - kappa[0]*kappa[0] + 2*I*n[0][0]*kappa[0] - n_bar*n_bar);
	for (int j = 1; j < N-1; j++){
		next_nmj_sqr_term[j] = k0*k0*(n[0][j]*n[0][j] - kappa[j]*kappa[j] + 2*I*n[0][j]*kappa[j] - n_bar*n_bar);
	}
	next_nmj_sqr_term[N-1] = k0*k0*(n[0][N-1]*n[0][N-1] - kappa[N-1]*kappa[N-1] + 2*I*n[0][N-1]*kappa[N-1] - n_bar*n_bar);

	double complex shared_r_b_term = 2*I*k0*n_bar/del_z;
	double complex later_part_of_r = -1*(2*(1-alpha)/(del_x*del_x)) + shared_r_b_term;
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
		
		solve_tridiag(N, m, u, gamma, a, b, r);	

		for (int j =0; j < N; j++){
			allFeild[m+1][j] = u[j];
		}

	}

	free(next_nmj_sqr_term);
}

void integrate_around_in_and_outvoid(int start_ind, double complex ** allFeild, double * x, double * z, int N, int Nzpts,
	double w_B, double s, double R, double w, double angle, double turn_start, double lambd,
	double * out_db, double * bus_db){

	//calc dist between the waveguides at the end
	double x_start_angle = w_B/2 + s + R + w - cos(angle)*(R+w);
	double z0 = turn_start + (R+w)*sin(angle);
	double dist_between = (x_start_angle + tan(angle)*(z[Nzpts] - z0)) - w_B/2;
	double tot_around_bus_start = 0;
	double tot_around_bus_end = 0;
	double tot_around_out = 0;

	for(int j = 0; j<N; j++){
		if((x[j] >= -1*w_B/2 - dist_between/2) && (x[j] <= w_B/2 + dist_between/2)){
			tot_around_bus_start += (double) (allFeild[start_ind][j]*conj(allFeild[start_ind][j]));
			tot_around_bus_end += (double) (allFeild[Nzpts][j]*conj(allFeild[Nzpts][j]));
		}
		if((x[j] >= x_start_angle + tan(angle)*(z[Nzpts] - z0) - dist_between/2) && (x[j] <= x_start_angle + w + tan(angle)*(z[Nzpts] - z0) + dist_between/2)){
			tot_around_out += (double) (allFeild[Nzpts][j]*conj(allFeild[Nzpts][j]));	
		}
	}

	*out_db = 10*log10l(tot_around_out/tot_around_bus_start);
	*bus_db = 10*log10l(tot_around_bus_end/tot_around_bus_start);
	//printf("%0.14f,%0.14f\n", 10*log10l(tot_around_out/tot_around_bus_start),10*log10l(tot_around_bus_end/tot_around_bus_start));

}

void find_max_in_and_out(int start_ind, double complex ** allFeild, double * x, double * z, int N, int Nzpts,
	double w_B, double s, double R, double w, double angle, double turn_start, double lambd){
	//this will be used for our evaluation of power transfer, and is ultimatly what we need
	//start_ind = where start looking (where assume bus reach steady state)

	//assumes the end happens on the z axis, which may not happen, need to make sure

	double bus_start_max = 0;
	double bus_end_max = 0;
	double output_end_max = 0;

	double x_start_angle = w_B/2 + s + R + w - cos(angle)*(R+w);
	double z0 = turn_start + (R+w)*sin(angle);
	for(int j = 0; j<N; j++){
		if((x[j] >= -1*w_B/2) && (x[j] <= w_B/2)){
			double calc_val_0 = (double) (allFeild[start_ind][j]*conj(allFeild[start_ind][j]));
			double calc_val_end = (double) (allFeild[Nzpts][j]*conj(allFeild[Nzpts][j]));
			//the bus waveguide
			if(calc_val_0 > bus_start_max){
				bus_start_max = calc_val_0;
			}
			if(calc_val_end > bus_end_max){
				bus_end_max = calc_val_end;
			}
		}
		if((x[j] >= x_start_angle + tan(angle)*(z[Nzpts] - z0)) && (x[j] <= x_start_angle + w + tan(angle)*(z[Nzpts] - z0))){
			double calc_val_out_end = (double) (allFeild[Nzpts][j]*conj(allFeild[Nzpts][j]));
			if(calc_val_out_end > output_end_max){
				output_end_max = calc_val_out_end;
			}
		}
	}

	double out_sum_db;
	double bus_sum_db;
	integrate_around_in_and_outvoid(start_ind, allFeild, x, z, N, Nzpts, w_B, s, R, w, angle, turn_start, lambd, &out_sum_db, &bus_sum_db);

	double transfer_power = 10*log10l(output_end_max/bus_start_max);
	double through_power = 10*log10l(bus_end_max/bus_start_max);
	if(SWEEP_MODE){
		printf("%0.14f,%0.14f,%0.14f,%0.14f,%0.14f,%0.14f,%0.14f,%0.14f,%0.14f\n",s,lambd,bus_start_max,bus_end_max,output_end_max,through_power,transfer_power, bus_sum_db, out_sum_db);
	}
	else{
		printf("bus_start_max=%0.14f\n", bus_start_max);
		printf("bus_end_max=%0.14f\n", bus_end_max);
		printf("output_end_max=%0.14f\n", output_end_max);
		printf("\n");
		printf("through_power=%0.14f\n",through_power);
		printf("transfer_power=%0.14f\n",transfer_power);
		printf("\n");
	}
	return;
}



int main(int argc, char **argv){
	

	if(SWEEP_MODE){
		if((argc != 3) && (argc != 4)){
			printf("Usage project_sim_2d_c <input file> s(um) lambda(um)(optional)\n Aborting...\n");
			exit(0);
		}
	}
	else{
		if((argc != 1) && (argc != 2) && (argc != 3)){
			printf("Usage project_sim_2d_c <input file (optional)> <output file(optional)>\n Aborting...\n");
			exit(0);
		}
	}

	int save_n = 0; //if you want to dump the n data, only use this to check the file, will not optimize that saving 

	struct timespec start_gen, start_processing, end_processing, end;
	
	clock_gettime(CLOCK_MONOTONIC, &start_gen);	
	
	//All of the parameters that you may want to adjust, should adjust them in the file

	double complex nCladding;	//cladding index
	double complex nCore;		//core (effective) index
	double lambd;				//wavelength (um)
	double widthDomain;		//width of total thing simulating (um)
	double lengthDomain; 		//length of total thing simulating (um)
	double w_B;				//width of the bus waveguide (um)
	double s;					//separation between the two waveguides (um)
	double w;					//width of the output waveguide (um)
	double taper_start;		//top of the taper for the output waveguide (um)
	double taper_end;			//bottom point of the taper for the output waveguide (um)
	double turn_start;			//startpoint for where the output waveguide turns
	double R;					//radius of the curve (um)
	double angle;				//angle the output waveguide curves to (degrees)
	double tip_width;			//the width of the tip of the taper of the output waveguide (um)
	double sig;				//width of the initial amplitude distribution (um)
	double del_x;				//discritization in x (um)
	double del_z;				//discritization in z (um)
	double alpha;				//mix between the two methods of simulation

	if(argc == 1 && !SWEEP_MODE){
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
		read_input_file(argc, argv, &nCladding, &nCore, &lambd, &widthDomain, &lengthDomain, &w_B, &s, &w, &taper_start, 
			&taper_end, &turn_start, &R, &angle, &tip_width, &sig, &del_x, &del_z, &alpha);
		if(SWEEP_MODE){
			//set s param based on arg
			s = (double) atof(argv[2])*1E-6;
			if(argc == 4)
				lambd = (double) atof(argv[3])*1E-6;
		}

	}

	//double del_z = 0.25*1E-6;
	double k0 = 2*M_PI/lambd;
	double n_bar = (nCladding + nCore)/2;

	int Nzpts = round(lengthDomain/del_z);

	//double del_x_2 = del_x*del_x;
	int N = round(widthDomain / del_x);	

	if(!SWEEP_MODE){
		printf("\nN=%d\nNzpts=%d\nnCladding=%0.9f\nnCore=%0.9f\n", N, Nzpts, creal(nCladding), creal(nCore));
		printf("lambd(um)=%0.9f\nwidthDomain(um)=%0.9f\nlengthDomain(um)=%0.9f\n", lambd*1E6,widthDomain*1E6,lengthDomain*1E6);
		printf("w_B(um)=%0.9f\ns(um)=%0.9f\nw(um)=%0.9f\ntaper_start(um)=%0.9f\n",w_B*1E6,s*1E6,w*1E6,taper_start*1E6);
		printf("taper_end(um)=%0.9f\nturn_start(um)=%0.9f\nR(um)=%0.9f\nangle(rad)=%0.9f\n",taper_end*1E6,turn_start*1E6,R*1E6,angle);
		printf("tip_width(um)=%0.9f\nsig(um)=%0.9f\ndel_x(um)=%0.9f\ndel_z(um)=%0.9f\nalpha=%0.9f\n\n",tip_width*1E6,sig*1E6,del_x*1E6,del_z*1E6,alpha);
	}

	//all the 1d arrays
	double complex * u = (double complex *) malloc (N * sizeof(double complex));
	double complex * b = (double complex *) malloc (N * sizeof(double complex));
	double complex * gamma = (double complex *) malloc (N * sizeof(double complex));
	double complex * r = (double complex *) malloc (N * sizeof(double complex));
	double * x = (double *) malloc (N * sizeof(double));
	double complex * kappa = (double complex *) malloc (N * sizeof(double complex));
	//linspace as well as a lot of allocation

	double widthAbsEdge = 3*1E-6;
	double kappa_max = -0.2;
		

	//all the 2d arrays
	double complex ** allFeild = (double complex **) malloc((Nzpts + 1)*sizeof(double complex *));
	double complex ** n = (double complex **) malloc((Nzpts + 1)*sizeof(double complex *));
	double * z = (double *) malloc ((Nzpts + 1) * sizeof(double));
	double z_step = del_z*((double)Nzpts)/((double)(Nzpts+1));

	for (int m = 0; m < Nzpts+1; m++){
		allFeild[m] = (double complex *) malloc(N * sizeof(double complex));
		n[m] = (double complex *) malloc(N * sizeof(double complex));
		z[m] = (double) m*z_step;
		for (int j =0; j < N; j++){
			n[m][j] = nCladding;
		}
	}
	
	init_x_u_kappa(N, widthDomain, sig, x, u, allFeild, widthAbsEdge, kappa_max, kappa);

	generate_n_grid(nCore, n, N, Nzpts, x, z, w_B, w, s, R, angle, taper_start, taper_end, turn_start, tip_width);

	if(save_n){
		//output the N data to check that it is right
		save_n_data(Nzpts, N, n);
	}

	clock_gettime(CLOCK_MONOTONIC, &start_processing);

	perform_beam_prop(N, Nzpts, alpha, del_x, k0, n, kappa, n_bar, del_z, u, b, r, gamma, allFeild);

	clock_gettime(CLOCK_MONOTONIC, &end_processing);
	
	if(!SWEEP_MODE){
		save_data(argc, argv, allFeild, Nzpts, N);
	}

	//find_max_in_and_out(50, allFeild, x, z, N, Nzpts, w_B, s, R, w, angle, turn_start, lambd);
	find_max_in_and_out(500, allFeild, x, z, N, Nzpts, w_B, s, R, w, angle, turn_start, lambd);
	clock_gettime(CLOCK_MONOTONIC, &end);
	
	if(!SWEEP_MODE){
		double time_from_init = BILLION *(end.tv_sec - start_gen.tv_sec) +(end.tv_nsec - start_gen.tv_nsec);
		time_from_init = time_from_init / BILLION;
		printf("Elapsed time from init: %lf seconds\n", time_from_init);

		double time_from_proc = BILLION *(end.tv_sec - start_processing.tv_sec) +(end.tv_nsec - start_processing.tv_nsec);
		time_from_proc = time_from_proc / BILLION;
		printf("Elapsed time from processing: %lf seconds\n", time_from_proc);	

		double time_proc = BILLION *(end_processing.tv_sec - start_processing.tv_sec) +(end_processing.tv_nsec - start_processing.tv_nsec);
		time_proc = time_proc / BILLION;
		printf("Elapsed time processing: %lf seconds\n", time_proc);

		double time_saving = BILLION *(end.tv_sec - end_processing.tv_sec) +(end.tv_nsec - end_processing.tv_nsec);
		time_saving = time_saving / BILLION;
		printf("Elapsed time saving: %lf seconds\n", time_saving);
	}

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


}




