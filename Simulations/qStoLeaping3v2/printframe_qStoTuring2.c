
//#include"convert_i2a.c"

//code to output a,b,U,V,C,N to a file indexed with f_index.
//void printframe(double *U, double *V, double *Iu, double *Iv, double *C, double *N, int f_index, int mov_flag)
void printframe(double *S, int f_index, int mov_flag)
{
//	void itoa_mine(int n, char s[]);

	char str[50]; // MUST be big enough to hold all the characters of your number!! 
	long a, b, ab;
	FILE *fp;

	//flag: whether needs to output files for making a movie.
	if (mov_flag!=1)
	{
		//printf("Don't make movie, no print\n");
		return ; 
	}

//	itoa(f_index, str, 10) ; //convert integer to string
	itoa_mine(f_index, str) ; //convert integer to string

	fp=fopen(str,"w");
	//print U, V, C, N
	for(a=0;a<Nx;a++)
	{
		for(b=0;b<Ny;b++)  
		{
			ab=a*Ny + b;
			fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf\n", a*dx, b*dy, S[ab*N_spg+0], S[ab*N_spg+1], 
				S[ab*N_spg+2], S[ab*N_spg+3], S[ab*N_spg+4], S[ab*N_spg+5]);
			fflush(fp);
		}
	}
	fclose(fp);
	

}//program end 

