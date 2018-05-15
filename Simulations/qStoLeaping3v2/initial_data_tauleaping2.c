//#include"convert_i2a.c"

extern double Us;
extern double Vs;
extern double Cs;
extern double Ns;

//parameters used:
// Us, Vs, Cs, Ns, Ius, Ivs

//ATT: For cell-growth and non-cell growth cases: growth: randomize initial density; no-growth: fix initial density be identical 
//ATT: iRndamp  (iRndamp = 0.01)


//void initial_data(double *U, double *V, double *Iu, double *Iv, double *C, double *N, int read)
//void initial_data(double *S, int read)
void initial_data(double *S, int read, int f_index)
{

//	void itoa_mine(int n, char s[]);
	char str[50]; // MUST be big enough to hold all the characters of your number!! 

	int a, b;
	long int ab;
	int unxx;
	double trash,temp,temp2,temp3,temp4,temp5,temp6;
	FILE *INPUT;

	//int read=0;
	
	//double iRndamp = 0.01; //random amplitude of initial values

	srand(time(0)*4321);

//	if (read==1) INPUT = fopen("FINAL", "r");

	itoa_mine(f_index, str) ; //convert integer to string
	if (read==1) 
	{
		INPUT = fopen(str, "r");
		printf("f_index=%d\n",f_index);
	}

	//INITIAL DATA
	for(a=0;a<Nx;a++)
	{
		for(b=0;b<Ny;b++)
		{
			ab=a*Ny+b;

			if (read==1)
			{
				fscanf(INPUT, "%lf %lf %lf %lf %lf %lf %lf %lf ",&trash,&trash,&temp,&temp2,&temp3,&temp4,&temp5,&temp6);
				S[ab*N_spg+0]  = floor(temp+0.49 ) ;
				S[ab*N_spg+1]  = floor(temp2+0.49 ) ;
				S[ab*N_spg+2]  = floor(temp3+0.49 ) ;
				S[ab*N_spg+3]  = floor(temp4+0.49 ) ;
				S[ab*N_spg+4]  = floor(temp5+0.49 ) ;
				S[ab*N_spg+5]  = floor(temp6) ;
				/*
				S[ab*N_spg+0]  = temp  ;
				S[ab*N_spg+1]  = temp2 ;
				S[ab*N_spg+2]  = temp3 ;
				S[ab*N_spg+3]  = temp4 ;
				S[ab*N_spg+4]  = temp5 ;
				S[ab*N_spg+5]  = temp6 ;
				*/
				printf("%lf %lf %lf %lf %lf %lf\n",temp,temp2,temp3,temp4,temp5,temp6);
			}
      
			if (read==0)
			{
				S[ab*N_spg+0]  = Us  ;
				S[ab*N_spg+1]  = Vs  ;
				S[ab*N_spg+2]  = Ius ;
				S[ab*N_spg+3]  = Ivs ;
				S[ab*N_spg+4]  = Cs  ;
				S[ab*N_spg+5]  = Ns  ;
				
				
				temp=(rand()*1.0/RAND_MAX-0.5);
				S[ab*N_spg+0] *= 1 + iRndamp*temp;
				
				temp=(rand()*1.0/RAND_MAX-0.5);
				S[ab*N_spg+1] *= 1 + iRndamp*temp;
				
				temp=(rand()*1.0/RAND_MAX-0.5);
				S[ab*N_spg+2] *= 1 + iRndamp*temp;
				
				temp=(rand()*1.0/RAND_MAX-0.5);
				S[ab*N_spg+3] *= 1 + iRndamp*temp;

				temp=(rand()*1.0/RAND_MAX-0.5);
				S[ab*N_spg+4] *= 1 + iRndamp*temp;

				//temp=(rand()*1.0/RAND_MAX-0.5);
				//S[ab*N_spg+5] *= 1 + iRndamp*temp;

//For no-cell growth case, cell density is set the same all around.

				//spherical initial condition
				//if ((a-Nx/2)*(a-Nx/2) + (b-Ny/2)*(b-Ny/2) < (Nx/4)*(Nx/4) )
				//	U[ab]=1.1103*(1.5);
			}

		}//b loop
	 }//a loop


	unxx = (Nx/2*Ny+Ny/2)*N_spg ;
//	S[unxx+0] = 100000;
//	S[unxx+1] = 100000;
//	S[unxx+5] = 1000;

	printf("\ncenter high! xx=%d\n",xx);


} //end program

