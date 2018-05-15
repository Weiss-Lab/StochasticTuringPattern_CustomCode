extern int Nx;
extern int Ny;

extern double Us;
extern double Vs;
extern double Cs;
extern double Ns;

//parameters used:
// Us, Vs, Cs, Ns, Ius, Ivs

//ATT: For cell-growth and non-cell growth cases: growth: randomize initial density; no-growth: fix initial density be identical 
//ATT: iRndamp  (iRndamp = 0.01)


void initial_data(double *U, double *V, double *Iu, double *Iv, double *C, double *N, int read)
{
	double trash,temp,temp2,temp3,temp4,temp5,temp6;
	FILE *INPUT;
	//int read=0;
	
	//double iRndamp = 0.01; //random amplitude of initial values

	srand(time(0)*4321);

	if (read==1) INPUT = fopen("FINAL", "r");

	//INITIAL DATA
	for(int a=0;a<Nx;a++)
	{
		for(int b=0;b<Ny;b++)
		{
			long int ab=a*Ny+b;

			if (read==1)
			{
				fscanf(INPUT, "%lf %lf %lf %lf %lf %lf %lf %lf ",&trash,&trash,&temp,&temp2,&temp3,&temp4,&temp5,&temp6);
				U[ab]  = temp;
				V[ab]  = temp2;
				Iu[ab] = temp3;
				Iv[ab] = temp4;
				C[ab]  = temp5;
				N[ab]  = temp6;
				printf("%lf %lf %lf %lf %lf %lf\n",temp,temp2,temp3,temp4,temp5,temp6);
			}
      
			if (read==0)
			{
				U[ab] = Us;
				V[ab] = Vs;
				Iu[ab]= Ius;
				Iv[ab]= Ivs;
				C[ab] = Cs;
				N[ab] = Ns;
				
				temp=(rand()*1.0/RAND_MAX-0.5);
				U[ab]+=iRndamp*temp*U[ab];
				
				temp=(rand()*1.0/RAND_MAX-0.5);
				V[ab]+=iRndamp*temp*V[ab];
				
				temp=(rand()*1.0/RAND_MAX-0.5);
				Iu[ab]+=iRndamp*temp*Iu[ab];
				
				temp=(rand()*1.0/RAND_MAX-0.5);
				Iv[ab]+=iRndamp*temp*Iv[ab];

				temp=(rand()*1.0/RAND_MAX-0.5);
				C[ab]+=iRndamp*temp*C[ab];

//For no-cell growth case, cell density is set the same all around.
				//temp=(rand()*1.0/RAND_MAX-0.5);
				//N[ab]+=iRndamp*temp*N[ab];

				//spherical initial condition
				//if ((a-Nx/2)*(a-Nx/2) + (b-Ny/2)*(b-Ny/2) < (Nx/4)*(Nx/4) )
				//	U[ab]=1.1103*(1.5);
			}

		}//b loop
	 }//a loop


	int unxx = Nx/2*Ny+Ny/2;
//	U[unxx] = 100000;
//	V[unxx] = 100000;
//	Iu[unxx] = 1000;

//	C[unxx] = 10;

	printf("\ncenter high! xx=%d\n",xx);


} //end program
