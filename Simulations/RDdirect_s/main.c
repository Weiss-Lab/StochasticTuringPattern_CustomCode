#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>


#include "allvaripat.c" // all variables and parameters
//#include "RDcurr.c"  //calculate the change of X:  dX for t--> t+dt. 
#include "RDcurr_s.c"  //calculate the change of X:  dX for t--> t+dt. simplified X1,X2 expression.

#include "printparams.c"	//
#include "initial_data.c"	//
#include "printframe.c"		//my code for making movie


int main()
{

	int a,b;
	long ab;
	FILE *DISini, *DISfin;

	//ARRAY DECLARATIONS
	double *U;
	double *V;
	double *Iu;
	double *Iv;
	double *C;
	double *N;

	double *U1;
	double *V1;
	double *Iu1;
	double *Iv1;
	double *C1;
	double *N1;

	double *U2;
	double *V2;
	double *Iu2;
	double *Iv2;
	double *C2;
	double *N2;

	double *dU1;
	double *dV1;
	double *dIu1;
	double *dIv1;
	double *dC1;
	double *dN1;

	double *dU2;
	double *dV2;
	double *dIu2;
	double *dIv2;
	double *dC2;
	double *dN2;

	//Dynamic allocation
	U  = (double *) malloc( Nx * Ny * sizeof(double) );
	V  = (double *) malloc( Nx * Ny * sizeof(double) );
	Iu = (double *) malloc( Nx * Ny * sizeof(double) );
	Iv = (double *) malloc( Nx * Ny * sizeof(double) );
	C  = (double *) malloc( Nx * Ny * sizeof(double) );
	N  = (double *) malloc( Nx * Ny * sizeof(double) );

	U1  = (double *) malloc( Nx * Ny * sizeof(double) );
	V1  = (double *) malloc( Nx * Ny * sizeof(double) );
	Iu1 = (double *) malloc( Nx * Ny * sizeof(double) );
	Iv1 = (double *) malloc( Nx * Ny * sizeof(double) );
	C1  = (double *) malloc( Nx * Ny * sizeof(double) );
	N1  = (double *) malloc( Nx * Ny * sizeof(double) );

	U2  = (double *) malloc( Nx * Ny * sizeof(double) );
	V2  = (double *) malloc( Nx * Ny * sizeof(double) );
	Iu2 = (double *) malloc( Nx * Ny * sizeof(double) );
	Iv2 = (double *) malloc( Nx * Ny * sizeof(double) );
	C2  = (double *) malloc( Nx * Ny * sizeof(double) );
	N2  = (double *) malloc( Nx * Ny * sizeof(double) );


	dU1  = (double *) malloc( Nx * Ny * sizeof(double) );
	dV1  = (double *) malloc( Nx * Ny * sizeof(double) );
	dIu1 = (double *) malloc( Nx * Ny * sizeof(double) );
	dIv1 = (double *) malloc( Nx * Ny * sizeof(double) );
	dC1  = (double *) malloc( Nx * Ny * sizeof(double) );
	dN1  = (double *) malloc( Nx * Ny * sizeof(double) );

	dU2  = (double *) malloc( Nx * Ny * sizeof(double) );
	dV2  = (double *) malloc( Nx * Ny * sizeof(double) );
	dIu2 = (double *) malloc( Nx * Ny * sizeof(double) );
	dIv2 = (double *) malloc( Nx * Ny * sizeof(double) );
	dC2  = (double *) malloc( Nx * Ny * sizeof(double) );
	dN2  = (double *) malloc( Nx * Ny * sizeof(double) );

	printf("all allocated\n");
	printparams();

	//RANDOM SEED
	seed=1163307032;//time(0);
	srand(seed);//time(0));
	//seed=time(0);
	//srand(time(0));

	initial_data(U, V, Iu, Iv, C, N, read);

	printframe(U, V, Iu, Iv, C, N, frame_ind, movie);

	printf("\nPrinting initial data to INITIAL file ... \n");
	DISini = fopen("INITIAL", "w");
	for(a=0;a<Nx;a++)	  
	{
		for(b=0;b<Ny;b++)
		{	
			ab=a*Ny+b;
			fprintf(DISini, "%lf %lf %lf %lf %lf %lf %lf %lf\n", a*dx, b*dy, U[ab],V[ab],Iu[ab],Iv[ab],C[ab],N[ab]);
		}
	}
	fclose(DISini);
	printf("\nInitial data Printed! \n");

	//TIME LOOP

	printf("Starting Time Loop..@time=%lf\n",t);


	//TIME EVOLUTION HERE

	for(t=1;t<=tsteps;t++)
	{
	    if (t%100==0)
		{
			printf("N_t=%d, t=%lf\n",t, t*dt);
		}
		//void printframe(fftw_real *U, fftw_real *V, fftw_real *C, fftw_real *N, int f_index, int mov_flag)
		if (t%300==0) 
		{	
			frame_ind +=1;
			printframe(U, V, Iu, Iv, C, N,frame_ind, movie);

			printf("\nt=%ld, time=%lf\n", t, t*dt);
			printf("U[%d]=%lf V[%d]=%lf \nIu[%d]=%lf Iv[%d]=%lf\n C[%d]=%lf N[%d]=%lf\n", 
				xx, U[xx],xx,V[xx],xx,Iu[xx],xx,Iv[xx],xx,C[xx],xx, N[xx]);


			double Ru = lam_u*Iu[xx] ;
			double Rv = lam_v*( 1+1/fod_5*pow(C[xx]/Kd_5,the5) )/( 1+pow(C[xx]/Kd_5,the5) ) ;							
			double Ds_u = 1 + Ru/Kc_1 + Rv/Kc_3 ;
			double Ds_v = 1 + Rv/Kc_2 ;


			//Assumption: molecular pools existed. molecules in cells are neglectiable. all cells have same density
			//double X1 = pow(Ru*U[xx],2)/pow(Kc_1,3) /Ds_u /( 1+U[xx]/Kc_1 ) ;			
			//double X2 = pow(Rv*V[xx],2)/pow(Kc_2,3) /Ds_v /( 1+V[xx]/Kc_2+U[xx]/Kc_3 ) ;
			double X1 = Ru*U[xx];			
			double X2 = Rv*V[xx]/( 1+U[xx]/Kc_3 ) ;

			printf("X1[%d]=%lf, X2[%d]=%lf \n",xx,X1,xx,X2);
			printf("Ru[%d]=%lf, Rv[%d]=%lf \n",xx, Ru,xx,Rv);
			//printf("Dif_u[%d]=%lf, Dif_v[%d]=%lf \n", xx, dif_u/Ds_u,xx, dif_v/Ds_v);

		}


		//ADVANCE TIME STEP, RETURNS U, U_NONLINEAR
		//void RDcurr(double *U,  double *V,  double *Iu,  double *Iv,  double *C,  double *N, 
		//			double *dUc, double *dVc, double *dIuc, double *dIvc, double *dCc, double *dNc )

//printf("1\n") ;

		RDcurr(U, V, Iu, Iv, C, N,  dU1, dV1, dIu1, dIv1, dC1, dN1); 

//printf("2\n") ;

		for (a=0; a<Nx; a++)
		{
			for (b=0; b<Ny; b++)
			{
				ab=a*Ny+b;

				 U1[ab] =  U[ab] + dU1 [ab] * dt ;
				 V1[ab] =  V[ab] + dV1 [ab] * dt ;

				Iu1[ab] = Iu[ab] + dIu1[ab] * dt ;
				Iv1[ab] = Iv[ab] + dIv1[ab] * dt ;

				 C1[ab] =  C[ab] + dC1 [ab] * dt ;
				 N1[ab] =  N[ab] + dN1 [ab] * dt ;

			}//b
		}//a


//printf("3\n") ;

		RDcurr(U1, V1, Iu1, Iv1, C1, N1,  dU2, dV2, dIu2, dIv2, dC2, dN2); 

//printf("4\n") ;

		for (a=0; a<Nx; a++)
		{
			for (b=0; b<Ny; b++)
			{
				ab=a*Ny+b;

				 U2[ab] =  U[ab] + 0.5*( dU1[ab] + dU2[ab])* dt ;
				 V2[ab] =  V[ab] + 0.5*( dV1[ab] + dV2[ab])* dt ;

				Iu2[ab] = Iu[ab] + 0.5*( dIu1[ab] + dIu2[ab])* dt ;
				Iv2[ab] = Iv[ab] + 0.5*( dIv1[ab] + dIv2[ab])* dt ;

				 C2[ab] =  C[ab] + 0.5*( dC1[ab] + dC2[ab])* dt ;
				 N2[ab] =  N[ab] + 0.5*( dN1[ab] + dN2[ab])* dt ;

			}//b
		}//a

//printf("5\n") ;

		// UPDATE
		for (a=0; a<Nx; a++)
		{
			for (b=0; b<Ny; b++)
			{
				ab=a*Ny+b;

				 U[ab] =  U2[ab] ;
				 V[ab] =  V2[ab] ;
				Iu[ab] = Iu2[ab] ;
				Iv[ab] = Iv2[ab] ;
				 C[ab] =  C2[ab] ;
				 N[ab] =  N2[ab] ;
/*
				 U[ab] =  U1[ab] ;
				 V[ab] =  V1[ab] ;
				Iu[ab] = Iu1[ab] ;
				Iv[ab] = Iv1[ab] ;
				 C[ab] =  C1[ab] ;
				 N[ab] =  N1[ab] ;
*/
				 if( (U[ab]<0) || (V[ab]<0) || (Iu[ab]<0) || (Iv[ab]<0) || (C[ab]<0) || (N[ab]<0) )
				 {
					int xx2=ab;
					printf("U[%d]=%lf V[%d]=%lf \nIu[%d]=%lf Iv[%d]=%lf\n C[%d]=%lf N[%d]=%lf\n", 
						xx2, U[xx2],xx2,V[xx2],xx2,Iu[xx2],xx2,Iv[xx2],xx2,C[xx2],xx2, N[xx2]);
					//xx2=xx2-1;
					//printf("U[%d]=%lf V[%d]=%lf \nIu[%d]=%lf Iv[%d]=%lf\n C[%d]=%lf N[%d]=%lf\n", 
					//	xx2, U[xx2],xx2,V[xx2],xx2,Iu[xx2],xx2,Iv[xx2],xx2,C[xx2],xx2, N[xx2]);
								 
					DISfin = fopen("Negative", "w");
					for(a=0;a<Nx;a++)	  
					{
						for(b=0;b<Ny;b++)
						{	
							ab=a*Ny+b;
							fprintf(DISfin, "%lf %lf %lf %lf %lf %lf %lf %lf\n", a*dx, b*dy, U[ab],V[ab],Iu[ab],Iv[ab],C[ab],N[ab]);
						}//b
					}//a
					exit(-2);
				 }

			}//b
		}//a		

//printf("6\n") ;

	}//end of time loop


	//PARAMETER LISTING
	printparams();

	DISfin = fopen("FINAL", "w");
	for(a=0;a<Nx;a++)	  
	{
		for(b=0;b<Ny;b++)
		{	
			ab=a*Ny+b;
			fprintf(DISfin, "%lf %lf %lf %lf %lf %lf %lf %lf\n", a*dx, b*dy, U[ab],V[ab],Iu[ab],Iv[ab],C[ab],N[ab]);
		}//b
	}//a


	//close file pointers
	fclose(DISini);
	fclose(DISfin);

	//DEALLOCATE 
	free(U);
	free(V);
	free(Iu);
	free(Iv);
	free(C);
	free(N);

	free(U1);
	free(V1);
	free(Iu1);
	free(Iv1);
	free(C1);
	free(N1);

	free(U2);
	free(V2);
	free(Iu2);
	free(Iv2);
	free(C2);
	free(N2);

	free(dU1);
	free(dV1);
	free(dIu1);
	free(dIv1);
	free(dC1);
	free(dN1);

	free(dU2);
	free(dV2);
	free(dIu2);
	free(dIv2);
	free(dC2);
	free(dN2);

	return 0;

}

