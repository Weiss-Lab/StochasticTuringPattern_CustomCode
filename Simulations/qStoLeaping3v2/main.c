
//Attention:	1. Updating of Propensity at each necessary time
//				2. zero/small rate coefficients. popultion growth rate, population death.
//				3. Model for pFNK-805 is different for pFNK-804


#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "time.h"

#include "UnirandMersenne.cpp" 
#include "Poisson.c" //after "math.h" and "UnirandMersenne.cpp"

#include "varipara.c"		

#define t_print	0.5			//time interval for print

#define Nx 64				//simulation size: X axis
#define Ny 64				//simulation size: Y axis
#define N_rpg 12//24		//reactions per grid
#define N_spg 6				//species per grid


#define N_sp	N_spg*Nx*Ny	//total speices number
#define N_rc	N_rpg*Nx*Ny	//total reaction channel number

//#define Nlgs	  1.0		//critical cell density for logistic growth 
#define Nlgs	  Ncri		//critical cell density for logistic growth 
#define n_crireac  3.0		//1.0 //10.0
#define epsilon	  0.03
#define Nstep_ssa 100		//Step numbers for SSA.


double	*S ;				//species, molecules
double  *rS ;				//molecular number changes due to reaction
double  *dS ;				//molecular number changes due to diffusion
double	*h ;				//propensity functions
int		*flag_cr ;			//flag indicating critical reactions
long	*KJ_value ;			//Kj value for each reaction

int		Q_ncr ;
double	k[N_rpg] ;			//rate coefficients for reaction in each signle grid.
double	g_value[N_spg] ;	//g value for each species 
double	t_now ;				//current time
double  h_total ;			//overall propensity function value

int read;
int frame_ind;

#include "convert_i2a.c"
//#include "initial_data_tauleaping.c"
#include "initial_data_tauleaping2.c" //auto-input a frame as an initial state
//#include "printframe_qStoTuring.c" 
#include "printframe_qStoTuring2.c" // for-input a frame as an initial state 
#include "printparams.c"

//Home-made subroutines
void PROPENSITY(void) ;
void CriticalReactions(void) ;
void PropensityAll(void) ;
double PropensityCriRea(void) ;
void gvalue_fun(void) ;
double tauprime(void) ;
double maxer(double x1, double x2) ;
double miner(double x1, double x2) ;
void SSA(void) ;
double NextTime(void) ;
int ChoosingChannel(void) ;
void UpdatingNum(int s) ;
double KJ_assign(double tau1, double tau2) ;
long Poisson(double mean) ;
long int Choosing1CRchan(void) ;
double tauprime2(void) ;
int Check_n_Update(int flag_CorU, double tau) ;
int AddingDiff(double tau) ;


int main(void)
//int main(int read, int frame_ind)
{
	long i, a, b ;
	long sab ;
	long num_count, pri_count ;
	double tau ;

	extern double *S, *h, *rS, *dS ;
	extern int *flag_cr ;
	extern long *KJ_value ;
	extern int Q_ncr ;
	extern double k[N_rpg], g_value[N_spg] ;
	extern double t_now ;
	extern double h_total ;


	double taup1, taup2 ;
	double t_now_old ;
	int nega_up ;
	int adif ;

	FILE *DISini ;

//	printf("read=%d,frame_ind=%d\n",read, frame_ind);

	printf("please input read and frame_ind..\n");
	scanf("%d %d", &read, &frame_ind);
	printf("read=%d,frame_ind=%d\n",read, frame_ind);

	//dynamic allocation
	S  = (double *) malloc( N_sp * sizeof(double) ) ;	//all molecular species
	dS = (double *) calloc( N_sp , sizeof(double) ) ;	//all molecular species
	rS = (double *) calloc( N_sp , sizeof(double) ) ;	//all molecular species
	h  = (double *) calloc( N_rc , sizeof(double) ) ;	//reactions propensity

	flag_cr = (int *) malloc( N_rc * sizeof(int) ) ;	//critical reaction indicator
	KJ_value  = (long *) calloc( N_rc , sizeof(long) ); //Kj value for each reaction

	printf("all allocated\n");

	//overall system time
	t_now = 0 ;

	//Print parameters
	printparams() ; 

	//RANDOM SEED
	init_genrand( (unsigned long)(time( NULL )*1234567) );

	//initialization of k[N_rpg]
	for(i=0;i<N_rpg;i++)
	{	k[i] = 1.0 ; }
//	k[10] = 0 ;	//population growth
//	k[11] = 0 ;	//population death


	//Input initial data
//	initial_data(S, read);
	initial_data(S, read, frame_ind);

	//print initial frame
	printframe(S, frame_ind, movie);

	//save initial frame
	DISini = fopen("INITIAL", "w");
	for(a=0;a<Nx;a++)	  
	{
		for(b=0;b<Ny;b++)
		{	
			sab = (a*Ny+b) * N_spg ;
			fprintf(DISini, "%lf %lf %lf %lf %lf %lf %lf %lf\n", a*dx, b*dy, 
				S[sab+0], S[sab+1], S[sab+2], S[sab+3], S[sab+4], S[sab+5] );
		}
	}
	fclose(DISini);
	printf("\nInitial data Printed! \n");


	////////////////////////////////////////////////
	//TIME LOOP
	printf("Starting Time Loop..@time=%lf\n",t);

	num_count = 0 ;	//iteration count
	pri_count = 0 ;	//printing count


	//calculate g_value. //void gvalue_fun(void)
	gvalue_fun() ;

	do
	{

STEP0:	

		//count the number of iteration
		num_count ++ ;
		if ( num_count%100==0 )
		{	printf("num_count=%d,\t t_now=%.10f\n", num_count, t_now); }
		
		//print frame for movie.
		if ( pri_count == (int)floor(t_now/t_print) )
		{	
			pri_count ++ ;
			frame_ind +=1 ;
			printframe(S, frame_ind, movie) ;
			printf("num_count=%d,\t t_now=%.15f\n", num_count, t_now); 
		}

		//step 0: Calculate propensity and propensityall
		PROPENSITY() ; //update h[i]

		PropensityAll() ; //return h_total = sum(h[i])

		//step 1: Identify critical reactions
		CriticalReactions() ; //to get: flag_cr[]

		//step 2: Calculate tau'
		taup1 = tauprime() ; //need gvalue_fun() and CriticalReactions() before tauprime().


STEP3:
		//step 3: if tau' < 1/a0(x), goto 100 steps of SSA. else goto step 4.
		if ( taup1 < 10*1.0/h_total )
		{
			t_now_old = t_now ;
			
			//run SSA.
			SSA();
						
			//update diffusion for SSA.
			adif = AddingDiff( t_now-t_now_old );
			if(adif==0)
			{	printf("Diffusion negative in SSA\n");
				exit(-2);
			}
			printf("SSA: t_now=%.10f\n", t_now) ;

			goto STEP0 ;
			//goto STEP3 ; //for SSA loop all the time.
		}
		else //taup1 > 10.0/h_total
		{
			//step 4: Calculate tau''
			//CriticalReactions() ; //to get: flag_cr[], has been done in Step 1. no chance since then. 

			taup2 = tauprime2() ;

			//step 5: compare tau' vs. tau'', Compute Kj and return tau.
			tau = KJ_assign(taup1, taup2) ;	

			//CHECK Kj
			/*for(i=0;i<N_rc;i++)		
			{
				if(KJ_value[i]<0)
				{
					printf("taup1-taup2=%.10f-%.10f=%.10f,tau=%.10f\n",taup1, taup2, taup1-taup2,tau);
					printf("KJ_value[%d]=%d\n",i,KJ_value[i]);
					exit(-2);
				}
			}
			*/
			//printf("taup1=%.10f, taup2=%.10f\n", taup1, taup2);

			//step 6: check and Update:
			nega_up = Check_n_Update(0, tau) ; //0--check

			// if any S_(t+tau)[i] negative, tau = tau/2. goto Step 3.
			if ( nega_up == 0 )
			{
				taup1 = taup1 * 0.5 ;
				printf("negative Langevin, tau halved! \n", tau) ; //1--update

				goto STEP3 ;
			}
			// else: t = t + tau, S_(t+tau) = S_t + Kj*vj, goto Step 1.
			else
			{
				//molecular number update
				Check_n_Update(1, tau) ;

				//time update
				t_now += tau ;

				if ( num_count%10==0 )
				{	printf("Langevin: t_now=%.10f\n", t_now) ;	}

			}//if-else

		}//if-else

		/*
		//Positive check
		for(i=0;i<N_rc;i++)		
		{
			if(h[i]<0)
			{
				printf("do end: h[%d]=%.3f\n",i,h[i]);
				exit(-2);
			}
		}
		for(i=0;i<N_sp;i++)		
		{
			if(S[i]<0)
			{
				printf("do end: S[%d]=%.3f\n",i,S[i]);
				exit(-2);
			}
		}
		*/

	}//do iteration
	while(t_now < totaltime ) ;
		
	////////////////////////////////////////////////
	free( S ) ;
	free(dS ) ;
	free(rS ) ;
	free( h ) ;
	free( flag_cr ) ;
	free( KJ_value ) ;

	printf("%d %d\n",N_sp, N_rc);

	return 0;

}

/////////////////////////////////////////////////////////////////////////

//calculate propensity functions
void PROPENSITY(void)
{
	int a, b, i;
//	int ag, as, bg, bs ;
	long ab, sab, rab;
//	long agb, asb, abg, abs;
	double U, V, Iu, Iv, C, N;
	double Ru, Rv, Leff, Ds_u, Ds_v, X1, X2, HFun_uv, HFun_c ; 

	extern double *S, *h ;
	extern double k[N_rpg] ;

	//printf("propensity starts\n");

	for(a=0; a<Nx; a++)
	{
		for(b=0; b<Ny; b++)
		{
			//printf("%d,%d\n",a,b);

			ab=a*Ny+b;
			sab = ab * N_spg ;
			rab = ab * N_rpg ;


			U  = S[sab+0] ;
			V  = S[sab+1] ;
			Iu = S[sab+2] ;
			Iv = S[sab+3] ;
			C  = S[sab+4] ;
			N  = S[sab+5] ;


			//prefactors
			Ru = lam_u * Iu ;
			Rv = lam_v * ( 1+1/fod_5*pow(C/Kd_5,the5) )/( 1+pow(C/Kd_5,the5) ) ; //Model 804
//			Rv = lam_v ; //Model 805

			Leff= Lac*( 1+1/fod_6*pow(IPTG/Kd_6, the6) )/( 1+pow(IPTG/Kd_6, the6) ); //effective Lac
							
			Ds_u = 1 + Ru/Kc_1 + Rv/Kc_3 ;
			Ds_v = 1 + Rv/Kc_2 ;

			//Assumption: molecular pools existed. molecules in cells are neglectiable. all cells have same density
//			X1 = pow(Ru*U[ab],2)/pow(Kc_1,3) /(1 + Ru/Kc_1 + Rv/Kc_3) /( 1+U[ab]/Kc_1 ) ;			
//			X2 = pow(Rv*V[ab],2)/pow(Kc_2,3) /(1 + Rv/Kc_2) /( 1+V[ab]/Kc_2+U[ab]/Kc_3 ) ;

			X1 = Ru*U ;			
			X2 = Rv*V/( 1+U/Kc_3 ) ;

			HFun_uv = ( 1+fod_1*pow(X1/Kd_1, the1) )/( 1+pow(X1/Kd_1, the1) )
					*( 1+1/fod_2*pow(C/Kd_2, the2) )/( 1+pow(C/Kd_2, the2) );

			HFun_c  = ( 1+fod_3*pow(X2/Kd_3, the3) )/( 1+pow(X2/Kd_3, the3) )
					*( 1+1/fod_4*pow( Leff/Kd_4, the4) )/( 1+pow( Leff/Kd_4, the4) ); 

			
			//Propensity functions
			h[rab +  0] =  afa_u * N * Iu ;			// k[ 0] * N * Iu	;
			h[rab +  1] =  gam_u * U;				// k[ 1] * U		;
			
			h[rab +  2] =  afa_v * N * Iv ;			// k[ 2] * N * Iv	;
			h[rab +  3] =  gam_v * V ;				// k[ 3] * V		;
			
			h[rab +  4] =   afa_iu * HFun_uv ;		// k[ 4] * HFun_uv	;
			h[rab +  5] =   gam_iu * Iu ;			//k[ 5] * Iu		;
			
			h[rab +  6] =   afa_iv * HFun_uv ;		//k[ 6] * HFun_uv;
			h[rab +  7] =   gam_iv * Iv ;			//k[ 7] * Iv		;

			h[rab +  8] =   afa_c * HFun_c	;		//k[ 8] * HFun_c	;
			h[rab +  9] =   gam_c * C ;				//k[ 9] * C		;

			h[rab + 10] =   afa_n * N * fabs(1-N*1.0/Nlgs) ;	//k[10] * N * fabs(1-N*1.0/Nlgs)	;
			h[rab + 11] =   gam_n * N ;				//k[11] * N		;
/*
			h[rab + 12] =   dif_u * U ;				//k[12] * U		;
			h[rab + 13] =   dif_u * U ;				//k[13] * U		;
			h[rab + 14] =   dif_u * U ;				//k[14] * U		;
			h[rab + 15] =   dif_u * U ;				//k[15] * U		;

			h[rab + 16] =   dif_v * V ;				//k[16] * V		;
			h[rab + 17] =   dif_v * V ;				//k[17] * V		;
			h[rab + 18] =   dif_v * V ;				//k[18] * V		;
			h[rab + 19] =   dif_v * V ;				//k[19] * V		;

			h[rab + 20] =   dif_n * N ;				//k[20] * N		;
			h[rab + 21] =   dif_n * N ;				//k[21] * N		;
			h[rab + 22] =   dif_n * N ;				//k[22] * N		;
			h[rab + 23] =   dif_n * N ;				//k[23] * N		;
*/
			//CHECK Negativity
			for(i=0;i<N_rpg;i++)
			{
				if(h[rab+i]<0)
				{	
					printf("in Pro: h[%d+%d]=%.3f\n",rab,i, h[rab+i]); 
					exit(-2);
				}
			}

			/*
			if( (HFun_uv<0)||(HFun_c<0) )
			{

				printf("in Pro: HFun_uv=%.10f,HFun_c=%.10f\n", HFun_uv, HFun_c) ;
				for(i=0;i<N_spg;i++)
				{
					printf("S[%d]=%.3f\n",sab+i, S[sab+i]);
				}
				printf("Ru=%.10f, Rv=%.10f\n",Ru,Rv);
				printf("X1=%.10f, X2=%.10f\n",X1,X2);

			}
			*/
			

		}//b
	}//a

	//printf("end of propensity\n") ;
	return ;
}


////////////////////////////////////////////////////////////////////////
//total propensity, return: h_total.
void PropensityAll(void)
{
	extern double *S, *h ;
	long i ;
	extern double h_total ;

	h_total = 0 ;
	for(i=0; i<N_rc; i++)
	{
		h_total += h[i] ;
	}

	//printf("h_total=%f \n", h_total);
	return ;

}


/////////////////////////////////////////////////////////////////////////
//determine critical reactions %flag_cr
void CriticalReactions(void)
//need *h, *S, return: *flag_cr
//flag_cr[]: 0----if non-critical, 1--if critical.
 
{
	extern double *S, *h ;
	extern int *flag_cr ;
	extern int Q_ncr ;
	
	int a, b ;
	long ab, sab, rab;
	long i ;
	double Lj, U, V, Iu, Iv, C, N ;

	//Yes. or Not. of non-critical reaction 
	Q_ncr = 1 ;  

	//set all reactions as 'non-critical' (0)
	for(i=0; i<N_rc; i++)
	{
		flag_cr[i] = 0 ;
	}

	//find and reset 'critical' reactions (1)
	for(a=0; a<Nx; a++)
	{
		for(b=0; b<Ny; b++)
		{
			ab=a*Ny+b;
			
			sab = ab*N_spg ;
			rab = ab*N_rpg ;

			U  = S[sab+0] ;
			V  = S[sab+1] ;
			Iu = S[sab+2] ;
			Iv = S[sab+3] ;
			C  = S[sab+4] ;
			N  = S[sab+5] ;

			//three conditions: 1. vij is negative, 2. h>0, 3. Lj<n_crireac.
			//0
			//if ( h[rab+0]>0 )
			{	flag_cr[rab + 0] = 0 ; }
			//1
			Lj = U/1.0 ;
			if ( ( Lj < n_crireac )  && ( h[rab+1]>0 ) )
			{	flag_cr[rab + 1] = 1 ; }
			//2
			//if ( h[rab+2]>0 )
			{ flag_cr[rab + 2] = 0 ; }
			//3			
			Lj = V/1.0;
			if ( ( Lj < n_crireac ) && ( h[rab+3]>0 ) )
			{	flag_cr[rab + 3] = 1 ; }
			//4
			//if ( h[rab+4]>0 )
			{ flag_cr[rab + 4] = 0 ; }
			//5			
			Lj = Iu/1.0;
			if ( ( Lj < n_crireac ) && ( h[rab+5]>0 ) )
			{	flag_cr[rab + 5] = 1 ; }
			//6
			//if ( h[rab+6]>0 )
			{   flag_cr[rab + 6] = 0 ; }
			//7			
			Lj = Iv/1.0;
			if ( ( Lj < n_crireac ) && ( h[rab+7]>0 ) )  				
			{	flag_cr[rab + 7] = 1 ; }
			//8
			//if ( h[rab+8]>0 )
			{   flag_cr[rab + 8] = 0 ; }
			//9
			Lj = C/1.0;
			if ( ( Lj < n_crireac ) && ( h[rab+9]>0 ) )	
			{	flag_cr[rab + 9] = 1 ; }
			//10
			Lj = N/1.0;
			if ( ( Lj<n_crireac ) && ( N>Nlgs ) && ( h[rab+10]>0 ) ) //N>Nlgs, this term goes to negative, i.e., cell death
			{   flag_cr[rab + 10] = 1 ; }
			//11
			Lj = N/1.0;
			if ( ( Lj<n_crireac ) && ( h[rab+11]>0 ) )	
			{   flag_cr[rab + 11] = 1 ; }

/*
			//12-15
			Lj = U/1.0;
			if ( ( Lj < n_crireac ) && ( h[rab+12]>0 ) )
			{
				flag_cr[rab + 12] = 1 ;
				flag_cr[rab + 13] = 1 ;
				flag_cr[rab + 14] = 1 ;
				flag_cr[rab + 15] = 1 ;
			}
			//16-19
			Lj = V/1.0;
			if ( ( Lj < n_crireac ) && ( h[rab+16]>0 ) )
			{
				flag_cr[rab + 16] = 1 ;
				flag_cr[rab + 17] = 1 ;
				flag_cr[rab + 18] = 1 ;
				flag_cr[rab + 19] = 1 ;
			}
			//20-23
			Lj = N/1.0;
			if ( ( Lj < n_crireac ) && ( h[rab+20]>0 ) )
			{
				flag_cr[rab + 20] = 1 ;
				flag_cr[rab + 21] = 1 ;
				flag_cr[rab + 22] = 1 ;
				flag_cr[rab + 23] = 1 ;
			}
*/

		}//b
	}//a

}


/////////////////////////////////////////////////////////////////////
//total propensity of critical reactions
double PropensityCriRea(void)
//need: *h, *flag_cr. return: critical h_total
//need: CriticalReactions() 
{
	extern double *S, *h ;
	extern int *flag_cr ;

	long i ;
	double propen_cr ;
	
	long ind1;

	ind1=0;
	propen_cr = 0 ;
	for(i=0; i<N_rc; i++)
	{
		if (flag_cr[i]==1) //flag_cr[]=1, ==>critical reaction
		{	
			propen_cr += h[i] ;			
			ind1 ++ ;
		}
	}
//	printf("CriRea2=%.2f (1/100)\n\n", 100*ind1/(1.0*N_rc)) ;

	return propen_cr ;

}



//g values: highest order 
void gvalue_fun(void)
{
	extern double g_value[N_spg] ;

	g_value[0] = 1 ;	//U
	g_value[1] = 1 ;	//V
	g_value[2] = 2 ;	//Iu
	g_value[3] = 2 ;	//Iv
	g_value[4] = 1 ;	//C
	g_value[5] = 2 ;	//N
	
	return ;
}

//calculate tau' 
double tauprime(void)
{
	void gvalue_fun(void) ;
	void CriticalReactions(void) ;
	int sign(double x) ;

	double maxer(double x1, double x2) ;
	double miner(double x1, double x2) ;

	int a, b, c; 
//	int ag, as, bg, bs ;
	long ab, sabc, sab, rab;
//	long agb, asb, abg, abs, ragb, rasb, rabg, rabs ;
	double W1, W2, max1, max2, min_c, taup ;

	extern double *S, *h ;
	extern int *flag_cr ;
	extern double g_value[N_spg] ;

	extern int Q_ncr ; 

	double *mu ;
	double *sig2 ;

	//calculate g_value.
//	gvalue_fun() ;	
//	CriticalReactions() ;

	//if 'non-critical reactions' pool is empty, return tau' = infinity.
	if ( Q_ncr == 0) //In this system, there are some non-critical reactions all the time. e.g., R[ab*N_rpg+0]: -->U . 
	{	
		return 1.0e10 ; 
	}
	//else: do the following
	else
	{
		mu		= (double *) calloc( N_sp , sizeof(double) ) ;
		sig2	= (double *) calloc( N_sp , sizeof(double) ) ;
		
		//calculate mu and sig2
		for(a=0; a<Nx; a++)
		{
			for(b=0; b<Ny; b++)
			{
				ab=a*Ny+b;
				sab = ab*N_spg ;
				rab = ab*N_rpg ;
/*
				//a+1 =[0, Nx-1]
				if (a<=Nx-2) {	ag = a+1; }
				else		 {	ag = 0;   }
				//a-1 =[0, Nx-1]
				if (a>=1)	 {	as = a-1;  }
				else		 {	as = Nx-1; }
				//b+1 =[0, Ny-1]
				if (b<=Ny-2) {	bg = b+1; }
				else		 {	bg = 0;   }
				//b-1 =[0, Ny-1]
				if (b>=1)	 {	bs = b-1;  }
				else		 {	bs = Ny-1; }

				agb = ag*Ny+b; 
				asb = as*Ny+b; 
				abg = a*Ny+bg; 
				abs = a*Ny+bs;
				ragb= agb*N_rpg;
				rasb= asb*N_rpg;
				rabg= abg*N_rpg;
				rabs= abs*N_rpg;

*/

				//CALCULATING mu 
				//U
				mu[sab+0] = + (1-flag_cr[rab+0]) * h[rab+0] - (1-flag_cr[rab+1]) * h[rab+1] ;
//								-  4 * (1-flag_cr[rab+12]) * h[rab+12] 
//								+ (1-flag_cr[ragb+12]) * h[ragb+12] + (1-flag_cr[rasb+13]) * h[rasb+13]   
//								+ (1-flag_cr[rabg+14]) * h[rabg+14] + (1-flag_cr[rabs+15]) * h[rabs+15] ;

				//V
				mu[sab+1] = + (1-flag_cr[rab+2]) * h[rab+2] - (1-flag_cr[rab+3]) * h[rab+3] ;
//								-  4 * (1-flag_cr[rab+16]) * h[rab+16] 
//								+ (1-flag_cr[ragb+16]) * h[ragb+16] + (1-flag_cr[rasb+17]) * h[rasb+17]   
//								+ (1-flag_cr[rabg+18]) * h[rabg+18] + (1-flag_cr[rabs+19]) * h[rabs+19] ;

				//Iu
				mu[sab+2] = + (1-flag_cr[rab+4]) * h[rab+4] - (1-flag_cr[rab+5]) * h[rab+5] ;

				//Iv
				mu[sab+3] = + (1-flag_cr[rab+6]) * h[rab+6] - (1-flag_cr[rab+7]) * h[rab+7] ;

				//C
				mu[sab+4] = + (1-flag_cr[rab+8]) * h[rab+8] - (1-flag_cr[rab+9]) * h[rab+9] ;

				//N = S[sab+5]
				mu[sab+5] = sign(1-S[sab+5]/Nlgs)*(1-flag_cr[rab+10]) * h[rab+10] - (1-flag_cr[rab+11]) * h[rab+11] ; 
//								-  4 * (1-flag_cr[rab+20]) * h[rab+20] 
//								+ (1-flag_cr[ragb+20]) * h[ragb+20] + (1-flag_cr[rasb+21]) * h[rasb+21]   
//								+ (1-flag_cr[rabg+22]) * h[rabg+22] + (1-flag_cr[rabs+23]) * h[rabs+23] ;


				//CALCULATE sig2
				//U
				sig2[sab+0] = + (1-flag_cr[rab+0]) * h[rab+0] + (1-flag_cr[rab+1]) * h[rab+1] ;
//								+  4 * (1-flag_cr[rab+12]) * h[rab+12] 
//								+ (1-flag_cr[ragb+12]) * h[ragb+12] + (1-flag_cr[rasb+13]) * h[rasb+13]   
//								+ (1-flag_cr[rabg+14]) * h[rabg+14] + (1-flag_cr[rabs+15]) * h[rabs+15] ;

				//V
				sig2[sab+1] = + (1-flag_cr[rab+2]) * h[rab+2] + (1-flag_cr[rab+3]) * h[rab+3] ;
//								+  4 * (1-flag_cr[rab+16]) * h[rab+16] 
//								+ (1-flag_cr[ragb+16]) * h[ragb+16] + (1-flag_cr[rasb+17]) * h[rasb+17]   
//								+ (1-flag_cr[rabg+18]) * h[rabg+18] + (1-flag_cr[rabs+19]) * h[rabs+19] ;

				//Iu
				sig2[sab+2] = + (1-flag_cr[rab+4]) * h[rab+4] + (1-flag_cr[rab+5]) * h[rab+5] ;

				//Iv
				sig2[sab+3] = + (1-flag_cr[rab+6]) * h[rab+6] + (1-flag_cr[rab+7]) * h[rab+7] ;

				//C
				sig2[sab+4] = + (1-flag_cr[rab+8]) * h[rab+8] + (1-flag_cr[rab+9]) * h[rab+9] ;

				//N
				sig2[sab+5] = + (1-flag_cr[rab+10]) * h[rab+10] + (1-flag_cr[rab+11]) * h[rab+11] ;
//								+  4 * (1-flag_cr[rab+20]) * h[rab+20] 
//								+ (1-flag_cr[ragb+20]) * h[ragb+20] + (1-flag_cr[rasb+21]) * h[rasb+21]   
//								+ (1-flag_cr[rabg+22]) * h[rabg+22] + (1-flag_cr[rabs+23]) * h[rabs+23] ;

/*				
				for(i=0; i<N_spg; i++)
				{
					if ( (mu[ab*N_spg+i]>0) && (sig2[ab*N_spg+i]>0) )
						printf("%d, %.8f, %.8f\n",ab*N_spg+i, mu[ab*N_spg+i], sig2[ab*N_spg+i]) ;
				}
*/

			}//b
		}//a
 

		//Find minimum from maximum.	
		
		taup = 1.0e10 ;
		//printf("searching begining: taup=%.10f\n",taup);
		for(a=0; a<Nx; a++)
		{
			for(b=0; b<Ny; b++)
			{
				ab=a*Ny+b;
				for(c=0; c<N_spg; c++)
				{
					sabc = ab*N_spg + c;
					

					W1 = epsilon*S[sabc]/g_value[c] ;
					W2 = W1*W1 ;


					max1 = maxer(W1, 1)*1.0/fabs(mu[sabc]) ;
					//avoid infinitiy
					if( (fabs(mu[sabc]) <= 1e-10) || (max1>=1e10) )
					{	max1=1.0e10 ; }
					
					max2 = maxer(W2, 1)*1.0/sig2[sabc] ;					
					//avoid infinitiy
					if( (sig2[sabc] <= 1e-10) || (max2>=1e10) )
					{	max2=1.0e10 ; }

					//find minimum
					min_c = miner(max1, max2) ;
					
					/*
					if( (mu[sabc]==0)||(sig2[sabc]==0) )
					{
						printf("S[%d]=%.3f, g_value[c]=%.3f\n", sabc, S[sabc], g_value[c]);
						printf("max1=%.8f, max2=%.8f\n", max1, max2) ;
						printf("mu=%.10f, sig2=%.10f\n",mu[sabc], sig2[sabc]) ;
						printf("min_c=%.10f\n",min_c);
						exit(-2);
					}
					*/

					if ( min_c < taup )
					{
						taup = min_c ;
					}

				}//c
			}//b
		}//a

		//printf("taup=%.10f\n", taup);

		free(mu) ;
		free(sig2) ;

		//printf("returned taup=%.10f\n", taup) ;

		return taup ;

	}//if-else 

}


////////////////////////////////////////////////////////////////////
int sign(double x)
{
	if( x>0 )
	{	return 1 ; }
	else
	{	return -1 ; }

}


///////////////////////////////////////////////////////////////////
//maximum of two inputs
double maxer(double x1, double x2)
{
	if (x1>x2)
	{	return x1; }
	else
	{	return x2; }
}


//////////////////////////////////////////////////////////////////
//minimum of two inputs
double miner(double x1, double x2)
{
	if (x1>x2)
	{	return x2; }
	else
	{	return x1; }
}



///////////////////////////////////////////////////////////////
//Gillespie simulation code
void SSA(void)
{
	void PROPENSITY(void) ;
	void PropensityAll(void) ;

	double NextTime(void) ;
	int ChoosingChannel(void) ;
	void UpdatingNum(int s) ;

	extern double *S, *h ;
	extern double t_now ;

	long ind_x ;
	long r ;
	double tou ;

	for (ind_x=0 ; ind_x < Nstep_ssa; ind_x++)
	{

		//Calculate propensity.
		PROPENSITY() ;

		//Calculate toatl propensity.
		PropensityAll() ;

		//Calculate next reaction time.
		tou = NextTime();

		//update time.
		t_now += tou;			
				
		//Choosing reaction channel in SSA
		r = ChoosingChannel();
		
		//Update values of every species.
		UpdatingNum(r);

	}//ind_x loop

	//printf("SSA out!\n");
}



//////////////////////////////////////////////////////////////////
double NextTime(void)
//need: h_total <---PropensityAll() <-- h[] <-- PROPENSITY() : PropensityAll(), PROPENSITY().
/*To produce the time interval (tou). */
{
	//void PROPENSITY(void);	
	extern double *S, *h ;
	extern double h_total ;

	long i;
	double tou,rand01,rand_01;
	
	//get random number
	rand01=UNIRND();
	rand_01=log(1.0/rand01);

	//get total propentisy: h_total

	//Choosing a random time.
	tou=1.0/h_total*rand_01;	

	if (tou<0)
	{
		printf("\ntou=%f\n WRONG in NextTime subroutine\n",tou);
		for(i=0;i<N_rc;i++)
		{	
			printf("i=%d\t\th[%d]=%f\n",i, i, h[i]);
		}
		exit(-2);
	}

    return tou;

}

int ChoosingChannel(void)
//need: h_total <---PropensityAll() <-- h[] <-- PROPENSITY(): PROPENSITY(), PropensityAll().
/*To decide the chemical channel.*/
{
	//void PROPENSITY(void);
	extern double *S, *h ;
	extern double h_total ;

	long i,u;
	double rand02,S_chan,h_cumulant,h_cumulant_pre;
	
	//rand sample
	rand02=UNIRND();
	S_chan=h_total*rand02;
					
	h_cumulant=0;
	h_cumulant_pre=0;
	for(i=0;i<N_rc;i++)
	{
		h_cumulant=h_cumulant_pre+h[i];
		if(h_cumulant< S_chan)
		{ h_cumulant_pre=h_cumulant; }
		else 
		{ u=i; break; }		
	}

//	printf("SSA: u=%d\n",u);
	return u;
}

void UpdatingNum(int s)
//need: *S
/*To Update the numbers of molecules*/
{
	//rate constants of reactions. 
	extern double k[N_rpg] ;
	extern double *S, *h ;
	extern double t_now ;
	//extern double h[];

	int z0, z1, z2, z3 ;
//	int ag, as, bg, bs;
	long ab, sab;
//	long agb, asb, abg, abs ;
//	long sagb, sasb, sabg, sabs;
	//total = N_rpg*Nx*Ny ;

	//s=( z1*Ny+z2 )*N_rpg + z3
	z0 = (int)floor(s/N_rpg);
	z3 =  s - z0*N_rpg ;  // 
	z1 = (int)floor(z0/Ny ) ;
	z2 = z0 - z1*Ny ;

	ab	= z1*Ny+z2 ;
	sab = ab*N_spg ;

/*
	//z1+1 =[0, Nx-1]
	if (z1<=Nx-2) {	ag = z1+1; }
	else		 {	ag = 0;   }
	//z1-1 =[0, Nx-1]
	if (z1>=1)	 {	as = z1-1;  }
	else		 {	as = Nx-1; }
	
	//z2+1 =[0, Ny-1]
	if (z2<=Ny-2) {	bg = z2+1; }
	else		 {	bg = 0;   }
	//z2-1 =[0, Ny-1]
	if (z2>=1)	 {	bs = z2-1;  }
	else		 {	bs = Ny-1; }

	agb = ag*Ny+z2; 
	asb = as*Ny+z2; 
	abg = z1*Ny+bg; 
	abs = z1*Ny+bs;

	sasb = asb*N_spg ;
	sagb = agb*N_spg ;
	sabs = abs*N_spg ;
	sabg = abg*N_spg ;
*/
	//U  = S[sab+0] ;
	//V  = S[sab+1] ;
	//Iu = S[sab+2] ;
	//Iv = S[sab+3] ;
	//C  = S[sab+4] ;
	//N  = S[sab+5] ;


	//update the numbers of cells.
	switch (z3)
	{

		/*case :	{ 
			if( (>=1) )
			{	;	}
			break; 
		}
		case :	{ 
			if( (>=1)&&(>=) )
			{	;	}
			break; 
		}*/


		//propensity functions

		//h[ab*N_rpg +  0] =   k[ 0] * N * Iu	;
		case 0:	{ 
			{	S[sab+0] ++ ;	}
			break; 
		}
		//h[ab*N_rpg +  1] =   k[ 1] * U		;
		case 1:	{ 
			if( (S[sab+0]>=1) )
			{	S[sab+0]--;	}
			break; 
		}
			
		//h[ab*N_rpg +  2] =   k[ 2] * N * Iv	;
		case 2:	{ 
			{	S[sab+1] ++ ;	}
			break; 
		}
		//h[ab*N_rpg +  3] =   k[ 3] * V		;
		case 3:	{ 
			if( (S[sab+1]>=1) )
			{	S[sab+1]--;	}
			break; 
		}
			
		//h[ab*N_rpg +  4] =   k[ 4] * HFun_uv;
		case 4:	{ 
			{	S[sab+2]++;	}
			break; 
		}
		//h[ab*N_rpg +  5] =   k[ 5] * Iu		;
		case 5:	{ 
			if( (S[sab+2]>=1) )
			{	S[sab+2]--;	}
			break; 
		}
			
		//h[ab*N_rpg +  6] =   k[ 6] * HFun_uv;
		case 6:	{ 
			{	S[sab+3]++;	}
			break; 
		}
		//h[ab*N_rpg +  7] =   k[ 7] * Iv		;
		case 7:	{ 
			if( (S[sab+3]>=1) )
			{	S[sab+3]--;	}
			break; 
		}

		//h[ab*N_rpg +  8] =   k[ 8] * HFun_c	;
		case 8:	{ 
			{	S[sab+4]++;	}
			break; 
		}
		//h[ab*N_rpg +  9] =   k[ 9] * C		;
		case 9:	{ 
			if( (S[sab+4]>=1) )
			{	S[sab+4]--;	}
			break; 
		}

		//h[ab*N_rpg + 10] =   k[10] * N * (1-N*1.0/Nlgs)	;
		case 10:{ 
			if (S[sab+5]<Nlgs)
			{	S[sab+5]++;	}
			else
			{
				if( (S[sab+5]>=1) )
				{	S[sab+5]--;	}
			}
			break; 
		}
		//h[ab*N_rpg + 11] =   k[11] * N		;
		case 11:	{ 
			if( (S[sab+5]>=1) )
			{	S[sab+5]--;	}
			break; 
		}

/*
		//h[ab*N_rpg + 12] =   k[12] * U		;
		case 12:	{ 
			if( (S[sab+0]>=1) )
			{	S[sab+0]--;	S[sasb+0]++; }
			break; 
		}
		//h[ab*N_rpg + 13] =   k[13] * U		;
		case 13:	{ 
			if( (S[sab+0]>=1) )
			{	S[sab+0]--;	S[sagb+0]++; }
			break; 
		}
		//h[ab*N_rpg + 14] =   k[14] * U		;
		case 14:	{ 
			if( (S[sab+0]>=1) )
			{	S[sab+0]--;	S[sabs+0]++; }
			break; 
		}
		//h[ab*N_rpg + 15] =   k[15] * U		;
		case 15:	{ 
			if( (S[sab+0]>=1) )
			{	S[sab+0]--;	S[sabg+0]++; }
			break; 
		}

		//h[ab*N_rpg + 16] =   k[16] * V		;
		case 16:	{ 
			if( (S[sab+1]>=1) )
			{	S[sab+1]--;	S[sasb+1]++; }
			break; 
		}
		//h[ab*N_rpg + 17] =   k[17] * V		;
		case 17:	{ 
			if( (S[sab+1]>=1) )
			{	S[sab+1]--;	S[sagb+1]++; }
			break; 
		}
		//h[ab*N_rpg + 18] =   k[18] * V		;
		case 18:	{ 
			if( (S[sab+1]>=1) )
			{	S[sab+1]--;	S[sabs+1]++; }
			break; 
		}
		//h[ab*N_rpg + 19] =   k[19] * V		;
		case 19:	{ 
			if( (S[sab+1]>=1) )
			{	S[sab+1]--;	S[sabg+1]++; }
			break; 
		}

		//h[ab*N_rpg + 20] =   k[20] * N		;
		case 20:	{ 
			if( (S[sab+5]>=1) )
			{	S[sab+5]--;	S[sasb+5]++; }
			break; 
		}
		//h[ab*N_rpg + 21] =   k[21] * N		;
		case 21:	{ 
			if( (S[sab+5]>=1) )
			{	S[sab+5]--;	S[sagb+5]++; }
			break; 
		}
		//h[ab*N_rpg + 22] =   k[22] * N		;
		case 22:	{ 
			if( (S[sab+5]>=1) )
			{	S[sab+5]--;	S[sabs+5]++; }
			break; 
		}
		//h[ab*N_rpg + 23] =   k[23] * N		;
		case 23:	{ 
			if( (S[sab+5]>=1) )
			{	S[sab+5]--;	S[sabg+5]++; }
			break; 
		}
*/


	}//end of switch
	/////////////end of the updating trajectory//////////
	

}//end of function




//calculate tau'' 
double tauprime2(void)
//need: PropensityCriRea()
{
	double PropensityCriRea(void) ;
	extern double *S, *h ;
	extern int *flag_cr ;

	double tau2 ;
	double h_cr ;
	double rand01 ;

	//calculate overall propensity for critical reactions.
	h_cr = PropensityCriRea();
	
	rand01 = UNIRND();
	
	tau2 = 1.0/h_cr * log(1.0/rand01) ;

//	printf("h_cr=%.8f, tau2=%.8f\n", h_cr, tau2);

	return tau2 ;

}



///////////////////////////////////////////////////////////////////
double KJ_assign(double tau1, double tau2)
//need: *flag_cr, *R, Choosing1CRchan(), Poisson(), 
{
	long Poisson(double mean) ;
	long Choosing1CRchan(void) ;

	extern double *S, *h ;
	extern int *flag_cr ;
	extern long *KJ_value ;
	
	double tau_rtn ;

	int a, b, c;
	long ab, rabc, j_cri ;

	//pick the smaller time
	if (tau1<tau2)
	{
		tau_rtn = tau1 ; 
	}
	else
	{
		tau_rtn = tau2 ;
	}

	//assign the Kj values.
	for(a=0; a<Nx; a++)
	{
		for(b=0; b<Ny; b++)
		{
			ab=a*Ny+b;
			for(c=0; c<N_rpg; c++)
			{
				rabc = ab*N_rpg + c;

				if (flag_cr[rabc]==1) //if critical....goto zero.
				{
					KJ_value[rabc] = 0 ; 
				}
				else //non-critical: a poisson random number.
				{
					KJ_value[rabc] = Poisson(h[rabc]*tau_rtn) ;

//					printf("h[(%d*Ny+%d)*N_rpg+%d]*tau_rtn=%.10f\nKJ=%d\n",a,b,c,h[rabc]*tau_rtn,KJ_value[rabc]);
//					if(KJ_value[rabc]!=0)	{	exit(-2); }

				}//end-if

			}//c
		}//b
	}//a


	//find K_jc for the case: tau2=<tau1 .
	if (tau2<=tau1)
	{
		j_cri = Choosing1CRchan() ;
		KJ_value[j_cri] = 1 ;

//		printf("j_cri=%d,KJ=%d\n",j_cri, KJ_value[j_cri]);
	}

	return tau_rtn ;

}


/////////////////////////////////////////////////////////////////
//total propensity of critical reactions
long int Choosing1CRchan(void)
//need: *flag_cr, PropensityCriRea() 
{
	//void PROPENSITY(void) ;
	double PropensityCriRea(void) ;

	extern double *S, *h ;
	extern int *flag_cr ;
	long i, u;	
	double hcr_total, h_cumulant, h_cumulant_pre ;
	double rand02, S_chan;

	//calculate critical propensity
	hcr_total = PropensityCriRea() ;

	//sample random number
	rand02=UNIRND();
	S_chan=hcr_total*rand02;
					
	//pick reaction channel
	h_cumulant=0;
	h_cumulant_pre=0;
	for(i=0;i<N_rc;i++)
	{
		if (flag_cr[i]==1)
		{
			h_cumulant=h_cumulant_pre+h[i];
			if(h_cumulant< S_chan)
			{ h_cumulant_pre=h_cumulant; }
			else 
			{ u=i; break; }	

		}//if
	}//for

	return u ;

}



/////////////////////////////////////////////////////////////
int Check_n_Update(int flag_CorU, double tau)
// flag_CorU: 0---Check, 1---Update
//need: *S, *KJ_value
{
	extern double *S, *h, *rS, *dS ;
	extern long *KJ_value ;

	int a, b, ag, as, bg, bs;
	long i, ab, agb, asb, abg, abs;
	long sab, rab;
	long sasb, sagb, sabs, sabg;

	int val_rtn;

	val_rtn = 1 ; //returning value

	//if it is used for check
	if( flag_CorU == 0 )
	{
		for(a=0; a<Nx; a++)
		{
			for(b=0; b<Ny; b++)
			{
				ab=a*Ny+b;

				sab=ab*N_spg;
				rab = ab*N_rpg;

				//a+1 =[0, Nx-1]
				if (a<=Nx-2) {	ag = a+1; }
				else		 {	ag = 0;   }
				//a-1 =[0, Nx-1]
				if (a>=1)	 {	as = a-1;  }
				else		 {	as = Nx-1; }
				//b+1 =[0, Ny-1]
				if (b<=Ny-2) {	bg = b+1; }
				else		 {	bg = 0;   }
				//b-1 =[0, Ny-1]
				if (b>=1)	 {	bs = b-1;  }
				else		 {	bs = Ny-1; }

				agb = ag*Ny+b; 
				asb = as*Ny+b; 
				abg = a*Ny+bg; 
				abs = a*Ny+bs;

				sasb = asb*N_spg;
				sagb = agb*N_spg;
				sabs = abs*N_spg;
				sabg = abg*N_spg;

				//calculate reaction part[]
				rS[sab+0] = 1*KJ_value[rab+0] - 1*KJ_value[rab+1] ;

				rS[sab+1] = 1*KJ_value[rab+2] - 1*KJ_value[rab+3] ;

				rS[sab+2] = 1*KJ_value[rab+4] - 1*KJ_value[rab+5] ;

				rS[sab+3] = 1*KJ_value[rab+6] - 1*KJ_value[rab+7] ;
					
				rS[sab+4] = 1*KJ_value[rab+8] - 1*KJ_value[rab+9] ;

				rS[sab+5] = 1*KJ_value[rab+10] - 1*KJ_value[rab+11] ; 

				//calculate diffusion part90
				//U diffusion
				dS[sab+0] = dif_u*( -4*S[sab+0] + S[sagb+0] + S[sasb+0] + S[sabg+0] + S[sabs+0] ) * tau ;
				//V diffusion
				dS[sab+1] = dif_v*( -4*S[sab+1] + S[sagb+1] + S[sasb+1] + S[sabg+1] + S[sabs+1] ) * tau ;
				//Iu difussion
				//Iv difussion
				//C difussion
				//N diffusion
				dS[sab+5] = dif_n*( -4*S[sab+5] + S[sagb+5] + S[sasb+5] + S[sabg+5] + S[sabs+5] ) * tau ;

				//Negativity check.
				for(i=0; i<N_spg; i++)
				{	
					if( (S[sab+i]+rS[sab+i]+dS[sab+i])<0 )
					{
						val_rtn = 0 ;
						goto CorU_END ;
					}
				}//i loop

			}//b
		}//a

	}//end-if
	else //update
	{
		for(a=0;a<Nx;a++)
		{
			for(b=0;b<Ny;b++)
			{
				sab = (a*Ny+b) * N_spg ;
				
				S[sab+0] += rS[sab+0] + dS[sab+0] ;
				S[sab+1] += rS[sab+1] + dS[sab+1] ;
				S[sab+2] += rS[sab+2] ;
				S[sab+3] += rS[sab+3] ;
				S[sab+4] += rS[sab+4] ;
				S[sab+5] += rS[sab+5] + dS[sab+5] ;

			}//b
		}//a

	}

	//if loop all succeed, return all positive.
	
CorU_END:

	return val_rtn ; 
	
}




///////////////////////////////////////////////////////

int AddingDiff(double tau)
//return: 0 if negative appeared. 1 if all positive.
{
	int a, b, ag, as, bg, bs;
	long ab, agb, asb, abg, abs;
	long sab, sagb, sasb, sabg, sabs;

	extern double *S, *dS;

	int rtn=1;

	for(a=0;a<Nx;a++)
	{
		for(b=0;b<Ny;b++)
		{
			ab=a*Ny+b;
			sab=ab*N_spg;

			//a+1 =[0, Nx-1]
			if (a<=Nx-2) {	ag = a+1; }
			else		 {	ag = 0;   }
			//a-1 =[0, Nx-1]
			if (a>=1)	 {	as = a-1;  }
			else		 {	as = Nx-1; }
			//b+1 =[0, Ny-1]
			if (b<=Ny-2) {	bg = b+1; }
			else		 {	bg = 0;   }
			//b-1 =[0, Ny-1]
			if (b>=1)	 {	bs = b-1;  }
			else		 {	bs = Ny-1; }

			agb = ag*Ny+b; 
			asb = as*Ny+b; 
			abg = a*Ny+bg; 
			abs = a*Ny+bs;

			sasb = asb*N_spg;
			sagb = agb*N_spg;
			sabs = abs*N_spg;
			sabg = abg*N_spg;

			//U diffusion
			dS[sab+0] = dif_u*( -4*S[sab+0] + S[sagb+0] + S[sasb+0] + S[sabg+0] + S[sabs+0] ) * tau ;
			//V diffusion
			dS[sab+1] = dif_v*( -4*S[sab+1] + S[sagb+1] + S[sasb+1] + S[sabg+1] + S[sabs+1] ) * tau ;
			//N diffusion
			dS[sab+5] = dif_n*( -4*S[sab+5] + S[sagb+5] + S[sasb+5] + S[sabg+5] + S[sabs+5] ) * tau ;

			//negativity check 
			if( ( (S[sab+0]+dS[sab+0])<0 ) || ( (S[sab+1]+dS[sab+1])<0 ) || ( (S[sab+5]+dS[sab+5])<0 ) ) 
			{
				printf("Diffusion negatived!\n");
				rtn = 0 ;
				goto AddD_END ;
			}
			

		}//b
	}//a

AddD_END:

	return rtn;


}




/////////////////////////////////////////////////////
/*
int Check_n_Update(int flag_CorU, double tau)
// flag_CorU: 0---Check, 1---Update
//need: *S, *KJ_value
{
	int a, b, ag, as, bg, bs;
	long i, ab, agb, asb, abg, abs;
	long sab, rab, ragb, rasb, rabg, rabs;
	double temp_ab[N_spg] ;
	extern long *KJ_value ;
	extern double *S, *h, *dS ;

	int val_rtn;

	val_rtn = 1 ; //returning value

	for(a=0; a<Nx; a++)
	{
		for(b=0; b<Ny; b++)
		{
			ab=a*Ny+b;

			sab=ab*N_spg;
			rab = ab*N_rpg;


			//calculate temp_ab[]
			temp_ab[0] = 1*KJ_value[rab+0] - 1*KJ_value[rab+1] ;

			temp_ab[1] = 1*KJ_value[rab+2] - 1*KJ_value[rab+3] ;

			temp_ab[2] = 1*KJ_value[rab+4] - 1*KJ_value[rab+5] ;

			temp_ab[3] = 1*KJ_value[rab+6] - 1*KJ_value[rab+7] ;
				
			temp_ab[4] = 1*KJ_value[rab+8] - 1*KJ_value[rab+9] ;

			temp_ab[5] = 1*KJ_value[rab+10] - 1*KJ_value[rab+11] ; 


//			//calculate temp_ab[]
//			temp_ab[0] = 1*KJ_value[rab+0] - 1*KJ_value[rab+1] 
//				- KJ_value[rab+12]  - KJ_value[rab+13]  - KJ_value[rab+14]  - KJ_value[rab+15] 
//				+ KJ_value[ragb+12] + KJ_value[rasb+13] + KJ_value[rabg+14] + KJ_value[rabs+15]  ;

//			temp_ab[1] = 1*KJ_value[rab+2] - 1*KJ_value[rab+3] 
//				- KJ_value[rab+16]  - KJ_value[rab+17]  - KJ_value[rab+18]  - KJ_value[rab+19] 
//				+ KJ_value[ragb+16] + KJ_value[rasb+17] + KJ_value[rabg+18] + KJ_value[rabs+19]  ;


//			temp_ab[2] = 1*KJ_value[rab+4] - 1*KJ_value[rab+5] ;

//			temp_ab[3] = 1*KJ_value[rab+6] - 1*KJ_value[rab+7] ;
				
//			temp_ab[4] = 1*KJ_value[rab+8] - 1*KJ_value[rab+9] ;

//			temp_ab[5] = 1*KJ_value[rab+10] - 1*KJ_value[rab+11]  
//				- KJ_value[rab+20]  - KJ_value[rab+21]  - KJ_value[rab+22]  - KJ_value[rab+23] 
//				+ KJ_value[ragb+20] + KJ_value[rasb+21] + KJ_value[rabg+22] + KJ_value[rabs+23]  ;
 


			//if it is used for check
			if( flag_CorU == 0 )
			{

				for(i=0; i<N_spg; i++)
				{	
					if( (S[sab+i]+rS[sab+i]+dS[sab+i])<0 )
					//if( (S[sab+i]+rS[sab+i]+dS[sab+i])<0 )
					{
						val_rtn = 0 ;
						goto CorU_END ;
					}
				}//i loop

			}
			//if it is used for update
			else
			{

				S[sab+0] += temp_ab[0] ;
				S[sab+1] += temp_ab[1] ;
				S[sab+2] += temp_ab[2] ;
				S[sab+3] += temp_ab[3] ;
				S[sab+4] += temp_ab[4] ;
				S[sab+5] += temp_ab[5] ;

			}//if-else


		}//b
	}//a

	//if loop all succeed, return all positive.
	
CorU_END:

	return val_rtn ; 
	
}
*/

/////////////////////////////////////////////////////////////////

