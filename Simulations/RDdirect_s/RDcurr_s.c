//Note:
/*
	1. 2008.07.22: bugs in "	dIuc[ab] = afa_iu*HFun_uv - gam_iu*Iv[ab] ;
								dIvc[ab] = afa_iv*HFun_uv - gam_iv*Iu[ab] ; "
			have been corrected: dIuc[ab] ... gam_iu*  Iu  [ab].

	2. 2008.07.24 Dif_u[ab]=dif_u, Dif_v[ab]=dif_v for simplicity with 1:100 fold 


  */


void RDcurr(double *U,  double *V,  double *Iu,  double *Iv,  double *C,  double *N, 
			double *dUc, double *dVc, double *dIuc, double *dIvc, double *dCc, double *dNc )
{

	int HIGG(int X, int Y);
	int LOWG(int X, int Y);

	int a, b, ag, as, bg, bs;
	int ab, agb, asb, abg, abs;
	double Ru, Rv, Leff, Ds_u, Ds_v, X1, X2, HFun_uv, HFun_c ;
	double Diff1, Diff2, Diff3;

	//allocation of difx[ab] and Reac[ab].
	double *Dif_u, *Dif_v;
	Dif_u = (double *) malloc( Nx * Ny * sizeof(double) );
	Dif_v = (double *) malloc( Nx * Ny * sizeof(double) );

	if ( (Dif_u==NULL) || (Dif_v==NULL) )  {	printf("RDcurr.c: \n");	exit(-2); }

/*
printf("rd-1\n") ;

				a=0, b=0;
			if (a<=Nx-2) {	ag = a+1; }
			else		 {	ag = 0;   }

			if (a>=1)	 {	as = a-1;  }
			else		 {	as = Nx-1; }


			if (b<=Ny-2) {	bg = b+1; }
			else		 {	bg = 0;   }

			if (b>=1)	 {	bs = b-1;  }
			else		 {	bs = Ny-1; }

			printf("a*Ny+Ny-1=%d\n", (Nx-1)*Ny+Ny-1);

			printf("a=%d,b=%d\n",a,b);
			ab=a*Ny+b;
			agb = ag*Ny+b; printf("agb=%d\n",agb);
			asb = as*Ny+b; printf("asb=%d\n",asb);
			abg = a*Ny+bg; printf("abg=%d\n",abg);
			abs = a*Ny+bs; printf("abs=%d\n",abs);

				printf("a=%d,b=%d\n",a,b);
				printf("dead\n");

				printf("%d\n",agb);
				printf("U[]: %f\n", U[agb]);
				
				printf("asb=%d\n",asb);
				printf("U[]=%f\n", U[asb]);
				
				printf("abg=%d\n",abg);
				printf("U[]=%f\n", U[abg]);

				printf("abs=%d\n",abs);
				printf("U[]=%f\n", U[abs]);

				printf("ab=%d\n",ab);
				printf("U[]=%f\n", U[ab]);
				printf("%lf %lf %lf %lf %lf\n", U[agb], U[asb], U[abg], U[abs], U[ab]);
*/				

	//calculate diffusion coefficients and the reaction part
	for(a=0;a<Nx;a++)
	{
		for(b=0;b<Ny;b++)
		{
			ab=a*Ny+b;

			Ru = lam_u*Iu[ab] ;
			Rv = lam_v*( 1+1/fod_5*pow(C[ab]/Kd_5,the5) )/( 1+pow(C[ab]/Kd_5,the5) ) ;

			Leff= Lac*( 1+1/fod_6*pow(IPTG/Kd_6, the6) )/( 1+pow(IPTG/Kd_6, the6) ); //effective Lac
							
			Ds_u = 1 + Ru/Kc_1 + Rv/Kc_3 ;
			Ds_v = 1 + Rv/Kc_2 ;

			//Assumption: molecular pools existed. molecules in cells are neglectiable. all cells have same density
//			X1 = pow(Ru*U[ab],2)/pow(Kc_1,3) /(1 + Ru/Kc_1 + Rv/Kc_3) /( 1+U[ab]/Kc_1 ) ;			
//			X2 = pow(Rv*V[ab],2)/pow(Kc_2,3) /(1 + Rv/Kc_2) /( 1+V[ab]/Kc_2+U[ab]/Kc_3 ) ;

			X1 = Ru*U[ab] ;			
			X2 = Rv*V[ab]/( 1+U[ab]/Kc_3 ) ;

			HFun_uv = ( 1+fod_1*pow(X1/Kd_1, the1) )/( 1+pow(X1/Kd_1, the1) )
					*( 1+1/fod_2*pow(C[ab]/Kd_2, the2) )/( 1+pow(C[ab]/Kd_2, the2) );

			HFun_c  = ( 1+fod_3*pow(X2/Kd_3, the3) )/( 1+pow(X2/Kd_3, the3) )
					*( 1+1/fod_4*pow( Leff/Kd_4, the4) )/( 1+pow( Leff/Kd_4, the4) ); 

		
			//effective diffusion coefficient
//			Dif_u[ab] = dif_u/Ds_u;
//			Dif_v[ab] = dif_v/Ds_v;

			Dif_u[ab] = dif_u ;
			Dif_v[ab] = dif_v ;

			//adding reaction part to dXc.
			dUc[ab]	 = afa_u*N[ab]*Iu[ab] - gam_u*U[ab] ;
			dVc[ab]	 = afa_v*N[ab]*Iv[ab] - gam_v*V[ab] ;

			dIuc[ab] = afa_iu*HFun_uv - gam_iu*Iu[ab] ;
			dIvc[ab] = afa_iv*HFun_uv - gam_iv*Iv[ab] ;
			dCc[ab]	 = afa_c *HFun_c  - gam_c * C[ab] ; 

			dNc[ab]	 = afa_n * N[ab] * ( 1-N[ab]/Ncri ) ;

			if( (Ru<0) || (Rv<0) || (X1<0) || (X2<0) )
			{
				printf("negative value: ab=%d, Ru=%lf, Rv=%lf, X1=%lf, X2=%lf\n", ab, Ru, Rv, X1, X2) ;
			}

		}//b
	}//a
	
//printf("rd-2\n") ;

	//adding diffusion part into dXc.
	for(a=0;a<Nx;a++)
	{
		for(b=0;b<Ny;b++)
		{

			ab= a*Ny+b;

//			ag = HIGG(a+1, Nx) ; 
//			as = LOWG(a-1, Nx) ; 
//			bg = HIGG(b+1, Ny) ;
//			bs = LOWG(b-1, Ny) ;

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

//				printf("U: %lf %lf %lf %lf %lf", 
//					Dif_u[agb], Dif_u[asb], Dif_u[abg], Dif_u[abs], Dif_u[ab]);
//				printf("U: %lf %lf %lf %lf %lf", 
//					U[agb], U[asb], U[abg], U[abs], U[ab]);

			Diff1 = Dif_u[agb]*U[agb] + Dif_u[asb]*U[asb] + Dif_u[abg]*U[abg] + Dif_u[abs]*U[abs] - 4*Dif_u[ab]*U[ab];
			Diff2 = Dif_v[agb]*V[agb] + Dif_v[asb]*V[asb] + Dif_v[abg]*V[abg] + Dif_v[abs]*V[abs] - 4*Dif_v[ab]*V[ab];

			Diff3 = dif_n*( N[agb] + N[asb] + N[abg] + N[abs] - 4*N[ab]) ;

			dUc[ab]	+= Diff1 ;
			dVc[ab]	+= Diff2 ;

			dNc[ab]	+= Diff3 ; //ATT: chemotatic term is not inlcuded.

/*			if ( (Diff1<0) || (Diff2<0) )
			{
				printf("a=%d,b=%d,Diff1=%lf,Diff2=%lf\n",a,b,Diff1,Diff2);
				printf("Dif_u[%d]=%f, Dif_v[%d]=%f\n",ab, Dif_u[ab], ab, Dif_v[ab]);
				printf("Dif_x: %lf %lf %lf %lf %lf\n", 
					Dif_u[agb], Dif_u[asb], Dif_u[abg], Dif_u[abs], Dif_u[ab]);
				printf("X: %lf %lf %lf %lf %lf\n", 
					U[agb], U[asb], U[abg], U[abs], U[ab]);
				
				exit(-2);
			}

			if ( ( (a==Nx-1) && (b==Ny-1) ) || ( (a==0) && (b==0) ) )
			{
				printf("Nx=%d, Ny=%d\n",Nx, Ny);
				printf("a=%d, b=%d, ag=%d, as=%d, bg=%d, bs=%d\n", a, b, ag, as, bg, bs);
				printf("a=%d,b=%d,Diff1=%lf,Diff2=%lf\n",a,b,Diff1,Diff2);
			}
*/

		}//a
	}//a

//printf("rd-3\n") ;

	//free
	free(Dif_u);
	free(Dif_v);

}//end

//(a+1,Nx),(b+1,Ny)
int HIGG(int X, int Y) //return a(b)+1
{
	if (X<=Y-1)
		return X;
	else
		return 0;

}

//(a-1,Nx), (b-1,Ny)
int LOWG(int X, int Y) //return a(b)-1
{
	if (X>=0)
		return X;
	else
		return Y;
}

