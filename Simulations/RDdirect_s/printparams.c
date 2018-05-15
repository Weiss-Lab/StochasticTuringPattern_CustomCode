

extern int t;
extern int seed;
extern double Us;
extern double Vs;
extern double Cs;
extern double Ns;


void printparams()
{

	FILE *PARAMS;

	printf("\n--- PARAMETERS ---\n");
	printf("t=%ld\t, seed=%ld\n", t,seed);
	printf("totaltime=%g,  dt=%g, outputstep=%ld\n",totaltime, dt, outputstep);
	printf("Nx=%d, Ny=%d,  dx=%g,  dy=%g\n",Nx, Ny, dx, dy);
	printf("Us=%g, Vs=%g, Ius=%g, Ivs=%g, Cs=%g, Ns=%g\n", Us, Vs, Ius, Ivs, Cs, Ns);
	printf("GAMX=%g\n", GAMX);
	
	printf("gam_u=%g, gam_v=%g, gam_iu=%g, gam_iv=%g, gam_c=%g\n",gam_u, gam_v, gam_iu, gam_iv, gam_c);
	printf("afa_u=%g, afa_v=%g, afa_iu=%g, afa_iv=%g, afa_c=%g, afa_n=%g\n", afa_u, afa_v, afa_iu, afa_iv, afa_c, afa_n);
	printf("dif_u=%g, dif_v=%g, dif_iu=%g, dif_iv=%g, dif_c=%g, dif_n=%g\n", dif_u, dif_v, dif_iu, dif_iv, dif_c, dif_n);
	printf("the1=%d, the1=%d, the3=%d, the4=%d, the5=%d, the6=%d\n", the1, the2, the3, the4, the5, the6);
	printf("Kd_1=%g, Kd_2=%g, Kd_3=%g, Kd_4=%g, Kd_5=%g, Kd_6=%g\n", Kd_1, Kd_2, Kd_3, Kd_4, Kd_5, Kd_6);
	printf("Kc_1=%g, Kc_2=%g, Kc_3=%g\n", Kc_1, Kc_2, Kc_3);
	printf("fod_1=%g, fod_2=%g, fod_3=%g, fod_4=%g fod_5=%g, fod_6=%g\n", fod_1, fod_2, fod_3, fod_4, fod_5, fod_6);
	printf("lam_u=%g, lam_v=%g\n", lam_u, lam_v);
	printf("Ncri=%g, CHI=%g, Lac=%g, IPTG=%g\n", Ncri, CHI, Lac, IPTG);
	


	PARAMS = fopen("PARAMS_RUN","w");
	
	fprintf(PARAMS, "t=%ld\t, seed=%ld\n", t,seed);
	fprintf(PARAMS, "totaltime=%g, dt=%g, outputstep=%ld\n",totaltime, dt, outputstep);
	fprintf(PARAMS, "Nx=%d, Ny=%d, dx=%g, dy=%g\n",Nx, Ny, dx, dy);
	fprintf(PARAMS, "Us=%g, Vs=%g, Ius=%g, Ivs=%g, Cs=%g, Ns=%g\n", Us, Vs, Ius, Ivs, Cs, Ns);
	fprintf(PARAMS, "GAMX=%g\n", GAMX);
	
	fprintf(PARAMS, "gam_u=%g, gam_v=%g, gam_iu=%g, gam_iv=%g, gam_c=%g\n",gam_u, gam_v, gam_iu, gam_iv, gam_c);
	fprintf(PARAMS, "afa_u=%g, afa_v=%g, afa_iu=%g, afa_iv=%g, afa_c=%g, afa_n=%g\n", afa_u, afa_v, afa_iu, afa_iv, afa_c, afa_n);
	fprintf(PARAMS, "dif_u=%g, dif_v=%g, dif_iu=%g, dif_iv=%g, dif_c=%g, dif_n=%g\n", dif_u, dif_v, dif_iu, dif_iv, dif_c, dif_n);
	fprintf(PARAMS, "the1=%d, the1=%d, the3=%d, the4=%d, the5=%d, the6=%d\n", the1, the2, the3, the4, the5, the6);
	fprintf(PARAMS, "Kd_1=%g, Kd_2=%g, Kd_3=%g, Kd_4=%g, Kd_5=%g, Kd_6=%g\n", Kd_1, Kd_2, Kd_3, Kd_4, Kd_5, Kd_6);
	fprintf(PARAMS, "Kc_1=%g, Kc_2=%g, Kc_3=%g\n", Kc_1, Kc_2, Kc_3);
	fprintf(PARAMS, "fod_1=%g, fod_2=%g, fod_3=%g, fod_4=%g fod_5=%g, fod_6=%g\n", fod_1, fod_2, fod_3, fod_4, fod_5, fod_6);
	fprintf(PARAMS, "lam_u=%g, lam_v=%g\n", lam_u, lam_v);
	fprintf(PARAMS, "Ncri=%g, CHI=%g, Lac=%g, IPTG=%g\n", Ncri, CHI, Lac, IPTG);


	fclose(PARAMS);


}

