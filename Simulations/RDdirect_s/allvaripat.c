



//PARAMETERS
int Nx=64;
int Ny=64;

double dt=2e-3 ;
double totaltime=1000000*dt; 
int tsteps= int(totaltime/dt) ;
int outputstep = (tsteps/5); 

double dx=0.05;//0.005//1.0/Nx;
double dy=0.05;//0.005//1.0/Ny;
 

int xx=Nx/2*Ny+Ny/2;

double Us=1.0;//100.0;
double Vs=1.0;
double Ns=1.0;   //tried: 1000-patt, 5e-4-nopatt, 1.0--patt.
double Ius=1.0; //100.0;
double Ivs=1.0; //100.0;
double Cs=1.0;


double iRndamp=0.1;
int read = 0; //0--from given state; 1- from FINAL file
int frame_ind=100 ;
int movie=1 ;
int t=0;
int seed;

double GAMX=1.0;

double gam_u =GAMX*1.0; //0.01; //
double gam_v =GAMX*1.0; //0.005; //
double gam_c =GAMX*1.0;
double gam_iu=GAMX*1.0;
double gam_iv=GAMX*1.0;

double dif_u = 5e-1;  //1e-4/(dx*dy);//0.1;
double dif_v = 100*dif_u; //1e-5/(dx*dy);//0.005;
double dif_iu=0.0;
double dif_iv=0.0;
double dif_c =0.0;
double dif_n =0.0;


double afa_u =GAMX*3e1; //0.1;//
double afa_v =GAMX*3e1; //0.2;//
double afa_iu=GAMX*1e0;
double afa_iv=GAMX*0.3e0;
double afa_c =GAMX*1.0e0; //GAMX*1e0;
//double afa_n =GAMX*1.0;
double afa_n =GAMX*0.0;
//double afa_n=1.0;

double lam_u=1.0;
double lam_v=1e3;

double fod_1 =1e3;
double fod_2 =1e5;
double fod_3 =1e3;
double fod_4 =1e5;

double fod_5 =1e5; //1e5;

double fod_6 =1e5;

double Kd_1 =1e1;
double Kd_2 =1e1;  //Kd_2*Kd_2 
double Kd_3 =1e4;

double Kd_5 =1e2;

double Kc_3 =1e0;



double Kd_4 =1e1;
double Kd_6 =1e1;


double Kc_1 =1e1;
double Kc_2 =1e3;

double Lac = 0.0;
double IPTG = 0.0;

double Ncri = Ns;
double CHI=0.0;


#define the1 1
#define the2 2
#define the3 1
#define the4 4
#define the5 2
#define the6 1






