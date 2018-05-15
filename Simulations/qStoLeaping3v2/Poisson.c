// This is an ANSI C code for Poissonian random number generator.
// It calls "UnirandMersenne.cpp" for providing uniformly distributed random numbers (0,1).
// This is modified from the original code developed by Agner Fog(2008) in C++ library.
// Developed by Ting Lu @ 2008.11.19


//#include "UnirandMersenne.cpp" 
//#include "math.h"


const double SHAT1 = 2.943035529371538573;    // 8/e
const double SHAT2 = 0.8989161620588987408;   // 3-sqrt(12/e)
#define FAK_LEN 1024       // length of factorial table

/*
   // Variables used by Poisson distribution
   double pois_L_last;                           // previous value of L
   double pois_f0;                               // value at x=0 or at mode
   double pois_a;                                // hat center
   double pois_h;                                // hat width
   double pois_g;                                // ln(L)
   int32  pois_bound;                            // upper bound
*/

/***********************************************************************
Poisson distribution
***********************************************************************/
long Poisson (double L) 
{

	long PoissonLow(double L) ;
	long PoissonInver(double L) ;
	long PoissonRatioUniforms(double L) ;

   /*
   This function generates a random variate with the poisson distribution.

   Uses inversion by chop-down method for L < 17, and ratio-of-uniforms
   method for L >= 17.

   For L < 1.E-6 numerical inaccuracy is avoided by direct calculation.
   */

   //------------------------------------------------------------------
   //                 choose method
   //------------------------------------------------------------------
   if (L < 17) {
      if (L < 1.E-6) {
         if (L == 0) return 0;
         if (L < 0)  
		 {
			 printf("Parameter negative in poisson function\n") ;
			 exit(-2);
		 }

         //--------------------------------------------------------------
         // calculate probabilities
         //--------------------------------------------------------------
         // For extremely small L we calculate the probabilities of x = 1
         // and x = 2 (ignoring higher x). The reason for using this 
         // method is to prevent numerical inaccuracies in other methods.
         //--------------------------------------------------------------
         return PoissonLow(L);
      }    
      else {
         //--------------------------------------------------------------
         // inversion method
         //--------------------------------------------------------------
         // The computation time for this method grows with L.
         // Gives overflow for L > 80
         //--------------------------------------------------------------
         return PoissonInver(L);
      }
   }      
   else {
      if (L > 2.E9) 
		 {
			 printf("Parameter too big in poisson function\n") ;
			 exit(-2);
		 }

      //----------------------------------------------------------------
      // ratio-of-uniforms method
      //----------------------------------------------------------------
      // The computation time for this method does not depend on L.
      // Use where other methods would be slower.
      //----------------------------------------------------------------
      return PoissonRatioUniforms(L);
   }
}


/***********************************************************************
Subfunctions used by poisson
***********************************************************************/
//int32 StochasticLib1::PoissonLow(double L) {
long PoissonLow(double L) 
{
   /*
   This subfunction generates a random variate with the poisson 
   distribution for extremely low values of L.

   The method is a simple calculation of the probabilities of x = 1
   and x = 2. Higher values are ignored.

   The reason for using this method is to avoid the numerical inaccuracies 
   in other methods.
   */   
   double d, r;
   d = sqrt(L);
   if (UNIRND() >= d) return 0;
   r = UNIRND() * d;
   if (r > L * (1.-L)) return 0;
   if (r > 0.5 * L*L * (1.-L)) return 1;
   return 2;
}


//int32 StochasticLib1::PoissonInver(double L) {
long PoissonInver(double L) 
{
   /*
   This subfunction generates a random variate with the poisson 
   distribution using inversion by the chop down method (PIN).

   Execution time grows with L. Gives overflow for L > 80.

   The value of bound must be adjusted to the maximal value of L.
   */   
   const int bound = 130;              // safety bound. Must be > L + 8*sqrt(L).
   double r;                           // uniform random number
   double f;                           // function value
   long x;                            // return value
   double pois_L_last ;
   double pois_f0 ;

   pois_L_last = -1.;                         // Last values of Poisson parameters

   if (L != pois_L_last) {             // set up
      pois_L_last = L;
      pois_f0 = exp(-L);               // f(0) = probability of x=0
   }
   while (1) {  
      r = UNIRND();  x = 0;  f = pois_f0;
      do {                             // recursive calculation: f(x) = f(x-1) * L / x
         r -= f;
         if (r <= 0) return x;
         x++;
         f *= L;
         r *= x;                       // instead of f /= x
      }
      while (x <= bound);
   }
}  


//int32 StochasticLib1::PoissonRatioUniforms(double L) {
long PoissonRatioUniforms(double L) 
{
   /*
   This subfunction generates a random variate with the poisson 
   distribution using the ratio-of-uniforms rejection method (PRUAt).

   Execution time does not depend on L, except that it matters whether L
   is within the range where ln(n!) is tabulated.

   Reference: E. Stadlober: "The ratio of uniforms approach for generating
   discrete random variates". Journal of Computational and Applied Mathematics,
   vol. 31, no. 1, 1990, pp. 181-189.
   */

	double LnFac(long n) ;

   double u;                                          // uniform random
   double lf;                                         // ln(f(x))
   double x;                                          // real sample
   long k;                                           // integer sample
   long pois_bound ;
   double pois_a, pois_g, pois_h ;
   long mode ;
   double pois_f0 ;
   double pois_L_last ;

   pois_L_last = -1.;                         // Last values of Poisson parameters

   if (pois_L_last != L) {
      pois_L_last = L;                                // Set-up
      pois_a = L + 0.5;                               // hat center
      mode = (long)L;                          // mode
      pois_g  = log(L);
      pois_f0 = mode * pois_g - LnFac(mode);          // value at mode
      pois_h = sqrt(SHAT1 * (L+0.5)) + SHAT2;         // hat width
      pois_bound = (long)(pois_a + 6.0 * pois_h);    // safety-bound
   }
   while(1) {
      u = UNIRND();
      if (u == 0) continue;                           // avoid division by 0
      x = pois_a + pois_h * (UNIRND() - 0.5) / u;
      if (x < 0 || x >= pois_bound) continue;         // reject if outside valid range
      k = (long)(x);
      lf = k * pois_g - LnFac(k) - pois_f0;
      if (lf >= u * (4.0 - u) - 3.0) break;           // quick acceptance
      if (u * (u - lf) > 1.0) continue;               // quick rejection
      if (2.0 * log(u) <= lf) break;                  // final acceptance
   }
   return(k);
}


/***********************************************************************
Log factorial function
***********************************************************************/
double LnFac(long n) 
{
   // log factorial function. gives natural logarithm of n!

   // define constants
   static const double                 // coefficients in Stirling approximation     
      C0 =  0.918938533204672722,      // ln(sqrt(2*pi))
      C1 =  1./12., 
      C3 = -1./360.;
   // C5 =  1./1260.,                  // use r^5 term if FAK_LEN < 50
   // C7 = -1./1680.;                  // use r^7 term if FAK_LEN < 20
   // static variables
   static double fac_table[FAK_LEN];   // table of ln(n!):
   static int initialized = 0;         // remember if fac_table has been initialized

   double sum;
   int i;
   double n1, r ;

   if (n < FAK_LEN) {
      if (n <= 1) {
		  if (n < 0) { printf("Parameter negative in LnFac function\n"); exit(-2); }  
         return 0;
      }
      if (!initialized) {              // first time. Must initialize table
         // make table of ln(n!)
         sum = fac_table[0] = 0.;
         for (i=1; i<FAK_LEN; i++) {
            sum += log((double)(i));
            fac_table[i] = sum;
         }
         initialized = 1;
      }
      return fac_table[n];
   }
   // not found in table. use Stirling approximation
   n1 = n;  r  = 1. / n1;
   return (n1 + 0.5)*log(n1) - n1 + C0 + r*(C1 + r*r*C3);
}

