#include <string.h>

//code to output a,b,U,V,C,N to a file indexed with f_index.
void printframe(double *U, double *V, double *Iu, double *Iv, double *C, double *N, int f_index, int mov_flag)
{
	void itoa(int n, char s[]);

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
	itoa(f_index, str) ; //convert integer to string

	fp=fopen(str,"w");
	//print U, V, C, N
	for(a=0;a<Nx;a++)
	{
		for(b=0;b<Ny;b++)  
		{
			ab=a*Ny + b;
			fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf\n", a*dx, b*dy, U[ab], V[ab], Iu[ab], Iv[ab], C[ab], N[ab]);
			fflush(fp);
		}
	}
	fclose(fp);
	

}//program end 



/* itoa:  convert n to characters in s */
void itoa(int n, char s[])
{
	void reverse(char s[]);
    int i, sign;

    if ((sign = n) < 0)  /* record sign */
        n = -n;          /* make n positive */
    i = 0;
    do {       /* generate digits in reverse order */
        s[i++] = n % 10 + '0';   /* get next digit */
    } while ((n /= 10) > 0);     /* delete it */
    if (sign < 0)
        s[i++] = '-';
    s[i] = '\0';
    reverse(s);
} 

/* reverse:  reverse string s in place */
void reverse(char s[])
{
    int c, i, j;

    for (i = 0, j = strlen(s)-1; i<j; i++, j--) {
        c = s[i];
        s[i] = s[j];
        s[j] = c;
    }
}


