#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 3

double f1(double x[N])
{	return sqrt(7*x[1]*x[2]);
}

double f2(double x[N])
{	return 2*x[2] - 2*pow(x[1],2) - pow(x[0],2);
}

double f3(double x[N])
{	return pow(x[0],2)*0.1 +pow(x[1],2)*0.8;
}

int main(int argc, char **argv)
{	
	double xa[N]={1,2,3};
	double eps=1e-5;
	double norm, norma;
	double (*equacao[N])()={f1,f2,f3};
	int i;

	do
	{	norma = norm = 0;
		for(i=0; i<N; i++)
		{	
			norm += pow(xa[i], 2);
			xa[i] = equacao[i](xa);
			norma += pow(xa[i], 2);
		}
	}while(eps < sqrt(fabs(norma-norm)));
	
	for(i=0; i<N; i++)
	{
		printf("%lf\n", xa[i]);
	}
}
