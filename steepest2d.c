
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 3
#define eps 5e-2
#define dt 1e-6

double f1(double x[N])
{	
	return x[0] + cos(x[0]*x[1]*x[3]) - 1;
}

double f2(double x[N])
{	
	return pow(1-x[0],1/4) + x[1] + 0.05*pow(x[2],2) - 0.15*x[2];
}

double f3(double x[N])
{	
	return -pow(x[0],2) - 0.1*pow(x[1],2) + 0.01*x[1] + x[2] - 1;
}

double G(double (*f[N])(), double x[N])
{	
	int i;
	double soma=0;
	
	for(i=0; i<N; i++)
	{
		soma += pow(f[i](x), 2);
	}
	return soma;
}

double *gradiente(double **J, double (*F[N])(), double x[N], double *NormaGrad)
{	
	int i, j;
	double *grad;

	grad = malloc(N*sizeof(double));
	
	for(i=0 ; i<N; i++)
	{	
		for(j=0; j<N; j++)
		{
			grad[i] += 2*J[j][i]*F[j](x);
		}
		*NormaGrad += pow(grad[i],2);
	}
	
	*NormaGrad = sqrt(*NormaGrad);
	
	return grad;
}

double df(double f(), double x[N], int i)
{	
	double xant, dx, dy;

	xant = x[i];
	dx = 1e-6;
	
	x[i]+=dx;
	dy=f(x);
	
	x[i]-=2*dx;
	dy-=f(x);
	
	x[i] = xant;
	
	return dy/(2*dx);
}

double h(double x[N], double *grad, double a, double (*F[N])())
{	
	int i;
	double aux[N]={0};

	for(i=0; i<N; i++)
	{
		aux[i] = x[i] - a*grad[i];
	}
	
	return G(F, aux);
}

void imprime( double **M, int NL, int NC)
{
	int i, j;
	for(i=0; i<NL; i++)
	{
		for(j=0; j<NC; j++)
		{
			printf("%.2lf\t ", M[i][j]);
		}	  
		puts("");
	}	
}

int main(int argc, char **argv)
{	
	double x0[N]={0,0,0};
	double a[N], gn[N], k[N];
	double norm, norma, alpha, g, NormaGradg;
	double *y, *gradg;
	double (*F[N])()={f1,f2,f3};
	double **J, **L;
	int i, j;

	y = malloc(N*sizeof(double));
	L = malloc(N*sizeof(double));
	
	for(i=0; i<N; i++)
	{
		L[i] = malloc((N+1)*sizeof(double));
	}

	J = malloc(N*sizeof(double));
	
	for(i=0; i<N; i++)
	{
		J[i] = malloc(N*sizeof(double));
	}

	do
	{	
		norma = norm = 0;
		NormaGradg = 0;
		
		for(i=0; i<N; i++)
		{
			gn[i] = 0;
		}
		
		for(i=0; i<N; i++)
		{	
			for(j=0; j<N; j++)
			{
				J[i][j] = df(F[i], x0, j);
			}
		}
		
		for(i=0; i<N; i++)
		{	
			for(j=0; j<N; j++)
			{
				L[i][j] = J[i][j];
			}
			L[i][N] = -F[i](x0);
		}
		
		g = G(F, x0);
		gradg = gradiente(J, F, x0, &NormaGradg);
		
		for(i=0; i<N; i++)
		{
			gradg[i] /= NormaGradg;
		}
		
		a[0] = 0;
		a[2] = 1;
		
		while(h(x0, gradg, a[0], F) < h(x0, gradg, a[2], F))
		{	
			a[2] /= 2.;
			
			if( a[2] < 1e-6)
			{
				break;
			}
		}
		
		if(h(x0, gradg, a[0], F) < h(x0, gradg, a[2], F) )
		{	
			a[1] = a[2];
			a[2] = a[0];
			a[0] = a[1];
			a[1] = a[0]/2;
		}
		else
		{			
			a[1] = a[2]/ 2.;
		}
		
		for(i=0; i<N; i++)
		{
			gn[i] = h(x0, gradg, a[i], F);
		}
		
		for(i=0; i<N-1; i++)
		{  
			k[i] = (gn[i+1] - gn[i])/(a[i+1] - a[i]);
		}
		
		k[2] = (k[1] - k[0])/(a[1] - a[0]);

		alpha = (k[2]* a[1] - k[0])/(2* k[2]);

		for(i=0; i<N; i++)
		{	
			norma += pow(x0[i],2);
			x0[i] -= alpha*gradg[i];
			norm += pow(x0[i],2);
		}
	
		printf("%f \t%f \t%f \t%f\n",x0[0], x0[1], x0[2], sqrt(fabs(norma-norm)));	
	}while(eps < sqrt(fabs(norma-norm)));
}

