#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 2

double f1(double x[N])
{	
	return 4*pow(x[0],2) - 20*x[0] + pow(x[1],2)/4. + 8;
}

double f2(double x[N])
{	
	return x[0]*pow(x[1],2)*0.5 + 2*x[0] - 5*x[1] + 8;
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

void **pivoteamento(double **M, int NL, int NC)
{
        double l, pivot, maior, aux;
        int i, j, k, m, n;
        
        for(j=0; j<NL-1; j++)
        {	
			pivot = M[j][j];
			maior = pivot;

			for(k=j; k<NL; k++)
			{	
				if(fabs(maior) < fabs(M[k][j]))
				{	
					maior = M[k][j];
					m = k;
				}
			}

			if(maior != pivot)
			{	
				for(n=0; n<NC; n++)
				{
					aux = M[m][n];
					M[m][n] = M[j][n];
					M[j][n] = aux;
				}
			}

			for(i=j+1; i<NL; i++)
			{
				l = M[i][j]/M[j][j];       
				for(k=0; k<NC; k++)
				{
					M[i][k] -= l*M[j][k];
				} 
			}
		}       
}

void subsreversa(double **M, double *x, int dim)
{
	int i,j;
	double soma;
	
	for(i=dim-1; i>=0; i--)
	{
		soma = 0;
		
		for(j=i+1; j<dim; j++)
		{
			soma += M[i][j]*x[j];
		}
		
		x[i] = (M[i][dim] - soma)/M[i][i];
		//printf("a[%d]: %lf\n",i,x[i]);
	}
}

int main(int argc, char **argv)
{	
	double x0[N] = {0.1,0.1};
	double norm, norma;
	double eps=1e-5;
	double *y;
	double (*F[N])()={f1,f2};
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
	{	norma = norm = 0;
		
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
		
		pivoteamento(L, N, N+1);
		subsreversa(L, y, N);
		
		for(i=0; i<N; i++)
		{	
			norma += pow(x0[i],2);
			x0[i] = x0[i] + y[i];
			norm += pow(x0[i],2);
		}
		printf("%f \t%f \t%f\n",x0[0], x0[1], sqrt(fabs(norma-norm)));	
	}while(eps < sqrt(fabs(norma-norm)));
}
