#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 3

double f1(double x[N])
{	
	return 10*x[0] - 2*pow(x[1],2) + x[1] - 2*x[2] - 5;
}

double f2(double x[N])
{
	return 8*pow(x[1],2) + 4*pow(x[2],2) - 9;
}

double f3(double x[N])
{	
	return 8*x[1]*x[2] + 4;
}

double g1(double x[N]) //df1/dx0
{
	return 10;
}

double j2(double x[N]) //df1/dx1
{
	return -4*x[1] + 1;
}

double j3(double x[N]) //df1/dx2
{
	return -2;
}

double j4(double x[N]) //df2/dx0
{
	return 0;
}

double j5(double x[N]) //df2/dx1
{
	return -16*x[1];
}

double j6(double x[N]) //df2/dx2
{
	return 8*x[2];
}

double j7(double x[N]) //df3/dx0
{
	return 0;
}

double j8(double x[N]) //df3/dx1
{
	return 8*x[2];
}

double j9(double x[N]) //df3/dx2
{
	return 8*x[1];
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
	double x[N]={0}, x0[N]={0.1,0.1,-0.1};
	double eps = 1e-5;
	double norm, norma;
	double (*F[N])()={f1,f2,f3};
	double (*J[N][N])() = {{g1,j2,j3},{j4,j5,j6},{j7,j8,j9}};
	double **L, *y;
	int i, j;

	y = malloc(N*sizeof(double));
	
	L = malloc(N*sizeof(double *));
	for(i=0; i<N; i++)
	{
		L[i] = malloc((N+1)*sizeof(double));
	}
	
	do
	{	norma = norm = 0;
		for(i=0; i<N ; i++)
		{	
			for(j=0; j<N; j++)
			{
				L[i][j] = J[i][j](x0);
			}
			L[i][N] = -F[i](x0);
			norma += x0[i];
		}
		
		pivoteamento(L, N, N+1);
		
		
		subsreversa(L, y, N);
		
		for(i=0; i<N; i++ )
		{	
			x0[i] = x0[i] + y[i];
			norm += x0[i];
		}	
	}while(eps < sqrt(fabs(norma-norm)));
	
	for(i=0; i<N; i++)
	{
		printf("x[%d]:%f\n",i, x0[i]);
	}
	
}


