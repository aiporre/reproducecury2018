#include "setulb.cpp"
#include <stdio.h>

static double evaluate(const uvector &x,
		       uvector &g,
		       int N)
{
    int i;
    double fx = 0.0;
    
    for (i = 0; i < N; i += 2)
	{
	double t1 = 1.0 - x[i];
	double t2 = 10.0 * (x[i+1] - x[i] * x[i]);
	g[i+1] = 20.0 * t2;
	g[i] = -2.0 * (x[i] * g[i+1] + t1);
	fx += t1 * t1 + t2 * t2;
      }
    return fx;
}

#define N   100

int main(int argc, char *argv[])
{
    int i, ret = 0;
    double fx;
    uvector x(N, 0.0), g(N, 0.0);

    for (i = 0;i < N;i += 2)
      {
        x[i] = -1.2;
        x[i+1] = 1.0;
      }

    Lbfgs lb (N);    
    int k = 0;
    while (1)
      {
	int r = lb.iterate(x, fx, g);
	if (r == LBFGS_FG)
	    fx = evaluate(x, g, N);
	else if (r == LBFGS_NEW_X)
	  {
		k++;
		printf("Iteration %d:\n", k);
	    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
	    printf("\n");
	  }
	else
	  break;

      };

    return 0;
}
