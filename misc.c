#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "allvars.h"

const gsl_rng_type *rgtype, *rutype;
gsl_rng *rg, *ru;

void init_random_number(void)
{
    gsl_rng_env_setup();
    gsl_rng_default_seed = Random_number_seed;

    rgtype = gsl_rng_default;
    rg = gsl_rng_alloc (rgtype);

    rutype = gsl_rng_default;
    ru = gsl_rng_alloc (rutype);
}

void free_random_number(void)
{
    gsl_rng_free(rg);
    gsl_rng_free(ru);
}

double random_number(int flag)
{
    double rn;

    if(flag == 0) rn = gsl_rng_uniform(ru);

    if(flag == 1) rn = gsl_ran_ugaussian(rg);
//printf("flag=%d rn=%g\n", flag, rn);
    return rn;
}


double interpolate_bipoint(double *x, double *y, int n, double x0, int flag)
{
	int i;
	double x1, x2, dx, y0;

	x1 = x[0]; x2 = x[n-1];
	dx = (x2 - x1)/(n-1);

	i = (int) ((x0 - x1)/dx);
	if(i<0 || i>n-1)
	{
printf("x1=%g x2=%g x0=%g\n",x1, x2, x0);
		fprintf(stderr,"the interpolation point is out of range! %d\n", flag);
	}
	y0 = (y[i+1] - y[i])/(x[i+1] - x[i]) * (x0 - x[i]) + y[i];
	return y0;
}
