#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#include "cosmo.h"
#include "allvars.h"
#include "proto.h"

void readparameters()
{
	char base[200];

	Mass_min = 1.e7;
	Mass_max = 1.e15;
	Num_Bin_LnM = 300;

	Redshift_max = 14.5;
	Redshift_min = 0.0;
	Num_Bin_Z = 300;
	Num_Bin_Ju = 300;

	Redshift = 0.0;

	Mass_res = 2.726e9 ;//* 0.1;
	//Mass_res = pow(10.,7.0)*0.7;
/*
    Model.g0 = 1.0;
    Model.g1 = 0.0;
    Model.g2 = 0.0;
*/

	Model.g0 = 0.57;
	Model.g1 = 0.38;
	Model.g2 = -0.01;

	Model.epsilon1 = Model.epsilon2 = 0.1;
}


void init(int mode)
{
	int i;

	init_cosmo();

    init_tabs();

    init_random_number();

    //read_catalogue(Halo_catalog_name, 0);
	read_catalogue(Halo_catalog_name, mode); 

	make_snapshot_list( Snapshot_list_name );

	TreeNHalos = malloc( sizeof(int) * NTrees );

	for(i = 0; i < NTrees; i++) TreeNHalos[i] = 0.0;
}

double fun_j(double x, void *params)
{
	double u2;

	u2 = x*x;
    return pow(1.+1./u2, 0.5*Model.g1);
}

double sigma(double m)
{
	return var_Mo(m);
	//return var_numerical(m);
}

void init_tabs(void)
{
	int i;
	double lnm_min, lnm_max, lnm, mass, dlnm;
	double m1, m2, lnm1, lnm2, lns1, lns2;
	double z_min, z_max, dz, z, z1, z2, d1, d2;
	double lnu_min, lnu_max, dlnu, u_res, result, error, gamma1;
	double sig_2, sig_res;

	gsl_integration_workspace *w
        = gsl_integration_workspace_alloc(1000);

// making table for ln(mass), ln(sigma), and alpha=-dln\sigma/dlnM
	lnm_min = log(Mass_min);
	lnm_max = log(Mass_max);
	dlnm = (lnm_max - lnm_min)/(Num_Bin_LnM-1);

	Tab_M_S.lnM = malloc(Num_Bin_LnM * sizeof(double));
	Tab_M_S.lnS = malloc(Num_Bin_LnM * sizeof(double));
	Tab_M_S.alpha = malloc(Num_Bin_LnM * sizeof(double));

	for(i=0; i<Num_Bin_LnM; i++)
	{
		Tab_M_S.lnM[i] = lnm = lnm_min + i*dlnm;
		mass = exp(lnm);
		Tab_M_S.lnS[i] = log(sigma(mass));

		lnm1 = lnm-0.01; lnm2 = lnm+0.01;
		m1 = exp(lnm1);  m2 = exp(lnm2);
		lns1 = log(sigma(m1)); lns2 = log(sigma(m2));
		Tab_M_S.alpha[i] = -(lns2-lns1)/(lnm2-lnm1);

	}

// making table for z, delta, d\delta/dz
	z_min = Redshift_min;
	z_max = Redshift_max+0.1;
	dz = (z_max - z_min)/(Num_Bin_Z - 1);

	Tab_Z_D.z = malloc(Num_Bin_Z * sizeof(double));
	Tab_Z_D.delta = malloc(Num_Bin_Z * sizeof(double));
	Tab_Z_D.dddz = malloc(Num_Bin_Z * sizeof(double));

	for(i=0; i<Num_Bin_Z; i++)
	{
		Tab_Z_D.z[i] = z = z_min + i*dz;
		Tab_Z_D.delta[i] = delta_c(z);
		z1 = z - 0.01; z2 = z+0.01;
		d1 = delta_c(z1); d2 = delta_c(z2);
		Tab_Z_D.dddz[i] = (d2-d1)/(z2-z1);
	}

// making table for J(u_res), see astro-ph/0708.1382, I notice that J(u)~u is good approximation when u>100
	gamma1 = Model.g1;
	sig_2 = sigma(Mass_max);
	Sigma_res = sig_res = sigma(Mass_res);
	lnu_min = log(sig_2/sqrt(sig_res*sig_res - sig_2*sig_2));
	lnu_max = log(1000.);
	dlnu = (lnu_max - lnu_min)/(Num_Bin_Ju - 1);

	Tab_U_J.lnures = malloc(Num_Bin_Ju * sizeof(double));
	Tab_U_J.lnJ = malloc(Num_Bin_Ju * sizeof(double));

	gsl_function F;
    F.function = &fun_j;
    F.params = &gamma1;

	for(i=0; i<Num_Bin_Ju; i++)
	{
		u_res = exp(lnu_min + i*dlnu);
//printf("u=%g\n", u_res);
		gsl_integration_qags(&F, 0.0, u_res, 0.0, 1.e-7, 1000, w, &result, &error);
		Tab_U_J.lnures[i] = log(u_res);
		Tab_U_J.lnJ[i] = log(result);
	}

	gsl_integration_workspace_free (w);

}


void make_snapshot_list( char *fname)
{
    int i;
    double redshift1, amin, amax, a, da, z;
    FILE *fd;

    Redshift_snap = malloc(Nsnapshot*sizeof(float));
    if( fd=fopen(fname, "r"))
    {
        for(i=0; i<Nsnapshot; i++)
        {
            fscanf(fd, "%lf", &a);
            z = 1./a - 1.0;
            Redshift_snap[i] = (float) z;
        }
    }
    else
    {
        if(!(fd=fopen(fname,"w")))
        {
            printf("can't open snapshot file `%s`\n", fname );
            exit(1);
        }

        redshift1 = Redshift_max;
        amin = 1./(redshift1+1); amax = 1.;
        da = pow(amax/amin, 1./(Nsnapshot-1));
        for(i=0, a=amin; i<Nsnapshot; i++)
        {
            z = 1./a -1.0;
            Redshift_snap[i] = (float) z;
            //fprintf(fd, "%d %f %f\n", i, Redshift_snap[i], hubble_time(z));
            fprintf(fd, "%f\n", (float) a);
            a *= da;
        }
    }
    fclose(fd);
}


int snapshot_number_from_z(float z)
{
	int index;
	float a, amin, amax, lnaratio;

/*
	amin = 1./(1.+Redshift_max);
    amax = 1.0;
    lnaratio = log(amax/amin);
	a = 1./(z + 1);
	index = (int) ceil((Nsnapshot - 1) * log(a/amin) / lnaratio);
*/
	index = 0;
	while(index < Nsnapshot && z < Redshift_snap[index])
	{
		index++;
	}	
	return index;
}

void free_tabs(void)
{
}
