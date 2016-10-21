#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "cosmo.h"
#include "proto.h"

double fun_F(double ures, void *params)
{
	double d2, s2, dddz, dz;
	double x, y, j, f;

	d2 = Halo.del2;
	s2 = Halo.sig2;
	dddz = Halo.dddz2;
	dz = Halo.dz;
	//dw = Halo.dw;

	x = log(ures);
//printf("fun_F: ures=%g, x=%g\n", ures, x);
	y = interpolate_bipoint(Tab_U_J.lnures, Tab_U_J.lnJ, Num_Bin_Ju, x, 4);
	j = exp(y);

	f = j * Model.g0 * pow(d2/s2, Model.g2) * dz;

	f /= s2;
	f *= sqrt(2./PI); // corrected by luyu !!!!!

//printf("fun_F: ures=%g, x=%g y=%g j=%g dz=%g f=%g\n", ures, x, y, j, dz, f);
	return f;
}


double fun_R(double q)
{
	double b, beta;
	double t2, t3, t4;

	b = Halo.B;
	beta = Halo.beta;

	t2 = Halo.alp1/Halo.alph;
	t3 = fun_V(1)/b/pow(q, beta);
	t4 = pow(pow(q, Halo.mu)*Halo.sig1/Halo.sigh, Model.g1);
//printf("alph1=%g alph=%g %g %g %g\n", Halo.alp1, Halo.alph, t2, t3, t4);
	return t2*t3*t4;
}

double fun_V(int mode)
{
	double t1, t2;

	switch(mode)
	{
		case 1:
			t1 = Halo.sig1*Halo.sig1;
			break;
		case 2:
			t1 = Sigma_res*Sigma_res;
			break;
		case 3:
			t1 = Halo.sigh*Halo.sigh;
			break;
		default:
			printf("You did not choose an appropriate mode for fun_V(q)!\n");
			exit(3);
	}
	t2 = t1 - Halo.sig2*Halo.sig2;
//printf("fun_V: t1=%g t2=%g %g %g\n", t1, t2, Halo.sig1, Halo.sig2);
	return t1/(t2*sqrt(t2));
}


void params_durham( void )
{
	double t1, t2, t3, t4, lnb;

	if(Model.g1 > 0.0 ) Halo.mu = Halo.alph;
	else
	{
		t1 = log(Sigma_res/Halo.sigh);
		t2 = log(2.*Halo.qres);
		Halo.mu = -t1/t2;
	}

	t1 = log(fun_V(2)); t2 = log(fun_V(3));
	t3 = log(Halo.qres); t4 = log(0.5);
	Halo.beta = (t1 - t2)/(t3 - t4);

	lnb = (t2*t3 - t1*t4)/(t3 - t4);
	//Halo.B = exp(lnb);
	Halo.B = fun_V(3) * pow(2., Halo.beta);

	Halo.eta = Halo.beta - 1. - Model.g1* Halo.mu;
//printf("params_durham: t1=%g t2=%g t3=%g t4=%g\n", *b, *eta, *beta, *mu);
}


void progenitor_mass_function(void)
{
    int i, j;
    double xmass0,xmass, lnm1, lnm2, lns1, lns2;
    double xlogm, sig, sig1, sig2, ds, dc1, dc2, dww, fac, m1, m2, dsdm, xn;
    double fac_sim, alpha;
	float fMass_res = 1./10.0;
    float zwrite[7];
    FILE *fz[7];

    fz[0]=fopen("r.2.dat","w");
    fz[1]=fopen("r.5.dat","w");
    fz[2]=fopen("r1.dat","w");
    fz[3]=fopen("r2.dat","w");
    fz[4]=fopen("r3.dat","w");
    fz[5]=fopen("r4.dat","w");
    fz[6]=fopen("r6.dat","w");

    zwrite[0] = 0.2; zwrite[1] = 0.5; zwrite[2] = 1.0; zwrite[3] = 2.0; zwrite[4]= 3.0; zwrite[5]=4.0; zwrite[6]=6.0;

    xmass0 = 1.e12; //fMass_res * 1.e12;

	for(j=0; j<7; j++)
	{
    	for(i=0;i<1000;i++)
    	{
        	xlogm = log10(fMass_res*0.01) + (double) i /999.9 * (2.-log10(fMass_res));
        	xlogm = log10(xmass0) + xlogm;
        	xmass = pow(10.,xlogm);

			m1 = xmass;
            m2 = xmass0;

            lnm1 = log(m1);
            lnm2 = log(m2);
            lns1 = interpolate_bipoint(Tab_M_S.lnM, Tab_M_S.lnS, Num_Bin_LnM, lnm1, 1);
            lns2 = interpolate_bipoint(Tab_M_S.lnM, Tab_M_S.lnS, Num_Bin_LnM, lnm2, 1);
            sig1 = sig = exp(lns1);
            sig2 = exp(lns2);
            sig1 = sig1*sig1;
            sig2 = sig2*sig2;

        	ds = sig1 - sig2;

        	dc1 = delta_c(0.0);
        	dc2 = delta_c(zwrite[j]);
        	dww = dc2 - dc1;

        	fac = (1./sqrt(2.*M_PI)) * (dww/pow(ds,1.5) * exp(-dww*dww/(2.*ds)));

			alpha = interpolate_bipoint(Tab_M_S.lnM, Tab_M_S.alpha, Num_Bin_LnM, lnm1, 3);

			dsdm = 2.*sig*sig/m1*alpha;
        	xn = (xmass0/xmass) * fac * dsdm;
fac_sim = 0.4*pow(dww/sqrt(ds),0.75)*exp(-pow(dww/sqrt(ds),3)/10.)*dsdm/ds*0.5;
//fac_sim = sqrt(2./PI)*dww/sqrt(ds)*exp(-pow(dww/sqrt(ds),2)/2.)*dsdm/ds*0.5;
        	fprintf(fz[j],"%g %g %g %g\n",xlogm,dsdm,fac_sim,xn);
		}
    }
	for(j=0;j<7;j++) fclose(fz[j]);
}


