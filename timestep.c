#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#include "allvars.h"
#include "proto.h"


double timestep(double qres)
// maximum timestep needs to be considered!!!!
{
	int isnap;
	double result, dz1, dz;

	if(fabs(Halo.eta) > SMALLEPSILON)
		result = sqrt(2./PI)*Halo.alph*Halo.B/Halo.eta*(pow(2., -Halo.eta) - pow(Halo.qres, Halo.eta)) 
				* Model.g0 * pow(Halo.del2/Halo.sig2, Model.g2) * pow(Halo.sigh/Halo.sig2, Model.g1);
	else 
		result = -sqrt(2./PI)*Halo.alph*Halo.B*log(2.*Halo.qres) 
                * Model.g0 * pow(Halo.del2/Halo.sig2, Model.g2) * pow(Halo.sigh/Halo.sig2, Model.g1);
	dz = Model.epsilon2/result;

	dz1 = Model.epsilon1*sqrt(2.* (Halo.sigh * Halo.sigh - Halo.sig2 * Halo.sig2));
//printf("dz=%g dz1=%g eta=%g result=%g the#=%g\n", dz, dz1, Halo.eta, result, pow(2., -Halo.eta) - pow(Halo.qres, Halo.eta));
	if(dz1 < dz) dz = dz1;

	isnap = snapshot_number_from_z(Halo.z2);
	if(isnap > 0)
	{
		dz1 = Redshift_snap[isnap-1] - Redshift_snap[isnap];
		dz1 *= Halo.dddz2;
		if(dz1 < dz) dz = dz1;
	}

	N_Upper = dz * result;
//printf("timestep: %g\n", pow(Halo.sigh/Halo.sig2, Model.g1));
	return dz;
}
