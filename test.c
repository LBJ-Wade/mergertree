#include <stdio.h>
#include <math.h>

#include "allvars.h"

void test_tabs(void)
{
	int i;

	printf("I am in test_tabs!\n");

	for(i=0; i<Num_Bin_LnM; i++)
	{
//		printf("%d %g %g %g\n", i, Tab_M_S.lnM[i], Tab_M_S.lnS[i], Tab_M_S.alpha[i]);
		printf("%g %g\n", log10(exp(Tab_M_S.lnM[i])), exp(Tab_M_S.lnS[i]));
	}

/*

	for(i=0; i<Num_Bin_Z; i++)
	{
		printf("%d %g %g %g\n", i, Tab_Z_D.z[i], Tab_Z_D.delta[i], Tab_Z_D.dddz[i]);
	}


	for(i=0; i<Num_Bin_Ju; i++)
	{
		printf("%d %g %g\n", i, Tab_U_J.lnures[i], Tab_U_J.lnJ[i]);
	}
*/
}
