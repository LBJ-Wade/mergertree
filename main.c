#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_sort.h>

#include "cosmo.h"
#include "allvars.h"
#include "proto.h"

int main(int argc, char *argv[])
{
	int i, mode;

	if(argc != 4)
    {
		printf("\n  usage: mergertree <catalogfile> <mode> <Random_number_seed>\n\n");
      	exit(1);
    }

	strcpy(Halo_catalog_name, argv[1]);

	mode = atoi(argv[2]);
	
	Random_number_seed = atoi(argv[3]);

	readparameters();

	init(mode);

//	test_tabs();

	write_tree_table( Tree_file_name, 0 );

	TotNHalos = 0;
    // loop over all trees
    for (i=0; i<NTrees; i++) 
	{

		printf("working on the %dth tree!\n", i+1);
		Nodes_On_Tree = 0;

        // initialize tree array & allocate memory
        init_tree(i);
//printf("OK init!\n");


        // generate tree recursively
		//Nodes_On_Tree = construct_tree_recursive(0) + 1;
		// generate tree iteratively
		Nodes_On_Tree = construct_tree_iterative(0) + 1;
//printf("OK construct!\n");

        // save the tree, or add to a histogram
//        save_snapshot( 1.0 ) ;

		reduce_tree(i);
//printf("OK reduce!\n");
		save_reduced_tree(i, Tree_file_name);
//printf("OK save!\n");
        // free memory
        free_tree( ) ;
//printf("OK free!\n");
		TotNHalos += TreeNHalos[i];
//printf("Nhalo=%d tothalo=%d\n", TreeNHalos[i], TotNHalos);
    }

	write_tree_table( Tree_file_name, 1 );
printf("after write_tree_table!\n");
//	progenitor_mass_function();

	free_random_number();

	//free_catalog();

	free_tabs();
	return 1;
}


void init_tree(int itree)
{
	int i;
    // This function initialize the tree array

    // allocate memory for the tree array
    if(!(tree = malloc(MaxNumNode*sizeof(struct tree_node))) )
    {
        printf("memory not allocated!\n");
        exit(1);
    }

	if(!(reduced_tree = malloc(MaxNumNode*sizeof(struct tree_output))) )
    {
        printf("memory not allocated!\n");
        exit(1);
    }

	for(i=0; i<MaxNumNode; i++)
	{
		reduced_tree[i].CountProgenitor = 0;
		reduced_tree[i].FirstProgenitor = -1;
		reduced_tree[i].NextProgenitor = -1;
/*
		reduced_tree[i].FirstHaloInFOFgroup = -1;
		reduced_tree[i].NextHaloInFOFgroup = -1;
*/
	}

    // initial values for the tree trunk
	Nodes_On_Tree = 0;
    tree[0].descendant = 0 ;
	tree[0].nextprogenitor = 0;
    tree[0].mass =  Mvir[itree];
    tree[0].redshift = Redshift ;

	tree[0].newid = 0; //test
}

void free_tree(void)
{
	free(tree);
	free(reduced_tree);
}


void build_node(int this_node, int ind, float redshift, float mp1, float z2)
{
	tree[this_node].mass = mp1;
    tree[this_node].redshift = z2;
    tree[this_node].dredshift = redshift;
    tree[this_node].descendant = ind;
}

int construct_tree_recursive(int ind)
{
	int ind_next;
	float mass, redshift;
	float mp1, mp2, z2;

	mass = tree[ind].mass;
	redshift = tree[ind].redshift;
	ind_next = ind;

	if(ind >= MaxNumNode)
	{
		printf("MaxNumNode (%d) has been reached!\n To run the code, make it larger!\n", MaxNumNode);
		exit(0);
	}
	if(ind >= 0)
	{
		split_node(mass, redshift, &mp1, &mp2, &z2);

		if(z2 <= Redshift_max)
		{
			if(mp1 >= 2*Mass_res)
			{
				build_node(ind+1, ind, redshift, mp1, z2);
            	ind_next = construct_tree_recursive(ind+1);
			}
			if(mp2 >= 2*Mass_res)
			{
				build_node(ind_next+1, ind, redshift, mp2, z2);
                ind_next = construct_tree_recursive(ind_next+1);
			}
		}
	}
	return ind_next;
}

int construct_tree_iterative(int ind)
{
	int totalnode;
	float mass, redshift, mp1, mp2, z2;

	totalnode = ind;
	while(totalnode >= ind)
	{
		if(totalnode >= MaxNumNode)
	    {
    	    printf("MaxNumNode (%d) has been reached!\n To run the code, make it larger!\n", MaxNumNode);
        	exit(0);
    	}

		mass = tree[ind].mass;
		redshift = tree[ind].redshift;
		split_node(mass, redshift, &mp1, &mp2, &z2);

		if(z2 <= Redshift_max)
		{
			if(mp1 >= 2*Mass_res)
			{
				totalnode++;
				build_node(totalnode, ind, redshift, mp1, z2);
			}
			if(mp2 >= 2*Mass_res)
			{
				totalnode++;
				build_node(totalnode, ind, redshift, mp2, z2);
			}
		}
		ind ++;
	}
	return totalnode;
}
		
// This function generate TWO progenitors of the halo mass: tree[ind].M
void split_node(float mass2, float redshift2, float *mp1, float *mp2, float *zp)
{
	float idle;
	double q, r_q, f, m2, m_child1, m_child2, mtemp;
	double z1, z2, dz, r1, r2, r3;
	double lnm1, lnm2, lnmh, lns1, lns2, lnsh;
	
    // prepare variables
	Halo.m2 = m2 = mass2;
	lnm2 = log(Halo.m2);
	lnmh = lnm2 - log(2.);
	lns2 = interpolate_bipoint(Tab_M_S.lnM, Tab_M_S.lnS, Num_Bin_LnM, lnm2, 1);
	lnsh = interpolate_bipoint(Tab_M_S.lnM, Tab_M_S.lnS, Num_Bin_LnM, lnmh, 1);
	Halo.sig2 = exp(lns2);
	Halo.sigh = exp(lnsh);
	Halo.alp2 = interpolate_bipoint(Tab_M_S.lnM, Tab_M_S.alpha, Num_Bin_LnM, lnm2, 1);
	Halo.alph = interpolate_bipoint(Tab_M_S.lnM, Tab_M_S.alpha, Num_Bin_LnM, lnmh, 1);

	Halo.qres = Mass_res/Halo.m2;
//printf("here:%d %g %g", ind, Halo.m2, Mass_res);
//printf("qres= %g\n", Halo.qres);
	Halo.ures = Halo.sig2/sqrt(Sigma_res*Sigma_res - Halo.sig2*Halo.sig2);
//printf("z2=%g\n", redshift2);
	Halo.z2 = z2 = redshift2;
	Halo.del2 = interpolate_bipoint(Tab_Z_D.z, Tab_Z_D.delta, Num_Bin_Z, z2, 2);
	Halo.dddz2 = interpolate_bipoint(Tab_Z_D.z, Tab_Z_D.dddz, Num_Bin_Z, z2, 2);

	params_durham( );
//printf("params: %g %g %g %g\n", Halo.B, Halo.eta, Halo.beta, Halo.mu);
//printf("params2: %g %g %g %g\n", Halo.alph, Halo.sig2, Halo.del2, Halo.sigh);

	// repeat recursive for all progenitors
	Halo.dz = timestep( Halo.qres );
    dz = Halo.dz / Halo.dddz2;
	z1 = z2 + dz;
//printf("Halo.dz=%g Halo.dddz2=%g\n", Halo.dz, Halo.dddz2);
//printf("construct_tree: qres=%g ures=%g Sigma_res=%g sig2=%g\n", Halo.qres, Halo.ures,Sigma_res,Halo.sig2);
//if(fabs(Halo.eta) <=1.e-6) printf("B=%g eta=%g beta=%g mu=%g m=%g z=%g\n", Halo.B, Halo.eta, Halo.beta, Halo.mu, m2, z2);
//printf("qres=%g dz=%g\n", Halo.qres, dz);
	f = fun_F(Halo.ures, &idle);
//printf("f=%g\n", f);
	r1 = random_number( 0 );
	if(r1 < N_Upper)
	{
		r2 = random_number( 0 );
		if(fabs(Halo.eta) > SMALLEPSILON)
			q = pow( pow(Halo.qres, Halo.eta) + (pow(2., -Halo.eta) - pow(Halo.qres, Halo.eta))*r2, 1./Halo.eta);
		else
			q = Halo.qres * pow(2.*Halo.qres, -r2);
//printf("eta=%g q=%g\n", Halo.eta, q);
		lnm1 = log(q * Halo.m2);
		lns1 = interpolate_bipoint(Tab_M_S.lnM, Tab_M_S.lnS, Num_Bin_LnM, lnm1, 1);
		Halo.sig1 = exp(lns1);
		Halo.alp1 = interpolate_bipoint(Tab_M_S.lnM, Tab_M_S.alpha, Num_Bin_LnM, lnm1, 1);
//printf("before fun_R\n");
		r_q = fun_R(q);
		if(r_q >= 1)
		{
			printf("R(q) is even greater then 1! q=%g R(q)=%g\n", q, r_q);
			exit(0);
		}

		r3 = random_number( 0 );
		if(r3 < r_q)
		{
/*
			p = 1. - q;
//printf("f=%g q=%g p=%g qres=%g\n", f, q, p, Halo.qres);
			if(p < q) p = q;
			m = m2 * (p - f);
*/
			m_child1= m2*(1.-q -f);
			m_child2= m2*q;
		}
		else
		{
			m_child1 = m2 * (1.-f);
			m_child2 = m2 * f;
		}
	}
	else
	{
		m_child1 = m2 * (1.-f);
		m_child2 = m2 * f;
	}

	if (m_child2 > m_child1)
	{
		mtemp = m_child1;
		m_child1 = m_child2;
		m_child2 = mtemp;
	}

	*mp1 = m_child1;
	*mp2 = m_child2;
	*zp = z1;
}

//sort progenitors into ascending order according to mass
void sort_progenitor_mass(int itree)
{
    int i, j, np, p, *prog_buf;
    size_t *ind;
    double *mass_buf;

    for(i=0; i<TreeNHalos[itree]; i++)
    {
        np = reduced_tree[i].CountProgenitor;
        if( np > 0)
        {
            ind = malloc(np*sizeof(size_t));
            prog_buf = malloc((np)*sizeof(int));
            mass_buf = malloc(np*sizeof(double));

            for(j=0, p = reduced_tree[i].FirstProgenitor; j<np; j++)
            {
                prog_buf[j] = p;
                mass_buf[j] = reduced_tree[p].Mvir;
//printf("**before: %d %d %f %d\n", j, p, reduced_tree[p].Mvir, reduced_tree[p].NextProgenitor);
                p = reduced_tree[p].NextProgenitor;
            }
            gsl_sort_index( ind, mass_buf, 1, np);

            p = reduced_tree[i].FirstProgenitor = prog_buf[ind[np-1]];
            for(j=0; j<np-1; j++)
            {
                reduced_tree[p].NextProgenitor = prog_buf[ind[np-2-j]];
//printf("ordered: %d %d %f %d\n", j, p, reduced_tree[p].Mvir, reduced_tree[p].NextProgenitor);
                p = reduced_tree[p].NextProgenitor;
            }
            reduced_tree[p].NextProgenitor = -1;
//printf("ordered: %d %d %f %d\n", j, p, reduced_tree[p].Mvir, reduced_tree[p].NextProgenitor);

            free(ind);
            free(prog_buf);
            free(mass_buf);
        }
    }
}

void reduce_tree( int itree )
{
	int i, j, k, icur, ides, istart;
	int descendant, iprog;
	double z, mvir, rvir, rho;

// for the root halo
	//reduced_tree[0].SnapNum = Nsnapshot-1;
	reduced_tree[0].Descendant = -1;
    mvir  = tree[0].mass;
	reduced_tree[0].Mvir = mvir/1.e10;
	z = tree[0].redshift;
	reduced_tree[0].SnapNum = istart = snapshot_number_from_z(z);
	rho = Delta_vir(z) * rho_crit(z) * 1.e-9;
	rvir = pow(mvir/(4.0*PI/3.0 * rho), 1.0/3);
    reduced_tree[0].Rvir = (float) rvir;
	tree[0].newid = 0;
	//reduced_tree[0].FirstHaloInFOFgroup = 0;
//printf("%d %d %d %g %g\n", 0, reduced_tree[0].SnapNum, reduced_tree[0].Descendant, reduced_tree[0].Mvir, tree[0].redshift);

// go up the tree
	for(i = 1, j=1; i < Nodes_On_Tree; i++)
	{
		icur = snapshot_number_from_z(tree[i].redshift);
		ides = snapshot_number_from_z(tree[tree[i].descendant].redshift);
		if(icur > 0 && icur < istart && ides - icur > 0) // halos in the first snapshot are not saved for the offset of 1 in the zlist
		{
			tree[i].newid = j;
			reduced_tree[j].SnapNum = icur;
			mvir = tree[i].mass;
			reduced_tree[j].Mvir = mvir/1.e10;
			z = tree[i].redshift;
			rho = Delta_vir(z) * rho_crit(z) * 1.e-9;
			rvir = pow(mvir/(4.0*PI/3.0 * rho), 1.0/3);
			reduced_tree[j].Rvir = (float) rvir;
			k = tree[i].descendant;
			reduced_tree[j].Descendant = descendant = tree[k].newid;
/*
	// circled nextprogenitor pointers, pointer initializer needs to be changed too
			if(reduced_tree[descendant].CountProgenitor == 0) {reduced_tree[descendant].FirstProgenitor = j;}
			for(pcount=1, iprog=reduced_tree[descendant].FirstProgenitor; pcount<reduced_tree[descendant].CountProgenitor; 
								pcount++, iprog = reduced_tree[iprog].NextProgenitor);
			reduced_tree[iprog].NextProgenitor = j;
			reduced_tree[j].NextProgenitor = reduced_tree[descendant].FirstProgenitor;
	//
*/
	// straight nextprogenitor pointers, pointer initializer needs to be changed too
			if(reduced_tree[descendant].CountProgenitor < 1)
			{
				reduced_tree[descendant].FirstProgenitor = j;
			}
			else
			{
				iprog = reduced_tree[descendant].FirstProgenitor;
				while(reduced_tree[iprog].NextProgenitor >= 0)
				{
					iprog = reduced_tree[iprog].NextProgenitor;
				}
				reduced_tree[iprog].NextProgenitor = j;
			}
	//
			reduced_tree[descendant].CountProgenitor++ ;

			//reduced_tree[j].FirstHaloInFOFgroup = j;
//if(descendant==2) printf("i=%d j=%d\n",i, j);
//printf("i=%d scur=%f sdes=%f %d %d\n", j, scur, sdes, icur, ides);
//printf("%d %d %d %g %g\n", j, reduced_tree[j].SnapNum, reduced_tree[j].Descendant, reduced_tree[j].Mvir, tree[i].redshift);
			j++;
		}
		else
		{
			tree[i].newid = tree[tree[i].descendant].newid;
		}
	}

	TreeNHalos[itree] = j;  // as the total number of halos in the reduced tree

	//sort progenitors into ascending order according to mass
	sort_progenitor_mass(itree);
}
