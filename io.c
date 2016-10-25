#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

int compare (const void * a, const void * b)
{
	int c;

	if((*(float*)b) - (*(float*)a) > 0.0) c = 1;
	if((*(float*)b) - (*(float*)a) == 0.0) c = 0;
	if((*(float*)b) - (*(float*)a) < 0.0) c = -1;
  	return c;
}


void read_catalogue(char *fname, int mode)
{
	FILE  *fd;
  	int   i, j, check, nsub;
	long int libuf;
  	float z, fbuf;

	if(!(fd=fopen(fname,"r")))
    {
      	printf("can't open file `%s`\n", fname );
      	exit(0);
    }

    fscanf(fd, "%d", &Nsnapshot);
    fscanf(fd, "%s", Snapshot_list_name);
    printf("%d snapshots in the Snapshot_alist, %s.\n", Nsnapshot, Snapshot_list_name);

    fscanf(fd, "%s", Simu_tree_file_name);
    printf("Simulation trees are in %s.\n", Simu_tree_file_name);

    fscanf(fd, "%s", Tree_file_name);
    printf("Trees will be written in the file, %s.\n", Tree_file_name);

    if(mode == 0 || mode == 1) fscanf(fd, "%d", &NTrees);
    else fread(&NTrees, sizeof(int), 1, fd);
    printf("The catalog contains %d halo(s)!\n", NTrees);

    if(mode==2) fread(&z, sizeof(float), 1, fd); // gadget postprocessing

    Mvir = malloc(NTrees*sizeof(float));
    Rvir = malloc(NTrees*sizeof(float));
    //Fpos = malloc(NTrees*sizeof(long int));

	switch(mode)
	{
		case 0:
        	for(i=0; i<NTrees; i++)
        	{
            	//fscanf(fd, "%e", &Mvir[i]);
				fscanf(fd, "%f %ld %d", &Mvir[i], &libuf, &nsub);

				//Mvir[i] *= 1.e10;
            	//printf("%f\n", Mvir[i]);
        	}
			break;

		case 1:
			fscanf(fd, "%f %ld %d", &fbuf, &libuf, &nsub);
			for(i=0; i<NTrees; i++)
			{
				Mvir[i] = fbuf;
			}
			break;

		case 2:
        	for(i=0; i<NTrees; i++)
        	{
     	//       fscanf(fd, "%f %ld", &Mvir[i], &Fpos[i]);
				fscanf(fd, "%f %ld %d", &Mvir[i], &libuf, &nsub);
				for(j=0; j<nsub; j++)    fscanf(fd, "%ld ", &libuf);
        	}
			break;

		case 3:
        	fread(&Mvir[0], sizeof(float), NTrees, fd);
        	fread(&Rvir[0], sizeof(float), NTrees, fd);
			break;
		default:
            printf("give an appropriate choice!\n");
            exit(2);
    }

    fclose(fd);

//	qsort(Mvir, NTrees, sizeof(float), compare);

	printf("%d\n", NTrees);
	for(i=0; i<NTrees; i++) 
	{
		printf("%e\n", Mvir[i]);
	}
}


void write_tree_table( char *fname, int mode)
{
	int i;
	FILE *fd;

	if(mode)
	{
		if(!(fd=fopen(fname,"r+b")))
    	{
			printf("can't open file `%s`\n", fname );
        	exit(1);
    	}
	}
	else
	{
		if(!(fd=fopen(fname,"w")))
        {
            printf("can't open file `%s`\n", fname );
            exit(1);
        }
	}
/*
	printf("%d\n", NTrees);
    printf("%d\n", TotNHalos);
    for(i=0; i<NTrees; i++)
    {
        printf("%d \n", TreeNHalos[i]);
    }
*/
	fwrite(&NTrees, 1, sizeof(int), fd);
	fwrite(&TotNHalos, 1, sizeof(unsigned long), fd);
	fwrite(TreeNHalos, NTrees, sizeof(int), fd);

	fclose(fd);
}


void save_reduced_tree( int itree, char *fname )
{
	int i;
	FILE *fd;

	if(!(fd=fopen(fname,"ab")))
    {
        printf("can't open file `%s`\n", fname );
        exit(1);
    }

	for(i=0; i<TreeNHalos[itree]; i++)
	{
		//printf("%d %d %d %d %d %d %g\n", i, reduced_tree[i].SnapNum, reduced_tree[i].Descendant, reduced_tree[i].CountProgenitor, reduced_tree[i].FirstProgenitor, reduced_tree[i].NextProgenitor, reduced_tree[i].Mvir);

		fwrite(&reduced_tree[i], sizeof(struct tree_output), 1, fd);
    }

	fclose(fd);
}


void save_snapshot(double z)
{
    int i,j;
	float zsave[7];
	FILE *fz[7];

    fz[0]=fopen("hz.2.dat","a");
    fz[1]=fopen("hz.5.dat","a");
    fz[2]=fopen("hz1.dat","a");
    fz[3]=fopen("hz2.dat","a");
    fz[4]=fopen("hz3.dat","a");
    fz[5]=fopen("hz4.dat","a");
    fz[6]=fopen("hz6.dat","a");

    zsave[0] = 0.2; zsave[1] = 0.5; zsave[2] = 1.0; zsave[3] = 2.0; zsave[4]= 3.0; zsave[5]=4.0; zsave[6]=6.0;
    for(i=0;i<Nodes_On_Tree;i++)
    {
            for(j=0;j<7;j++)
            {
                    if(tree[i].redshift >= zsave[j] && tree[i].dredshift < zsave[j])
                    {
                        fprintf(fz[j],"%ld %g %g %d %d\n", i, 1+ tree[i].redshift, tree[i].mass, 1, 1);
                    }
            }
    }

	for(i=0;i<7;i++) fclose(fz[i]);
}
