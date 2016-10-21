#include "allvars.h"

char Halo_catalog_name[500];
char Snapshot_list_name[500];
char Tree_table_name[500];
char Simu_tree_file_name[500];
char Tree_file_name[500];
int Random_number_seed;
double Redshift;

int Nodes_On_Tree;
float N_Upper;

int Nsnapshot;
float *Redshift_snap;
int NTrees;
int *TreeNHalos;
unsigned long TotNHalos;
//unsigned long long TotNHalos;
struct tree_node *tree;
struct tree_output *reduced_tree;

float *Mvir;
float *Rvir;

int Num_Bin_LnM;
double Mass_min;
double Mass_max;
struct tab_mass_sigma Tab_M_S;

int Num_Bin_Z;
double Redshift_max;
double Redshift_min;
struct tab_z_delta Tab_Z_D;

int Num_Bin_Ju;
struct tab_u_j Tab_U_J;

double Mass_res;
double Sigma_res;

struct model_paramsters Model;

struct halo_properties Halo;

