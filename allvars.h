//#define int unsigned long
#define PI 3.14159265358979
#define SMALLEPSILON  1.e-6

#define MaxNumNode 90000000

extern char Halo_catalog_name[500];
extern char Snapshot_list_name[500];
extern char Tree_table_name[500];
extern char Simu_tree_file_name[500];
extern char Tree_file_name[500];
extern int Random_number_seed;
extern double Redshift;

extern int Nodes_On_Tree;
extern float N_Upper;

extern int Nsnapshot;
extern float *Redshift_snap;
extern int NTrees;
extern int *TreeNHalos;
extern unsigned long TotNHalos;

extern struct tree_node
{
	int descendant;
	int firstprogenitor;
    int nextprogenitor;
	int newid; // id in reduced tree
    float mass;
    float redshift;
	float dredshift;
} *tree;

extern struct tree_output
{
	int SnapNum;
	int Descendant;
	int CountProgenitor;
	int FirstProgenitor;
	int NextProgenitor;
/*
	int FirstHaloInFOFgroup;
	int NextHaloInFOFgroup;
*/
	float Mvir;
	float Rvir;
} *reduced_tree;

extern float *Mvir;
extern float *Rvir;

extern int Num_Bin_LnM;
extern double Mass_min;
extern double Mass_max;
extern struct tab_mass_sigma
{
    double *lnM, *lnS, *alpha;
} Tab_M_S;


extern int Num_Bin_Z;
extern double Redshift_max;
extern double Redshift_min;
extern struct tab_z_delta
{
	double *z, *delta, *dddz;
} Tab_Z_D;


extern int Num_Bin_Ju;
extern struct tab_u_j 
{
	double *lnures, *lnJ;
} Tab_U_J;

extern double Mass_res;
extern double Sigma_res;

extern struct model_paramsters
{
	float g0, g1, g2;
	float epsilon1, epsilon2;
} Model;

extern struct halo_properties 
{
	double m2;
	double sig2, sigh, sig1;
	double alp2, alph, alp1;
	double z2, del2, dddz2, dz;
	double ures, qres;
	double B, eta, beta, mu;
} Halo;



