// main.c
void init(int mode);
void init_tree(int itree);
void free_tree(void);
//void construct_tree( struct tree_node node );
int construct_tree_recursive(int ind);
int construct_tree_iterative(int ind);
//void construct_tree( int ind, double m2, double z2 );
int construct_tree02(double m2, double z2);
void split_node(float mass2, float redshift2, float *mp1, float *mp2, float *zp);
void reduce_tree( int itree );


//funcs_durham.c
double fun_F(double ures, void *params);
double fun_S(double q, void *params);
double fun_R(double q);
double fun_V(int mode);
void params_durham( void );
void progenitor_mass_function( void );

//init.c
void readparameters();
double fun_j(double x, void *params);
double sigma(double m);
void init_tabs(void);
void make_snapshot_list( char *fname);
int snapshot_number_from_z(float z);
void free_tabs(void);

//misc.c
void init_random_number(void);
void free_random_number(void);
double random_number(int flag);
double interpolate_bipoint(double *x, double *y, int n, double x0, int flag);

//timestep.c
double timestep(double qres);

// io.c
void read_catalogue(char *fname, int mode);
void write_tree_table( char *fname, int mode);
void save_reduced_tree( int itree, char *fname );
void save_snapshot(double z);

//test.c
void test_tabs(void);
