#include "BB/std.h"
#include "BB/calc.h"
#include "BB/vtk.h"
#include "BBFE/sys/read.h"
#include "BBFE/sys/memory.h"
#include "BBFE/sys/write.h"
#include "BBFE/sys/FE_dataset.h"
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

static const char* CODENAME = "inflow_sups_poiseuille >";
static const char* VOIDNAME = "                       >";
static const int BUFFER_SIZE = 10000;

const char* ID_PENALTY_VELOCITY_2D_POISEUILLE     = "#velocity_2d_poiseuille";

static const char* DEF_DIRECTORY     = ".";
static const char* DEF_FILENAME_NODE = "node.dat";
static const char* DEF_FILENAME_D_BC  = "D_bc.dat";
 
static const char* OPTION_DIRECTORY           = "-d";
static const char* OPTION_INFILE_NODE         = "-in";
static const char* OPTION_OUTFILE_D_BC        = "-od";


typedef struct {
    const char* infile_node;
    const char* outfile_dbc;
    const char* directory;
} SETTINGS;

typedef struct {
    double x0;  double y0;  double z0;
    double l_x; double l_y; double l_z;
    double zc;   double radius;
} CONDITION;

void args_manager(
    SETTINGS*   set,
    int         argc,
    char*       argv[])
{
    if(argc < 8) {
        printf("%s ERROR: Please specify parameters for 2d_poiseuille\n", CODENAME);
        printf("%s Format: \n", VOIDNAME);
        printf("%s ./inflow_sups_poiseuille [1: refarence point x(x0)][2: refarence point(y0)][3: refarence point(z0)][4: global x length(l_x)][5: global y length(l_y)][6: global z length(l_z)][7: center of hight(z_c)] [8: radius (radius)] \n\n", VOIDNAME);
        printf("%s -d  [directory] (Default: ./) \n", VOIDNAME);
        printf("%s -in [input nodes] (Default: node.dat) \n", VOIDNAME);
        printf("%s -od [output D_bc] (Default: D_bc.dat) \n", VOIDNAME);

        exit(0);
    }

    int num;
    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_DIRECTORY);
    if(num == -1) {
        set->directory = DEF_DIRECTORY;
    }
    else {
        set->directory = argv[num+1];
    }

    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_INFILE_NODE);
    if(num == -1) {
        set->infile_node = DEF_FILENAME_NODE;
    }
    else {
        set->infile_node = argv[num+1];
    }

    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_OUTFILE_D_BC);
    if(num == -1) {
        set->outfile_dbc = DEF_FILENAME_D_BC;
    }
    else {
        set->outfile_dbc = argv[num+1];
    }

    printf("%s Input/output directory   : %s\n", CODENAME, set->directory);
    printf("%s Input filename (node)    : %s\n", CODENAME, set->infile_node);
    printf("%s Output filename (D_bc)   : %s\n", CODENAME, set->outfile_dbc);
}

void set_condition(
        char* argv[],
        CONDITION* cond)
{
    cond->x0      = atof( argv[1] );
    cond->y0      = atof( argv[2] );
    cond->z0      = atof( argv[3] );
    cond->l_x     = atof( argv[4] );
    cond->l_y     = atof( argv[5] );
    cond->l_z     = atof( argv[6] );
    cond->zc      = atof( argv[7] );
    cond->radius  = atof( argv[8] );

    printf("%s Reference point (%e %e %e)\n", CODENAME,
            cond->x0, cond->y0 , cond->z0);
    printf("%s xy-plane area length (%e %e %e)\n", CODENAME,
            cond->l_x, cond->l_y , cond->l_z);
    printf("%s hight of center :%e\n", CODENAME,
            cond->zc);
    printf("%s radius  :%e\n", CODENAME, cond->radius);
}


static bool node_is_inside_plane(
        const double 	x,
        const double 	val,
        const double 	epsilon) 
{
    double min_val = val - epsilon;
    double max_val = val + epsilon;

    if ( min_val < x && x < max_val ) {
        return true;
    }
    else {
        return false;
    }
}


int detect_Dirichlet_bc(
        BBFE_DATA*      fe,
        CONDITION*      cond)
{
    double x_min = cond->x0;
    double y_min = cond->y0;
    double z_min = cond->z0;
    // double x_max = cond->x0 + cond->l_x;
    double y_max = cond->y0 + cond->l_y;
    double z_max = cond->z0 + cond->l_z;

    double epsilon = 1.0e-05;

    // Count the num. of Dirichlet B.C.
    int num_bc_nodes = 0;
    for(int i=0; i<(fe->total_num_nodes); i++) {
        // inflow
        if( node_is_inside_plane(fe->x[i][0], x_min, epsilon) ) { num_bc_nodes += 3; }
              
        // non-slip
        else if( node_is_inside_plane(fe->x[i][2], z_max, epsilon) ||
            node_is_inside_plane(fe->x[i][2], z_min, epsilon) ) { num_bc_nodes += 3; }
            
        // slip
        else if( node_is_inside_plane(fe->x[i][1], y_min, epsilon) ||
            node_is_inside_plane(fe->x[i][1], y_max, epsilon) ) { num_bc_nodes += 1; }

    }

    return num_bc_nodes;
}


double calc_velocity_theory(
    double zc,
    double z
)
{
    double val;
    double dz;
    double tmp;

    dz = z - zc;
    tmp = dz*dz;

    val = 1-4*tmp;

    return val;
}


void output_Dirichlet_bc_2d_poiseuille(
    BBFE_DATA* fe,
    SETTINGS*  set,
    CONDITION* cond,
    const int  num_bc_nodes)
{
    double x_min = cond->x0;
    double y_min = cond->y0;
    double z_min = cond->z0;
    // double x_max = cond->x0 + cond->l_x;
    double y_max = cond->y0 + cond->l_y;
    double z_max = cond->z0 + cond->l_z;

    double epsilon = 1.0e-05;

    // Outout Dirichlet B.C. file
    char fname_dbc[BUFFER_SIZE];
    snprintf(fname_dbc, BUFFER_SIZE, "%s/%s", set->directory, set->outfile_dbc);
    FILE* fp_dbc;
    fp_dbc = fopen(fname_dbc, "w");
    if( fp_dbc == NULL ) {
        printf("%s ERROR: File \"%s\" cannot be opened.\n", 
                CODENAME, fname_dbc);
    }

    int block_length = 4;

    fprintf(fp_dbc, "%d %d\n", num_bc_nodes, block_length);
    for(int i=0; i<(fe->total_num_nodes); i++) {
        // 流入境界条件
        if( node_is_inside_plane(fe->x[i][0], x_min, epsilon) ){
            double vx = 0.0;
            vx = calc_velocity_theory(cond->zc, fe->x[i][2]);
            fprintf(fp_dbc, "%d 0 %lf \n", i, vx);
            fprintf(fp_dbc, "%d 1 0.0\n", i);
            fprintf(fp_dbc, "%d 2 0.0\n", i);
        }

        // non-slip
        else if( node_is_inside_plane(fe->x[i][2], z_min, epsilon) ||
            node_is_inside_plane(fe->x[i][2], z_max, epsilon) ) {
                fprintf(fp_dbc, "%d 0 0.0\n", i);
                fprintf(fp_dbc, "%d 1 0.0\n", i);
                fprintf(fp_dbc, "%d 2 0.0\n", i);
            }

        // slip
        else if( node_is_inside_plane(fe->x[i][1], y_min, epsilon) ||
            node_is_inside_plane(fe->x[i][1], y_max, epsilon) ) {
                fprintf(fp_dbc, "%d 1 0.0\n", i);
            }
    }
    fclose(fp_dbc);
}



int main(
    int     argc,
    char*   argv[])
{
    printf("\n");

    BBFE_DATA	    fe;
    SETTINGS        set;
    CONDITION       cond;

    //BBFE_sys_read_void();
    BB_vtk_void();

    args_manager(&set, argc, argv);
    set_condition(argv, &cond);

    BBFE_sys_read_node(&fe, set.infile_node, set.directory);

    printf("%s start write D_bc file.\n", CODENAME);
    int num_bc_nodes = detect_Dirichlet_bc(&fe, &cond);
    output_Dirichlet_bc_2d_poiseuille(&fe, &set, &cond, num_bc_nodes);
    printf("%s Finish writing Dirichlet b.c. file.\n", CODENAME);

    return 0;
}