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
#include <math.h>

static const char* CODENAME = "inflow_sups_pipe_poiseuille >";
static const char* VOIDNAME = "                            >";
static const int BUFFER_SIZE = 10000;

const char* ID_PENALTY_VELOCITY_3D_PIPE_POISEUILLE = "#velocity_3d_pipe_poiseuille";

static const char* DEF_DIRECTORY      = ".";
static const char* DEF_FILENAME_NODE  = "node.dat";
static const char* DEF_FILENAME_D_BC  = "D_bc.dat";

static const char* OPTION_DIRECTORY    = "-d";
static const char* OPTION_INFILE_NODE  = "-in";
static const char* OPTION_OUTFILE_D_BC = "-od";


typedef struct {
    const char* infile_node;
    const char* outfile_dbc;
    const char* directory;
} SETTINGS;


typedef struct {
    double x0;      // inflow plane position
    double l_x;     // pipe length
    double yc;      // center of pipe in y direction
    double zc;      // center of pipe in z direction
    double radius;  // pipe radius
    double umax;    // maximum velocity at pipe center
} CONDITION;


void args_manager(
    SETTINGS* set,
    int       argc,
    char*     argv[])
{
    /*
        Required arguments:
        argv[1] = x0
        argv[2] = l_x
        argv[3] = yc
        argv[4] = zc
        argv[5] = radius
        argv[6] = umax

        Options should be written after the required arguments.
    */
    if(argc < 7) {
        printf("%s ERROR: Please specify parameters for 3d_pipe_poiseuille\n", CODENAME);
        printf("%s Format:\n", VOIDNAME);
        printf("%s ./inflow_sups_poiseuille "
               "[1: inflow x-plane (x0)] "
               "[2: pipe length (l_x)] "
               "[3: center y (yc)] "
               "[4: center z (zc)] "
               "[5: radius] "
               "[6: maximum velocity (umax)]\n\n", VOIDNAME);

        printf("%s Options:\n", VOIDNAME);
        printf("%s -d  [directory]    Default: ./\n", VOIDNAME);
        printf("%s -in [input nodes]  Default: node.dat\n", VOIDNAME);
        printf("%s -od [output D_bc]  Default: D_bc.dat\n", VOIDNAME);
        printf("\n");

        printf("%s Example:\n", VOIDNAME);
        printf("%s ./inflow_sups_poiseuille 0.0 10.0 0.0 0.0 0.5 1.0 "
               "-d ./mesh -in node.dat -od D_bc.dat\n", VOIDNAME);

        exit(0);
    }

    int num;

    num = BB_std_read_args_return_char_num(argc, argv, OPTION_DIRECTORY);
    if(num == -1) {
        set->directory = DEF_DIRECTORY;
    }
    else {
        set->directory = argv[num + 1];
    }

    num = BB_std_read_args_return_char_num(argc, argv, OPTION_INFILE_NODE);
    if(num == -1) {
        set->infile_node = DEF_FILENAME_NODE;
    }
    else {
        set->infile_node = argv[num + 1];
    }

    num = BB_std_read_args_return_char_num(argc, argv, OPTION_OUTFILE_D_BC);
    if(num == -1) {
        set->outfile_dbc = DEF_FILENAME_D_BC;
    }
    else {
        set->outfile_dbc = argv[num + 1];
    }

    printf("%s Input/output directory : %s\n", CODENAME, set->directory);
    printf("%s Input filename node    : %s\n", CODENAME, set->infile_node);
    printf("%s Output filename D_bc   : %s\n", CODENAME, set->outfile_dbc);
}


void set_condition(
    char*      argv[],
    CONDITION* cond)
{
    cond->x0     = atof(argv[1]);
    cond->l_x    = atof(argv[2]);
    cond->yc     = atof(argv[3]);
    cond->zc     = atof(argv[4]);
    cond->radius = atof(argv[5]);
    cond->umax   = atof(argv[6]);

    if(cond->radius <= 0.0) {
        printf("%s ERROR: radius must be positive.\n", CODENAME);
        exit(1);
    }

    if(cond->l_x <= 0.0) {
        printf("%s ERROR: l_x must be positive.\n", CODENAME);
        exit(1);
    }

    printf("%s Inflow x-plane : %e\n", CODENAME, cond->x0);
    printf("%s Pipe length    : %e\n", CODENAME, cond->l_x);
    printf("%s Center y,z     : %e %e\n", CODENAME, cond->yc, cond->zc);
    printf("%s Radius         : %e\n", CODENAME, cond->radius);
    printf("%s Umax           : %e\n", CODENAME, cond->umax);
}


static bool node_is_inside_plane(
    const double x,
    const double val,
    const double epsilon)
{
    if(fabs(x - val) <= epsilon) {
        return true;
    }
    else {
        return false;
    }
}


static bool node_is_inside_range(
    const double x,
    const double min_val,
    const double max_val,
    const double epsilon)
{
    if((min_val - epsilon <= x) && (x <= max_val + epsilon)) {
        return true;
    }
    else {
        return false;
    }
}


static bool node_is_inside_circle_yz(
    const double y,
    const double z,
    const double yc,
    const double zc,
    const double radius,
    const double epsilon)
{
    double dy = y - yc;
    double dz = z - zc;

    double r2 = dy * dy + dz * dz;

    double r_max = radius + epsilon;

    if(r2 <= r_max * r_max) {
        return true;
    }
    else {
        return false;
    }
}


static bool node_is_on_cylinder_wall(
    const double y,
    const double z,
    const double yc,
    const double zc,
    const double radius,
    const double epsilon)
{
    double dy = y - yc;
    double dz = z - zc;

    double r2 = dy * dy + dz * dz;
    double R2 = radius * radius;

    /*
        判定条件:
            r = R

        sqrt を使わずに r^2 で判定する。
        |r^2 - R^2| <= 2 R eps + eps^2
    */
    double tolerance = 2.0 * radius * epsilon + epsilon * epsilon;

    if(fabs(r2 - R2) <= tolerance) {
        return true;
    }
    else {
        return false;
    }
}


double calc_velocity_pipe_poiseuille(
    const double y,
    const double z,
    const double yc,
    const double zc,
    const double radius,
    const double umax)
{
    double dy = y - yc;
    double dz = z - zc;

    double r2 = dy * dy + dz * dz;
    double R2 = radius * radius;

    /*
        3D circular pipe Poiseuille profile:

            u_x(r) = umax * (1 - r^2 / R^2)
            u_y    = 0
            u_z    = 0
    */
    double val = umax * (1.0 - r2 / R2);

    /*
        丸め誤差で壁面近傍が微小負になる場合を防ぐ
    */
    if((val < 0.0) && (val > -1.0e-10)) {
        val = 0.0;
    }

    return val;
}


int detect_Dirichlet_bc(
    BBFE_DATA* fe,
    CONDITION* cond)
{
    double x_min = cond->x0;
    double x_max = cond->x0 + cond->l_x;

    double yc = cond->yc;
    double zc = cond->zc;
    double R  = cond->radius;

    double epsilon = 1.0e-05;

    int num_bc_nodes = 0;

    for(int i = 0; i < fe->total_num_nodes; i++) {
        double x = fe->x[i][0];
        double y = fe->x[i][1];
        double z = fe->x[i][2];

        /*
            inflow:
                x = x0 かつ yz断面内の円内部
        */
        bool is_inflow =
            node_is_inside_plane(x, x_min, epsilon) &&
            node_is_inside_circle_yz(y, z, yc, zc, R, epsilon);

        /*
            cylinder wall:
                x方向は管長範囲内
                yz断面では r = R
        */
        bool is_wall =
            node_is_inside_range(x, x_min, x_max, epsilon) &&
            node_is_on_cylinder_wall(y, z, yc, zc, R, epsilon);

        if(is_inflow) {
            num_bc_nodes += 3;
        }
        else if(is_wall) {
            num_bc_nodes += 3;
        }
    }

    return num_bc_nodes;
}


void output_Dirichlet_bc_3d_pipe_poiseuille(
    BBFE_DATA* fe,
    SETTINGS*  set,
    CONDITION* cond,
    const int  num_bc_nodes)
{
    double x_min = cond->x0;
    double x_max = cond->x0 + cond->l_x;

    double yc   = cond->yc;
    double zc   = cond->zc;
    double R    = cond->radius;
    double umax = cond->umax;

    double epsilon = 1.0e-05;

    char fname_dbc[BUFFER_SIZE];
    snprintf(fname_dbc, BUFFER_SIZE, "%s/%s", set->directory, set->outfile_dbc);

    FILE* fp_dbc = fopen(fname_dbc, "w");
    if(fp_dbc == NULL) {
        printf("%s ERROR: File \"%s\" cannot be opened.\n",
                CODENAME, fname_dbc);
        exit(1);
    }

    int block_length = 4;

    fprintf(fp_dbc, "%d %d\n", num_bc_nodes, block_length);

    for(int i = 0; i < fe->total_num_nodes; i++) {
        double x = fe->x[i][0];
        double y = fe->x[i][1];
        double z = fe->x[i][2];

        bool is_inflow =
            node_is_inside_plane(x, x_min, epsilon) &&
            node_is_inside_circle_yz(y, z, yc, zc, R, epsilon);

        bool is_wall =
            node_is_inside_range(x, x_min, x_max, epsilon) &&
            node_is_on_cylinder_wall(y, z, yc, zc, R, epsilon);

        /*
            inflow boundary condition
            u_x = Poiseuille profile
            u_y = 0
            u_z = 0
        */
        if(is_inflow) {
            double vx = calc_velocity_pipe_poiseuille(
                            y, z, yc, zc, R, umax);

            fprintf(fp_dbc, "%d 0 %lf\n", i, vx);
            fprintf(fp_dbc, "%d 1 0.0\n", i);
            fprintf(fp_dbc, "%d 2 0.0\n", i);
        }

        /*
            no-slip wall boundary condition
            u_x = 0
            u_y = 0
            u_z = 0
        */
        else if(is_wall) {
            fprintf(fp_dbc, "%d 0 0.0\n", i);
            fprintf(fp_dbc, "%d 1 0.0\n", i);
            fprintf(fp_dbc, "%d 2 0.0\n", i);
        }
    }

    fclose(fp_dbc);
}


int main(
    int   argc,
    char* argv[])
{
    printf("\n");

    BBFE_DATA fe;
    SETTINGS  set;
    CONDITION cond;

    BB_vtk_void();

    args_manager(&set, argc, argv);
    set_condition(argv, &cond);

    BBFE_sys_read_node(&fe, set.infile_node, set.directory);

    printf("%s Start writing D_bc file.\n", CODENAME);

    int num_bc_nodes = detect_Dirichlet_bc(&fe, &cond);

    output_Dirichlet_bc_3d_pipe_poiseuille(
        &fe,
        &set,
        &cond,
        num_bc_nodes);

    printf("%s Finish writing Dirichlet b.c. file.\n", CODENAME);

    return 0;
}