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

static const char* CODENAME = "inlet_bfs_only >";
static const char* VOIDNAME = "                >";
static const int BUFFER_SIZE = 10000;

static const char* DEF_DIRECTORY      = ".";
static const char* DEF_FILENAME_NODE  = "node.dat";
static const char* DEF_FILENAME_D_BC  = "D_bc_inlet.dat";

static const char* OPTION_DIRECTORY    = "-d";
static const char* OPTION_INFILE_NODE  = "-in";
static const char* OPTION_OUTFILE_D_BC = "-od";

/*
    false:
      流入面の上下端ノードを出力しない。
      壁 no-slip 条件を別ファイルで与える場合に推奨。

    true:
      流入面の上下端ノードも出力する。
*/
static const bool INCLUDE_INLET_WALL_NODES = false;

typedef struct {
    const char* infile_node;
    const char* outfile_dbc;
    const char* directory;
} SETTINGS;

typedef struct {
    double U_max;
    double delta;
} CONDITION;

typedef struct {
    double x_min;
    double y_min;
    double y_max;
    double epsilon;
} INLET_BOUNDARY;


void args_manager(
        SETTINGS* set,
        int       argc,
        char*     argv[])
{
    if(argc < 3) {
        printf("%s ERROR: Please specify parameters for BFS inlet condition.\n", CODENAME);
        printf("%s Format:\n", VOIDNAME);
        printf("%s ./inlet_bfs_only [1: U_max] [2: delta]\n\n", VOIDNAME);

        printf("%s Example:\n", VOIDNAME);
        printf("%s ./inlet_bfs_only 1.0 0.3 -d . -in node.dat -od D_bc_inlet.dat\n\n", VOIDNAME);

        printf("%s -d  [directory]      Default: ./\n", VOIDNAME);
        printf("%s -in [input nodes]    Default: node.dat\n", VOIDNAME);
        printf("%s -od [output D_bc]    Default: D_bc_inlet.dat\n", VOIDNAME);

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
    cond->U_max = atof(argv[1]);
    cond->delta = atof(argv[2]);

    printf("%s Maximum inlet velocity : %e\n", CODENAME, cond->U_max);
    printf("%s Wall ramp thickness    : %e\n", CODENAME, cond->delta);

    if(INCLUDE_INLET_WALL_NODES) {
        printf("%s Inlet wall-edge nodes : included\n", CODENAME);
    }
    else {
        printf("%s Inlet wall-edge nodes : excluded\n", CODENAME);
    }
}


static double max3(
        const double a,
        const double b,
        const double c)
{
    double m = a;

    if(b > m) {
        m = b;
    }

    if(c > m) {
        m = c;
    }

    return m;
}


/*
    node.dat から、最も左側にある流入面を自動検出する。

    判定:
      x が全ノード中の最小値 x_min に近いノードを流入面ノードとする。

    さらに、その流入面上の y_min, y_max も自動検出する。
*/
static void detect_leftmost_inlet_boundary(
        BBFE_DATA*      fe,
        INLET_BOUNDARY* inlet)
{
    double x_min =  1.0e+100;
    double x_max = -1.0e+100;
    double y_min =  1.0e+100;
    double y_max = -1.0e+100;
    double z_min =  1.0e+100;
    double z_max = -1.0e+100;

    for(int i = 0; i < fe->total_num_nodes; i++) {
        const double x = fe->x[i][0];
        const double y = fe->x[i][1];
        const double z = fe->x[i][2];

        if(x < x_min) {
            x_min = x;
        }

        if(x > x_max) {
            x_max = x;
        }

        if(y < y_min) {
            y_min = y;
        }

        if(y > y_max) {
            y_max = y;
        }

        if(z < z_min) {
            z_min = z;
        }

        if(z > z_max) {
            z_max = z;
        }
    }

    const double lx = x_max - x_min;
    const double ly = y_max - y_min;
    const double lz = z_max - z_min;

    const double length_scale = max3(lx, ly, lz);

    /*
        座標の丸め誤差や微小なずれを吸収するための許容幅。
        大きすぎると隣の x 面まで拾うので注意。
    */
    double epsilon = 1.0e-6 * length_scale;

    if(epsilon < 1.0e-10) {
        epsilon = 1.0e-10;
    }

    /*
        最左端面上の y 範囲を検出する。
    */
    double inlet_y_min =  1.0e+100;
    double inlet_y_max = -1.0e+100;

    int count = 0;

    for(int i = 0; i < fe->total_num_nodes; i++) {
        const double x = fe->x[i][0];
        const double y = fe->x[i][1];

        if(x <= x_min + epsilon) {
            if(y < inlet_y_min) {
                inlet_y_min = y;
            }

            if(y > inlet_y_max) {
                inlet_y_max = y;
            }

            count++;
        }
    }

    if(count == 0) {
        printf("%s ERROR: No nodes found on leftmost inlet boundary.\n", CODENAME);
        exit(1);
    }

    inlet->x_min   = x_min;
    inlet->y_min   = inlet_y_min;
    inlet->y_max   = inlet_y_max;
    inlet->epsilon = epsilon;

    printf("%s Detected inlet boundary by leftmost x.\n", CODENAME);
    printf("%s inlet x_min : %.15e\n", CODENAME, inlet->x_min);
    printf("%s inlet y_min : %.15e\n", CODENAME, inlet->y_min);
    printf("%s inlet y_max : %.15e\n", CODENAME, inlet->y_max);
    printf("%s inlet eps   : %.15e\n", CODENAME, inlet->epsilon);
    printf("%s inlet nodes detected on x-min side, including edge nodes: %d\n",
            CODENAME, count);
}


/*
    最小 x 側の面を流入面として判定する。

    INCLUDE_INLET_WALL_NODES = false のため、
    y_min, y_max の上下端ノードは除外する。
*/
static bool node_is_leftmost_inlet(
        const double x,
        const double y,
        const INLET_BOUNDARY* inlet)
{
    const double eps = inlet->epsilon;

    if(!(x <= inlet->x_min + eps)) {
        return false;
    }

    if(INCLUDE_INLET_WALL_NODES) {
        return (inlet->y_min - eps <= y && y <= inlet->y_max + eps);
    }
    else {
        return (inlet->y_min + eps < y && y < inlet->y_max - eps);
    }
}


/*
    C2 continuous smoothstep:
      s = 0 で 0
      s = 1 で 1
      s = 0, 1 で勾配も滑らか
*/
static double smoothstep5(
        const double s_in)
{
    double s = s_in;

    if(s <= 0.0) {
        return 0.0;
    }

    if(s >= 1.0) {
        return 1.0;
    }

    return s*s*s*(10.0 - 15.0*s + 6.0*s*s);
}


/*
    流入速度:
      最大流速は U_max。
      壁近傍 delta の範囲だけ smoothstep で 0 -> U_max に立ち上げる。
      壁近傍以外では常に U_max。

    重要:
      ここでは断面平均速度の補正は行わない。
      したがって U_max = 1.0 を指定した場合、最大値は必ず 1.0。
*/
static double calc_bfs_inlet_velocity(
        const double y,
        const double y_min_inlet,
        const double y_max_inlet,
        const double U_max,
        const double delta)
{
    const double H = y_max_inlet - y_min_inlet;

    if(H <= 0.0) {
        return 0.0;
    }

    /*
        delta <= 0 の場合:
          上下端以外は完全一様流速 U_max。
    */
    if(delta <= 0.0) {
        if(y <= y_min_inlet || y >= y_max_inlet) {
            return 0.0;
        }

        return U_max;
    }

    double delta_eff = delta;

    /*
        上下の ramp が重ならないように制限する。
    */
    if(delta_eff >= 0.5 * H) {
        delta_eff = 0.499 * H;
    }

    const double d_lower = y - y_min_inlet;
    const double d_upper = y_max_inlet - y;

    double d_wall = d_lower;

    if(d_upper < d_wall) {
        d_wall = d_upper;
    }

    /*
        壁上では 0。
        ただし INCLUDE_INLET_WALL_NODES = false なら、
        通常この点は出力対象にはならない。
    */
    if(d_wall <= 0.0) {
        return 0.0;
    }

    /*
        壁から delta 以上離れていれば、常に最大流速。
        U_max = 1.0 なら、ここは常に 1.0。
    */
    if(d_wall >= delta_eff) {
        return U_max;
    }

    /*
        壁近傍のみ、0 -> U_max に滑らかに立ち上げる。
    */
    const double s = d_wall / delta_eff;

    return U_max * smoothstep5(s);
}


/*
    D_bc 対象の入口節点数を数える。

    注意:
      この関数が返すのは「入口節点数」。
      出力では1節点につき ux, uy, uz の3行を書くため、
      D_bc のヘッダ先頭には

        num_bc_nodes * 3

      を出力する。
*/
int detect_inlet_Dirichlet_bc(
        BBFE_DATA*             fe,
        const INLET_BOUNDARY*  inlet)
{
    int num_bc_nodes = 0;

    for(int i = 0; i < fe->total_num_nodes; i++) {
        const double x = fe->x[i][0];
        const double y = fe->x[i][1];

        if(node_is_leftmost_inlet(x, y, inlet)) {
            num_bc_nodes += 1;
        }
    }

    return num_bc_nodes;
}


void output_inlet_Dirichlet_bc_bfs(
        BBFE_DATA*             fe,
        SETTINGS*              set,
        CONDITION*             cond,
        const INLET_BOUNDARY*  inlet,
        const int              num_bc_nodes)
{
    char fname_dbc[BUFFER_SIZE];

    snprintf(fname_dbc, BUFFER_SIZE, "%s/%s", set->directory, set->outfile_dbc);

    FILE* fp_dbc = fopen(fname_dbc, "w");

    if(fp_dbc == NULL) {
        printf("%s ERROR: File \"%s\" cannot be opened.\n",
                CODENAME, fname_dbc);
        exit(1);
    }

    /*
        D_bc 形式:
          1行目: 条件行数 3

        1節点につき
          node_id 0 ux
          node_id 1 uy
          node_id 2 uz

        の3行を出力するので、
          条件行数 = num_bc_nodes * 3
        とする。

        例:
          num_bc_nodes = 254 の場合
          header = 762 3
    */
    const int block_length = 3;
    const int num_bc_rows  = num_bc_nodes * block_length;

    fprintf(fp_dbc, "%d %d\n", num_bc_rows, block_length);

    for(int i = 0; i < fe->total_num_nodes; i++) {
        const double x = fe->x[i][0];
        const double y = fe->x[i][1];

        if(node_is_leftmost_inlet(x, y, inlet)) {
            const double vx = calc_bfs_inlet_velocity(
                                    y,
                                    inlet->y_min,
                                    inlet->y_max,
                                    cond->U_max,
                                    cond->delta);

            fprintf(fp_dbc, "%d 0 %.15e\n", i, vx);
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
    INLET_BOUNDARY inlet;

    BB_vtk_void();

    args_manager(&set, argc, argv);
    set_condition(argv, &cond);

    BBFE_sys_read_node(&fe, set.infile_node, set.directory);

    /*
        node.dat から、一番左端の入口面を自動検出する。
    */
    detect_leftmost_inlet_boundary(&fe, &inlet);

    printf("%s Start detecting inlet Dirichlet boundary conditions.\n", CODENAME);

    const int num_bc_nodes = detect_inlet_Dirichlet_bc(&fe, &inlet);
    const int num_bc_rows  = num_bc_nodes * 3;

    printf("%s Number of inlet Dirichlet nodes: %d\n", CODENAME, num_bc_nodes);
    printf("%s Number of inlet Dirichlet rows : %d\n", CODENAME, num_bc_rows);
    printf("%s D_bc header will be: %d 3\n", CODENAME, num_bc_rows);

    printf("%s Start writing inlet D_bc file.\n", CODENAME);

    output_inlet_Dirichlet_bc_bfs(&fe, &set, &cond, &inlet, num_bc_nodes);

    printf("%s Finish writing inlet Dirichlet b.c. file.\n", CODENAME);

    return 0;
}
