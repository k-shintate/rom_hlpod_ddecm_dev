#include <stdio.h>
#include "BB/std.h"
#include "BB/calc.h"
#include "BB/vtk.h"
#include "BBFE/sys/FE_dataset.h"
#include "BBFE/sys/read.h"
#include "BBFE/sys/memory.h"
#include "BBFE/sys/write.h"

typedef struct
{
    const char* directory;
    const char* infile_elem;
    const char* outfile_graph;
} SETTINGS;

static const char* CODENAME = "util/converter/elem2graph >";

static const char* OPTION_DIRECTORY = "-d";
static const char* OPTION_ELEM      = "-ie";
static const char* OPTION_GRAPH     = "-og";

static const char* DEF_DIRECTORY     = ".";
static const char* DEF_INFILE_ELEM   = "elem.dat";
static const char* DEF_OUTFILE_GRAPH = "graph.dat";

void args_manager_elem2graph(
    SETTINGS*   set,
    int         argc,
    char*       argv[])
{
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
                argc, argv, OPTION_ELEM);
    if(num == -1) {
        set->infile_elem = DEF_INFILE_ELEM;
    }
    else {
        set->infile_elem = argv[num+1];
    }

    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_GRAPH);
    if(num == -1) {
        set->outfile_graph = DEF_OUTFILE_GRAPH;
    }
    else {
        set->outfile_graph = argv[num+1];
    }

    printf("%s Input/output directory  : %s\n", CODENAME, set->directory);
    printf("%s Input filename (elem)   : %s\n", CODENAME, set->infile_elem);
    printf("%s Output filename (graph) : %s\n", CODENAME, set->outfile_graph);
}


void BBFE_conv_elem2graph_hex(
    BBFE_DATA*   fe,
    const char*  filename,
    const char*  directory)
{

    FILE* fp;
    fp = BBFE_sys_write_fopen(fp, filename, directory);

    fprintf(fp, "%d\n", fe->total_num_elems);

    for(int e=0; e<(fe->total_num_elems); e++) {
        fprintf(fp, "%d %d ", e, fe->local_num_nodes);

        for(int i=0; i<(fe->local_num_nodes); i++) {
            fprintf(fp, "%d ", fe->conn[e][i]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}


int main(
    int     argc,
    char*   argv[])
{
    printf("\n");
    BB_calc_void();
    BB_vtk_void();

    BBFE_DATA   fe;
    SETTINGS    set;

    args_manager_elem2graph(&set, argc, argv);
    BBFE_sys_read_elem(&fe, set.infile_elem, set.directory, 1);
    BBFE_conv_elem2graph_hex(&fe, set.outfile_graph, set.directory);

    return 0;
}