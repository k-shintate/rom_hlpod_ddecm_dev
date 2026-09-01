#include "core_FOM_ST_spaceblock.h"
#include "st_fom_project_adapter.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void st_accuracy_path(
    const ST_FOM_PROJECT* project,
    char* path,
    size_t path_size)
{
    const char* env = getenv("ST_FOM_ACCURACY_CSV");

    if(env != NULL && env[0] != '\0'){
        snprintf(path, path_size, "%s", env);
        return;
    }

    if(project->directory[0] == '\0' ||
       strcmp(project->directory, ".") == 0 ||
       strcmp(project->directory, "./") == 0)
    {
        snprintf(path, path_size, "st_accuracy.csv");
    }
    else{
        snprintf(path, path_size, "%s/st_accuracy.csv", project->directory);
    }
}

int main(int argc, char* argv[])
{
    ST_FOM_PROJECT project;
    ST_DIFFUSION_SPACEBLOCK solver;
    char accuracy_csv[4096];
    int num_windows;
    int status;
    double time_start;
    double time_end;

    MPI_Init(&argc, &argv);
    time_start = MPI_Wtime();

    memset(&project, 0, sizeof(project));
    memset(&solver, 0, sizeof(solver));

    status = ST_fom_project_initialize(&project, argc, argv);
    if(status != 0){
        if(monolis_mpi_get_global_my_rank() == 0){
            fprintf(stderr, "ERROR: ST-FOM project initialization failed.\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }

    num_windows = ST_fom_project_num_windows(&project);
    if(num_windows <= 0){
        ST_fom_project_finalize(&project);
        MPI_Abort(MPI_COMM_WORLD, 2);
        return 2;
    }

    st_accuracy_path(&project, accuracy_csv, sizeof(accuracy_csv));

    status = ST_diffusion_spaceblock_initialize(
        &solver,
        &project.fom,
        project.time_degree,
        project.slabs_per_window,
        accuracy_csv);
    if(status != ST_FOM_SUCCESS){
        if(monolis_mpi_get_global_my_rank() == 0){
            fprintf(
                stderr,
                "ERROR: ST-FOM spaceblock initialization failed: %d\n",
                status);
        }
        ST_fom_project_finalize(&project);
        MPI_Abort(MPI_COMM_WORLD, 3);
        return 3;
    }

    status = ST_diffusion_spaceblock_run(&solver, num_windows);
    if(status != ST_FOM_SUCCESS){
        if(monolis_mpi_get_global_my_rank() == 0){
            fprintf(stderr, "ERROR: ST-FOM solve failed: %d\n", status);
        }
        ST_diffusion_spaceblock_finalize(&solver);
        ST_fom_project_finalize(&project);
        MPI_Abort(MPI_COMM_WORLD, 4);
        return 4;
    }

    ST_diffusion_spaceblock_finalize(&solver);
    ST_fom_project_finalize(&project);

    time_end = MPI_Wtime();
    if(monolis_mpi_get_global_my_rank() == 0){
        printf("\nNormal ST-FOM completed.\n");
        printf("Accuracy CSV: %s\n", accuracy_csv);
        printf("** Total time: %.15e\n", time_end - time_start);
    }

    MPI_Finalize();
    return 0;
}
