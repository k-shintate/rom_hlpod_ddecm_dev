
#include "write_BB.h"
#include "write_std.h"
#include "read_BB.h"

static const char* CODENAME = "HLPOD_std/write >";

static const char* OUTPUT_FILENAME_SVD          = "singularvalue.dat";
static const char* OUTPUT_FILENAME_NNLS_RESIDUALS    = "NNLS_residuals.dat";
static const char* OUTPUT_FILENAME_NNLS_NUM_ELEMS    = "NNLS_num_elems.dat";

static const int BUFFER_SIZE = 10000;


void ROM_std_hlpod_write_num_modes(
    const int       num_modes,
    const int		subdomain_id,
    const char*		label,
    const char*		directory)
{
    const char* filename;
    char fname[BUFFER_SIZE];
    FILE* fp;

    snprintf(fname, BUFFER_SIZE,"%s/num_modes.dat.%d", label, subdomain_id);

    fp = ROM_BB_write_fopen(fp, fname, directory);
    fprintf(fp,"%d\n", num_modes);
    fclose(fp);
}

void ROM_std_hlpod_write_pod_modes(
    double**		basis,
    const int 		num_modes,
    const int 		n_internal_vertex,
    const int		subdomain_id,
    const int 		dof,
    const char*		label,
    const char*		directory)
{
    char fname[BUFFER_SIZE];

    for(int k = 0; k < num_modes; k++){
        snprintf(fname, BUFFER_SIZE,"%s/subdomain%d/lb_pod_modes%d.dat", label, subdomain_id, k);

        FILE* fp;
        fp = ROM_BB_write_fopen(fp, fname, directory);

        fprintf(fp,"#%s\n", label);
        fprintf(fp,"%d %d\n", n_internal_vertex, dof);
    
        for(int i = 0; i < n_internal_vertex; i++){
            for(int j = 0; j < dof; j++){
                fprintf(fp,"%.30e ", basis[i*dof + j][k]);
            }
            fprintf(fp,"\n");
        }

        fclose(fp);
    }
}


void ROM_std_hlpod_write_time_svd(
    double			time_svd,
    const int 		subdomain_id,
    const char*		label,
    const char*		directory)
{
    FILE* fp;
    char fname[BUFFER_SIZE];

    snprintf(fname, BUFFER_SIZE, "%s/subdomain%d/time_SVD.txt", label, subdomain_id);

    fp = ROM_BB_write_fopen(fp, fname, directory);
    fprintf(fp,"%e\n",time_svd);
    fclose(fp);
}


void ROM_std_hlpod_write_singular_values(
    double*			singular_value,
    const int 		num_modes,
    const int 		subdomain_id,
    const char*		label,
    const char*		directory)
{
    FILE* fp;
    char fname[BUFFER_SIZE];

    snprintf(fname, BUFFER_SIZE, "%s/subdomain%d/%s", label, subdomain_id, OUTPUT_FILENAME_SVD);

    fp = ROM_BB_write_fopen(fp, fname, directory);
    for(int i = 0; i < num_modes; i++){
        fprintf(fp,"%.16f\n", singular_value[i]);
    }
    fclose(fp);
}

void hr_write_NNLS_residual(
    const double 	residual,
    const int		myrank,
    const int		num,
    const char*		directory)
{
    const char* filename;
    char fname[BUFFER_SIZE];
    snprintf(fname, BUFFER_SIZE,"hr_prm/%s.%d.%d", OUTPUT_FILENAME_NNLS_RESIDUALS, myrank, num);

    FILE* fp;

    fp = ROM_BB_write_fopen(fp, fname, directory);
    fprintf(fp,"%.30e\n", residual);
    fclose(fp);
}

void hr_write_NNLS_num_elems(
    const int 		num_elems,
    const int		myrank,
    const int		num,
    const char*		directory)
{
    const char* filename;
    char fname[BUFFER_SIZE];
    snprintf(fname, BUFFER_SIZE,"hr_prm/%s.%d.%d", OUTPUT_FILENAME_NNLS_NUM_ELEMS, myrank, num);

    FILE* fp;

    fp = ROM_BB_write_fopen(fp, fname, directory);
    fprintf(fp,"%d\n",num_elems);
    fclose(fp);
}

void ddhr_lb_write_reduced_mat(
    double**		reduced_mat,
    const int		m,
    const int 		n,
    const char*		directory)
{
    const char* filename;
    char fname[BUFFER_SIZE];
    snprintf(fname, BUFFER_SIZE, "DDECM/%s", "redued_mat.dat");
    filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname);
    FILE* fp;

    fp = ROM_BB_write_fopen(fp, filename, directory);
    for(int i = 0; i < m; i++){
        for (int j = 0; j < n; j++)
        {
            fprintf(fp,"%.16f\n",reduced_mat[i][j]);
        }		
    }
    fclose(fp);
}

void ROM_std_hlpod_output_calc_time_fopen(
    const char*  	directory)
{
    if(monolis_mpi_get_global_my_rank() == 0){
        FILE* fp;
        fp = ROM_BB_write_fopen(fp, "calctime/time_fem_total.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "calctime/time_lpod_total.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "calctime/time_hr_total.txt", directory);
        fclose(fp);
    }
}

void ROM_std_hlpod_output_calc_time(
    double			calctime,
    double 			t,
    const char*  	fname,
    const char*  	directory)
{ 
    if(monolis_mpi_get_global_my_rank() == 0){
        FILE* fp;
        fp = ROM_BB_write_fopen(fp, fname, directory);
        fprintf(fp, "%e %e\n", t, calctime);
        fclose(fp);
    }
}


void ROM_std_hlpod_output_add_calc_time(
    double			calctime,
    double 			t,
    const char*  	fname,
    const char*  	directory)
{
    if(monolis_mpi_get_global_my_rank() == 0){
        FILE* fp;
        fp = ROM_BB_write_add_fopen(fp, fname, directory);
        fprintf(fp, "%e %e\n", t, calctime);
        fclose(fp);
    }
}


void ROM_std_hlpod_output_fem_solver_prm_fopen(
    const char*  	directory)
{
    if(monolis_mpi_get_global_my_rank() == 0){
        FILE* fp;
        fp = ROM_BB_write_fopen(fp, "fem_solver_prm/converge_iter.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "fem_solver_prm/converge_reidual.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "fem_solver_prm/time_solver.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "fem_solver_prm/time_preparing.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "fem_solver_prm/time_spmv.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "fem_solver_prm/time_inner_product.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "fem_solver_prm/time_precondition.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "fem_solver_prm/time_comm_inner_product.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "fem_solver_prm/time_comm_spmv.txt", directory);
        fclose(fp);
    }
}

void ROM_std_hlpod_write_solver_prm(
    MONOLIS*  		monolis,
    double 			t,
    const char*     fname,
    const char*  	directory)
{
    char filename[BUFFER_SIZE];

    //if(monolis_mpi_get_global_my_rank() == 0){
        FILE* fp;
        int int_val;
        double val;
    
        snprintf(filename, BUFFER_SIZE, "%s/converge_iter.txt", fname);
        fp = ROM_BB_write_add_fopen(fp, filename, directory);
        monolis_get_converge_iter(monolis, &int_val);
        fprintf(fp, "%e %d\n", t, int_val);
        fclose(fp);

        snprintf(filename, BUFFER_SIZE,  "%s/converge_reidual.txt", fname);
        fp = ROM_BB_write_add_fopen(fp, filename, directory);
        monolis_get_converge_residual(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        snprintf(filename, BUFFER_SIZE,  "%s/time_solver.txt", fname);
        fp = ROM_BB_write_add_fopen(fp, filename, directory);
        monolis_get_time_solver(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        snprintf(filename, BUFFER_SIZE,  "%s/time_preparing.txt", fname);
        fp = ROM_BB_write_add_fopen(fp, filename, directory);
        monolis_get_time_preparing(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        snprintf(filename, BUFFER_SIZE,  "%s/time_spmv.txt", fname);
        fp = ROM_BB_write_add_fopen(fp, filename, directory);
        monolis_get_time_spmv(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        snprintf(filename, BUFFER_SIZE,  "%s/time_inner_product.txt", fname);
        fp = ROM_BB_write_add_fopen(fp, filename, directory);
        monolis_get_time_inner_product(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);


        snprintf(filename, BUFFER_SIZE,  "%s/time_precondition.txt", fname);
        fp = ROM_BB_write_add_fopen(fp, filename, directory);
        monolis_get_time_precondition(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        snprintf(filename, BUFFER_SIZE,  "%s/time_comm_inner_product.txt", fname);
        fp = ROM_BB_write_add_fopen(fp, filename, directory);
        monolis_get_time_comm_inner_product(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        snprintf(filename, BUFFER_SIZE,  "%s/time_comm_spmv.txt", fname);
        fp = ROM_BB_write_add_fopen(fp, filename, directory);
        monolis_get_time_comm_spmv(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);
    //}
}

void ROM_std_hlpod_write_solver_prm_fopen(
    const char*     output_directory,
    const char*  	directory)
{
    char filename[BUFFER_SIZE];
    FILE* fp;


        snprintf(filename, BUFFER_SIZE, "%s/converge_iter.txt", output_directory);
        fp = ROM_BB_write_fopen(fp, filename, directory);
        fclose(fp);

        snprintf(filename, BUFFER_SIZE,  "%s/time_spmv.txt", output_directory);
        fp = ROM_BB_write_fopen(fp, filename, directory);
        fclose(fp);

        snprintf(filename, BUFFER_SIZE,  "%s/time_solver.txt", output_directory);
        fp = ROM_BB_write_fopen(fp, filename, directory);
        fclose(fp);

        snprintf(filename, BUFFER_SIZE,  "%s/time_preparing.txt", output_directory);
        fp = ROM_BB_write_fopen(fp, filename, directory);
        fclose(fp);

        snprintf(filename, BUFFER_SIZE,  "%s/time_spmv.txt", output_directory);
        fp = ROM_BB_write_fopen(fp, filename, directory);
        fclose(fp);

        snprintf(filename, BUFFER_SIZE,  "%s/time_inner_product.txt", output_directory);
        fp = ROM_BB_write_fopen(fp, filename, directory);
        fclose(fp);

        snprintf(filename, BUFFER_SIZE,  "%s/time_precondition.txt", output_directory);
        fp = ROM_BB_write_fopen(fp, filename, directory);
        fclose(fp);

        snprintf(filename, BUFFER_SIZE,  "%s/time_comm_inner_product.txt", output_directory);
        fp = ROM_BB_write_fopen(fp, filename, directory);
        fclose(fp);

        snprintf(filename, BUFFER_SIZE,  "%s/time_comm_spmv.txt", output_directory);
        fp = ROM_BB_write_fopen(fp, filename, directory);
        fclose(fp);
}

void ROM_std_hlpod_fem_solver_prm(
    MONOLIS*  		monolis,
    const char*  	directory,
    double 			t)
{
    if(monolis_mpi_get_global_my_rank() == 0){
        FILE* fp;
        int int_val;
        double val;

        fp = ROM_BB_write_add_fopen(fp, "fem_solver_prm/converge_iter.txt", directory);
        monolis_get_converge_iter(monolis, &int_val);
        fprintf(fp, "%e %d\n", t, int_val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "fem_solver_prm/converge_reidual.txt", directory);
        monolis_get_converge_residual(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "fem_solver_prm/time_solver.txt", directory);
        monolis_get_time_solver(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "fem_solver_prm/time_preparing.txt", directory);
        monolis_get_time_preparing(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "fem_solver_prm/time_spmv.txt", directory);
        monolis_get_time_spmv(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "fem_solver_prm/time_inner_product.txt", directory);
        monolis_get_time_inner_product(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);


        fp = ROM_BB_write_add_fopen(fp, "fem_solver_prm/time_precondition.txt", directory);
        monolis_get_time_precondition(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "fem_solver_prm/time_comm_inner_product.txt", directory);
        monolis_get_time_comm_inner_product(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "fem_solver_prm/time_comm_spmv.txt", directory);
        monolis_get_time_comm_spmv(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);
    }
}

void ROM_std_hlpod_output_rom_solver_prm_fopen(
    const char*  	directory)
{
    if(monolis_mpi_get_global_my_rank() == 0){
        FILE* fp;
        fp = ROM_BB_write_fopen(fp, "pod_solver_prm/converge_iter.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "pod_solver_prm/converge_reidual.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "pod_solver_prm/time_solver.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "pod_solver_prm/time_preparing.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "pod_solver_prm/time_spmv.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "pod_solver_prm/time_inner_product.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "pod_solver_prm/time_precondition.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "pod_solver_prm/time_comm_inner_product.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "pod_solver_prm/time_comm_spmv.txt", directory);
        fclose(fp);
    }
}

void ROM_std_hlpod_rom_solver_prm(
    MONOLIS*  		monolis,
    const char*  	directory,
    double 			t)
{
    if(monolis_mpi_get_global_my_rank() == 0){
        FILE* fp;
        int int_val;
        double val;

        fp = ROM_BB_write_add_fopen(fp, "pod_solver_prm/converge_iter.txt", directory);
        monolis_get_converge_iter(monolis, &int_val);
        fprintf(fp, "%e %d\n", t, int_val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "pod_solver_prm/converge_reidual.txt", directory);
        monolis_get_converge_residual(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "pod_solver_prm/time_solver.txt", directory);
        monolis_get_time_solver(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "pod_solver_prm/time_preparing.txt", directory);
        monolis_get_time_preparing(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "pod_solver_prm/time_spmv.txt", directory);
        monolis_get_time_spmv(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "pod_solver_prm/time_inner_product.txt", directory);
        monolis_get_time_inner_product(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);


        fp = ROM_BB_write_add_fopen(fp, "pod_solver_prm/time_precondition.txt", directory);
        monolis_get_time_precondition(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "pod_solver_prm/time_comm_inner_product.txt", directory);
        monolis_get_time_comm_inner_product(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "pod_solver_prm/time_comm_spmv.txt", directory);
        monolis_get_time_comm_spmv(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);
    }
}

void ROM_std_hlpod_hrom_solver_prm_fopen(
    const char*  	directory)
{
    if(monolis_mpi_get_global_my_rank() == 0){
        FILE* fp;
        fp = ROM_BB_write_fopen(fp, "hr_solver_prm/converge_iter.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "hr_solver_prm/converge_reidual.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "hr_solver_prm/time_solver.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "hr_solver_prm/time_preparing.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "hr_solver_prm/time_spmv.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "hr_solver_prm/time_inner_product.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "hr_solver_prm/time_precondition.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "hr_solver_prm/time_comm_inner_product.txt", directory);
        fclose(fp);

        fp = ROM_BB_write_fopen(fp, "hr_solver_prm/time_comm_spmv.txt", directory);
        fclose(fp);
    }
}

void ROM_std_hlpod_hrom_solver_prm(
    MONOLIS*  		monolis,
    const char*  	directory,
    double 			t)
{
    if(monolis_mpi_get_global_my_rank() == 0){
        FILE* fp;
        int int_val;
        double val;

        fp = ROM_BB_write_add_fopen(fp, "hr_solver_prm/converge_iter.txt", directory);
        monolis_get_converge_iter(monolis, &int_val);
        fprintf(fp, "%e %d\n", t, int_val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "hr_solver_prm/converge_reidual.txt", directory);
        monolis_get_converge_residual(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "hr_solver_prm/time_solver.txt", directory);
        monolis_get_time_solver(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "hr_solver_prm/time_preparing.txt", directory);
        monolis_get_time_preparing(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "hr_solver_prm/time_spmv.txt", directory);
        monolis_get_time_spmv(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "hr_solver_prm/time_inner_product.txt", directory);
        monolis_get_time_inner_product(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "hr_solver_prm/time_precondition.txt", directory);
        monolis_get_time_precondition(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "hr_solver_prm/time_comm_inner_product.txt", directory);
        monolis_get_time_comm_inner_product(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);

        fp = ROM_BB_write_add_fopen(fp, "hr_solver_prm/time_comm_spmv.txt", directory);
        monolis_get_time_comm_spmv(monolis, &val);
        fprintf(fp, "%e %e\n", t, val);
        fclose(fp);
    }
}

