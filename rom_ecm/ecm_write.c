
#include "ecm_write.h"

void output_hr_monolis_solver_prm(
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

