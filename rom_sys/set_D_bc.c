
#include "set_D_bc.h"

/*要変更*/
void ROM_sys_hlpod_fe_set_bc_id(
        BBFE_BC*   	bc,
        const int   total_num_nodes,
        const int   num_dofs_on_node,
        ROM_BC*	    rom_bc)
{
    int j = 0;
    rom_bc->D_bc_node_id = BB_std_calloc_1d_int(rom_bc->D_bc_node_id, bc->num_D_bcs);

    for(int i = 0; i < total_num_nodes; i++) {
        if( bc->D_bc_exists[num_dofs_on_node*i] ) {
            rom_bc->D_bc_node_id[j] = i;
            j++;
        }
    }
    
}

void ROM_sys_hlpod_fe_manusol_set_bc_id(
        BBFE_DATA* 	fe,
        BBFE_BC*   	bc,
        const int   num_dofs_on_node,
        double     	t,
        double      (*func)(double, double, double, double), // scalar function(x, y, z, t)
        ROM_BC*	    rom_bc)
{
    int index = 0;

    for(int i=0; i<bc->num_D_bcs; i++) {
        index = rom_bc->D_bc_node_id[i];
        bc->imposed_D_val[index] = func(fe->x[index][0], fe->x[index][1], fe->x[index][2], t);
    }
}



void ROM_sys_hlpod_fe_monowrap_set_Dirichlet_bc(
        MONOLIS*    monolis,
        const int   num_dofs_on_node,
        BBFE_BC*    bc,
        double*     g_rhs,
        ROM_BC*	    rom_bc)
{
    int index =0;
    for(int i=0; i<(bc->num_D_bcs/num_dofs_on_node); i++) {
        for(int k=0; k<num_dofs_on_node; k++) {
            index = rom_bc->D_bc_node_id[i];
            monolis_set_Dirichlet_bc_R(
                    monolis,
                    g_rhs,
                    index,
                    k,
                    bc->imposed_D_val[ num_dofs_on_node*index + k]);
        }
    }
}




void HROM_manusol_set_bc(
		BBFE_DATA* 	fe,
		BBFE_BC*   	bc,
        const int   num_dofs_on_node,
		double     	t,
		double      (*func)(double, double, double, double), // scalar function(x, y, z, t)
		HLPOD_MAT*	hlpod_mat)
{
    int index = 0;

    for(int i=0; i < hlpod_mat->num_hr_D_bc_nodes; i++) {
        index = hlpod_mat->hr_D_bc_node_id[i];
		bc->imposed_D_val[index] = func(fe->x[index][0], fe->x[index][1], fe->x[index][2], t);
    }

}
