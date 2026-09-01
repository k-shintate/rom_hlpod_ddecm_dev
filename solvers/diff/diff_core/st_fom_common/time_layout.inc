/*
 * Time-element and space-time layout helpers for st_fom_common.c.
 *
 * The dG matrices below intentionally match the formulas used by the supplied
 * core_FOM_ST implementation for dG(0), dG(1), and dG(2).
 */

static void st_fom_time_element_clear(ST_FOM_TIME_ELEMENT* time)
{
    int i;
    int j;
    int q;

    memset(time, 0, sizeof(*time));
    time->family = ST_FOM_TIME_FAMILY_DG;

    for(i = 0; i < ST_FOM_MAX_TIME_DOF; i++){
        time->dof_tau[i] = 0.0;
        time->e_minus[i] = 0.0;
        time->e_plus[i] = 0.0;
        for(j = 0; j < ST_FOM_MAX_TIME_DOF; j++){
            time->Mt[i][j] = 0.0;
            time->Dt[i][j] = 0.0;
            time->E_self[i][j] = 0.0;
            time->E_prev[i][j] = 0.0;
        }
        for(q = 0; q < ST_FOM_MAX_TIME_IP; q++){
            time->psi[i][q] = 0.0;
        }
    }

    for(q = 0; q < ST_FOM_MAX_TIME_IP; q++){
        time->tau[q] = 0.0;
        time->weight[q] = 0.0;
    }
}

static void st_fom_fill_endpoint_flux(ST_FOM_TIME_ELEMENT* time)
{
    int a;
    int b;

    for(a = 0; a < time->n_dof; a++){
        for(b = 0; b < time->n_dof; b++){
            time->E_self[a][b] =
                time->e_minus[a] * time->e_minus[b];
            time->E_prev[a][b] =
                time->e_minus[a] * time->e_plus[b];
        }
    }
}

int ST_fom_time_element_init_dg(
    ST_FOM_TIME_ELEMENT* time,
    int degree)
{
    int q;

    if(time == NULL){
        return ST_FOM_ERR_ARGUMENT;
    }
    if(degree < 0 || degree > 2){
        return ST_FOM_ERR_UNSUPPORTED;
    }

    st_fom_time_element_clear(time);
    time->family = ST_FOM_TIME_FAMILY_DG;
    time->degree = degree;
    time->n_dof = degree + 1;
    time->is_nodal = 1;

    if(degree == 0){
        /* dG(0): one constant coefficient; use the right-endpoint label. */
        time->dof_tau[0] = 1.0;
        time->Mt[0][0] = 1.0;
        time->Dt[0][0] = 0.0;
        time->e_minus[0] = 1.0;
        time->e_plus[0] = 1.0;
        st_fom_fill_endpoint_flux(time);

        /* Right-endpoint quadrature: backward-Euler-compatible F(t_{n+1}). */
        time->num_time_ip = 1;
        time->tau[0] = 1.0;
        time->weight[0] = 1.0;
        time->psi[0][0] = 1.0;
        return ST_FOM_SUCCESS;
    }

    if(degree == 1){
        const double s = 0.57735026918962576451;

        time->dof_tau[0] = 0.0;
        time->dof_tau[1] = 1.0;

        time->Mt[0][0] = 1.0 / 3.0;
        time->Mt[0][1] = 1.0 / 6.0;
        time->Mt[1][0] = 1.0 / 6.0;
        time->Mt[1][1] = 1.0 / 3.0;

        time->Dt[0][0] = -0.5;
        time->Dt[0][1] =  0.5;
        time->Dt[1][0] = -0.5;
        time->Dt[1][1] =  0.5;

        time->e_minus[0] = 1.0;
        time->e_minus[1] = 0.0;
        time->e_plus[0] = 0.0;
        time->e_plus[1] = 1.0;
        st_fom_fill_endpoint_flux(time);

        time->num_time_ip = 2;
        time->tau[0] = 0.5 * (1.0 - s);
        time->tau[1] = 0.5 * (1.0 + s);
        time->weight[0] = 0.5;
        time->weight[1] = 0.5;

        for(q = 0; q < time->num_time_ip; q++){
            const double tau = time->tau[q];
            time->psi[0][q] = 1.0 - tau;
            time->psi[1][q] = tau;
        }
        return ST_FOM_SUCCESS;
    }

    /* degree == 2: nodal quadratic basis at tau = 0, 1/2, 1. */
    {
        const double s = 0.77459666924148337704;

        time->dof_tau[0] = 0.0;
        time->dof_tau[1] = 0.5;
        time->dof_tau[2] = 1.0;

        time->Mt[0][0] =  2.0 / 15.0;
        time->Mt[0][1] =  1.0 / 15.0;
        time->Mt[0][2] = -1.0 / 30.0;
        time->Mt[1][0] =  1.0 / 15.0;
        time->Mt[1][1] =  8.0 / 15.0;
        time->Mt[1][2] =  1.0 / 15.0;
        time->Mt[2][0] = -1.0 / 30.0;
        time->Mt[2][1] =  1.0 / 15.0;
        time->Mt[2][2] =  2.0 / 15.0;

        time->Dt[0][0] = -1.0 / 2.0;
        time->Dt[0][1] =  2.0 / 3.0;
        time->Dt[0][2] = -1.0 / 6.0;
        time->Dt[1][0] = -2.0 / 3.0;
        time->Dt[1][1] =  0.0;
        time->Dt[1][2] =  2.0 / 3.0;
        time->Dt[2][0] =  1.0 / 6.0;
        time->Dt[2][1] = -2.0 / 3.0;
        time->Dt[2][2] =  1.0 / 2.0;

        time->e_minus[0] = 1.0;
        time->e_minus[1] = 0.0;
        time->e_minus[2] = 0.0;
        time->e_plus[0] = 0.0;
        time->e_plus[1] = 0.0;
        time->e_plus[2] = 1.0;
        st_fom_fill_endpoint_flux(time);

        time->num_time_ip = 3;
        time->tau[0] = 0.5 * (1.0 - s);
        time->tau[1] = 0.5;
        time->tau[2] = 0.5 * (1.0 + s);
        time->weight[0] = 5.0 / 18.0;
        time->weight[1] = 4.0 / 9.0;
        time->weight[2] = 5.0 / 18.0;

        for(q = 0; q < time->num_time_ip; q++){
            const double tau = time->tau[q];
            time->psi[0][q] =  2.0 * tau * tau - 3.0 * tau + 1.0;
            time->psi[1][q] = -4.0 * tau * tau + 4.0 * tau;
            time->psi[2][q] =  2.0 * tau * tau - tau;
        }
    }

    return ST_FOM_SUCCESS;
}

int ST_fom_layout_initialize(
    ST_FOM_WINDOW_LAYOUT* layout,
    int num_space_points,
    int num_owned_space_points,
    int physical_dof,
    int num_slabs,
    int degree)
{
    ST_FOM_TIME_ELEMENT time;
    int status;

    if(layout == NULL || num_space_points <= 0 ||
       num_owned_space_points <= 0 ||
       num_owned_space_points > num_space_points ||
       physical_dof <= 0 || num_slabs <= 0)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    status = ST_fom_time_element_init_dg(&time, degree);
    if(status != ST_FOM_SUCCESS){
        return status;
    }

    layout->num_space_points = num_space_points;
    layout->num_owned_space_points = num_owned_space_points;
    layout->physical_dof = physical_dof;
    layout->num_slabs = num_slabs;
    layout->time = time;

    if(ST_fom_window_dof_per_space_point(layout) <= 0 ||
       ST_fom_window_vector_length(layout) == 0)
    {
        memset(layout, 0, sizeof(*layout));
        return ST_FOM_ERR_DIMENSION;
    }

    return ST_FOM_SUCCESS;
}

int ST_fom_resolve_owned_space_points(
    int num_local_space_points,
    int communicator_size,
    int communication_owned_space_points,
    int* num_owned_space_points)
{
    if(num_owned_space_points == NULL ||
       num_local_space_points <= 0 || communicator_size <= 0 ||
       communication_owned_space_points < 0)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    if(communicator_size == 1){
        *num_owned_space_points = num_local_space_points;
        return ST_FOM_SUCCESS;
    }

    if(communication_owned_space_points <= 0 ||
       communication_owned_space_points > num_local_space_points)
    {
        return ST_FOM_ERR_DIMENSION;
    }

    *num_owned_space_points = communication_owned_space_points;
    return ST_FOM_SUCCESS;
}

int ST_fom_window_dof_per_space_point(
    const ST_FOM_WINDOW_LAYOUT* layout)
{
    size_t value;

    if(layout == NULL || layout->physical_dof <= 0 ||
       layout->num_slabs <= 0 || layout->time.n_dof <= 0)
    {
        return -1;
    }

    value = (size_t)layout->physical_dof
        * (size_t)layout->num_slabs
        * (size_t)layout->time.n_dof;

    if(value > (size_t)INT_MAX){
        return -1;
    }

    return (int)value;
}

size_t ST_fom_window_vector_length(
    const ST_FOM_WINDOW_LAYOUT* layout)
{
    const int dof_per_space_point =
        ST_fom_window_dof_per_space_point(layout);

    if(layout == NULL || layout->num_space_points <= 0 ||
       dof_per_space_point <= 0)
    {
        return 0;
    }

    if((size_t)layout->num_space_points >
       ((size_t)-1) / (size_t)dof_per_space_point)
    {
        return 0;
    }

    return (size_t)layout->num_space_points
        * (size_t)dof_per_space_point;
}

int ST_fom_block_dof(
    const ST_FOM_WINDOW_LAYOUT* layout,
    int slab,
    int time_dof,
    int component)
{
    size_t value;

    if(layout == NULL || slab < 0 || slab >= layout->num_slabs ||
       time_dof < 0 || time_dof >= layout->time.n_dof ||
       component < 0 || component >= layout->physical_dof)
    {
        return -1;
    }

    value =
        (((size_t)slab * (size_t)layout->time.n_dof
          + (size_t)time_dof)
         * (size_t)layout->physical_dof)
        + (size_t)component;

    if(value > (size_t)INT_MAX){
        return -1;
    }

    return (int)value;
}

size_t ST_fom_vector_index(
    const ST_FOM_WINDOW_LAYOUT* layout,
    int space_point,
    int slab,
    int time_dof,
    int component)
{
    const int block_dof =
        ST_fom_block_dof(layout, slab, time_dof, component);
    const int dof_per_space_point =
        ST_fom_window_dof_per_space_point(layout);

    if(layout == NULL || space_point < 0 ||
       space_point >= layout->num_space_points ||
       block_dof < 0 || dof_per_space_point <= 0)
    {
        return (size_t)-1;
    }

    return (size_t)space_point * (size_t)dof_per_space_point
        + (size_t)block_dof;
}

int ST_fom_pack_dg0_slab(
    const ST_FOM_WINDOW_LAYOUT* layout,
    int slab,
    const double* state,
    double* window_state)
{
    int node;
    int component;

    if(layout == NULL || state == NULL || window_state == NULL ||
       layout->time.family != ST_FOM_TIME_FAMILY_DG ||
       layout->time.degree != 0 || layout->time.n_dof != 1 ||
       slab < 0 || slab >= layout->num_slabs)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    for(node = 0; node < layout->num_space_points; node++){
        for(component = 0; component < layout->physical_dof; component++){
            const size_t source =
                (size_t)node * (size_t)layout->physical_dof
                + (size_t)component;
            const size_t destination =
                ST_fom_vector_index(layout, node, slab, 0, component);
            window_state[destination] = state[source];
        }
    }

    return ST_FOM_SUCCESS;
}

int ST_fom_unpack_dg0_slab(
    const ST_FOM_WINDOW_LAYOUT* layout,
    int slab,
    const double* window_state,
    double* state)
{
    int node;
    int component;

    if(layout == NULL || state == NULL || window_state == NULL ||
       layout->time.family != ST_FOM_TIME_FAMILY_DG ||
       layout->time.degree != 0 || layout->time.n_dof != 1 ||
       slab < 0 || slab >= layout->num_slabs)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    for(node = 0; node < layout->num_space_points; node++){
        for(component = 0; component < layout->physical_dof; component++){
            const size_t source =
                ST_fom_vector_index(layout, node, slab, 0, component);
            const size_t destination =
                (size_t)node * (size_t)layout->physical_dof
                + (size_t)component;
            state[destination] = window_state[source];
        }
    }

    return ST_FOM_SUCCESS;
}

int ST_fom_extract_slab_right_trace(
    const ST_FOM_WINDOW_LAYOUT* layout,
    int slab,
    const double* window_state,
    double* trace)
{
    int node;
    int component;
    int alpha;

    if(layout == NULL || window_state == NULL || trace == NULL ||
       layout->time.family != ST_FOM_TIME_FAMILY_DG ||
       layout->time.n_dof <= 0 ||
       layout->time.n_dof > ST_FOM_MAX_TIME_DOF ||
       slab < 0 || slab >= layout->num_slabs)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    for(node = 0; node < layout->num_space_points; node++){
        for(component = 0; component < layout->physical_dof; component++){
            double value = 0.0;
            for(alpha = 0; alpha < layout->time.n_dof; alpha++){
                value += layout->time.e_plus[alpha]
                    * window_state[
                        ST_fom_vector_index(
                            layout, node, slab, alpha, component)];
            }
            trace[(size_t)node * (size_t)layout->physical_dof
                + (size_t)component] = value;
        }
    }

    return ST_FOM_SUCCESS;
}

int ST_fom_extract_last_right_trace(
    const ST_FOM_WINDOW_LAYOUT* layout,
    const double* window_state,
    double* trace)
{
    if(layout == NULL || layout->num_slabs <= 0){
        return ST_FOM_ERR_ARGUMENT;
    }

    return ST_fom_extract_slab_right_trace(
        layout,
        layout->num_slabs - 1,
        window_state,
        trace);
}
