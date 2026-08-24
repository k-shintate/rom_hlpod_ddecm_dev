static void ST_fom_time_clear(ST_FOM_TIME_ELEMENT* time)
{
    if(time != NULL){
        memset(time, 0, sizeof(*time));
        time->family = ST_FOM_TIME_FAMILY_DG;
    }
}

static void ST_fom_set_flux_matrices(ST_FOM_TIME_ELEMENT* time)
{
    int alpha;
    int beta;

    for(alpha = 0; alpha < time->n_dof; alpha++){
        for(beta = 0; beta < time->n_dof; beta++){
            time->E_self[alpha][beta] =
                time->e_minus[alpha] * time->e_minus[beta];
            time->E_prev[alpha][beta] =
                time->e_minus[alpha] * time->e_plus[beta];
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

    ST_fom_time_clear(time);
    time->degree = degree;
    time->is_nodal = 1;

    if(degree == 0){
        time->n_dof = 1;
        time->dof_tau[0] = 1.0;
        time->Mt[0][0] = 1.0;
        time->Dt[0][0] = 0.0;
        time->e_minus[0] = 1.0;
        time->e_plus[0] = 1.0;
        time->num_time_ip = 1;
        time->tau[0] = 1.0;
        time->weight[0] = 1.0;
        time->psi[0][0] = 1.0;
    }
    else if(degree == 1){
        const double s = 0.57735026918962576451;

        time->n_dof = 2;
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

        time->num_time_ip = 2;
        time->tau[0] = 0.5 * (1.0 - s);
        time->tau[1] = 0.5 * (1.0 + s);
        time->weight[0] = 0.5;
        time->weight[1] = 0.5;
        for(q = 0; q < time->num_time_ip; q++){
            time->psi[0][q] = 1.0 - time->tau[q];
            time->psi[1][q] = time->tau[q];
        }
    }
    else if(degree == 2){
        const double s = 0.77459666924148337704;

        time->n_dof = 3;
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

        time->num_time_ip = 3;
        time->tau[0] = 0.5 * (1.0 - s);
        time->tau[1] = 0.5;
        time->tau[2] = 0.5 * (1.0 + s);
        time->weight[0] = 5.0 / 18.0;
        time->weight[1] = 4.0 / 9.0;
        time->weight[2] = 5.0 / 18.0;
        for(q = 0; q < time->num_time_ip; q++){
            const double tau = time->tau[q];
            time->psi[0][q] = 2.0 * tau * tau - 3.0 * tau + 1.0;
            time->psi[1][q] = -4.0 * tau * tau + 4.0 * tau;
            time->psi[2][q] = 2.0 * tau * tau - tau;
        }
    }
    else{
        ST_fom_time_clear(time);
        return ST_FOM_ERR_UNSUPPORTED;
    }

    ST_fom_set_flux_matrices(time);
    return ST_FOM_SUCCESS;
}

int ST_fom_resolve_owned_space_points(
    int num_local_space_points,
    int communicator_size,
    int communication_owned_space_points,
    int* num_owned_space_points)
{
    if(num_owned_space_points == NULL ||
       num_local_space_points <= 0 || communicator_size <= 0)
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

int ST_fom_layout_initialize(
    ST_FOM_WINDOW_LAYOUT* layout,
    int num_space_points,
    int num_owned_space_points,
    int physical_dof,
    int num_slabs,
    int degree)
{
    int status;

    if(layout == NULL || num_space_points <= 0 ||
       num_owned_space_points <= 0 ||
       num_owned_space_points > num_space_points ||
       physical_dof <= 0 || num_slabs <= 0)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    memset(layout, 0, sizeof(*layout));
    layout->num_space_points = num_space_points;
    layout->num_owned_space_points = num_owned_space_points;
    layout->physical_dof = physical_dof;
    layout->num_slabs = num_slabs;

    status = ST_fom_time_element_init_dg(&layout->time, degree);
    if(status != ST_FOM_SUCCESS){
        memset(layout, 0, sizeof(*layout));
        return status;
    }

    if(ST_fom_window_vector_length(layout) == 0u){
        memset(layout, 0, sizeof(*layout));
        return ST_FOM_ERR_DIMENSION;
    }

    return ST_FOM_SUCCESS;
}

int ST_fom_window_dof_per_space_point(
    const ST_FOM_WINDOW_LAYOUT* layout)
{
    long long value;

    if(layout == NULL || layout->physical_dof <= 0 ||
       layout->num_slabs <= 0 || layout->time.n_dof <= 0)
    {
        return 0;
    }

    value = (long long)layout->physical_dof
        * (long long)layout->num_slabs
        * (long long)layout->time.n_dof;
    if(value <= 0 || value > INT_MAX){
        return 0;
    }
    return (int)value;
}

size_t ST_fom_window_vector_length(
    const ST_FOM_WINDOW_LAYOUT* layout)
{
    const int block_dof = ST_fom_window_dof_per_space_point(layout);

    if(layout == NULL || block_dof <= 0 || layout->num_space_points <= 0){
        return 0u;
    }
    if((size_t)layout->num_space_points > SIZE_MAX / (size_t)block_dof){
        return 0u;
    }
    return (size_t)layout->num_space_points * (size_t)block_dof;
}

int ST_fom_block_dof(
    const ST_FOM_WINDOW_LAYOUT* layout,
    int slab,
    int time_dof,
    int component)
{
    long long value;

    if(layout == NULL || slab < 0 || slab >= layout->num_slabs ||
       time_dof < 0 || time_dof >= layout->time.n_dof ||
       component < 0 || component >= layout->physical_dof)
    {
        return -1;
    }

    value = ((long long)slab * (long long)layout->time.n_dof
        + (long long)time_dof) * (long long)layout->physical_dof
        + (long long)component;
    if(value < 0 || value > INT_MAX){
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
    const int block_dof = ST_fom_window_dof_per_space_point(layout);
    const int local_dof = ST_fom_block_dof(
        layout, slab, time_dof, component);

    if(layout == NULL || block_dof <= 0 || local_dof < 0 ||
       space_point < 0 || space_point >= layout->num_space_points)
    {
        return SIZE_MAX;
    }

    if((size_t)space_point >
       (SIZE_MAX - (size_t)local_dof) / (size_t)block_dof)
    {
        return SIZE_MAX;
    }

    return (size_t)space_point * (size_t)block_dof
        + (size_t)local_dof;
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
       layout->time.degree != 0 || layout->time.n_dof != 1 ||
       slab < 0 || slab >= layout->num_slabs)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    for(node = 0; node < layout->num_space_points; node++){
        for(component = 0; component < layout->physical_dof; component++){
            const size_t index = ST_fom_vector_index(
                layout, node, slab, 0, component);
            if(index == SIZE_MAX){
                return ST_FOM_ERR_DIMENSION;
            }
            window_state[index] = state[
                (size_t)node * (size_t)layout->physical_dof
                + (size_t)component];
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
       layout->time.degree != 0 || layout->time.n_dof != 1 ||
       slab < 0 || slab >= layout->num_slabs)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    for(node = 0; node < layout->num_space_points; node++){
        for(component = 0; component < layout->physical_dof; component++){
            const size_t index = ST_fom_vector_index(
                layout, node, slab, 0, component);
            if(index == SIZE_MAX){
                return ST_FOM_ERR_DIMENSION;
            }
            state[(size_t)node * (size_t)layout->physical_dof
                + (size_t)component] = window_state[index];
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
       slab < 0 || slab >= layout->num_slabs ||
       layout->time.n_dof <= 0)
    {
        return ST_FOM_ERR_ARGUMENT;
    }

    for(node = 0; node < layout->num_space_points; node++){
        for(component = 0; component < layout->physical_dof; component++){
            double value = 0.0;
            for(alpha = 0; alpha < layout->time.n_dof; alpha++){
                const size_t index = ST_fom_vector_index(
                    layout, node, slab, alpha, component);
                if(index == SIZE_MAX){
                    return ST_FOM_ERR_DIMENSION;
                }
                value += layout->time.e_plus[alpha] * window_state[index];
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
    if(layout == NULL){
        return ST_FOM_ERR_ARGUMENT;
    }
    return ST_fom_extract_slab_right_trace(
        layout, layout->num_slabs - 1, window_state, trace);
}

