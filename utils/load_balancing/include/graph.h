#pragma once

#include "lb_dataset.h"
#include "std.h"
#include "read_BB.h"
#include "write_BB.h"

const char* get_directory_name(
	int         argc,
	char*       argv[],
	const char* codename);

void lb_set_node_weight(
    const int       num_subdomains,
    const char*     directory);

void lb_set_node_weight_for_hyperreduction(
    const int       num_subdomains,
    const char*     directory);

const char* hlpod_get_directory_name(
    const char* directory);

void hlpod_set_subdomain_graph(
    const int       num_subdomains,
    const char*     directory);

void open_gedatsu_external_file(
    const int       nd1,
    const char*  	directory);

void lb_read_meta_graph_mpi(
    const int       nd1,
    const char*     directory,
    LPOD_LB*        lpod_lb);