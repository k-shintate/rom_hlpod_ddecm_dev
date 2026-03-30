#pragma once

#include <stdio.h>
#include <stdlib.h>

int*** ROM_BB_calloc_3d_int(
	int***  array,
	const int  size1,
	const int  size2,
	const int  size3);

bool*** std_calloc_3d_bool(
	bool***  array,
	const int  size1,
	const int  size2,
	const int  size3);

void std_bubble_sort(
    int* vec,
    int total_num);

void std_bubble_sort_remove_duplicate_node(
    int* vec,
    int total_num);

void std_bubble_sort_remove_duplicate_node_return(
    int* vec,
    int total_num,
    int* return_val);

void std_bubble_sort_remove_duplicate_elem(
    int* vec,
    int* savings,
    int* return_val,
    int total_num);

void ROM_BB_bubble_sort_with_id(
    int* vec,
    int* id,
    int total_num);

void std_bubble_sort_remove_duplicate_elem_with_id(
    int* vec,
    int* savings_vec,
    int* return_val,
    int total_num,
    int* id,
    int* savings_id);

void std_count_negative_val(
    int* vec,
    int total_num,
    int* return_val);

int ROM_BB_binarySearch(
    int* a,
    int x, 
    const int N);