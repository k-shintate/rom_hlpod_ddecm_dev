#include "std.h"

//3d int型calloc
int*** ROM_BB_calloc_3d_int(
    int***  array,
    const int  size1,
    const int  size2,
    const int  size3)
{
	array = (int***)calloc(size1, sizeof(int**));
	for(int i=0; i<size1; i++) {
		array[i] = (int**)calloc(size2, sizeof(int*));
		
		for(int j=0; j<size2; j++) {
			array[i][j] = (int*)calloc(size3, sizeof(int));
		}
	}

	return array;
}

//3d int型calloc
bool*** std_calloc_3d_bool(
	bool***  array,
	const int  size1,
    const int  size2,
	const int  size3)
{
	array = (bool***)calloc(size1, sizeof(bool**));
	for(int i=0; i<size1; i++) {
		array[i] = (bool**)calloc(size2, sizeof(bool*));
		
		for(int j=0; j<size2; j++) {
			array[i][j] = (bool*)calloc(size3, sizeof(bool));
		}
	}

	return array;
}


//バブルソート

void std_bubble_sort(
    int*    vec,
    int     total_num)
{
    int tmp;
    for (int i = 0; i < total_num; ++i) {
        for (int j = i + 1; j < total_num; ++j) {
            if (vec[i] > vec[j]) {
                tmp =  vec[i];
                vec[i] = vec[j];
                vec[j] = tmp;
            }
        }
    }
}

//逆引き用ソート関数
void ROM_BB_bubble_sort_with_id(
    int*    vec,
    int*    id,
    int     total_num)
{
    int tmp_vec;
    int tmp_id;
    for (int i = 0; i < total_num; ++i) {
        for (int j = i + 1; j < total_num; ++j) {
            if (vec[i] > vec[j]) {
                tmp_vec =  vec[i];
                vec[i] = vec[j];
                vec[j] = tmp_vec;

                tmp_id =  id[i];
                id[i] = id[j];
                id[j] = tmp_id;
            }
        }
    }
}

//バブルソート後に重複を削除する関数
//安定動作のために使用
void std_bubble_sort_remove_duplicate_node(
    int*    vec,
    int     total_num)
{
    int tmp;
    for (int i = 0; i < total_num; ++i) {
        for (int j = i + 1; j < total_num; ++j) {
            if (vec[i] == vec[j] && vec[i] != -1) {
                while(vec[i] == vec[j]){
                    vec[j] = -1;
                    j++;
                }
            }
        }
    }
}

//削除した合計値を返す
void std_bubble_sort_remove_duplicate_node_return(
    int*    vec,
    int     total_num,
    int*    return_val)
{
    int tmp;
    int sum = 0;
    for (int i = 0; i < total_num; ++i) {
        for (int j = i + 1; j < total_num; ++j) {
            if (vec[i] == vec[j] && vec[i] != -1) {
                while(vec[i] == vec[j]){
                    vec[j] = -1;
                    j++;
                    sum++;
                }
            }
        }
    }

    *return_val = sum;
}

//削除した値を数える
void std_count_negative_val(
    int*    vec,
    int     total_num,
    int*    return_val)
{
    int sum = 0;
    for (int i = 0; i < total_num; ++i) {
        if (vec[i] == -1){
            sum++;
        }
    }

    *return_val = sum;
}


/*バブルソート後に重複を削除し、
その値を除き、内部要素に差し込むために保存する関数*/
void std_bubble_sort_remove_duplicate_elem(
    int*    vec,
    int*    savings,
    int*    return_val,
    int     total_num)
{
    int tmp;
    int index = 0;
    for (int i = 0; i < total_num; ++i) {
        for (int j = i + 1; j < total_num; ++j) {
            if (vec[i] == vec[j] && vec[i] != -1) {
                while(vec[i] == vec[j]){
                    vec[j] = -1;
                    j++;
                }
            
                savings[index] = vec[i];
                vec[i] = -1;
                index++;
            }
        }
    }
    *return_val = index;
}

//逆引き用関数
void std_bubble_sort_remove_duplicate_elem_with_id(
    int*    vec,
    int*    savings_vec,
    int*    return_val,
    int     total_num,
    int*    id,
    int*    savings_id)
{
    int tmp;
    int index = 0;
    for (int i = 0; i < total_num; ++i) {
        for (int j = i + 1; j < total_num; ++j) {
            if (vec[i] == vec[j] && vec[i] != -1) {
                while(vec[i] == vec[j]){
                    vec[j] = -1;
                    j++;
                }

                savings_vec[index] = vec[i];
                vec[i] = -1;

                savings_id[index] = id[i];
                index++;
            }
        }
    }

    *return_val = index;
}



//二部探索でサーチするプログラム
int ROM_BB_binarySearch(
    int*    a,
    int     x, 
    const int N)
{
    int left = 0;
    int right = N - 1;

    int mid;

    while(left <= right) {
        mid = (left + right) / 2;
        if(a[mid] == x) {
            return mid;
        }
        else if(a[mid] < x) {
            left = mid + 1;
        }
        else {
            right = mid - 1;
        }
    }
    
    //見つからなければ抜ける x = 0のときにバグ
    if(x == 0){
        return 0;
    }
    printf("\nBinary search error\n value = %d\n", x);
    exit(1);
}