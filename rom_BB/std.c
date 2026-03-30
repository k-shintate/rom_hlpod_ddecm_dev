#include "std.h"

double ROM_BB_sum(
    const int start,
    const int end,
    double* x)
{
    double sum = 0.0;
    for (int i = start; i < end; i++){
        sum += x[i];
    }
    return sum;
}


void ROM_BB_matvec(
    double** A,
    double* x,
    double* y,
    const int rows,
    const int cols)
{
    for (int i = 0; i < rows; i++) {
        double sum = 0.0;
        for (int j = 0; j < cols; j++) {
            sum += A[i][j] * x[j];
        }
        y[i] = sum;
    }
}

void ROM_BB_matmat(
    double** A,
    double** B,
    double** C,
    const int m,
    const int k,
    const int n)
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int p = 0; p < k; p++) {
                sum += A[i][p] * B[p][j];
            }
            C[i][j] = sum;
        }
    }
}

void ROM_BB_transposemat_mat(
    double** A,  // A: m x k
    double** B,  // B: m x n
    double** C,  // C: k x n
    const int m,
    const int k,
    const int n)
{
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int p = 0; p < m; p++) {
                sum += A[p][i] * B[p][j]; // (A^T * B)[i][j] = Σ_p A[p][i] * B[p][j]
            }
            C[i][j] = sum;
        }
    }
}

void ROM_BB_transposemat_vec(
    double** A,  // A: m x n
    double* b,   // b: m
    double* y,   // y: n
    const int m,
    const int n)
{
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < m; j++) {
            sum += A[j][i] * b[j]; // y[i] = Σ_j A[j][i] * b[j]
        }
        y[i] = sum;
    }
}

void ROM_BB_vec_copy(
		double*  in,	//input
		double*  out,	//output
		const int num)
{
    for (int i = 0; i < num; i++) {
       out[i] = in[i];
    }
}

void ROM_BB_vec_copy_2d_to_1d(
		double**  in,
		double*   out,
		const int num)
{
	for(int i = 0; i < num; i++) {
		for(int d = 0; d < 3; d++) {
			out[3*i + d] = in[i][d] ;
		}
	}
}

// 配列の最大値を見つける関数
int ROM_BB_findMax(
    int array[],
    const int size) {
    int max = array[0];  // 最初の要素を最大値とする

    for (int i = 1; i < size; i++) {
        if (array[i] > max) {
            max = array[i];  // 新しい最大値を更新
        }
    }

    return max;
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

//二部探索でサーチ
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
    
    printf("\nBinary search error\n value = %d\n", x);
    exit(1);
}

int ROM_BB_estimate_num_pod_modes(
    double* V,
    const int num_modes_max,
    const double rom_epsilon)
{
    double si;
    int val;
    double ep = 1 - rom_epsilon;
    double sum = ROM_BB_sum(0,num_modes_max,V);
    for(int i = 0; i < num_modes_max; i++){
        si = ROM_BB_sum(0,i,V) / sum;
        if(si > ep){
		val = i;
		printf("val = %d 1 - rom_epsilon = %lf", val, ep);
            break;
        }
    }
    if(val==0){
        return 1;
    }
    else{
        return val;
    }
}


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


int ROM_BB_gauss_elimination(
        int n,
        double **A,
        double *b,
        double *x)
{
    for(int i = 0; i < n; i++) {
		// Pivot selection (partial pivoting)
		int max_row = i;
		double max_val = fabs(A[i][i]);
        for(int k = i+1; k < n; k++) {
            if(fabs(A[k][i]) > max_val) {
                max_val = fabs(A[k][i]);
                max_row = k;
            }
        }

		// If pivot is close to 0, solution does not exist or infinite solutions exist
		if(max_val < 1e-12) {
			return -1; // No solution or infinite solutions
        }

		// Row exchange
		if(max_row != i) {
            double *temp_row = A[i];
            A[i] = A[max_row];
            A[max_row] = temp_row;

            double temp_b = b[i];
            b[i] = b[max_row];
            b[max_row] = temp_b;
        }

		// Normalize pivot row
		double pivot = A[i][i];
        for(int j = i; j < n; j++) {
            A[i][j] /= pivot;
        }
        b[i] /= pivot;

		// Eliminate pivot element from other rows
		for(int k = 0; k < n; k++) {
            if(k != i) {
                double factor = A[k][i];
                for(int j = i; j < n; j++) {
                    A[k][j] -= factor * A[i][j];
                }
                b[k] -= factor * b[i];
            }
        }
    }

	// Extract solution (matrix is already diagonal)
	for(int i = 0; i < n; i++) {
        x[i] = b[i];
    }
	
	printf("Gauss elimination succeeded\n");

    return 0;
}
