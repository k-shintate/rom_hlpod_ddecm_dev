
#include "std.h"

void ROM_BB_vec_copy_1d(
        double*  in,   //input
        double*  out,  //output
        const int num1,
        const int num2)
{
    for (int i = 0; i < num1; i++) {
        out[i] = in[i];
    }
}

void ROM_BB_vec_copy_2d(
        double**  in,   //input
        double**  out,  //output
        const int num1,
        const int num2)
{
    for (int i = 0; i < num1; i++) {
        for(int j = 0; j < num2; j++){
            out[i][j] = in[i][j];
        }
    }
}

double calc_internal_norm_1d(
        double*  in,       //input
        const int num1,
        const int num2)
{
    double norm = 0.0;

    for (int i = 0; i < num1; i++) {    
        norm += in[i] * in[i];
    }

    return norm;
}


double _Complex BBFE_std_integ_calc_C(
    const int num_integ_points,
    const double _Complex *value,
    const double *weight,
    const double *Jacobian)
{
    double _Complex val = 0.0;

    for (int i = 0; i < num_integ_points; i++)
    {
        val += value[i] * weight[i] * Jacobian[i];
    }

    return val;
}


double calc_nodal_error_vector(double** edge_result,
    double** source,
    double** edge_error,
    int N)
{
    double sum_sq_diff = 0.0;
    double sum_sq_src  = 0.0;

    for (int i = 0; i < N; ++i) {
        for (int d = 0; d < 3; ++d) {
            double diff = edge_result[i][d] - source[i][d];
            edge_error[i][d] = sqrt(diff*diff);
            sum_sq_diff += diff * diff;
            sum_sq_src  += source[i][d] * source[i][d];
        }
    }

    if (sum_sq_src == 0.0) {
        return (sum_sq_diff == 0.0) ? 0.0 : INFINITY;
    }

    return sqrt(sum_sq_diff) / sqrt(sum_sq_src);
}


void copy_Aphi_to_V_phi(
    double _Complex * Aphi,
    double * V,
    double * phi,
    const int total_num_nodes,
    const int total_num_edges)
{

    for(int i = 0; i < total_num_edges; i++){
        V[i] = creal(Aphi[i + total_num_nodes]);
    }

    for(int i = 0; i < total_num_nodes; i++){
        phi[i] = creal(Aphi[i]);
    }

}

void copy_Aphi_to_V_phi_time(
    double * Aphi,
    double * V,
    double * phi,
    const int total_num_nodes,
    const int total_num_edges)
{

    for(int i = 0; i < total_num_edges; i++){
        V[i] = Aphi[i + total_num_nodes];
    }

    for(int i = 0; i < total_num_nodes; i++){
        phi[i] = Aphi[i];
    }

}


void compute_Js(double x, double y, double z,
    double nu, double sigma, double omega,
    double _Complex Js[3])
{
    double pi = M_PI;
    double pi3 = pi * pi * pi; // π^3

    // 三角関数の値（引数は πx, πy, πz）
    double sin_x = sin(pi * x);
    double sin_y = sin(pi * y);
    double sin_z = sin(pi * z);
    double cos_x = cos(pi * x);
    double cos_y = cos(pi * y);
    double cos_z = cos(pi * z);

    // 各係数項を計算
    // coeff1 = 6*nu*pi^3 + 2*j*omega*sigma*pi
    double _Complex coeff1 = 6.0 * nu * pi3 + 2.0 * I * omega * sigma * pi;
    // coeff2 = 3*nu*pi^3 + j*omega*sigma*pi
    double _Complex coeff2 = 3.0 * nu * pi3 + I * omega * sigma * pi;

    // J_s の各成分の計算
    // Jx = coeff1 * sin(πx) * cos(πy) * cos(πz)
    Js[0] = 2 * x * I * omega * sigma;
    // Jy = - coeff2 * cos(πx) * sin(πy) * cos(πz)
    Js[1] = -2 * y * I * omega * sigma;
    // Jz = - coeff2 * cos(πx) * cos(πy) * sin(πz)
    Js[2] = 2*M_PI*M_PI*nu*sin_x*sin_y + sin_x*sin_y*I*omega*sigma;
}

void compute_Js_time(
  double x,double y,double z,double t,
  double lambda,double Nu,double Sigma,
  double out[3])
{
  (void)z;
  const double Az = sin(M_PI*x)*sin(M_PI*y)*exp(-lambda*t);

  // σ ∇ϕ⋆ の寄与（x,y 成分）
  out[0] =  2.0 * Sigma * x;   // J_s,x
  out[1] = -2.0 * Sigma * y;   // J_s,y

  // z 成分： σ∂t A⋆ + ν curl curl A⋆
  out[2] = Sigma*(-lambda)*Az + Nu*(2.0*M_PI*M_PI)*Az; // J_s,z
}

double**** BB_std_calloc_4d_double(
		double****  array,
		const int  size1,
		const int  size2,
		const int  size3,
        const int  size4)
{
	array = (double****)calloc(size1, sizeof(double***));
	for(int i=0; i<size1; i++) {
		array[i] = (double***)calloc(size2, sizeof(double**));
		
		for(int j=0; j<size2; j++) {
			array[i][j] = (double**)calloc(size3, sizeof(double*));
		
        	for(int k=0; k<size3; k++) {
		    	array[i][j][k] = (double*)calloc(size4, sizeof(double));
		    }
        }
	}

	return array;
}
