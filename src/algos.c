#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#pragma STDC CX_LIMITED_RANGE on



// Функции двойной точности с действительными числами ------------------
double* dallocation_1d(int size) {
	double* array_1d = (double*) malloc(size * sizeof(double));
	return array_1d;
}

double** dallocation_2d(int nrow, int ncol) {
	int row = 0;
	double** array_2d = (double**) malloc(nrow * sizeof(double));
	for (row; row < nrow; row++) {
		array_2d[row] = dallocation_1d(ncol);
	}
	return array_2d;
}


double dvec_norm_L2(double* vec, int size) {
	int index = 0;
	double result_L2 = 0.0;
	for (index; index < size; index++) {
		result_L2 += vec[index] * vec[index];
	}
	result_L2 = sqrt(result_L2);
	return result_L2;
}


// При теплицевой матрице -- 3 * N * log_{2}(N) + 2 * N  = N * (3 log_{2}(N) + 2)
// Сложность вычисления 2 * N^{2}
void dmatrix_dvec_permutation(double** matrix, double* vector, 
							  double* result, int nrow, int ncol, 
							  char transpose) {
	int row = 0, col = 0;
	double dot = 0.0;
	if (transpose == 0 && nrow == ncol) {
		for (row; row < nrow; row++) {
			dot = 0.0;
			for (col = 0; col < ncol; col++) {
				dot += matrix[col][row] * vector[col];
			}
			result[row] = dot; 
		}
	} else {
		for (row; row < nrow; row++) {
			dot = 0.0;
			for (col = 0; col < ncol; col++) {
				dot += matrix[row][col] * vector[col];
			} 
			result[row] = dot;
		}
	}
}


void dcopy_vec(double* donor_array, double* result_array, int size) {
	int index = 0;
	for (index; index < size; index++) {
		result_array[index] = donor_array[index];
	}
}


void dfree_1d(double* arg) {
	free(arg);
	arg = NULL;
}



void dfree_2d(double** arg, int rows) {
	int index = 0;
	for (index; index < rows; index++) {
		dfree_1d(arg[index]);
	}
	free(arg);
	arg = NULL;
}


double ddot(double* array_1, double* array_2, int size) {
	int index = 0;
	double result = 0.0;
	for (index; index < size; index++) {
		result += array_1[index] * array_2[index];
	}
	return result;
}



void dfunciton_array(double (*func) (double), double* darray, double* result, int size){
	int index = 0;	
	for (index; index < size; index++) {
		result[index] = func(darray[index]);
	}
}


void dlinspace(double* result, double low, double high, int size) {
	int index = 0;
	double step = (high - low) / (size - 1);
	for (index; index < size; index++) {
		result[index] = low + index * step;
	}
}


void ddiff_1d(double* vec1, double* vec2, double* result, int size) {
	int index = 0;
	double diff = 0.0;
	for (index; index < size; index++) {
		diff = vec1[index] - vec2[index];
		result[index] = diff;
	}
}


void dsum_1d(double* vec1, double* vec2, double* result, int size) {
	int index = 0;
	double diff = 0.0;
	for (index; index < size; index++) {
		diff = vec1[index] - vec2[index];
		result[index] = diff;
	}
}


void dprod_1d(double* vec1, double* vec2, double* result, int size) {
	int index = 0;
	double diff = 0.0;
	for (index; index < size; index++) {
		diff = vec1[index] * vec2[index];
		result[index] = diff;
	}
}


void ddiv_1d(double* vec1, double* vec2, double* result, int size) {
	int index = 0;
	double diff = 0.0;
	for (index; index < size; index++) {
		diff = vec1[index] / vec2[index];
		result[index] = diff;
	}
}


void dmul_coeff_1d(double* vec, double coeff, double* result, int size) {
	int index = 0;
	for (index; index < size; index++) {
		result[index] = vec[index] * coeff;
	}
}


void drandom_array_1d(double* result, double low, double high, int size, int seed) {
	int index = 0;
	srand(seed);
	for (index; index < size; index++) {
		result[index] = (double) rand() / RAND_MAX * (high - low) + low;
	}
}


void dprint_1d(double* array, int size) {
	int index = 1;
	printf("array: [%lf", array[0]);
	for (index; index < size; index++) {
		printf(", %lf", array[index]);
	}
	printf("]\n");
}


// Интерфейс для модифицированного метода градиентного спуска
void dtwo_step_sgd(double** matrix, double* f_vector, double* u0_vector, 
				  double* result, int size, double eps, int max_iter) {
	int iter = 0, permutes = 0;
	
	double alpha_coeff = 0.0, gamma_coeff = 0.0;	// iter coeffs alpha, gamma
	double a = 0.0, b = 0.0, c = 0.0, det = 0.0; 				// SLE coeffs for alpha, gamma	
	double f_norm = dvec_norm_L2(f_vector, size);
	double* residuals0 = (double *) malloc(size * sizeof(double));
	double* residuals1 = (double *) malloc(size * sizeof(double));
	double* asr = (double *) malloc(size * sizeof(double));
	double* aasr = (double *) malloc(size * sizeof(double));
	double* u1_vector = (double *) malloc(size * sizeof(double));
	double* deltau = (double *) malloc(size * sizeof(double));

	dmatrix_dvec_permutation(matrix, u0_vector, residuals0, size, size, 1); // r_0 := Au_0
	ddiff_1d(residuals0, f_vector, residuals0, size); // r_0 := Au_0 - f
	dmatrix_dvec_permutation(matrix, residuals0, asr, size, size, 0); // Asr := A*r_0
	dmatrix_dvec_permutation(matrix, asr, aasr, size, size, 1); // AAsr := AA*r_0

	alpha_coeff = ddot(asr, asr, size) / ddot(aasr, aasr, size); // alpha := (||A*r||^2) / (||AA*r||^2)
	dmul_coeff_1d(asr, alpha_coeff, asr, size); // A*r := alpha * (A*r)
	ddiff_1d(u0_vector, asr, u1_vector, size); // u1 := u0 - alpha * (A*r)
	permutes += 3;

	for (iter; iter < max_iter; iter++) {
		dmatrix_dvec_permutation(matrix, u1_vector, residuals1, size, size, 1); // r_n := Au_n
		ddiff_1d(residuals1, f_vector, residuals1, size); // r_n := Au_n - f
		ddiff_1d(residuals1, residuals0, residuals0, size);	// deltar = r_0 := r_{n} - r_{n-1}
		dmatrix_dvec_permutation(matrix, residuals1, asr, size, size, 0); // Asr := A*r_n
		dmatrix_dvec_permutation(matrix, asr, aasr, size, size, 1); // AAsr := AA*r_n
		
		a = ddot(residuals0, residuals0, size);
		b = ddot(asr, asr, size);
		c = ddot(aasr, aasr, size);
		det = a * c - b * b;
		alpha_coeff = -(b * b) / det;
		gamma_coeff = a * b / det;

		ddiff_1d(u1_vector, u0_vector, u0_vector, size); // u_0 := u_1 - u_0
		dmul_coeff_1d(u0_vector, alpha_coeff, u0_vector, size); // u_0 := alpha * u_o
		dmul_coeff_1d(asr, gamma_coeff, asr, size); // Asr := gamma * Asr
		ddiff_1d(u1_vector, u0_vector, result, size); // u_2 := u1 - alpha * (u_1 - u_0)
		ddiff_1d(result, asr, result, size); // u_2 := u_1 - alpha * (u_1 - u_0) - gamma * A*r_0
		ddiff_1d(result, u1_vector, deltau, size); // diffu = u_1 := u_2 - u_1
		
		if ((dvec_norm_L2(deltau, size) / f_norm) < eps) {
			break;
		} else {
			dcopy_vec(u1_vector, u0_vector, size);
			dcopy_vec(result, u1_vector, size);
			dcopy_vec(residuals1, residuals0, size);		
		}
	}
	dfree_1d(residuals0);
	dfree_1d(residuals1);
	dfree_1d(asr);
	dfree_1d(aasr);
	dfree_1d(u1_vector);
	dfree_1d(deltau);
}


// Функции теста -----------------------------------------------------
void test_double_functions() {
	int size = 10;
	double min = -5.3, max = 6.7;
	double* random_array = dallocation_1d(size);
	drandom_array_1d(random_array, min, max, size, 124);
	double* linspace_array = dallocation_1d(size);
	dlinspace(linspace_array, min, max, size);
	
	dprint_1d(random_array, size);
	dprint_1d(linspace_array, size);
	dfree_1d(random_array);
	dfree_1d(linspace_array);
	random_array = NULL;
	linspace_array = NULL;
}


void test_twosgd() {
	int size = 3;
	double** matrix = dallocation_2d(size, size);
	matrix[0][0] = 3.4; matrix[0][1] = 1.8; matrix[0][2] = 0.8;
	matrix[1][0] = -5.4; matrix[1][1] = 0.87; matrix[1][2] = 1.87;
	matrix[2][0] = 0.9; matrix[2][1] = 4.56; matrix[2][2] = 9.821;

	double u0_vector[3] = {1., 1., 1.};
	double vec[3] = {1.0, -2.0, 3.0};
	double result_vec_tsgd[3], result_vec[3], result_permutation[3];
	
	dtwo_step_sgd(matrix, vec, u0_vector, result_vec_tsgd, size, 0.00000001, 10000);
	dmatrix_dvec_permutation(matrix, result_vec_tsgd, result_vec, size, size, 1);
	

	dprint_1d(vec, size);
	dprint_1d(result_vec, size);
}


void memtest() {
	int index  = 100000;
	int size = 90000000;
	while(index > 2) {
		double* vec1 = dallocation_1d(size);
		double* vec2 = dallocation_1d(size);
		drandom_array_1d(vec1, -5.5, 5.5, size, index);
		drandom_array_1d(vec2, -5.5, 5.5, size, index + 1);
		printf("iter %d done, dot = %lf\n", index, ddot(vec1, vec2, size));
		dfree_1d(vec1);
		dfree_1d(vec2);
		index -= 1;
	}
}


int main() {
	//test_twosgd();
	memtest();
	return 0;
}
