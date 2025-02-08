#include "mex.h"
#include "half.hpp"
#include <complex>

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // Check input arguments
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:complex_matrix_multiplication:nrhs", "Two inputs required.");
    }

    // Check output arguments
    if (nlhs != 3) {
        mexErrMsgIdAndTxt("MyToolbox:complex_matrix_multiplication:nlhs", "Three outputs required.");
    }

    // Create input matrices from MATLAB data
    double* input_matrix_a_real = mxGetPr(prhs[0]);
    double* input_matrix_a_imag = mxGetPi(prhs[0]);
    double* input_matrix_b_real = mxGetPr(prhs[1]);
    double* input_matrix_b_imag = mxGetPi(prhs[1]);
    mwSize m = mxGetM(prhs[0]);  // Number of rows of A
    mwSize n = mxGetN(prhs[0]);  // Number of columns of A
    mwSize k = mxGetN(prhs[1]);  // Number of columns of B

    // Check if the matrix dimensions are compatible for multiplication
    if (n != mxGetM(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:complex_matrix_multiplication:dim", "Matrix dimensions are not compatible for multiplication.");
    }

    // Create output matrices
    plhs[2] = mxCreateDoubleMatrix(m, k, mxCOMPLEX);
    double* output_matrix_product_real = mxGetPr(plhs[2]);
    double* output_matrix_product_imag = mxGetPi(plhs[2]);

    plhs[0] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
    double* a_halfprecision_real = mxGetPr(plhs[0]);
    double* a_halfprecision_imag = mxGetPi(plhs[0]);

    plhs[1] = mxCreateDoubleMatrix(n, k, mxCOMPLEX);
    double* b_halfprecision_real = mxGetPr(plhs[1]);
    double* b_halfprecision_imag = mxGetPi(plhs[1]);

    // Perform complex matrix multiplication with half-precision conversion
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            half_float::half sum_real = static_cast<half_float::half>(0.0);
            half_float::half sum_imag = static_cast<half_float::half>(0.0);
            for (int l = 0; l < n; l++) {
                half_float::half a_real = static_cast<half_float::half>(input_matrix_a_real[i + m * l]);
                half_float::half a_imag = static_cast<half_float::half>(input_matrix_a_imag[i + m * l]);
                half_float::half b_real = static_cast<half_float::half>(input_matrix_b_real[l + n * j]);
                half_float::half b_imag = static_cast<half_float::half>(input_matrix_b_imag[l + n * j]);

                // Complex multiplication
                half_float::half product_real = half_float::operator-(half_float::operator*(a_real,b_real), half_float::operator*(a_imag,b_imag));
                half_float::half product_imag = half_float::operator+(half_float::operator*(a_real,b_imag), half_float::operator*(a_imag,b_real));

                sum_real = half_float::operator+(sum_real, product_real);
                sum_imag = half_float::operator+(sum_imag, product_imag);
                // Store intermediate results
                a_halfprecision_real[i + m * l] = static_cast<double>(a_real);
                a_halfprecision_imag[i + m * l] = static_cast<double>(a_imag);
                b_halfprecision_real[l + n * j] = static_cast<double>(b_real);
                b_halfprecision_imag[l + n * j] = static_cast<double>(b_imag);
            }
            output_matrix_product_real[i + m * j] = static_cast<double>(sum_real);
            output_matrix_product_imag[i + m * j] = static_cast<double>(sum_imag);
        }
    }
}
