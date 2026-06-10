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
//     if (nlhs != 3) {
//         mexErrMsgIdAndTxt("MyToolbox:complex_matrix_multiplication:nlhs", "Three outputs required.");
//     }

    // Create input matrices from MATLAB data
    double* input_matrix_a_real = mxGetPr(prhs[0]);
    double* input_matrix_a_imag = mxGetPi(prhs[0]);
    double* input_b_real = mxGetPr(prhs[1]);
    double* input_b_imag = mxGetPi(prhs[1]);
    
    
    mwSize m = mxGetM(prhs[0]);  // Number of rows of A
    mwSize n = mxGetN(prhs[0]);  // Number of columns of A

    // Create output matrices
    
        if (nlhs == 1) {
plhs[0] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
    
    double* output_matrix_product_real = mxGetPr(plhs[0]);
    double* output_matrix_product_imag = mxGetPi(plhs[0]);
    
    half_float::half b_real = static_cast<half_float::half>(*input_b_real);
    half_float::half b_imag = static_cast<half_float::half>(*input_b_imag);
    
    // Perform complex matrix multiplication with half-precision conversion
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
                half_float::half a_real = static_cast<half_float::half>(input_matrix_a_real[i + m * j]);
                half_float::half a_imag = static_cast<half_float::half>(input_matrix_a_imag[i + m * j]);

                // Complex multiplication
                half_float::half product_real = half_float::operator-(half_float::operator*(a_real,b_real), half_float::operator*(a_imag,b_imag));
                half_float::half product_imag = half_float::operator+(half_float::operator*(a_real,b_imag), half_float::operator*(a_imag,b_real));

                // Store intermediate results
                output_matrix_product_real[i + m * j] = static_cast<double>(product_real);
                output_matrix_product_imag[i + m * j] = static_cast<double>(product_imag);
        }
    }
        }
        else{
    plhs[0] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
    plhs[2] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
    
    double* output_matrix_product_real = mxGetPr(plhs[2]);
    double* output_matrix_product_imag = mxGetPi(plhs[2]);
    double* a_halfprecision_real = mxGetPr(plhs[0]);
    double* a_halfprecision_imag = mxGetPi(plhs[0]);
    double* b_halfprecision_real = mxGetPr(plhs[1]);
    double* b_halfprecision_imag = mxGetPi(plhs[1]);
    
    half_float::half b_real = static_cast<half_float::half>(*input_b_real);
    half_float::half b_imag = static_cast<half_float::half>(*input_b_imag);
    
    *b_halfprecision_real = b_real;
    *b_halfprecision_imag = b_imag;
    
    // Perform complex matrix multiplication with half-precision conversion
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
                half_float::half a_real = static_cast<half_float::half>(input_matrix_a_real[i + m * j]);
                half_float::half a_imag = static_cast<half_float::half>(input_matrix_a_imag[i + m * j]);

                // Complex multiplication
                half_float::half product_real = half_float::operator-(half_float::operator*(a_real,b_real), half_float::operator*(a_imag,b_imag));
                half_float::half product_imag = half_float::operator+(half_float::operator*(a_real,b_imag), half_float::operator*(a_imag,b_real));

                // Store intermediate results
                a_halfprecision_real[i + m * j] = static_cast<double>(a_real);
                a_halfprecision_imag[i + m * j] = static_cast<double>(a_imag);
                output_matrix_product_real[i + m * j] = static_cast<double>(product_real);
                output_matrix_product_imag[i + m * j] = static_cast<double>(product_imag);
        }
    }
        }
    
}
