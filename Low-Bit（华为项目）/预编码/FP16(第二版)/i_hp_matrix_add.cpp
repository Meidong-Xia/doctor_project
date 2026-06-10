#include "mex.h"
#include "half.hpp"
#include <complex>

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // Check input arguments
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:matrix_addition:nrhs", "Two inputs required.");
    }

    // Check output arguments
    if (nlhs != 3) {
        mexErrMsgIdAndTxt("MyToolbox:matrix_addition:nlhs", "Three outputs required.");
    }

    // Create input variables from MATLAB data
    double* input_matrix_a_real = mxGetPr(prhs[0]);
    double* input_matrix_b_real = mxGetPr(prhs[1]);
    double* input_matrix_a_imag = mxGetPi(prhs[0]);
    double* input_matrix_b_imag  = mxGetPi(prhs[1]);
    
    mwSize m = mxGetM(prhs[0]);  // Number of rows
    mwSize n = mxGetN(prhs[0]);  // Number of columns

    // Create output matrices
    plhs[0] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
    plhs[1] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
    plhs[2] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);

    double* a_real = mxGetPr(plhs[0]);
    double* b_real = mxGetPr(plhs[1]);
    double* output_real = mxGetPr(plhs[2]);
    double* a_imag = mxGetPi(plhs[0]);
    double* b_imag = mxGetPi(plhs[1]);
    double* output_imag = mxGetPi(plhs[2]);
    
    // Perform element-wise addition for complex numbers
            half_float::half a_real_half = static_cast<half_float::half>(0.0);
            half_float::half a_imag_half = static_cast<half_float::half>(0.0);
            half_float::half b_real_half = static_cast<half_float::half>(0.0);
            half_float::half b_imag_half = static_cast<half_float::half>(0.0);

            half_float::half sub_real = static_cast<half_float::half>(0.0);
            half_float::half sub_imag = static_cast<half_float::half>(0.0);
    // Perform element-wise subtraction for complex numbers
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {

            a_real_half = static_cast<half_float::half>(input_matrix_a_real[i + m * j]);
            a_imag_half = static_cast<half_float::half>(input_matrix_a_imag[i + m * j]);
            b_real_half = static_cast<half_float::half>(input_matrix_b_real[i + m * j]);
            b_imag_half = static_cast<half_float::half>(input_matrix_b_imag[i + m * j]);

            sub_real = half_float::operator+(a_real_half,b_real_half);
            sub_imag = half_float::operator+(a_imag_half,b_imag_half); 

            // Set the real and imaginary parts of the output
            a_real[i + m * j] = a_real_half;
            b_real[i + m * j] = b_real_half;
            output_real[i + m * j] = sub_real;

            a_imag[i + m * j] = a_imag_half;
            b_imag[i + m * j] = b_imag_half;
            output_imag[i + m * j] = sub_imag;
    }
}

}
