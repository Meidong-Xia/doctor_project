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
    
    mwSize m = mxGetM(prhs[0]);  // Number of rows
    mwSize n = mxGetN(prhs[0]);  // Number of columns

    // Create output matrices
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(m, n, mxREAL);

    double* a_real = mxGetPr(plhs[0]);
    double* b_real = mxGetPr(plhs[1]);
    double* output_real = mxGetPr(plhs[2]);
    
    // Perform element-wise addition for complex numbers
    // Perform element-wise subtraction for complex numbers
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {

            half_float::half a_real_half = static_cast<half_float::half>(input_matrix_a_real[i + m * j]);
            half_float::half b_real_half = static_cast<half_float::half>(input_matrix_b_real[i + m * j]);

            half_float::half sub_real = half_float::operator-(a_real_half,b_real_half);

            // Set the real and imaginary parts of the output
            a_real[i + m * j] = a_real_half;
            b_real[i + m * j] = b_real_half;
            output_real[i + m * j] = sub_real;
    }
}

}
