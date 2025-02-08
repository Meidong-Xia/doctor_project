#include "mex.h"
#include "half.hpp"
#include <complex>

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // Check input arguments
    if (nrhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:matrix_addition:nrhs", "one inputs required.");
    }

    // Check output arguments
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:matrix_addition:nlhs", "one outputs required.");
    }

    // Create input variables from MATLAB data
    double* input_matrix_a_real = mxGetPr(prhs[0]);
    double* input_matrix_a_imag = mxGetPi(prhs[0]);
    
    mwSize m = mxGetM(prhs[0]);  // Number of rows
    mwSize n = mxGetN(prhs[0]);  // Number of columns

    // Create output matrices
    plhs[0] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);

    double* a_real = mxGetPr(plhs[0]);
    double* a_imag = mxGetPi(plhs[0]);
    
    // Perform element-wise addition for complex numbers
            half_float::half a_real_half = static_cast<half_float::half>(0.0);
            half_float::half a_imag_half = static_cast<half_float::half>(0.0);
    // Perform element-wise subtraction for complex numbers
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {

            a_real_half = static_cast<half_float::half>(input_matrix_a_real[i + m * j]);
            a_imag_half = static_cast<half_float::half>(input_matrix_a_imag[i + m * j]);

            // Set the real and imaginary parts of the output
            a_real[i + m * j] = a_real_half;
            a_imag[i + m * j] = a_imag_half;
    }
}

}
