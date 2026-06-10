#include "mex.h"
#include "half.hpp"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // Check input arguments
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:matrix_division:nrhs", "Two inputs required.");
    }

    // Check output arguments
    if (nlhs != 3) {
        mexErrMsgIdAndTxt("MyToolbox:matrix_division:nlhs", "Three outputs required.");
    }

    // Create input matrices from MATLAB data
    double* input_matrix = mxGetPr(prhs[0]);
    double divisor = mxGetScalar(prhs[1]);  // Get the divisor as a double

    mwSize m = mxGetM(prhs[0]);  // Number of rows
    mwSize n = mxGetN(prhs[0]);  // Number of columns

    // Create output matrices
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    plhs[1] = mxCreateDoubleScalar(static_cast<double>(divisor));
    plhs[2] = mxCreateDoubleMatrix(m, n, mxREAL);

    double* result_matrix = mxGetPr(plhs[0]);
    double* divisor_halfprecision = mxGetPr(plhs[1]);
    double* result_halfprecision = mxGetPr(plhs[2]);
    *divisor_halfprecision = static_cast<double>(half_float::half(divisor));
    
    // Perform matrix division with half-precision conversion
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            half_float::half a = static_cast<half_float::half>(input_matrix[i + m * j]);
            half_float::half b = static_cast<half_float::half>(divisor);
            half_float::half value=half_float::operator*(a ,b)	;
            result_matrix[i + m * j] = static_cast<double>(a);
            result_halfprecision[i + m * j] = static_cast<double>(value);
        }
    }
}
