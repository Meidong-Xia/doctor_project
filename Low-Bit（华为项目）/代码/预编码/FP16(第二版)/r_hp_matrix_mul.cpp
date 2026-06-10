#include "mex.h"
#include "half.hpp"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // Check input arguments
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:matrix_multiplication:nrhs", "Two inputs required.");
    }

    // Check output arguments
    if (nlhs != 3) {
        mexErrMsgIdAndTxt("MyToolbox:matrix_multiplication:nlhs", "Three outputs required.");
    }

    // Create input matrices from MATLAB data
    double* input_matrix_a = mxGetPr(prhs[0]);
    double* input_matrix_b = mxGetPr(prhs[1]);
    mwSize m = mxGetM(prhs[0]);  // Number of rows of A
    mwSize n = mxGetN(prhs[0]);  // Number of columns of A
    mwSize k = mxGetN(prhs[1]);  // Number of columns of B

    // Check if the matrix dimensions are compatible for multiplication
    if (n != mxGetM(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:matrix_multiplication:dim", "Matrix dimensions are not compatible for multiplication.");
    }

    // Create output matrices
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n, k, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(m, k, mxREAL);

    double* a_halfprecision = mxGetPr(plhs[0]);
    double* b_halfprecision = mxGetPr(plhs[1]);
    double* output_matrix_product = mxGetPr(plhs[2]);

    // Perform matrix multiplication with half-precision conversion
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < k; j++) {
            half_float::half sum = static_cast<half_float::half>(0.0);
            for (int l = 0; l < n; l++) {
                half_float::half a = static_cast<half_float::half>(input_matrix_a[i + m * l]);
                half_float::half b = static_cast<half_float::half>(input_matrix_b[l + n * j]);
                sum = half_float::operator+(half_float::operator*(a,b),sum);
//                 sum += a * b;            
                a_halfprecision[i + m * l] = static_cast<double>(a);
                b_halfprecision[l + n * j] = static_cast<double>(b);
            }
            output_matrix_product[i + m * j] = static_cast<double>(sum);
        }
    }
}
