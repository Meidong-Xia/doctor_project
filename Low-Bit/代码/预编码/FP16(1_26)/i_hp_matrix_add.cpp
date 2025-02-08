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
//     if (!(nlhs == 1 || nlhs == 3)) {
//         mexErrMsgIdAndTxt("MyToolbox:half_conversion:nlhs", "Either one or three outputs required.");
//     }

    // Create input variables from MATLAB data
    double* input_matrix_a_real = mxGetPr(prhs[0]);
    double* input_matrix_b_real = mxGetPr(prhs[1]);
    double* input_matrix_a_imag = mxGetPi(prhs[0]);
    double* input_matrix_b_imag  = mxGetPi(prhs[1]);
    
    mwSize m = mxGetM(prhs[0]);  // Number of rows
    mwSize n = mxGetN(prhs[0]);  // Number of columns

    // 创建输出矩阵
    if (nlhs == 1) {
        // 如果只有一个输出参数，则创建一个输出矩阵
        plhs[0] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
    } else {
        // 如果有三个输出参数，则创建三个输出矩阵
        plhs[0] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
        plhs[1] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
        plhs[2] = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
    }

    double* a_real =  (nlhs == 3) ? mxGetPr(plhs[0]) : nullptr;
    double* b_real = (nlhs == 3) ? mxGetPr(plhs[1]) : nullptr;
    double* output_real = (nlhs == 3) ? mxGetPr(plhs[2]) : mxGetPr(plhs[0]);
    double* a_imag = (nlhs == 3) ? mxGetPi(plhs[0]) : nullptr;
    double* b_imag = (nlhs == 3) ? mxGetPi(plhs[1]) : nullptr;
    double* output_imag = (nlhs == 3) ? mxGetPi(plhs[2]) : mxGetPi(plhs[0]);
    
    // Perform element-wise addition for complex numbers
            half_float::half a_real_half = static_cast<half_float::half>(0.0);
            half_float::half a_imag_half = static_cast<half_float::half>(0.0);
            half_float::half b_real_half = static_cast<half_float::half>(0.0);
            half_float::half b_imag_half = static_cast<half_float::half>(0.0);

            half_float::half sub_real = static_cast<half_float::half>(0.0);
            half_float::half sub_imag = static_cast<half_float::half>(0.0);
    if (nlhs == 3) {
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
    else {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {

            a_real_half = static_cast<half_float::half>(input_matrix_a_real[i + m * j]);
            a_imag_half = static_cast<half_float::half>(input_matrix_a_imag[i + m * j]);
            b_real_half = static_cast<half_float::half>(input_matrix_b_real[i + m * j]);
            b_imag_half = static_cast<half_float::half>(input_matrix_b_imag[i + m * j]);

            sub_real = half_float::operator+(a_real_half,b_real_half);
            sub_imag = half_float::operator+(a_imag_half,b_imag_half); 

            // Set the real and imaginary parts of the output
            output_real[i + m * j] = sub_real;
            output_imag[i + m * j] = sub_imag;
            }
          }   
    }
}
