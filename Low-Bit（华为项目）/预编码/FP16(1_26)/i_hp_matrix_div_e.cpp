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
    
    half_float::half b_real_half = static_cast<half_float::half>(*input_b_real);
    half_float::half b_imag_half = static_cast<half_float::half>(0.0);

    double panduan = 0;
    // 处理第二个输入
    if (!mxIsComplex(prhs[1])) {
        // 第二个输入是实数
        b_imag_half = static_cast<half_float::half>(0.0);
        panduan = 1;
        
    } else {
        // 第二个输入是复数
        b_imag_half = static_cast<half_float::half>(*input_b_imag);
    }
    
    half_float::half result_real = static_cast<half_float::half>(0.0);
    half_float::half result_imag = static_cast<half_float::half>(0.0);
    // Perform complex matrix multiplication with half-precision conversion
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
                half_float::half a_real_half = static_cast<half_float::half>(input_matrix_a_real[i + m * j]);
                half_float::half a_imag_half = static_cast<half_float::half>(input_matrix_a_imag[i + m * j]);

                // Complex multiplication
    // Perform complex division
                if(panduan==0){
    half_float::half temp5=half_float::operator*(b_real_half,b_real_half);
    half_float::half temp6=half_float::operator*(b_imag_half,b_imag_half);
    half_float::half denominator = half_float::operator+(temp5,temp6);
    half_float::half temp1=half_float::operator*(a_real_half,b_real_half);
    half_float::half temp2=half_float::operator*(a_imag_half,b_imag_half);
    half_float::half temp3=half_float::operator*(a_imag_half,b_real_half);
    half_float::half temp4=half_float::operator*(a_real_half,b_imag_half);
    half_float::half molecule1 = half_float::operator+(temp1,temp2);
    half_float::half molecule2 = half_float::operator-(temp3,temp4);
    result_real = half_float::operator/(molecule1,denominator);
    result_imag = half_float::operator/(molecule2,denominator);
                }
                else {
    result_real = half_float::operator/(a_real_half,b_real_half);
    result_imag = half_float::operator/(a_imag_half,b_real_half);
                }
    
                // Store intermediate results
                output_matrix_product_real[i + m * j] = static_cast<double>(result_real);
                output_matrix_product_imag[i + m * j] = static_cast<double>(result_imag);
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
    
    half_float::half b_real_half = static_cast<half_float::half>(*input_b_real);
    half_float::half b_imag_half = static_cast<half_float::half>(0.0);

    
    double panduan = 0;
    // 处理第二个输入
    if (!mxIsComplex(prhs[1])) {
        // 第二个输入是实数
        b_imag_half = static_cast<half_float::half>(0.0);
        panduan = 1;
        
    } else {
        // 第二个输入是复数
        b_imag_half = static_cast<half_float::half>(*input_b_imag);
    }
    
    *b_halfprecision_real = b_real_half;
    *b_halfprecision_imag = b_imag_half;
    half_float::half result_real = static_cast<half_float::half>(0.0);
    half_float::half result_imag = static_cast<half_float::half>(0.0);
    // Perform complex matrix multiplication with half-precision conversion
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
                half_float::half a_real_half = static_cast<half_float::half>(input_matrix_a_real[i + m * j]);
                half_float::half a_imag_half = static_cast<half_float::half>(input_matrix_a_imag[i + m * j]);

                // Complex multiplication
    // Perform complex division
                if(panduan==0){
    half_float::half temp5=half_float::operator*(b_real_half,b_real_half);
    half_float::half temp6=half_float::operator*(b_imag_half,b_imag_half);
    half_float::half denominator = half_float::operator+(temp5,temp6);
    half_float::half temp1=half_float::operator*(a_real_half,b_real_half);
    half_float::half temp2=half_float::operator*(a_imag_half,b_imag_half);
    half_float::half temp3=half_float::operator*(a_imag_half,b_real_half);
    half_float::half temp4=half_float::operator*(a_real_half,b_imag_half);
    half_float::half molecule1 = half_float::operator+(temp1,temp2);
    half_float::half molecule2 = half_float::operator-(temp3,temp4);
    result_real = half_float::operator/(molecule1,denominator);
    result_imag = half_float::operator/(molecule2,denominator);
                }
                else {
    result_real = half_float::operator/(a_real_half,b_real_half);
    result_imag = half_float::operator/(a_imag_half,b_real_half);
                }
    
                // Store intermediate results
                a_halfprecision_real[i + m * j] = static_cast<double>(a_real_half);
                a_halfprecision_imag[i + m * j] = static_cast<double>(a_imag_half);
                output_matrix_product_real[i + m * j] = static_cast<double>(result_real);
                output_matrix_product_imag[i + m * j] = static_cast<double>(result_imag);
        }
    }

    }
}
