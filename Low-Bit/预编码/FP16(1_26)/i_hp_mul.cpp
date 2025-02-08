#include "mex.h"
#include "half.hpp"
#include <complex>

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // Check input arguments
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:complex_multiplication:nrhs", "Two inputs required.");
    }

    // Check output arguments
    if (!(nlhs == 1 || nlhs == 3)) {
        mexErrMsgIdAndTxt("MyToolbox:half_conversion:nlhs", "Either one or three outputs required.");
    }

    // Create input variables from MATLAB data
    double* input_a_real = mxGetPr(prhs[0]);
    double* input_a_imag = mxGetPi(prhs[0]);
    double* input_b_real = mxGetPr(prhs[1]);
    double* input_b_imag = mxGetPi(prhs[1]);
    
   // Create output matrices
        // 创建输出矩阵
    if (nlhs == 1) {
        // 如果只有一个输出参数，则创建一个输出矩阵
        plhs[0] = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
    } else {
        // 如果有三个输出参数，则创建三个输出矩阵
        plhs[0] = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
        plhs[1] = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
        plhs[2] = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
    }
    
    double* a_real =  (nlhs == 3) ? mxGetPr(plhs[0]) : nullptr;
    double* b_real = (nlhs == 3) ? mxGetPr(plhs[1]) : nullptr;
    double* output_real = (nlhs == 3) ? mxGetPr(plhs[2]) : mxGetPr(plhs[0]);
    double* a_imag = (nlhs == 3) ? mxGetPi(plhs[0]) : nullptr;
    double* b_imag = (nlhs == 3) ? mxGetPi(plhs[1]) : nullptr;
    double* output_imag = (nlhs == 3) ? mxGetPi(plhs[2]) : mxGetPi(plhs[0]);
    
    // Perform half-precision conversion
    half_float::half a_imag_half = static_cast<half_float::half>(0.0);
    half_float::half b_imag_half = static_cast<half_float::half>(0.0);    
    half_float::half a_real_half = static_cast<half_float::half>(*input_a_real);
    half_float::half b_real_half = static_cast<half_float::half>(*input_b_real);
    // 处理第一个输入
    if (!mxIsComplex(prhs[0])) {
        // 第一个输入是实数
        a_imag_half = static_cast<half_float::half>(0.0);
    } else {
        // 第一个输入是复数
        a_imag_half = static_cast<half_float::half>(*input_a_imag);
    }

    // 处理第二个输入
    if (!mxIsComplex(prhs[1])) {
        // 第二个输入是实数
        b_imag_half = static_cast<half_float::half>(0.0);
    } else {
        // 第二个输入是复数
        b_imag_half = static_cast<half_float::half>(*input_b_imag);
    }

    // Perform complex multiplication
    half_float::half temp1=half_float::operator*(a_real_half,b_real_half);
    half_float::half temp2=half_float::operator*(a_imag_half,b_imag_half);
    half_float::half temp3=half_float::operator*(a_real_half,b_imag_half);
    half_float::half temp4=half_float::operator*(a_imag_half,b_real_half);
    
    half_float::half result_real = half_float::operator-(temp1,temp2);
    half_float::half result_imag = half_float::operator+(temp3,temp4);

        // 计算
    if (nlhs == 3) {
    *a_real = a_real_half;
    *b_real = b_real_half;
    *output_real = result_real;
    *a_imag = a_imag_half;
    *b_imag = b_imag_half;
    *output_imag = result_imag;
    } else {
        // 如果只有一个输出参数，则将结果存储在单个输出矩阵中
        *output_real = result_real;
        *output_imag = result_imag;
    }
    
}
