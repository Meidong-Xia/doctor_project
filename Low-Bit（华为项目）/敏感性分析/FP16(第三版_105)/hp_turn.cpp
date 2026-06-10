#include "mex.h"
#include "half.hpp"
#include <complex>

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // 检查输入参数
    if (nrhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:half_conversion:nrhs", "Two inputs required.");
    }
    // Check output arguments
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:half_conversion:nlhs", "Three outputs required.");
    }
    // 从MATLAB数据创建输入变量
    double* input_a_real = mxGetPr(prhs[0]);
    double* input_a_imag = mxGetPi(prhs[0]);

    // 创建输出矩阵
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
    double* a_real = mxGetPr(plhs[0]);
    double* a_imag = mxGetPi(plhs[0]);

    half_float::half a_imag_half = static_cast<half_float::half>(0.0);   
    half_float::half a_real_half = static_cast<half_float::half>(*input_a_real);
    // 处理第一个输入
    if (!mxIsComplex(prhs[0])) {
        // 第一个输入是实数
        a_imag_half = static_cast<half_float::half>(0.0);
    } else {
        // 第一个输入是复数
        a_imag_half = static_cast<half_float::half>(*input_a_imag);
    }
    // 计算
    *a_real = a_real_half;
    *a_imag = a_imag_half;
}
