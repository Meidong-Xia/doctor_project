#include "mex.h"
#include "half.hpp"
#include <complex>

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // Check input arguments
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:complex_division:nrhs", "Two inputs required.");
    }

    // Check output arguments
    if (nlhs != 3) {
        mexErrMsgIdAndTxt("MyToolbox:complex_division:nlhs", "Three outputs required.");
    }

    // Create input variables from MATLAB data
    double* input_a_real = mxGetPr(prhs[0]);
    double* input_a_imag = mxGetPi(prhs[0]);
    double* input_b_real = mxGetPr(prhs[1]);
    double* input_b_imag = mxGetPi(prhs[1]);
    
    // Create output matrices
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
    
    double* a_real = mxGetPr(plhs[0]);
    double* b_real = mxGetPr(plhs[1]);
    double* output_real = mxGetPr(plhs[2]);
    double* a_imag = mxGetPi(plhs[0]);
    double* b_imag = mxGetPi(plhs[1]);
    double* output_imag = mxGetPi(plhs[2]);
    
    // Perform half-precision conversion
    half_float::half a_imag_half = static_cast<half_float::half>(0.0);
    half_float::half b_imag_half = static_cast<half_float::half>(0.0);    
    half_float::half a_real_half = static_cast<half_float::half>(*input_a_real);
    half_float::half b_real_half = static_cast<half_float::half>(*input_b_real);
    half_float::half result_real = static_cast<half_float::half>(0.0);   
    half_float::half result_imag = static_cast<half_float::half>(0.0);   
    double panduan = 0;
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
        panduan = 1;
        
    } else {
        // 第二个输入是复数
        b_imag_half = static_cast<half_float::half>(*input_b_imag);
    }
    
    if(panduan==0){
    // Perform complex division
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
    else  {
    result_real = static_cast<half_float::half>(half_float::operator/(a_real_half,b_real_half));
    result_imag = static_cast<half_float::half>(half_float::operator/(a_imag_half,b_real_half));
    }
    *a_real = a_real_half;
    *b_real = b_real_half;
    *output_real = result_real;
    *a_imag = a_imag_half;
    *b_imag = b_imag_half;
    *output_imag = result_imag;
}
