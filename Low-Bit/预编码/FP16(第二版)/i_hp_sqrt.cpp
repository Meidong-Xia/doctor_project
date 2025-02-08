#include "mex.h"
#include "half.hpp"
#include <complex>

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // Check input arguments
    if (nrhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:complex_sqrt:nrhs", "One input required.");
    }

    // Check output arguments
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("MyToolbox:complex_sqrt:nlhs", "One output required.");
    }

    // Create input variable from MATLAB data
    double* input_real = mxGetPr(prhs[0]);
    double* input_imag = mxGetPi(prhs[0]);
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
    
    double* real = mxGetPr(plhs[0]);
    double* image = mxGetPr(plhs[0]);
    double* output_real = mxGetPr(plhs[1]);
    double* output_imag = mxGetPi(plhs[1]);
    
    // Perform half-precision conversion
    half_float::half real_half = static_cast<half_float::half>(*input_real);
    half_float::half imag_half = static_cast<half_float::half>(*input_imag);
    
    // Perform complex square root
    half_float::half temp1=half_float::operator*(real_half,real_half);
    half_float::half temp2=half_float::operator*(imag_half,imag_half);
    half_float::half temp3=half_float::operator+(temp1,temp2);
    half_float::half abs_real_half = half_float::sqrt(temp3);
    half_float::half temp4=half_float::operator+(abs_real_half,real_half);
    half_float::half temp5=half_float::operator-(abs_real_half,real_half);
    half_float::half two = static_cast<half_float::half>(2);
    half_float::half temp6=half_float::operator/(temp4,two);
    half_float::half temp7=half_float::operator/(temp5,two);
    half_float::half sqrt_real = half_float::sqrt(temp6);
    half_float::half sqrt_imag = half_float::sqrt(temp7);

    // Create output variable
    *real = real_half;
    *image = imag_half;
    *output_real = static_cast<double>(sqrt_real);
    *output_imag = static_cast<double>(sqrt_imag);
}
