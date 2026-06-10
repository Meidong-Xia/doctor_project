#include "mex.h"
#include "half.hpp"
//Half Float Library：这是一个用于 C 语言的库，用于处理 IEEE 754 半精度浮点数。它提供了半精度浮点数的基本算术运算和转换。
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // Check input arguments
    if (nrhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:half_conversion:nrhs", "One inputs required.");
    }
    // Create input variables from MATLAB data
    double* input_a = mxGetPr(prhs[0]);

    // Perform half-precision conversion
    half_float::half a = static_cast<half_float::half>(*input_a);
    half_float::half sqrt=half_float::sqrt(a)	;
    // half_float::half sum = a + b;

    // Create output variables
    plhs[0] = mxCreateDoubleScalar(static_cast<double>(a));
    plhs[1] = mxCreateDoubleScalar(static_cast<double>(sqrt));
}
