function complexMatrix = fp8Mul_e5m2(b,c)
    [fp8RealMatrix_1, fp8ImagMatrix_1] = complexToFP8Matrix_e5m2(b);
    [fp8RealMatrix_2, fp8ImagMatrix_2] = complexToFP8Matrix_e5m2(c);
    [fp8result_real, fp8result_imag] = i_fp8Mul_e5m2(fp8RealMatrix_1{1}, fp8ImagMatrix_1{1}, fp8RealMatrix_2{1}, fp8ImagMatrix_2{1});
    complexMatrix = fp8MatrixToComplex_e5m2({fp8result_real}, {fp8result_imag});
end