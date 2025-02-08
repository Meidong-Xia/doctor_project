function complexMatrix = fp8Mul_e4m3(b,c)
    [fp8RealMatrix_1, fp8ImagMatrix_1] = complexToFP8Matrix_e4m3(b);
    [fp8RealMatrix_2, fp8ImagMatrix_2] = complexToFP8Matrix_e4m3(c);
    
    [fp8result_real, fp8result_imag] = i_fp8Mul_e4m3(fp8RealMatrix_1{1}, fp8ImagMatrix_1{1}, fp8RealMatrix_2{1}, fp8ImagMatrix_2{1});
    complexMatrix = fp8MatrixToComplex_e4m3({fp8result_real}, {fp8result_imag});
end