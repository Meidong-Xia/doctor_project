function complexMatrix = fp8Mul_e4mX(b,c,m)
    [fp8RealMatrix_1, fp8ImagMatrix_1] = complexToFP8Matrix_e4mX(b,m);
    [fp8RealMatrix_2, fp8ImagMatrix_2] = complexToFP8Matrix_e4mX(c,m);
    decimalValue_real_1 = fp8Todecimal_e4mX(fp8RealMatrix_1{1},m);
    decimalValue_real_2 = fp8Todecimal_e4mX(fp8RealMatrix_2{1},m);
    decimalValue_imag_1 = fp8Todecimal_e4mX(fp8ImagMatrix_1{1},m);
    decimalValue_imag_2 = fp8Todecimal_e4mX(fp8ImagMatrix_2{1},m);
    temp1 = fp8Todecimal_e4mX(decimalTofp8_e4mX(decimalValue_real_1*decimalValue_real_2,m),m);
    temp2 = fp8Todecimal_e4mX(decimalTofp8_e4mX(decimalValue_imag_1*decimalValue_imag_2,m),m);
    fp8result_real = decimalTofp8_e4mX(temp1-temp2,m);
    temp3 = fp8Todecimal_e4mX(decimalTofp8_e4mX(decimalValue_real_1*decimalValue_imag_2,m),m);
    temp4 = fp8Todecimal_e4mX(decimalTofp8_e4mX(decimalValue_real_2*decimalValue_imag_1,m),m);
    fp8result_imag = decimalTofp8_e4mX(temp3+temp4,m);
    complexMatrix = fp8MatrixToComplex_e4mX({fp8result_real}, {fp8result_imag},m);
end