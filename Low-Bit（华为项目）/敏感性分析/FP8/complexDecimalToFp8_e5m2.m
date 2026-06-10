function [fp8Binary] = complexDecimalToFp8_e5m2(decimalValue)
    decimalValue_real = real(decimalValue);
    decimalValue_imag = imag(decimalValue);
    fp8Binary_real = decimalTofp8_e5m2(decimalValue_real);
    fp8Binary_imag = decimalTofp8_e5m2(decimalValue_imag);
    fp8Binary = [fp8Binary_real;fp8Binary_imag];
end


