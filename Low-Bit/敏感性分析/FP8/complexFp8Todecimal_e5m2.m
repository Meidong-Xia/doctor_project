function decimalValue = complexFp8Todecimal_e5m2(fp8Binary)
    fp8Binary_real = fp8Binary(1,:);
    fp8Binary_imag = fp8Binary(2,:);
    decimalValue_real = fp8Todecimal_e5m2(fp8Binary_real);
    decimalValue_imag = fp8Todecimal_e5m2(fp8Binary_imag);
    decimalValue = decimalValue_real + 1i*decimalValue_imag;
end
