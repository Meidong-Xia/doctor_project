function result_fp8 = i_fp8_minus(a,b)
a_fp8 = complexFp8Todecimal_e5m2(complexDecimalToFp8_e5m2(a));
b_fp8 = complexFp8Todecimal_e5m2(complexDecimalToFp8_e5m2(b));
result = a_fp8-b_fp8;
result_fp8 = complexFp8Todecimal_e5m2(complexDecimalToFp8_e5m2(result));
end