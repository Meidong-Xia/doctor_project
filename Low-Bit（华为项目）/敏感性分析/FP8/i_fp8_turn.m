function a_fp8 = i_fp8_turn(a)
a_fp8 = complexFp8Todecimal_e5m2(complexDecimalToFp8_e5m2(a));
end