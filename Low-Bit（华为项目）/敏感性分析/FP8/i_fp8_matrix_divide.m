function result_fp8 = i_fp8_matrix_divide(A)
A_fp8 = complexFp8ToMatrix_e5m2(complexMatrixToFp8_e5m2(A));
result = inv(A_fp8);
result_fp8 = complexFp8ToMatrix_e5m2(complexMatrixToFp8_e5m2(result));
end