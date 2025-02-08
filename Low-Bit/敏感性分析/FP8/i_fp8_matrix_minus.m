function result_fp8 = i_fp8_matrix_minus(A,B)
A_fp8 = complexFp8ToMatrix_e5m2(complexMatrixToFp8_e5m2(A));
B_fp8 = complexFp8ToMatrix_e5m2(complexMatrixToFp8_e5m2(B));
result = A_fp8-B_fp8;
result_fp8 = complexFp8ToMatrix_e5m2(complexMatrixToFp8_e5m2(result));
end