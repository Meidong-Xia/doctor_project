function A_fp8 = i_fp8_matrix_turn(A)
A_fp8 = complexFp8ToMatrix_e5m2(complexMatrixToFp8_e5m2(A));
end