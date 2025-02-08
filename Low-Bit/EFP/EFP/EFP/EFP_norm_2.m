function [result] = EFP_norm_2(fp81,m_bit,fraction_tables_par,table_par)
    [M,N] = size(fp81);
    [matrix_efp1, matrix_config1] = decToEFP_auto(fp81,m_bit,base, fraction_tables_par);
    [matrix_efp2, matrix_config2] = hermitian_transpose_2(fp81,m_bit,fraction_tables_par,table_par);
    result = 0;
    for m = 1:M
        [result_efp1, result_config1] = EFP_mul_matrix(matrix_efp1(m,:),matrix_config1(m,:),matrix_efp2(:,m), ...
                matrix_config2(:,m), base, fraction_tables_par,table_par);
        result_efp = EFPTodec(result_efp1, result_config1,base,fraction_tables_par);
        result = result_efp+result;
    end
    
end