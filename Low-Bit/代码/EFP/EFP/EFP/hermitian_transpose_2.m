function [result_efp, result_config] = hermitian_transpose_2(fp81,m_bit,fraction_tables_par,table_par)
    [matrix_efp1, matrix_config1] = decToEFP_auto(fp81,m_bit,base, fraction_tables_par);
    [result_efp, result_config] = hermitian_transpose(matrix_efp1,matrix_config1);
end

