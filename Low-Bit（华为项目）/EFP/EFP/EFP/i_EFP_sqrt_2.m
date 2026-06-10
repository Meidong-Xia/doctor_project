function [result] = i_EFP_sqrt_2(fp81, m_bit, base, fraction_tables_par,table_par)
    [matrix_efp1, matrix_config1] = decToEFP_auto(fp81,m_bit,base, fraction_tables_par);
    [result_efp,result_config] = i_EFP_sqrt(matrix_efp1,matrix_config1, base,fraction_tables_par,table_par);
    result = EFPTodec({result_efp}, {result_config},base,fraction_tables_par);
end