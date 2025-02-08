function [result] = EFP_add_matrix_2(fp81, fp82, m_bit, base, fraction_tables_par,table_par)
    [matrix_efp1, matrix_config1] = decToEFP_auto(fp81,m_bit,base, fraction_tables_par);
    [matrix_efp2, matrix_config2] = decToEFP_auto(fp82,m_bit,base, fraction_tables_par);

    [result_efp, result_config] = EFP_add_matrix(matrix_efp1,matrix_config1,matrix_efp2, ...
        matrix_config2, base, fraction_tables_par,table_par);

    result = EFPTodec(result_efp, result_config,base,fraction_tables_par);
end