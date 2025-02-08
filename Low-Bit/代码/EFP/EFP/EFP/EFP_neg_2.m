function [result] = EFP_neg_2(fp81,m_bit,fraction_tables_par,table_par)
    [matrix_efp1, matrix_config1] = decToEFP_auto(fp81,m_bit,base, fraction_tables_par);
    [result_efp, result_config] = EFP_neg(matrix_efp1,matrix_config1);
    result = EFPTodec(result_efp, result_config,base,fraction_tables_par);
end

