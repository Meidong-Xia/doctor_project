function [result_efp,result_config] = nx_efp(v_efp,v_config, base, fraction_tables,table)
    % v_efp = r_efp(:,j);
    % v_config = r_config(:,j);
    % global m;
    [her,her_config] = hermitian_transpose(v_efp,v_config);
    [result_efp,result_config] = EFP_mul_matrix(her,her_config,v_efp,v_config, base, fraction_tables,table);
    [result_efp,result_config] = EFP_real(result_efp,result_config);
    % result_efp{1}(m+4:end) = dec2bin(0,m+3);
    % result_config{1}(6:end) = [2 m 2 0 0];
    % i_EFP_sqrt(result_efp,result_config)

    % result = real(EFPTodec(result_efp,result_config)));
    % [result_efp_real,result_config_real] = decimalTonew8_auto(result,result_config(1:3));
    % [result_efp,result_config] = decToEFP(result,result_config);

    [result_efp,result_config] = i_EFP_sqrt(result_efp,result_config, base, fraction_tables,table);
    [result_efp,result_config] = EFP_real({result_efp},{result_config});
    % result_efp = {result_efp};
    % result_config = {result_config};
end