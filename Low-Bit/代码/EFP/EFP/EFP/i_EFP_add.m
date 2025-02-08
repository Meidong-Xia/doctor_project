function [result_fp8Binary, result_config] = i_EFP_add(fp8Binary1, config1, fp8Binary2, config2, base, fraction_tables,table)
    % 做复数二进制的拆解工作
    config1 = config1{1};
    config2 = config2{1};
    fp8Binary1 = fp8Binary1{1};
    fp8Binary2 = fp8Binary2{1};

    fp8Binary1_real = fp8Binary1(1:1+config1(1)+config1(2));
    fp8Binary2_real = fp8Binary2(1:1+config2(1)+config2(2));
    fp8Binary1_imag = fp8Binary1(1+config1(1)+config1(2)+1:end);
    fp8Binary2_imag = fp8Binary2(1+config2(1)+config2(2)+1:end);

    % EFP_add
    [result_real, real_config] = EFP_add(fp8Binary1_real, config1(1:3), fp8Binary2_real, config2(1:3), base, fraction_tables,table);
    [result_imag, imag_config] = EFP_add(fp8Binary1_imag, config1(6:8), fp8Binary2_imag, config2(6:8), base, fraction_tables,table);
    
    result_fp8Binary = [result_real result_imag];
    result_config = [real_config 0 0 imag_config 0 0];
    
end