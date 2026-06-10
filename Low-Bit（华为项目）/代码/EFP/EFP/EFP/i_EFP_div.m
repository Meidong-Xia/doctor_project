function [result_fp8Binary, result_config] = i_EFP_div(fp8Binary1, config1, fp8Binary2, config2, base, fraction_tables,table)
    % % 做复数二进制的拆解工作
    %     decimalValue1 = rand() + rand() * 1i;
    % decimalValue2 = rand() + rand() * 1i;
    % % % % % 生成一个 1×6 的随机矩阵
    % config_1 = [2 10 2  0 0 2 10 2 0 0];
    % % % % % 将其复制到整个 4×4 矩阵中
    % config = repmat({config_1}, 1, 1);
    % [fp8Binary1, config1] = arrayfun(@i_decimalTonew8_auto, decimalValue1, config,'UniformOutput',false);
    % [fp8Binary2, config2] = arrayfun(@i_decimalTonew8_auto, decimalValue2, config,'UniformOutput',false);

    config1 = config1{1};
    config2 = config2{1};
    fp8Binary1 = fp8Binary1{1};
    fp8Binary2 = fp8Binary2{1};

    fp8Binary1_real = fp8Binary1(1:1+config1(1)+config1(2));
    fp8Binary2_real = fp8Binary2(1:1+config2(1)+config2(2));
    fp8Binary1_imag = fp8Binary1(1+config1(1)+config1(2)+1:end);
    fp8Binary2_imag = fp8Binary2(1+config2(1)+config2(2)+1:end);


    % EFP_add
    % a1*a2
    [temp1_real, temp1_real_config] = EFP_mul(fp8Binary1_real, config1(1:3), fp8Binary2_real, config2(1:3));
    % a1*b2 i
    [temp1_imag, temp1_imag_config] = EFP_mul(fp8Binary1_real, config1(1:3), fp8Binary2_imag, config2(6:8));
    % b1*a2 i
    [temp2_imag, temp2_imag_config] = EFP_mul(fp8Binary1_imag, config1(6:8), fp8Binary2_real, config2(1:3));
    % b1*b2 
    [temp2_real, temp2_real_config] = EFP_mul(fp8Binary1_imag, config1(6:8), fp8Binary2_imag, config2(6:8));

    % a2*a2
    [a22, a22_config] = EFP_mul(fp8Binary2_real, config2(1:3), fp8Binary2_real, config2(1:3));
    % b2* b2
    [b22, b22_config] = EFP_mul(fp8Binary2_imag, config2(6:8), fp8Binary2_imag, config2(6:8));
    
    [denominator, denominator_config] = EFP_add(a22, a22_config, b22, b22_config, base, fraction_tables,table);
   
    [numerator_real, numerator_real_config] = EFP_add(temp1_real, temp1_real_config, temp2_real, temp2_real_config, base, fraction_tables,table);
    [numerator_imag, numerator_imag_config] = EFP_sub(temp2_imag, temp2_imag_config, temp1_imag, temp1_imag_config, base, fraction_tables,table);
 
    [result_real, real_config] = EFP_div(numerator_real, numerator_real_config, denominator, denominator_config);
    [result_imag, imag_config] = EFP_div(numerator_imag, numerator_imag_config, denominator, denominator_config);
   
    result_fp8Binary = [result_real result_imag];
    result_config = [real_config 0 0 imag_config 0 0];

    %     decimalValue = EFPTodec({result_fp8Binary},{result_config});
    % decimalValue_sum = decimalValue1 / decimalValue2;
end