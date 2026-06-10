function [result_fp8Binary, result_config] = i_EFP_sqrt(fp8Binary1, config1, base, fraction_tables,table)
    % 做复数二进制的拆解工作
    % 做复数二进制的拆解工作
    %     decimalValue1 = randn() + randn() * 1i;
    % % % % % 生成一个 1×6 的随机矩阵
    % config_1 = [2 5 2 0 0 2 5 2 0 0];
    % % % % % 将其复制到整个 4×4 矩阵中
    % config = repmat({config_1}, 1, 1);
    % global base;
    % base = 10;
    % [fp8Binary1, config1] = arrayfun(@i_decimalTonew8_auto, decimalValue1, config,'UniformOutput',false);

    config1 = config1{1};
    fp8Binary1 = fp8Binary1{1};

    fp8Binary1_real = fp8Binary1(1:1+config1(1)+config1(2));
    fp8Binary1_imag = fp8Binary1(1+config1(1)+config1(2)+1:end);

    % EFP_add
    % double temp1= real_half*real_half;
    [temp1, temp1_config] = EFP_mul(fp8Binary1_real, config1(1:3), fp8Binary1_real, config1(1:3));
    % double temp2= imag_half*imag_half;
    [temp2, temp2_config] = EFP_mul(fp8Binary1_imag, config1(6:8), fp8Binary1_imag, config1(6:8));
    % double temp3= temp1+temp2;
    [temp3, temp3_config] = EFP_add(temp1, temp1_config, temp2, temp2_config, base, fraction_tables,table);
    % double abs_real_half = sqrt(temp3);
    [abs_real_half, abs_real_half_config] = EFP_sqrt(temp3, temp3_config);
    % double temp4= abs_real_half+real_half;
    [temp4, temp4_config] = EFP_add(abs_real_half, abs_real_half_config, fp8Binary1_real, config1(1:3), base, fraction_tables,table);
    % double temp5= abs_real_half-real_half;
    [temp5, temp5_config] = EFP_sub(abs_real_half, abs_real_half_config, fp8Binary1_real, config1(1:3), base, fraction_tables,table);
    % 
    % double temp6=temp4/2;
    [efp2,efp2_config] = decToEFP(2,{[config1(1:3) 0 0 config1(1:3) 0 0]}, base, fraction_tables);
    efp2 = efp2{1};
    efp2_config = efp2_config{1};

    [temp6, temp6_config] = EFP_div(temp4, temp4_config, efp2(1:1+efp2_config(1)+efp2_config(2)), efp2_config(1:3));
    % double temp7=temp5/2;
    [temp7, temp7_config] = EFP_div(temp5, temp5_config, efp2(1:1+efp2_config(1)+efp2_config(2)), efp2_config(1:3));
    % double sqrt_real_1 = sqrt(temp6);
    [sqrt_real, sqrt_real_config] = EFP_sqrt(temp6, temp6_config);
    % double sqrt_imag_1 = sqrt(temp7);
    % temp7_dec = EFPTodec({[temp7 '0000']},{[temp7_config 0 0 2 1 2 0 0]})
    [sqrt_imag, sqrt_imag_config] = EFP_sqrt(temp7, temp7_config);

    if fp8Binary1_imag(1) == '1' && bin2dec(sqrt_imag) ~= 0
        sqrt_imag(1) = dec2bin(1-bin2dec(sqrt_imag(1)));
    end

    result_fp8Binary = [sqrt_real sqrt_imag];
    result_config = [sqrt_real_config 0 0 sqrt_imag_config 0 0];

    % decimalValue = EFPTodec({result_fp8Binary},{result_config});
    % decimalValue_sum = sqrt(decimalValue1);
     % error = abs(decimalValue - decimalValue_sum) / abs(decimalValue_sum)
end