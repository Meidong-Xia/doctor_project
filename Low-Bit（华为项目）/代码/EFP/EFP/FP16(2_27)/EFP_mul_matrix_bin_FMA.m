function [result_fp8Binary, result_config] = EFP_mul_matrix_bin_FMA(fp8Binary1, config1, fp8Binary2, config2)
        % 生成一个 4x4 复数域矩阵
    % decimalValue_1 = rand(4) + rand(4) * 1i;
    % decimalValue_2 = rand(4) + rand(4) * 1i;
    % % % % 生成一个 1×6 的随机矩阵
    % config_1 = [2 10 2 0 0 2 10 2 0 0];
    % % % % 将其复制到整个 4×4 矩阵中
    % config = repmat({config_1}, 4, 4);
    % global base;
    % base = 2;
    % % base = repmat(10, 4, 4);
    % % % s = 0;
    % % % 
    % [fp8Binary1, config1] = arrayfun(@i_decimalTonew8_auto, decimalValue_1, config,'UniformOutput',false);
    % [fp8Binary2, config2] = arrayfun(@i_decimalTonew8_auto, decimalValue_2, config,'UniformOutput',false);

    % 检查矩阵的尺寸是否兼容
    [m1, n1] = size(fp8Binary1);
    [m2, n2] = size(fp8Binary2);
    if n1 ~= m2
        error('矩阵尺寸不兼容，无法相乘。');
    end
    
    % 初始化结果矩阵

    result_fp8Binary = cell(m1, n2);
    result_config = cell(m1, n2);
    % 计算矩阵乘法
    for i = 1:m1
        for j = 1:n2
            % 计算结果矩阵的每个元素
            sum = 0;
            for k = 1:n1 
                if k == 1
                    [result_fp8Binary{i, j}, result_config{i, j}] = i_EFP_mul_bin(fp8Binary1(i, k), config1(i, k), fp8Binary2(k, j), config2(k, j));
                    sum = EFPTodec(result_fp8Binary(i, j), result_config(i, j));
                    % decimalValue1 = EFPTodec(fp8Binary1(i, k), config1(i, k));
                    % decimalValue2 = EFPTodec(fp8Binary2(k, j), config2(k, j));
                    % decimalValue = decimalValue1*decimalValue2;
                    % [result_fp8Binary(i, j), result_config(i, j)] = decToEFP(decimalValue,config2(k, j));
                else
                    [temp, temp_config] = i_EFP_mul_bin(fp8Binary1(i, k), config1(i, k), fp8Binary2(k, j), config2(k, j));
                    % decimalValue1 = EFPTodec(fp8Binary1(i, k), config1(i, k));
                    % decimalValue2 = EFPTodec(fp8Binary2(k, j), config2(k, j));
                    % decimalValue = decimalValue1*decimalValue2;
                    % [temp, temp_config] = decToEFP(decimalValue,config2(k, j));

                    % [result_fp8Binary{i, j}, result_config{i, j}] = i_EFP_add(result_fp8Binary(i, j), result_config(i, j), {temp}, {temp_config});
                    % decimalValue1 = EFPTodec(result_fp8Binary(i, j), result_config(i, j));
                    decimalValue2 = EFPTodec({temp}, {temp_config});
                    sum = sum + decimalValue2;
                end
            end
            [result_fp8Binary(i, j), result_config(i, j)] = decToEFP(sum,{temp_config});
        end
    end
    % decimalValue = EFPTodec(result_fp8Binary,result_config);
    % decimalValue_sum = decimalValue_1 * decimalValue_2;
end

