function [result_fp8Binary, result_config] = EFP_mul_matrix(fp8Binary1, config1, fp8Binary2, config2, base, fraction_tables,table)
    % %     生成一个 4x4 复数域矩阵
    % decimalValue1 = rand(4) + rand(4) * 1i;
    % decimalValue2 = rand(4) + rand(4) * 1i;
    % % % % 生成一个 1×6 的随机矩阵
    % config_1 = [2 10 2 0 0 2 10 2 0 0];
    % % % % 将其复制到整个 4×4 矩阵中
    % config = repmat({config_1}, 4, 4);
    % % base = repmat(10, 4, 4);
    % % % s = 0;
    % % % 
    % [fp8Binary1, config1] = arrayfun(@i_decimalTonew8_auto, decimalValue1, config,'UniformOutput',false);
    % [fp8Binary2, config2] = arrayfun(@i_decimalTonew8_auto, decimalValue2, config,'UniformOutput',false);

    % 检查矩阵的尺寸是否兼容
    [m1, n1] = size(fp8Binary1);
    [m2, n2] = size(fp8Binary2);
    if n1 ~= m2
        error('矩阵尺寸不兼容，无法相乘。');
    end
    
    % % 初始化结果矩阵
    % global m;
    m = config1{1}(2);
    config_1 = [2 m 2 0 0 2 m 2 0 0];
    config = repmat({config_1},m1, n2);

    % result_fp8Binary = cell(m1, n2);
    % result_config = cell(m1, n2);
    % 计算矩阵乘法
    result = zeros(m1, n2);
    for i = 1:m1
        for j = 1:n2
            % 计算结果矩阵的每个元素
            for k = 1:n1
                if k == 1
                    % [result_fp8Binary_temp{i, j}, result_config_temp{i, j}] = i_EFP_mul(fp8Binary1(i, k), config1(i, k), fp8Binary2(k, j), config2(k, j));
                    [result_fp8Binary_temp, result_config_temp] = i_EFP_mul(fp8Binary1(i, k), config1(i, k), fp8Binary2(k, j), config2(k, j), base, fraction_tables,table);
                    result(i,j) = EFPTodec({result_fp8Binary_temp}, {result_config_temp}, base, fraction_tables);
                else
                    [temp, temp_config] = i_EFP_mul(fp8Binary1(i, k), config1(i, k), fp8Binary2(k, j), config2(k, j), base, fraction_tables,table);
                    temp_dec = EFPTodec({temp}, {temp_config}, base, fraction_tables);
                    result(i,j) = temp_dec + result(i,j);
                    % [result_fp8Binary{i, j}, result_config{i, j}] = i_EFP_add(result_fp8Binary(i, j), result_config(i, j), {temp}, {temp_config});
                end
            end
        end
    end
    [result_fp8Binary,result_config] = decToEFP(result,config, base, fraction_tables);
    % decimalValue = EFPTodec(result_fp8Binary,result_config);
    % decimalValue_sum = decimalValue1 * decimalValue2;
end

