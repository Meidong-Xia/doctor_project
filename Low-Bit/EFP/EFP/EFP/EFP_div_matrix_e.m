function [result_fp8Binary, result_config] = EFP_div_matrix_e(fp8Binary1, config1, fp8Binary2, config2, base, fraction_tables,table)
    % %     生成一个 4x4 复数域矩阵
    % decimalValue1 = rand(4) + rand(4) * 1i;
    % decimalValue2 = rand() + rand() * 1i;
    % % % % % 生成一个 1×6 的随机矩阵
    % config_1 = [2 10 2 0 0 2 10 2 0 0];
    % % % % % 将其复制到整个 4×4 矩阵中
    % config = repmat({config_1}, 4, 4);
    % config_2 = repmat({config_1}, 1, 1);
    % global base;
    % base = 10;
    % global mu_xishu;
    % mu_xishu = 0;
    % % base = repmat(10, 4, 4);
    % % % % s = 0;
    % % % % 
    % [fp8Binary1, config1] = arrayfun(@i_decimalTonew8_auto, decimalValue1, config,'UniformOutput',false);
    % [fp8Binary2, config2] = arrayfun(@i_decimalTonew8_auto, decimalValue2, config_2,'UniformOutput',false);
    % 
    % % 检查矩阵的尺寸是否兼容
    % [m1, n1] = size(fp8Binary1);
    % fp8Binary2 = repmat(fp8Binary2, m1, n1);
    % config2 = repmat(config2, m1, n1);
    % 初始化结果矩阵
    % [result_fp8Binary, result_config] = arrayfun(@i_EFP_div,fp8Binary1,config1,fp8Binary2,config2,'UniformOutput',false);
    % 获取矩阵的尺寸
    [numRows, numCols] = size(fp8Binary1);

    % 将矩阵转换为线性索引
    numElements = numel(fp8Binary1);
    result_fp8Binary = cell(numRows, numCols);
    result_config = cell(numRows, numCols);
    % 并行循环逐元素处理
    for idx = 1:numElements
        % [row, col] = ind2sub([numRows, numCols], idx);
        [result_fp8Binary{idx}, result_config{idx}] = i_EFP_div(fp8Binary1(idx), config1(idx),fp8Binary2,config2, base, fraction_tables,table);
    end

    % 将线性索引结果重构为矩阵
    result_fp8Binary = reshape(result_fp8Binary, [numRows, numCols]);
    result_config = reshape(result_config, [numRows, numCols]);
end

