function [fp8Binary,result_config] = decToEFP(decimalValue,config, base, fraction_tables)
    % 生成一个 4x4 复数域矩阵
    % decimalValue = rand(4) + rand(4) * 1i;
    % % 生成一个 1×6 的随机矩阵
    % config_1 = [2 5 2 2 5 2];
    % % 将其复制到整个 4×4 矩阵中
    % config = repmat({config_1}, 4, 4);


    % 初始化输出数组
    fp8Binary = cell(size(decimalValue));
    result_config = cell(size(decimalValue));

    % 获取矩阵的尺寸
    [numRows, numCols] = size(decimalValue);

    % 将矩阵转换为线性索引
    numElements = numel(decimalValue);

    % 并行循环逐元素处理
    for idx = 1:numElements
        % [row, col] = ind2sub([numRows, numCols], idx);
        [fp8Binary{idx}, result_config{idx}] = i_decimalTonew8_auto(decimalValue(idx), config(idx), base, fraction_tables);
    end

    % 将线性索引结果重构为矩阵
    fp8Binary = reshape(fp8Binary, [numRows, numCols]);
    result_config = reshape(result_config, [numRows, numCols]);
end