function decimalValue = EFPTodec(fp8Binary,config, base, fraction_tables)
    decimalValue = zeros(size(fp8Binary));
    % 获取矩阵的尺寸
    [numRows, numCols] = size(fp8Binary);

    % 将矩阵转换为线性索引
    numElements = numel(fp8Binary);

    % 并行循环逐元素处理
    for idx = 1:numElements
        % [row, col] = ind2sub([numRows, numCols], idx);
        decimalValue(idx) = i_new8Todecimal(fp8Binary{idx}, config{idx}, base, fraction_tables);
    end

    % 将线性索引结果重构为矩阵
    decimalValue = reshape(decimalValue, [numRows, numCols]);
end