function decimalMatrix = fp8MatrixToDecimal(fp8Matrix)
    % 输入参数 fp8Matrix 应为包含fp8二进制字符串的二维矩阵
    
    % 获取矩阵的大小
    [rows, cols] = size(fp8Matrix);

    % 初始化十进制矩阵
    decimalMatrix = zeros(rows, cols);

    % 循环遍历矩阵元素并进行转换
    for i = 1:rows
        for j = 1:cols
            % 调用之前定义的 fp8ToDecimal 函数
            decimalMatrix(i, j) = fp8ToDecimal(fp8Matrix(i, j));
        end
    end
end
