function complexMatrix = fp8MatrixToComplex_e3m4(fp8RealMatrix, fp8ImagMatrix,bias)
    % 输入参数 fp8RealMatrix 和 fp8ImagMatrix 应为包含fp8二进制字符串的二维矩阵
    % 获取矩阵的大小
    [rows, cols] = size(fp8RealMatrix);
    % 初始化复数矩阵
    complexMatrix = zeros(rows, cols);
    % 循环遍历矩阵元素并进行转换
    for i = 1:rows
        for j = 1:cols
            % 调用之前定义的 fp8ToDecimal 函数，将实部和虚部转换为十进制
            realPart = fp8Todecimal_e3m4(fp8RealMatrix{i, j},bias);
            imagPart = fp8Todecimal_e3m4(fp8ImagMatrix{i, j},bias);
            % 构建复数
            complexMatrix(i, j) = complex(realPart, imagPart);
        end
    end
end