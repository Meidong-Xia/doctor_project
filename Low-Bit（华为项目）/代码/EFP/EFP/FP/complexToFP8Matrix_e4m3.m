function [fp8RealMatrix, fp8ImagMatrix] = complexToFP8Matrix_e4m3(complexMatrix)
    % 输入参数 complexMatrix 应为包含十进制复数的二维矩阵
    % 获取矩阵的大小
    [rows, cols] = size(complexMatrix);
    % 初始化实部和虚部的 fp8 二进制矩阵
    fp8RealMatrix = cell(rows, cols);
    fp8ImagMatrix = cell(rows, cols);
    % 循环遍历矩阵元素并进行转换
    for i = 1:rows
        for j = 1:cols
            % 调用之前定义的 decimalTofp8_e5m2 函数，将实部和虚部转换为 fp8 二进制字符串
            realPartBinary = decimalTofp8_e4m3(real(complexMatrix(i, j)));
            imagPartBinary = decimalTofp8_e4m3(imag(complexMatrix(i, j)));

            % 存储 fp8 二进制字符串
            fp8RealMatrix{i, j} = realPartBinary;
            fp8ImagMatrix{i, j} = imagPartBinary;
        end
    end
end
