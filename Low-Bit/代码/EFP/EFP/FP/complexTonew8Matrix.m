function [fp8Matrix] = complexTonew8Matrix(complexMatrix,bias,m)
    % 输入参数 complexMatrix 应为包含十进制复数的二维矩阵
    % 获取矩阵的大小
    [rows, cols] = size(complexMatrix);
    fp8Matrix = zeros(rows,cols);
    % 初始化实部和虚部的 fp8 二进制矩阵
    % 循环遍历矩阵元素并进行转换
    for i = 1:rows
        for j = 1:cols
            % 调用之前定义的 decimalTofp8_e5m2 函数，将实部和虚部转换为 fp8 二进制字符串
            realPartBinary = new8Todecimal_1(decimalTonew8_1(real(complexMatrix(i, j)),bias,m),bias,m);
            imagPartBinary = new8Todecimal_1(decimalTonew8_1(imag(complexMatrix(i, j)),bias,m),bias,m);

            % 存储 fp8 二进制字符串
            fp8Matrix(i, j) = realPartBinary + 1i*imagPartBinary;
        end
    end
end