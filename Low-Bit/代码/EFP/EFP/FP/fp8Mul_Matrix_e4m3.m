function result = fp8Mul_Matrix_e4m3(matrix1, matrix2)
    % 检查输入矩阵的维度是否符合矩阵乘法的规则
    [m, n] = size(matrix1);
    [n1, k] = size(matrix2);
    
    if n1 ~= n
        error('矩阵维度不符合矩阵乘法的规则。');
    end

    % 初始化结果矩阵
    result = zeros(m, k);

    % 执行矩阵乘法
    for i = 1:m
        for j = 1:k
            for l = 1:n
                result(i, j) = result(i, j) + fp8Mul_e4m3(matrix1(i, l),matrix2(l, j));
                %result(i, j) = result(i, j) + matrix1(i, l)*matrix2(l, j);
            end
            result_real_bin = decimalTofp8_e4m3(real(result(i, j)));
            result_imag_bin = decimalTofp8_e4m3(imag(result(i, j)));
            result(i, j) = fp8Todecimal_e4m3(result_real_bin)+fp8Todecimal_e4m3(result_imag_bin) * 1i;
        end
    end
end
