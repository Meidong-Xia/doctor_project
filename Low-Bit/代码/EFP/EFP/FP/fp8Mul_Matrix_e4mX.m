function result = fp8Mul_Matrix_e4mX(matrix1, matrix2,X)
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
                result(i, j) = result(i, j) + fp8Mul_e4mX(matrix1(i, l),matrix2(l, j),X);
                result_real_bin = decimalTofp8_e4mX(real(result(i, j)),X);
                result_imag_bin = decimalTofp8_e4mX(imag(result(i, j)),X);
                result(i, j) = fp8Todecimal_e4mX(result_real_bin,X)+fp8Todecimal_e4mX(result_imag_bin,X) * 1i;
                %result(i, j) = result(i, j) + matrix1(i, l)*matrix2(l, j);
            end
            result_real_bin = decimalTofp8_e4mX(real(result(i, j)),X);
            result_imag_bin = decimalTofp8_e4mX(imag(result(i, j)),X);
            result(i, j) = fp8Todecimal_e4mX(result_real_bin,X)+fp8Todecimal_e4mX(result_imag_bin,X) * 1i;
        end
    end
end
