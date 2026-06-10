function result = fp32Mul_Matrix(matrix1, matrix2)
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
                result(i, j) = result(i, j) + single(single(matrix1(i, l))*single(matrix2(l, j)));
                %result(i, j) = result(i, j) + matrix1(i, l)*matrix2(l, j);
            end
            result(i, j) = single(result(i, j));
        end
    end
end
