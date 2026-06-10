function [result] = improve_Gaussian_inv_16(arr)
% 复数域求逆 部分高斯消元法
% A = i_random_symmetric_positive_definite_matrix(dimension);
[N, ~] = size(arr);
W = [arr, eye(N)]; % 扩增矩阵


for i = 1:N
    % 部分主元高斯消元法
    [~, pivot_row] = max(abs(W(i:N, i))); 
    pivot_row = pivot_row + i - 1;

    % 如果主元素为0，则交换行
    if W(pivot_row, i) == 0
        error('矩阵不可逆');
    end

    % 交换行
    W([i, pivot_row], :) = W([pivot_row, i], :);

    % 将主元素变为1
    divisor = W(i, i);
%     W(i, :) = W(i, :) / divisor;
    [~, ~, W(i, :)] = hp_matrix_div_e(complex(W(i, :)),complex(divisor));

    % 将其他行的元素变为0
    for j = 1:N
        if j ~= i
            factor = W(j, i);
            if factor ~= 0
%                 W(j, :) = W(j, :) - factor * W(i, :);
                [~, ~, temp] = hp_matrix_mul_e(complex(W(i, :)),complex(factor));
                [~, ~, W(j, :)] = hp_matrix_sub(complex(W(j, :)),complex(temp));
            end
        end
    end
end

% 提取逆矩阵
result = W(:, N+1:end);
end
