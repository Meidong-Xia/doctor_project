function [result] = GaussianElimination_single(arr)
    [N, ~] = size(arr);
    W = [arr, eye(N)]; % 扩增矩阵
    
    for i = 1:N
        % 如果主元素为0，则交换行
        if W(i, i) == 0
            [~, j] = max(abs(W(i+1:end, i))); % 找到绝对值最大的元素所在的行
            j = j + i; % 在原始矩阵中的行号
            if W(j, i) == 0
                error('矩阵不可逆');
            end
            % 交换行
            temp = W(i, :);
            W(i, :) = W(j, :);
            W(j, :) = temp;
        end

% 将主元素变为1
divisor = single(W(i, i));
W(i, :) = single(W(i, :)) / single(divisor);

        % 将其他行的元素变为0
        for j = 1:N
            if j ~= i
                factor = single(W(j, i));
                if factor ~= 0
                    W(j, :) = single(W(j, :)) - single(factor) * single(W(i, :));
                end
            end
        end
    end
    % 提取逆矩阵
    result = single(W(:, N+1:end));
end
