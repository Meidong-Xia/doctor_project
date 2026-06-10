function [result, count] = GaussianElimination(arr)
    [N, ~] = size(arr);
    W = [arr, eye(N)]; % 扩增矩阵
    
    % 初始化操作次数
    count = struct('multiplication', 0, 'addition', 0, 'reciprocal', 0);
    
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
divisor = W(i, i);
W(i, :) = W(i, :) / divisor;

% 更新操作次数
count.reciprocal = count.reciprocal + 1;
k = size(W, 2);
    % 检查分子中非零元素的数量，然后增加相应的乘法次数
    num_nonzero = sum(W(i, :) ~= 0);
    count.multiplication = count.multiplication + num_nonzero;

        
        % 将其他行的元素变为0
        for j = 1:N
            if j ~= i
                factor = W(j, i);
                if factor ~= 0
                    W(j, :) = W(j, :) - factor * W(i, :);
                    % 更新操作次数
                    k = size(W, 2);
                    num_nonzero = sum(W(i, :) ~= 0);
                    count.multiplication = count.multiplication + num_nonzero;
                    count.addition = count.addition + k;
                end
            end
        end
    end
    
    % 提取逆矩阵
    result = W(:, N+1:end);
    
%     fprintf('使用高斯消元法矩阵求逆的结果为：\n');
%     disp(result);
end
