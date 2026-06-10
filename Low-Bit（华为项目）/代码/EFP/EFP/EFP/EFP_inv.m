function [result_efp,result_config] = EFP_inv(efp,config1)

    % %     生成一个 4x4 复数域矩阵
    % decimalValue1 = rand(4) + rand(4) * 1i;
    % % % % 生成一个 1×6 的随机矩阵
    % config_1 = [2 10 2 0 0 2 10 2 0 0];
    % % % % 将其复制到整个 4×4 矩阵中
    % config = repmat({config_1}, 4, 4);
    % 
    % %
    % [efp, config1] = arrayfun(@i_decimalTonew8_auto, decimalValue1, config,'UniformOutput',false);
    % 
    
    dec = EFPTodec(efp, config1);
    [N, ~] = size(dec);
    W = [dec, eye(N)]; % 扩增矩阵
% 复数域求逆 部分高斯消元法
% A = i_random_symmetric_positive_definite_matrix(dimension);
    % config_N = repmat(config1(1), N, N);
    % eye_N = decToEFP(eye(N),config_N);
    % W_efp = [efp, eye_N]; % 扩增矩阵
    % W_config = [config1 config_N];
    config_N = repmat(config1(1), N, N);
    [eye_N,config_N_new] = decToEFP(eye(N),config_N);
    W_efp = [efp, eye_N]; % 扩增矩阵
    W_config = [config1,config_N_new];


for i = 1:N
    % 部分主元高斯消元法
    [~, pivot_row] = max(abs(W(i:N, i))); 
    pivot_row = pivot_row + i - 1;

    % 如果主元素为0，则交换行
    % if W(pivot_row, i) == 0
    %     error('矩阵不可逆');
    % end

    % 交换行
    W([i, pivot_row], :) = W([pivot_row, i], :);
    W_efp([i, pivot_row], :) = W_efp([pivot_row, i], :);
    W_config([i, pivot_row], :) = W_config([pivot_row, i], :);

    % 将主元素变为1
    % divisor = W(i, i);
    divisor_efp = W_efp(i, i);
    divisor_config = W_config(i, i);
    % W(i, :) = W(i, :) / divisor;
    [W_efp(i, :),W_config(i, :)] = EFP_div_matrix_e(W_efp(i, :),W_config(i, :),divisor_efp,divisor_config);
    W(i, :) = EFPTodec(W_efp(i, :),W_config(i, :));

    % 将其他行的元素变为0
    for j = 1:N
        if j ~= i
            factor = W(j, i);
            factor_efp = W_efp(j, i);
            factor_config = W_config(j, i);
            if factor ~= 0
%                 W(j, :) = W(j, :) - factor * W(i, :);
                % temp = i_hp_matrix_mul_e(complex(W(i, :)),complex(factor));
                % W(j, :) = i_hp_matrix_sub(complex(W(j, :)),complex(temp));
                [temp_efp,temp_config] = EFP_mul_matrix_e(W_efp(i, :),W_config(i,:),factor_efp,factor_config);
                [W_efp(j, :),W_config(j, :)] = EFP_sub_matrix(W_efp(j, :),W_config(j, :),temp_efp,temp_config);
                W(j, :) = EFPTodec(W_efp(j, :),W_config(j, :));
            end
        end
    end

end

% 提取逆矩阵
result_efp = W_efp(:, N+1:end);
result_config = W_config(:, N+1:end);

% result = EFPTodec(result_efp,result_config);
% true = inv(decimalValue1);

end