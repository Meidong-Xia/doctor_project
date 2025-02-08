% The implement of several precoding algorithms with different computation precisions
% 
% Based on:
%    
%
% log:
%   - initialized by Meidong Xia on 11/26/2023
%   - added ZF/MMSE precoding scheme by Meidong Xia on 11/30/2023
%   - added BD precoding scheme by Meidong Xia on 12/01/2023
%   - modified receiver combiner to use 64bit computation precision by
%     Meidong Xia on 01/22/2024
%   - added Rician channel by Meidong Xia on 01/22/2024
% 
%
% Rest of the code...

% clear
%% 参数设置
addpath('FP16(2_27)');
addpath('FP');
addpath('SVD');
Nl = 648; % 码长
R = 3/4; % 码率
blockSize = 27; 
Nb = Nl*R; % 每个流一个块的信息比特数
% Nf = 1e1; % 信道实现数
M = 64; % M-QAM
Ns = 8; % 每个用户流数
Nu = 1; % 用户数
maxNumIter = 10;
P = Nu*Ns; % 最大发送功率（每个流单位功率）
kappa = 10; % 莱斯因子(dB)
angles = rand(Nu,2)*pi; % 用户与基站之间主路径的到达角和离开角（[0,pi]内均匀分布）

rho = 0; % 计算误差

% resolution = Nf/10; % 调试信息颗粒度


Nt = 64; % 发送天线数>=用户数*流数
Nr = 8; % 接收天线数>=流数
 
snr_dB = -1:1; % SNR
snr = 10.^(snr_dB./10);

flg_p = 1; % 预编码方式，1：SVD(单用户)，2：ZF，3: BD(多用户)

flg_c = 1; % 信道类型，1：瑞利信道，2：莱斯信道，3：CDL信道

flg_b = 2; % 计算模式，1：FP  2: EFP  3： EFP_VPC

flg_a = 2; % 功率分配，1：平均功率分配，2：最小BER，3：最小化最大MSE，4：最小化最大SINR

% 生成LDPC相关Object
Hb = [16 17 22 24  9  3 14 -1  4  2  7 -1 26 -1  2 -1 21 -1  1  0 -1 -1 -1 -1
     25 12 12  3  3 26  6 21 -1 15 22 -1 15 -1  4 -1 -1 16 -1  0  0 -1 -1 -1
     25 18 26 16 22 23  9 -1  0 -1  4 -1  4 -1  8 23 11 -1 -1 -1  0  0 -1 -1
      9  7  0  1 17 -1 -1  7  3 -1  3 23 -1 16 -1 -1 21 -1  0 -1 -1  0  0 -1
     24  5 26  7  1 -1 -1 15 24 15 -1  8 -1 13 -1 13 -1 11 -1 -1 -1 -1  0  0
      2  2 19 14 24  1 15 19 -1 21 -1  2 -1 24 -1  3 -1  2  1 -1 -1 -1 -1  0
    ];
pcmatrix = ldpcQuasiCyclicMatrix(blockSize,Hb);

cfgLDPCEnc = ldpcEncoderConfig(pcmatrix);

cfgLDPCDec = ldpcDecoderConfig(pcmatrix);

    % global K;
    % global F0;
    % global sigma_u;
    % global x_c_spark;
    % global alpha;
    % global mu_xishu;
    % m = 5;   %平均精度与m无关
    % f0 = m-1;    %决定平均精度的下限
    % mu_a = 0;  %决定K
    % % sigma_a = 0.00017;
    % mu_xishu = 3.5;   %  误差变小，精度边大
    % sigma_xishu = 10; %
    % sigma_a = 4^(-f0)* sigma_xishu;
    % alpha = 2;
    % sigma_u = 0.01;   %
    % F0 = alpha^(-f0);
    % K = (F0 - abs(mu_a)) / sqrt(sigma_a);
    % x_c_spark = 10;
    %%
    
    %%
% global base;
% snr_dB = [-2.2918 -1.442 -0.864]; % 4流
% snr_dB = [-3.22 -2.47 -1.93]; % 2流
% snr_dB = [-3.679 -2.99 -2.465]; % 1流
% snr_dB = [5 6]; % 4流
% snr = 10.^(snr_dB./10);
% snr = [1.0533 1.4874 1.8144];

% flg_p = 1; % 预编码方式，1：SVD(单用户)，2：ZF，3: BD(多用户)130
% flg_b = 1;  % 计算模式，1：FP  2: EFP固定精度  3： EFP_VPC
m_bit_values = [10];  % 指定当前模型下的尾数位宽  若FP  则指定 16 32 64

    % P = Nu * Ns;
    Nf = 1e3;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
    Nblock = 1; % 每个信道实现发送的数据块数
    resolution = Nf / 100; % 调试信息颗粒度

% 创建CDL信道
% num_channels         = Nf;
% H_set = genChannels(num_channels, Nt, Nr);
% load("H_set_B4.mat")
% H_set = sqrt(5)*H_set;
bler_results = zeros(length(snr_dB), length(m_bit_values));
ber_results = zeros(length(snr_dB), length(m_bit_values));
% delete(gcp('nocreate'));
% 创建并行池常量

delete(gcp('nocreate'));
c = parcluster('local');
c.NumWorkers = 14;
parpool(c, c.NumWorkers);


run_cnt = 0;
% fraction_tables_par = fraction_tables;
% table_par = table;
% 主循环 tic;
% val = zeros(Nt,Ns,Nf);
% run = 1;
for m_index = 1:length(m_bit_values)
    m_bit = m_bit_values(m_index);
    fprintf('当前 m_bit 值：%d\n', m_bit);
    

    bler = zeros(length(snr_dB), 1);
    ber = zeros(length(snr_dB), 1);
    total_run_cnt = 0;

    for i = 1:length(snr_dB) % SNR
        % error_cnt = 0;
        % error_bits = 0;
        % run_cnt = 0;
        % run_bits = 0;
            par_error_cnt = zeros(1,Nf);
            par_error_bits = zeros(1,Nf);
            par_run_cnt = zeros(1,Nf);
            par_run_bits = zeros(1,Nf);
        parfor j = 1:Nf % 信道实现
            local_error_cnt = 0;
            local_error_bits = 0;
            local_run_cnt = 0;
            local_run_bits = 0;

            m_bit = 4;  %输入8bit
            if mod(j, resolution) == 0
                fprintf('当前SNR：%d\n 当前信道实现值：%d\n', i,j);
                fprintf('当前进度：%G%%\n 耗时：%G小时\n\n', ...
                    ((i - 1) * Nf + j) / (length(snr_dB) * Nf) * 100, toc / 3600);
                % fprintf('当前进度：%G%%\n 耗时：%G小时\n 总块个数：%G\n 有效块数：%G\n\n', ...
                %     ((i - 1) * Nf + j) / (length(snr_dB) * Nf) * 100, toc / 3600, ...
                %     ((i - 1) * Nf + j - 1) * Nblock,(total_run_cnt+run_cnt));
            end

            % 生成信道
            if flg_c == 1 % 瑞利信道
                H = (randn(Nu * Nr, Nt) + 1i * randn(Nu * Nr, Nt)) / sqrt(2);
            elseif flg_c == 2 % 莱斯信道
                H = f_channel_generator(Nu, Nr, Nt, kappa, angles);
            elseif flg_c ==3 % CDL信道
                % H = squeeze(H_set(j,:,:));
            end

            % 生成预编码矩阵
            if flg_p == 1 % SVD
                if flg_b == 2 % 计算模式，1：FP  2: EFP固定精度  3： EFP_VPC
                        m_bit = 4;  %输入8bit
                        % config_1 = [2 m_bit 2 0 0 2 m_bit 2 0 0];
                        % config = repmat({config_1}, Nu*Nr, Nt);
                        % [H_efp, H_config] = decToEFP(H, config,base,fraction_tables_par);
                        [H_efp, H_config] = decToEFP_auto(H,m_bit,base, fraction_tables_par);
                        H_input = EFPTodec(H_efp, H_config,base,fraction_tables_par);
                        
                        m_bit_jisuan = m_bit_values(m_index);   %计算指定
                        [H_efp, H_config] = decToEFP_auto(H_input,m_bit_jisuan, base, fraction_tables_par);

                        [U,Sigma,V,Plist,stop] = SVD_EFP(H,H_efp, H_config,P,Ns,base,fraction_tables_par, ...
                            table_par,flg_a,1/snr(i),rho);
                        % [U,Sigma,V,Plist,stop] = ...
                        % VPC_subiter_ADSSVD(H,H_efp, H_config,2,15, ...
                        %     m_bit_jisuan,base,fraction_tables_par,table_par,P,Ns,flg_a,1/snr(i),rho);

                 elseif flg_b == 1 % 计算模式，1：FP64
                        
                        [U,Sigma,V,Plist,stop] = SVD_FP(H,P,Ns,flg_b,1/snr(i),flg_a,rho);
                        % acc = norm(H-U*Sigma*V','fro')/norm(H,'fro');
                 end
                    if stop
                        continue;
                    end
                    eta = 1;
            elseif flg_p == 2 % ZF
                if flg_b == 2 % 计算模式，1：FP  2: EFP固定精度  3： EFP_VPC
                    config_1 = [2 m_bit 2 0 0 2 m_bit 2 0 0];
                    config = repmat({config_1}, Nu * Nr, Nt);
                    [H_efp, H_config] = decToEFP(H, config);
                    [eta, F] = ZF_EFP(H_efp, H_config, P, Nu, Ns, Nr);
                elseif flg_b == 1
                    [eta, F] = ZF_FP(H, P, Nu, Ns, Nr, m_bit);
                elseif flg_b == 3
                    config_1 = [2 m_bit 2 0 0 2 m_bit 2 0 0];
                    config = repmat({config_1}, Nu * Nr, Nt);
                    [H_efp, H_config] = decToEFP(H, config);
                    [eta, F] = ZF_vpc(H_efp, H_config, P, Nu, Ns, Nr);
                end

                if any(any(isnan(eta))) || any(any(isnan(F)))
                    continue;
                end
            end

                % 生成接收端combiner（64bit计算精度）   
                % 处理16bit   接收64bit  误差会大
                % 处理16bit   接收16bit  误差会小
                if m_bit~=64 && flg_p~=2  %不是64bit的ZF
                    if flg_p == 1         % SVD   
                        % sym = sign(real(U(1,:)));
                        % [U,Sigma,~] = lansvd(H);
                        % [U,Sigma,~] = cfp32_svd_lanczos(H);
                        U = H*V*inv(Sigma);
                        U = U(:,1:Ns);
                        Sigma = Sigma(1:Ns,1:Ns);
                        % V_2 = inv(Sigma)*U'*H;
                        % sym = sym.*sign(real(U(1,:)));
                        % U = U*diag(sym); % 符号对齐
                        % val(:,:,run) = V-V_2';
                        % run = run+1;
                    elseif flg_p == 3 % BD todo
                        % sym = sign(real(U(1,:)));
                        [U,Sigma] = f_combiner_bd(Ns,Nu,Nr,H);
                        % sym = sym.*sign(real(U(1,:)));
                        % U = U*diag(sym);
                    end
                end

            % 每个信道实现发送多个块的数据
            for b = 1:Nblock % 块数
                data = randi([0, 1], Nb, Ns * Nu); % 生成信号
                ldpcEncodeData = ldpcEncode(data, cfgLDPCEnc); % 信道编码
                moduledData = qammod(ldpcEncodeData, M, 'InputType', 'bit'); % 调制
                moduledData = moduledData.'; % 转置
                sigPower = sum(abs(moduledData) .^ 2, 2) / size(moduledData, 2); % 测量信号功率
                moduledData = diag(sqrt(1 ./ sigPower)) * moduledData; % 信号功率归一化

                % 预编码
                if flg_p == 2 % ZF
                    tmp = zeros(Nu * Nr, size(moduledData, 2));
                    for k = 1:Nu
                        tmp((k - 1) * Nr + 1:(k - 1) * Nr + Ns, :) = moduledData((k - 1) * Ns + 1:k * Ns, :);
                    end
                    precodedData = F * tmp;
                else % SVD & BD
                    precodedData = eta * V* Plist * moduledData;
                end

                receivedData = H * precodedData; % 经过信道
                noise = sqrt(1 / (2 * snr(i))) * (randn(size(receivedData)) + 1i * randn(size(receivedData))); % 加噪声
                receivedDataPlusNoise = receivedData + noise; % 添加噪声

                % 恢复信号功率
                if flg_p == 1 || flg_p == 3 % SVD || BD
                    receivedDataPlusNoise = inv(Sigma(1:Nu * Ns, 1:Nu * Ns) * Plist(1:Nu * Ns, 1:Nu * Ns)) * ...
                    U' * receivedDataPlusNoise;
                else % ZF
                    for k = Nu:-1:1
                        receivedDataPlusNoise(((k - 1) * Nr + Ns + 1):(k * Nr), :) = [];
                    end
                end

                receivedDataPlusNoise = receivedDataPlusNoise / eta;
                receivedDataPlusNoise = diag(sqrt(sigPower)) * receivedDataPlusNoise;
                receivedDataPlusNoise = receivedDataPlusNoise.';
                demoduledData = qamdemod(receivedDataPlusNoise, M, 'OutputType', 'approxllr'); %
                ldpcDecodeData = ldpcDecode(demoduledData, cfgLDPCDec, maxNumIter); % LDPC解码

                local_error_cnt = local_error_cnt + sum(any(ldpcDecodeData ~= data)); % 统计错误块和运行块个数
                local_run_cnt = local_run_cnt + 1;
                local_error_bits = local_error_bits + sum(sum(ldpcDecodeData ~= data)); % 统计错误比特数和总比特数
                local_run_bits = local_run_bits + numel(data);
            end
            par_error_cnt(j) = local_error_cnt;
            par_error_bits(j) = local_error_bits;
            par_run_cnt(j) = local_run_cnt;
            par_run_bits(j) = local_run_bits;
        end
        error_cnt = sum(par_error_cnt);
        error_bits = sum(par_error_bits);
        run_cnt = sum(par_run_cnt);
        run_bits = sum(par_run_bits);
        total_run_cnt = run_cnt + total_run_cnt;
        bler(i) = error_cnt / (run_cnt * Ns * Nu); % 计算bler
        ber(i) = error_bits / run_bits; % 计算BER
    end

    % 将当前 m_bit 值的结果存储到数组中
    bler_results(:, m_index) = bler;
    ber_results(:, m_index) = ber;
end

% val = val(:,:,1:run-1);
% 
% me = mean(val,3);
% 
% va = var(val,0,3);

%%
% 存储 bler 和 ber 到文件中
storage_folder = 'SNR结果';
if ~exist(storage_folder, 'dir')
    mkdir(storage_folder);
end
save(fullfile(storage_folder, 'SVD_EFPM11_all_2(64_8,Nf=10000,Ns=8,SNR=5_6).mat'), 'bler_results', 'ber_results');


%% 画图

% m_bit_values = [16 32 64];  % 指定当前模型下的尾数位宽  若FP  则指定 16 32 64

% 绘制曲线
% figure(1);
% hold on;
% for m_index = 1:length(m_bit_values)
%     plot(snr_dB, ber_results(:, m_index), 'LineWidth', 1.5); % m_bit = 6 对应 ber_results 的第一列
% end

% 添加标签和图例
% xlabel('信噪比 (SNR) (dB)');
% ylabel('误比特率 (BER)');
% if flg_b == 1 % 计算模式，1：FP  2: EFP  3： EFP_VPC
%     title('SVD预编码BER FP');
%     % 构造标签字符串数组
%     labels = cell(1, length(m_bit_values));
%     for i = 1:length(m_bit_values)
%         labels{i} = ['FP' num2str(m_bit_values(i))];
%     end
% elseif flg_b == 2
%     title('SVD预编码 BER EFP固定精度');
%     % 构造标签字符串数组
%     labels = cell(1, length(m_bit_values));
%     for i = 1:length(m_bit_values)
%         labels{i} = ['EFP尾数位宽' num2str(m_bit_values(i))];
%     end
% elseif flg_b == 3
%     title('SVD预编码 BER EFP可变精度');
%     % 构造标签字符串数组
%     labels = cell(1, length(m_bit_values));
%     for i = 1:length(m_bit_values)
%         labels{i} = ['EFP尾数平均位宽' num2str(m_bit_values(1))];
%     end
% end
% legend(labels, 'Location', 'best');
% set(gca, 'YScale', 'log');
% % 显示网格和保持图形可视化
% grid on;
% hold off;
% 
% %% 可变位宽
%     x_add = find(result_add(:, 2) ~= 0);
%     x_sub = find(result_sub(:, 2) ~= 0);
%     x_mul = find(result_mul(:, 2) ~= 0);
%     x_div = find(result_div(:, 2) ~= 0);
%     x_sqrt = find(result_sqrt(:, 2) ~= 0);
%     y_add = result_add(x_add, 2);
%     y_sub = result_sub(x_sub, 2);
%     y_mul = result_mul(x_mul, 2);
%     y_div = result_div(x_div, 2);
%     y_sqrt = result_sqrt(x_div, 2);
% 
%     % 绘制图形
%     figure();
%     hold on;
%     plot(x_add, y_add, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 5); % 添加加法的数据点（蓝色）
%     plot(x_sub, y_sub, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5); % 添加减法的数据点（红色）
%     plot(x_mul, y_mul, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 5); % 添加乘法的数据点（绿色）
%     plot(x_div, y_div, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5); % 添加乘法的数据点（绿色）
%     plot(x_sqrt, y_sqrt, 'autoo', 'MarkerFaceColor', 'auto', 'MarkerSize', 5); % 添加乘法的数据点（绿色）
%     % 标记 f0 和 x_c_spark 的值
% 
%     % 添加标签和标题
%     xlabel('顺序');
%     ylabel('EFP尾数位宽');
%     title('ZF中的位宽变化');
%     legend('加法', '减法', '乘法','除法','开方');
% 
%     % 显示图例和网格
%     legend('Location', 'best');
%     grid on;
%     hold off;
% 
% %% 
% m_bit_values_temp = [6,8,10];
% 
% % 绘制曲线
% figure(1);
% hold on;
% for m_index = 1:length(m_bit_values_temp)
%     plot(snr_dB, bler_results(:, m_index), 'LineWidth', 1.5); % m_bit = 6 对应 ber_results 的第一列
% end
% 
% % 突出显示 FP16 曲线
% plot(snr_dB, bler_results(:,4), 'LineWidth', 2, 'LineStyle', '--', 'Color', 'red'); 
% plot(snr_dB, bler_results(:,5), 'LineWidth', 2, 'LineStyle', '--', 'Color', 'blue'); 
% plot(snr_dB, bler_results(:,6), 'LineWidth', 2, 'LineStyle', '--', 'Color', 'black'); 
% 
% % 构造标签字符串数组
% labels = cell(1, length(m_bit_values_temp) + 1);
% for i = 1:length(m_bit_values_temp)
%     labels{i} = ['尾数位宽' num2str(m_bit_values_temp(i))];
% end
% labels{length(m_bit_values_temp)+1} = 'FP16';
% labels{length(m_bit_values_temp)+2} = 'FP32';
% labels{length(m_bit_values_temp)+3} = 'FP64';
% 
% 
% % 添加标签和图例
% xlabel('信噪比 (SNR) (dB)');
% ylabel('误块率 (BLER)');
% title('ZF预编码 BLER EFP固定精度');
% legend(labels, 'Location', 'best');
% set(gca, 'YScale', 'log');
% % 显示网格和保持图形可视化
% grid on;
% hold off;