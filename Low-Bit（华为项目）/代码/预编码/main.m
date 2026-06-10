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

clear
%% 参数设置
Nl = 648; % 码长
R = 3/4; % 码率
blockSize = 27; 
Nb = Nl*R; % 每个流一个块的信息比特数
Nf = 1e3; % 信道实现数
M = 64; % M-QAM
Ns = 8; % 每个用户流数
Nu = 4; % 用户数
maxNumIter = 10;
P = Nu*Ns; % 最大发送功率（每个流单位功率）
kappa = 10; % 莱斯因子(dB)
angles = rand(Nu,2)*pi; % 用户与基站之间主路径的到达角和离开角（[0,pi]内均匀分布）

rho = 0; % 计算误差

resolution = Nf/10; % 调试信息颗粒度


Nt = 64; % 发送天线数>=用户数*流数
Nr = 8; % 接收天线数>=流数
 
snr_dB = -4; % SNR
snr = 10.^(snr_dB./10);

flg_p = 4; % 预编码方式，1：SVD(单用户)，2：ZF，3: BD(多用户)，4：EZF(多用户)

flg_c = 1; % 信道类型，1：瑞利信道，2：莱斯信道

flg_b = 3; % 计算精度，1：64位，2：32位，3：16位

flg_a = 1; % 功率分配，1：平均功率分配，2：注水功率分配，3：数值方法，4：二分查找


%% 生成LDPC相关Object
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
 


%% 主循环
tic
for Ns = [2]
P = Nu*Ns;
Nf = 2*1e4; % 信道实现数
Nblock = 1e2; % 每个信道实现发送的数据块数
resolution = Nf/100; % 调试信息颗粒度
bler = zeros(length(snr_dB),1);
total_run_cnt = 0;
for i=1:length(snr_dB) % SNR
    error_cnt = 0;
    run_cnt = 0;
    for j=1:Nf % 信道实现
        
        % 输出调试信息
        if mod(j,resolution)==0 
            fprintf('当前进度：%G%%\n 耗时：%G小时\n 总块个数：%G\n 有效块个数：%G\n\n', ...
                ((i-1)*Nf+j)/(length(snr_dB)*Nf)*100,toc/3600,((i-1)*Nf+j-1)*Nblock,(total_run_cnt+run_cnt));
        end

        % 生成信道
        if flg_c == 1 % 瑞利信道
            H = (randn(Nu*Nr,Nt) + 1i*randn(Nu*Nr,Nt)) / sqrt(2);
        elseif flg_c == 2 % 莱斯信道
            H = f_channel_generator(Nu,Nr,Nt,kappa,angles);
        end

        % 生成预编码矩阵（对于固定的信道实现，预编码矩阵是一样的）
        if flg_p == 1 % SVD
            [U,Sigma,V,Plist,stop] = f_svd_precoding(H,P,Ns,flg_b,1/snr(i),flg_a,rho);
            
            if stop
                continue;
            end
            eta = 1;
            
        elseif flg_p == 2 % ZF
            [eta,F] = f_zf_precoding(H,P,Nu,Ns,Nr,flg_b);
            

        elseif flg_p == 3 % BD
            [U,Sigma,V,Plist,stop] = f_bd_precoding(H,P,Nu,Ns,Nr,flg_b,1/snr(i),flg_a); % todo
            
            if stop
                continue;
            end
            eta = 1;
        elseif flg_p == 4 % EZF
            [eta,F,U,Sigma,stop] = f_ezf_precoding(H,P,Nu,Ns,Nr,flg_b,1/snr(i),flg_a); 

            Plist = eye(Nu*Ns,Nu*Ns);
            
            if stop
                continue;
            end
             
        end

        % 生成接收端combiner（64bit计算精度）
        if flg_b~=1 && flg_p~=2 
            if flg_p == 1 % SVD
                [U,Sigma,~] = pre_mysvd_lanczos(H);
                U = U(:,1:Ns);
            elseif flg_p == 3 % BD
                % sym = sign(real(U(1,:)));
                [U,Sigma] = f_combiner_bd(Ns,Nu,Nr,H);
                % sym = sym.*sign(real(U(1,:)));
                % U = U*diag(sym);
            elseif flg_p == 4 % EZF
                [U,Sigma] = f_combiner_ezf(Ns,Nu,Nr,H);
            end
        end

        % 每个信道实现发送多个块的数据
        for b = 1:Nblock % 块数
            % 生成信号
            data = randi([0,1],Nb,Ns*Nu);
    
            % 信道编码
            ldpcEncodeData = ldpcEncode(data,cfgLDPCEnc);
    
            % 调制
            moduledData = qammod(ldpcEncodeData,M,'InputType','bit');
            moduledData = moduledData.';
    
            % 测量信号功率、信号功率归一化（单个流）
            sigPower = sum(abs(moduledData).^2,2)/size(moduledData,2);
            moduledData = diag(sqrt(1./sigPower))*moduledData;
    
            
            % 预编码
            if flg_p == 2 % ZF
                tmp = zeros(Nu*Nr,size(moduledData,2));
                for k=1:Nu
                    tmp((k-1)*Nr+1:(k-1)*Nr+Ns,:)=moduledData((k-1)*Ns+1:k*Ns,:);
                end
                precodedData = F*tmp;
            elseif flg_p==4 % EZF
                precodedData = F*Plist*moduledData;
            else % SVD & BD
                precodedData = eta*V*moduledData;
            end
            
     
            % 经过信道
            receivedData = H*precodedData;
    
            % 加噪声
            noise = sqrt(1/(2*snr(i)))*(randn(size(receivedData)) + 1i*randn(size(receivedData)));
            receivedDataPlusNoise = receivedData+noise;
    
            % 恢复信号功率
            if flg_p==1 || flg_p == 3 || flg_p == 4 % SVD || BD || EZF
                receivedDataPlusNoise = inv(Sigma(1:Nu*Ns,1:Nu*Ns)*Plist(1:Nu*Ns,1:Nu*Ns))*U'*receivedDataPlusNoise;
            else % ZF
                for k=Nu:-1:1
                    receivedDataPlusNoise(((k-1)*Nr+Ns+1):(k*Nr),:)=[];
                end
            end
            receivedDataPlusNoise = receivedDataPlusNoise/eta;
            
            receivedDataPlusNoise = diag(sqrt(sigPower))*receivedDataPlusNoise;
    
            % 解调
            receivedDataPlusNoise = receivedDataPlusNoise.';
            demoduledData = qamdemod(receivedDataPlusNoise,M,'OutputType','approxllr');
            
            % 信道解码
            ldpcDecodeData = ldpcDecode(demoduledData,cfgLDPCDec,maxNumIter); 
    
            % 统计错误块和运行块个数
            error_cnt = error_cnt + sum(any(ldpcDecodeData~=data));
            run_cnt = run_cnt+1;
        end
 
    end
    % 统计总运行块
    total_run_cnt = total_run_cnt+run_cnt;
    % 计算bler
    bler(i) = error_cnt/(run_cnt*Ns*Nu);
end

%% 可视化 || 存储
% figure
% semilogy(snr_dB,bler);
% writematrix(bler.','bler.xlsx','WriteMode','append');
bler.'

end


function [U,Sigma] = f_combiner_ezf(Ns,Nu,Nr,H)
U = zeros(Nr*Nu,Ns*Nu);
Sigma = zeros(Ns*Nu,Ns*Nu);
for i = 1:Nu
    H_tmp = H((i-1)*Nr+1:i*Nr,:);
    [U_i,Sigma_i,~] = pre_mysvd_household(H_tmp);
    U((i-1)*Nr+1:i*Nr,(i-1)*Ns+1:i*Ns) = U_i(:,1:Ns);
    Sigma((i-1)*Ns+1:i*Ns,(i-1)*Ns+1:i*Ns) = Sigma_i(1:Ns,1:Ns);
end
end


function [U,Sigma] = f_combiner_bd(Ns,Nu,Nr,H)
    Sigma = zeros(Ns,Ns,Nu);
    U = zeros(Nr,Ns,Nu);
    for k = 1:Nu
        H_tmp = H;
        H_tmp((k-1)*Nr+1:k*Nr,:)=[];
        H_j = H((k-1)*Nr+1:k*Nr,:);
        [~,Sigma_tilde,V_tilde]=pre_mysvd_household(H_tmp);
        V_tilde = V_tilde(:,rank(Sigma_tilde)+1:end);
        [U_j,Sigma_j,~]=pre_mysvd_household(H_j*V_tilde);
        U_j = U_j(:,1:Ns);
        Sigma_j = Sigma_j(1:Ns,1:Ns);
        Sigma(:,:,k)=Sigma_j;
        U(:,:,k) = U_j;
    end
    B=Sigma;
    A=U;
    U = zeros(size(U_j)*Nu);
    Sigma = zeros(size(Sigma_j)*Nu);            
    for k=1:Nu
        U((k-1)*Nr+1:k*Nr,(k-1)*Ns+1:k*Ns)=A(:,:,k);
        Sigma((k-1)*Ns+1:k*Ns,(k-1)*Ns+1:k*Ns)=B(:,:,k);
    end
end
    