clear
%% 参数设置
Nl = 648; % 码长
R = 3/4; % 码率
Z = 27; % 提升因子
Nb = Nl*R; % 信息比特数量
Nf = 1e4; % 帧数
M = 64; % M-QAM
maxNumIter = 10;


N_t = 4; % 发送天线数
N_r = N_t; % 接收天线数  

snr_dB = 8:2:20; % SNR
snr = 10.^(snr_dB./10);

delta_dB = 20; % 信道功率对和信道加的总噪声功率的比值（dB）

flg = 5; % 加噪声方式，1：不加；2：全部加；3：加主对角线；4：加第一行；5：加一个元素


%% 生成LDPC相关Object
Hb = [16 17 22 24  9  3 14 -1  4  2  7 -1 26 -1  2 -1 21 -1  1  0 -1 -1 -1 -1
     25 12 12  3  3 26  6 21 -1 15 22 -1 15 -1  4 -1 -1 16 -1  0  0 -1 -1 -1
     25 18 26 16 22 23  9 -1  0 -1  4 -1  4 -1  8 23 11 -1 -1 -1  0  0 -1 -1
      9  7  0  1 17 -1 -1  7  3 -1  3 23 -1 16 -1 -1 21 -1  0 -1 -1  0  0 -1
     24  5 26  7  1 -1 -1 15 24 15 -1  8 -1 13 -1 13 -1 11 -1 -1 -1 -1  0  0
      2  2 19 14 24  1 15 19 -1 21 -1  2 -1 24 -1  3 -1  2  1 -1 -1 -1 -1  0
    ];
pcmatrix = ldpcQuasiCyclicMatrix(Z,Hb);

cfgLDPCEnc = ldpcEncoderConfig(pcmatrix);

cfgLDPCDec = ldpcDecoderConfig(pcmatrix);



% 生成随机信道
% H = (randn(N_r,N_t) + 1i*randn(N_r,N_t)) / sqrt(2);
% save("H1.mat","H");
H = load("H1.mat","H").H;

% norm(H(:,1),2)^2
% norm(H(:,2),2)^2
% norm(H(:,3),2)^2
% norm(H(:,4),2)^2
% 
% norm(H(1,:),2)^2
% norm(H(2,:),2)^2
% norm(H(3,:),2)^2
% norm(H(4,:),2)^2

for m=1:N_t
    for n=1:N_r
        %% 主循环
        error_rate = zeros(length(snr_dB),1);
        
        for i=1:length(snr_dB)
            error_cnt = 0;
            for j=1:Nf
        
                % 生成信号
                data = randi([0,1],Nb,1,'int8');
        
                % 信道编码
                ldpcEncodeData = ldpcEncode(data,cfgLDPCEnc);
        
                % 调制
                moduledData = qammod(ldpcEncodeData,M,'InputType','bit');
        
                % Padding
                paddLength = mod(length(moduledData),N_t);
                moduledDataPadding = [moduledData;moduledData(1:paddLength)];
        
                % 串并转换
                transmittedData = reshape(moduledDataPadding,N_t,[]);
        
                % 测量信号功率
                sigPower = norm(transmittedData,2)^2/numel(transmittedData);
        
                % 预编码
                if flg == 1 % 不加噪声
                    F = H'*inv(H*H');
                    
                elseif flg == 2 % 全部加
                    delta = 1/(10^(delta_dB/10))/numel(H);
                    tmp = (randn(N_r,N_t) + 1i*randn(N_r,N_t)) * sqrt(delta/2);
                    H_head = H + tmp;
                    F = H_head'*inv(H_head*H_head');
                elseif flg == 3 % 加对角线
                    delta = 1/(10^(delta_dB/10))/N_t; 
                    tmp = (randn(N_t,1) + 1i*randn(N_t,1)) * sqrt(delta/2);
                    H_head = H + diag(tmp);
                    F = H_head'*inv(H_head*H_head');
                elseif flg == 4 % 加第一行
                    delta = 1/(10^(delta_dB/10))/N_t; 
                    tmp = (randn(N_t,1) + 1i*randn(N_t,1)) * sqrt(delta/2);
                    H_head = H + [tmp.';zeros(N_r-1,N_t)];
                    F = H_head'*inv(H_head*H_head');
                elseif flg == 5 % 加第一个元素
                    delta = norm(H(m,n))^2/(10^(delta_dB/10));
                    tmp = (randn(1,1) + 1i*randn(1,1)) * sqrt(delta/2);
                    tmp2 = zeros(N_r,N_t);
                    tmp2(m,n) = tmp;
                    H_head = H + tmp2;
                    F = H_head'*inv(H_head*H_head');
                end
                precodedTransmittedData = F*transmittedData;
         
                % 经过信道
                receivedData = H*precodedTransmittedData;
        
                % 串并转换
                receivedData= reshape(receivedData,[],1);
        
        
                % 加噪声
                noiPower = sigPower/snr(i);
                noise = sqrt(noiPower/2)*(randn(numel(receivedData),1) + 1i*randn(numel(receivedData),1));
                receivedDataPlusNoise = receivedData+noise;
        
                % 去Padding
                receivedDataPlusNoise(end-paddLength+1:end) = [];
        
                % 解调
                demoduledData = qamdemod(receivedDataPlusNoise,M,'OutputType','approxllr','NoiseVariance',noiPower);
                
                % 信道解码
                ldpcDecodeData = ldpcDecode(demoduledData,cfgLDPCDec,maxNumIter); 
        
                % 统计错误比特数
                error_cnt = error_cnt + length(find(data~=ldpcDecodeData));
        
            end
            % 计算误比特率
            error_rate(i) = error_cnt/(Nb*Nf);
        end
        
        %% 可视化
        % figure
        % semilogy(snr_dB,error_rate);
        % xlabel('snr/dB');
        % ylabel('ber')
        
        writematrix(error_rate.','mat1.xls','WriteMode','append');

    end
end
    