clear all
%% 参数设置
Nl = 648; % 码长
R = 3/4; % 码率
blockSize = 27; 
Nb = Nl*R; % 每个流一个块的信息比特数
Nf = 1e4; % 单个流的块数
M = 64; % M-QAM
Ns = 8; % 每个用户流数
Nu = 1; % 用户数
maxNumIter = 10;
P = Nu*Ns; % 最大发送功率（每个流单位功率）
kappa = 5; % 莱斯因子(dB)
angles = rand(Nu,2)*pi; % 用户与基站之间主路径的到达角和离开角（[0,pi]内均匀分布）



resolution = Nf/10; % 调试信息颗粒度


Nt = 64; % 发送天线数>=用户数*流数
Nr = 8; % 接收天线数>=流数
 
snr_dB = -5:10; % SNR
snr = 10.^(snr_dB./10);

flg_p = 1; % 预编码方式，1：SVD(单用户)，2：ZF，3: BD(多用户)

flg_c = 1; % 信道类型，1：瑞利信道，2：莱斯信道

flg_b = 1; % 计算精度，1：64位，2：32位，3：16位

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
 

% save("Sigma",'Sigma');

% load Sigma.mat
rho = 0.001; % 计算误差(方差)

%% 主循环
tic
for Ns = [4]
P = Ns;
Nf = 1e4;
NBlock = 0;
Nb = 486;
ber_sim = zeros(length(snr_dB),1);
ber = zeros(length(snr_dB),1);
ber2 = zeros(length(snr_dB),1);
ber_opt = zeros(length(snr_dB),1);
ber_opt_2 = zeros(length(snr_dB),1);
ber_opt_3 = zeros(length(snr_dB),1);
ber_opt_4 = zeros(length(snr_dB),1);
resolution = ceil(length(snr_dB) /10);
for i= 1:length(snr_dB) 
    % 输出调试信息
    if mod(i,resolution)==0 
        fprintf('当前进度：%G%%\n 耗时：%G小时\n \n\n', ...
            i/length(snr_dB)*100,toc/3600);
    end
    error_cnt = 0;
    run_cnt = 0;
    for c = 1:Nf % 信道实现
        % 生成信道
        if flg_c == 1 % 瑞利信道
            H = (randn(Nu*Nr,Nt) + 1i*randn(Nu*Nr,Nt)) / sqrt(2);
        elseif flg_c == 2 % 莱斯信道
            H = f_channel_generator(Nu,Nr,Nt,kappa,angles);
        end

        % 生成预编码矩阵（对于固定的信道实现，预编码矩阵是一样的）
        [U,Sigma,V]=svd(H);
        V = V(:,1:Ns);
        U = U(:,1:Ns);
        Sigma = Sigma(1:Ns,1:Ns);

        

        % p_ave = ones(Ns,1);
        p_opt = f_power_mse2(M,Ns,Sigma,0,1/snr(i),P);
        p_opt2 = f_power_mse2(M,Ns,Sigma,rho,1/snr(i),P);
        % p_opt_2 = f_power_mse(M,Ns,Sigma,rho,1/snr(i),P);
        % p_opt_3 = f_power_mse2(M,Ns,Sigma,rho,1/snr(i),P);
        % p_opt_4 = f_power_bler(M,Ns,Sigma,rho,1/snr(i),P);
        % p_opt_3 = f_power_blerNP(M,Ns,Sigma,rho,1/snr(i),P);

        % p_opt_3 = ones(Ns,1);
            
        % Plist = diag(p_opt);

        % tmp = 0;

        for j = 1:Ns % 理论BER
            % 计算bler/ber
            % ber(i) = ber(i)+f_calber(Sigma,j,M,rho,P,1/snr(i),p_ave);

            ber2(i) = ber2(i)+f_calber2(Sigma,j,M,rho,P,1/snr(i),p_opt);
    
            % ber_opt(i) = ber_opt(i)+f_calber(Sigma,j,M,rho,P,1/snr(i),p_opt);

            % ber_opt_2(i) = ber_opt_2(i)+f_calber_mse(Sigma,j,M,rho,P,1/snr(i),p_opt_2);

            % ber_opt_3(i) = ber_opt_3(i)+f_calber(Sigma,j,M,rho,P,1/snr(i),p_opt_3);

            % tmp = max(tmp,f_calmse(Sigma,j,M,rho,P,1/snr(i),p_opt_3));
        end

        % ber_opt_3(i) = ber_opt_3(i) + tmp;

        for k = 1:NBlock
            % 生成信号
            data = randi([0,1],Nb,Ns*Nu);
    
            % 信道编码
            % ldpcEncodeData = ldpcEncode(data,cfgLDPCEnc);
            ldpcEncodeData = data;
    
            % 调制
            moduledData = qammod(ldpcEncodeData,M,'InputType','bit');
            moduledData = moduledData.';
    
            % 测量信号功率、信号功率归一化（单个流）
            sigPower = sum(abs(moduledData).^2,2)/size(moduledData,2);
            moduledData = diag(sqrt(1./sigPower))*moduledData;
    
            
            % 预编码
            % E = sqrt(1/2*rho)*(randn(size(V)) + 1i*randn(size(V)));
            % V1 = V + E;
            V1 = V;
            V1 = V1*Plist;
            precodedData = V1*moduledData;
            
     
            % 经过信道
            receivedData = H*precodedData;
    
            % 加噪声
            noise = sqrt(1/(2*snr(i)))*(randn(size(receivedData)) + 1i*randn(size(receivedData)));
            receivedDataPlusNoise = receivedData+noise;
    
            % 恢复信号功率
            % Help = diag(ones(Ns,1));
            % for kk = 1:Ns
            %     Help(kk,kk) = 1+V(:,kk)'*E(:,kk);
            % end
            receivedDataPlusNoise = inv(Sigma(1:Nu*Ns,1:Nu*Ns)*Plist(1:Nu*Ns,1:Nu*Ns))*U'*receivedDataPlusNoise;
            
            receivedDataPlusNoise = diag(sqrt(sigPower))*receivedDataPlusNoise;
    
            % 解调
            receivedDataPlusNoise = receivedDataPlusNoise.';
            demoduledData = qamdemod(receivedDataPlusNoise,M,'OutputType','bit');
            
            % 信道解码
            % ldpcDecodeData = ldpcDecode(demoduledData,cfgLDPCDec,maxNumIter); 
            ldpcDecodeData = demoduledData;
    
            % 统计错误bit和运行块个数
            error_cnt = error_cnt + sum(sum(ldpcDecodeData~=data));
            run_cnt = run_cnt+1;
        end
    end
% ber_sim(i) = error_cnt/(run_cnt*Nb*Ns);
    
end

% ber = ber./Ns./Nf;
ber2 = ber2./Ns./Nf;
% ber_opt = ber_opt./Ns./Nf;
% ber_opt_2 = ber_opt_2./Ns./Nf;
% ber_opt_3 = ber_opt_3./Nf;

%% 可视化 || 存储
% figure
% semilogy(snr_dB,bler);
% writematrix(bler.','bler.xlsx','WriteMode','append');
% bler.'
% ber_sim.'
% ber.'
% ber_opt.'
% ber_opt_2.'
% ber_opt_3.'
end

function p = f_power_blerNP(M,Ns,Sigma,rho,noise,P)
iter = 10;
R = ones(Ns,Ns)*rho/(1+rho);
for i = 1:Ns
    R(i,i) = 0;
end
h = zeros(Ns,1);
for i = 1:Ns
    h(i) = noise/(Sigma(i,i)^2*(1+rho));
end
C = [eye(Ns),zeros(Ns,1);ones(1,Ns),-P];
B = [R,h;zeros(1,Ns),0];
D = inv(C)*B;
[y,lambda] = f_power_iteration(D,iter);
idx = y(end);
y = y/idx;
p = y(1:end-1);
end

function p = f_power_bler(M,Ns,Sigma,rho,noise,P)
p = ones(Ns,1);
f_pre = Inf;
iter = 100;
for i = 1:iter
    y = f_power_bler_1(M,Ns,Sigma,rho,noise,P,p);
    [f,p] = f_power_bler_2(M,Ns,Sigma,rho,noise,P,y);
    if abs(f-f_pre) < 1e-3
        break
    end
    f_pre = f;
end
end

function y = f_power_bler_1(M,Ns,Sigma,rho,noise,P,p)
y = zeros(Ns,1);
for i = 1:Ns
    y(i) = sqrt(p(i)*Sigma(i,i)^2*(1+rho))/((P-p(i))*Sigma(i,i)^2*rho+noise);
end
end

function [f,p] = f_power_bler_2(M,Ns,Sigma,rho,noise,P,y)
cvx_begin quiet
    cvx_solver mosek
    variable p(Ns) nonnegative
    variable f
    maximize(f)
    subject to
        sum(p) == P
        for i = 1:Ns
            2*y(i)*sqrt(p(i)*Sigma(i,i)^2*(1+rho))-y(i)^2*((P-p(i))*Sigma(i,i)^2 ...
                *rho+noise)>=f
        end
cvx_end
end


function p = f_power_mse(M,Ns,Sigma,rho,noise,P)
cvx_begin quiet
    cvx_solver mosek
    variable p(Ns) nonnegative
    variable t
    minimize(t)
    subject to
        sum(p) == P
        for i = 1:Ns
            inv_pos(p(i))*(P*rho+noise/Sigma(i,i)^2)<=t
        end
cvx_end
end

function p = f_power_mse2(M,Ns,Sigma,rho,noise,P) % 二分
u_min = 0;
u_max = 10*P;
eposilon = 1e-5;
p = zeros(Ns,1);
while 1
    u = (u_max+u_min)/2;
    f = 0;
    for i = 1:Ns
        f = f+max((P*rho+noise/Sigma(i,i)^2)/u*P,eposilon);
    end
    if abs(f-P) <= 1e-5
        break
    end
    if f>P
        u_min = u;
    elseif f<P
        u_max = u;
    end
end
for i = 1:Ns
    p(i) = max((P*rho+noise/Sigma(i,i)^2)/u*P,eposilon);
end
end

function p = f_power_ber(M,Ns,Sigma,rho,noise,P)
p = ones(Ns,1);
f_pre = Inf;
iter = 100;
for i = 1:iter
    y = f_power_1(M,Ns,Sigma,rho,noise,P,p);
    [f,p] = f_power_2(M,Ns,Sigma,rho,noise,P,y);
    if abs(f-f_pre) < 1e-5
        break
    end
    f_pre = f;
end
end

function [f,p] = f_power_2(M,Ns,Sigma,rho,noise,P,y)
tmp1 = 2*(1-1/sqrt(M))/log2(M);
cvx_begin quiet
    cvx_solver mosek
    variable p(Ns) nonnegative
    expression a(Ns)
    expression b(Ns)
    expression f
    for i = 1:Ns
        a(i) = p(i)*Sigma(i,i)^2*(1+rho);
        b(i) = (M-1)*(Sigma(i,i)^2*rho*(P-p(i))+noise);
    end
    for i = 1:Ns
        f =  f + tmp1*0.5*(exp(y(i,1)^2-2*y(i,1)*sqrt(3*a(i)))+exp(2*y(i,2)^2 ...
            -2*y(i,2)*sqrt(3*a(i))));
    end
    minimize(f)
    subject to
        sum(p) == P
cvx_end
end

function y = f_power_1(M,Ns,Sigma,rho,noise,P,p)
a = zeros(Ns,1);
b = zeros(Ns,1);
y = zeros(Ns,2);
for i = 1:Ns
    a(i) = p(i)*Sigma(i,i)^2*(1+rho);
    b(i) = (M-1)*(Sigma(i,i)^2*rho*(P-p(i))+noise);
end
for k = 1:Ns
    y(k,1) = sqrt(3*a(k))/b(k);
    y(k,2) = sqrt(3*a(k))/(2*b(k));
end
end

function ber = f_calber(Sigma,j,M,rho,P,noise,p)
    ttt = (3*Sigma(j,j)^2*p(j)*(1+rho))/(2*(M-1)*(Sigma(j,j)^2*rho*(P-p(j))+noise));
    tttt1 = erfc(sqrt(ttt));
    tttt2 = erfc(3*sqrt(ttt));
    ber = (sqrt(M)-1)/(sqrt(M)*log2(sqrt(M)))*tttt1+(sqrt(M)-2)/(sqrt(M)*log2(sqrt(M)))*tttt2;
end


function ber = f_calber2(Sigma,j,M,rho,P,noise,p)
    ttt = (3*Sigma(j,j)^2*p(j)*(1+rho))/(2*(M-1)*(Sigma(j,j)^2*rho*(P-p(j))+noise));
    tttt1 = erfc(sqrt(ttt));
    ber = 2*(1-1/sqrt(M))/log2(M)*tttt1;
end

function ber = f_calber_mse(Sigma,j,M,rho,P,noise,p)
    ttt = (3*Sigma(j,j)^2*p(j))/(2*(M-1)*(Sigma(j,j)^2*rho*P+noise));
    tttt1 = erfc(sqrt(ttt));
    tttt2 = erfc(3*sqrt(ttt));
    ber = (sqrt(M)-1)/(sqrt(M)*log2(sqrt(M)))*tttt1+(sqrt(M)-2)/(sqrt(M)*log2(sqrt(M)))*tttt2;
end
function f = f_calsinr(Sigma,j,M,rho,P,noise,p)
f = (Sigma(j,j)^2*p(j)*(1+rho))/(Sigma(j,j)^2*rho*(P-p(j))+noise);
end

function f = f_calmse(Sigma,j,M,rho,P,noise,p)
f = (P*rho+noise/Sigma(j,j)^2)/p(j);
end


function [u,lambda] = f_power_iteration(H,iter)
[~,n]=size(H);
u = ones(n,1)/sqrt(n);
for i=1:iter
    v = H*u; 
    u = v/(norm(v,2));
end
lambda = (u'*H*u)/(u'*u);
end




